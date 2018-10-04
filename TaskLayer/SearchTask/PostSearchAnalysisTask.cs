﻿using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Localization;
using EngineLayer.ModificationAnalysis;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class PostSearchAnalysisTask : MetaMorpheusTask
    {
        public PostSearchAnalysisParameters Parameters { get; set; }
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        private IEnumerable<IGrouping<string, PeptideSpectralMatch>> PsmsGroupedByFile { get; set; }

        public PostSearchAnalysisTask()
            : base(MyTask.Search)
        {
        }

        public MyTaskResults Run()
        {
            // Stop loop if canceled
            if (GlobalVariables.StopLoops) { return Parameters.SearchTaskResults; }

            if (Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ModOpen
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Open
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Custom
                )
            {
                // This only makes sense if there is a mass difference that you want to localize. No use for exact and missed monoisotopic mass searches.
                Parameters.SearchParameters.DoLocalizationAnalysis = true;
            }
            else
            {
                Parameters.SearchParameters.DoLocalizationAnalysis = false;
            }

            CalculatePsmFdr();
            DoMassDifferenceLocalizationAnalysis();
            ProteinAnalysis();

            if (Parameters.SearchParameters.DoQuantification)
            {
                QuantificationAnalysis();
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files" }));
            HistogramAnalysis();
            WritePsmResults();
            WriteProteinResults();
            WriteQuantificationResults();
            WritePrunedDatabase();

            return Parameters.SearchTaskResults;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        /// <summary>
        /// Calculate estimated false-discovery rate (FDR) for peptide spectral matches (PSMs)
        /// </summary>
        private void CalculatePsmFdr()
        {
            // TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein, there may be conflicts between how it's handled in parsimony and the FDR engine here
            // for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined by the parsimony algorithm which is agnostic of target/decoy assignments
            // this could cause weird PSM FDR issues

            Status("Estimating PSM FDR...", Parameters.SearchTaskId);
            
            int massDiffAcceptorNumNotches = Parameters.NumNotches;
            var fdrAnalysisResults = (FdrAnalysisResults)(new FdrAnalysisEngine(Parameters.AllPsms, massDiffAcceptorNumNotches, Parameters.CommonParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // sort the list of psms by the score (for FDR analysis)
            // TODO: move this into FdrAnalysisEngine
            if (fdrAnalysisResults.DeltaScoreImprovement)
            {
                Parameters.AllPsms = Parameters.AllPsms.Where(b => b != null)
                    .OrderByDescending(b => b.DeltaScore)
                    .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                    .GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
            }
            else
            {
                Parameters.AllPsms = Parameters.AllPsms.Where(b => b != null)
                    .OrderByDescending(b => b.Score)
                    .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                    .GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
            }

            Status("Done estimating PSM FDR!", Parameters.SearchTaskId);
        }

        private void ProteinAnalysis()
        {
            if (!Parameters.SearchParameters.DoParsimony)
            {
                return;
            }

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(Parameters.AllPsms, Parameters.SearchParameters.ModPeptidesAreDifferent, Parameters.CommonParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, Parameters.AllPsms,
                Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent, true, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            foreach (PeptideSpectralMatch psm in Parameters.AllPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
        }

        private void DoMassDifferenceLocalizationAnalysis()
        {
            if (Parameters.SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running mass-difference localization analysis...", Parameters.SearchTaskId);
                for (int spectraFileIndex = 0; spectraFileIndex < Parameters.CurrentRawFileList.Count; spectraFileIndex++)
                {
                    CommonParameters combinedParams = SetAllFileSpecificCommonParams(Parameters.CommonParameters, Parameters.FileSettingsList[spectraFileIndex]);

                    var origDataFile = Parameters.CurrentRawFileList[spectraFileIndex];
                    Status("Running mass-difference localization analysis...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile });
                    MsDataFile myMsDataFile = Parameters.MyFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                    new LocalizationEngine(Parameters.AllPsms.Where(b => b.FullFilePath.Equals(origDataFile)).ToList(), Parameters.IonTypes,
                        myMsDataFile, combinedParams, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }).Run();
                    Parameters.MyFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }));
                }
            }

            // count different modifications observed
            new ModificationAnalysisEngine(Parameters.AllPsms, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();
        }

        private void QuantificationAnalysis()
        {
            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

            // construct file info for FlashLFQ
            var spectraFileInfo = new List<SpectraFileInfo>();

            // get experimental design info for normalization
            if (Parameters.SearchParameters.Normalize)
            {
                string assumedExperimentalDesignPath = Directory.GetParent(Parameters.CurrentRawFileList.First()).FullName;
                assumedExperimentalDesignPath = Path.Combine(assumedExperimentalDesignPath, GlobalVariables.ExperimentalDesignFileName);

                if (File.Exists(assumedExperimentalDesignPath))
                {
                    var experimentalDesign = File.ReadAllLines(assumedExperimentalDesignPath)
                        .ToDictionary(p => p.Split('\t')[0], p => p);

                    foreach (var file in Parameters.CurrentRawFileList)
                    {
                        string filename = Path.GetFileNameWithoutExtension(file);

                        var expDesignForThisFile = experimentalDesign[filename];
                        var split = expDesignForThisFile.Split('\t');

                        string condition = split[1];
                        int biorep = int.Parse(split[2]);
                        int fraction = int.Parse(split[3]);
                        int techrep = int.Parse(split[4]);

                        // experimental design info passed in here for each spectra file
                        spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file,
                            condition: condition,
                            biorep: biorep - 1,
                            fraction: fraction - 1,
                            techrep: techrep - 1));

                        Parameters.MyFileManager.DoneWithFile(file);
                    }
                }
                else
                {
                    throw new MetaMorpheusException("Could not find experimental design file at location:\n" + assumedExperimentalDesignPath);
                }
            }
            else
            {
                foreach (var file in Parameters.CurrentRawFileList)
                {
                    // experimental design info passed in here for each spectra file
                    spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file, condition: "", biorep: 0, fraction: 0, techrep: 0));
                    Parameters.MyFileManager.DoneWithFile(file);
                }
            }

            // get PSMs to pass to FlashLFQ
            var unambiguousPsmsBelowOnePercentFdr = Parameters.AllPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList();

            var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

            // pass protein group info for each PSM
            var psmToProteinGroups = new Dictionary<PeptideSpectralMatch, List<FlashLFQ.ProteinGroup>>();
            if (ProteinGroups != null)
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    var proteinsOrderedByAccession = proteinGroup.Proteins.OrderBy(p => p.Accession);

                    var flashLfqProteinGroup = new FlashLFQ.ProteinGroup(proteinGroup.ProteinGroupName,
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault())),
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.Organism).Distinct()));

                    foreach (var psm in proteinGroup.AllPsmsBelowOnePercentFDR.Where(v => v.FullSequence != null))
                    {
                        if (psmToProteinGroups.TryGetValue(psm, out var flashLfqProteinGroups))
                        {
                            flashLfqProteinGroups.Add(flashLfqProteinGroup);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { flashLfqProteinGroup });
                        }
                    }
                }
            }
            else
            {
                // if protein groups were not constructed, just use accession numbers
                var accessionToPg = new Dictionary<string, FlashLFQ.ProteinGroup>();
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    var proteins = psm.BestMatchingPeptideWithSetMods.Select(b => b.Pwsm.Protein).Distinct();

                    foreach (var protein in proteins)
                    {
                        if (!accessionToPg.ContainsKey(protein.Accession))
                        {
                            accessionToPg.Add(protein.Accession, new FlashLFQ.ProteinGroup(protein.Accession, string.Join("|", protein.GeneNames.Select(p => p.Item2).Distinct()), protein.Organism));
                        }

                        if (psmToProteinGroups.TryGetValue(psm, out var proteinGroups))
                        {
                            proteinGroups.Add(accessionToPg[protein.Accession]);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { accessionToPg[protein.Accession] });
                        }
                    }
                }
            }

            // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
            // the peptides should still be quantified but not considered for protein quantification
            var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
            foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
            {
                if (!psmToProteinGroups.ContainsKey(psm))
                {
                    psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { undefinedPg });
                }
            }

            // pass PSM info to FlashLFQ
            var flashLFQIdentifications = new List<Identification>();
            foreach (var spectraFile in psmsGroupedByFile)
            {
                var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

                foreach (var psm in spectraFile)
                {
                    flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
                        psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
                }
            }

            // run FlashLFQ
            var FlashLfqEngine = new FlashLFQEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: Parameters.SearchParameters.Normalize,
                ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
                matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
                silent: true,
                optionalPeriodicTablePath: GlobalVariables.ElementsLocation);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = FlashLfqEngine.Run();
            }

            // get protein intensity back from FlashLFQ
            if (ProteinGroups != null && Parameters.FlashLfqResults != null)
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    proteinGroup.FilesForQuantification = spectraFileInfo;
                    proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                    foreach (var spectraFile in proteinGroup.FilesForQuantification)
                    {
                        if (Parameters.FlashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
                        {
                            proteinGroup.IntensitiesByFile.Add(spectraFile, flashLfqProteinGroup.GetIntensity(spectraFile));
                        }
                        else
                        {
                            if (proteinGroup.AllPsmsBelowOnePercentFDR.Count(v => v.FullSequence != null) > 0)
                            {
                                //TODO This if clase didn't exist before.
                                if (!proteinGroup.IsDecoy)
                                    throw new MetaMorpheusException("Could not get intensity back from FlashLFQ! " + proteinGroup.ProteinGroupName);
                            }
                            else
                            {
                                // this is pretty rare... PSMs below 1% FDR were observed for this protein in this file, but they were all ambiguous PSMs
                                // in this case, set the protein quantity to 0 for this file (FlashLFQ didn't quantify it because there were no valid PSMs)
                                proteinGroup.IntensitiesByFile.Add(spectraFile, 0);
                            }
                        }
                    }
                }
            }
        }

        private void HistogramAnalysis()
        {
            if (Parameters.SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = Parameters.AllPsms.Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                {
                    Status("Running histogram analysis...", new List<string> { Parameters.SearchTaskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr, Parameters.SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(Parameters.OutputFolder, "aggregate.mytsv");
                    WriteTree(myTreeStructure, writtenFile);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        private void WritePsmResults()
        {
            Status("Writing results...", Parameters.SearchTaskId);
            if (Parameters.CurrentRawFileList.Count > 1)
            {
                var writtenFile = Path.Combine(Parameters.OutputFolder, "aggregatePSMs.psmtsv");
                WritePsmsToTsv(Parameters.AllPsms, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                SucessfullyFinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

                var writtenFileForPercolator = Path.Combine(Parameters.OutputFolder, "forPercolator.tsv");
                WritePsmsForPercolator(Parameters.AllPsms, writtenFileForPercolator);
                SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { Parameters.SearchTaskId });
            }
            Parameters.SearchTaskResults.AddNiceText("All target PSMS within 1% FDR: " + Parameters.AllPsms.Count(a => a.FdrInfo.QValue < .01 && !a.IsDecoy));

            var peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            if (Parameters.CurrentRawFileList.Count > 1)
            {
                var writtenFile = Path.Combine(Parameters.OutputFolder, "aggregatePeptides.psmtsv");
                WritePsmsToTsv(peptides, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                SucessfullyFinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
            }
            Parameters.SearchTaskResults.AddNiceText("All target peptides within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue < 0.01 && !a.IsDecoy));

            if (Parameters.SearchParameters.DoParsimony)
            {
                Parameters.SearchTaskResults.AddNiceText("All target protein groups within 1% FDR: " + ProteinGroups.Count(b => b.QValue < 0.01 && !b.IsDecoy) + Environment.NewLine);
            }

            PsmsGroupedByFile = Parameters.AllPsms.GroupBy(p => p.FullFilePath);

            // individual psm files (with global psm fdr, global parsimony)
            foreach (var group in PsmsGroupedByFile) //just spectra
            {
                var psmsForThisFile = group.ToList();
                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);
                Parameters.SearchTaskResults.AddNiceText("MS2 spectra in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][0]);
            }
            foreach (var group in PsmsGroupedByFile) //just fragmented precursors
            {
                var psmsForThisFile = group.ToList();
                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);
                Parameters.SearchTaskResults.AddNiceText("Precursors fragmented in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][1]);
            }
            foreach (var group in PsmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                var writtenFile = Path.Combine(Parameters.OutputFolder, strippedFileName + "_PSMs.psmtsv");
                WritePsmsToTsv(psmsForThisFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                SucessfullyFinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", group.First().FullFilePath });
                Parameters.SearchTaskResults.AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));

                var writtenFileForPercolator = Path.Combine(Parameters.OutputFolder, strippedFileName + "_forPercolator.tsv");
                WritePsmsForPercolator(psmsForThisFile, writtenFileForPercolator);
                SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", group.First().FullFilePath });
            }
            foreach (var group in PsmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                var writtenFile = Path.Combine(Parameters.OutputFolder, strippedFileName + "_Peptides.psmtsv");
                WritePsmsToTsv(peptidesForFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                SucessfullyFinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", group.First().FullFilePath });
                Parameters.SearchTaskResults.AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));
            }
        }

        private void WriteProteinResults()
        {
            if (Parameters.SearchParameters.DoParsimony)
            {
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    WriteProteinGroupsToTsv(ProteinGroups, Parameters.OutputFolder, "aggregateProteinGroups", new List<string> { Parameters.SearchTaskId }, PsmsGroupedByFile.Select(b => b.Key).ToList());
                }

                // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                foreach (var fullFilePath in Parameters.CurrentRawFileList)
                {
                    List<PeptideSpectralMatch> psmsForThisFile = PsmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();

                    var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                    var subsetProteinGroupsForThisFile = new List<EngineLayer.ProteinGroup>();
                    foreach (var pg in ProteinGroups)
                    {
                        subsetProteinGroupsForThisFile.Add(pg.ConstructSubsetProteinGroup(fullFilePath));
                    }

                    var subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                        Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent,
                        false, CommonParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();
                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                    Parameters.SearchTaskResults.AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue < 0.01 && !b.IsDecoy));

                    WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, Parameters.OutputFolder, strippedFileName + "_ProteinGroups", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }, new List<string> { fullFilePath });
                    if (Parameters.SearchParameters.OutMzId)
                    {
                        Status("Writing mzid...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        var mzidFilePath = Path.Combine(Parameters.OutputFolder, strippedFileName + ".mzid");
                        MzIdentMLWriter.WriteMzidentml(psmsForThisFile, subsetProteinGroupsForThisFile, Parameters.VariableModifications, Parameters.FixedModifications, new List<Protease> { Parameters.CommonParameters.DigestionParams.Protease }, 0.01, Parameters.CommonParameters.ProductMassTolerance, Parameters.CommonParameters.PrecursorMassTolerance, Parameters.CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);
                        SucessfullyFinishedWritingFile(mzidFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                    }
                    if (Parameters.SearchParameters.OutPepXML)
                    {
                        Status("Writing pep.XML...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        var pepXMLFilePath = Path.Combine(Parameters.OutputFolder, strippedFileName + ".pep.XML");
                        PepXMLWriter.WritePepXML(psmsForThisFile, Parameters.DatabaseFilenameList, Parameters.VariableModifications, Parameters.FixedModifications, Parameters.CommonParameters, pepXMLFilePath);
                        SucessfullyFinishedWritingFile(pepXMLFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                    }
                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
                }
            }
        }

        private void WriteQuantificationResults()
        {
            if (Parameters.SearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                foreach (var file in Parameters.FlashLfqResults.Peaks)
                {
                    WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                }

                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "aggregateQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });
                }

                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "aggregateQuantifiedPeptides", new List<string> { Parameters.SearchTaskId });
            }
        }

        private void WritePrunedDatabase()
        {
            if (Parameters.SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { Parameters.SearchTaskId });
                List<Modification> modificationsToWriteIfBoth = new List<Modification>();
                List<Modification> modificationsToWriteIfInDatabase = new List<Modification>();
                List<Modification> modificationsToWriteIfObserved = new List<Modification>();

                var confidentPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.BaseSequence != null).ToList();
                var proteinToConfidentBaseSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
                foreach (PeptideSpectralMatch psm in confidentPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptideWithSetMods.Select(p => p.Pwsm);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentBaseSequences.TryGetValue(peptide.Protein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentBaseSequences.Add(peptide.Protein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                // Add user mod selection behavours to Pruned DB
                foreach (var modType in Parameters.SearchParameters.ModsToWriteSelection)
                {
                    if (modType.Value == 1) // Write if observed and in database
                    {
                        modificationsToWriteIfBoth.AddRange(GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)));
                    }
                    if (modType.Value == 2) // Write if in database
                    {
                        modificationsToWriteIfInDatabase.AddRange(GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)));
                    }
                    if (modType.Value == 3) // Write if observed
                    {
                        modificationsToWriteIfObserved.AddRange(GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)));
                    }
                }

                //generates dictionary of proteins with only localized modifications
                var ModPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();
                var proteinToConfidentModifiedSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                foreach (PeptideSpectralMatch psm in ModPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptideWithSetMods.Select(p => p.Pwsm);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentModifiedSequences.TryGetValue(peptide.Protein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentModifiedSequences.Add(peptide.Protein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
                foreach (var protein in Parameters.ProteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, Modification>> modsObservedOnThisProtein = new HashSet<Tuple<int, Modification>>();
                        if (proteinToConfidentModifiedSequences.ContainsKey(protein))
                        {
                            modsObservedOnThisProtein = new HashSet<Tuple<int, Modification>>(proteinToConfidentModifiedSequences[protein].SelectMany(b => b.AllModsOneIsNterminus.Select(c => new Tuple<int, Modification>(GetOneBasedIndexInProtein(c.Key, b), c.Value))));
                        }

                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();

                        foreach (var observedMod in modsObservedOnThisProtein)
                        {
                            //Add if observed (regardless if in database)
                            var tempMod = observedMod.Item2;

                            if (modificationsToWriteIfObserved.Contains(tempMod as Modification))
                            {
                                if (!modsToWrite.ContainsKey(observedMod.Item1))
                                {
                                    modsToWrite.Add(observedMod.Item1, new List<Modification> { observedMod.Item2 as Modification });
                                }
                                else
                                {
                                    modsToWrite[observedMod.Item1].Add(observedMod.Item2 as Modification);
                                }
                                continue;
                            }
                        }

                        // Add if in database (two cases: always or if observed)
                        foreach (var modkv in protein.OneBasedPossibleLocalizedModifications)
                        {
                            foreach (var mod in modkv.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod as Modification) ||
                                    (modsObservedOnThisProtein.Contains(new Tuple<int, Modification>(modkv.Key, mod as Modification))
                                        && modificationsToWriteIfBoth.Contains(mod as Modification)))
                                {
                                    if (!modsToWrite.ContainsKey(modkv.Key))
                                    {
                                        modsToWrite.Add(modkv.Key, new List<Modification> { mod });
                                    }
                                    else
                                    {
                                        modsToWrite[modkv.Key].Add(mod);
                                    }
                                }
                            }
                        }

                        if (proteinToConfidentBaseSequences.TryGetValue(protein, out var peptideSequences))
                        {
                            // removes all annotated mods on proteins
                            if (protein.Accession == protein.Accession)
                            {
                                protein.OneBasedPossibleLocalizedModifications.Clear();
                                // adds confidently localized and identified mods
                                foreach (var kvp in modsToWrite)
                                {
                                    protein.OneBasedPossibleLocalizedModifications.Add(kvp);
                                }
                            }
                        }
                    }
                }

                //writes all proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }

                //writes only detected proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
            {
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            }
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
            {
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            }
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
        }

        private static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.Combos
                        + "\t" + string.Join(",", bin.ModsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.ResidueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.UniprotID);
                }
            }
        }

        private void WritePsmsForPercolator(List<PeptideSpectralMatch> psmList, string writtenFileForPercolator)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
                output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
                for (int i = 0; i < psmList.Count; i++)
                {
                    PeptideSpectralMatch psm = psmList[i];

                    output.Write(i.ToString());
                    output.Write('\t' + (psm.IsDecoy ? -1 : 1).ToString());
                    output.Write('\t' + psm.ScanNumber.ToString());

                    // Features
                    output.Write('\t' + string.Join("\t", psm.Features));

                    // HACKY: Ignores all ambiguity
                    var pwsm = psm.BestMatchingPeptideWithSetMods.First().Pwsm;

                    output.Write('\t' + (pwsm.PreviousAminoAcid + "." + pwsm.FullSequence + "." + pwsm.NextAminoAcid).ToString());
                    output.Write('\t' + (pwsm.Protein.Accession).ToString());
                    output.WriteLine();
                }
            }
        }

        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> items, string outputFolder, string strippedFileName, List<string> nestedIds, List<string> FileNames)
        {
            if (items != null && items.Any())
            {
                var writtenFile = Path.Combine(outputFolder, strippedFileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(items.First().GetTabSeparatedHeader());
                    for (int i = 0; i < items.Count; i++)
                    {
                        output.WriteLine(items[i]);
                    }
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeptideQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var baseSeqPath = Path.Combine(outputFolder, fileName + "BaseSequences.tsv");
            var fullSeqPath = Path.Combine(outputFolder, fileName + "FullSequences.tsv");

            flashLFQResults.WriteResults(null, fullSeqPath, baseSeqPath, null);

            SucessfullyFinishedWritingFile(baseSeqPath, nestedIds);
            SucessfullyFinishedWritingFile(fullSeqPath, nestedIds);
        }

        private void WritePeakQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(peaksPath, null, null, null);

            SucessfullyFinishedWritingFile(peaksPath, nestedIds);
        }
    }
}