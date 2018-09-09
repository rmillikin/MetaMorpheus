using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Text;
using EngineLayer.FdrAnalysis;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected static readonly Dictionary<DissociationType, double> complementaryIonConversionDictionary = new Dictionary<DissociationType, double>
        {
            { DissociationType.HCD, Constants.ProtonMass },
            { DissociationType.ETD, 2*Constants.ProtonMass }
        };

        protected readonly CommonParameters commonParameters;

        protected readonly List<string> nestedIds;

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<string> nestedIds)
        {
            this.commonParameters = commonParameters;
            this.nestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        private static Dictionary<string, Modification> fds = new Dictionary<string, Modification>();
        public static Dictionary<int, IsotopicEnvelope[]> scanToEnvelopes = new Dictionary<int, IsotopicEnvelope[]>();

        public static void DeconvoluteAndStoreAllMs2(Ms2ScanWithSpecificMass[] ms2Scans, CommonParameters commonParameters)
        {
            double ms2DeconvolutionPpmTolerance = 5.0;
            int minZ = 1;
            int maxZ = 10;

            foreach (var scan in ms2Scans)
            {
                // deconvolute the scan
                if (!scanToEnvelopes.ContainsKey(scan.OneBasedScanNumber))
                {
                    var isotopicEnvelopes = scan.TheScan.MassSpectrum.Deconvolute(scan.TheScan.MassSpectrum.Range, minZ, maxZ,
                        ms2DeconvolutionPpmTolerance, commonParameters.DeconvolutionIntensityRatio).ToArray();

                    scanToEnvelopes.Add(scan.OneBasedScanNumber, isotopicEnvelopes);
                }
            }
        }

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
        {
            // scoring if some fragments get doubled for scoring purposes
            if (maximumMassThatFragmentIonScoreIsDoubled > 0)
            {
                double score = 0;

                foreach (var fragment in matchedFragmentIons)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);

                    if (fragment.NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled) // TODO: not sure if this is supposed to be neutral mass or mz
                    {
                        score += fragmentScore * 2;
                    }
                    else
                    {
                        score += fragmentScore;
                    }
                }

                return score;
            }

            // normal scoring
            return matchedFragmentIons.Count + (matchedFragmentIons.Sum(v => v.Intensity) / thisScan.TotalIonCurrent);
        }

        public static List<MatchedFragmentIon> MatchFragmentIons(MsDataScan spectrum, List<Product> theoreticalProducts, CommonParameters commonParameters, bool deconvolutionSearch = false)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var alreadyCountedMzs = new HashSet<double>();

            // if the spectrum has no peaks
            if (spectrum.MassSpectrum.Size == 0)
            {
                return matchedFragmentIons;
            }

            if (deconvolutionSearch)
            {
                return MatchFragmentIonsDeconvolution(spectrum, theoreticalProducts, commonParameters);
            }

            //search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak assuming z=1
                int matchedPeakIndex = spectrum.MassSpectrum.GetClosestPeakIndex(product.NeutralMass.ToMz(1)).Value;

                double mz = spectrum.MassSpectrum.XArray[matchedPeakIndex];
                double intensity = spectrum.MassSpectrum.YArray[matchedPeakIndex];

                // is the mass error acceptable and has it been counted already?
                if (commonParameters.ProductMassTolerance.Within(mz, product.NeutralMass.ToMz(1)) && !alreadyCountedMzs.Contains(mz))
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, intensity, 1));
                    alreadyCountedMzs.Add(mz);
                }
            }

            return matchedFragmentIons;
        }

        public static MzSpectrum GenerateComplementarySpectrum(MzSpectrum spectrum, double precursorMass, DissociationType dissociationType)
        {
            double protonMassShift = complementaryIonConversionDictionary[dissociationType].ToMass(1);

            double[] newMzSpectrum = new double[spectrum.Size];
            double[] intensity = new double[spectrum.Size];

            for (int i = spectrum.Size - 1; i >= 0; i--)
            {
                int j = spectrum.Size - i - 1;

                double mz = spectrum.XArray[i];
                double compFragmentMass = (precursorMass + protonMassShift) - mz.ToMass(1); //FIXME, not valid for all fragmentation (b+y+H = precursor, but c+zdot+2H = precursor)

                newMzSpectrum[j] = compFragmentMass.ToMz(1);
                intensity[j] = spectrum.YArray[i];
            }

            return new MzSpectrum(newMzSpectrum, intensity, false);
        }

        public MetaMorpheusEngineResults Run()
        {
            StartingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            FinishedSingleEngine(myResults);
            return myResults;
        }

        public string GetId()
        {
            return string.Join(",", nestedIds);
        }

        public static List<DissociationType> DetermineDissociationType(List<ProductType> lp)
        {
            List<DissociationType> dissociationTypes = new List<DissociationType>();

            if (lp.Contains(ProductType.b) || lp.Contains(ProductType.y))
            {
                dissociationTypes.Add(DissociationType.HCD);
            }

            if (lp.Contains(ProductType.c) || lp.Contains(ProductType.zPlusOne))
            {
                dissociationTypes.Add(DissociationType.ETD);
            }

            return dissociationTypes;
        }

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MetaMorpheusEngineResults RunSpecific();

        public static void DoNewFdr(PeptideSpectralMatch psm, MsDataScan scan, CommonParameters commonParameters)
        {
            int numDatabases = 20;
            var peptide = psm.BestMatchingPeptideWithSetMods.First().Pwsm;
            bool deconScoring = commonParameters.DigestionParams.Protease.Name == "top-down";

            double numDecoys = 0;
            for (int i = 0; i < numDatabases; i++)
            {
                string randomizedSequence = ShuffleString(peptide.BaseSequence);

                while (randomizedSequence == peptide.BaseSequence)
                {
                    randomizedSequence = ShuffleString(peptide.BaseSequence);
                }

                PeptideWithSetModifications randomizedPeptide = new PeptideWithSetModifications(randomizedSequence, fds, peptide.NumFixedMods, peptide.DigestionParams, peptide.Protein, int.MinValue, int.MinValue, peptide.MissedCleavages, peptide.PeptideDescription);
                var theorFragments = randomizedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
                var matchedFragments = MatchFragmentIons(scan, theorFragments, commonParameters, deconScoring);

                double decoyScore = CalculatePeptideScore(scan, matchedFragments, 0);
                if (decoyScore >= psm.Score)
                {
                    numDecoys++;
                }
                //if (commonParameters.CalculateEValue)
                //{
                psm.AllScores.Add(decoyScore);
                //}
            }

            double percentDecoy = numDecoys / numDatabases;
            psm.percentDecoy = percentDecoy;
        }

        private static string ShuffleString(string start)
        {
            Random r = new Random();
            int startLength = start.Length;
            StringBuilder s = new StringBuilder();
            s.Append(start.First());
            var shuffleablePart = start.Substring(1, start.Length - 2).ToArray();

            int ind = r.Next(0, shuffleablePart.Length - 1);

            while (shuffleablePart.Any(p => p != '-'))
            {
                ind = r.Next(0, shuffleablePart.Length);

                if (shuffleablePart[ind] != '-')
                {
                    s.Append(shuffleablePart[ind]);
                }
                shuffleablePart[ind] = '-';
            }
            s.Append(start.Last());

            int endLength = s.Length;
            return s.ToString();
        }

        private static void fff()
        {
            foreach (var mod in GlobalVariables.AllModsKnown)
            {
                if (!fds.ContainsKey(mod.IdWithMotif))
                {
                    fds.Add(mod.IdWithMotif, mod);
                }
            }
        }


        private static List<MatchedFragmentIon> MatchFragmentIonsDeconvolution(MsDataScan spectrum, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var isotopicEnvelopes = scanToEnvelopes[spectrum.OneBasedScanNumber];

            // return empty list of matched fragments if there are no deconvoluted isotopic envelopes
            if (!isotopicEnvelopes.Any())
            {
                return matchedFragmentIons;
            }

            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the first isotopic envelope within the ppm tolerance
                // if there is no envelope within the desired ppm tolerance, "bestEnvelope" will be null
                IsotopicEnvelope bestEnvelope = isotopicEnvelopes.FirstOrDefault(p =>
                    commonParameters.ProductMassTolerance.Within(p.monoisotopicMass, product.NeutralMass));

                // add the matched fragment
                if (bestEnvelope != null)
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product,
                        bestEnvelope.monoisotopicMass.ToMz(bestEnvelope.charge),
                        bestEnvelope.totalIntensity, bestEnvelope.charge));
                }
            }

            return matchedFragmentIons;
        }

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }
    }
}