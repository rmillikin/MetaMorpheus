﻿using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

/// <summary>
/// The mass difference between an experimental precursor and the theoretical mass of the assigned peptide is determined. The LocalizationEngine attempts
/// to localize this mass to one of the residues. It does this by adding the mass difference to each theoretical ion mass and looking for additional matches
/// in the experimental spectrum. This engine should only be run for open, notch or custom searches. It should not be run for exact mass or missed
/// monoisopic searches.
/// </summary>

namespace EngineLayer.Localization
{
    public class LocalizationEngine : MetaMorpheusEngine
    {
        private readonly IEnumerable<PeptideSpectralMatch> AllResultingIdentifications;
        private readonly List<ProductType> ProductTypes;
        private readonly MsDataFile MyMsDataFile;


        public LocalizationEngine(IEnumerable<PeptideSpectralMatch> allResultingIdentifications, List<ProductType> lp, MsDataFile myMsDataFile, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            AllResultingIdentifications = allResultingIdentifications;
            ProductTypes = lp;
            MyMsDataFile = myMsDataFile;

        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            foreach (PeptideSpectralMatch psm in AllResultingIdentifications.Where(b => b.NumDifferentMatchingPeptides == 1))
            {
                // Stop loop if canceled
                if (GlobalVariables.StopLoops)
                {
                    break;
                }
                
                // don't try to localize mass-differences for ambiguous sequences
                if (psm.FullSequence == null)
                {
                    continue;
                }

                MsDataScan scan = MyMsDataFile.GetOneBasedScan(psm.ScanNumber);
                PeptideWithSetModifications peptide = psm.BestMatchingPeptideWithSetMods.First().Pwsm;
                double massDifference = psm.ScanPrecursorMass - peptide.MonoisotopicMass;

                // this section will iterate through all residues of the peptide and try to localize the mass-diff at each residue and report a score for each residue
                var localizedScores = new List<double>();
                for (int r = 0; r < peptide.Length; r++)
                {
                    // create new PeptideWithSetMods with unidentified mass difference at the given residue
                    PeptideWithSetModifications peptideWithLocalizedMassDiff = peptide.Localize(r, massDifference);

                    // this is the list of theoretical products for this peptide with mass-difference on this residue
                    List<Product> productsWithLocalizedMassDiff = peptideWithLocalizedMassDiff.Fragment(commonParameters.DissociationType, commonParameters.FragmentationTerminus).ToList();

                    var matchedIons = MatchFragmentIons(scan, productsWithLocalizedMassDiff, commonParameters);

                    if (commonParameters.AddCompIons)
                    {
                        MzSpectrum complementarySpectrum = GenerateComplementarySpectrum(scan.MassSpectrum, psm.ScanPrecursorMass, commonParameters.DissociationType);
                        //matchedIons.AddRange(MatchFragmentIons(complementarySpectrum, productsWithLocalizedMassDiff, commonParameters));
                    }

                    // score when the mass-diff is on this residue
                    double localizedScore = CalculatePeptideScore(scan, matchedIons, 0);
                    
                    localizedScores.Add(localizedScore);
                }

                psm.LocalizedScores = localizedScores;
            }

            return new LocalizationEngineResults(this);
        }
    }
}