﻿using Nett;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        private static List<Modification> _AllModsKnown = new List<Modification>();
        private static HashSet<string> _AllModTypesKnown = new HashSet<string>();
        public static Dictionary<string, Modification> AllModsKnownDictionary = new Dictionary<string, Modification>();

        static GlobalVariables()
        {
            MetaMorpheusVersion = typeof(GlobalVariables).Assembly.GetName().Version.ToString();

            if (MetaMorpheusVersion.Equals("1.0.0.0"))
            {
#if DEBUG
                MetaMorpheusVersion = "Not a release version. DEBUG.";
#else
                MetaMorpheusVersion = "Not a release version.";
#endif
            }
            else
            {
                // as of 0.0.277, AppVeyor appends the build number
                // this is intentional; it's to avoid conflicting AppVeyor build numbers
                // trim the build number off the version number for displaying/checking versions, etc
                var foundIndexes = new List<int>();
                for (int i = 0; i < MetaMorpheusVersion.Length; i++)
                {
                    if (MetaMorpheusVersion[i] == '.')
                        foundIndexes.Add(i);
                }
                MetaMorpheusVersion = MetaMorpheusVersion.Substring(0, foundIndexes.Last());
            }

            {
                var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
                if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles) && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
                    DataDir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
                else
                    DataDir = AppDomain.CurrentDomain.BaseDirectory;
            }

            ElementsLocation = Path.Combine(DataDir, @"Data", @"elements.dat");
            UsefulProteomicsDatabases.Loaders.LoadElements(ElementsLocation);

            ExperimentalDesignFileName = "ExperimentalDesign.tsv";

            UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(Path.Combine(DataDir, @"Data", @"unimod.xml")).ToList();
            PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(Path.Combine(DataDir, @"Data", @"PSI-MOD.obo.xml"));
            var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
            UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(Path.Combine(DataDir, @"Data", @"ptmlist.txt"), formalChargesDictionary).ToList();

            foreach (var modFile in Directory.GetFiles(Path.Combine(DataDir, @"Mods")))
                AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(modFile));

            // TODO: need to add motif to Unimod/UniProt ID
            //AddMods(UnimodDeserialized.OfType<Modification>());
            //AddMods(UniprotDeseralized.OfType<Modification>());

            foreach (var mod in AllModsKnown)
            {
                if (!AllModsKnownDictionary.ContainsKey(mod.IdWithMotif))
                {
                    AllModsKnownDictionary.Add(mod.IdWithMotif, mod);
                }
            }

            GlobalSettings = Toml.ReadFile<GlobalSettings>(Path.Combine(DataDir, @"settings.toml"));
        }

        // File locations
        public static string DataDir { get; }
        public static bool StopLoops { get; set; }
        public static string ElementsLocation { get; }
        public static string MetaMorpheusVersion { get; }
        public static IGlobalSettings GlobalSettings { get; }
        public static IEnumerable<Modification> UnimodDeserialized { get; }
        public static IEnumerable<Modification> UniprotDeseralized { get; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }
        public static IEnumerable<Modification> AllModsKnown { get { return _AllModsKnown.AsEnumerable(); } }
        public static IEnumerable<string> AllModTypesKnown { get { return _AllModTypesKnown.AsEnumerable(); } }
        public static string ExperimentalDesignFileName { get; }

        public static void AddMods(IEnumerable<Modification> enumerable)
        {
            foreach (var ye in enumerable)
            {
                if (string.IsNullOrEmpty(ye.ModificationType) || string.IsNullOrEmpty(ye.IdWithMotif))
                    throw new MetaMorpheusException(ye.ToString() + Environment.NewLine + " has null or empty modification type");
                if (AllModsKnown.Any(b => b.IdWithMotif.Equals(ye.IdWithMotif) && b.ModificationType.Equals(ye.ModificationType) && !b.Equals(ye)))
                {
                    throw new MetaMorpheusException("Modification id and type are equal, but some fields are not! Please modify/remove one of the modifications: " + Environment.NewLine + Environment.NewLine + ye.ToString() + Environment.NewLine + Environment.NewLine + " has same and id and modification type as " + Environment.NewLine + Environment.NewLine + AllModsKnown.First(b => b.IdWithMotif.Equals(ye.IdWithMotif) && b.ModificationType.Equals(ye.ModificationType)) + Environment.NewLine + Environment.NewLine);
                }
                else if (AllModsKnown.Any(b => b.IdWithMotif.Equals(ye.IdWithMotif) && b.ModificationType.Equals(ye.ModificationType)))
                    continue;
                else
                {
                    _AllModsKnown.Add(ye);
                    _AllModTypesKnown.Add(ye.ModificationType);
                }
            }
        }

        public static string CheckLengthOfOutput(string psmString)
        {
            if (psmString.Length > 32000 && GlobalSettings.WriteExcelCompatibleTSVs)
            {
                return "Output too long for Excel";
            }
            else
            {
                return psmString;
            }
        }

        public static Dictionary<string, Protease> LoadProteaseDictionary(string proteasesLocation)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            using (StreamReader proteases = new StreamReader(proteasesLocation))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[][] fields = line.Split('\t').Select(x => x.Split('|')).ToArray();
                    string name = fields[0][0];
                    string[] preventing;
                    List<Tuple<string, FragmentationTerminus>> sequences_inducing_cleavage = new List<Tuple<string, FragmentationTerminus>>();
                    List<Tuple<string, FragmentationTerminus>> sequences_preventing_cleavage = new List<Tuple<string, FragmentationTerminus>>();
                    for (int i = 0; i < fields[1].Length; i++)
                    {
                        if (!fields[1][i].Equals(""))
                        {
                            sequences_inducing_cleavage.Add(new Tuple<string, FragmentationTerminus>(fields[1][i], ((FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), fields[3][i], true))));
                            if (!fields[2].Contains(""))
                            {
                                preventing = (fields[2][i].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries));
                                for (int j = 0; j < preventing.Length; j++)
                                {
                                    sequences_preventing_cleavage.Add(new Tuple<string, FragmentationTerminus>(preventing[j], (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), fields[3][i], true)));
                                }
                            }
                        }
                    }
                    var cleavage_specificity = ((CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4][0], true));
                    string psi_ms_accession_number = fields[5][0];
                    string psi_ms_name = fields[6][0];
                    string site_regexp = fields[7][0];
                    var protease = new Protease(name, sequences_inducing_cleavage, sequences_preventing_cleavage, cleavage_specificity, psi_ms_accession_number, psi_ms_name, site_regexp);
                    dict.Add(protease.Name, protease);
                }
            }
            return dict;
        }
    }
}