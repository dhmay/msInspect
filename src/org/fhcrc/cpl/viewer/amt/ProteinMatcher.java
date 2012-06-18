/*
 * Copyright (c) 2003-2012 Fred Hutchinson Cancer Research Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.fhcrc.cpl.viewer.amt;

import java.util.*;
import java.io.*;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.apache.log4j.Logger;

public class ProteinMatcher
{
    private static Logger _log = Logger.getLogger(ProteinMatcher.class);

    //minimum peptide length we want to try to match
    protected static final int MIN_PEPTIDE_LENGTH=5;
    //maximum number of missed cleavages to account for when matching
    protected static final int MAX_MISSED_CLEAVAGES=1;

    /**
     * This is a bit complicated.  The point is to take an MS1 feature set that has already been matched
     * with MS2 (or with an AMT database), and mine it for more matches, against all the proteins that
     * we think might be there.
     *
     * Why might we think proteins are there?  Well, they might be matched in the MS1 file itself.  Or
     * they might be in the embedded MS2 (ms2Features).  Or they might be everything in the AMT database,
     * which presumably represents related runs.
     * @param matchedMs1Features
     * @param ms2Features
     * @param amtDatabase
     * @param allProteins
     * @param minPeptideMass
     * @param maxPeptideMass
     * @param featureSetMatcher
     * @param modifications
     * @return
     */
    public static Map<Protein,List<Feature>> matchAdditionalPeptidesForMatchedProteins(
            Feature[] matchedMs1Features, Feature[] ms2Features, AmtDatabase amtDatabase, Protein[] allProteins,
            float minPeptideMass, float maxPeptideMass,
            FeatureSetMatcher featureSetMatcher,
            MS2Modification[] modifications)
    {
        List<Feature> ms1FeaturesWithoutMatches = new ArrayList<Feature>();
        Set<String> initiallyMatchedPeptides = new HashSet<String>();
        for (Feature matchedMs1Feature : matchedMs1Features)
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(matchedMs1Feature);
            if (peptide == null)
                ms1FeaturesWithoutMatches.add(matchedMs1Feature);
            else
                initiallyMatchedPeptides.add(peptide);
        }
        if (ms2Features != null)
        {
            for (Feature ms2Feature : ms2Features)
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(ms2Feature);
                if (peptide == null)
                    ms1FeaturesWithoutMatches.add(ms2Feature);
                else
                    initiallyMatchedPeptides.add(peptide);
            }
        }


        if (amtDatabase != null)
            for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
                initiallyMatchedPeptides.add(peptideEntry.getPeptideSequence());

        Set<Protein> matchedProteins = findProteinsForPeptides(initiallyMatchedPeptides, allProteins);
        _log.debug("matched proteins: " + matchedProteins.size());

        Map<Protein,List<Feature>> proteinNewMatchMap =
                new HashMap<Protein,List<Feature>>();
        Map<Feature, Protein> peptideFeatureProteinMap =
                new HashMap<Feature, Protein>();
        for (Protein matchedProtein : matchedProteins)
        {
            List<Feature> unmatchedPeptideFeatures =
                    createFeaturesForProteinPeptides(matchedProtein, initiallyMatchedPeptides,
                                                     minPeptideMass, maxPeptideMass, modifications);
            for (Feature feature : unmatchedPeptideFeatures)
            {
                if (!initiallyMatchedPeptides.contains(MS2ExtraInfoDef.getFirstPeptide(feature)))
                    peptideFeatureProteinMap.put(feature, matchedProtein);
            }
//            _log.debug("protein: created peptide features: " + unmatchedPeptideFeatures.size());
        }
        _log.debug("Created " + peptideFeatureProteinMap.size() + " features representing unmatched peptides");
        FeatureSet unmatchedPeptideFeatureSet =
                new FeatureSet(peptideFeatureProteinMap.keySet().toArray(new Feature[peptideFeatureProteinMap.size()]));
        _log.debug("Matching against " + ms1FeaturesWithoutMatches.size() + " unmatched MS1 features");
        FeatureSet initiallyUnmatchedFeatureSet =
                new FeatureSet(ms1FeaturesWithoutMatches.toArray(new Feature[ms1FeaturesWithoutMatches.size()]));

        FeatureSetMatcher.FeatureMatchingResult matchingResult =
                featureSetMatcher.matchFeatures(initiallyUnmatchedFeatureSet, unmatchedPeptideFeatureSet);

        for (Feature feature : matchingResult.getMasterSetFeatures())
        {
            Feature bestPeptideFeatureMatch = matchingResult.getBestMatch(feature);
            Protein sourceProtein = peptideFeatureProteinMap.get(bestPeptideFeatureMatch);
//System.err.println("protein: " + sourceProtein.getName());
            List<Feature> thisProteinMatchedFeatures =
                    proteinNewMatchMap.get(sourceProtein);
            if (thisProteinMatchedFeatures == null)
            {
                thisProteinMatchedFeatures = new ArrayList<Feature>();
                proteinNewMatchMap.put(sourceProtein, thisProteinMatchedFeatures);
            }
            MS2ExtraInfoDef.setSinglePeptide(feature,
                                             MS2ExtraInfoDef.getFirstPeptide(bestPeptideFeatureMatch));

            thisProteinMatchedFeatures.add(feature);
        }

        _log.debug("Matched " + proteinNewMatchMap.size() + " proteins out of " +
                   matchedProteins.size() + " initially matched");
        _log.debug("Peptides matched: " + matchingResult.getMasterSetFeatures().size());
        return proteinNewMatchMap;
    }

    public static Set<Protein> findProteinsForPeptides(Set<String> peptides, Protein[] proteins)
    {
        Set<Protein> result = new HashSet<Protein>();
        for (Protein protein : proteins)
        {
            String proteinSequence = protein.getSequenceAsString();
            for (String peptide : peptides)
            {
                if (proteinSequence.contains(peptide))
                    result.add(protein);
            }
        }
        return result;
    }

    /**
     * TODO: modifications
     * @param protein
     * @param peptidesToExclude
     * @param minPeptideMass
     * @param maxPeptideMass
     * @return
     */
    public static List<Feature> createFeaturesForProteinPeptides(Protein protein,
                                                             Set<String> peptidesToExclude,
                                                             float minPeptideMass,
                                                             float maxPeptideMass,
                                                             MS2Modification[] modifications)
    {
        List<Feature> result = new ArrayList<Feature>();
        List<Peptide> proteinPeptides =
                generatePeptidesFromProtein(protein, 0, MIN_PEPTIDE_LENGTH);
//        _log.debug("createFeatures: peptides: " + proteinPeptides.size());
        for (Peptide peptide : proteinPeptides)
        {
            if (peptide.getMonoisotopicMass() >= minPeptideMass &&
                peptide.getMonoisotopicMass() <= maxPeptideMass &&
                !peptidesToExclude.contains(peptide.getChars()))
            {
                AmtDatabaseFeatureSetGenerator featureGen =
                        new AmtDatabaseFeatureSetGenerator(null);
        List<MS2Modification> varModList = new ArrayList<MS2Modification>();
        List<MS2Modification> staticModList = new ArrayList<MS2Modification>();
        if (modifications != null)
            for (MS2Modification modification : modifications)
            {
                if (modification.getVariable())
                    varModList.add(modification);
                else
                    staticModList.add(modification);
            }
                List<Feature> peptideFeatures = null;
//                        featureGen.generateModFeaturesForPeptide(
//                        new String(peptide.getChars()),
//                        HydrophobicityNormalizer.normalize(peptide.getHydrophobicity3(),"krokhin",3),
//                        staticModList,
//                        varModList);
                for (Feature feature : peptideFeatures)
                    result.add(feature);

//                Feature peptideFeature = new Feature();
//                peptideFeature.setScan(1);
//                peptideFeature.setPeaks(1);
//                AmtExtraInfoDef.setObservedHydrophobicity(peptideFeature,
//                        HydrophobicityNormalizer.normalize(peptide.getHydrophobicity3(),"krokhin",3));
//                MS2ExtraInfoDef.setSinglePeptide(peptideFeature, new String(peptide.getChars()));
//                peptideFeature.setCharge(1);
//                peptideFeature.setMass((float) peptide.getMonoisotopicMass());
//                peptideFeature.updateMz();
////System.err.println("Adding feature with mass " + peptideFeature.getMass() + ", hydro " + AmtExtraInfoDef.getObservedHydrophobicity(peptideFeature));
//
//                result.add(peptideFeature);
            }
        }
        return result;
    }


    /**
     * Returns a HashMap containing, for each Protein, all the features that match peptides found
     * by tryptic digestion.  Each feature in the FeatureSet will be marked with its peptide
     * @param ms1Features
     * @param proteins
     * @param minProteinMatchFeatures minimum number of features to match in order for
     * protein to be included in results
     * @return
     */
    public static HashMap<Protein, Map<String, Feature>> findMatchesForAllProteins(FeatureSet ms1Features,
                                                  Protein[] proteins,
                                                  MS2Modification[] modifications,
                                                  AmtDatabase amtDatabase,
                                                  int minProteinMatchFeatures)
    {
        Feature[] ms1FeatureArray = ms1Features.getFeatures().clone();
        //make sure none of the features are excluded, to start with
        for (Feature feature : ms1FeatureArray)
        {
            feature.excluded = false;
        }

        Feature.MassAscComparator comp = new Feature.MassAscComparator();
        Arrays.sort(ms1FeatureArray, comp);

//        double[] lowMassTable = PeptideMatcher.generateMassTableWithStaticMods(
//                PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES,
//                modifications);
        HashMap<Protein, Map<String, Feature>> result = new HashMap<Protein, Map<String, Feature>>();

        for (int i=0; i<proteins.length; i++)
        {
            Protein currentProtein = proteins[i];
            //status message
            if (i%(proteins.length/10) == 0)
                _log.debug("Processed " + i + " out of " + proteins.length + " proteins");

            Map<String, Feature> foundFeatures = findFeaturesForProtein(ms1FeatureArray, currentProtein,
                                                    modifications, amtDatabase);
            if (foundFeatures.size() >= minProteinMatchFeatures)
                result.put(currentProtein, foundFeatures);
        }

        return result;
    }



    /**
     * Runs digestProtein to break stuff up, setting the digest type, cleavages and peptide
     * length appropriately
     * @param protein
     * @return
     */
    public static ArrayList<Peptide> generatePeptidesFromProtein(Protein protein,
                                                                 int maxMissedCleavages,
                                                                 int minPeptideLength)
    {
        PeptideGenerator peptideGenerator = new PeptideGenerator();
        //use tryptic digest
        peptideGenerator.setDigest(PeptideGenerator.DIGEST_TRYPTIC);
        //allow one missed cleavage
        peptideGenerator.setMaxMissedCleavages(maxMissedCleavages);
        //set minimum residues
        peptideGenerator.setMinResidues(minPeptideLength);

        Peptide[] peptidesFromProtein = peptideGenerator.digestProtein(protein);

        ArrayList<Peptide> peptidesFromProteinList = new ArrayList<Peptide>();

        for (Peptide peptide : peptidesFromProtein)
            peptidesFromProteinList.add(peptide);

        return peptidesFromProteinList;
    }

    /**
     * Runs digestProtein to break stuff up, setting the digest type, cleavages and peptide
     * length appropriately
     * @param protein
     * @return
     */
    public static ArrayList<Peptide> generatePeptidesFromProtein(Protein protein)
    {
        return generatePeptidesFromProtein(protein, MAX_MISSED_CLEAVAGES, MIN_PEPTIDE_LENGTH);
    }


    /**
     * returns an ArrayList of features, the matches between the ms1Features and the peptides
     * in the protein that are also in the amt database.  First copies the database and pares out
     * the entries that don't exist as peptides in this protein, then matches the remaining database
     * against the MS1 features
     * @param ms1Features
     * @param protein
     * @return
     */
    public static Map<String,Feature> findFeaturesForProtein(Feature[] ms1Features, Protein protein,
                                                  MS2Modification[] modifications,
                                                  AmtDatabase amtDatabase)
    {
        ArrayList<Peptide> peptidesFromProtein = generatePeptidesFromProtein(protein);
//System.err.println("protein peptides: " + peptidesFromProtein.size());
        Set<String> peptideSequencesInAmtFromProtein = new HashSet<String>();
        for (Peptide peptide : peptidesFromProtein)
        {
            String peptideSequence = new String(peptide.getChars());
            if (amtDatabase.contains(peptideSequence))
            {
                peptideSequencesInAmtFromProtein.add(peptideSequence);
            }
//for (String amtPeptide : amtDatabase.getPeptides())
//{
//    if (peptideSequence.equals(amtPeptide)) System.err.println("A MATCH!!!! entry?" + amtDatabase.getEntry(peptideSequence) + ",***contains: " + amtDatabase.contains(peptideSequence));
//}
        }

        //if the protein and the db don't even overlap in peptide IDs, give up
        if (peptideSequencesInAmtFromProtein.isEmpty())
        {
            _log.debug("No AMT database entries with peptides contained in this protein");
            return new HashMap<String, Feature>();
        }

//        AmtDatabase amtDatabaseCopy = (AmtDatabase) amtDatabase.waistDeepCopy();
//        for (AmtPeptideEntry peptideEntry : amtDatabaseCopy.getEntries())
//        {
//            if (!peptideSequencesInAmtFromProtein.contains(peptideEntry.getPeptideSequence()))
//            {
//                amtDatabaseCopy.removeEntry(peptideEntry.getPeptideSequence());
//            }
//        }
        AmtDatabase amtDatabaseCopy = new AmtDatabase();
        for (AmtRunEntry runEntry : amtDatabase.getRuns())
            amtDatabaseCopy.addRunEntry(runEntry);
        for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
        {
            if (peptideSequencesInAmtFromProtein.contains(peptideEntry.getPeptideSequence()))
                amtDatabaseCopy.addOrOverrideEntry(peptideEntry);
        }

        //TODO: some way to do something other than the defaults

        //TODO: use alignment
        AmtDatabaseMatcher amtDatabaseMatcher = new AmtDatabaseMatcher();
        Map<Feature, List<AmtPeptideEntry>> featureEntryListMap = null;
//                amtDatabaseMatcher.matchAmtDatabaseWithMS1Features(amtDatabaseCopy, ms1Features,
//                                                                modifications, null, 1, false, false);
        //TODO: this is so suspect, it's not even funny
        Map<String, Feature> result =
                new HashMap<String, Feature>();
        for (Feature feature : featureEntryListMap.keySet())
        {
            List<AmtPeptideEntry> entryList = featureEntryListMap.get(feature);
            result.put(entryList.get(0).getPeptideSequence(), feature);
        }

        return result;
    }


    //methods for comparison of protein/peptide files


    /**
     * Just handles file IO
     * @param file
     * @return
     */
    public static BufferedReader openProteinPeptideFile(File file)
    {
        BufferedReader reader = null;
        try
        {
            reader = new BufferedReader(new FileReader(file));
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        return reader;
    }

    /**
     * skips past any peptide lines to the next protein, in my funky file format
     * @param reader
     * @return
     */
    public static Protein getNextProtein(BufferedReader reader)
    {
            Protein result = null;
        try{
        String line;
            while ((line = reader.readLine()) != null)
            {
                if (line.length() > 0 && line.charAt(0) == '>')
                {
                    String proteinHeader = line.substring(1);

                    //process next line as protein sequence
                    ByteArrayOutputStream aaStream = new ByteArrayOutputStream(2048);
                    line = reader.readLine();
                    byte[] bytes = line.getBytes();
                    for (int i = 0; i < bytes.length; i++)
                        if ((bytes[i] >= 'A') && (bytes[i] <= 'Z'))
                        {
                            //_aaCounts[bytes[i] - 'A'] ++;
                            aaStream.write(bytes[i]);
                        }

                    result = new Protein(proteinHeader, aaStream.toByteArray());
//System.err.println("Writing protein " + currentProtein.getSequenceAsString());
                    break;
                }
            }
        }
            catch (Exception e) { e.printStackTrace(System.err); }
        return result;
    }

    /**
     * In my funky file format, gobbles a list of peptides, stopping when it hits eof or
     * a protein definition
     * @param reader
     * @return
     */
    public static ArrayList<String> getPeptides(BufferedReader reader)
    {
        ArrayList<String> result = new ArrayList<String>();
        try
        {
            String line;
            reader.mark(2000);
            while ((line = reader.readLine()) != null)
            {
                if (line.length() > 0 && line.charAt(0) == '>')
                {
                    //set the mark back to the last protein
                    reader.reset();
                    break;
                }
                else
                {
                        result.add(line);
                    reader.mark(2000);
                }
            }
        }
                catch (Exception e)
        {
//TODO: better handling
            e.printStackTrace(System.err);
        }
        return result;
    }

    public static int compareProteinsBySequence(Protein prot1, Protein prot2)
    {
        return new Protein.SequenceComparator().compare(prot1, prot2);
    }


}
