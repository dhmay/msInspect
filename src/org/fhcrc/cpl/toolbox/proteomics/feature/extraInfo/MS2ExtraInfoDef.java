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
package org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.BrowserController;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;

import javax.swing.*;
import java.util.*;
import java.lang.reflect.InvocationTargetException;
import java.awt.event.ActionEvent;

/**
 * FeatureExtraInformationDef related to MS2.  Defines a class that handles the conversion
 * of peptide and protein lists to and from their file-storage format, and allows you
 * to work with them in memory.
 *
 * Ideally, rather than two separate lists, we'd store a single list of a combined
 * peptide-protein association.  Unfortunately, that's quite difficult, as we need
 * for the columns in the .tsv to correspond to single objects.  If we ever change
 * the format to be a combined peptide-protein column, we can revisit this.
 *
 * Some effort is made to keep peptide and protein lists in synch.  That is, if you
 * add a peptide without an associated protein, a null entry will be put on the
 * protein list.  In a situation where NO features have proteins, this is relatively
 * harmless, and in a situation where some features have proteins, this allows you
 * to maintain the association between peptides and proteins -- the peptides and
 * proteins will be in corresponding array indices.
 */
public class MS2ExtraInfoDef extends FeatureExtraInformationDef
{
    static Logger _log = Logger.getLogger(MS2ExtraInfoDef.class);

    public static final int NO_ENZYMATIC_ENDS_SPECIFIED = -1;


    public static final float defaultPeptideProphet = 0;
    public static final float defaultDeltaMass = 0;
    public static final float defaultFDR = -1;

    public static final float defaultFval = 0;
    public static final int defaultNumEnzymaticEnds = NO_ENZYMATIC_ENDS_SPECIFIED;




    /**
     * initialize the name of the extra information, the column list, and the datatypes
     * of the columns
     */
    public MS2ExtraInfoDef()
    {
        super(
                    "MS2",
                    new String[]{
                            "peptide",
                            "protein",
                            "peptideprophet",
                            "deltamass",
                            "modifiedaminoacids",
                            "search_scores",
                            "fval",
                            "num_enzymatic_ends",
                            "prev_aa",
                            "next_aa",
                            "all_ntt_prob",
                            "alt_protein_ntts",
                            "nterm_modmass",
                            "cterm_modmass",
                            "fdr"
                    },
                    new Class[]{
                            List.class, List.class,
                            Double.class, Double.class, Map.class, Map.class, Double.class,
                            Integer.class, String.class, String.class, String.class, List.class,
                            Double.class, Double.class, Double.class
                    },
                    new String[]{
                            "modifications",
                            "search_database_path",
                            "search_constraint_max_int_cleavages",
                            "search_constraint_min_termini",
                            "base_name"
                    }

            );
    }

    protected static MS2ExtraInfoDef singletonInstance = null;

    public static MS2ExtraInfoDef getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new MS2ExtraInfoDef();
        return singletonInstance;
    }

    /**
     * Does this feature have PeptideProphet information associated with it?
     * @param feature
     * @return
     */
    public static boolean hasPeptideProphet(Feature feature)
    {
        return feature.hasProperty("peptideprophet");
    }


    public static boolean hasAllNttProb(Feature feature)
    {
        return feature.hasProperty("all_ntt_prob");
    }

    //getters and setters for the various columns.
    //These are here so that a caller doesn't have to type the column name,
    //which wouldn't get validated by the compiler.

    public static double getPeptideProphet(Feature feature)
    {
        return feature.getDoubleProperty("peptideprophet", defaultPeptideProphet);
    }

    public static void setPeptideProphet(Feature feature, double peptideprophet)
    {
        feature.setProperty("peptideprophet", peptideprophet);
    }

    public static double getFval(Feature feature)
    {
        return feature.getDoubleProperty("fval", defaultFval);
    }

    public static void setFval(Feature feature, double fval)
    {
        feature.setProperty("fval", fval);
    }

    public static float getDeltaMass(Feature feature)
    {
        try
        {
            return (float) feature.getDoubleProperty("deltamass", defaultDeltaMass);
        }
        catch (ClassCastException cce)
        {
            return feature.getFloatProperty("deltamass", defaultDeltaMass);
        }
    }

    public static void setDeltaMass(Feature feature, float deltaMass)
    {
        feature.setProperty("deltamass", deltaMass);
    }

    public static float getFalseDiscoveryRate(Feature feature)
    {
        try
        {
            return (float) feature.getDoubleProperty("fdr", defaultFDR);
        }
        catch (ClassCastException cce)
        {
            return feature.getFloatProperty("fdr", defaultFDR);
        }
    }

    public static boolean hasFalseDiscoveryRate(Feature feature)
    {
        return getFalseDiscoveryRate(feature) != defaultFDR;
    }

    public static void setFalseDiscoveryRate(Feature feature, float fdr)
    {
        feature.setProperty("fdr", fdr);
    }

    public static int getNumEnzymaticEnds(Feature feature)
    {
        return feature.getIntProperty("num_enzymatic_ends", defaultNumEnzymaticEnds);
    }

    public static void setNumEnzymaticEnds(Feature feature, int numEnzymaticEnds)
    {
        feature.setProperty("num_enzymatic_ends", numEnzymaticEnds);
    }

    public static boolean hasNumEnzymaticEnds(Feature feature)
    {
        return feature.hasProperty("num_enzymatic_ends");
    }

    public static Character getPrevAminoAcid(Feature feature)
    {
        String resultAsString = feature.getStringProperty("prev_aa", null);
        if (resultAsString == null)
            return null;
        return resultAsString.charAt(0);
    }

    public static void setPrevAminoAcid(Feature feature, char prevAminoAcid)
    {
        feature.setProperty("prev_aa", "" + prevAminoAcid);
    }

    public static Character getNextAminoAcid(Feature feature)
    {
        String resultAsString = feature.getStringProperty("next_aa", null);
        if (resultAsString == null)
            return null;
        return resultAsString.charAt(0);
    }

    public static void setNextAminoAcid(Feature feature, char nextAminoAcid)
    {
        feature.setProperty("next_aa", "" + nextAminoAcid);
    }

    public static String getAllNttProb(Feature feature)
    {
        return feature.getStringProperty("all_ntt_prob", null);
    }

    public static void setAllNttProb(Feature feature, String allNttProb)
    {
        feature.setProperty("all_ntt_prob", allNttProb);
    }

    public static List<Integer> getAltProteinNTTs(Feature feature)
    {
        try
        {
            return (List<Integer>) feature.getProperty("alt_protein_ntts", null);
        }
        catch (ClassCastException e)
        {
            ApplicationContext.infoMessage("getAltProteinNTTs, ClassCastException processing value " +
                    feature.getProperty("alt_protein_ntts", null) + ", class " +
                    feature.getProperty("alt_protein_ntts", null).getClass().getName());
            System.err.println(feature);                        
            throw e;
        }
    }

    public static void setAltProteinNTTs(Feature feature, List<Integer> altProteinNTTs)
    {
        feature.setProperty("alt_protein_ntts", altProteinNTTs);
    }

    public static float getNtermModMass(Feature feature)
    {
        try
        {
            return (float) feature.getDoubleProperty("nterm_mod_mass", 0);
        }
        catch (ClassCastException cce)
        {
            return feature.getFloatProperty("nterm_mod_mass", 0);
        }
    }

    public static void setNtermModMass(Feature feature, float mass)
    {
        feature.setProperty("nterm_mod_mass", mass);
    }

    public static float getCtermModMass(Feature feature)
    {
        try
        {
            return (float) feature.getDoubleProperty("cterm_mod_mass", 0);
        }
        catch (ClassCastException cce)
        {
            return feature.getFloatProperty("cterm_mod_mass", 0);
        }
    }

    public static void setCtermModMass(Feature feature, float mass)
    {
        feature.setProperty("cterm_mod_mass", mass);
    }





    public static Feature createMS2Feature(int scan, float mass, int charge,
                                           String peptideListString, String proteinListString,
                                           List<ModifiedAminoAcid>[] modifiedAminoAcids)

    {
        List<String> peptideList = parseStringListString(peptideListString);
        List<String> proteinList = parseStringListString(proteinListString);

        return createMS2Feature(scan, mass, charge, peptideList, proteinList, modifiedAminoAcids);
    }


    public static Feature createMS2Feature(int scan, float mass, int charge,
                                           List<String> peptideList, List<String> proteinList,
                                           List<ModifiedAminoAcid>[] modifiedAminoAcids)

    {
        Feature result = new Feature();
        result.setScan(scan);
        result.setScanFirst(scan);
        result.setScanLast(scan);
        result.setMass(mass);
        result.setCharge(charge);
        result.setScanCount(1);
        result.updateMz();
        setPeptideList(result, peptideList);
        setProteinList(result, proteinList);
        setModifiedAminoAcids(result,modifiedAminoAcids);
        return result;
    }

    public Object convertStringValue(String columnName, String value)
            throws InstantiationException, IllegalAccessException,
            InvocationTargetException, NoSuchMethodException
    {
        if ("modifiedaminoacids".equalsIgnoreCase(columnName))
            return parsePositionModifiedAminoAcidListMapString(value);
        else if ("peptide".equalsIgnoreCase(columnName) ||
                 "protein".equalsIgnoreCase(columnName))
            return parseStringListString(value);
        else if ("search_scores".equalsIgnoreCase(columnName))
            return parseStringDoubleMapString(value);
        else if ("alt_protein_ntts".equalsIgnoreCase(columnName))
            return parseIntListString(value);
        else
            return super.convertStringValue(columnName, value);
    }

    public static String convertModifiedAminoAcidsMapToString(Map<Integer, List<ModifiedAminoAcid>> modifiedAminoAcids)
    {
        if (modifiedAminoAcids == null)
            return "";
        List<String> positionsList = new ArrayList<String>();
        for (int position :  modifiedAminoAcids.keySet())
        {
            if (modifiedAminoAcids.get(position) == null)
                continue;
            List<String> modStringList = new ArrayList<String>();
            for (ModifiedAminoAcid mod : modifiedAminoAcids.get(position))
                modStringList.add(mod.toString());
            //convert zero-based to one-based
            positionsList.add((position + 1) + "(" +
                    convertStringListToString(modStringList, ",") + ")");
        }
        return convertStringListToString(positionsList);
    }

    /**
     * Create a modified sequence string, e.g., MC[139.02]EMK
     * @param peptideSequence
     * @param modifiedAminoAcids
     * @return
     */
    public static String createModifiedSequenceString(String peptideSequence,
                                                      List<ModifiedAminoAcid>[] modifiedAminoAcids)
    {
        StringBuffer resultBuf = new StringBuffer();
        for (int i=0; i<peptideSequence.length(); i++)
        {
            resultBuf.append(peptideSequence.charAt(i));
            if (modifiedAminoAcids != null && modifiedAminoAcids[i] != null)
            {
                double massDiff = 0;
                for (ModifiedAminoAcid mod : modifiedAminoAcids[i])
                {
                    massDiff += mod.getMass() - PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)];
                }
//if (peptideSequence.charAt(i) == 'C' && massDiff > 0 && massDiff < 57)
//{
//    System.err.println(i + ": " + MS2ExtraInfoDef.convertModifiedAminoAcidsMapToString(modifiedAminoAcids));
//                for (ModifiedAminoAcid mod : modifiedAminoAcids.get(i))
//                {
//                    System.err.println(mod.getMass() + "-" + PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)]);
//                }
//}
                if (massDiff != 0)
                    resultBuf.append("[" + Rounder.round(PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)] + massDiff, 2) + "]");
            }
        }
        return resultBuf.toString();
    }

    /**
     * Create a modified sequence string, e.g., MC[139.02]EMK
     * @param peptideSequence
     * @param modifiedAminoAcids
     * @return
     */
    public static String createModifiedSequenceString(String peptideSequence, 
                                                      Map<Integer, List<ModifiedAminoAcid>> modifiedAminoAcids)
    {
        StringBuffer resultBuf = new StringBuffer();
        for (int i=0; i<peptideSequence.length(); i++)
        {
            resultBuf.append(peptideSequence.charAt(i));
            if (modifiedAminoAcids != null && modifiedAminoAcids.containsKey(i) && modifiedAminoAcids.get(i) != null)
            {
                double massDiff = 0;
                for (ModifiedAminoAcid mod : modifiedAminoAcids.get(i))
                {
                    massDiff += mod.getMass() - PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)];
                }
//if (peptideSequence.charAt(i) == 'C' && massDiff > 0 && massDiff < 57)
//{
//    System.err.println(i + ": " + MS2ExtraInfoDef.convertModifiedAminoAcidsMapToString(modifiedAminoAcids));
//                for (ModifiedAminoAcid mod : modifiedAminoAcids.get(i))
//                {
//                    System.err.println(mod.getMass() + "-" + PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)]);
//                }
//}
                if (massDiff != 0)
                    resultBuf.append("[" + Rounder.round(PeptideGenerator.getMasses(true)[peptideSequence.charAt(i)] + massDiff, 2) + "]");
            }
        }
        return resultBuf.toString();
    }

    /**
     * Special handling for peptide, protein columns
     *
     * TODO: special handling for modifiedaminoacids
     * @param columnName
     * @param value
     * @return
     */
    public String convertToString(String columnName, Object value)
    {
//System.err.println(columnName);   
        if ("modifiedaminoacids".equalsIgnoreCase(columnName))
        {
            Map<Integer, List<ModifiedAminoAcid>> positionModListMap =
                 (Map<Integer, List<ModifiedAminoAcid>>) value;
            return convertModifiedAminoAcidsMapToString(positionModListMap);
        }
        else if ("peptide".equalsIgnoreCase(columnName) || "protein".equalsIgnoreCase(columnName))
            return convertStringListToString((List<String>) value);
        else if ("search_scores".equalsIgnoreCase(columnName))
            return convertStringDoubleMapToString((Map<String,String>) value);
        else if ("alt_protein_ntts".equalsIgnoreCase(columnName))
        {
            if (value != null && !(value instanceof List))
            {
                new RuntimeException("artifical exception for stack trace").printStackTrace(System.err);
                ApplicationContext.infoMessage("WARNING!!! Non-List value found for alt_protein_ntts.  " +
                        "Setting to null.  Class: " + value.getClass().getName() + ", value: " + value);
                value = null;
            }
            return convertIntListToString((List<Integer>) value);
        }
        else return super.convertToString(columnName, value);
    }

    /**
     * Same as convertToString, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public String convertFeatureSetPropertyToString(String propertyName, Object value)
    {
        if (propertyName.equals("modifications"))
        {
            MS2Modification[] modifications = (MS2Modification[]) value;
            if (modifications == null || modifications.length == 0)
                return "";
            List<String> modificationsAsStrings = new ArrayList<String>(modifications.length);

            for (MS2Modification modification : modifications)
            {
                modificationsAsStrings.add(modification.toString());
            }

            return convertStringListToString(modificationsAsStrings);
        }
        else if (propertyName.equals("search_database_path") ||
                 propertyName.equals("search_constraint_max_int_cleavages") ||
                 propertyName.equals("search_constraint_min_termini") ||
                 propertyName.equals("base_name"))

        {
            return value.toString();
        }
        else throw new IllegalArgumentException(
                "MS2ExtraInfoDef doesn't know about a feature set property named " +
                propertyName);
    }

    /**
     * If either mass or massDiff is missing, set it correctly
     * @param ms2Mod
     */
    public static void updateMS2ModMassOrDiff(MS2Modification ms2Mod)
    {
        if (ms2Mod.getSymbol() != null)
        {
            float symbolMass =
                    (float) PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[ms2Mod.getAminoAcid().charAt(0)];


            if (ms2Mod.getMass() == 0f && ms2Mod.getMassDiff() != 0f)
            {
                ms2Mod.setMass(symbolMass + ms2Mod.getMassDiff());
            }
            else if (ms2Mod.getMass() != 0f && ms2Mod.getMassDiff() == 0f)
            {
                ms2Mod.setMassDiff(ms2Mod.getMass() - symbolMass);
            }
        }
    }

    /**
     * Once all modifications are loaded, we can correct the "mass" values for variable
     * modifications for residues that also have static modifications
     * @param modifications
     */
    public static void correctMS2ModMasses(MS2Modification[] modifications)
    {
        Map<String, List<MS2Modification>> staticModMap = new HashMap<String, List<MS2Modification>>();
        Map<String, List<MS2Modification>> varModMap = new HashMap<String, List<MS2Modification>>();

        for (MS2Modification mod : modifications)
        {
            String residue = mod.getAminoAcid();
            if (mod.getVariable())
            {
                List<MS2Modification> varMods = varModMap.get(residue);
                if (varMods == null)
                {
                    varMods = new ArrayList<MS2Modification>();
                    varModMap.put(residue, varMods);
                }
                varMods.add(mod);
            }
            else
            {
                List<MS2Modification> staticMods = staticModMap.get(residue);
                if (staticMods == null)
                {
                    staticMods = new ArrayList<MS2Modification>();
                    staticModMap.put(residue,staticMods);
                }
                staticMods.add(mod);
            }
        }

        for (String residue : varModMap.keySet())
        {
            List<MS2Modification> staticModsThisResidue = staticModMap.get(residue);
            if (staticModsThisResidue != null && staticModsThisResidue.size() == 1)
            {
                for (MS2Modification varMod : varModMap.get(residue))
                {
                    varMod.setMass(staticModsThisResidue.get(0).getMass() + varMod.getMassDiff());
                }
            }
        }
    }

    /**
     * Save as convertStringValue, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public Object convertFeatureSetPropertyStringValue(String propertyName, String value)
    {
        if (propertyName.equals("modifications"))
        {
            _log.debug("Parsing MS2 modifications");
            List<String> modificationStrings =
                    parseStringListString(value);
            List<MS2Modification> modificationsList =
                    new ArrayList<MS2Modification>(modificationStrings.size());
            for (String modificationString : modificationStrings)
            {
                if (modificationString != null && modificationString.length() > 0)
                {
                    MS2Modification ms2Mod = MS2Modification.parseString(modificationString);
                    updateMS2ModMassOrDiff(ms2Mod);

                    _log.debug("\tadding modification " + ms2Mod);
                    modificationsList.add(ms2Mod);
                }
            }

            correctMS2ModMasses(modificationsList.toArray(new MS2Modification[modificationsList.size()]));

            MS2Modification[] result =
                    modificationsList.toArray(new MS2Modification[modificationsList.size()]);
            return result;
        }
        else if (propertyName.equals("search_database_path") ||
                 propertyName.equals("base_name"))
        {
            return value;
        }
        else if (propertyName.equals("search_constraint_max_int_cleavages") ||
                 propertyName.equals("search_constraint_min_termini"))
        {
            return Integer.parseInt(value);
        }
        else throw new IllegalArgumentException(
                "MS2ExtraInfoDef doesn't know about a feature set property named " +
                        propertyName);
    }



    public static MS2Modification[] getFeatureSetModifications(FeatureSet featureSet)
    {
        return (MS2Modification[])
                getSingletonInstance().getFeatureSetProperty(featureSet, "modifications");
    }

    public static void setFeatureSetModifications(FeatureSet featureSet,
                                                  MS2Modification[] modifications)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "modifications", modifications);
    }

    public static String getFeatureSetSearchDatabasePath(FeatureSet featureSet)
    {
        return (String)
                getSingletonInstance().getFeatureSetProperty(featureSet, "search_database_path");
    }

    public static void setFeatureSetSearchDatabasePath(FeatureSet featureSet,
                                                       String databasePath)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "search_database_path", databasePath);
    }

    public static int getFeatureSetSearchConstraintMaxIntCleavages(FeatureSet featureSet)
    {
        Integer result =
                (Integer) getSingletonInstance().getFeatureSetProperty(featureSet,
                        "search_constraint_max_int_cleavages");
        if (result == null)
            result = 0;
        return result;

    }

    public static void setFeatureSetSearchConstraintMaxIntCleavages(FeatureSet featureSet,
                                                       int maxCleavages)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "search_constraint_max_int_cleavages", maxCleavages);
    }

    public static int getFeatureSetSearchConstraintMinTermini(FeatureSet featureSet)
    {
        Integer result = (Integer)
                getSingletonInstance().getFeatureSetProperty(featureSet, "search_constraint_min_termini");
        if (result == null)
            result = 0;
        return result;
    }

    public static void setFeatureSetSearchConstraintMinTermini(FeatureSet featureSet,
                                                       int maxTermini)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "search_constraint_min_termini", maxTermini);
    }

    public static String getFeatureSetBaseName(FeatureSet featureSet)
    {
        return (String) getSingletonInstance().getFeatureSetProperty(featureSet, "base_name");
    }

    public static void setFeatureSetBaseName(FeatureSet featureSet,
                                                       String baseName)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "base_name", baseName);
    }





    /**
     * This is simple, rather than safe.  If the score name has "=" or "," in it,
     * it will fail, badly
     * @param mapToConvert
     * @return
     */
    protected String convertStringDoubleMapToString(Map<String,String> mapToConvert)
    {
        StringBuffer resultBuf = new StringBuffer();
        if (mapToConvert == null || mapToConvert.size() == 0)
            return "";
        for (String key : mapToConvert.keySet())
        {
            resultBuf.append(MULTI_VALUE_LIST_SEPARATOR + key + "=" +
                             mapToConvert.get(key));
        }
        String result = resultBuf.toString();
        if (result.startsWith(MULTI_VALUE_LIST_SEPARATOR))
            result = result.substring(MULTI_VALUE_LIST_SEPARATOR.length());
        return result;
    }

    public static void setSearchScores(Feature feature,
                                       Map<String,String> searchScoreMap)
    {
        feature.setProperty("search_scores", searchScoreMap);
    }

    /**
     * Convenience method, since most scores are doubles
     * @param feature
     * @param scoreName
     * @param scoreValue
     */
    public static void addSearchScore(Feature feature, String scoreName, double scoreValue)
    {
        addSearchScore(feature, scoreName, "" + scoreValue);
    }

    /**
     * Add a new search score.  Makes no effort to prevent overriding existing scores
     * @param feature
     * @param scoreName
     * @param scoreValue
     */
    public static void addSearchScore(Feature feature, String scoreName, String scoreValue)
    {
        Map<String,String> searchScoreMap = getSearchScores(feature);
        boolean previouslyNull = false;
        if (searchScoreMap == null)
        {
            previouslyNull = true;
            searchScoreMap = new HashMap<String,String>(1);
        }
        searchScoreMap.put(scoreName, scoreValue);

        if (previouslyNull)
            setSearchScores(feature, searchScoreMap);
    }

    public static Map<String,String> getSearchScores(Feature feature)
    {
        return (Map<String,String>) feature.getProperty("search_scores");
    }

    /**
     * Find an individual search score for this feature.  Null if doesn't exist
     * @param feature
     * @param scoreName
     * @return
     */
    public static String getSearchScore(Feature feature, String scoreName)
    {
        Map<String,String> searchScores = getSearchScores(feature);
        if (searchScores == null)
            return null;
        return searchScores.get(scoreName);
    }

    /**
     * Brittle.  If the score doesn't exist for this feature, will throw a RTE
     * @param feature
     * @param scoreName
     * @return
     */
    public static double getDoubleSearchScore(Feature feature, String scoreName)
    {
        return Double.parseDouble(getSearchScore(feature, scoreName));
    }


    public static void setModifiedAminoAcids(Feature feature,
                                             List<ModifiedAminoAcid>[] modifiedAminoAcids)
    {
        Map<Integer,List<ModifiedAminoAcid>> positionModListMap = null;
        if (modifiedAminoAcids != null)
        {
            positionModListMap = new HashMap<Integer,List<ModifiedAminoAcid>>();
            for (int i=0; i<modifiedAminoAcids.length; i++)
                positionModListMap.put(i, modifiedAminoAcids[i]);
        }

        feature.setProperty("modifiedaminoacids", positionModListMap);
    }

    public static void setModifiedAminoAcids(Feature feature, Map<Integer, List<ModifiedAminoAcid>> positionModListMap)
    {
        feature.setProperty("modifiedaminoacids", positionModListMap);
    }

    public static List<ModifiedAminoAcid>[] getModifiedAminoAcids(Feature feature)
    {
        Map<Integer, List<ModifiedAminoAcid>> positionModListMap =
            getModifiedAminoAcidsMap(feature);
        if (positionModListMap == null)
            return null;
        List<ModifiedAminoAcid>[] result =
                (List<ModifiedAminoAcid>[]) new List[getFirstPeptide(feature).length()];
        for (int i=0; i<result.length; i++)
            result[i] = positionModListMap.get(i);
        return result;
    }

    public static Map<Integer, List<ModifiedAminoAcid>> getModifiedAminoAcidsMap(Feature feature)
    {
        return (Map<Integer, List<ModifiedAminoAcid>>)
                feature.getProperty("modifiedaminoacids");
    }

    /**
     * Returns the first peptide in the peptide list.  If the list is
     * empty, returns null.  If the list has one item, that's what you get.
     *
     * The first item will generally be the most likely peptide
     * @param feature
     * @return
     */
    public static String getFirstPeptide(Feature feature)
    {
        List<String> peptideList = getPeptideList(feature);
        if (peptideList == null)
            return null;
        return peptideList.get(0);
    }

    public static List<String> getPeptideList(Feature feature)
    {
        return (List<String>) feature.getProperty("peptide");
    }

    public static void setPeptideList(Feature feature, String peptideListString)
    {
        setPeptideList(feature, parseStringListString(peptideListString));
    }

    public static void setPeptideList(Feature feature, String[] peptideArray)
    {
        List<String> peptideList = new ArrayList<String>();
        for (String peptide : peptideArray)
            peptideList.add(peptide);
        setPeptideList(feature, peptideList);
    }

    public static void setPeptideList(Feature feature,
                                      List<String> peptideOrProteinList)
    {
        feature.setProperty("peptide", peptideOrProteinList);
    }

    public static void removeAllPeptides(Feature feature)
    {
        feature.setProperty("peptide", null);
    }

    /**
     * Adds a peptide.  Also adds an associated null protein, to maintain
     * correspondence between positions in the two lists
     * @param feature
     * @param peptide
     */
    public static void addPeptide(Feature feature, String peptide)
    {
        addPeptideWithProtein(feature, peptide, null);
    }

    /**
     *
     * @param feature
     * @param peptide
     */
    public static void setSinglePeptide(Feature feature, String peptide)
    {
        if (peptide != null &&
            peptide.contains(MULTI_VALUE_LIST_SEPARATOR))
            throw new RuntimeException("MS2ExtraInfoDef.setSinglePeptide: Attempted to set a list of peptides as a single peptide");
        setPeptideList(feature, peptide);
    }

    /**
     * Adds a peptide with an associated Protein
     * @param feature
     * @param peptide
     * @param protein
     */
    public static void addPeptideWithProtein(Feature feature, String peptide,
                                            String protein)
    {
        List<String> peptideList = getPeptideList(feature);
        if (peptideList == null)
        {
            peptideList = new ArrayList<String>();
            setPeptideList(feature, peptideList);
        }
        peptideList.add(peptide);

        addProtein(feature, protein);
//System.err.println("  last pep: " + getPeptideList(feature).get(getPeptideList(feature).size()-1) +
//                   ", last prot: " + getProteinList(feature).get(getProteinList(feature).size()-1));
    }




    /**
     * Returns the first peptide in the peptide list.  If the list is
     * empty, returns null.  If the list has one item, that's what you get.
     *
     * The first item will generally be the most likely peptide
     * @param feature
     * @return
     */
    public static String getFirstProtein(Feature feature)
    {
        List<String> proteinList = getProteinList(feature);
        if (proteinList == null)
            return null;
        return proteinList.get(0);
    }

    public static List<String> getProteinList(Feature feature)
    {
        return (List<String>) feature.getProperty("protein");
    }

    public static void setProteinList(Feature feature, String proteinListString)
    {
        //if we've got a non-null list with something in it, go ahead and
        //set the property normally
        if (proteinListString != null && proteinListString.length() > 0)
        {
            setProteinList(feature, parseStringListString(proteinListString));
        }
        else
        {
            //we've got nothing.  If the Feature already has a list of peptides
            //defined, give it a list of null proteins of the same length.
            List<String> featurePeptideList =
                    getPeptideList(feature);
            if (featurePeptideList != null)
            {
                for (int i=0; i<featurePeptideList.size(); i++)
                    addProtein(feature, null);
            }
        }

    }

    public static void setProteinList(Feature feature,
                                      List<String> peptideOrProteinList)
    {
        feature.setProperty("protein", peptideOrProteinList);
    }

    public static void addProtein(Feature feature, String protein)
    {
        List<String> proteinList = getProteinList(feature);
        if (proteinList == null)
        {
            proteinList = new ArrayList<String>();
            setProteinList(feature, proteinList);
        }
        proteinList.add(protein);
    }

    /**
     * parse a string into a map from sequence position to list of modifications.  Never return null --
     * empty map for empty string
     * @param stringListString
     * @return
     * @throws IllegalArgumentException
     */
    public static Map<Integer, List<ModifiedAminoAcid>>
         parsePositionModifiedAminoAcidListMapString(String stringListString)
            throws IllegalArgumentException
    {
        if (stringListString == null || stringListString.length() == 0)
            return new HashMap<Integer, List<ModifiedAminoAcid>>(0);

        List<String> positionList = parseStringListString(stringListString);
        Map<Integer, List<ModifiedAminoAcid>> result =
                new HashMap<Integer, List<ModifiedAminoAcid>>(positionList.size());
        for (String positionEntryString : positionList)
        {
            
            if (!positionEntryString.contains("(") ||
                      positionEntryString.charAt(positionEntryString.length()-1) != ')')
                          throw new IllegalArgumentException("Badly formatted modified amino acid column: " + stringListString);

            int positionOfOpenParen = positionEntryString.indexOf('(');
            int position = Integer.parseInt(positionEntryString.substring(0,positionOfOpenParen));
            //convert one-based to zero-based
            position--;
             String modListString =
                    positionEntryString.substring(positionOfOpenParen+1,positionEntryString.length()-2);

            List<String> modStringList = parseStringListString(modListString, ",");
            List<ModifiedAminoAcid> modsThisPosition =
                    new ArrayList<ModifiedAminoAcid>(modStringList.size());
            for (String modString : modStringList)
            {
                //first character is residue
                char residue = modString.charAt(0);
                //second character is space.  Rest is mass
                double mass = Double.parseDouble(modString.substring(2));
                ModifiedAminoAcid mod = new ModifiedAminoAcid(residue, mass);
                modsThisPosition.add(mod);
            }
            result.put(position, modsThisPosition);
        }
        return result;
    }






    /**
     * Parse the value of a map of score types to score values.
     * This is intentionally brittle:  If there are extra "="s or ","s, or
     * if the score value doesn't parse as a Double, it'll throw a runtime exception
     * @param stringDoubleMapString
     * @return
     */
    public static Map<String,String>
            parseStringDoubleMapString(String stringDoubleMapString)
    {
        Map<String,String> result = new HashMap<String,String>();
        String[] stringArray =
                stringDoubleMapString.split(MULTI_VALUE_LIST_SEPARATOR);
        for (String string : stringArray)
        {
            String[] scoreNameAndValue = string.split("=");
            result.put(scoreNameAndValue[0],scoreNameAndValue[1]);
        }
        return result;
    }


    /**
     * If there are duplicate peptides next to each other, collapse them
     * @param feature
     */
    public static void collapsePeptideSequentialDuplicates(Feature feature)
    {
        collapseSequentialDuplicates(getPeptideList(feature));
    }

    /**
     * If there are duplicate proteins next to each other, collapse them
     * @param feature
     */
    public static void collapseProteinSequentialDuplicates(Feature feature)
    {
        collapseSequentialDuplicates(getProteinList(feature));
    }

    /**
     * If there are duplicate Strings next to each other, collapse them
     * @param peptideOrProteinList
     */
    public static void collapseSequentialDuplicates(List<String> peptideOrProteinList)
    {
        String previous = null;
        for (int i = peptideOrProteinList.size() - 1; i >= 0; i--)
        {
            String peptideOrProtein = peptideOrProteinList.get(i);
            if (peptideOrProtein != null && previous != null &&
                    peptideOrProtein.equals(previous))
            {
                peptideOrProteinList.remove(i);
            }
            else
            {
                previous = peptideOrProtein;
            }
        }
    }


    /**
     * Cover method, convenience
     * @param featureSet
     */
    public static void removeAllButFirstFeatureForEachPeptide(FeatureSet featureSet)
    {
        removeAllButFirstFeatureForEachPeptide(featureSet, false);
    }

    /**
     * As it says... removes all but the first feature (by scan) for each peptide.  If
     * we fail to find a peptide for ANY of this featureset's features, nothing happens.
     *
     * Also, the "firstscan" and "lastscan" values reflect a summary of all features for the peptide.
     * "firstscan" will be the same as "scan", since we're returning the first occurrence.  But
     * "lastscan" will be the scan of the last observation.
     *
     * Also, if perModificationState is specified, we'll try to keep one feature per
     * modification state, potentially multiple features per peptide sequence.  We do this
     * by checking if the modifications on each are _exactly_ the same.
     * @param featureSet
     * @param perModState If true, keep one feature per distinct modification state
     */
    public static void removeAllButFirstFeatureForEachPeptide(FeatureSet featureSet,
                                                              boolean perModState)
    {
        //map peptides to features
        HashMap<String,List<Feature>> peptideFeatureListMap =
                new HashMap<String,List<Feature>>();

        Feature[] features = featureSet.getFeatures();

        _log.debug("removeAllButFirstFeatureForEachPeptide, perModState=" + perModState + ", #features=" + features.length);

        for (Feature feature : features)
        {
            String peptide =
                    MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            List<Feature> thisPeptideFeatureList = peptideFeatureListMap.get(peptide);
            if (thisPeptideFeatureList == null)
            {
                thisPeptideFeatureList = new ArrayList<Feature>();
                peptideFeatureListMap.put(peptide, thisPeptideFeatureList);
            }
            thisPeptideFeatureList.add(feature);
        }

        List<Feature> resultFeatureList = new ArrayList<Feature>();

        //for debugging
        int sumNumModStatesPerPeptide = 0;

        for (String peptide : peptideFeatureListMap.keySet())
        {
            List<Feature> featuresThisPeptide = peptideFeatureListMap.get(peptide);

            if (perModState)
            {
                //only exact mod state matches fold together
                Map<String, List<Feature>> modStateFeatureListMap =
                        new HashMap<String, List<Feature>>();
                for (Feature featureThisPeptide : featuresThisPeptide)
                {
                    MS2ExtraInfoDef.getModifiedAminoAcids(featureThisPeptide);
                    Map<Integer, List<ModifiedAminoAcid>> positionModListMap =
                            (Map<Integer, List<ModifiedAminoAcid>>)
                                    featureThisPeptide.getProperty("modifiedaminoacids");
                    //special handling for no modifications
                    String modStateString = "";
                    if (positionModListMap != null && positionModListMap.size() > 0)
                        modStateString = MS2ExtraInfoDef.getSingletonInstance().convertToString("modifiedaminoacids", positionModListMap);
                    List<Feature> featuresThisModState = modStateFeatureListMap.get(modStateString);
                    if (featuresThisModState == null)
                    {
                        featuresThisModState = new ArrayList<Feature>();
                        modStateFeatureListMap.put(modStateString, featuresThisModState);
                    }
                    featuresThisModState.add(featureThisPeptide);
                }

                if (_log.isDebugEnabled())
                    sumNumModStatesPerPeptide += modStateFeatureListMap.size();
                //add the first for each mod state
                for (String modStateString : modStateFeatureListMap.keySet())
                {
                    List<Feature> featuresThisModState =
                            modStateFeatureListMap.get(modStateString);

                    int lastScan = 0;
                    Feature firstFeature = null;

                    for (Feature feature : featuresThisModState)
                    {
                        if (feature.getScan() > lastScan)
                            lastScan = feature.getScan();
                        if (firstFeature == null || feature.getScan() < firstFeature.getScan())
                            firstFeature = feature;
                    }
                    firstFeature.setScanLast(lastScan);
                    resultFeatureList.add(firstFeature);
                }
            }
            else
            {
                int lastScan = 0;
                Feature firstFeature = null;

                for (Feature feature : featuresThisPeptide)
                {
                    if (feature.getScan() > lastScan)
                        lastScan = feature.getScan();
                    if (firstFeature == null || feature.getScan() < firstFeature.getScan())
                        firstFeature = feature;
                }
                firstFeature.setScanLast(lastScan);
                resultFeatureList.add(firstFeature);
            }
        }

        if (perModState && _log.isDebugEnabled())
            _log.debug("\tmean mod states per peptide: " +
                   ((double) sumNumModStatesPerPeptide) / ((double) peptideFeatureListMap.size()));

        //sort all the earliest features for each peptide & set them as the feature array
        //for this featureset
        Collections.sort(resultFeatureList, new Feature.MzScanAscComparator());
        _log.debug("removeAllButFirstFeatureForEachPeptide, resulting features: " + resultFeatureList.size());

        featureSet.setFeatures(resultFeatureList.toArray(new Feature[0]));


    }

    /**
     * handy
     */
    public static class PeptideProphetDescComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            double ppDiff = getPeptideProphet(o1) - getPeptideProphet(o2);

             return ppDiff > 0 ? -1 :
                     ppDiff < 0 ? 1 : 0;
        }
    }

    /**
     * handy
     */
    public static class PeptideProphetAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            double ppDiff = getPeptideProphet(o2) - getPeptideProphet(o1);

             return ppDiff > 0 ? -1 :
                     ppDiff < 0 ? 1 : 0;
        }
    }




    public List<JMenuItem> createPopupMenuItems(Feature feature)
    {
        List<JMenuItem> result = new ArrayList<JMenuItem>();
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
        if (peptide != null)
        {
            JMenuItem googleSearchMenuItem =
                    new JMenuItem(TextProvider.getText("GOOGLE_PEPTIDE_SEARCH"));
            googleSearchMenuItem.setActionCommand(peptide);
            googleSearchMenuItem.addActionListener(
                    new AbstractAction ()
                    {
                        public void actionPerformed(ActionEvent e)
                        {
                            try
                            {
                                BrowserController.navigate("http://www.google.com/search?q=" + e.getActionCommand());
                            }
                            catch (Exception ex) {}
                        }
                    }
            );
            result.add(googleSearchMenuItem);
        }
        return result;
    }

}
