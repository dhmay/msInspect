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
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.PeptideUtilities;
import org.fhcrc.cpl.toolbox.proteomics.QuantitationUtilities;

import java.util.List;

/**
 * Contains column name and datatype information about each column.
 * Each column datatype must be a class with a constructor that accepts
 * one String argument.
 */
public class IsotopicLabelExtraInfoDef extends FeatureExtraInformationDef
{
    static Logger _log = Logger.getLogger(IsotopicLabelExtraInfoDef.class);


    public static final float NO_RATIO_FOR_FEATURE = -1;

    public IsotopicLabelExtraInfoDef()
    {
        super(
                "ISOTOPIC_LABEL",
                new String[]{
                        "lightIntensity",
                        "heavyIntensity",
                        "ratio",
                        "labelCount",
                        "label",
                        "lightMass",
                        "heavyMass",
                        "lightFirstScan",
                        "lightLastScan",
                        "heavyFirstScan",
                        "heavyLastScan",
                },
                new Class[]{
                        Float.class, Float.class, Float.class,
                        Integer.class, AnalyzeICAT.IsotopicLabel.class, Float.class, Float.class,
                        Integer.class, Integer.class, Integer.class, Integer.class
                },
                new String[]{
                        "algorithm"
                }
            );
    }

    /**
     * Same as convertToString, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public String convertFeatureSetPropertyToString(String propertyName, Object value)
    {
        if (propertyName.equals("algorithm"))
        {
            return (String) value;
        }
        else throw new IllegalArgumentException(
                "IsotopicLabelExtraInfoDef doesn't know about a feature set property named " +
                propertyName);
    }

    /**
     * Save as convertStringValue, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public Object convertFeatureSetPropertyStringValue(String propertyName, String value)
    {
        if (propertyName.equals("algorithm"))
        {
            if (QuantitationUtilities.ALGORITHM_Q3.equals(value) || QuantitationUtilities.ALGORITHM_XPRESS.equals(value))
                return value;
            else throw new IllegalArgumentException(
                "IsotopicLabelExtraInfoDef doesn't know about a quantitation algorithm named  " +
                propertyName);
        }
        else throw new IllegalArgumentException(
                "IsotopicLabelExtraInfoDef doesn't know about a feature set property named " +
                propertyName);
    }

    protected static IsotopicLabelExtraInfoDef singletonInstance = null;

    public static IsotopicLabelExtraInfoDef getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new IsotopicLabelExtraInfoDef();
        return singletonInstance;
    }

    public static boolean hasLabel(Feature feature)
    {
        return (feature.hasProperty("label") && (feature.getProperty("label") != null));
    }

    public static int getLabelCount(Feature feature)
    {
        return feature.getIntProperty("labelCount",0);
    }

    public static AnalyzeICAT.IsotopicLabel getLabel(Feature feature)
    {
        return (AnalyzeICAT.IsotopicLabel) feature.getProperty("label");
    }

    public static void setLabel(Feature feature, AnalyzeICAT.IsotopicLabel label)
    {
        feature.setProperty("label", label);
    }

    public static void setHeavyIntensity(Feature feature, double intensity)
    {
        feature.setProperty("heavyIntensity", intensity);
    }

    public static void setLightIntensity(Feature feature, double intensity)
    {
        feature.setProperty("lightIntensity", intensity);
    }

    public static double getHeavyIntensity(Feature feature)
    {
        return feature.getDoubleProperty("heavyIntensity", 0);
    }

    public static double getLightIntensity(Feature feature)
    {
        return feature.getDoubleProperty("lightIntensity", 0);
    }

    public static void setHeavyMass(Feature feature, double mass)
    {
        feature.setProperty("heavyMass", mass);
    }

    public static double getHeavyMass(Feature feature)
    {
        return feature.getDoubleProperty("heavyMass", 0);
    }

    public static void setLightMass(Feature feature, double mass)
    {
        feature.setProperty("lightMass", mass);
    }

    public static double getLightMass(Feature feature)
    {
        return feature.getDoubleProperty("lightMass", 0);
    }

    public static void setLightFirstScan(Feature feature, int scan)
    {
        feature.setProperty("lightFirstScan", scan);
    }

    public static int getLightFirstScan(Feature feature)
    {
        return feature.getIntProperty("lightFirstScan", 0);
    }

    public static void setLightLastScan(Feature feature, int scan)
    {
        feature.setProperty("lightLastScan", scan);
    }

    public static int getLightLastScan(Feature feature)
    {
        return feature.getIntProperty("lightLastScan", 0);
    }

    public static void setHeavyFirstScan(Feature feature, int scan)
    {
        feature.setProperty("heavyFirstScan", scan);
    }

    public static int getHeavyFirstScan(Feature feature)
    {
        return feature.getIntProperty("heavyFirstScan", 0);
    }

    public static void setHeavyLastScan(Feature feature, int scan)
    {
        feature.setProperty("heavyLastScan", scan);
    }

    public static int getHeavyLastScan(Feature feature)
    {
        return feature.getIntProperty("heavyLastScan", 0);
    }



    public static void setRatio(Feature feature, double ratio)
    {
        feature.setProperty("ratio", ratio);
    }

    public static void removeRatio(Feature feature)
    {
        feature.unsetProperty("ratio");
    }

    public static double getRatio(Feature feature)
    {
        try
        {
            return feature.getFloatProperty("ratio",NO_RATIO_FOR_FEATURE);
        }
        catch (Exception e)
        {
            return feature.getDoubleProperty("ratio",NO_RATIO_FOR_FEATURE);
        }
    }

    public static boolean hasRatio(Feature feature)
    {
        return getRatio(feature) != NO_RATIO_FOR_FEATURE;
    }

    public static void setLabelCount(Feature feature, int labelCount)
    {
        feature.setProperty("labelCount", labelCount);
    }

    public static String getFeatureSetAlgorithm(FeatureSet featureSet)
    {
        return (String) getSingletonInstance().getFeatureSetProperty(featureSet, "algorithm");
    }

    public static void setFeatureSetAlgorithm(FeatureSet featureSet,
                                                       String baseName)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "algorithm", baseName);
    }

    public static boolean isLightLabeled(Feature feature, int labelType)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean lightLabeled = false;
        switch(labelType)
        {
            case QuantitationUtilities.LABEL_ACRYLAMIDE:
                lightLabeled =  PeptideUtilities.checkForModAllResidues(peptide, mods, 'C', QuantitationUtilities.ACRYLAMIDE_LABEL_LIGHTMASS);
                break;
            case QuantitationUtilities.LABEL_LYCINE:
                lightLabeled =  PeptideUtilities.checkForModNoResidues(peptide, mods, 'K', QuantitationUtilities.SILAC_LABEL_MASS);
                break;
            case QuantitationUtilities.LABEL_LYCINE_ARGININE:
                lightLabeled = (!peptide.contains("K") ||
                        PeptideUtilities.checkForModNoResidues(peptide, mods, 'K', QuantitationUtilities.SILAC_LABEL_MASS)) &&
                        (!peptide.contains("R") ||
                        PeptideUtilities.checkForModNoResidues(peptide, mods, 'R', QuantitationUtilities.SILAC_LABEL_MASS));
                break;
        }
        return lightLabeled;
    }

    public static boolean isHeavyLabeled(Feature feature, int labelType)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean heavyLabeled = false;
        switch(labelType)
        {
            case QuantitationUtilities.LABEL_ACRYLAMIDE:
                heavyLabeled = PeptideUtilities.checkForModAllResidues(peptide, mods, 'C', QuantitationUtilities.ACRYLAMIDE_LABEL_HEAVYMASS);
                break;
            case QuantitationUtilities.LABEL_LYCINE:
                heavyLabeled = PeptideUtilities.checkForModAllResidues(peptide, mods, 'K', QuantitationUtilities.SILAC_LABEL_MASS);
                break;
            case QuantitationUtilities.LABEL_LYCINE_ARGININE:
                heavyLabeled = (!peptide.contains("K") ||
                        PeptideUtilities.checkForModAllResidues(peptide, mods, 'K', QuantitationUtilities.SILAC_LABEL_MASS)) &&
                        (!peptide.contains("$") ||
                        PeptideUtilities.checkForModAllResidues(peptide, mods, 'R', QuantitationUtilities.SILAC_LABEL_MASS));
                break;
        }
        return heavyLabeled;
    }

}
