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
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;

import java.util.List;

/**
 * Contains column name and datatype information about each column.
 * Each column datatype must be a class with a constructor that accepts
 * one String argument.
 */
public class AmtExtraInfoDef extends FeatureExtraInformationDef
{
    static Logger _log = Logger.getLogger(AmtExtraInfoDef.class);

    public static final double NO_OBSERVED_HYDROPHOBICITY = -999999;
    public static final double NO_MATCH_PROBABILITY = -999999;
    public static final double NO_MATCH_FDR = -999999;



    public AmtExtraInfoDef()
    {
        super(
                    "AMT",
                    new String[]{
                            "observedhydrophobicity",
                            "match_probability",
                            "match_fdr"
                    },
                    new Class[]{
                            Double.class,
                            Double.class,
                            Double.class,
                    },
                    new String[]{
                            "db_runs_matched",
                            "estimated_true_positives",
                            "false_assignment_rate",
                            "matched_database_name"
                    }
            );

    }

    protected static AmtExtraInfoDef singletonInstance = null;

    public static AmtExtraInfoDef getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new AmtExtraInfoDef();
        return singletonInstance;
    }

    public static double getObservedHydrophobicity(Feature feature)
    {       
        return feature.getDoubleProperty("observedhydrophobicity",NO_OBSERVED_HYDROPHOBICITY);
    }

    public static void setObservedHydrophobicity(Feature feature, double hydrophobicity)
    {
        feature.setProperty("observedhydrophobicity", hydrophobicity);
    }

    public static boolean hasObservedHydrophobicity(Feature feature)
    {
        return getObservedHydrophobicity(feature) != NO_OBSERVED_HYDROPHOBICITY;
    }

    public static double getMatchProbability(Feature feature)
    {
        return feature.getDoubleProperty("match_probability",NO_MATCH_PROBABILITY);
    }

    public static void setMatchProbability(Feature feature, double matchProbability)
    {
        feature.setProperty("match_probability", matchProbability);
    }

    public static boolean hasMatchProbability(Feature feature)
    {
        return getMatchProbability(feature) != NO_MATCH_PROBABILITY;
    }

    public static double getMatchFDR(Feature feature)
    {
        return feature.getDoubleProperty("match_fdr",NO_MATCH_FDR);
    }

    public static void setMatchFDR(Feature feature, double matchFDR)
    {
        feature.setProperty("match_fdr", matchFDR);
    }

    public static boolean hasMatchFDR(Feature feature)
    {
        return getMatchFDR(feature) != NO_MATCH_FDR;
    }




    /**
     * Same as convertToString, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public String convertFeatureSetPropertyToString(String propertyName, Object value)
    {
        if (propertyName.equals("db_runs_matched"))
        {
            return convertStringListToString((List<String>) value);
        }
        else if (propertyName.equals("estimated_true_positives") ||
                 propertyName.equals("false_assignment_rate") ||
                 propertyName.equals("matched_database_name"))
        {
            return value.toString();
        }
        else throw new IllegalArgumentException(
                "AmtExtraInfoDef doesn't know about a feature set property named " +
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
        if (propertyName.equals("db_runs_matched"))
        {
            List<String> runsMatchedStrings =
                    MS2ExtraInfoDef.parseStringListString(value);
            return runsMatchedStrings;
        }
        else if (propertyName.equals("estimated_true_positives") ||
                 propertyName.equals("false_assignment_rate"))
        {
            return Double.parseDouble(value);
        }
        else if (propertyName.equals("matched_database_name"))
        {
            return value;
        }
        else throw new IllegalArgumentException(
                "AmtExtraInfoDef doesn't know about a feature set property named " +
                        propertyName);
    }

    public static List<String> getFeatureSetRunsMatched(FeatureSet featureSet)
    {
        return (List<String>)
                getSingletonInstance().getFeatureSetProperty(featureSet, "db_runs_matched");
    }

    public static void setFeatureSetRunsMatched(FeatureSet featureSet,
                                                List<String> runsMatched)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "db_runs_matched", runsMatched);
    }

    public static double getFeatureSetEstimatedTruePositives(FeatureSet featureSet)
    {
        return (Double) getSingletonInstance().getFeatureSetProperty(featureSet,
                "estimated_true_positives");
    }

    public static void setFeatureSetEstimatedTruePositives(FeatureSet featureSet,
                                                             double estTruePos)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "estimated_true_positives",
                estTruePos);
    }

    public static double getFeatureSetFalseAssignmentRate(FeatureSet featureSet)
    {
        return (Double) getSingletonInstance().getFeatureSetProperty(featureSet,
                "false_assignment_rate");
    }

    public static void setFeatureSetFalseAssignmentRate(FeatureSet featureSet,
                                                             double estTruePos)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "false_assignment_rate",
                estTruePos);
    }

    public static String getFeatureSetMatchedDatabaseName(FeatureSet featureSet)
    {
        return (String) getSingletonInstance().getFeatureSetProperty(featureSet,
                "matched_database_name");
    }

    public static void setFeatureSetMatchedDatabaseName(FeatureSet featureSet,
                                                        String matchedDatabaseName)
    {
        getSingletonInstance().setFeatureSetProperty(featureSet, "matched_database_name",
                matchedDatabaseName);
    }

}
