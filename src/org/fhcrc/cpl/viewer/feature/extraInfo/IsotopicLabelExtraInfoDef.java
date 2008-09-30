/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.feature.extraInfo;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.feature.AnalyzeICAT;
import org.fhcrc.cpl.viewer.feature.Feature;

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
                        "heavyMass"
                },
                new Class[]{
                        Float.class, Float.class, Float.class,
                        Integer.class, AnalyzeICAT.IsotopicLabel.class, Float.class, Float.class
                }
            );
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
}
