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
package org.fhcrc.cpl.toolbox.proteomics.feature.matching;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Attempts to match one feature set against another using a window of minimum and maximum
 * mass and elution time differences, not necessarily centered at 0.
 */
public class Window2DFeatureSetMatcher
{
    private static Logger _log = Logger.getLogger(Window2DFeatureSetMatcher.class);

    protected float minMassDiff, maxMassDiff;
    protected float minElutionDiff, maxElutionDiff;

    //Allowed values: FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE and FeatureSetMatcher.DELTA_MASS_TYPE_PPM
    protected int massDiffType = FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE;

    //Define different elution modes
    public static final int ELUTION_MODE_HYDROPHOBICITY=0;
    public static final int ELUTION_MODE_TIME=1;
    public static final int ELUTION_MODE_SCAN=2;

    //default elution mode is hydrophobicity
    public static final int DEFAULT_ELUTION_MODE=ELUTION_MODE_HYDROPHOBICITY;

    //determines whether this FSM works in terms of scan, time, or hydrophobicity
    protected int elutionMode = DEFAULT_ELUTION_MODE;

    //Should we consider elution differences based on a point (e.g., apex intensity time) or
    //a range (e.g., minScan-maxScan)?
    public static final int ELUTION_RANGE_MODE_POINT=0;
    public static final int ELUTION_RANGE_MODE_RANGE=1;
    //NOTE: because we only store point values for time and hydrophobicity, MODE_RANGE only makes sense for MODE_SCAN

    protected int elutionRangeMode = ELUTION_RANGE_MODE_POINT;

    //match features only if they have the same charge?
    protected boolean matchWithinChargeOnly = false;


    public Window2DFeatureSetMatcher()
    {

    }

    public String toString()
    {
        return "Window2DFeatureSetMatcher: massDiffRange=(" + minMassDiff + "," + maxMassDiff +
               "), elutionDiffRange=(" + minElutionDiff + ", " + maxElutionDiff + ")";
    }

    /**
     * Convenience method to set a bunch of parameters at once

     */
    public void setMatchingParameters(float minMassDiff, float maxMassDiff,
                                      float minElutionDiff, float maxElutionDiff,
                                      int massDiffType)
    {
        setMinMassDiff(minMassDiff);
        setMaxMassDiff(maxMassDiff);
        setMinElutionDiff(minElutionDiff);
        setMaxElutionDiff(maxElutionDiff);

        setMassDiffType(massDiffType);
    }

    /**
     * Match the features in ms2Features against ms1Features.  Populate the passed-in unmatchedMS2Features
     * ArrayList with all unmatched features.
     * Can also be used to match ms1 features against other ms1 features
     * @return
     */
    public FeatureSetMatcher.FeatureMatchingResult matchFeatures(FeatureSet masterFeatures,
                                               FeatureSet slaveFeatures)
    {
        FeatureSetMatcher.FeatureMatchingResult result =
                new FeatureSetMatcher.FeatureMatchingResult();

        //clone the feature array, because we'll be sorting it differently from how
        //FeatureSet expects it to be sorted
        Feature[] masterFeatureArray = masterFeatures.getFeatures();

        Feature.MassAscComparator comp = new Feature.MassAscComparator();

        Feature[] slaveFeatureArray = slaveFeatures.getFeatures().clone();

        Arrays.sort(slaveFeatureArray, comp);

        for (Feature masterFeature : masterFeatureArray)
        {

            List<Feature> matchedSlaveFeatures =
                findMatchingFeatures(slaveFeatureArray, masterFeature, comp);

            if (matchedSlaveFeatures != null &&
                !matchedSlaveFeatures.isEmpty())
            {
                //order slave set features in the order determined by
                //variable settings at the superclass level
                result.put(masterFeature, matchedSlaveFeatures);
            }

        }
        return result;
    }

    /**
     * Find the single slave features that most closely match this Master feature
     * @return the index of the nearest feature, or -1 if none found in the range
     */
    public List<Feature> findMatchingFeatures(Feature[] slaveFeatures, Feature masterFeature,
                                                       Comparator<Feature> featureComparator)
    {
        int index = Arrays.binarySearch(slaveFeatures, masterFeature, featureComparator);
        int centerPos = index;
        //if no exact match, restore the would-have-been position of the feature in the array,
        //so we can work outward from there in our search
        if (index < 0)
            centerPos = -index - 1;
//System.err.println("fMF, centerPos=" + centerPos + ", slaveAtCP: " + slaveFeatures[centerPos].getMass());

        //calculate the absolute min and max mass diff allowed at this feature mass
        float absoluteMinMassDiff =
                MassUtilities.calculateAbsoluteDeltaMass(masterFeature.mass,
                                                        minMassDiff,
                                                        massDiffType);
        //calculate the absolute delta mass allowed at this feature mass
        float absoluteMaxMassDiff =
                MassUtilities.calculateAbsoluteDeltaMass(masterFeature.mass,
                                                        maxMassDiff,
                                                        massDiffType);
//System.err.println("Slave features length: " + slaveFeatures.length + ", centerPos = " + centerPos);
        boolean doneUp = false;
        boolean doneDown = false;
        boolean searchingUp = false;
        int position=-1;

        List<Feature> result = new ArrayList<Feature>();

        for (int i=0; !(doneUp && doneDown); i++)
        {
            int offset = i/2;

            if (i % 2 == 1)
            {
                //search with decreasing index, below the starting mass
                offset = -offset - 1;
                position = centerPos + offset;
                searchingUp = false;
                if (position < 0)
                    doneDown = true;
                if (doneDown)
                    continue;
            }
            else
            {
                //search with increasing index, above the starting mass
                position = centerPos + offset;
                searchingUp = true;
                if (position >= slaveFeatures.length)
                    doneUp = true;
                if (doneUp)
                    continue;
            }
            Feature slaveFeature = slaveFeatures[position];

            float massDiff = masterFeature.mass - slaveFeature.mass;

            //since we've sorted primarily by mass, if we've gone out of range
            //we can stop searching in this direction.
            //searching up means negative massDiff values, searching down means positive
            if (!searchingUp && massDiff > absoluteMaxMassDiff)
            {
                doneDown=true;
                continue;
            }
            else if(searchingUp && massDiff < absoluteMinMassDiff)
            {
                doneUp=true;
                continue;
            }

            float elutionDiff = calcElutionDiff(masterFeature, slaveFeature);

            if (elutionDiff >= minElutionDiff &&
                elutionDiff <= maxElutionDiff)
            {
//System.err.println("Adding match: searchingup? " + searchingUp + ", master feature mass=" + masterFeature.getMass() + ", slave mass=" + slaveFeature.getMass() + ", massDiff = " + massDiff + ", minMD=" + absoluteMinMassDiff + ", maxMD=" + absoluteMaxMassDiff);
//System.err.println("\tcp+1: " + slaveFeatures[Math.min(centerPos+1,slaveFeatures.length-1)].getMass() +", cp: " + slaveFeatures[centerPos].getMass() + ", cp-1: " + slaveFeatures[Math.max(0,centerPos-1)].getMass());

                //if we're matching only within charge, check charge
                if (!matchWithinChargeOnly ||
                    masterFeature.getCharge() == slaveFeature.getCharge())
                {
                    result.add(slaveFeature);
                }
            }
        }

        return result;
    }

    protected float calcElutionDiff(Feature feature1, Feature feature2)
    {
        switch (elutionRangeMode)
        {
            case ELUTION_RANGE_MODE_POINT:
                return (float) (getPointElutionValue(feature1) - getPointElutionValue(feature2));
            case ELUTION_RANGE_MODE_RANGE:
                float f1Min = (float) getMinElutionRangeValue(feature1);
                float f1Max = (float) getMaxElutionRangeValue(feature1);

                float f2Min = (float) getMinElutionRangeValue(feature2);
                float f2Max = (float) getMaxElutionRangeValue(feature2);

                if (f1Min > f2Min)
                    return Math.max(0, f1Min - f2Max);
                else if (f2Max > f1Max)
                    return Math.max(0, f2Min - f1Max);
                return 0;
            default:
                throw new IllegalArgumentException("Unknown elution range mode");
        }
    }

    /**
     * switch on elutionMode to return the appropriate elution value for this feature
     * @param feature
     * @return
     */
    protected double getMinElutionRangeValue(Feature feature)
    {
        switch (elutionMode)
        {
            case ELUTION_MODE_HYDROPHOBICITY:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
            case ELUTION_MODE_TIME:
                return feature.getTime();
            case ELUTION_MODE_SCAN:
                return feature.getScanFirst();
            default:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
        }
    }

    /**
     * switch on elutionMode to return the appropriate elution value for this feature
     * @param feature
     * @return
     */
    protected double getMaxElutionRangeValue(Feature feature)
    {
        switch (elutionMode)
        {
            case ELUTION_MODE_HYDROPHOBICITY:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
            case ELUTION_MODE_TIME:
                return feature.getTime();
            case ELUTION_MODE_SCAN:
                return feature.getScanLast();
            default:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
        }
    }


    /**
     * switch on elutionMode to return the appropriate elution value for this feature
     * @param feature
     * @return
     */
    protected double getPointElutionValue(Feature feature)
    {
        switch (elutionMode)
        {
            case ELUTION_MODE_HYDROPHOBICITY:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
            case ELUTION_MODE_TIME:
                return feature.getTime();
            case ELUTION_MODE_SCAN:
                return feature.getScan();
            default:
                return AmtExtraInfoDef.getObservedHydrophobicity(feature);
        }
    }

    public float getMinMassDiff()
    {
        return minMassDiff;
    }

    public void setMinMassDiff(float minMassDiff)
    {
        this.minMassDiff = minMassDiff;
    }

    public float getMaxMassDiff()
    {
        return maxMassDiff;
    }

    public void setMaxMassDiff(float maxMassDiff)
    {
        this.maxMassDiff = maxMassDiff;
    }

    public float getMinElutionDiff()
    {
        return minElutionDiff;
    }

    public void setMinElutionDiff(float minElutionDiff)
    {
        this.minElutionDiff = minElutionDiff;
    }

    public float getMaxElutionDiff()
    {
        return maxElutionDiff;
    }

    public void setMaxElutionDiff(float maxElutionDiff)
    {
        this.maxElutionDiff = maxElutionDiff;
    }


    public int getMassDiffType()
    {
        return massDiffType;
    }

    public void setMassDiffType(int massDiffType)
    {
        this.massDiffType = massDiffType;
    }


    public int getElutionRangeMode()
    {
        return elutionRangeMode;
    }

    public void setElutionRangeMode(int elutionRangeMode)
    {
        this.elutionRangeMode = elutionRangeMode;
    }


    public int getElutionMode()
    {
        return elutionMode;
    }

    public void setElutionMode(int elutionMode)
    {
        this.elutionMode = elutionMode;
    }


    public boolean isMatchWithinChargeOnly()
    {
        return matchWithinChargeOnly;
    }

    public void setMatchWithinChargeOnly(boolean matchWithinChargeOnly)
    {
        this.matchWithinChargeOnly = matchWithinChargeOnly;
    }
}
