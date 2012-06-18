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
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;

import java.util.List;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Base implementing class for FeatureSetMatcher classes.  FeatureSetMatchers don't
 * have to extend from this class, it just provides some helpful methods.
 *
 * In particular, it makes it easy for implementing classes to order the feature
 * matches, and to use different elution modes (scan, time, hydrophobicity)
 */
public abstract class BaseFeatureSetMatcherImpl
    implements FeatureSetMatcher
{
    //modes for determining the ordering of slave set features.
    //feature "quality"
    public static final int ORDER_BY_FEATURE_QUALITY_MODE = 0;
    //closeness in observed hydrophobicity to the master feature
    public static final int ORDER_BY_ELUTION_CLOSENESS_TO_MASTER_MODE = 1;

    //default ordering of slave set features is by closeness in elution time
    protected int slaveSetOrderingMode = ORDER_BY_ELUTION_CLOSENESS_TO_MASTER_MODE;

    //Define different elution modes
    public static final int ELUTION_MODE_HYDROPHOBICITY=0;
    public static final int ELUTION_MODE_TIME=1;
    public static final int ELUTION_MODE_SCAN=2;

    //default elution mode is hydrophobicity
    public static final int DEFAULT_ELUTION_MODE=ELUTION_MODE_HYDROPHOBICITY;

    //determines whether this FSM works in terms of scan, time, or hydrophobicity
    protected int elutionMode = DEFAULT_ELUTION_MODE;

    //tolerances
    protected float deltaMass = DEFAULT_DELTA_MASS_PPM;
    protected int deltaMassType = DEFAULT_DELTA_MASS_TYPE;
    protected float deltaElution = (float) DEFAULT_DELTA_HYDROPHOBICITY;

    public BaseFeatureSetMatcherImpl()
    {

    }

    /**
     * Set the parameters that are common to all FeatureSetMatchers
     * @param deltaMass
     * @param deltaMassType
     * @param deltaHydrophobicity
     */
    public void init(float deltaMass, int deltaMassType,
                     float deltaHydrophobicity)
    {
        setDeltaMass(deltaMass);
        setDeltaMassType(deltaMassType);
        setDeltaElution(deltaHydrophobicity);
    }


    /**
     * Order the slave set features appropriately, in descending order of "goodness"
     *
     * @param slaveSetFeatures
     * @param masterSetFeature
     * @return a DISTINCT LIST from the input list, ordered
     */
    protected List<Feature> orderSlaveSetFeatures(List<Feature> slaveSetFeatures,
                                        Feature masterSetFeature)
    {
        Comparator<Feature> slaveOrderingComp;
//System.err.println("Master scan: " + masterSetFeature.getScan());
//System.err.println("sorting mode: " + slaveSetOrderingMode);
//System.err.println(" Unsorted slaves: " );
//for (Feature feature : slaveSetFeatures) System.err.println("  " + feature.getScan());

        switch(slaveSetOrderingMode)
        {
            //if we're ordering by elution closeness to master, then the appropriate
            //statistic to sort on is the same elution statistic we're using in
            //the overall matching
            case ORDER_BY_ELUTION_CLOSENESS_TO_MASTER_MODE:
                slaveOrderingComp =
                        new FeatureElutionClosenessToBaseFeatureComparatorDesc(masterSetFeature);
                break;
            case ORDER_BY_FEATURE_QUALITY_MODE:
                slaveOrderingComp = Feature.PeaksKLComparatorDesc.getSingletonInstance();
                break;
            //default to feature quality, but we shouldn't really hit this.
            //default needed for compilation.
            default:
                slaveOrderingComp = Feature.PeaksKLComparatorDesc.getSingletonInstance();
                break;
        }

        //copy the list, then sort
        List<Feature> result = new ArrayList<Feature>(slaveSetFeatures);
        Collections.sort(result, slaveOrderingComp);
//System.err.println( "Sorted slaves: " );
//for (Feature feature : result) System.err.println("  ***" + feature.getScan() + "***");
        return result;
    }

    /**
     * switch on elutionMode to return the appropriate elution value for this feature
     * @param feature
     * @return
     */
    protected double getElutionValue(Feature feature)
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

    /**
     * Compare features based on closeness in hydrophobicity to a base feature
     */
    protected class FeatureElutionClosenessToBaseFeatureComparatorDesc implements Comparator<Feature>
    {
        protected double baseFeatureValue = -999999;

        //be sure to update check in constructor whenever adding new types

        public FeatureElutionClosenessToBaseFeatureComparatorDesc(Feature baseFeature)
        {
             baseFeatureValue = getElutionValue(baseFeature);
        }

        public int compare(Feature o1, Feature o2)
        {
            double statDiff1 =
                    Math.abs(getElutionValue(o1) -
                             baseFeatureValue);
            double statDiff2 =
                    Math.abs(getElutionValue(o2) -
                             baseFeatureValue);

            return (statDiff1 == statDiff2 ? 0 : statDiff1 < statDiff2 ? -1 : 1);
        }
    }

    //getters and setters

    public float getDeltaMass()
    {
        return deltaMass;
    }

    public void setDeltaMass(float deltaMass)
    {
        this.deltaMass = deltaMass;
    }

    public int getDeltaMassType()
    {
        return deltaMassType;
    }

    public void setDeltaMassType(int deltaMassType)
    {
        this.deltaMassType = deltaMassType;
    }

    public float getDeltaElution()
    {
        return deltaElution;
    }

    public void setDeltaElution(float deltaElution)
    {
        this.deltaElution = deltaElution;
    }

    public int getSlaveSetOrderingMode()
    {
        return slaveSetOrderingMode;
    }

    public void setSlaveSetOrderingMode(int slaveSetOrderingMode)
    {
        this.slaveSetOrderingMode = slaveSetOrderingMode;
    }

    public void setElutionMode(int elutionMode)
    {
        this.elutionMode = elutionMode;
    }

    public int getElutionMode()
    {
        return elutionMode;
    }
}
