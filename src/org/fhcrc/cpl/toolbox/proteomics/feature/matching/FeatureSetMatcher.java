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

import java.util.*;


/**
 * Interface that must be implemented by any FeatureSetMatcher class.  Defines some
 * useful constants
 */
public interface FeatureSetMatcher
{
    //default delta masses for both types of mass tolerance
    public static final float DEFAULT_DELTA_MASS_ABSOLUTE = .2f;
    public static final float DEFAULT_DELTA_MASS_PPM = 5f;
    //constants used to specify how mass tolerance should be calculated
    public static final int DELTA_MASS_TYPE_ABSOLUTE = 0;
    public static final int DELTA_MASS_TYPE_PPM = 1;

    //default to PPM
    public static final int DEFAULT_DELTA_MASS_TYPE = DELTA_MASS_TYPE_PPM;

    //Default values for scan, time and hydrophobicity tolerance.  These values are
    //awfully meaningless... it's important for these to be defined appropriately for
    //the individual analysis
    public static final int DEFAULT_DELTA_SCAN = 3;
    public static final double DEFAULT_DELTA_TIME = 20;
    public static final double DEFAULT_DELTA_HYDROPHOBICITY = 0.05;

    /**
     * init method must define the parameters common to all FeatureSetMatchers
     * @param deltaMass
     * @param deltaMassType
     * @param deltaHydrophobicity
     */
    public void init(float deltaMass, int deltaMassType,
                     float deltaHydrophobicity);

    /**
     * do the actual featureset matching.
     * @return
     */
    public FeatureMatchingResult
            matchFeatures(FeatureSet masterSet, FeatureSet slaveSet);

    /**
     * This class represents the result of running
     * FeatureSetMatcher.matchFeatures().
     *
     * Pulling this out into its own class to make it utterly obvious what goes
     * where, and to make it easier to change in the future.
     *
     * This has already had to undergo one revision.  Initially, I returned
     * a list of pairs of features.  Now I'm returning (essentially) a Map.
     *
     * Each matched feature in the master feature set is associated with a list
     * of features in the slave feature set.  By convention, if a featuresetmatcher
     * is able to assign an order of likelihood to the matches, the list should be
     * in descending order of likelihood, i.e., best one first.
     */
    public static class FeatureMatchingResult
    {
        Map<Feature, List<Feature>> masterSlaveFeatureMap;

        public FeatureMatchingResult()
        {
            init();
        }

        protected void init()
        {
            masterSlaveFeatureMap = new HashMap<Feature, List<Feature>>();
        }

        public void put(Feature masterSetFeature, List<Feature> slaveSetFeatures)
        {
            masterSlaveFeatureMap.put(masterSetFeature, slaveSetFeatures);
        }

        public void add(Feature masterSetFeature, Feature slaveSetFeature)
        {
            List<Feature> existingSlaveSetList =
                    masterSlaveFeatureMap.get(masterSetFeature);

            if (existingSlaveSetList == null)
            {
                existingSlaveSetList = new ArrayList<Feature>();
                put(masterSetFeature, existingSlaveSetList);
            }

            existingSlaveSetList.add(slaveSetFeature);
        }

        public Set<Feature> getMasterSetFeatures()
        {
            return masterSlaveFeatureMap.keySet();
        }

        public Set<Feature> getSlaveSetFeatures()
        {
            Set<Feature> slaveSetFeatures = new HashSet<Feature>();

            for (List<Feature> slaveList : masterSlaveFeatureMap.values())
            {
                for (Feature slaveFeature : slaveList)
                    slaveSetFeatures.add(slaveFeature);
            }

            return slaveSetFeatures;
        }

        public List<Feature> getSlaveSetFeatures(Feature masterSetFeature)
        {
            return masterSlaveFeatureMap.get(masterSetFeature);
        }


        


        public List<Feature> get(Feature masterSetFeature)
        {
            return masterSlaveFeatureMap.get(masterSetFeature);
        }

        public Feature getBestMatch(Feature masterSetFeature)
        {
            List<Feature> matches = get(masterSetFeature);
            if (matches == null || matches.isEmpty())
                return null;
            //by convention, best match is first
            return matches.get(0);
        }

        public int size()
        {
            return masterSlaveFeatureMap.size();
        }

        public Map<Feature, List<Feature>> getMasterSlaveFeatureMap()
        {
            return masterSlaveFeatureMap;
        }
    }




}
