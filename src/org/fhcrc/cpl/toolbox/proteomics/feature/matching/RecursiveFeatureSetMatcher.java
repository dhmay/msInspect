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
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.BaseFeatureSetMatcherImpl;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Attempts to match one feature set against another.
 * This version is the closest to what CS types would call hierarchical clustering:
 * we don't attempt to guess the global "best" values for every bucket.  Instead, we
 * start loose, make what matches we can, and drill down into the smaller buckets, breaking
 * them up until either everything's cool or we hit a defined minimum in both dimensions.
 *
 * Strictly speaking, though, this isn't exactly hierarchical clustering, either.  Instead of
 * breaking clusters in half each time, I'm decreasing the cluster size slightly and re-clustering.
 *
 * In order to do true hierarchical clustering appropriately, it would be necessary to rewrite
 * Clusterer2D.  I think this is a close enough approximation that that's not a high priority.
 * Since it would affect array creation, too, I'm reluctant to do it.
 */
public class RecursiveFeatureSetMatcher extends BaseFeatureSetMatcherImpl
        implements FeatureSetMatcher
{
    private static Logger _log = Logger.getLogger(RecursiveFeatureSetMatcher.class);

    public static final double DEFAULT_HYDRO_ELUTION_BUCKET_INCREMENT = 0.01;
    public static final double DEFAULT_SCAN_ELUTION_BUCKET_INCREMENT = 5;
    public static final double DEFAULT_TIME_ELUTION_BUCKET_INCREMENT = 10;

    protected static final double ABSOLUTE_MASS_BUCKET_INCREMENT = 0.05;
    protected static final double PPM_MASS_BUCKET_INCREMENT = 2;
    protected double massBucketIncrement = ABSOLUTE_MASS_BUCKET_INCREMENT;
    protected double elutionBucketIncrement = DEFAULT_HYDRO_ELUTION_BUCKET_INCREMENT;

    public static final double DEFAULT_MIN_DELTA_MASS_PPM = 1;
    public static final double DEFAULT_MIN_DELTA_MASS_DA = .05;

    public static final double DEFAULT_MIN_DELTA_SCAN = 5;
    public static final double DEFAULT_MIN_DELTA_TIME = 10;
    public static final double DEFAULT_MIN_DELTA_H = .005;
    protected double minDeltaMass = DEFAULT_MIN_DELTA_MASS_PPM;
    protected double minDeltaElution = DEFAULT_MIN_DELTA_SCAN;

    protected int[] matchedAtDepth;

    protected boolean useMassInsteadOfMz = true;



    public RecursiveFeatureSetMatcher()
    {
    }

    public RecursiveFeatureSetMatcher(float deltaMass, int deltaMassType,
                                      float deltaElution)
    {
        super();
        init(deltaMass, deltaMassType, deltaElution);
    }

    public void init(float deltaMass, int deltaMassType, float deltaElution)
    {
        super.init(deltaMass, deltaMassType, deltaElution);
        if (deltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_PPM)
            massBucketIncrement = PPM_MASS_BUCKET_INCREMENT;
    }


    public String toString()
    {
        return "RecursiveFeatureSetMatcher: deltaMass=" + deltaMass +
                ", deltaElution=" + deltaElution;
    }

    /**
     *
     * //      In the PPM case, we have to convert from PPM to Da before figuring out whether we can split.
     //      This is done by converting maxBucket based on this bucket's getMin() value.  This will
     //      ensure that a bucket will be split if any subregion of the bucket would be split for the
     //      given PPM value.
     //      This should be reasonable behavior, perhaps splitting slightly more than necessary.
     * @param clusterer
     */
    protected void setDimensionSplitCalculator(FeatureClusterer clusterer)
    {
        if (deltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_PPM)
        {
            clusterer.setDimensionSplitCalculator(new Clusterer2D.ClusterDimensionSplitCalculator()
            {
                public double calculateDimension1ForSplit(double referenceValue, double dimensionValue)
                {
                    return MassUtilities.calculateAbsoluteDeltaMass((float) referenceValue,
                            (float) dimensionValue,
                            FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
                }

                public double calculateDimension2ForSplit(double referenceValue, double dimensionValue)
                {
                    return dimensionValue;
                }
            });
        }
    }


    /**
     * Match two feature sets using clustering
     * @param featureSet1
     * @param featureSet2
     * @return
     */
    public FeatureMatchingResult
    matchFeatures(FeatureSet featureSet1, FeatureSet featureSet2)
    {
        _log.debug("matchFeatures 1");

        List<Feature> set1Features =
                new ArrayList<Feature>(featureSet1.getFeatures().length);
        List<Feature> set2Features =
                new ArrayList<Feature>(featureSet2.getFeatures().length);

        for (Feature feature : featureSet1.getFeatures())
            set1Features.add(feature);
        for (Feature feature : featureSet2.getFeatures())
            set2Features.add(feature);

        matchedAtDepth = new int[100];

        FeatureMatchingResult result = recursivelyMatch(set1Features,set2Features,
                deltaMass, deltaElution,
                1);
        _log.debug("Done.  Features matched at depth...");
        for (int i=1; i<matchedAtDepth.length; i++)
        {
            if (matchedAtDepth[i] == 0)
                break;
            _log.debug("\t" + i + "  :  " + matchedAtDepth[i]);
        }

        return result;
    }

    protected FeatureMatchingResult recursivelyMatch(List<Feature> set1Features,
                                                     List<Feature> set2Features,
                                                     double currentDeltaMass,
                                                     double currentDeltaElution,
                                                     int depth)
    {
        FeatureMatchingResult result = new FeatureMatchingResult();

        boolean adjustedMass=false;
        boolean adjustedElution=false;
        if (currentDeltaMass < minDeltaMass)
        {
            currentDeltaMass = minDeltaMass;
            adjustedMass=true;
        }
        if (currentDeltaElution < minDeltaElution)
        {
            currentDeltaElution = minDeltaElution;
            adjustedElution=true;
        }
        if (adjustedMass && adjustedElution)
            return result;

        FeatureSet currentFeatureSet1 =
                new FeatureSet(set1Features.toArray(
                        new Feature[set1Features.size()]));
        FeatureSet currentFeatureSet2 =
                new FeatureSet(set2Features.toArray(
                        new Feature[set2Features.size()]));

        FeatureClusterer clusterer =
                new FeatureClusterer(FeatureClusterer.MASS_MZ_MODE_MASS,
                        elutionMode, currentFeatureSet1);
        clusterer.addSet(currentFeatureSet2);
        setDimensionSplitCalculator(clusterer);

        clusterer.split2D(currentDeltaMass, currentDeltaElution);
        Clusterer2D.BucketSummary[] bucketSummaries = clusterer.summarize();

        int numConflictBuckets = 0;
        int matchedThisIteration = 0;
        for (Clusterer2D.BucketSummary bucketSummary : bucketSummaries)
        {
            //for each bucket, if it has at least one pair, do something with it
            if (bucketSummary.setCount == 2)
            {
                if (bucketSummary.entries().length == 2)
                {
                    //exactly one from each, no-brainer
                    Feature masterSetFeature = ((FeatureClusterer.FeatureClusterable)
                            bucketSummary.getParentListForSetIndex(0).get(0)).getParentFeature();
                    Feature slaveSetFeature = ((FeatureClusterer.FeatureClusterable)
                            bucketSummary.getParentListForSetIndex(1).get(0)).getParentFeature();
                    result.add(masterSetFeature, slaveSetFeature);
                    matchedAtDepth[depth]++;
                }
                else
                {
                    //Conflict bucket.  Try to resolve recursively.  If failure, add as is
                    List<Feature> set2FeaturesThisBucket =
                            FeatureGrouper.getFeaturesFromSet(1, bucketSummary);
                    List<Feature> set1FeaturesThisBucket =
                            FeatureGrouper.getFeaturesFromSet(0, bucketSummary);

                    double nextDeltaMass = currentDeltaMass;
                    double nextDeltaElution = currentDeltaElution;
                    nextDeltaMass -= massBucketIncrement;
                    nextDeltaElution -= elutionBucketIncrement;

                    FeatureMatchingResult subResult =
                            recursivelyMatch(set1FeaturesThisBucket,
                                    set2FeaturesThisBucket,
                                    nextDeltaMass,
                                    nextDeltaElution,
                                    depth+1);

                    if (subResult.size() == 0)
                    {
                        //if we got nothing from the next level, make the matches at this level

                        for (Feature masterSetFeature : set1FeaturesThisBucket)
                        {
                            //order slave set features in the order determined by
                            //variable settings at the superclass level
                            List<Feature> orderedSlaveSetFeatures =
                                    orderSlaveSetFeatures(set2FeaturesThisBucket, masterSetFeature);
                            result.put(masterSetFeature, orderedSlaveSetFeatures);
                            matchedAtDepth[depth]++;
                        }
                    }
                    else
                    {
                        //if we got some stuff at the next level, add that stuff, then
                        //go back and pick up the spares at this level

                        Set<Feature> subResultMasterSetFeatures = subResult.getMasterSetFeatures();
                        Set<Feature> subResultSlaveSetFeatures = subResult.getSlaveSetFeatures();
                        for (Feature feature : subResultMasterSetFeatures)
                        {
                            result.put(feature, subResult.get(feature));
                        }

                        Set<Feature> remainingMasterSetFeatures = new HashSet<Feature>();

                        for (Feature feature : set1FeaturesThisBucket)
                        {
                            if (!subResultMasterSetFeatures.contains(feature))
                               remainingMasterSetFeatures.add(feature);
                        }

                        if (!remainingMasterSetFeatures.isEmpty())
                        {
                            List<Feature> remainingSlaveSetFeatures = new ArrayList<Feature>();
                            for (Feature feature : set2FeaturesThisBucket)
                                if (!subResultSlaveSetFeatures.contains(feature))
                                    remainingSlaveSetFeatures.add(feature);
                            if (!remainingSlaveSetFeatures.isEmpty())
                            {
                                for (Feature masterSetFeature : remainingMasterSetFeatures)
                                {
                                    List<Feature> orderedSlaveSetFeatures =
                                            orderSlaveSetFeatures(set2FeaturesThisBucket,
                                                                  masterSetFeature);
                                    result.put(masterSetFeature, orderedSlaveSetFeatures);
                                    matchedAtDepth[depth]++;
                                }
                            }
                        }
                    }
                    numConflictBuckets++;
                }
            }
            else
            {
                //if bucket contains only features from one set, those features are
                //unmatchable.
            }

        }

//            _log.debug("\t\tDone. matched: " + matchedThisIteration +
//                    ", Conflict buckets: " + numConflictBuckets);


        return result;
    }




    /**
     * Convert a list of Clusterables, assumed to be FeatureMassHydroClusterables, to a list of
     * Features
     * @param clusterableList
     * @return
     */
    protected List<Feature>
    getFeatureListForClusterableList(List<Clusterer2D.Clusterable> clusterableList)
    {
        List<Feature> result = new ArrayList<Feature>(clusterableList.size());
        for (Clusterer2D.Clusterable clusterable : clusterableList)
        {
            result.add(((FeatureClusterer.FeatureClusterable) clusterable).getParentFeature());
        }
        return result;
    }



    public double getMassBucketIncrement()
    {
        return massBucketIncrement;
    }

    public void setMassBucketIncrement(double massBucketIncrement)
    {
        this.massBucketIncrement = massBucketIncrement;
    }

    public double getElutionBucketIncrement()
    {
        return elutionBucketIncrement;
    }

    public void setElutionBucketIncrement(double elutionBucketIncrement)
    {
        this.elutionBucketIncrement = elutionBucketIncrement;
    }

    public void setUseMassInsteadOfMz(boolean value)
    {
        useMassInsteadOfMz = value;
    }


    public double getMinDeltaMass()
    {
        return minDeltaMass;
    }

    public void setMinDeltaMass(double minDeltaMass)
    {
        this.minDeltaMass = minDeltaMass;
    }

    public double getMinDeltaElution()
    {
        return minDeltaElution;
    }

    public void setMinDeltaElution(double minDeltaElution)
    {
        this.minDeltaElution = minDeltaElution;
    }

    public void setElutionMode(int newElutionMode)
    {
        super.setElutionMode(newElutionMode);

        switch (newElutionMode)
        {
            case FeatureClusterer.ELUTION_MODE_SCAN:
                elutionBucketIncrement = (float) DEFAULT_SCAN_ELUTION_BUCKET_INCREMENT;
                minDeltaElution = DEFAULT_MIN_DELTA_SCAN;
                break;
            case FeatureClusterer.ELUTION_MODE_TIME:
                elutionBucketIncrement = (float) DEFAULT_TIME_ELUTION_BUCKET_INCREMENT;
                minDeltaElution = DEFAULT_MIN_DELTA_TIME;
                break;
        }
    }
}
