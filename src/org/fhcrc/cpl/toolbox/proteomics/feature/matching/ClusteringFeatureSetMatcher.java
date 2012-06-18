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

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.feature.*;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Attempts to match one feature set against another using clustering
 */
public class ClusteringFeatureSetMatcher extends BaseFeatureSetMatcherImpl
        implements FeatureSetMatcher
{
    private static Logger _log = Logger.getLogger(ClusteringFeatureSetMatcher.class);

    public static final double DEFAULT_HYDRO_ELUTION_BUCKET_INCREMENT = 0.05;

    protected static final double ABSOLUTE_MASS_BUCKET_INCREMENT = 0.05;
    protected static final double PPM_MASS_BUCKET_INCREMENT = 1;
    protected double massBucketIncrement = ABSOLUTE_MASS_BUCKET_INCREMENT;
    protected double elutionBucketIncrement = DEFAULT_HYDRO_ELUTION_BUCKET_INCREMENT;

    //Up these values only with caution!  The complexity of this feature-matching grows with
    //numMassBuckets * numElutionBuckets
    protected int numMassBuckets = 4;
    protected int numElutionBuckets = 4;

    protected double bestMassBucketSize = 0;
    protected double bestElutionBucketSize = 0;

   protected boolean useMassInsteadOfMz = true;



    public ClusteringFeatureSetMatcher()
    {
    }

    public ClusteringFeatureSetMatcher(float deltaMass, int deltaMassType,
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
        return "ClusteringFeatureSetMatcher: deltaMass=" + deltaMass +
                ", deltaElution=" + deltaElution;
    }

    /**
     * Calculate the buckets to use for either mass or hydro clustering
     * @param numBuckets
     * @param bucketIncrement
     * @param maxBucketSize
     * @return
     */
    protected double[] calculateBuckets(int numBuckets, double bucketIncrement,
                                        double maxBucketSize)
    {
        _log.debug("Buckets:");
        ArrayList<Double> resultList = new ArrayList<Double>();
        double lowBucket = maxBucketSize - (bucketIncrement * (numBuckets - 1));

        for (int i=0; i<numBuckets; i++)
        {
            double bucketSize = lowBucket + (i * bucketIncrement);
            if (bucketSize > 0)
            {
                resultList.add(bucketSize);
                _log.debug("  " + resultList.get(resultList.size()-1));
            }
        }
        double[] result = new double[resultList.size()];
        for (int i=0; i<resultList.size(); i++)
            result[i] = resultList.get(i);
        return result;
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
//for (Feature feature : featureSet1.getFeatures())
//{
//    System.err.println(AmtExtraInfoDef.getObservedHydrophobicity(feature));
//}
        FeatureClusterer clusterer =
                new FeatureClusterer(useMassInsteadOfMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                          FeatureClusterer.MASS_MZ_MODE_MZ,
                                     elutionMode, featureSet1);
        clusterer.addSet(featureSet2);
//System.err.println("matcF 1, ms1Features: " + featureSet1.getFeatures().length + ", ms2Features: " + featureSet2.getFeatures().length);
//System.err.println("mass? " + useMassInsteadOfMz + ", elut: " + elutionMode);
//System.err.println("deltaMass: " + deltaMass + ", deltaElut: " + deltaElution);

//for (Feature feature : featureSet1.getFeatures()) System.err.println(" H: " + AmtExtraInfoDef.getObservedHydrophobicity(feature));
//System.err.println("set 1 features: " + featureSet1.getFeatures().length);

//      In the PPM case, we have to convert from PPM to Da before figuring out whether we can split.
//      This is done by converting maxBucket based on this bucket's getMin() value.  This will
//      ensure that a bucket will be split if any subregion of the bucket would be split for the
//      given PPM value.
//      This should be reasonable behavior, perhaps splitting slightly more than necessary.
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

        FeatureMatchingResult result = new FeatureMatchingResult();

        //maximum number of perfect matches we could possibly have
//        int maxPossiblePerfectMatches = Math.min(featureSet1.getFeatures().length,
//                                                 featureSet2.getFeatures().length);



        double[] massBuckets = calculateBuckets(numMassBuckets, massBucketIncrement, getDeltaMass());
        double[] elutionBuckets = calculateBuckets(numElutionBuckets,
                                                          elutionBucketIncrement,
                                                          getDeltaElution());

        setBestElutionBucketSize(0);
        setBestMassBucketSize(0);

        Pair<Double,Double> bestBucketSizes = clusterer.calculateBestBuckets(massBuckets, elutionBuckets, false);
        clusterer.split2D(bestBucketSizes.first, bestBucketSizes.second);
        _log.debug("Splitting on buckets " + bestBucketSizes.first + " and " +
                   bestBucketSizes.second);
        Clusterer2D.BucketSummary[] bucketSummaries = clusterer.summarize();

int numConflictBuckets = 0;
        for (Clusterer2D.BucketSummary bucketSummary : bucketSummaries)
        {
            //for each bucket, if it has at least one pair, do something with it
            if (bucketSummary.setCount == 2)
            {
                if (bucketSummary.entries().length > 2)
                {
                    List<Feature> slaveSetFeatures =
                            FeatureGrouper.getFeaturesFromSet(1, bucketSummary);
                    for (Feature masterSetFeature :
                            FeatureGrouper.getFeaturesFromSet(0, bucketSummary))
                    {
                        //order slave set features in the order determined by
                        //variable settings at the superclass level
                        List<Feature> orderedSlaveSetFeatures =
                                orderSlaveSetFeatures(slaveSetFeatures, masterSetFeature);
                        result.put(masterSetFeature, orderedSlaveSetFeatures);
                    }
                    numConflictBuckets++;
//                    Feature[] oneFromEachSet =
//                            FeatureGrouper.resolveClusterConflicts(bucketSummary,
//                                    new FeatureSet[] {featureSet1, featureSet2},
//                                    FeatureGrouper.CONFLICT_RESOLVER_BEST);
                }
                else
                {
                    //exactly one from each, no-brainer
                    Feature masterSetFeature = ((FeatureClusterer.FeatureClusterable)
                            bucketSummary.getParentListForSetIndex(0).get(0)).getParentFeature();
                    Feature slaveSetFeature = ((FeatureClusterer.FeatureClusterable)
                            bucketSummary.getParentListForSetIndex(1).get(0)).getParentFeature();
                    result.add(masterSetFeature, slaveSetFeature);
//System.err.println("first scan: " + firstOne.getParentFeature().scan + ", second scan: " + secondOne.getParentFeature().scan);                    
                }
/*
                List<Clusterer2D.Clusterable> set2Clusterables = bucketSummary.getParentListForSetIndex(1);
                List<Clusterer2D.Clusterable> set1Clusterables = bucketSummary.getParentListForSetIndex(0);
                List<Feature> set2Features = getFeatureListForClusterableList(set2Clusterables);
                List<Feature> set1Features = getFeatureListForClusterableList(set1Clusterables);
                //if there are only two features in the bucket, then obviously we've got one from each
                if (bucketSummary.featureCount == 2)
                {
                    result.add(new Pair<Feature, Feature>(set2Features.get(0),
                                                          set1Features.get(0)));
                }
                else if (bucketSummary.featureCount > 2)
                {
                    numConflictBuckets++;
                    result.add(resolvePairingConflicts(set2Features, set1Features));
                }
*/
            }
        }
//System.err.println("Perfect buckets: " + clusterer.rowsWithOneFromEach(bucketSummaries) + ", conflicts: " + numConflictBuckets);

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



    /**
     * Pick the closest feature to a given feature by mass, from a list.
     * @param baseFeature
     * @param featuresToMatch
     * @return
     */
    protected Feature pickClosestByMass(Feature baseFeature, List<Feature> featuresToMatch)
    {
        Feature result = null;
        //worst value possible
        double bestMassDistance = Double.MAX_VALUE;
        for (Feature featureToMatch : featuresToMatch)
        {
            double currentMassDistance = Math.abs(baseFeature.getMass() - featureToMatch.getMass());
            if (currentMassDistance < bestMassDistance)
            {
                bestMassDistance = currentMassDistance;
                result = featureToMatch;
            }
        }
        return result;
    }


    public double getBestMassBucketSize()
    {
        return bestMassBucketSize;
    }

    public void setBestMassBucketSize(double bestMassBucketSize)
    {
        this.bestMassBucketSize = bestMassBucketSize;
    }

    public double getBestElutionBucketSize()
    {
        return bestElutionBucketSize;
    }

    public void setBestElutionBucketSize(double bestElutionBucketSize)
    {
        this.bestElutionBucketSize = bestElutionBucketSize;
    }


    public int getNumMassBuckets()
    {
        return numMassBuckets;
    }

    public void setNumMassBuckets(int numMassBuckets)
    {
        this.numMassBuckets = numMassBuckets;
    }

    public int getNumElutionBuckets()
    {
        return numElutionBuckets;
    }

    public void setNumElutionBuckets(int numElutionBuckets)
    {
        this.numElutionBuckets = numElutionBuckets;
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
}
