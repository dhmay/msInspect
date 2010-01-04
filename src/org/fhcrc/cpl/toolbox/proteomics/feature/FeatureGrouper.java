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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.normalize.Normalizer;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;

import org.apache.log4j.Logger;

import java.io.PrintWriter;
import java.util.*;

/**
 * Uses clustering (Clusterer2D) to align multiple feature sets
 *
 * dhmay changing 20100104, for peptide array changes.
 * Big changes: Clustering may optionally be done separately per charge state
 */
public class FeatureGrouper
{
    private static Logger _log = Logger.getLogger(FeatureGrouper.class);

    private List<FeatureSet> _featureSets = new ArrayList<FeatureSet>();

    //handles the actual clustering
    protected FeatureClusterer _featureClusterer;
    protected Map<Integer, FeatureClusterer> _chargeClustererMap = new HashMap<Integer, FeatureClusterer>();

    protected boolean _useMassInsteadOfMz = false;

    protected boolean _shouldGroupByCharge = false;
    protected Set<Integer> _allObservedCharges = new HashSet<Integer>();

    protected int _conflictResolver = CONFLICT_RESOLVER_SUM;

    //modes for determining which bucket settings are the best
    public static final int BUCKET_EVALUATION_MODE_ONE_FROM_EACH = 0;
    public static final int BUCKET_EVALUATION_MODE_PEPTIDE_AGREEMENT = 1;

    //These parameters determine the scoring used to decide what the best cluster
    //diameters are, when factoring in peptide agreement

    //minimum number of features in a bucket -- one per set -- to consider
    //the bucket a hit
    protected double bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple = 1;
    //minimum number of matching peptide IDs in the bucket to consider it a hit
    protected double bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple = 0.5;
    //score for a matched bucket
    protected int bucketEvaluationPeptideAgreementMatchScore =
            DEFAULT_PEPTIDE_AGREEMENT_MATCH_SCORE;
    //penalty for a mismatched bucket
    protected int bucketEvaluationPeptideAgreementMismatchPenalty =
            DEFAULT_PEPTIDE_AGREEMENT_MISMATCH_PENALTY;

    public static final int DEFAULT_PEPTIDE_AGREEMENT_MATCH_SCORE=1;
    public static final int DEFAULT_PEPTIDE_AGREEMENT_MISMATCH_PENALTY=1;


    //how to resolve conflicts? I.e., features from the same set in the same
    //cluster
    public static final int CONFLICT_RESOLVER_SUM=0;
    public static final int CONFLICT_RESOLVER_BEST=1;

    //by default, sum the intensities of conflicting features
    public static final int DEFAULT_CONFLICT_RESOLVER=CONFLICT_RESOLVER_SUM;





    /**
     * Initialize the clusterer
     */
    public FeatureGrouper()
    {
        _featureClusterer = createFeatureClusterer();
    }


    protected FeatureClusterer createFeatureClusterer()
    {
        return new FeatureClusterer(_useMassInsteadOfMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                    FeatureClusterer.MASS_MZ_MODE_MZ,
                                     FeatureClusterer.ELUTION_MODE_SCAN);
    }

    /**
     * By default, simply evaluate the bucket settings based on the number of rows with
     * exactly one feature from each bucket
     * @param massOrMzBuckets
     * @param scanBuckets
     * @return
     */
    public Pair<Double, Integer> calculateBestBuckets(double[] massOrMzBuckets,
                                                      int[] scanBuckets)
    {
        return calculateBestBuckets(massOrMzBuckets, scanBuckets,
                                    BUCKET_EVALUATION_MODE_ONE_FROM_EACH);
    }

    public Pair<Double, Integer> calculateBestBuckets(double[] massOrMzBuckets,
                                                      int[] scanBuckets,
                                                      int bucketEvaluationMode)
    {
        Pair<Double, Integer> bestBuckets = null;
        switch (bucketEvaluationMode)
        {
            case BUCKET_EVALUATION_MODE_ONE_FROM_EACH:

                double[] scanBucketsDouble = new double[scanBuckets.length];
                for (int i = 0; i < scanBuckets.length; i++)
                    scanBucketsDouble[i] = scanBuckets[i];
                Pair<Double, Double> bestGenericBuckets =
                        _featureClusterer.calculateBestBuckets(massOrMzBuckets, scanBucketsDouble);

                bestBuckets = new Pair<Double, Integer>(bestGenericBuckets.first,
                                                        bestGenericBuckets.second.intValue());
                break;
            case BUCKET_EVALUATION_MODE_PEPTIDE_AGREEMENT:
                //evaluate bucket breakdowns based on the number of rows with minAligned
                //features, one from each represented set, minPeptideIds of which match
                //to the same peptide ID. Each passing row scores that breakdown matchScore
                //points, and each peptide mismatch loses it mismatchPenalty points
                int numSets = _featureClusterer.getClusterableArrays().size();



                //round up
                int bucketEvaluationPeptideAgreementMinPeptideIds =
                     (int) Math.round(bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple *
                             numSets);
                //round up
                int bucketEvaluationPeptideAgreementMinAligned =
                     (int) Math.round(bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple *
                             numSets);


                double bestMassOrMzBucket = 0;
                int bestScanBucket = 0;

                int[][] agreeingPeptides = new int[massOrMzBuckets.length][scanBuckets.length];
                int bestNumAgreeingPeptides = -1;

                //evaluate all combinations of mass and hydrophobicity buckets
                for (int iMass = 0; iMass < massOrMzBuckets.length; iMass++)
                {
                    for (int iHydrophobicity = 0; iHydrophobicity < scanBuckets.length;
                         iHydrophobicity++)
                    {
                        int scanBucketsize = scanBuckets[iHydrophobicity];
                        double massOrMzBucketsize = massOrMzBuckets[iMass];
                        if (scanBucketsize <= 0 || massOrMzBucketsize <= 0)
                            continue;
                        split2D(massOrMzBucketsize, scanBucketsize);

                        int bucketScore = 0;
                        Clusterer2D.BucketSummary[] bucketSummaries = summarize();



                        for (Clusterer2D.BucketSummary bucketSummary : bucketSummaries)
                        {
                            //see if there's one feature from each represented set in the bucket
                            if (bucketSummary.featureCount >= bucketEvaluationPeptideAgreementMinAligned &&
                                    bucketSummary.setCount == bucketSummary.featureCount)
                            {
                                String commonPeptide = null;
                                boolean failed = false;

                                int numMatched = 0;
                                for (Clusterer2D.Clusterable clusterable : bucketSummary.getParentList())
                                {
                                    Feature parentFeature =
                                            ((FeatureClusterer.FeatureClusterable) clusterable).getParentFeature();
                                    String peptide = MS2ExtraInfoDef.getFirstPeptide(parentFeature);
                                    if (peptide == null)
                                        continue;
                                    if (commonPeptide == null)
                                    {
                                        commonPeptide = peptide;
                                        numMatched++;
                                    }
                                    else
                                    {
                                        if (commonPeptide.equals(peptide))
                                            numMatched++;
                                        else
                                        {
                                            failed=true;
                                            break;
                                        }
                                    }
                                }
                                if (failed)
                                    bucketScore -= bucketEvaluationPeptideAgreementMismatchPenalty;
                                else if (numMatched >= bucketEvaluationPeptideAgreementMinPeptideIds)
                                {
                                    bucketScore += bucketEvaluationPeptideAgreementMatchScore;
                                }
                            }

                        }
                        _log.debug("Evaluating bucket size: mass " + massOrMzBucketsize + ", scan " + scanBucketsize);
                        _log.debug("    Bucket score: " + bucketScore);

                        agreeingPeptides[iMass][iHydrophobicity] = bucketScore;

                        if (agreeingPeptides[iMass][iHydrophobicity] >
                                bestNumAgreeingPeptides)
                        {
                            bestMassOrMzBucket = massOrMzBucketsize;
                            bestScanBucket = scanBucketsize;
                            bestNumAgreeingPeptides = agreeingPeptides[iMass][iHydrophobicity];
                        }
                    }
                }
                ApplicationContext.setMessage("Best number of agreeing peptides: " + bestNumAgreeingPeptides);
                bestBuckets = new Pair<Double, Integer>(bestMassOrMzBucket,
                                                        bestScanBucket);
                break;
        }
        return bestBuckets;
    }


    /**
     * Old school entry point; never does normalization
     */
    public void writePeptideArray(PrintWriter writer)
    {
        writePeptideArray(writer, false, false);
    }

    /**
     * by default, don't sum intensities
     * @param writer
     * @param normalize
     */
    public void writePeptideArray(PrintWriter writer, boolean normalize)
    {
        writePeptideArray(writer, normalize, false);
    }

    public static Feature[] getFeatures(Clusterer2D.BucketSummary bucketSummary)
    {
        List<Feature> resultList = new ArrayList<Feature>();
        for (Clusterer2D.TreeEntry treeEntry : bucketSummary.entries())
            resultList.add(((FeatureClusterer.FeatureClusterable) treeEntry.parent).getParentFeature());

        return resultList.toArray(new Feature[0]);
    }

    public static List<Feature> getFeaturesFromSet(int setIndex,
                                                   Clusterer2D.BucketSummary bucketSummary)
    {
        List<Feature> resultList = new ArrayList<Feature>();
        for (Clusterer2D.Clusterable clusterable : bucketSummary.getParentListForSetIndex(setIndex))
            resultList.add(((FeatureClusterer.FeatureClusterable) clusterable).getParentFeature());

        return resultList;
    }

    /**
     * Write out a peptide array. Optionally normalize.
     * If sumIntensities is true, all intensities from a given FeatureSet in the bucket will
     * be summed.  Otherwise, the intensity included will be from the feature that was picked
     * to represent the bucket
     * @param writer
     * @param normalize
     */
    public void writePeptideArray(PrintWriter writer, boolean normalize,
                                  boolean sumIntensities)
    {
        //should we be writing out peptides and proteins?  Only if any featuresets
        //have those fields
        boolean writePeptidesAndProteins = false;

        //Write peptides and proteins iff any FeatureSet has peptide
        //and protein info
        for (int i = 0; i < _featureSets.size() && (!writePeptidesAndProteins); i++)
        {
            writePeptidesAndProteins =
                 _featureSets.get(i).hasExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());
        }

        List<String> columns = new ArrayList<String>();
        columns.add("id");
        if (this._shouldGroupByCharge)
            columns.add("charge");
        if (_useMassInsteadOfMz)
        {
                    columns.add("minMass");
            columns.add("maxMass");
        }
        else
        {
            columns.add("minMz");
            columns.add("maxMz");
        }
        columns.add("minScan"); columns.add("maxScan");
        columns.add("featureCount"); columns.add("setCount");

        for (int i = 0; i < _featureSets.size(); i++)
        {
            FeatureSet fs = _featureSets.get(i);
            String featureSetColumnSuffix =
                    buildFeatureSetColumnSuffix(fs, i+1);
            columns.add("intensity" + featureSetColumnSuffix);
            if (writePeptidesAndProteins)
            {
                columns.add("peptide" + featureSetColumnSuffix);
                columns.add("protein" + featureSetColumnSuffix);
            }
        }

        for (int i=0; i<columns.size(); i++)
        {
            if (i>0)
                writer.print("\t");
            writer.print(columns.get(i));
        }
        writer.println();





        //TODO: probably this should be a parameter that gets passed in.
        //True, though, it only makes sense if writePeptidesAndProteins is true
        boolean showMultiplePeptideProteinMatches = writePeptidesAndProteins;

        List<FeatureClusterer> clusterersToSummarize = new ArrayList<FeatureClusterer>();
        List<Integer> chargesOfClusterers = new ArrayList<Integer>();
        if (this._shouldGroupByCharge)
        {
            for (int charge : _chargeClustererMap.keySet())
            {
                clusterersToSummarize.add(_chargeClustererMap.get(charge));
                chargesOfClusterers.add(charge);
            }
        }
        else
        {
            clusterersToSummarize.add(_featureClusterer);
        }

        int numTotalRows = 0;
        List<Integer> firstRowIndexesEachSummarizer = new ArrayList<Integer>();


        List<Clusterer2D.BucketSummary[]> summariesAllClusterers = new ArrayList<Clusterer2D.BucketSummary[]>();
        for (int i=0; i<clusterersToSummarize.size(); i++)
        {
            FeatureClusterer clusterer = clusterersToSummarize.get(i);
            firstRowIndexesEachSummarizer.add(numTotalRows);
            if (_shouldGroupByCharge)
                ApplicationContext.infoMessage("\tSummarizing charge " + chargesOfClusterers.get(i) + "...");

            Clusterer2D.BucketSummary[] summaries = clusterer.summarize();
            summariesAllClusterers.add(summaries);
            if (_shouldGroupByCharge)
                ApplicationContext.infoMessage("\t" + summaries.length + " clusters this charge");

            numTotalRows += summaries.length;
        }
        ApplicationContext.infoMessage("Done clustering.");




        Clusterer2D.BucketSummary[] allSummaries = new Clusterer2D.BucketSummary[numTotalRows];
        for (int i=0; i<summariesAllClusterers.size(); i++)
        {
            Clusterer2D.BucketSummary[] summaries = summariesAllClusterers.get(i);
            System.arraycopy(summaries, 0, allSummaries, firstRowIndexesEachSummarizer.get(i), summaries.length);
        }
        summariesAllClusterers = null;
        System.gc();

        ApplicationContext.infoMessage("Calculating intensities for " + allSummaries.length + " cells....");

        //resolve conflicts, determine intensities/features to represent rows

        //stupid java 1.5 can't handle generic array creation
        List<Feature>[][] allSummariesFeatures =
                (List<Feature>[][]) (new List[numTotalRows][_featureSets.size()]);
        float[][] allSummariesIntensities =
                new float[numTotalRows][_featureSets.size()];

        int arrayRowIndex = -1;
        FeatureSet[] featureSetsArray = new FeatureSet[_featureSets.size()];
        for (int i=0; i<_featureSets.size(); i++)
            featureSetsArray[i] = _featureSets.get(i);
        for (Clusterer2D.BucketSummary summary : allSummaries)
        {
            arrayRowIndex++;
            //order features.  If the resolver is currently "sum",
            //order them by "best" anyway; might as well do that so that the
            //peptides and proteins have a known ordering.
            //Order is descending order of goodness, whatever goodness is
            allSummariesFeatures[arrayRowIndex] =
                    orderClusterFeatures(summary, featureSetsArray,
                            _conflictResolver == CONFLICT_RESOLVER_SUM ?
                                    CONFLICT_RESOLVER_BEST : _conflictResolver);

            //Now, assign intensities.
            for (int j = 0; j < _featureSets.size(); j++)
            {
                allSummariesIntensities[arrayRowIndex][j] = 0;
                List<Feature> thisFeatureList = allSummariesFeatures[arrayRowIndex][j];
                if (thisFeatureList != null &&
                        !thisFeatureList.isEmpty())
                {
                    //sum all intensities
                    if (_conflictResolver == CONFLICT_RESOLVER_SUM)
                    {
                        for (Feature feature : thisFeatureList)
                        {
                            allSummariesIntensities[arrayRowIndex][j] += feature.getIntensity();
                        }
                    }
                    else
                    {
                        //pick just one intensity.  take the intensity of the
                        //first one, because they're in descending order of goodness.
                        allSummariesIntensities[arrayRowIndex][j] =
                                thisFeatureList.get(0).getIntensity();
                    }

                }
            }
        }
        ApplicationContext.infoMessage("Done calculating intensities.");


        

        //throw away return value of normalizeSummaryIntensities
        if (normalize)
            normalizeSummaryIntensities(allSummariesIntensities);
        int currentCharge = 0;
        int nextSummarizerFirstRowIndex = 0;
        int currentSummarizerIndex = -1;
        if (_shouldGroupByCharge)
        {
            currentCharge = chargesOfClusterers.get(0);
        }
        ApplicationContext.infoMessage("Writing array...");
        for (int i = 0; i < numTotalRows; i++)
        {
            writer.print(i + 1 + "\t");
            if (_shouldGroupByCharge)
            {
                if (i == nextSummarizerFirstRowIndex)
                {
                    currentSummarizerIndex++;
                    currentCharge = chargesOfClusterers.get(currentSummarizerIndex);
                    ApplicationContext.infoMessage("\tWriting charge " + currentCharge);
                    if (currentSummarizerIndex < chargesOfClusterers.size()-1)
                        nextSummarizerFirstRowIndex = firstRowIndexesEachSummarizer.get(currentSummarizerIndex+1);
                }
                writer.print(currentCharge + "\t");
            }
            writer.println(
                    createArrayRow(allSummaries[i],allSummariesFeatures[i],
                                   allSummariesIntensities[i],
                                   writePeptidesAndProteins,
                                   showMultiplePeptideProteinMatches));
        }
    }

    /**
     * build the column suffix for this featureset.
     *
     * Separating this to make it easy to find and mess with.
     * @param featureSet
     * @return
     */
    protected String buildFeatureSetColumnSuffix(FeatureSet featureSet,
                                                 int featureSetNumber)
    {
        StringBuffer resultBuf = new StringBuffer("_");
        if (featureSet.getSourceFile() == null)
        {
            resultBuf.append(featureSetNumber);
            return resultBuf.toString();
        }
        //add characters from the filename until you run into the first "."
        //turn whitespace into _
        for (byte charByte : featureSet.getSourceFile().getName().getBytes())
        {
            char thisChar = (char) charByte;
            if (charByte == '.')
                break;
            if (Character.isSpaceChar(thisChar))
                thisChar = '_';
            resultBuf.append(thisChar);
        }

        return resultBuf.toString();
    }

    /**
     * Filter each FeatureSet to include only those features that align across
     * at least minAligned sets. Intended for use when runs are grouped by
     * some tag (e.g. "case" or "control")
     */
    public FeatureSet[] filterByGroupedAlignment(int minAligned)
    {
        boolean useTree = false;
        Clusterer2D.BucketSummary[] clusters = _featureClusterer.summarize();
        FeatureSet[] featureSets = new FeatureSet[_featureSets.size()];
        for (int i = 0; i < _featureSets.size(); i++)
        {
            Tree2D tree = null;
            if (useTree)
            {
                tree = new Tree2D();
                for (Feature f :_featureSets.get(i).getFeatures())
                    tree.add(f.getMass(), f.getScan(), f);
            }

            ArrayList<Feature> filteredFeatures = new ArrayList<Feature>();

            for (int j = 0; j < clusters.length; j++)
            {
                Clusterer2D.BucketSummary cluster = clusters[j];
                if (cluster.setCount >= minAligned)
                {
                    float minMass = (float) cluster.minDimension1;
                    float maxMass = (float) cluster.maxDimension1;
                    int minScan = (int) cluster.minDimension2;
                    int maxScan = (int) cluster.maxDimension2;

                    if (useTree)
                    {
                        // use Tree2D to pull out all features within the bucket
                        ArrayList<Feature> bucketFeatures = 
                            (ArrayList<Feature>) tree.getPoints(minMass, minScan, maxMass, maxScan);
                        for (Feature f : bucketFeatures)
                            filteredFeatures.add(f);
                    }
                    else
                    {
                        // Tree2D is not inclusive...
                        for (Feature f :_featureSets.get(i).getFeatures())
                        {
                            if (f.getMass() >= minMass && f.getScan() >= minScan &&
                                f.getMass() <= maxMass && f.getScan() <= maxScan)
                                filteredFeatures.add(f);
                        }
                    }
                }
            }
            Feature[] features = filteredFeatures.toArray(new Feature[0]);
            featureSets[i] = (FeatureSet) _featureSets.get(i).clone();
            featureSets[i].setFeatures(features);
        }
        return featureSets;
    }

    /**
     * Cover method
     * @param cluster
     * @param featureSets
     * @param conflictResolver
     * @return
     */
    public static List<Feature>[]
            orderClusterFeatures(Clusterer2D.BucketSummary cluster,
                                 List<FeatureSet> featureSets,
                                 int conflictResolver)
    {
        return orderClusterFeatures(cluster,
                featureSets.toArray(new FeatureSet[0]), conflictResolver);
    }

    /**
     * Resolve conflicts between features in a cluster by ordering them along
     * some metric
     *
     * Result is a list of Features in descending order of goodness
     * @param cluster
     * @param conflictResolver
     * @return an array of Features, one representing each FeatureSet from
     * featureSets if at least one feature from that FeatureSet was present.
     */
    public static List<Feature>[]
            orderClusterFeatures(Clusterer2D.BucketSummary cluster,
                                 FeatureSet[] featureSetArray,
                                 int conflictResolver)
    {
        //stupid Java can't create arrays of generics
        List<Feature>[] result = (List<Feature>[]) new List[featureSetArray.length];
        switch (conflictResolver)
        {
//This is probably not useful, and difficult to implement right now
//            case CONFLICT_RESOLVER_CLOSEST:
//                Clusterer2D.Clusterable[] oneClusterableFromEachSet =
//                        cluster.pickOneFromEachSet(1);
//
//                for (int j = 0; j < featureSetArray.length; j++)
//                {
//                    if (oneClusterableFromEachSet[j] != null)
//                    {
//                        result[j] =
//                                ((FeatureClusterer.FeatureClusterable)
//                                        oneClusterableFromEachSet[j]).getParentFeature();
//                    }
//
//                }
//                break;
            case CONFLICT_RESOLVER_BEST:
                for (int i = 0; i < featureSetArray.length; i++)
                {
                    List<Clusterer2D.Clusterable> thisSetClusterables =
                            cluster.getParentListForSetIndex(i);
                    if (!thisSetClusterables.isEmpty())
                    {
                        ArrayList<Feature> featureList =
                                new ArrayList<Feature>(thisSetClusterables.size());
                        for (Clusterer2D.Clusterable featureClusterable : thisSetClusterables)
                        {
                            featureList.add(((FeatureClusterer.FeatureClusterable)
                                    featureClusterable).getParentFeature());
                        }
                        Collections.sort(featureList,
                                Feature.PeaksKLComparatorDesc.getSingletonInstance());
                        result[i] = featureList;
                    }
                }
                break;
        }
        return result;
    }

    public Clusterer2D.BucketSummary[] summarize()
    {
        return _featureClusterer.summarize();
    }

    /**
     * Create a single row for a peptide array.
     * @param summary
     * @param bucketFeatures
     * @param bucketIntensities
     * @return
     */
    protected static String createArrayRow(Clusterer2D.BucketSummary summary,
                                           List<Feature>[] bucketFeatures,
                                           float[] bucketIntensities,
                                           boolean writePeptidesAndProteins,
                                           boolean writeMultipleFeatures)
    {
        StringBuffer resultBuf = new StringBuffer(summary.toString());
        for (int i = 0; i < bucketIntensities.length; i++)
        {
            resultBuf.append("\t");
            if (bucketIntensities[i] > 0)
                resultBuf.append(bucketIntensities[i]);

            List<Feature> featureList = bucketFeatures[i];

            //make sure there are actually features to write out.
            boolean hasFeatures = (featureList != null &&
                                   !featureList.isEmpty());

            //if we're writing peptides and proteins, we've got a lot of work to do
            //to determine exactly how to do it
            if (writePeptidesAndProteins)
            {
                resultBuf.append("\t");

                if (writeMultipleFeatures)
                {
                    StringBuffer peptideBuf = new StringBuffer("");
                    StringBuffer proteinBuf = new StringBuffer("");

                    //If no features to write out, buffers remain empty
                    if (hasFeatures)
                    {
                        List<String> peptideListStrings =
                                new ArrayList<String>();
                        List<String> proteinListStrings =
                                new ArrayList<String>();

                        //For each feature, write out its peptide and protein lists,
                        //bracket-separated if necessary.
                        //Do not add duplicate peptide or protein lists
                        for (Feature feature : featureList)
                        {
                            List<String> peptideList =
                                    MS2ExtraInfoDef.getPeptideList(feature);
                            List<String> proteinList =
                                    MS2ExtraInfoDef.getProteinList(feature);
                            //if this peptide list contains duplicate peptides in
                            //sequence, only keep the first of each series.
                            //maintain the protein list in parallel.
                            if (peptideList != null && !peptideList.isEmpty())
                            {
                                MS2ExtraInfoDef.collapsePeptideSequentialDuplicates(feature);
                                String peptideListString =
                                        MS2ExtraInfoDef.convertStringListToString(peptideList);
                                //even if we have a peptide list with something in it,
                                //only add it if we don't already have an identical
                                //list
                                if (!peptideListStrings.contains(peptideListString))
                                {
                                    peptideListStrings.add(peptideListString);
                                    String proteinListString = null;
                                    if (proteinList != null &&
                                        !proteinList.isEmpty())
                                    {
                                        MS2ExtraInfoDef.collapseProteinSequentialDuplicates(
                                                feature);
                                        proteinListString =
                                                MS2ExtraInfoDef.convertStringListToString(
                                                    proteinList);
                                    }
                                    proteinListStrings.add(proteinListString);
                                }
                            }
                        }


                        //if we've got multiple peptides to write, have to
                        //separate them
                        boolean multiplePeptidesToWrite =
                                (peptideListStrings.size() > 1);

                        for (int j=0; j<peptideListStrings.size(); j++)
                        {
                            if (multiplePeptidesToWrite)
                                peptideBuf.append("[");
                            peptideBuf.append(peptideListStrings.get(j));
                            if (multiplePeptidesToWrite)
                                peptideBuf.append("]");
//System.err.println("writing eptides and proteins.  Peptides: " + peptideListStrings.size() + ", proteins: " + proteinListStrings.size());
                            if (proteinListStrings.get(j) != null)
                            {
                                if (multiplePeptidesToWrite)
                                    proteinBuf.append("[");
//System.err.println("   Writing protein " + proteinListStrings.get(j));
                                proteinBuf.append(proteinListStrings.get(j));
                                if (multiplePeptidesToWrite)
                                    proteinBuf.append("]");
                            }
                        }

                    }
                    resultBuf.append(peptideBuf);
                    resultBuf.append("\t");
                    resultBuf.append(proteinBuf);
                }
                else
                {
                    //not writing multiple features, just take the first
                    if (hasFeatures &&
                        MS2ExtraInfoDef.getFirstPeptide(featureList.get(0)) != null)
                        resultBuf.append(MS2ExtraInfoDef.getFirstPeptide(featureList.get(0)));
                    resultBuf.append("\t");
                    if (hasFeatures &&
                        MS2ExtraInfoDef.getFirstProtein(featureList.get(0)) != null)
                        resultBuf.append(MS2ExtraInfoDef.getFirstProtein(featureList.get(0)));
                }
            }
        }
        return resultBuf.toString();
    }


    /**
     * Write out all feature details (except peptide, protein)
     * @param writer
     * If true, writes the original, un-deconvoluted features
     */
    public void writeArrayDetails(PrintWriter writer, boolean writeMS2)
    {
        writer.print("id\t");
        writer.print("file\t");
        FeatureExtraInformationDef[] extraInfos = writeMS2 ?
                new FeatureExtraInformationDef[] {MS2ExtraInfoDef.getSingletonInstance()} : null;
        writer.println(Feature.getFeatureHeader(extraInfos));

        List<FeatureClusterer> clusterersToSummarize = new ArrayList<FeatureClusterer>();
        if (this._shouldGroupByCharge)
        {
            for (int charge : _chargeClustererMap.keySet())
            {
                clusterersToSummarize.add(_chargeClustererMap.get(charge));
            }
        }
        else
        {
            clusterersToSummarize.add(_featureClusterer);
        }

        int curRow = 0;
        for (FeatureClusterer clusterer : clusterersToSummarize)
        {
            Clusterer2D.BucketSummary[] summaries = clusterer.summarize();
            for (int i = 0; i < summaries.length; i++)
            {
                curRow++;
                Clusterer2D.TreeEntry[] entries= summaries[i].entries();
                for (int j = 0; j < entries.length; j++)
                {
                    String setName = getSetName(entries[j].iSet);
                    Feature deconvolutedFeature =
                            ((FeatureClusterer.FeatureClusterable) (entries[j].parent)).getParentFeature();
                    writeDetailsFeatureRow(writer, curRow, setName, deconvolutedFeature, extraInfos);
                }
                writer.flush();
            }
        }
    }

    /**
     * Helper method to write a single row in a details file representing one feature
     * @param writer
     * @param rowId
     * @param setName
     * @param feature
     */
    protected void writeDetailsFeatureRow(PrintWriter writer, int rowId, String setName, Feature feature,
                                          FeatureExtraInformationDef[] extraInfos)
    {
        writer.println((rowId) + "\t" + setName + "\t" + feature.toString(extraInfos));
    }

    /**
     * get the name of a FeatureSet, for the details file
     * @param i
     * @return
     */
    private String getSetName(int i)
    {
        FeatureSet fs = _featureSets.get(i);
        if (fs.getSourceFile() != null)
            return fs.getSourceFile().getName();
        else
            return String.valueOf(i);
    }


    /**
     * Normalize the intensities across all bucket summaries. Note that we don't
     * attempt to reflect normalization in the feature "details" file.
     */
    private static boolean normalizeSummaryIntensities(float[][] allSummariesIntensities)
    {
        List<float[]> rows = new ArrayList<float[]>();
        for (int i = 0; i < allSummariesIntensities.length; i++)
            rows.add(allSummariesIntensities[i]);
        return Normalizer.normalize(rows);
    }


    //cover methods on FeatureClusterer
    /**
     * Split based on a set of bucket sizes.  If not splitting by charge, just passes the command along to the clusterer.
     * If splitting by charge, splits each of the per-charge clusterers
     * @param dimension1Bucket
     * @param dimension2Bucket
     */
    public void split2D(double dimension1Bucket, double dimension2Bucket)
    {
        if (this._shouldGroupByCharge)
        {
            for (FeatureClusterer featureClustererThisCharge : _chargeClustererMap.values())
            {
                featureClustererThisCharge.split2D(dimension1Bucket, dimension2Bucket);
            }
        }
        else
        {
            _featureClusterer.split2D(dimension1Bucket, dimension2Bucket);
        }
    }

    public int numBuckets()
    {
        return _featureClusterer.numBuckets();
    }

    public int rowsWithOneFromEach()
    {
        if (this._shouldGroupByCharge)
        {
            int oneFromEachRowCount = 0;
            for (int charge : _allObservedCharges)
            {
                oneFromEachRowCount += _chargeClustererMap.get(charge).rowsWithOneFromEach();
            }
            return oneFromEachRowCount;
        }
        else
            return _featureClusterer.rowsWithOneFromEach();
    }




    //getters, setters

    public boolean getGroupByMass()
    {
        return _useMassInsteadOfMz;
    }

    public void setGroupByMass(boolean groupByMass)
    {
        this._useMassInsteadOfMz = groupByMass;
        if (_shouldGroupByCharge)
        {
            for (FeatureClusterer fc : _chargeClustererMap.values())
                fc.setMassMzMode(_useMassInsteadOfMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                              FeatureClusterer.MASS_MZ_MODE_MZ);
        }
        else
        {
        _featureClusterer.setMassMzMode(_useMassInsteadOfMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                              FeatureClusterer.MASS_MZ_MODE_MZ);
        }
    }

    /**
     * Add a featureset.  If we're not separating by charge, this is trivial.  If we are, then we need to
     * separate out features by charge, then create a featureset for each charge, then add each of those to
     * the featureclusterer for that charge
     * @param fs
     */
    public void addSet(FeatureSet fs)
    {
        if (this._shouldGroupByCharge)
        {

            Map<Integer,List<Feature>> chargeFeatureMap = new HashMap<Integer,List<Feature>>();
            {
                for (Feature feature : fs.getFeatures())
                {
                    if (feature.getCharge() == 0)
                        continue;
                    List<Feature> featuresThisCharge = chargeFeatureMap.get(feature.getCharge());
                    if (featuresThisCharge == null)
                    {
                        featuresThisCharge = new ArrayList<Feature>();
                        chargeFeatureMap.put(feature.getCharge(), featuresThisCharge);
                        _allObservedCharges.add(feature.getCharge());
                    }
                    featuresThisCharge.add(feature);
                }
            }
            for (int charge=1; charge<=6; charge++)
            {
                FeatureSet fsThisCharge = (FeatureSet) fs.clone();
                if (chargeFeatureMap.containsKey(charge))
                    fsThisCharge.setFeatures(chargeFeatureMap.get(charge).toArray(new Feature[chargeFeatureMap.get(charge).size()]));
                else fsThisCharge.setFeatures(new Feature[0]);

                FeatureClusterer featureClustererThisCharge = _chargeClustererMap.get(charge);
                if (featureClustererThisCharge == null)
                {
                    featureClustererThisCharge = createFeatureClusterer();
                    _chargeClustererMap.put(charge, featureClustererThisCharge);
                }
                featureClustererThisCharge.addSet(fsThisCharge);
            }
        }
        else
        {
            _featureClusterer.addSet(fs);
        }
        _featureSets.add(fs);
    }

    public FeatureSet getSet(int i)
    {
        return _featureClusterer.getSet(i);
    }

    public int getConflictResolver()
    {
        return _conflictResolver;
    }

    public void setConflictResolver(int conflictResolver)
    {
        this._conflictResolver = conflictResolver;
    }

    public double getBucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple()
    {
        return bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple;
    }

    public void setBucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple(double bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple)
    {
        this.bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple = bucketEvaluationPeptideAgreementMinAlignedNumSetsMultiple;
    }

    public double getBucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple()
    {
        return bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple;
    }

    public void setBucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple(double bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple)
    {
        this.bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple = bucketEvaluationPeptideAgreementMinPeptideIdsNumSetsMultiple;
    }

    public int getBucketEvaluationPeptideAgreementMatchScore()
    {
        return bucketEvaluationPeptideAgreementMatchScore;
    }

    public void setBucketEvaluationPeptideAgreementMatchScore(int bucketEvaluationPeptideAgreementMatchScore)
    {
        this.bucketEvaluationPeptideAgreementMatchScore = bucketEvaluationPeptideAgreementMatchScore;
    }

    public int getBucketEvaluationPeptideAgreementMismatchPenalty()
    {
        return bucketEvaluationPeptideAgreementMismatchPenalty;
    }

    public void setBucketEvaluationPeptideAgreementMismatchPenalty(int bucketEvaluationPeptideAgreementMismatchPenalty)
    {
        this.bucketEvaluationPeptideAgreementMismatchPenalty = bucketEvaluationPeptideAgreementMismatchPenalty;
    }

    public int getMassType()
    {
        return _featureClusterer.getMassType();
    }

    public void setMassType(int _massType)
    {
        _featureClusterer.setMassType(_massType);
    }

    public int getElutionMode()
    {
        return _featureClusterer.getElutionMode();
    }

    public void setElutionMode(int _elutionMode)       
    {
        _featureClusterer.setElutionMode(_elutionMode);
    }

    public boolean isGroupByCharge()
    {
        return _shouldGroupByCharge;
    }

    public void setGroupByCharge(boolean _shouldGroupByCharge)
    {
        this._shouldGroupByCharge = _shouldGroupByCharge;
    }
}
