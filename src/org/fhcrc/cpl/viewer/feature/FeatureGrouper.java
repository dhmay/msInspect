/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.feature;

import org.labkey.common.util.Pair;
import org.labkey.common.tools.Tree2D;
import org.labkey.common.tools.ApplicationContext;

import org.fhcrc.cpl.viewer.util.Clusterer2D;
import org.fhcrc.cpl.viewer.normalize.Normalizer;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.apache.log4j.Logger;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Uses clustering (Clusterer2D) to align multiple feature sets
 */
public class FeatureGrouper
{
    private static Logger _log = Logger.getLogger(FeatureGrouper.class);

    private List<FeatureSet> _featureSets = new ArrayList<FeatureSet>();

    //handles the actual clustering
    protected FeatureClusterer _featureClusterer;

    protected boolean _useMassInsteadOfMz = false;

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
        _featureClusterer = new FeatureClusterer();
        _featureClusterer.setMassMzMode(FeatureClusterer.MASS_MZ_MODE_MZ);
        _featureClusterer.setElutionMode(FeatureClusterer.ELUTION_MODE_SCAN);
    }

    /**
     * Initialize, set usemassformz, add a featureset
     * @param fs
     * @param useMassForMz
     */
    public FeatureGrouper(FeatureSet fs, boolean useMassForMz)
    {
        _useMassInsteadOfMz = useMassForMz;
        _featureClusterer =
                new FeatureClusterer(useMassForMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                    FeatureClusterer.MASS_MZ_MODE_MZ,
                                     FeatureClusterer.ELUTION_MODE_SCAN, fs);
        this.setGroupByMass(useMassForMz);
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


        //build the header row
        if (_useMassInsteadOfMz)
            writer.print("id\tminMass\tmaxMass\tminScan\tmaxScan\tfeatureCount\tsetCount");
        else
            writer.print("id\tminMz\tmaxMz\tminScan\tmaxScan\tfeatureCount\tsetCount");

        //dhmay changing -- we used to use column names like intensity1, intensity2, etc.
        //Now we try to use the source filename
        for (int i = 0; i < _featureSets.size(); i++)
        {
            FeatureSet fs = _featureSets.get(i);
            String featureSetColumnSuffix =
                    buildFeatureSetColumnSuffix(_featureSets.get(i), i+1);
            writer.print("\tintensity" + featureSetColumnSuffix);
            if (writePeptidesAndProteins)
                writer.print("\tpeptide" + featureSetColumnSuffix +
                             "\tprotein" + featureSetColumnSuffix);
        }
        writer.println();

        //build the row summaries
        Clusterer2D.BucketSummary[] summaries = _featureClusterer.summarize();

        //resolve conflicts, determine intensities/features to represent rows


        //TODO: probably this should be a parameter that gets passed in.
        //True, though, it only makes sense if writePeptidesAndProteins is true
        boolean showMultiplePeptideProteinMatches = writePeptidesAndProteins;

        //stupid java 1.5 can't handle generic array creation
        List<Feature>[][] allSummariesFeatures =
               (List<Feature>[][]) (new List[summaries.length][_featureSets.size()]);
        float[][] allSummariesIntensities =
                new float[summaries.length][_featureSets.size()];

        for (int i = 0; i < summaries.length; i++)
        {
            Clusterer2D.BucketSummary summary = summaries[i];

            //order features.  If the resolver is currently "sum",
            //order them by "best" anyway; might as well do that so that the
            //peptides and proteins have a known ordering.
            //Order is descending order of goodness, whatever goodness is
            allSummariesFeatures[i] =
                    orderClusterFeatures(summary, _featureSets,
                                         _conflictResolver == CONFLICT_RESOLVER_SUM ?
                                         CONFLICT_RESOLVER_BEST : _conflictResolver);

            //Now, assign intensities.
            for (int j = 0; j < _featureSets.size(); j++)
            {
                allSummariesIntensities[i][j] = 0;
                List<Feature> thisFeatureList = allSummariesFeatures[i][j];
                if (thisFeatureList != null &&
                        !thisFeatureList.isEmpty())
                {
                    //sum all intensities
                    if (_conflictResolver == CONFLICT_RESOLVER_SUM)
                    {
                        for (Feature feature : thisFeatureList)
                        {
                            allSummariesIntensities[i][j] += feature.getIntensity();
                        }
                    }
                    else
                    {
                        //pick just one intensity.  take the intensity of the
                        //first one, because they're in descending order of goodness.
                        allSummariesIntensities[i][j] =
                                thisFeatureList.get(0).getIntensity();
                    }

                }
            }
        }

        //throw away return value of normalizeSummaryIntensities
        if (normalize)
            normalizeSummaryIntensities(allSummariesIntensities);

        for (int i = 0; i < summaries.length; i++)
        {
            writer.print(i + 1 + "\t");
            writer.println(
                    createArrayRow(summaries[i],allSummariesFeatures[i],
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


/*
//a product of another time
                            MS2ExtraInfoDef.PeptideOrProteinList peptideList =
                                MS2ExtraInfoDef.getPeptideList(feature);
                            if (peptideList != null && !peptideList.isEmpty())
                            {
                                if (hasMultipleFeaturesWithPeptides)
                                    peptideBuf.append("[");
                                peptideBuf.append(peptideList);
                                if (hasMultipleFeaturesWithPeptides)
                                    peptideBuf.append("]");
                            }
                            MS2ExtraInfoDef.PeptideOrProteinList proteinList =
                                MS2ExtraInfoDef.getProteinList(feature);
                            if (proteinList != null && !proteinList.isEmpty())
                            {
                                if (hasMultipleFeaturesWithPeptides)
                                    proteinBuf.append("[");
                                proteinBuf.append(proteinList);
                                if (hasMultipleFeaturesWithPeptides)
                                    proteinBuf.append("]");
                            }
*/
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
     */
    public void writeArrayDetails(PrintWriter writer)
    {
        writer.print("id\t");
        writer.print("file\t");
        writer.println(Feature.getFeatureHeader(null));
        Clusterer2D.BucketSummary[] summaries = _featureClusterer.summarize();
        for (int i = 0; i < summaries.length; i++)
        {
            Clusterer2D.TreeEntry[] entries= summaries[i].entries();
            for (int j = 0; j < entries.length; j++)
            {
                writer.print((i + 1) + "\t");
                writer.print(getSetName(entries[j].iSet) + "\t");
                writer.println(((FeatureClusterer.FeatureClusterable) (entries[j].parent)).getParentFeature().toString(null));
            }
            writer.flush();
        }
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
     * Split based on a set of bucket sizes.  Just passes the command along to the clusterer
     * @param dimension1Bucket
     * @param dimension2Bucket
     */
    public void split2D(double dimension1Bucket, double dimension2Bucket)
    {
        _featureClusterer.split2D(dimension1Bucket, dimension2Bucket);
    }

    public int numBuckets()
    {
        return _featureClusterer.numBuckets();
    }

    public int rowsWithOneFromEach()
    {
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
        _featureClusterer.setMassMzMode(_useMassInsteadOfMz ? FeatureClusterer.MASS_MZ_MODE_MASS :
                                                              FeatureClusterer.MASS_MZ_MODE_MZ);
    }

    public void addSet(FeatureSet fs)
    {
        _featureClusterer.addSet(fs);
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

}
