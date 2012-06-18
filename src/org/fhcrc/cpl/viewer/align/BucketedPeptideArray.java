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
package org.fhcrc.cpl.viewer.align;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.viewer.amt.AmtMatchProbabilityAssigner;

import java.awt.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * User: migra
 * Date: Mar 9, 2005
 * Time: 11:06:52 AM
 *
 * dhmay changing 20100104.  Big changes: separating features by charge before clustering.  Alignment is still done based 
 * on a deconvoluted featureset, and normalization is done on a deconvoluted featureset in which all feature charges
 * are assigned 1.  But clustering is done separately for each charge state.  That way, feature intensity comparison
 * actually means something.
 *
 * dhmay changing 20100311.  More big changes.  Two kinds of optimization now:
 * 1.  "perfectbuckets".  This is the old-school optimization in which we try a bunch of combinations of parameters
 * and pick the one that gives us the most 'perfect buckets', i.e., rows with one feature from each run.  This
 * generally works pretty well, but the downside is that it can only evaluate the bucket sizes you give it.
 * 2.  "em".  This is a mixed-model approach, pretty much identical to what we do in calculating AMT match probability.
 * We perform a loose clustering and then calculate the distances, in both dimensions, between all pairs of features
 * that are clustered together (up to some maximum).  Then we assume that this 2D error distribution is a mixed
 * distribution of true matches (bivariate normal) and false (uniform), and calculate the parameters of the distributions
 * accordingly.  Cluster dimensions are then figured out based on the highest possible error for a match with probability
 * greater than some threashold of being in the true distribution.  The bounding box set by these values (note: not
 * equivalent to the area actually containing matches with those probabilities, which is not rectangular) is used
 * for clustering.  This isn't as good at maximizing "perfect buckets" as the "perfectbuckets" strategy is if you give
 * "perfectbuckets" a good set of tolerances, but it can be better if you're not sure about your mass accuracy / RT
 * reproducibility.
 */
public class BucketedPeptideArray implements Runnable
{
    private static Logger _log = Logger.getLogger(BucketedPeptideArray.class);

    private java.util.List _sets;
    private FeatureSet.FeatureSelector _sel;
    boolean _align = true;
    boolean _normalize = false;

    double _elutionBucket = 50;
    double _massBucket = .2;

    protected Aligner _aligner = null;

    //different bucket sizes used for optimization
    protected double[] _elutionBuckets = {20, 30, 50, 75, 100, 150, 200, 300, 400};
    public static final double[] DEFAULT_MASS_BUCKETS_DA = {.025, .05, .1, .15, .2};
    public static final double[] DEFAULT_MASS_BUCKETS_PPM = {3, 5, 7, 10, 20};
    protected double[] _massBuckets = DEFAULT_MASS_BUCKETS_DA;

    protected FeatureGrouper _featureGrouper;

    protected Aligner.FeaturePairSelector _featurePairSelector =
            Aligner.DEFAULT_FEATURE_PAIR_SELECTOR;

    private String _outFileName;

    protected int _conflictResolver = FeatureGrouper.DEFAULT_CONFLICT_RESOLVER;

    //parameters related to (optional) deconvolution
    public static final int DEFAULT_DECONVOLUTE_SCAN_WINDOW = 6;
    public static final double DEFAULT_DECONVOLUTE_MASS_WINDOW = 0.2;
    protected int _deconvoluteScanWindow = DEFAULT_DECONVOLUTE_SCAN_WINDOW;
    protected double _deconvoluteMassWindow = DEFAULT_DECONVOLUTE_MASS_WINDOW;
    protected boolean _shouldDeconvolute = false;


    //should we optimize based on the distribution of mass and rt error?
    public static final int OPTIMIZE_MODE_ERRORDIST = 0;
    public static final int OPTIMIZE_MODE_PERFECTBUCKETS = 1;

    public static final int DEFAULT_OPTIMIZATION_MODE = OPTIMIZE_MODE_PERFECTBUCKETS;

    protected int optimizationMode = DEFAULT_OPTIMIZATION_MODE;


    //For calculating the initial, wide tolerances for generating the dataset that's fed to the EM algorithm
    //during distribution optimization
    public static final float DEFAULT_MASS_WIDE_TOLERANCE_PPM = 40;
    public static final float DEFAULT_MASS_WIDE_TOLERANCE_DA = 0.2f;
    public static final float DEFAULT_RT_WIDE_TOLERANCE = 1000;

    protected float sdMultipleForWideToleranceCalc = 4f;
    protected float optimizationRTWideTolerance = DEFAULT_RT_WIDE_TOLERANCE;
    protected float optimizationMassWideTolerance = DEFAULT_MASS_WIDE_TOLERANCE_PPM;


    //The EM algorithm performance bogs down hugely if you give it too many datapoints.  Max datapoints to use.
    protected int maxDatapointsForEMAlgorithm = 10000;

    //In calculating a bounding box for optimized tolerances, minimum match probability for locating features
    //to be included in the box
    public static final float DEFAULT_MAX_MATCH_FDR_EM_OPT = 0.05f;
    protected float maxMatchFDRForToleranceBoxCalc = DEFAULT_MAX_MATCH_FDR_EM_OPT;
    //Initial assumption about the proportion of matches (with wide tolerances) that are true matches.  This is
    //likely to be an extremely high proportion if there's any real signal here.  If not, the EM algorithm will
    //correct this
    protected float initialDistAssumptionTrue = 0.95f;

    public BucketedPeptideArray(List<?> sets)
    {
        _sets = sets;
        _featureGrouper = new FeatureGrouper();
        _aligner = new SplineAligner();
    }

    public BucketedPeptideArray(List<?> sets, FeatureSet.FeatureSelector sel)
    {
        this(sets);
        _sel = sel;
    }

    public BucketedPeptideArray(List<?> sets, FeatureSet.FeatureSelector sel, double elutionBucket, double massBucket,
                                String outFileName, boolean align)
    {
        this(sets,sel);
        _elutionBucket = elutionBucket;
        _massBucket = massBucket;
        _outFileName = outFileName;
        _align = align;
    }

    public void run()
    {
        run(false);
    }

    public void run(boolean optimize)
    {
        run(optimize, FeatureGrouper.BUCKET_EVALUATION_MODE_ONE_FROM_EACH);
    }

    public void run(boolean optimize, int optimizationMode)
    {
        run(optimize, optimizationMode, false);
    }

    public void run(boolean optimize, int optimizationPerfectBucketMode, boolean showCharts)
    {
        String detailsFileName = null;
//System.err.println(_massBucket); System.err.println(_featureGrouper.getMassType());
        PrintWriter out = null;
        if (_sel == null)
        {
            _log.debug("*************CREATING DUMMY SELECTOR");
            _sel = new FeatureSet.FeatureSelector();
        }

        try
        {
            ApplicationContext.setMessage("Loading");
            java.util.List<FeatureSet> sourceFeatureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < _sets.size(); i++)
            {
                FeatureSet fs;
                if (_sets.get(i) instanceof File)
                    fs = new FeatureSet((File) _sets.get(i), Color.RED);
                else if (_sets.get(i) instanceof FeatureSet)
                    fs = (FeatureSet) _sets.get(i);
                else if (_sets.get(i) instanceof String)
                    fs = new FeatureSet(new File((String) _sets.get(i)), Color.RED);
                else
                {
                    ApplicationContext.errorMessage("Couldn't load feature set due to bad object: " + _sets.get(i).toString(), null);
                    return;
                }

                sourceFeatureSets.add(fs);
            }

            //Filter
            ApplicationContext.setMessage("Filtering");
            _log.debug("sel: " + _sel);
            List<FeatureSet> featureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < sourceFeatureSets.size(); i++)
            {
                FeatureSet fs = (sourceFeatureSets.get(i)).filter(_sel);
                _log.debug("\tbefore: " + sourceFeatureSets.get(i).getFeatures().length + ", after: " +
                        fs.getFeatures().length);                
                featureSets.add(fs);
            }

            //Align
            if (_align)
            {
                ApplicationContext.setMessage("Aligning");
                _aligner.setFeaturePairSelector(_featurePairSelector);
                featureSets = _aligner.alignFeatureSets(featureSets, showCharts);
            }

            //Deconvolute
            if (_shouldDeconvolute)
            {
                ApplicationContext.setMessage("Deconvoluting");
                for (int i=0; i<featureSets.size(); i++)
                {
                    FeatureSet origSet = featureSets.get(i);
                    FeatureSet deconvolutedSet = origSet.deconvolute(_deconvoluteScanWindow, _deconvoluteMassWindow,
                                                                     true);
                    ApplicationContext.setMessage("\tCollapsed " + origSet.getFeatures().length + " features into " +
                            deconvolutedSet.getFeatures().length);
                    featureSets.set(i, deconvolutedSet);
                }
            }

            _featureGrouper.setShowCharts(showCharts);
            _featureGrouper.setGroupByMass(true);
            _featureGrouper.setGroupByCharge(true);
            _featureGrouper.setConflictResolver(_conflictResolver);

            if (optimize)
            {
                switch (optimizationMode)
                {
                    case OPTIMIZE_MODE_PERFECTBUCKETS:
                        optimizePerfectBuckets(featureSets, optimizationPerfectBucketMode, showCharts);
                        break;
                    case OPTIMIZE_MODE_ERRORDIST:
                        optimizeDist(featureSets, showCharts);
                        break;
                }
            }


            for (FeatureSet fs : featureSets)
            {
                _featureGrouper.addSet(fs);
            }

            _featureGrouper.split2D(_massBucket, _elutionBucket);
            ApplicationContext.setMessage("Perfect buckets: " + _featureGrouper.rowsWithOneFromEach());

            if (_outFileName != null)
                out = new PrintWriter(new FileOutputStream(_outFileName));
            else
                out = new PrintWriter(System.out);

            writeHeader(out);
            float[][] intensitiesAfterProcessing =
                    _featureGrouper.writePeptideArray(out, _normalize); // ???? check for successful write?

            out.flush();

            //Write details file if outfilename != null
            if (_outFileName != null)
            {
                out.close();
                detailsFileName = calcDetailsFilepath(_outFileName);
                out = new PrintWriter(new FileOutputStream(detailsFileName));
                writeHeader(out);
                //todo: figure out dynamically if we need to write out MS2 extrainfo
                ApplicationContext.setMessage("Writing details array");
                _featureGrouper.writeArrayDetails(out, true);
                out.flush();
            }

            if (showCharts && featureSets.size() == 2)
            {
                List<Float> logInt1 = new ArrayList<Float>();
                List<Float> logInt2 = new ArrayList<Float>();


                for (Clusterer2D.BucketSummary summary :_featureGrouper.summarize())
                {
                    if ((summary.featureCount == summary.setCount) && (summary.setCount == 2))
                    {
                        logInt1.add((float)Math.log(
                                ((FeatureClusterer.FeatureClusterable)
                                        summary.getParentListForSetIndex(0).get(0)).getParentFeature().getTotalIntensity()));
                        logInt2.add((float)Math.log(
                                ((FeatureClusterer.FeatureClusterable)
                                        summary.getParentListForSetIndex(1).get(0)).getParentFeature().getTotalIntensity()));
                    }
                }
                if (logInt1.size() > 1)
                {
                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(logInt1, logInt2, "PerfectBucket_logint_nonorm");
                    pwsp.setPointSize(2);
                    pwsp.displayInTab();
                    ApplicationContext.setMessage("Perfect Bucket log-intensity correlation (without normalization): " +
                            BasicStatistics.correlationCoefficient(logInt1, logInt2));
                }

                if (_normalize)
                {
                    logInt1 = new ArrayList<Float>();
                    logInt2 = new ArrayList<Float>();

                    for (int i=0; i<intensitiesAfterProcessing.length; i++)
                    {
                         if (intensitiesAfterProcessing[i][0] > 0 && intensitiesAfterProcessing[i][1] > 0)
                         {
                             logInt1.add((float)Math.log(intensitiesAfterProcessing[i][0]));
                             logInt2.add((float)Math.log(intensitiesAfterProcessing[i][1]));
                         }
                    }
                    if (logInt1.size() > 1)
                    {
                        PanelWithScatterPlot pwsp = new PanelWithScatterPlot(logInt1, logInt2, "AllBucket_logint_norm");
                        pwsp.setPointSize(2);
                        pwsp.displayInTab();
                        ApplicationContext.setMessage("All Bucket log-intensity correlation (with normalization): " +
                                BasicStatistics.correlationCoefficient(logInt1, logInt2));
                    }

                }
            }
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Error building peptide array", x);
        }
        finally
        {
            ApplicationContext.setMessage("Peptide array complete." + (_outFileName == null ? "" : " See " + _outFileName + " and " + detailsFileName));
            if (out != null && _outFileName == null)
            {
                out.close();
//System.err.println("Closing out");
            }
        }
    }

    protected void optimizePerfectBuckets(List<FeatureSet> featureSets, int optimizationPerfectBucketMode, boolean showCharts)
    {
        ApplicationContext.setMessage("Optimizing");
        FeatureGrouper optimizeFeatureGrouper = new FeatureGrouper();
        optimizeFeatureGrouper.setGroupByCharge(false);
        optimizeFeatureGrouper.setGroupByMass(true);
        optimizeFeatureGrouper.setMassType(_featureGrouper.getMassType());
        optimizeFeatureGrouper.setConflictResolver(_conflictResolver);
        for (FeatureSet featureSet : featureSets)
        {
            optimizeFeatureGrouper.addSet(featureSet);
        }
        _log.debug("optimizing... added all sets");
        StringBuffer massBucketsToPrint = new StringBuffer();
        for (double massBucket : _massBuckets)
            massBucketsToPrint.append(massBucket + ", ");
        StringBuffer scanBucketsToPrint = new StringBuffer();
        for (double elutionBucket : _elutionBuckets)
            scanBucketsToPrint.append(elutionBucket + ", ");
        ApplicationContext.setMessage("Optimizing: mass buckets " + massBucketsToPrint +
                " elution buckets " + scanBucketsToPrint);

        Pair<Double, Double> bestBuckets =
                optimizeFeatureGrouper.calculateBestBuckets(_massBuckets, _elutionBuckets, optimizationPerfectBucketMode, showCharts);
        _massBucket = bestBuckets.first;
        _elutionBucket = bestBuckets.second;
        ApplicationContext.setMessage("Using mass and elution buckets " +
                _massBucket + " and " + _elutionBucket);
    }

    /**
     * 
     * @param featureSets
     * @param showCharts
     * @throws IOException
     */
    protected void optimizeDist(List<FeatureSet> featureSets, boolean showCharts)  throws IOException
    {
        ApplicationContext.setMessage("Beginning distribution optimization....");
        FeatureGrouper optimizeFeatureGrouper = new FeatureGrouper();
        optimizeFeatureGrouper.setGroupByCharge(true);

        FeatureSet.FeatureSelector accMassSel = new FeatureSet.FeatureSelector();
        accMassSel.setAccurateMzOnly(true);
        for (FeatureSet featureSet : featureSets)
            optimizeFeatureGrouper.addSet(featureSet);
        optimizeFeatureGrouper.setGroupByMass(true);
        optimizeFeatureGrouper.setMassType(_featureGrouper.getMassType());
        optimizeFeatureGrouper.setConflictResolver(_conflictResolver);

        optimizeFeatureGrouper.split2D(optimizationMassWideTolerance, optimizationRTWideTolerance);
        Clusterer2D.BucketSummary[] summaries = optimizeFeatureGrouper.summarize();

        float[] clusterRTDiameters = new float[summaries.length];
        float[] clusterMassDiameters = new float[summaries.length];
        for (int i=0; i<summaries.length; i++)
        {
            Clusterer2D.BucketSummary summary = summaries[i];
            if ((summary.featureCount < 2) || (summary.setCount != summary.featureCount))
                continue;
            clusterMassDiameters[i] = (float) (summary.maxDimension1 - summary.minDimension1);
            if (_featureGrouper.getMassType() == FeatureClusterer.DELTA_MASS_TYPE_PPM)
                    clusterMassDiameters[i] = MassUtilities.calculatePPMDeltaMass(
                            (float) ((summary.maxDimension1 + summary.minDimension1) / 2),
                            clusterMassDiameters[i], FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);

            clusterRTDiameters[i] =  (float) (summary.maxDimension2 - summary.minDimension2);
        }
        float massLimit = sdMultipleForWideToleranceCalc * BasicStatistics.standardDeviation(clusterMassDiameters);
        float rtLimit = sdMultipleForWideToleranceCalc * BasicStatistics.standardDeviation(clusterRTDiameters);


        ApplicationContext.setMessage(sdMultipleForWideToleranceCalc + "*SD mass diameter: " + massLimit + ", RT diameter = " + rtLimit);
        optimizeFeatureGrouper.split2D(massLimit, rtLimit);
        summaries = optimizeFeatureGrouper.summarize();
        int numPotentialDatapoints = 0;
        for (Clusterer2D.BucketSummary summary : summaries)
        {
            int numBucketPoints = summary.entries().length;
            if (summary.setCount >1 && (summary.featureCount == summary.setCount))
            {
                 numPotentialDatapoints += numBucketPoints * (numBucketPoints-1);
            }
        }

        float fractionDatapointsToKeep = 1f;
        if (numPotentialDatapoints > maxDatapointsForEMAlgorithm)
            fractionDatapointsToKeep = (float) maxDatapointsForEMAlgorithm / (float) numPotentialDatapoints;
        ApplicationContext.setMessage("Potential Datapoints: " + numPotentialDatapoints + ".  Keeping " +
                (fractionDatapointsToKeep * 100f) + "%");

        List<Float> massDistances = new ArrayList<Float>();
        List<Float> massesOfBuckets = new ArrayList<Float>();

        List<Float> rtDistances = new ArrayList<Float>();
        ApplicationContext.infoMessage("Calculating distances for " + summaries.length + " buckets");
        for (Clusterer2D.BucketSummary summary : summaries)
        {
            Feature[] bucketFeatures = optimizeFeatureGrouper.getFeatures(summary);


            if ((summary.setCount > 1) && (summary.featureCount == summary.setCount))
            {
                double[] masses = new double[bucketFeatures.length];
                double[] mzs = new double[bucketFeatures.length];
double[] rts = new double[bucketFeatures.length];


                for (int i=0; i<bucketFeatures.length; i++)
                {
                    masses[i] = bucketFeatures[i].mass;
                    mzs[i] = bucketFeatures[i].mz;
rts[i] = bucketFeatures[i].time;
                }
                double medianMass = BasicStatistics.median(masses);
                List<Float> massDistancesThisBucket = new ArrayList<Float>();
                List<Float> rtDistancesThisBucket = new ArrayList<Float>();
                for (int i=0; i<bucketFeatures.length-1; i++)
                {
                    for (int j=i+1; j<bucketFeatures.length; j++)
                    {
                        if (fractionDatapointsToKeep == 1 || Math.random() < fractionDatapointsToKeep)
                        {
                            massesOfBuckets.add((float) medianMass);
                            float deltaMass = bucketFeatures[i].getMass() - bucketFeatures[j].getMass();
                            if (_featureGrouper.getMassType() == FeatureClusterer.DELTA_MASS_TYPE_PPM)
                                deltaMass = MassUtilities.calculatePPMDeltaMass((float)medianMass,
                                    deltaMass, FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);
                            float deltaTime = _featureGrouper.getElution(bucketFeatures[i]) -
                                    _featureGrouper.getElution(bucketFeatures[j]);
                            //these features are ordered by time, secondarily by mass.  Need to randomize
                            //the ordering so that we don't end up with all negative or all pos values
                            boolean flipMassAndRT = Math.random() < 0.5f;
                            if (flipMassAndRT)
                            {
                                deltaMass = -deltaMass;
                                deltaTime = -deltaTime;
                            }

                            massDistancesThisBucket.add(deltaMass);
                            rtDistancesThisBucket.add(deltaTime);
                        }
                    }
                }

                massDistances.addAll(massDistancesThisBucket);
                rtDistances.addAll(rtDistancesThisBucket);

//if (rtDistancesThisBucket.get(0) < 10)
//    System.err.println(massDistancesThisBucket.get(0) + "\t" + masses[0] + "\t" + masses[1] + "\t" + mzs[0] + "\t" + mzs[1] + "\t" + rts[0] + "\t" + rts[1]);
            }
        }

        ApplicationContext.infoMessage("Calculating mixed model distribution for probability inference...");
        float minMassD = (float) BasicStatistics.min(massDistances);
        float maxMassD = (float) BasicStatistics.max(massDistances);
        float minRTD = (float) BasicStatistics.min(rtDistances);
        float maxRTD = (float) BasicStatistics.max(rtDistances);
//System.err.println("****" + minMassD + ", " + maxMassD + ", " + minRTD + ", " + maxRTD);        
        AmtMatchProbabilityAssigner probabilityAssigner = new AmtMatchProbabilityAssigner(
                minMassD, maxMassD, minRTD, maxRTD, 0.001f, maxMatchFDRForToleranceBoxCalc);
        float[] probabilities =
                probabilityAssigner.calculateProbabilitiesEM(massDistances, rtDistances, initialDistAssumptionTrue,
                        showCharts);
        float highestMassDiffWithHighProb = Float.MIN_VALUE;
        float highestRTDiffWithHighProb = Float.MIN_VALUE;
        for (int i=0; i<probabilities.length; i++)
        {
            if (probabilities[i] > maxMatchFDRForToleranceBoxCalc)
            {
                highestMassDiffWithHighProb = Math.max(highestMassDiffWithHighProb, Math.abs(massDistances.get(i)));
                highestRTDiffWithHighProb = Math.max(highestRTDiffWithHighProb, Math.abs(rtDistances.get(i)));
            }
        }
        ApplicationContext.infoMessage("Bounding box with passing features: mass=+-" +
                highestMassDiffWithHighProb + ", RT=+-" + highestRTDiffWithHighProb);
        _massBucket = highestMassDiffWithHighProb;
        _elutionBucket = highestRTDiffWithHighProb;
        ApplicationContext.setMessage("Using mass and elution buckets " +
                _massBucket + " and " + _elutionBucket);
        ApplicationContext.infoMessage("Distance dist datapoints: " + massDistances.size());
        if (showCharts)
        {
            PanelWithScatterPlot pwspDist = new PanelWithScatterPlot(rtDistances, massDistances, "Error distances");
            pwspDist.setPointSize(2);
            pwspDist.addVerticalLine(-highestRTDiffWithHighProb, -highestMassDiffWithHighProb, highestMassDiffWithHighProb);
            pwspDist.addVerticalLine(highestRTDiffWithHighProb, -highestMassDiffWithHighProb, highestMassDiffWithHighProb);
            pwspDist.addHorizontalLine(-highestMassDiffWithHighProb, -highestRTDiffWithHighProb, highestRTDiffWithHighProb);
            pwspDist.addHorizontalLine(highestMassDiffWithHighProb, -highestRTDiffWithHighProb, highestRTDiffWithHighProb);
            pwspDist.setSeriesColor(0, Color.RED);
            for (int i=1; i<6; i++)
                pwspDist.setSeriesColor(i, Color.BLUE);

            pwspDist.displayInTab();

            new PanelWithScatterPlot(massesOfBuckets, massDistances, "mass vs deltaMass").displayInTab();
        }
    }



    /**
     * Figure out the filepath for the details file related to an array file
     * @param arrayFilepath
     * @return
     */
    public static String calcDetailsFilepath(String arrayFilepath)
    {
        int dotPos = arrayFilepath.lastIndexOf('.');
        String detailsFileName = null;
        if (dotPos < 0)
            detailsFileName = arrayFilepath + ".details.tsv";
        else
            detailsFileName = arrayFilepath.substring(0, dotPos) + ".details" + arrayFilepath.substring(dotPos);
        return detailsFileName;
    }

    /**
     * Look up the tag for the given feature set. Tags a are supplied as key/value pairs,
     * where the key is the prefix of the last element of the file path (everything from
     * the first "." on is stripped off).
     */
    private static String lookupTag(FeatureSet fs, Map<String, String> tags)
    {
        if (null == fs || null == fs.getSourceFile())
            return null;

        String filename = fs.getSourceFile().getName();

        if (tags.containsKey(filename))
            return tags.get(filename);

        return null;
    }

    /**
     * Group feature sets by tag. 
     * Align within group to generate a representative FeatureSet for each.
     * Then align the representatives.
     */
    public void alignByGroup(Map<String,String> tags, boolean strict)
    {
        String detailsFileName = null;

        PrintWriter out = null;
        if (_sel == null)
            _sel = new FeatureSet.FeatureSelector();

        try
        {
            ApplicationContext.setMessage("Loading");
            List<FeatureSet> sourceFeatureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < _sets.size(); i++)
            {
                FeatureSet fs;
                if (_sets.get(i) instanceof File)
                    fs = new FeatureSet((File) _sets.get(i), Color.RED);
                else if (_sets.get(i) instanceof FeatureSet)
                    fs = (FeatureSet) _sets.get(i);
                else if (_sets.get(i) instanceof String)
                    fs = new FeatureSet(new File((String) _sets.get(i)), Color.RED);
                else
                {
                    ApplicationContext.errorMessage("Couldn't load feature set due to bad object: " + _sets.get(i).toString(), null);
                    return;
                }

                sourceFeatureSets.add(fs);
            }

            // ???? At this point, if there is only one tag, we could revert to regular alignment

            // Filter all source feature sets and add each to the appropriate tagged list
            ApplicationContext.setMessage("Filtering");

            HashMap<String, List<FeatureSet>> taggedFeatureSets = new HashMap<String, List<FeatureSet>>();
            ArrayList<String> tagList = new ArrayList<String>(); // Try to keep order of big alignment consistent
            for (int i = 0; i < sourceFeatureSets.size(); i++)
            {
                FeatureSet fs = ((FeatureSet) sourceFeatureSets.get(i)).filter(_sel);
                fs.setTag(lookupTag(fs, tags));
                if (!taggedFeatureSets.containsKey(fs.getTag()))
                    taggedFeatureSets.put(fs.getTag(), new ArrayList<FeatureSet>());
                taggedFeatureSets.get(fs.getTag()).add(fs);
                if (!tagList.contains(fs.getTag()))
                    tagList.add(fs.getTag());
            }

            // Loop through the tags and align all associated FeatureSets
            ArrayList<FeatureSet> filteredFeatureSets = new ArrayList<FeatureSet>();
            for (String tag : tagList)
            {
                List<FeatureSet> featureSets = taggedFeatureSets.get(tag);

                // If only one FeatureSet in this group just add it to the list for the next round
                if (featureSets.size() == 1)
                {
                    filteredFeatureSets.add(featureSets.get(0));
                }
                else
                {
                    //Align
                    ApplicationContext.setMessage("Aligning group " + tag);
                    _aligner.setFeaturePairSelector(_featurePairSelector);
                    featureSets = _aligner.alignFeatureSets(featureSets, false);

                    //Deconvolve
                    if (_shouldDeconvolute)
                    {
                        for (int i = 0; i < featureSets.size(); i++)
                        {

                            featureSets.set(i, featureSets.get(i).deconvolute(
                                    _deconvoluteScanWindow, _deconvoluteMassWindow, true));
                        }
                    }
            
                    FeatureGrouper grouper = new FeatureGrouper();
            
                    grouper.setGroupByMass(true);
                    grouper.setConflictResolver(_conflictResolver);
                    for (FeatureSet fs : featureSets)
                    {
                        grouper.addSet(fs);
                    }
            
                    ApplicationContext.setMessage("Creating array " + tag);
                    grouper.split2D(_massBucket, _elutionBucket);

                    // If strict, a feature must match across all runs to be included.
                    // Otherwise ~75% of the runs is good enough
                    int minAligned = strict ? featureSets.size() :(int) (featureSets.size() * 3.0/4.0 + 0.5);
                    for (FeatureSet fs : grouper.filterByGroupedAlignment(minAligned))
                        filteredFeatureSets.add(fs);
                }
            }

            // Now create a new BucketedPeptideArray to align the filtered feature sets.
            // Note that we've already applied any other filters during the within group
            // alignment, so there should be no need for a selector here.
            BucketedPeptideArray arr = new BucketedPeptideArray(filteredFeatureSets);
            arr.setElutionBucket(_elutionBucket);
            arr.setMassBucket(_massBucket);
            arr.setOutFileName(_outFileName);
            arr.setAlign(true);
            arr.setNormalize(_normalize);
            arr.setConflictResolver(_conflictResolver);
            arr.run();
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Error building peptide array", x);
        }
    }

    private void writeHeader(PrintWriter out)
    {
        String revision = (String) ApplicationContext.getProperty("REVISION");
        out.println("# Revision: " + revision);
        out.print("# Files:");
        for (int i = 0; i < _sets.size(); i++)
            out.print(" " + (i + 1) + ": " + getSetName(_sets.get(i)));
        out.println();
        out.print("# Params: " + _sel.toString());
        out.println(" --massWindow=" + _massBucket + " --scanWindow=" + _elutionBucket + (_normalize ? " --normalize" : ""));
        // ???? If normalize is true, it would be good to write the individual scales here
    }

    private String getSetName(Object set)
    {
        if (set instanceof File)
            return ((File) set).getName();
        else if (set instanceof FeatureSet)
            return ((FeatureSet) set).getSourceFile() != null ? ((FeatureSet) set).getSourceFile().getName() : set.toString();
        else
            return set.toString();
    }

    public java.util.List<?> getFiles()
    {
        return _sets;
    }

    public void setFiles(java.util.List<?> files)
    {
        this._sets = files;
    }

    public FeatureSet.FeatureSelector get_sel()
    {
        return _sel;
    }

    public void set_sel(FeatureSet.FeatureSelector _sel)
    {
        this._sel = _sel;
    }

    public boolean isAlign()
    {
        return _align;
    }

    public void setAlign(boolean align)
    {
        this._align = align;
    }

    public boolean getNormalize()
    {
        return _normalize;
    }

    public void setNormalize(boolean normalize)
    {
        this._normalize = normalize;
    }

    public int getDeconvoluteScanWindow()
    {
        return _deconvoluteScanWindow;
    }

    public void setDeconvoluteScanWindow(int deconvoluteScanWindow)
    {
        this._deconvoluteScanWindow = deconvoluteScanWindow;
    }

    public double getDeconvoluteMassWindow()
    {
        return _deconvoluteMassWindow;
    }

    public void setDeconvoluteMassWindow(double deconvoluteMassWindow)
    {
        this._deconvoluteMassWindow = deconvoluteMassWindow;
    }

    public double getMassBucket()
    {
        return _massBucket;
    }

    public void setMassBucket(double massBucket)
    {
        this._massBucket = massBucket;
    }

    public double getElutionBucket()
    {
        return _elutionBucket;
    }

    public void setElutionBucket(double scanBucket)
    {
        this._elutionBucket = scanBucket;
    }

    public String getOutFileName()
    {
        return _outFileName;
    }

    public void setOutFileName(String outFileName)
    {
        this._outFileName = outFileName;
    }

    public double[] getElutionBuckets()
    {
        return _elutionBuckets;
    }

    public void setScanBuckets(double[] _scanBuckets)
    {
        this._elutionBuckets = _scanBuckets;
    }

    public double[] getMassBuckets()
    {
        return _massBuckets;
    }

    public void setMassBuckets(double[] _massBuckets)
    {
        this._massBuckets = _massBuckets;
    }

    public int getConflictResolver()
    {
        return _conflictResolver;
    }

    public void setConflictResolver(int conflictResolver)
    {
        this._conflictResolver = conflictResolver;
    }

    public FeatureGrouper getFeatureGrouper()
    {
        return _featureGrouper;
    }

    public void setFeatureGrouper(FeatureGrouper _featureGrouper)
    {
        this._featureGrouper = _featureGrouper;
    }

    public Aligner.FeaturePairSelector getFeaturePairSelector()
    {
        return _featurePairSelector;
    }

    public void setFeaturePairSelector(Aligner.FeaturePairSelector featurePairSelector)
    {
        this._featurePairSelector = featurePairSelector;
    }

    public Aligner getAligner()
    {
        return _aligner;
    }

    public void setAligner(Aligner _aligner)
    {
        this._aligner = _aligner;
    }

    public boolean isShouldDeconvolute()
    {
        return _shouldDeconvolute;
    }

    public void setShouldDeconvolute(boolean _shouldDeconvolute)
    {
        this._shouldDeconvolute = _shouldDeconvolute;
    }
    
    public void setElutionMode(int _elutionMode)
    {
        _featureGrouper.setElutionMode(_elutionMode);
    }

    public int getOptimizationMode()
    {
        return optimizationMode;
    }

    public void setOptimizationMode(int optimizationMode)
    {
        this.optimizationMode = optimizationMode;
    }

    /**
     * Also sets the wide mass tolerance for EM dist optimization to a default.  So if you're overriding, call this,
     * then call setOptimizationMassWideTolerance
     * @param massType
     */
    public void setMassType(int massType)
    {
        _featureGrouper.setMassType(massType);
        switch (massType)
        {
            case FeatureClusterer.DELTA_MASS_TYPE_ABSOLUTE:
                optimizationMassWideTolerance = DEFAULT_MASS_WIDE_TOLERANCE_DA;
                _massBuckets = DEFAULT_MASS_BUCKETS_DA;
                break;
            case FeatureClusterer.DELTA_MASS_TYPE_PPM:
                optimizationMassWideTolerance = DEFAULT_MASS_WIDE_TOLERANCE_PPM;
                _massBuckets = DEFAULT_MASS_BUCKETS_PPM;
                break;
        }
    }

    public float getOptimizationRTWideTolerance()
    {
        return optimizationRTWideTolerance;
    }

    public void setOptimizationRTWideTolerance(float optimizationRTWideTolerance)
    {
        this.optimizationRTWideTolerance = optimizationRTWideTolerance;
    }

    public float getOptimizationMassWideTolerance()
    {
        return optimizationMassWideTolerance;
    }

    public void setOptimizationMassWideTolerance(float optimizationMassWideTolerance)
    {
        this.optimizationMassWideTolerance = optimizationMassWideTolerance;
    }

    public int getMaxDatapointsForEMAlgorithm()
    {
        return maxDatapointsForEMAlgorithm;
    }

    public void setMaxDatapointsForEMAlgorithm(int maxDatapointsForEMAlgorithm)
    {
        this.maxDatapointsForEMAlgorithm = maxDatapointsForEMAlgorithm;
    }

    public float getMaxMatchFDRForToleranceBoxCalc()
    {
        return maxMatchFDRForToleranceBoxCalc;
    }

    public void setMaxMatchFDRForToleranceBoxCalc(float minMatchProbForToleranceBoxCalc)
    {
        this.maxMatchFDRForToleranceBoxCalc = minMatchProbForToleranceBoxCalc;
    }
}
