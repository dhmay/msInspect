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

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.viewer.amt.AmtUtilities;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseMatcher;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.plot.PlotOrientation;

import java.io.*;
import java.util.*;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: migra
 * Date: Mar 1, 2005
 * Time: 10:10:36 AM
 *
 * dhmay, 5/7/08.  Making this class abstract.  Alignment can be performed via two different algorithms: a spline-based
 * regression, or a quantile regression similar to what's used in AMT.  The two algorithms perform very similarly
 * on data that are tightly filtered, but the quantile regression performs much better on noisy data.
 *
 * dhmay, 20100310.  Abstracting "elution" so that mapping can be from scan to scan or from time to time.
 */
public abstract class Aligner
{
    private static Logger _log = Logger.getLogger(Aligner.class);

    public static final FeaturePairSelector DEFAULT_FEATURE_PAIR_SELECTOR =
            new MzFeaturePairSelector();


    protected FeaturePairSelector featurePairSelector = DEFAULT_FEATURE_PAIR_SELECTOR;

    protected List<int[]> warpingMaps = null;

    protected boolean buildCharts = false;

    protected PanelWithChart alignmentChart = null;

    double[] mappingSourceValues = null;
    double[] mappingMappedValues = null;

    protected double maxStudRes = AmtDatabaseMatcher.DEFAULT_MAX_STUDENTIZED_RESIDUAL;
    protected double maxLeverageNumerator = AmtDatabaseMatcher.DEFAULT_LEVERAGE_NUMERATOR;

    //dhmay adding 20100310.  Different sources of data for alignment
    public static final int ALIGNMENT_DATASOURCE_TIME = 0;
    public static final int ALIGNMENT_DATASOURCE_SCAN = 1;

    public static final int ALIGNMENT_DATASOURCE_DEFAULT = ALIGNMENT_DATASOURCE_TIME;

    protected int alignmentDatasource = ALIGNMENT_DATASOURCE_DEFAULT;



    //dhmay adding 20100201.  Defining different ways of ordering the runs for alignment.

    //align all runs to the first run.  "Legacy" behavior
    public static final int ALIGNMENT_ORDER_MODE_ALIGNTOFIRST = 0;
    //daisychain the alignment.  Align run 1 to run 0, run 2 to run 1....
    public static final int ALIGNMENT_ORDER_MODE_DAISYCHAIN = 1;
    //accumulate a 'cumulative' dataset that we align to
    public static final int ALIGNMENT_ORDER_MODE_CUMULATIVE = 2;

    public static final int[] alignmentOrderModes = new int[]
            {
                ALIGNMENT_ORDER_MODE_ALIGNTOFIRST,
                ALIGNMENT_ORDER_MODE_DAISYCHAIN
            };
    public static final int ALIGNMENT_ORDER_MODE_DEFAULT = ALIGNMENT_ORDER_MODE_ALIGNTOFIRST;

    public static final String[] alignmentOrderModeDescs = new String[]
            {
                    "alltofirst",
                    "daisychain",
                    "cumulative",
            };

    protected int alignmentOrderMode = ALIGNMENT_ORDER_MODE_DEFAULT;

    public Aligner()
    {
    }

    /**
     * Calculate the maximum time of any of these features
     * @param features
     * @return
     */
    public double getMaxElution(Feature[] features)
    {
        double maxElution = 0;
        for (Feature feature : features)
        {
            if (getFeatureElutionValue(feature) > maxElution)
                maxElution = getFeatureElutionValue(feature);
        }
        return maxElution;
    }

    protected double getFeatureElutionValue(Feature feature)
    {
        switch (alignmentDatasource)
        {
            case ALIGNMENT_DATASOURCE_SCAN:
                return feature.getScan();
            default:
                return feature.getTime();
        }
    }

    /**
     * If scan mode, get time; if time mode, get scan
     * @param feature
     * @return
     */
    protected double getFeatureOppositeElutionValue(Feature feature)
    {
        switch (alignmentDatasource)
        {
            case ALIGNMENT_DATASOURCE_TIME:
                return feature.getScan();
            default:
                return feature.getTime();
        }
    }

    /**
     * Build a map from scan to RT or RT to scan, based on an array of features. In the resulting array,
     * index 0 represents minValueForResult, last index = maxValueForResult. minScanForResult can be negative.
     * Fill in values before first and after last.
     *
     * Between scan and time, which is the index and which is the value depends on alignmentDatasource:
     * SCAN: scan is index, time is value
     * TIME: time is index, scan is value
     * @param features
     * @return
     */
    protected float[] createElutionTypeMap(Feature[] features, int minIndexValForResult, int maxIndexValForResult, boolean showCharts)
    {
        Feature[] featuresCopy = new Feature[features.length];
        System.arraycopy(features, 0, featuresCopy, 0, featuresCopy.length);
        
        Arrays.sort(featuresCopy, new Feature.ScanAscComparator());

        //Calculate mean amount of time between scans.  We will use this to fill in gaps
        List<Float> elutionValuesBetweenTicks = new ArrayList<Float>();

        for (int i=1; i<featuresCopy.length; i++)
        {
            int numTicksBetween = (int) (getFeatureElutionValue(featuresCopy[i]) - getFeatureElutionValue(featuresCopy[i-1]));
            if (numTicksBetween > 0)
            {
                float valueBetweenTicks =  (float)
                        (getFeatureOppositeElutionValue(featuresCopy[i]) - getFeatureOppositeElutionValue(featuresCopy[i-1]));
                elutionValuesBetweenTicks.add(valueBetweenTicks / numTicksBetween);
            }
        }
        float meanValueBetweenTicks = (float) BasicStatistics.mean(elutionValuesBetweenTicks);

        int numScansInResult = maxIndexValForResult - minIndexValForResult + 1;

        float[] result = new float[numScansInResult];
        int featureIndex = 0;

        List<Integer> xValsForChart = new ArrayList<Integer>();
        List<Float> yValsForChart = new ArrayList<Float>();
        for (int i=0; i<numScansInResult; i++)
        {
            int integerizedValue = i + minIndexValForResult;

            //advance the feature array pointer to the first feature that has a higher or equivalent scan number to i
            while (integerizedValue>getFeatureElutionValue(featuresCopy[featureIndex]) && featureIndex < featuresCopy.length-1)
                featureIndex++;
            Feature firstSameOrHigherScanFeature = featuresCopy[featureIndex];
            //Assign the RT for i the value RT value for this feature.  But if the feature scan != i, adjust
            //by adding/subtracting the appropriate number of meanTimeBetweenScans
            float mappedVal = (float)getFeatureOppositeElutionValue(firstSameOrHigherScanFeature);
            result[i] =  mappedVal +
                   (float) ((integerizedValue - getFeatureElutionValue(featuresCopy[featureIndex])) * meanValueBetweenTicks);
            //For really nonlinear relationships, we might have created a situation in which one scan's value is lower
            //than the previous.  Bump up.
            if (i>0) result[i] =
                    Math.max(result[i-1],result[i]);
//System.err.println(i + ", " + result[i]);            

            if (showCharts)
                xValsForChart.add(integerizedValue); yValsForChart.add(result[i]);
        }
        if (showCharts)
        {
            PanelWithScatterPlot pwsp = new PanelWithScatterPlot(xValsForChart, yValsForChart, "scan vs rt, base run");
            pwsp.setAxisLabels("scan","time");
            if (alignmentDatasource == ALIGNMENT_DATASOURCE_TIME)
                pwsp.setAxisLabels("time","scan");
            pwsp.displayInTab();
        }
        return result;
    }

    /**
     * Perform the non-linear mapping between feature sets
     */
    public List<FeatureSet> alignFeatureSets(List<FeatureSet> sets, boolean showCharts)
    {
        warpingMaps = new ArrayList<int[]>(sets.size()-1);


        if (showCharts)
            buildCharts=true;

        List<FeatureSet> alignedSets = new ArrayList<FeatureSet>();
        alignedSets.add(sets.get(0));

        _log.debug("Align: base featureset  has " + sets.get(0).getFeatures().length + " features");

        double minElutionAllSetsPostAlign = Double.MAX_VALUE;
        double maxElutionAllSetsPostAlign = 0;
        for (Feature feature : sets.get(0).getFeatures())
        {
            double elutionVal = getFeatureElutionValue(feature);
            if (elutionVal > maxElutionAllSetsPostAlign)
                maxElutionAllSetsPostAlign = elutionVal;
            if (elutionVal < minElutionAllSetsPostAlign)
                minElutionAllSetsPostAlign = elutionVal;
        }

        //dhmay adding 20100201
        //keep track of the FeatureSet we're aligning to.  For align-to-first mode, always the first set.
        //For daisychain mode, the most recent (post-alignment) set
        //todo: if adding a new alignment mode, need to mess with this
        FeatureSet alignToSet = sets.get(0);
        if (alignmentOrderMode == ALIGNMENT_ORDER_MODE_CUMULATIVE)
            alignToSet = alignToSet.deepCopy();



        for (int i=1; i<sets.size(); i++)
        {
            String sourceFilename = "???";
            if (sets.get(i).getSourceFile() != null)
                sourceFilename = sets.get(i).getSourceFile().getAbsolutePath();
            _log.debug("Aligning Feature File " + sourceFilename + ", features: " + sets.get(i).getFeatures().length);
            Pair<Feature,Feature>[] pairedFeatures =
                    featurePairSelector.selectPairs(sets.get(i),alignToSet);
            _log.debug("Aligner.alignFeatureSets: selected " + pairedFeatures.length +
                       " pairs");

            double maxElution = 0;
            for (Feature feature : sets.get(i).getFeatures())
            {
                if (getFeatureElutionValue(feature) > maxElution)
                    maxElution = getFeatureElutionValue(feature);
            }


            //necessary cast
            Pair<Double,Double>[] pairedElutions = (Pair<Double,Double>[])
                    new Pair[pairedFeatures.length];
            for (int j=0; j<pairedElutions.length; j++)
            {
                pairedElutions[j] =
                        new Pair<Double,Double>(getFeatureElutionValue(pairedFeatures[j].first),
                                                getFeatureElutionValue(pairedFeatures[j].second));
            }

            Pair<Double,Double>[] restrictedPairs = restrictPairsForRegression(pairedElutions);
            _log.debug("Aligner.alignFeatureSets: restricted to " +
                    restrictedPairs.length + " pairs");

            String mapFileNameStart = (i+1) + "_onto_1";
            //this is where the work is done
            double[] elutionMap = alignPairs(restrictedPairs, maxElution, mapFileNameStart);

            //It's lame to take up extra memory this way, but it's nice to be
            //able to persist this scan map for later investigation
            int[] intElutionMap = new int[elutionMap.length];
            for (int j=0; j<elutionMap.length; j++)
                intElutionMap[j] = (int) Math.round(elutionMap[j]);
            warpingMaps.add(intElutionMap);

            FeatureSet aligned = sets.get(i).deepCopy();

            //indicate that this featureset has been aligned by putting "align"
            //in the source filename.  OK, so this is a bit of a hack for backward
            //compatibility of the detail files
            aligned.setSourceFile(new File(getAlignedFileName(sets, i)));


            for (Feature feature : aligned.getFeatures())
            {
                double newElution = intElutionMap[(int) getFeatureElutionValue(feature)];
                switch (alignmentDatasource)
                {
                    case ALIGNMENT_DATASOURCE_SCAN:
                        feature.setScan((int) newElution);
                        break;
                    default:
                        feature.setTime((float) newElution);
                        break;
                }
                if (newElution > maxElutionAllSetsPostAlign)
                    maxElutionAllSetsPostAlign = newElution;
                if (newElution < minElutionAllSetsPostAlign)
                {
                    minElutionAllSetsPostAlign = newElution;
                }
            }
            alignedSets.add(aligned);

            //dhmay adding 20100201
            //If doing daisy-chain alignment, save the current aligned set as the one to align the next one to
            //todo: if adding a new alignment mode, need to mess with this
            if (alignmentOrderMode == ALIGNMENT_ORDER_MODE_DAISYCHAIN)
                alignToSet = aligned;
            else if (alignmentOrderMode == ALIGNMENT_ORDER_MODE_CUMULATIVE)
            {
                Feature[] cumFeatures = new Feature[alignToSet.getFeatures().length + aligned.getFeatures().length];
                System.arraycopy(alignToSet.getFeatures(), 0, cumFeatures,0, alignToSet.getFeatures().length);
                System.arraycopy(aligned.getFeatures(), 0, cumFeatures, alignToSet.getFeatures().length, aligned.getFeatures().length);

                alignToSet.setFeatures(cumFeatures);
            }

            if (buildCharts)
            {
                createChart(pairedElutions, restrictedPairs, elutionMap);

                if (showCharts)
                {
                    alignmentChart.setName("Alignment " + i);
                    alignmentChart.displayInTab();

                    double[] massesPPM = new double[pairedFeatures.length];
                    double[] deltaMassesPPM = new double[pairedFeatures.length];
                    for (int j=0; j<deltaMassesPPM.length; j++)
                    {
                        deltaMassesPPM[j] = (pairedFeatures[j].second.getMass() -
                                pairedFeatures[j].first.getMass()) *
                                1000000 / pairedFeatures[j].first.getMass();
                        massesPPM[j] = pairedFeatures[j].first.getMass();
                    }
                    //todo: fix these so extreme values are gone, or just get rid of them permanently
//                    PanelWithHistogram deltaMassHist = new PanelWithHistogram(deltaMassesPPM);
//                    deltaMassHist.setName("massMatchDeltaPPM");
//                    deltaMassHist.displayInTab();
//
//
//                    PanelWithScatterPlot massMassDevPlot = new PanelWithScatterPlot(
//                             massesPPM, deltaMassesPPM, "Mass (x) vs. deltaMass (y), ppm");
//                    massMassDevPlot.displayInTab();
                }
            }

        }

        //maps from scan to time, or from time to scan, depending on alignmentDatasource.
        // SCAN: scan is index, time is value
        // TIME: time is index, scan is value
        //This is so that both values will be in sync.
        float[] baseElutionTypeMap = createElutionTypeMap(sets.get(0).getFeatures(),
                (int) minElutionAllSetsPostAlign, (int) maxElutionAllSetsPostAlign, showCharts);
        for (FeatureSet featureSet : alignedSets)
        {
            for (Feature feature : featureSet.getFeatures())
            {
                switch (alignmentDatasource)
                {
                    case ALIGNMENT_DATASOURCE_SCAN:
                        feature.setTime(baseElutionTypeMap[(int) (feature.getScan()-minElutionAllSetsPostAlign)]);
                        break;
                    default:
                        feature.setScan((int) baseElutionTypeMap[(int)(feature.getTime()-minElutionAllSetsPostAlign)]);
                        break;
                }
            }
        }


        //clean up
        TempFileManager.deleteTempFiles(this);

        return alignedSets;
    }

    protected Pair<Double,Double>[] restrictPairsForRegression(Pair<Double,Double>[] allPairs)
    {
        Feature[] featuresForRegression = new Feature[allPairs.length];

        double[] baseTimes = new double[allPairs.length];
        double[] toAlignTimes = new double[allPairs.length];
        //this is a bit silly -- we are modeling this as an AMT regression, so we're setting Time and
        //Observed Hydrophobicity appropriately in a dummy feature
        for (int j=0; j<allPairs.length; j++)
        {
            Pair<Double,Double> matchedPair = allPairs[j];
            Feature dummyFeature = new Feature();
            baseTimes[j] = matchedPair.first;
            toAlignTimes[j] = matchedPair.second;
            dummyFeature.setTime((float) baseTimes[j]);
            AmtExtraInfoDef.setObservedHydrophobicity(dummyFeature,
                    toAlignTimes[j]);
            featuresForRegression[j] = dummyFeature;
        }
        double[] baseTimesForRegression = new double[featuresForRegression.length];
        double[] toAlignTimesForRegression = new double[featuresForRegression.length];

        for (int j=0; j<featuresForRegression.length; j++)
        {
            //this call to getTime() is OK -- this is the dummy value for regression, which could be time or scan
            baseTimesForRegression[j] = featuresForRegression[j].getTime();
            toAlignTimesForRegression[j] =
                    AmtExtraInfoDef.getObservedHydrophobicity(featuresForRegression[j]);
        }

        featuresForRegression =
                ProteomicsRegressionUtilities.selectFeaturesWithLowLeverageAndStudentizedResidual(
                        featuresForRegression,
                        baseTimes, toAlignTimes,
                        maxLeverageNumerator,
                        maxStudRes,
                        false, 0, false, true);
        Pair<Double,Double>[] result = (Pair<Double,Double>[])
                new Pair[featuresForRegression.length];


        _log.debug(featuresForRegression.length + " pairs left after exclusion");

        for (int j=0; j<featuresForRegression.length; j++)
        {
            //this call to getTime() is OK -- this is the dummy value for regression, which could be time or scan            
            result[j] =  new Pair<Double,Double>((double) featuresForRegression[j].getTime(),
                    AmtExtraInfoDef.getObservedHydrophobicity(featuresForRegression[j]));
        }
        ApplicationContext.setMessage("Using " + featuresForRegression.length +
                " features (out of " + allPairs.length +
                " mass-matched) for regression");
        return result;
    }


    /**
     * Cover method. In the absence of filename information, use a filename that's
     * likely to be unique
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public double[] alignPairs(Pair<Double,Double>[] pairs, int maxValueToWarp)
    {
        return alignPairs(pairs, maxValueToWarp, "length" + pairs.length);
    }

    /**
     * Generic method for creating a warping based on pairs of values.
     * Create temp files with names based on tempFileNameStart, in case the user
     * wants to peruse them
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public abstract double[] alignPairs(Pair<Double,Double>[] pairs, double maxValueToWarp,
                               String tempFileNameStart);



    protected PanelWithChart createChart(Pair<Double,Double>[] allPairedScans,
                                         Pair<Double,Double>[] restrictedPairs,
                                         double[] warpedValues)
    {
        PanelWithScatterPlot pwsp = new PanelWithScatterPlot();

        double[] restrictedScans1 = new double[restrictedPairs.length];
        double[] restrictedScans2 = new double[restrictedPairs.length];

        for (int i=0; i<restrictedScans1.length; i++)
        {
            restrictedScans1[i] = restrictedPairs[i].first;
            restrictedScans2[i] = restrictedPairs[i].second;
        }
        pwsp.addData(restrictedScans1, restrictedScans2, "Used in regression");


        double[] scans1 = new double[allPairedScans.length];
        double[] scans2 = new double[allPairedScans.length];

        double minScan1 = Double.MAX_VALUE;
        double maxScan1 = Double.MIN_VALUE;

        for (int i=0; i<scans1.length; i++)
        {
            scans1[i] = allPairedScans[i].first;
            scans2[i] = allPairedScans[i].second;

            if (scans1[i] < minScan1)
                minScan1 = scans1[i];
            if (scans1[i] > maxScan1)
                maxScan1 = scans1[i];
        }
        pwsp.addData(scans1, scans2, "all");

        int firstBaseVal = Math.max((int) minScan1 - 100, 0);
        int lastBaseVal = Math.min((int) maxScan1 + 100, warpedValues.length-1);

        double[] plottedWarpedValues = new double[lastBaseVal - firstBaseVal + 1];
        System.arraycopy(warpedValues, firstBaseVal, plottedWarpedValues, 0, plottedWarpedValues.length);

        double[] unwarpedValues = new double[plottedWarpedValues.length];
        for (int i=0; i<plottedWarpedValues.length; i++)
            unwarpedValues[i] = i + firstBaseVal;



        pwsp.addData(unwarpedValues, plottedWarpedValues, "alignment map");

        pwsp.setAxisLabels("Time in run 1","Time in run 2");

        alignmentChart = pwsp;
        return pwsp;
    }

    /**
     * Write a file containing pairs
     * @param fileName
     * @throws FileNotFoundException
     */
    protected File writePairsFile(Pair<Double,Double>[] pairs, String fileName)
            throws FileNotFoundException
    {
        File currentPairsFile =
                    TempFileManager.createTempFile(fileName,
                                                   this);

        PrintWriter currentPairsPW = new PrintWriter(currentPairsFile);
        currentPairsPW.println("source\tdest");
        for (Pair<Double,Double> pair : pairs)
            currentPairsPW.println(pair.first + "\t" + pair.second);

        currentPairsPW.flush();
        _log.debug("wrote feature pairs to file " + currentPairsFile.getAbsolutePath());
        return currentPairsFile;
    }

    /**
     * Return a filename with an "align" in it.  If there's a source file for the
     * featureset, use that name to start with.  If not, index.
     * @param featureSets
     * @param i
     * @return
     */
    protected String getAlignedFileName(List featureSets, int i)
    {
        FeatureSet fs = (FeatureSet) featureSets.get(i);
        String baseFileName;
        String fileName = fs.getSourceFile() == null ? "set" + i + ".tsv" : fs.getSourceFile().getName();

        int dotPos = fileName.lastIndexOf('.');
        if (dotPos < 0)
            baseFileName = fileName;
        else
            baseFileName = fileName.substring(0, dotPos);

        return baseFileName + ".align.tsv";
    }

    /**
     * Provide access to the scan maps created for individual feature sets.
     * Obviously, only applicable if you've just done alignment on a bunch of
     * feature sets
     * @return
     */
    public List<int[]> getWarpingMaps()
    {
        return warpingMaps;
    }

    public void plotWarpings()
    {
        if (warpingMaps == null)
            return;
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (int i=0; i< warpingMaps.size(); i++)
        {
            XYSeries series = new XYSeries("run " + (i+1) + (i==0? " (reference)" : ""));
            int[] setScanMap = warpingMaps.get(i);
            for (int j=0; j<setScanMap.length; j++)
                series.add(j,setScanMap[j]);
            dataset.addSeries(series);
        }
        JFreeChart chart = ChartFactory.createXYLineChart(null, null, null, dataset,
                PlotOrientation.VERTICAL, true, false, false);

        ChartDialog chartDialog = new ChartDialog(chart);
        chartDialog.setSize(800,600);

        chartDialog.setVisible(true);
    }

    public void plotWarpings(double[][] warpings)
    {
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (double[] warping : warpings)
        {
            XYSeries series = new XYSeries("");
            for (int j=0; j<warping.length; j++)
                series.add(j,warping[j]);
            dataset.addSeries(series);
        }
        JFreeChart chart = ChartFactory.createXYLineChart(null, null, null, dataset,
                PlotOrientation.VERTICAL, true, false, false);

        ChartDialog chartDialog = new ChartDialog(chart);
        chartDialog.setSize(800,600);

        chartDialog.setVisible(true);
    }

    public void plotWarping(double[] warping)
    {
        double[][] warpings = new double[1][warping.length];
        warpings[0] = warping;
        plotWarpings(warpings);
    }






    /**
     * Interface for classes that provide a way of selecting pairs of features.
     * Each implementing class will define a different strategy for selecing pairs.
     */
    public static interface FeaturePairSelector
    {
        public abstract Pair<Feature,Feature>[] selectPairs(FeatureSet sourceFeatureSet,
                                                            FeatureSet destFeatureSet);
    }

    /**
     * Selecte pairs of features based on peptide agreement
     */
    public static class PeptideFeaturePairSelector implements FeaturePairSelector
    {
        public Pair<Feature,Feature>[] selectPairs(FeatureSet sourceFeatureSet,
                                                   FeatureSet destFeatureSet)
        {
            List<Pair<Feature,Feature>> resultList =
                new ArrayList<Pair<Feature,Feature>>();

            for (Feature sourceFeature : sourceFeatureSet.getFeatures())
            {
                String sourcePeptide = MS2ExtraInfoDef.getFirstPeptide(sourceFeature);
                if (sourcePeptide == null)
                    continue;

                for (Feature destFeature : destFeatureSet.getFeatures())
                {
                    String destPeptide = MS2ExtraInfoDef.getFirstPeptide(destFeature);
                    if (destPeptide != null && destPeptide.equals(sourcePeptide))
                       resultList.add(new Pair<Feature,Feature>(sourceFeature, destFeature));
                }
            }

            //have to cast, because Java is not good
            return (Pair<Feature,Feature>[]) resultList.toArray(new Pair[0]);
        }
    }


    /**
     * Select pairs of features based on peptide agreement, and also mass tolerance
     */
    public static class HybridFeaturePairSelector implements FeaturePairSelector
    {
        float massTolerance = 0.1f;

        public Pair<Feature,Feature>[] selectPairs(FeatureSet sourceFeatureSet,
                                                   FeatureSet destFeatureSet)
        {
            List<Pair<Feature,Feature>> resultList =
                new ArrayList<Pair<Feature,Feature>>();

            for (Feature sourceFeature : sourceFeatureSet.getFeatures())
            {
                String sourcePeptide = MS2ExtraInfoDef.getFirstPeptide(sourceFeature);
                if (sourcePeptide == null)
                    continue;

                for (Feature destFeature : destFeatureSet.getFeatures())
                {
                    String destPeptide = MS2ExtraInfoDef.getFirstPeptide(destFeature);
                    if (destPeptide != null && destPeptide.equals(sourcePeptide))
                       resultList.add(new Pair<Feature,Feature>(sourceFeature, destFeature));
                }
            }

            List<Pair<Feature,Feature>> forRemoval =
                    new ArrayList<Pair<Feature,Feature>>();

            for (Pair<Feature,Feature> resultPair : resultList)
            {
                if (Math.abs(resultPair.first.getMass() -
                             resultPair.second.getMass()) > massTolerance)
                    forRemoval.add(resultPair);
            }
            for (Pair<Feature,Feature> pairToRemove: forRemoval)
                resultList.remove(pairToRemove);
            _log.debug("Removed " + forRemoval.size() + " pairs out of mass tolerance");

            //have to cast, because Java is not good
            return (Pair<Feature,Feature>[]) resultList.toArray(new Pair[0]);
        }
    }

    public static class MassFeaturePairSelector extends MassOrMzFeaturePairSelector
        implements FeaturePairSelector
    {
        public static final float DEFAULT_DELTA_MASS = 10f;

        public MassFeaturePairSelector()
        {
            super();
            setMassOrMzMode(MODE_MASS);
            setDeltaMass(DEFAULT_DELTA_MASS);
        }

        public void setDeltaMass(float deltaMass)
        {
            setDeltaMassOrMz(deltaMass);
        }

        public float getDeltaMass()
        {
            return getDeltaMassOrMz();
        }
    }

    public static class MzFeaturePairSelector extends MassOrMzFeaturePairSelector
        implements FeaturePairSelector
    {
        public static final float DEFAULT_DELTA_MZ =0.1f;

        public MzFeaturePairSelector()
        {
            super();
            setMassOrMzMode(MODE_MZ);
            setDeltaMz(DEFAULT_DELTA_MZ);
        }

        public void setDeltaMz(float deltaMz)
        {
            setDeltaMassOrMz(deltaMz);
        }

        public float getDeltaMz()
        {
            return getDeltaMassOrMz();
        }
    }

    /**
     * This class selects pairs of features purely based on m/z -- it pairs up features
     * that are within a given m/z tolerance.  This emulates the old behavior, which was
     * implemented in R
     */
    public static class MassOrMzFeaturePairSelector implements FeaturePairSelector
    {
        public static final float DEFAULT_DELTA_MASS_OR_MZ = 10f;
        public static final float DEFAULT_MIN_INTENSITY=100;

        protected int topN = -1;
        protected float deltaMassOrMz = DEFAULT_DELTA_MASS_OR_MZ;
        protected float minIntensity = DEFAULT_MIN_INTENSITY;

        public static final int MODE_MASS = 0;
        public static final int MODE_MZ = 1;

        protected int massOrMzMode = MODE_MZ;

        protected int deltaMassType = DEFAULT_DELTA_MASS_TYPE;
        public static final int DEFAULT_DELTA_MASS_TYPE =
                FeatureSetMatcher.DELTA_MASS_TYPE_PPM;

        public void setTopN(int newTopN)
        {
            topN = newTopN;
        }

        public int getTopN()
        {
            return topN;
        }

        public float getDeltaMassOrMz()
        {
            return deltaMassOrMz;
        }

        public void setDeltaMassOrMz(float deltaMz)
        {
            this.deltaMassOrMz = deltaMz;
        }

        public float getMinIntensity()
        {
            return minIntensity;
        }

        public void setMinIntensity(float minIntensity)
        {
            this.minIntensity = minIntensity;
        }

        /**
         * This will return more than N entries when more than one entry shares the
         * _exact_ value of the boundry element; we do this to avoid making random exclusions
         * of "equivalent" values.
         *
         * Note that, in the degenerate case, this will return
         * all features if all members of the input list have the same value.
         */
        protected FeatureSet stripAllButTopNFeaturesByIntensity(FeatureSet featureSet)
        {
            Feature[] features = featureSet.getFeatures();
            Arrays.sort(features, new Feature.IntensityDescComparator());
            if (features.length == 0)
                ApplicationContext.infoMessage("ERROR! 0 features in this FeatureSet");
            float minIntensity = features[Math.min(features.length-1,topN-1)].getIntensity();

            FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
            sel.setMinIntensity(minIntensity);

            return featureSet.filter(sel);
        }

        /**
         *
         * Select pairs based on m/z or
         * @return
         */
        public Pair<Feature,Feature>[] selectPairs(FeatureSet sourceFeatureSet,
                                                   FeatureSet destFeatureSet)
        {
            _log.debug("selectPairs: topN=" + topN + ", mode (mass==0): " +
                    massOrMzMode + ", minIntensity: " + minIntensity + " sourceFeatures: " +
                    sourceFeatureSet.getFeatures().length + ", destFeatures: " +
                    destFeatureSet.getFeatures().length + ", deltaMassOrMz: " + deltaMassOrMz +
                    ", deltaMassType (ppm==1): " + deltaMassType);
            FeatureSet workingSourceFeatureSet = (FeatureSet) sourceFeatureSet.clone();
            FeatureSet workingDestFeatureSet = (FeatureSet) destFeatureSet.clone();

            if (topN <= 0)
            {
                FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
                sel.setMinIntensity(minIntensity);
                workingSourceFeatureSet = workingSourceFeatureSet.filter(sel);
                workingDestFeatureSet = workingDestFeatureSet.filter(sel);
            }
            else
            {
                workingSourceFeatureSet =
                        stripAllButTopNFeaturesByIntensity(workingSourceFeatureSet);
                workingDestFeatureSet =
                        stripAllButTopNFeaturesByIntensity(workingDestFeatureSet);
            }

            Feature[] sourceFeatures = workingSourceFeatureSet.getFeatures();
            Feature[] destFeatures = workingDestFeatureSet.getFeatures();

            Comparator<Feature> featureMassOrMzAscComparator = null;

            switch (massOrMzMode)
            {
                 case MODE_MASS:
                     featureMassOrMzAscComparator = new Feature.MassAscComparator();
                     break;
                 case MODE_MZ:
                     featureMassOrMzAscComparator = new Feature.MzAscComparator();
                     break;
            }

            Arrays.sort(sourceFeatures, featureMassOrMzAscComparator);
            Arrays.sort(destFeatures, featureMassOrMzAscComparator);

            List<Pair<Feature,Feature>> resultList =
                    new ArrayList<Pair<Feature,Feature>>();

            _log.debug("selectPairs: pairing sets of size " +
                       sourceFeatures.length + " and " + destFeatures.length);
            switch (massOrMzMode)
            {
                case MODE_MASS:
                    _log.debug("    deltaMass: " + deltaMassOrMz);
                    break;
                case MODE_MZ:
                    _log.debug("    deltaMz: " + deltaMassOrMz);
                    break;
            }


            int i2top=0;
            int i2;
            for (Feature sourceFeature : sourceFeatures)
            {
                float sourceMassOrMz = sourceFeature.getMz();
                switch (massOrMzMode)
                {
                    case MODE_MASS:
                        sourceMassOrMz = sourceFeature.getMass();
                        break;
                    case MODE_MZ:
                        sourceMassOrMz = sourceFeature.getMz();
                        break;
                }

                i2=i2top;

                while (i2<destFeatures.length)
                {
                    float destMassOrMz = destFeatures[i2].getMz();
                    float effectiveDeltaMassOrMz = deltaMassOrMz;

                    switch (massOrMzMode)
                    {
                        case MODE_MASS:
                            destMassOrMz = destFeatures[i2].getMass();
                            effectiveDeltaMassOrMz =
                                    MassUtilities.calculateAbsoluteDeltaMass(
                                            sourceMassOrMz, deltaMassOrMz,
                                            deltaMassType);
                            break;
                        case MODE_MZ:
                            destMassOrMz = destFeatures[i2].getMz();
                            effectiveDeltaMassOrMz = deltaMassOrMz;
                            break;
                    }

                    if (destMassOrMz < sourceMassOrMz - effectiveDeltaMassOrMz)
                    {
                        i2++;
                        i2top=i2;
                        continue;
                    }
                    if (destMassOrMz > sourceMassOrMz + effectiveDeltaMassOrMz)
                        break;
                    //this 'if' should not be necessary.  However, for some reason,
                    //there seems to be some rounding going on up above, and things close to
                    //the boundary don't always get assigned the right way.
                    //For some reason the 'if' below works.
                    if (Math.abs(destMassOrMz - sourceMassOrMz) <= effectiveDeltaMassOrMz)
                        resultList.add(new Pair<Feature,Feature>(sourceFeature,
                                                                 destFeatures[i2]));
                    i2++;
                }
            }


            //have to cast, because Java is not good
            return (Pair<Feature,Feature>[]) resultList.toArray(new Pair[0]);
        }







        public int getMassOrMzMode()
        {
            return massOrMzMode;
        }

        public void setMassOrMzMode(int massOrMzMode)
        {
            this.massOrMzMode = massOrMzMode;
        }

        public int getDeltaMassType()
        {
            return deltaMassType;
        }

        public void setDeltaMassType(int deltaMassType)
        {
            this.deltaMassType = deltaMassType;
        }
    }

    public FeaturePairSelector getFeaturePairSelector()
    {
        return featurePairSelector;
    }

    public void setFeaturePairSelector(FeaturePairSelector featurePairSelector)
    {
        this.featurePairSelector = featurePairSelector;
    }


    public boolean isBuildCharts()
    {
        return buildCharts;
    }

    public void setBuildCharts(boolean buildCharts)
    {
        this.buildCharts = buildCharts;
    }


    public PanelWithChart getAlignmentChart()
    {
        return alignmentChart;
    }

    public void setAlignmentChart(PanelWithChart alignmentChart)
    {
        this.alignmentChart = alignmentChart;
    }

    public double getMaxStudRes()
    {
        return maxStudRes;
    }

    public void setMaxStudRes(double maxStudRes)
    {
        this.maxStudRes = maxStudRes;
    }

    public double getMaxLeverageNumerator()
    {
        return maxLeverageNumerator;
    }

    public void setMaxLeverageNumerator(double maxLeverageNumerator)
    {
        this.maxLeverageNumerator = maxLeverageNumerator;
    }

    public int getAlignmentOrderMode()
    {
        return alignmentOrderMode;
    }

    public void setAlignmentOrderMode(int alignmentOrderMode)
    {
        this.alignmentOrderMode = alignmentOrderMode;
    }

    public int getAlignmentDatasource()
    {
        return alignmentDatasource;
    }

    public void setAlignmentDatasource(int alignmentDatasource)
    {
        this.alignmentDatasource = alignmentDatasource;
    }
}
