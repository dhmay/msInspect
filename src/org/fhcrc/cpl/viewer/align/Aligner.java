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
package org.fhcrc.cpl.viewer.align;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.toolbox.RegressionUtilities;
import org.fhcrc.cpl.viewer.amt.AmtFeatureSetMatcher;
import org.fhcrc.cpl.viewer.amt.AmtUtilities;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseMatcher;
import org.fhcrc.cpl.viewer.gui.util.ChartDialog;
import org.fhcrc.cpl.viewer.gui.util.PanelWithChart;
import org.fhcrc.cpl.viewer.gui.util.PanelWithScatterPlot;
import org.fhcrc.cpl.viewer.gui.util.PanelWithHistogram;
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

    public Aligner()
    {
    }

    /**
     * Calculate the maximum scan number of any of these features
     * @param features
     * @return
     */
    public static int getMaxScan(Feature[] features)
    {
        int maxScan = 0;
        for (Feature feature : features)
        {
            if (feature.getScan() > maxScan)
                maxScan = feature.getScan();
        }
        return maxScan;
    }

    /**
     * Calculate the maximum time of any of these features
     * @param features
     * @return
     */
    public static float getMaxTime(Feature[] features)
    {
        float maxTime = 0;
        for (Feature feature : features)
        {
            if (feature.getTime() > maxTime)
                maxTime = feature.getTime();
        }
        return maxTime;
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

        for (int i=1; i<sets.size(); i++)
        {
            Pair<Feature,Feature>[] pairedFeatures =
                    featurePairSelector.selectPairs(sets.get(i),sets.get(0));
            _log.debug("Aligner.alignFeatureSets: selected " + pairedFeatures.length +
                       " pairs");

            int maxScan = 0;
            for (Feature feature : sets.get(i).getFeatures())
            {
                if (feature.getScan() > maxScan)
                    maxScan = feature.getScan();
            }

            //necessary cast
            Pair<Integer,Double>[] pairedScans = (Pair<Integer,Double>[])
                    new Pair[pairedFeatures.length];
            for (int j=0; j<pairedScans.length; j++)
            {
                pairedScans[j] =
                        new Pair<Integer,Double>(pairedFeatures[j].first.getScan(),
                                                 (double) pairedFeatures[j].second.getScan());
            }

            Pair<Integer,Double>[] restrictedPairs = restrictPairsForRegression(pairedScans);
            _log.debug("Aligner.alignFeatureSets: restricted to " +
                    pairedFeatures.length + " pairs");

            String mapFileNameStart = (i+1) + "_onto_1";

            //this is where the work is done
            double[] scanMap = alignPairs(restrictedPairs, maxScan, mapFileNameStart);

            //It's lame to take up extra memory this way, but it's nice to be
            //able to persist this scan map for later investigation
            int[] intScanMap = new int[scanMap.length];
            for (int j=0; j<scanMap.length; j++)
                intScanMap[j] = (int) Math.round(scanMap[j]);
            warpingMaps.add(intScanMap);

            FeatureSet aligned = sets.get(i).deepCopy();

            //indicate that this featureset has been aligned by putting "align"
            //in the source filename.  OK, so this is a bit of a hack for backward
            //compatibility of the detail files
            aligned.setSourceFile(new File(getAlignedFileName(sets, i)));

            for (Feature feature : aligned.getFeatures())
            {
                feature.setScan(intScanMap[feature.getScan()]);
            }
            alignedSets.add(aligned);


            if (buildCharts)
            {
                createChart(pairedScans, restrictedPairs, scanMap);

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
                    PanelWithHistogram deltaMassHist = new PanelWithHistogram(deltaMassesPPM);
                    deltaMassHist.setName("massMatchDeltaPPM");
                    deltaMassHist.displayInTab();


                    PanelWithScatterPlot massMassDevPlot = new PanelWithScatterPlot(
                             massesPPM, deltaMassesPPM, "Mass (x) vs. deltaMass (y), ppm");
                    massMassDevPlot.displayInTab();
                }
            }

        }

        //clean up
        TempFileManager.deleteTempFiles(this);

        return alignedSets;
    }

    protected Pair<Integer,Double>[] restrictPairsForRegression(Pair<Integer,Double>[] allPairs)
    {
        Feature[] featuresForRegression = new Feature[allPairs.length];

        double[] baseTimes = new double[allPairs.length];
        double[] toAlignTimes = new double[allPairs.length];
        for (int j=0; j<allPairs.length; j++)
        {
            Pair<Integer,Double> matchedPair = allPairs[j];
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
            baseTimesForRegression[j] = featuresForRegression[j].getTime();
            toAlignTimesForRegression[j] =
                    AmtExtraInfoDef.getObservedHydrophobicity(featuresForRegression[j]);
        }

        featuresForRegression =
                RegressionUtilities.selectFeaturesWithLowLeverageAndStudentizedResidual(
                        featuresForRegression,
                        baseTimes, toAlignTimes,
                        maxLeverageNumerator,
                        maxStudRes,
                        false, 0, false);
        Pair<Integer,Double>[] result = (Pair<Integer,Double>[])
                new Pair[featuresForRegression.length];


        _log.debug(featuresForRegression.length + " pairs left after exclusion");

        for (int j=0; j<featuresForRegression.length; j++)
        {
            result[j] =  new Pair<Integer,Double>((int) featuresForRegression[j].getTime(),
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
    public double[] alignPairs(Pair<Integer,Double>[] pairs, int maxValueToWarp)
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
    public abstract double[] alignPairs(Pair<Integer,Double>[] pairs, int maxValueToWarp,
                               String tempFileNameStart);



    protected PanelWithChart createChart(Pair<Integer,Double>[] allPairedScans,
                                         Pair<Integer,Double>[] restrictedPairs,
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
    protected File writePairsFile(Pair<Integer,Double>[] pairs, String fileName)
            throws FileNotFoundException
    {
        File currentPairsFile =
                    TempFileManager.createTempFile(fileName,
                                                   this);

        PrintWriter currentPairsPW = new PrintWriter(currentPairsFile);
        currentPairsPW.println("source\tdest");
        for (Pair<Integer,Double> pair : pairs)
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
                AmtFeatureSetMatcher.DELTA_MASS_TYPE_PPM;

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

            float minIntensity = features[topN-1].getIntensity();

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
                                    AmtUtilities.calculateAbsoluteDeltaMass(
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
}
