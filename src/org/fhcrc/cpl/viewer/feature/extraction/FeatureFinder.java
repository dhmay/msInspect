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
package org.fhcrc.cpl.viewer.feature.extraction;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategy;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.BaseFeatureStrategy;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategyPeakClusters;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategyWindow;

import java.util.*;


/**
 * 20110917: adding scan window size parameters
 */
public class FeatureFinder 
{
    static Logger _log = Logger.getLogger(FeatureFinder.class);

    public static final int DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS = 3;

    public static final Class<? extends FeatureStrategy> DEFAULT_FEATURE_FINDING_CLASS =
            FeatureStrategyPeakClusters.class;    

    protected MSRun _run;

    protected int _maxCharge;
    protected FloatRange _mzRange;
    protected int _startScan;
    protected int _scanCount;
    int _accurateMassAdjustmentScans = 0; // Number of scans around a feature to consider for accurate mass adjustment

    protected FeatureStrategy _featureStrategy;

    protected CPUTimer timerMassAdjustment = new CPUTimer("FeatureFinder.massAdjustment");

    protected int _dumpWindowSize = 0;

    protected boolean _peakRidgeWalkSmoothed =
            PeakExtractor.DEFAULT_PEAK_RIDGE_WALK_SMOOTHED;

    protected boolean _plotStatistics = false;

    protected int _scanWindowSize = FeatureStrategyWindow.DEFAULT_WINDOW_WIDTH;


    public FeatureFinder()
    {
        
    }

    public FeatureFinder(MSRun run, int startScan, int scanCount, int maxCharge,
                             FloatRange mzRange, Class<? extends FeatureStrategy> featureStrategyClass,
                             boolean plotStatistics) {
        this(run, startScan, scanCount, maxCharge, mzRange, featureStrategyClass, plotStatistics, FeatureStrategyWindow.DEFAULT_WINDOW_WIDTH);
    }

    public FeatureFinder(MSRun run, int startScan, int scanCount, int maxCharge,
                             FloatRange mzRange, Class<? extends FeatureStrategy> featureStrategyClass,
                             boolean plotStatistics, int scanWindowSize)
    {
        this();
        init(run, startScan, scanCount, maxCharge, mzRange,
             featureStrategyClass, plotStatistics, scanWindowSize);
    }
    
    public void init(MSRun run, int startScan, int scanCount, int maxCharge,
                             FloatRange mzRange, Class<? extends FeatureStrategy> featureStrategyClass,
                             boolean plotStatistics, int scanWindowSize)
    {
        assert(mzRange.min < mzRange.max);
        _run = run;
        _startScan = startScan;
        _scanCount = scanCount;
        _maxCharge = maxCharge;
        _mzRange = mzRange;
        _plotStatistics = plotStatistics;

        _featureStrategy =
                BaseFeatureStrategy.getInstance(run, startScan,
                                                scanCount,
                                                maxCharge, mzRange,
                                                featureStrategyClass, _plotStatistics);
        _featureStrategy.setDumpWindowSize(_dumpWindowSize);
        if (_featureStrategy instanceof FeatureStrategyWindow)
            ((FeatureStrategyWindow) _featureStrategy).setWindowWidth(scanWindowSize);
    }

    /**
     * Return the mz range for the given scan. If computed lowMz and highMz were
     * supplied, use those; otherwise try startMz and endMz.
     */
    public static FloatRange getMzExtractionRange(MSRun.MSScan scan)
    {
        if (scan.getLowMz() >= 0 && scan.getHighMz() >= 0)
            return new FloatRange(scan.getLowMz() - 1, scan.getHighMz() + 1);
        return  new FloatRange(scan.getStartMz() - 1, scan.getEndMz() + 1);
    }

    /**
     * Return the mz range for the given run. First, check scan zero.
     * If range is no good, use the mzRange computed when the run was
     * read.
     */
    public static FloatRange getMzExtractionRange(MSRun run)
    {
        FloatRange range = getMzExtractionRange(run.getScan(0));
        if (range.min >= 0 && range.max >= 0)
            return range;
        return new FloatRange(run.getMzRange().min - 1, run.getMzRange().max + 1);
    }









//    protected float[][] CombineScans(Scan[] scans, FloatRange range, int resample_freq)
//    {
//        float[][] zeroChargeSpectrum = null;
//        int countScans = 0;
//        for (int s = 0; s < scans.length; s++)
//        {
//            if (null == scans[s]) continue;
//            countScans++;
//            float[][] spectrumRaw = scans[s].getSpectrum();
//            float[][] t = Spectrum.ResampleSpectrum(spectrumRaw, range, resample_freq, true);
//            if (null == zeroChargeSpectrum)
//                zeroChargeSpectrum = t;
//            else
//            {
//                int len = t[1].length;
//                for (int i = 0; i < len; i++)
//                    zeroChargeSpectrum[1][i] += t[1][i];
//            }
//        }
//        int len = zeroChargeSpectrum[1].length;
//        if (countScans > 1)
//            for (int i = 0; i < len; i++)
//                zeroChargeSpectrum[1][i] /= countScans;
//        return zeroChargeSpectrum;
//    }




    // implement this in order to use the analyzeScanAtATime() helper
    protected Collection<Feature> analyze1D(Scan scan) throws InterruptedException
    {
        throw new java.lang.UnsupportedOperationException();
    }


    // implement this in order to use the analyzeWindow() helper
    protected Collection<Feature> analyze2D(Scan[] scans) throws InterruptedException
    {
        throw new java.lang.UnsupportedOperationException();
    }


    protected static void _logDebug(String s)
    {
        _log.debug(s);
    }



    public FeatureSet findPeptides() throws InterruptedException
    {
        Feature[] features = _featureStrategy.findPeptides();

        FeatureSet featureSet = new FeatureSet(features);


        assert timerMassAdjustment.start();
        // if data are centroided, or if requested explicitly for profile
        // mode data, attempt to get accurate masses
        if (getAccurateMassAdjustmentScans() > 0 ||
            _run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
        {
            //dhmay changing 20100316, to allow default-setting for the mass adjustment within the strategy            
            AccurateMassAdjuster massAdjuster = _featureStrategy.getAccurateMassAdjuster();
            if (massAdjuster == null)
                massAdjuster = new AccurateMassAdjuster();
            massAdjuster.setScanWindowSize(getAccurateMassAdjustmentScans());
            massAdjuster.adjustAllMasses(_run, features);
        }
        assert timerMassAdjustment.stop();

        addInfoToFeatureSet(featureSet);
        
        Arrays.sort(features, new Feature.IntensityDescComparator());
        return featureSet;
    }

    protected void addInfoToFeatureSet(FeatureSet featureSet)
    {
        String revision = (String)ApplicationContext.getProperty("REVISION");
        if (null == revision) revision = "";
        if (revision.toLowerCase().startsWith("revision:"))
            revision = revision.substring("Revision:".length()).trim();

        Map<String, Object> m = featureSet.getProperties();
        m.put("revision", revision);
        m.put("algorithm", this.getClass().getName());
        m.put("java.vendor", System.getProperty("java.vendor"));
        m.put("java.version", System.getProperty("java.version"));
        m.put("user.name", System.getProperty("user.name"));
        m.put("date", (new Date()).toString());
    }



    /**
     * Utility method
     * @param run
     * @param start
     * @param count
     * @return
     */
    public static Scan[] getScans(MSRun run, int start, int count)
    {
        //eh?
        boolean bSkipLockSpray = false;

        start = Math.max(0, start);
        int maxScanIndex = run.getScanCount();
        List<Scan> scanList = new ArrayList<Scan>();
        for (int i = start; i < maxScanIndex && scanList.size() < count; i++)
        {
            Scan scan = run.getScan(i);
            if (bSkipLockSpray && 1 == (scan.getNum() % 20))
                continue;
            scanList.add(scan);
        }
        return scanList.toArray(new Scan[scanList.size()]);
    }

    public void setStatusListener(BaseFeatureStrategy.StatusListener status)
    {
        _featureStrategy.setStatusListener(status);
    }

    public void setAccurateMassAdjustmentScans(int accurateMassAdjustmentScans)
    {
        _accurateMassAdjustmentScans = accurateMassAdjustmentScans;
    }

    public int getAccurateMassAdjustmentScans()
    {
        return _accurateMassAdjustmentScans;
    }

    public void setDumpWindowSize(int dumpWindowSize)
    {
        _dumpWindowSize = dumpWindowSize;
        if (_featureStrategy != null)
            _featureStrategy.setDumpWindowSize(dumpWindowSize);
    }

    public int getDumpWindowSize()
    {
        return _dumpWindowSize;
    }


    public boolean isPeakRidgeWalkSmoothed()
    {
        return _peakRidgeWalkSmoothed;
    }

    public void setPeakRidgeWalkSmoothed(boolean peakRidgeWalkSmoothed)
    {
        _peakRidgeWalkSmoothed = peakRidgeWalkSmoothed;
        if (_featureStrategy != null)
            _featureStrategy.setPeakRidgeWalkSmoothed(peakRidgeWalkSmoothed);
    }

    public void plotStatistics()
    {
        if (_plotStatistics && _featureStrategy != null)
            _featureStrategy.plotStatistics();            
    }
}
