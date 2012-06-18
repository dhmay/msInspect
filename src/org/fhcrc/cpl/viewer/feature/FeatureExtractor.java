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
package org.fhcrc.cpl.viewer.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;

import java.lang.reflect.Constructor;
import java.util.*;


/**
 * User: mbellew
 * Date: Sep 7, 2004
 * Time: 12:08:54 PM
 */
public abstract class FeatureExtractor
{
    static Logger _log = Logger.getLogger(FeatureExtractor.class);
    public static final String DEFAULT_EXTRACTOR_PROPERTYNAME = FeatureExtractor.class + ".defaultExtractor";
    public static final int TYPE_1D = 1;    // process scans one at a time, cross scan features need to be combined
    public static final int TYPE_2D = 2;    // returns individual features in a window

    protected MSRun _run;
    protected int _maxCharge;
    protected FloatRange _mzRange;
    protected double _sn;
    protected StatusListener _status = null;
    protected int _startScan;
    protected int _scanCount;
    protected int _dumpWindowSize = 0; // size of the intensity window to grab around each feature
    int _accurateMassAdjustmentScans = 0; // Number of scans around a feature to consider for accurate mass adjustment

    static
    {
        setDefault(org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategyPeakClusters.class);
    }

    public interface StatusListener
    {
        void progress(float percent);
    }

    public static void setDefault(Class c)
    {
        ApplicationContext.setProperty(DEFAULT_EXTRACTOR_PROPERTYNAME, c);
    }

    public static Class getDefaultClass()
    {
        Class defaultClass =
                (Class) ApplicationContext.getProperty(DEFAULT_EXTRACTOR_PROPERTYNAME);
        return defaultClass;
    }

    /**
     * Instantiate the default FeatureExtractor
     * @param run
     * @param scan
     * @param count
     * @param maxCharge
     * @param range
     * @param sn
     * @return
     */
    public static FeatureExtractor getDefault(MSRun run, int scan, int count, int maxCharge, FloatRange range, double sn)
    {
        Class c = getDefaultClass();

        try
        {
            Constructor cons = c.getConstructor(
                    MSRun.class,
                    Integer.TYPE,
                    Integer.TYPE,
                    Integer.TYPE,
                    FloatRange.class,
                    Double.TYPE
            );
            return (FeatureExtractor)cons.newInstance(
                    run, new Integer(scan), new Integer(count), new Integer(maxCharge),
                    range, new Double(sn)
            );
        }
        catch (Exception x)
        {
            x.printStackTrace();
            System.exit(1);
            return null; // make compiler happy
        }
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

    protected FeatureExtractor(MSRun run, int startScan, int scanCount, int maxCharge, FloatRange mzRange, double sn)
    {
        assert(mzRange.min < mzRange.max);
        _run = run;
        _startScan = startScan;
        _scanCount = scanCount;
        _maxCharge = maxCharge;
        _mzRange = mzRange;
        _sn = sn;
    }

    public void setStatusListener(StatusListener status)
    {
        _status = status;
    }

    public void setDumpWindowSize(int dumpWindowSize)
    {
        _dumpWindowSize = dumpWindowSize;
    }

    protected int getDumpWindowSize()
    {
        return _dumpWindowSize;
    }

    public void setAccurateMassAdjustmentScans(int accurateMassAdjustmentScans)
    {
        _accurateMassAdjustmentScans = accurateMassAdjustmentScans;
    }

    protected int getAccurateMassAdjustmentScans()
    {
        return _accurateMassAdjustmentScans;
    }

    protected Scan[] getScans(MSRun run, int start, int count)
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


    protected float[][] CombineScans(Scan[] scans, FloatRange range, int resample_freq)
    {
        float[][] zeroChargeSpectrum = null;
        int countScans = 0;
        for (int s = 0; s < scans.length; s++)
        {
            if (null == scans[s]) continue;
            countScans++;
            float[][] spectrumRaw = scans[s].getSpectrum();
            float[][] t = Spectrum.ResampleSpectrum(spectrumRaw, range, resample_freq, true);
            if (null == zeroChargeSpectrum)
                zeroChargeSpectrum = t;
            else
            {
                int len = t[1].length;
                for (int i = 0; i < len; i++)
                    zeroChargeSpectrum[1][i] += t[1][i];
            }
        }
        int len = zeroChargeSpectrum[1].length;
        if (countScans > 1)
            for (int i = 0; i < len; i++)
                zeroChargeSpectrum[1][i] /= countScans;
        return zeroChargeSpectrum;
    }


    /**
     * Helper function for one scan at a time analyzer (TYPE_1D)
     *
     * @param scans
     * @return
     */
    protected Feature[] analyzeScanAtATime(Scan[] scans) throws InterruptedException
    {
        Thread currentThread = Thread.currentThread();
        ArrayList<Feature> all = new ArrayList<Feature>();

        if (null != _status)
            _status.progress(0.0F);

        for (int i = 0; i < scans.length; i++)
        {
            Collection<Feature> byScan = analyze1D(scans[i]);
            all.addAll(byScan);
            if (null != _status)
                _status.progress(i*100.0F/scans.length);
            if (currentThread.isInterrupted())
                throw new InterruptedException();
        }

        if (null != _status)
            _status.progress(100.0F);

        return all.toArray(new Feature[all.size()]);
    }


    /**
     * Helper function for window analyzer (TYPE_2D)
     *
     * @param scans
     * @return
     */
    protected Feature[] analyzeWindow(Scan[] scans, int windowWidth, int windowMargin)
            throws InterruptedException
    {
        _logDebug("analyzeWindow " + scans[0].getNum() + "-" + scans[scans.length-1].getNum());
        Thread currentThread = Thread.currentThread();
        List<Feature> all = new ArrayList<Feature>();
        int scan = 0;
        int end = 0;

        if (null != _status)
            _status.progress(0.0F);

        do
        {
            end = Math.min(scans.length, scan+windowWidth);
            int start = Math.max(0, end - windowWidth);
            Scan[] scanWindow = new Scan[end-start];
            System.arraycopy(scans, start, scanWindow, 0, end-start);

            Collection<Feature> byScan = analyze2D(scanWindow);
            // only add features found between the margins (exclude overlapping areas)
            int s = (scan == 0) ? scan : scan + windowMargin;
            int e = (end == scans.length) ? end : end - windowMargin;
            int fromScanNum = scans[s].getNum();
            int toScanNum = scans[e-1].getNum();
            _logDebug("WINDOW  [" + scans[start].getNum() + "-" + scans[end-1].getNum() + "] " + fromScanNum + "-" + toScanNum + "*******");
            for (Feature feature : byScan)
            {
                if (feature.scan >= fromScanNum && feature.scan <= toScanNum)
                    all.add(feature);
            }

            if (null != _status)
                _status.progress((end-windowMargin)*100.0F/scans.length);
            if (currentThread.isInterrupted())
                throw new InterruptedException();

            scan += windowWidth - 2 * windowMargin;
        }
        while (end < scans.length);

        if (null != _status)
            _status.progress(100.0F);

        return all.toArray(new Feature[all.size()]);
    }


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


    public int getType()
    {
        return TYPE_1D;
    }


    protected static void _logDebug(String s)
    {
        _log.debug(s);
    }


    protected abstract Feature[] _analyze() throws InterruptedException;


    public FeatureSet analyze() throws InterruptedException
    {
        Feature[] features = _analyze();
        FeatureSet fs = new FeatureSet(features);

        String revision = (String) ApplicationContext.getProperty("REVISION");
        if (null == revision) revision = "";
        if (revision.toLowerCase().startsWith("revision:"))
            revision = revision.substring("Revision:".length()).trim();

        Map m = fs.getProperties();
        m.put("revision", revision);
        m.put("algorithm", this.getClass().getName());
        m.put("java.vendor", System.getProperty("java.vendor"));
        m.put("java.version", System.getProperty("java.version"));
        m.put("user.name", System.getProperty("user.name"));
        m.put("date", (new Date()).toString());

        return fs;
    }
}
