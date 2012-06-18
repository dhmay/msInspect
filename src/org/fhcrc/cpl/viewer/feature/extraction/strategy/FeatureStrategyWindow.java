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

package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.ArrayList;
import java.util.Collection;

/**
 * A set of feature finders that works on 2D subsets of the spectrum space
 * at a time, then stitches all of the features found in those windows
 * together.
 */
public abstract class FeatureStrategyWindow extends BaseFeatureStrategy
{
    static Logger _log = Logger.getLogger(FeatureStrategyWindow.class);

    //start and end scan numbers for this window, including a margin
    //around what the user specified
    protected int _startScanNum = 0;
	protected int _endScanNum = 0;

    //Keep all the scans around.  Careful not to alter once set.
    //This is done because the window strategy requires a fringe of
    //scans around the window specified by the user, if the user specifies
    //a subset
    protected Scan[] _scans = null;

    //scans to keep around the edge of the window, to make sure we don't lose
    //features at edges
    protected static int _WindowMargin = 64;

    //default window width
    public static final int DEFAULT_WINDOW_WIDTH = 256;

    //window size to be considered (unless there aren't enough
    //scans in the run
    protected int _windowWidth = DEFAULT_WINDOW_WIDTH;

    public FeatureStrategyWindow()
    {

    }

    /**
     * This is where subclasses do the real work
     * @param spectra
     * @param scanWindow
     * @return
     * @throws InterruptedException
     */
    protected abstract Collection<Feature> findPeptidesIn2DWindow(
            float[][] spectra, Scan[] scanWindow)
             throws InterruptedException;


    /**
     * @param run
     * @param startScanIndex
     * @param scanCount
     * @param maxCharge
     * @param range
     */
    public void init(MSRun run, int startScanIndex,
                     int scanCount, int maxCharge, FloatRange range, boolean plotStatistics)
    {
        super.init(run, startScanIndex, scanCount, maxCharge, range, plotStatistics);

        //Determine the scans to use, with a margin of scans around the edge
        //of what the user specified, if available
        int c2 = Math.max(_windowWidth, scanCount + 2 * _WindowMargin);
		startScanIndex = Math.max(0, startScanIndex - (c2 - scanCount) / 2);
		int scanMax = Math.min(startScanIndex + c2, run.getScanCount());
		int scanCountWithWindow = scanMax - startScanIndex;

        _startScanNum = run.getScan(_startScan).getNum();
        _endScanNum = run.getScan(_endScan).getNum();

        _scans = FeatureFinder.getScans(_run, startScanIndex,
                                        scanCountWithWindow);
    }


    /**
     * Find features in the window, then filter out ones that fall outside
     * @return
     * @throws InterruptedException
     */
    public Feature[] findPeptides() throws InterruptedException
    {
        Feature[] features =
                analyzeWindow(_scans, _windowWidth, _WindowMargin);
        List<Feature> filtered = new ArrayList<Feature>();
        for (Feature feature : features)
        {
            if (feature.scan >= _startScanNum &&
                feature.scan <= _endScanNum)
                filtered.add(feature);
        }
        return filtered.toArray(new Feature[filtered.size()]);
    }

    /**
     * Divide the scans into windows, find features in each window, and
     * stitch together.
     *
     * @param scans
     * @return
     */
    protected Feature[] analyzeWindow(Scan[] scans, int windowWidth, int windowMargin)
            throws InterruptedException
    {
        _log.debug("analyzeWindow " + scans[0].getNum() + "-" + scans[scans.length-1].getNum());
        Thread currentThread = Thread.currentThread();
        List<Feature> allFeatures = new ArrayList<Feature>();
        int scanNum = 0;
        int endWindowScan = 0;

        if (null != _status)
            _status.progress(0.0F);

        //main loop
        do
        {
            endWindowScan = Math.min(scans.length, scanNum+windowWidth);
            int startWindowScan = Math.max(0, endWindowScan - windowWidth);
            Scan[] currentScanWindow = new Scan[endWindowScan-startWindowScan];
            System.arraycopy(scans, startWindowScan, currentScanWindow,
                             0, endWindowScan-startWindowScan);

            //resample just the spectra in this window
            float[][] resampledSpectra =
                    resampleSpectra(currentScanWindow);

            //Do the actual work
            Collection<Feature> byScan =
                    findPeptidesIn2DWindow(resampledSpectra, currentScanWindow);

            //If we're dumping spectra, dump spectra, in resampled space
            if (getDumpWindowSize() > 0)
                dumpWindow(byScan, resampledSpectra);

            // fix up scan numbers:  translate from resampled space back
            // to regular space
            for (Feature f : byScan)
            {
                Scan featureScan = currentScanWindow[f.scan];
                f.setTime((float)featureScan.getDoubleRetentionTime());
                f.scan = featureScan.getNum();
                f.setScanFirst(currentScanWindow[f.getScanFirst()].getNum());
                f.setScanLast(currentScanWindow[f.getScanLast()].getNum());

                //dhmay fixing up comprised peak scan numbers
                if (f.comprised != null)
                    for (Spectrum.Peak peak : f.comprised)
                    {
                        if (peak != null)
                            peak.scan = _run.getScanNumForIndex(startWindowScan + peak.scan);
                    }
            }

            // only add features found within the window proper
            // (exclude margins)
            int s = (scanNum == 0) ? scanNum : scanNum + windowMargin;
            int e = (endWindowScan == scans.length) ?
                    endWindowScan : endWindowScan - windowMargin;
            int fromScanNum = scans[s].getNum();
            int toScanNum = scans[e-1].getNum();
            _log.debug("WINDOW  [" + scans[startWindowScan].getNum() + "-" +
                    scans[endWindowScan-1].getNum() + "] " + fromScanNum
                    + "-" + toScanNum + "*******");
            for (Feature feature : byScan)
            {
                if (feature.scan >= fromScanNum && feature.scan <= toScanNum)
                    allFeatures.add(feature);
            }

            //report progress. This is coarse-grained
            if (null != _status)
                _status.progress((endWindowScan-windowMargin)*100.0F/scans.length);
            
            if (currentThread.isInterrupted())
                throw new InterruptedException();

            scanNum += windowWidth - (2 * windowMargin);
        }
        while (endWindowScan < scans.length);



        if (null != _status)
            _status.progress(100.0F);

        return allFeatures.toArray(new Feature[allFeatures.size()]);
    }

    /**
     * Dump a window of spectra in the resampled space
     * @param features
     * @param resampledSpectra
     */
    protected void dumpWindow(Collection<Feature> features,
                              float[][] resampledSpectra)
    {
        //
        // Extract a window of intensities around each Feature
        //
        // this is may be used by downstream analsysis programs
        //
        int dumpWindowSize = getDumpWindowSize();
        if (dumpWindowSize > 0)
        {
            int nSamples = dumpWindowSize * SpectrumResampler.getResampleFrequency();
            for (Feature f : features)
            {
                int mzIndex = (int)( (f.mz - _mzRange.min) *
                                      SpectrumResampler.getResampleFrequency() );
                int scanIndex = f.scan;
                f.intensityLeadingPeaks = dumpWindowSize;
                f.intensityTrailingPeaks = dumpWindowSize;
                f.intensityWindow = new float[2 * nSamples];
                int k = 0;
                for (int j = (mzIndex - nSamples); j < (mzIndex + nSamples); j++, k++)
                {
                    if ( j < 0 || j >= resampledSpectra[scanIndex].length )
                        f.intensityWindow[k] = 0.f; // pad with zeros if we hit an edge
                    else
                        f.intensityWindow[k] = resampledSpectra[scanIndex][j];
                }
            }
        }
    }


    public int getWindowWidth()
    {
        return _windowWidth;
    }

    public void setWindowWidth(int windowWidth)
    {
        this._windowWidth = windowWidth;
    }
}
