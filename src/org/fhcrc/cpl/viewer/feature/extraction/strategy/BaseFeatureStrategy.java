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

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.viewer.feature.extraction.AccurateMassAdjuster;
import org.fhcrc.cpl.viewer.commandline.modules.FindPeptidesCommandLineModule;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;

import java.lang.reflect.Constructor;

/**
 *
 */
public abstract class BaseFeatureStrategy implements FeatureStrategy
{
    protected MSRun _run;
    protected int _maxCharge;
    protected FloatRange _mzRange;
    protected double _sn;
    //This is a scan INDEX, not a scan NUMBER.  Should really change the name, but worried about 3rd-party effects
    protected int _startScan;
    protected int _endScan;
    protected int _scanCount;
    protected StatusListener _status = null;

    protected int _dumpWindowSize = 0; // size of the intensity window to grab around each feature
    protected CPUTimer timerResample = new CPUTimer("BaseFeatureStrategy.Resample");

    boolean useMedianSmooth = false;

    protected boolean _keepStatistics = false;

    protected FeatureSet.FeatureSelector defaultFeatureSelector = null;

    public interface StatusListener
    {
        void progress(float percent);
    }

    public BaseFeatureStrategy()
    {
    }

    public void init(MSRun run, int startScanIndex,
                     int scanCount, int maxCharge,
                     FloatRange mzRange, boolean plotStatistics)
    {
        _run = run;
        _startScan = startScanIndex;
        _scanCount = scanCount;
        _endScan = Math.min(startScanIndex + scanCount, run.getScanCount() - 1);
        _maxCharge = maxCharge;
        _mzRange = mzRange;
        _keepStatistics = true;
    }

    public static FeatureStrategy getInstance(Class<? extends FeatureStrategy> featureStrategyClass)
    {
        try
        {
            Constructor<? extends FeatureStrategy> cons = featureStrategyClass.getConstructor();
            FeatureStrategy featureStrategy = cons.newInstance();
            return featureStrategy;
        }
        catch (Exception x)
        {
            x.printStackTrace();
            System.exit(1);
            return null; // make compiler happy
        }

    }

    /**
     * This is what we use to instantiate FeatureStrategy objects.
     * It's a bit crude. Could reimplement
     * @param run
     * @param startScan
     * @param scanCount
     * @param maxCharge
     * @param mzRange
     * @param featureStrategyClass
     * @return
     */
    public static FeatureStrategy getInstance(MSRun run,
                                              int startScan,
                                              int scanCount, int maxCharge,
                                              FloatRange mzRange,
                                              Class<? extends FeatureStrategy> featureStrategyClass,
                                              boolean plotStatistics)
    {
        FeatureStrategy featureStrategy = getInstance(featureStrategyClass);
        featureStrategy.init(run, startScan,
                scanCount, maxCharge, mzRange, plotStatistics);
        return featureStrategy;
    }

    /**
     * Utility method to resample spectra onto a standard grid
     * @param scans
     * @return
     * @throws InterruptedException
     */
    protected float[][] resampleSpectra(Scan[] scans)
            throws InterruptedException
    {
            //
        // Convert data into 2D matrix
        // we will do all processing on this data until the end and
        // then process back to "scan" space
        //
        assert timerResample.start();
        SpectrumResampler spectrumResampler = new SpectrumResampler(_mzRange);
        spectrumResampler.setUseMedianSmooth(useMedianSmooth);
        float[][] spectra =
                spectrumResampler.resampleSpectra(scans);
        assert timerResample.stop();
        return spectra;
    }

    public AccurateMassAdjuster getAccurateMassAdjuster()
    {
        return null;
    }

    public FeatureSet.FeatureSelector getDefaultFeatureSelector()
    {
        if (defaultFeatureSelector == null)
        {
            defaultFeatureSelector = new FeatureSet.FeatureSelector();
            defaultFeatureSelector.setMinPeaks(FindPeptidesCommandLineModule.DEFAULT_MIN_PEAKS);
            defaultFeatureSelector.setMaxKL(FindPeptidesCommandLineModule.DEFAULT_MAX_KL);
        }
        return defaultFeatureSelector;
    }






    public void setDumpWindowSize(int dumpWindowSize)
    {
        _dumpWindowSize = dumpWindowSize;
    }

    public int getDumpWindowSize()
    {
        return _dumpWindowSize;
    }    

    public void setStatusListener(BaseFeatureStrategy.StatusListener status)
    {
        _status = status;
    }

    public boolean isKeepStatistics()
    {
        return _keepStatistics;
    }

    public void setKeepStatistics(boolean keepStatistics)
    {
        this._keepStatistics = keepStatistics;
    }

    /**
     * Can be overridden
     */
    public void plotStatistics()
    {

    }
}
