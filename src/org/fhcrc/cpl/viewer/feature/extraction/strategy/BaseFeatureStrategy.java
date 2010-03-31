package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.viewer.feature.extraction.AccurateMassAdjuster;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

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
    protected static int _resamplingFrequency =
            FeatureStrategy.DEFAULT_RESAMPLING_FREQUENCY;
    protected CPUTimer timerResample = new CPUTimer("BaseFeatureStrategy.Resample");

    boolean useMedianSmooth = false;

    protected boolean _keepStatistics = false;

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
        try
        {
            Constructor<? extends FeatureStrategy> cons = featureStrategyClass.getConstructor();
            FeatureStrategy featureStrategy = cons.newInstance();
            featureStrategy.init(run, startScan,
                                 scanCount, maxCharge, mzRange, plotStatistics);
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
        spectrumResampler.setResamplingFrequency(_resamplingFrequency);
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





    public void setDumpWindowSize(int dumpWindowSize)
    {
        _dumpWindowSize = dumpWindowSize;
    }

    public int getDumpWindowSize()
    {
        return _dumpWindowSize;
    }    


    public void setResamplingFrequency(int resamplingFrequency)
    {
        _resamplingFrequency = resamplingFrequency;
    }

    public int getResamplingFrequency()
    {
        return _resamplingFrequency;
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
