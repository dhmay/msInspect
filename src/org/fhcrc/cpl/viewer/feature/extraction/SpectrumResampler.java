package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * Resample spectra onto a grid with the specified frequency.
 * This is done for reasons of necessity:  we are looking for things that persist over time,
 * and  so we need to be able to compare things  across scans.
 *
 * The data are smoothed right after resampling, which can knock down the intensity of peaks
 * that appear right on the boundary of a resampled bin. 
 */
public class SpectrumResampler
{
    //We use a resampling rate of 1/36Da.  In theory this could be changed, but it
    //would require changes to the smoothing filters
    public static final int DEFAULT_RESAMPLE_FREQ = 36;

    protected int _resamplingFrequency = DEFAULT_RESAMPLE_FREQ;
	protected FloatRange _mzRange;
    protected boolean _useMedianSmooth = false;

    public SpectrumResampler(FloatRange mzRange)
    {
        setMzRange(mzRange);
    }

    public SpectrumResampler(FloatRange mzRange, int resamplingFrequency)
    {
        this(mzRange);
        setResamplingFrequency(resamplingFrequency);
    }

    /**
     * Resample the spectra, within the M/Z range _mzRange,
     * onto a regular grid, with frequency RESAMPLE_FREQ
     * @param scans
     * @return resampled spectra
     * @throws InterruptedException
     */
    public float[][] resampleSpectra(Scan[] scans)
            throws InterruptedException
    {
        Thread currentThread = Thread.currentThread();

        float[][] resampledSpectra = new float[scans.length][];
        for (int i = 0; i < scans.length; i++)
        {
            float[][] raw = scans[i].getSpectrum();
            if (currentThread.isInterrupted())
                throw new InterruptedException();
            resampledSpectra[i] =
                    Spectrum.Resample(raw, _mzRange, _resamplingFrequency);
        }
        int height = resampledSpectra[0].length;
        {
            float[] row = null, s = null;
            for (int i = 0; i < height; i++)
            {
                row = Spectrum.getRow(resampledSpectra, i, row);
                if (_useMedianSmooth) // this will remove lockspray
                { // use median
                    s = Spectrum.MedianSmooth(row, row.length, s);
                    Spectrum.SmoothALittle(s);
                    Spectrum.setRow(resampledSpectra, i, s);
                }
                else
                {
                    Spectrum.SmoothALittle(row);
                    Spectrum.setRow(resampledSpectra, i, row);
                }
            }
        }
        return resampledSpectra;
    }


    public int getResamplingFrequency() {
        return _resamplingFrequency;
    }

    public void setResamplingFrequency(int _resamplingFrequency) {
        this._resamplingFrequency = _resamplingFrequency;
    }

    public FloatRange getMzRange() {
        return _mzRange;
    }

    public void setMzRange(FloatRange _mzRange) {
        this._mzRange = _mzRange;
    }

    public boolean getUseMedianSmooth()
    {
        return _useMedianSmooth;
    }

    public void setUseMedianSmooth(boolean useMedianSmooth)
    {
        this._useMedianSmooth = useMedianSmooth;
    }
}
