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

import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.viewer.feature.ExtractMaxima2D;

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
	protected FloatRange _mzRange;
    protected boolean _useMedianSmooth = false;

    static final int DEFAULT_RESAMPLE_FREQUENCY = 36;

    public static int resampleFrequency = DEFAULT_RESAMPLE_FREQUENCY;
    public static float resampleInterval = 1.0f / resampleFrequency;

    public SpectrumResampler(FloatRange mzRange)
    {
        setMzRange(mzRange);
    }

    public SpectrumResampler(FloatRange mzRange, int resamplingFrequency)
    {
        this(mzRange);
        setResampleFrequency(resamplingFrequency);
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
                    Spectrum.Resample(raw, _mzRange, getResampleFrequency());
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

    public static int getResampleFrequency()
    {
        return resampleFrequency;
    }

    /**
     * Side effect: also sets resample interval to 1 / resampleFrequency.  resample interval not settable directly
     * @param resampleFrequencyNew
     */
    public static void setResampleFrequency(int resampleFrequencyNew)
    {
        resampleFrequency = resampleFrequencyNew;
        resampleInterval = 1.0f / resampleFrequency;

    }

    public static float getResampleInterval()
    {
        return resampleInterval;
    }
}
