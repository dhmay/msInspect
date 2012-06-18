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
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;


/**
 * User: mbellew
 * Date: Sep 15, 2004
 * Time: 6:42:33 PM
 */
public class Smooth2D
    {
    static Logger _log = Logger.getLogger(Smooth2D.class);

    public static float smoothYfactor = 8.0F;
    public static float smoothXfactor = 12.0F;


    /**
     * scans.length should be 2^n for efficiency of some FFT implementations
     */
    public void smooth(float[][] spectra)
        {
        //
        // smooth elution profile
        //

        float[] elution = new float[spectra.length];
        for (int imz = 0; imz < spectra[0].length ; imz++)
            {
            // create elution profile array
            getRow(spectra, imz, elution);
            float[] smooth = SmoothElution(elution);
            if (null == smooth)
                {
                _log.error("smooth: null==smooth, isInterrupted = " + (Thread.currentThread().isInterrupted() ? "true" : "false"));
                return;
                }
            setRow(spectra, imz, smooth);
            }

        //
        // smooth spectrum
        //

        for (int i = 0; i < spectra.length; i++)
            {
            spectra[i] = SmoothSpectra(spectra[i]);
            }
        }


    protected void getRow(float[][] m, int r, float[] out)
        {
        for (int s = 0; s < m.length; s++)
            out[s] = m[s][r];
        }


    protected void setRow(float[][] m, int r, float[] in)
        {
        for (int s = 0; s < m.length; s++)
            m[s][r] = in[s];
        }


    protected void SubtractBackground(float[][] spectra)
        {
        }


    protected float[] SmoothSpectra(float[] spectrum)
        {
        return Spectrum.FFTsmooth(spectrum, smoothYfactor, false);
        }


    protected float[] SmoothElution(float[] elution)
        {
        return Spectrum.FFTsmooth(elution, smoothXfactor, false);
        }
    }
