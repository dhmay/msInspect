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

import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;

/**
 * Removes background noise from spectra
 */
public class BackgroundRemover
{

    public float[][] removeBackground(float[][] spectra)
    {
        return Spectrum.RemoveBackground(spectra);
    }

    /**
     * Calculate median intensity at each point on the grid
     * @param spectra
     * @return
     */
    public float[][] calculateMedian(float[][] spectra)
    {
        int numSpectra = spectra.length;
        int spectrumHeight = spectra[0].length;

        float[][] median = new float[numSpectra][];
        for (int i = 0; i < numSpectra; i++)
            median[i] = Spectrum.MedianWindow(spectra[i], spectrumHeight, 2 * SpectrumResampler.getResampleFrequency(),
                    false);
        float[] row = null;
        for (int r = 0; r < spectrumHeight; r++)
        {
            row = Spectrum.getRow(spectra, r, row);
            float[] m = Spectrum.MedianWindow(row, numSpectra, SpectrumResampler.getResampleFrequency(), false);
            for (int s = 0; s < m.length; s++)
                median[s][r] = Math.max(median[s][r], m[s]);
        }
        return median;
    }

}
