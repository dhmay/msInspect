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

import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;

import java.util.ArrayList;

/**
 * User: mbellew
 * Date: Sep 15, 2004
 * Time: 6:42:33 PM
 *
 * 20100802: dhmay making resample frequency into a variable that can be modified
 */
public class ExtractMaxima2D
{
    static float smoothYfactor = 4.0F;
    static float smoothXfactor = 4.0F;

    private ExtractMaxima2D()
    {
    }


    public static Spectrum.Peak[] analyze(Scan[] scans, FloatRange rangePeaks)
    {
        return analyze(scans, rangePeaks, new Smooth2D(), 1);
    }


    /** find local maxima after smoothing using smooth2d,
     *
     *  peaks returned with scan numbers
     */
    public static Spectrum.Peak[] analyze(Scan[] scans, FloatRange rangePeaks, Smooth2D smooth2d, float minPeak)
    {
        // use larger ranger for resample to avoid edge checking later
        FloatRange range = new FloatRange(rangePeaks.min - .5F, rangePeaks.max + .5F);

        // allocate two extra scans to avoid edge checking later
        float[][] spectra = new float[scans.length + 2][];
        for (int i = 0; i < scans.length; i++)
            spectra[i + 1] = Spectrum.Resample(scans[i].getSpectrum(), range, SpectrumResampler.resampleFrequency);
        spectra[0] = spectra[1];
        spectra[spectra.length - 1] = spectra[spectra.length - 2];

        Spectrum.Peak[] peaks = analyze(spectra, range.min, SpectrumResampler.resampleInterval, smooth2d, minPeak);

        // fix up peaks to use scan num
        for (int p=0 ; p<peaks.length ; p++)
        {
            peaks[p].setScan(scans[peaks[p].scan-1].getNum());
        }
        return peaks;
    }


    /**
     * scans.length should be 2^n for efficiency
     *
     * NOTE: will not return maxima in first or last scan, or bottom/top mz value
     * NOTE: returns array index, not scan numbers
     */
    public static Spectrum.Peak[] analyze(float[][] spectra, float startMz, double interval, Smooth2D smooth2d, float minPeak)
    {
        //
        // smooth
        //

        if (null != smooth2d)
            smooth2d.smooth(spectra);

        //
        // find maxima
        //

        ArrayList peaks = new ArrayList();

        float[] s0 = null;
        float[] s1 = spectra[0];
        float[] s2 = spectra[1];
        int imzMax = spectra[0].length;
        for (int ispectra=1 ; ispectra<spectra.length-1 ; ispectra++)
        {
            s0 = s1;
            s1 = s2;
            s2 = spectra[ispectra+1];
            int[] scanPeaks = Spectrum.PickPeakIndexes(s1, minPeak);
            for (int p = 0 ; p<scanPeaks.length ; p++)
            {
                int imz = scanPeaks[p];
                float mz = (float)(startMz + imz * interval);
                if (imz < 1 || imz >= imzMax-1 )
                    continue;
                float intensity = s1[imz];
                assert intensity >= s1[imz-1] && intensity >= s1[imz+1];
                if (    intensity < s0[imz - 1] ||
                        intensity < s0[imz] ||
                        intensity < s0[imz + 1] ||
                        intensity < s2[imz - 1] ||
                        intensity < s2[imz] ||
                        intensity < s2[imz + 1])
                    continue;
                //int scanNum =scans[ispectra-1].getNum();
                peaks.add(new Spectrum.Peak(ispectra, mz, intensity));
            }
        }

        Spectrum.Peak[] arr = (Spectrum.Peak[]) peaks.toArray(new Spectrum.Peak[peaks.size()]);
        return arr;
    }

}


