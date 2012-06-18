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

import org.fhcrc.cpl.viewer.feature.Smooth2D;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import modwt.Transform;
import modwt.Filter;

/**
 // Utility methods for creating different kinds of smoothers

 // So what's up with all these smoothers???
 //
 // Really this is a problem, the sample rate varies drastically
 // across experimental setups.  Most notably due to the number
 // of MS1 vs. MS2 scans, but also scan rate of the machine, and
 // the elution gradient.  I've settled on smoothTHRESHOLD for now.
 //
 // If the MS1 scan rate is too low, using this type of 2D
 // algorithm will not work well.
 //
 */
public class SmootherCreator
{
    public static Smooth2D getThresholdSmoother()
    {
        return new Smooth2D()
        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                return spectrum;
            }

            public float[] SmoothElution(float[] elution)
            {
                return _thresholdElution(elution);
            }
        };
    }

    /**
     * Not currently used
     * @return
     */
    public static Smooth2D getSmoothALot()
    {
        return new Smooth2D()
        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                //Spectrum.SmoothALittle(spectrum);
                return spectrum;
            }

            protected float[] SmoothElution(float[] elution)
            {
                return Spectrum.FFTsmooth(elution, 12, false);
            }
        };
    }

    /**
     * Not currently used
     * @return
     */
    public static Smooth2D getSmoothMedium()
    {
        return new Smooth2D()
        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                //Spectrum.SmoothALittle(spectrum);
                return spectrum;
            }

            protected float[] SmoothElution(float[] elution)
            {
                return Spectrum.FFTsmooth(elution, 6, false);
            }
        };
    }

    /**
     * Not currently used
     * @return
     */
    public static Smooth2D getSmoothD4()
    {
        return new Smooth2D()

        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                //Spectrum.SmoothALittle(spectrum);
                return spectrum;
            }

            protected float[] SmoothElution(float[] elution)
            {
                return Spectrum.WaveletD4(elution);
            }
        };
    }

    /**
     * Not currently used
     * @return
     */
    public static Smooth2D getSmoothALittle()
    {
        return new Smooth2D()
        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                //Spectrum.SmoothALittle(spectrum);
                return spectrum;
            }

            protected float[] SmoothElution(float[] elution)
            {
                Spectrum.SmoothALittle(elution);
                return elution;
            }
        };
    }

    /**
     * Not currently used
     * @return
     */    
    public static Smooth2D getSmoothNone()
    {
        return new Smooth2D()
        {
            protected float[] SmoothSpectra(float[] spectrum)
            {
                return spectrum;
            }

            protected float[] SmoothElution(float[] elution)
            {
                return elution;
            }
        };
    }

    /**
     * This actually gets used
     * @param elution
     * @return
     */
    public static float[] _thresholdElution(float[] elution)
    {
        // UNDONE: cache intermediate arrays
        int K = 3;
        int N = elution.length;
        Pair<float[][], float[][]> tmp = new Pair<float[][], float[][]>(null, null);

        tmp.first = Transform.decompose(elution, N, K, new Filter("haar"), "modwt",
                "periodic", tmp.first);
        tmp.second = Transform.multiresolution(tmp.first, N, K, new Filter("haar"), "modwt",
                "periodic", tmp.second);
        float[][] mra = tmp.second;
        float[] threshold;

        // even smoother
        threshold = mra[3];
        return threshold;
    }
}
