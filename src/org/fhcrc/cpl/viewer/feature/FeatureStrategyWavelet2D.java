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

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.Smooth2D;
import org.fhcrc.cpl.viewer.feature.ExtractMaxima2D;
import org.fhcrc.cpl.viewer.feature.FeatureStrategyUsingWindow2D;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: Oct 12, 2004
 * Time: 9:23:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class FeatureStrategyWavelet2D extends FeatureStrategyUsingWindow2D
    {
    public FeatureStrategyWavelet2D(MSRun run, int scan, int count, int maxCharge, FloatRange range, double sn)
        {
        super(run, scan, count, maxCharge, range, sn);
        }


    public Spectrum.Peak[] ExtractPeaks(Scan[] scans, Spectrum.Peak[] grossFeatures)
        {
        return ExtractMaxima2D.analyze(scans, _mzRange, new SmoothWavelet(), -Float.MAX_VALUE);
        }


    public static class SmoothWavelet extends Smooth2D
        {
        Pair tmp = new Pair(null, null);

        public float[] SmoothSpectra(float[] spectrum)
            {
            float[] d3 = Spectrum.WaveletD3(spectrum, tmp);
            return d3;
            }

        public float[] SmoothElution(float[] elution)
            {
            return Spectrum.FFTsmooth(elution, smoothXfactor, false);
            }
        }
    }
