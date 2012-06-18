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
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.FeatureStrategyUsingWindow;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * User: mbellew
 * Date: Oct 6, 2004
 * Time: 9:30:21 PM
 */
public class FeatureStrategyWavelet extends FeatureStrategyUsingWindow
    {
    public FeatureStrategyWavelet(MSRun run, int scan, int count, int maxCharge, FloatRange range, double sn)
        {
        super(run, scan, count, maxCharge, range, sn);
        }


    Pair tmp = new Pair(null, null);
    protected float[] _smooth(float[] s)
        {
        return Spectrum.WaveletD3(s, tmp);
        }
    }
