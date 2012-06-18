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

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

/**
 * Extract peaks from spectra. Does NOT do background subtraction -- the
 * spectra passed to this method should already have background subtracted 
 */
public interface PeakExtractor
{
    public static final boolean DEFAULT_PEAK_RIDGE_WALK_SMOOTHED = false;    

    /**
     * Main peak-extraction method.  Does NOT do background subtraction -- the
     * spectra passed to this method should already have background subtracted.
     * @param scans
     * @param spectra
     * @param mzRange
     * @return
     * @throws InterruptedException
     */
    public Feature[] extractPeakFeatures(Scan[] scans, float[][] spectra,
                                             FloatRange mzRange)
            throws InterruptedException;
    public boolean isPeakRidgeWalkSmoothed();
    public void setPeakRidgeWalkSmoothed(boolean ridgeWalkSmoothed);
}
