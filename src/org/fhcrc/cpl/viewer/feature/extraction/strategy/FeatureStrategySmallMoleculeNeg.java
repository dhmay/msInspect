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
package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.extraction.*;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * Feature-finding strategy for small molecules (mass < ~700), e.g., metabolites
 * Based on FeatureStrategyPeakClusters.  The only difference is SmallMoleculePeakCombiner.
 *
 * This strategy also assumes all ions are negatively charged, and assigns masses and charges accordingly.
 */
public class FeatureStrategySmallMoleculeNeg extends FeatureStrategySmallMolecule
{
    private static Logger _log = Logger.getLogger(FeatureStrategySmallMoleculeNeg.class);


    public FeatureStrategySmallMoleculeNeg()
    {
    }


    /**
     * Plug in the peak combiner and peak extractor
     * @param run
     * @param startScan
     * @param scanCount
     * @param maxCharge
     * @param mzRange
     */
    public void init(MSRun run, int startScan,
                     int scanCount, int maxCharge,
                     FloatRange mzRange, boolean plotStatistics)
    {
        super.init(run, startScan, scanCount, maxCharge, mzRange, plotStatistics);
//
//        WaveletPeakExtractor waveletPeakExtractor = new WaveletPeakExtractor();
//        setPeakExtractor(waveletPeakExtractor);
//
//        SmallMoleculePeakCombiner peakCombiner = new SmallMoleculePeakCombiner();
//        peakCombiner.setKeepStatistics(plotStatistics);
//
//        peakCombiner.setMaxCharge(_maxCharge);
//        peakCombiner.setResamplingFrequency(_resamplingFrequency);
//        setPeakCombiner(peakCombiner);

        ((SmallMoleculePeakCombiner) peakCombiner).setNegativeChargeMode(true);             
    }
}
