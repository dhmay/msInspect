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
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.extraction.*;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * Feature-finding strategy for small molecules (mass < ~700), e.g., metabolites
 * Based on FeatureStrategyPeakClusters.  The biggest difference is SmallMoleculePeakCombiner.
 * 
 */
public class FeatureStrategySmallMolecule extends BaseFeatureStrategyModular
{
    private static Logger _log = Logger.getLogger(FeatureStrategySmallMolecule.class);


    public static final double DEFAULT_MAX_KL = 5;
    public static final int DEFAULT_MIN_PEAKS = 1;
    public static final int DEFAULT_MAX_MASS = 1000;

    //todo: infer ppm error in some way
    //This is the amount of PPM error allowed in each direction around each feature, when recalculating
    //intensity based on accurate mass.  Features for which no value can be found within tolerance are
    //removed.
    protected float intensityRecalcPPM = 2;



    protected AccurateMassAdjuster accMassAdjuster;

    public FeatureStrategySmallMolecule()
    {
    }

    public FeatureSet.FeatureSelector getDefaultFeatureSelector()
    {
        defaultFeatureSelector = super.getDefaultFeatureSelector();
        defaultFeatureSelector.setMinPeaks(DEFAULT_MIN_PEAKS);
        defaultFeatureSelector.setMaxKL(DEFAULT_MAX_KL);
        defaultFeatureSelector.setMaxMass(DEFAULT_MAX_MASS);
        defaultFeatureSelector.setAccurateMzOnly(true);
        return defaultFeatureSelector;
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

        WaveletPeakExtractor waveletPeakExtractor = new WaveletPeakExtractor();
        setPeakExtractor(waveletPeakExtractor);

        SmallMoleculePeakCombiner peakCombiner = new SmallMoleculePeakCombiner();
        peakCombiner.setKeepStatistics(plotStatistics);

        peakCombiner.setMaxCharge(_maxCharge);
        setPeakCombiner(peakCombiner);

        accMassAdjuster = new AccurateMassAdjuster();
        accMassAdjuster.setProfileMassMode(AccurateMassAdjuster.PROFILE_MASS_MODE_MAX);
        //todo: try to come up with some justification for this number
        accMassAdjuster.setResamplingSizeProportion(.8);
        accMassAdjuster.setShouldAdjustComprisedMasses(true);
        accMassAdjuster.setShouldRecalculateIntensities(true);
        accMassAdjuster.setIntensityRecalcPPMError(intensityRecalcPPM);

        //Can't control this here, at the moment
        //todo: control this default in the FeatureStrategies 
//        accMassAdjuster.setScanWindowSize(3);
    }

    public void plotStatistics()
    {
        if (_keepStatistics)
            ((DefaultPeakCombiner) peakCombiner).plotStatistics();
    }

    public AccurateMassAdjuster getAccurateMassAdjuster()
    {
        return accMassAdjuster;
    }

}
