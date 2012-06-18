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
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraction.*;
import org.fhcrc.cpl.viewer.feature.*;

import java.util.*;

/**
 * Base class for a modular feature strategy that contains well-defined components
 */
public class BaseFeatureStrategyModular extends FeatureStrategyWindow
{
    private static Logger _log = Logger.getLogger(BaseFeatureStrategyModular.class);

    protected static final int WINDOW_MARGIN = 64;

    protected CPUTimer timerAnalyze = new CPUTimer("BaseFeatureStrategyModular.analyze");
    protected CPUTimer timerExtractPeaks = new CPUTimer("BaseFeatureStrategyModular.peaks");
    protected CPUTimer timerExtractPeptides = new CPUTimer("BaseFeatureStrategyModular.peptides");

    //these two variables define the behavior of this feature strategy
    protected PeakExtractor peakExtractor = null;
    protected PeakCombiner peakCombiner = null;

    protected boolean peakRidgeWalkSmoothed =
            PeakExtractor.DEFAULT_PEAK_RIDGE_WALK_SMOOTHED;


    /**
     * THIS IS THE MAIN FEATURE FINDING ROUTINE
     *
     * Structure:
     *   Extract peaks -- wavelet decomposition, 
     *   FeatureStrategyUsingWindow.ExtractPeptideFeatures(), to tie features together
     *   Change scan numbers, which are currently indexes, to the actual scan numbers
     *   If centroided, call AccurateMassCentroid() to fix mass
     *
     */
    protected Collection<Feature> findPeptidesIn2DWindow(float[][] spectra, Scan[] scans)
            throws InterruptedException
    {
        if (peakExtractor == null || peakCombiner == null)
            throw new IllegalArgumentException("A peak extractor and peak combiner must both be specified for this strategy");

        WaveletPeakExtractor waveletPeakExtractor = new WaveletPeakExtractor();
        setPeakExtractor(waveletPeakExtractor);
        peakExtractor.setPeakRidgeWalkSmoothed(peakRidgeWalkSmoothed);

        Thread currentThread = Thread.currentThread();

        _log.debug("analyze2D " + scans[0].getNum() + "-" + scans[scans.length - 1].getNum());
        assert timerAnalyze.start();

        int numSpectra = spectra.length;
        int spectrumHeight = spectra[0].length;
        _log.debug("analyze2D datasize = " + (numSpectra * spectrumHeight * 4));

        if (currentThread.isInterrupted())
            throw new InterruptedException();

        // Extract peaks
        timerExtractPeaks.start();

        Feature[] peaks = peakExtractor.extractPeakFeatures(scans, spectra, _mzRange);
        assert timerExtractPeaks.stop();

        if (currentThread.isInterrupted())
            throw new InterruptedException();

        _log.debug("kept " + peaks.length + " peaks after filtering");

        // combine peaks into features representing peptides
        assert timerExtractPeptides.start();
        Arrays.sort(peaks, Spectrum.comparePeakMzAsc);
        Feature[] allPeptides = peakCombiner.createFeaturesFromPeaks(_run, peaks);
        assert timerExtractPeptides.stop();

        assert timerAnalyze.stop();
        CPUTimer.dumpAllTimers();
        
        List<Feature> result = new ArrayList<Feature>();

        for (Feature feature : allPeptides)
            result.add(feature);

        return result;
    }


    public PeakExtractor getPeakExtractor()
    {
        return peakExtractor;
    }

    public void setPeakExtractor(PeakExtractor peakExtractor)
    {
        this.peakExtractor = peakExtractor;
    }

    public PeakCombiner getPeakCombiner()
    {
        return peakCombiner;
    }

    public void setPeakCombiner(PeakCombiner peakCombiner)
    {
        this.peakCombiner = peakCombiner;
    }

    public boolean isPeakRidgeWalkSmoothed()
    {
        return peakRidgeWalkSmoothed;
    }

    public void setPeakRidgeWalkSmoothed(boolean peakRidgeWalkSmoothed)
    {
        this.peakRidgeWalkSmoothed = peakRidgeWalkSmoothed;
    }



}
