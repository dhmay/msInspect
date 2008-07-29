package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.ExtractMaxima2D;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraction.DefaultPeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.PeakExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.SmootherCreator;
import org.fhcrc.cpl.viewer.feature.extraction.WaveletPeakExtractor;
import org.labkey.common.tools.FloatRange;
import org.labkey.common.tools.Scan;

public class FeatureStrategyMSe extends FeatureStrategyWindow {

    protected boolean peakRidgeWalkSmoothed =
        PeakExtractor.DEFAULT_PEAK_RIDGE_WALK_SMOOTHED;
    
    protected int msLevel = 1;
    
	protected WaveletPeakExtractor extractor = new WaveletPeakExtractor();
	protected DefaultPeakCombiner peakCombiner = new DefaultPeakCombiner();


    public FeatureStrategyMSe() {
    }
    
    public FeatureStrategyMSe(int msLevel) {
    	this();
    	this.msLevel = msLevel;
    }
    
    /**
     * @param run
     * @param startScanIndex
     * @param scanCount
     * @param maxCharge
     * @param range
     */
    @Override
    public void init(MSRun run, int startScanIndex,
                     int scanCount, int maxCharge, FloatRange range, boolean plotStatistics)
    {
    	super.init(run, startScanIndex, scanCount, maxCharge, range, plotStatistics);
    	if (msLevel == 2) {
    		_scans = _run.getMS2Scans();
    		if (_scans.length > scanCount) {
    			int endScanIndex = Math.max(startScanIndex + scanCount, _scans.length);
    			_scans = Arrays.copyOfRange(_scans, startScanIndex, endScanIndex);
    		}
    	}
    	extractor.setPeakRidgeWalkSmoothed(peakRidgeWalkSmoothed);
    	extractor.setShortPeak(2);
    	extractor.setUseIntensityCutoff(false);
        peakCombiner.setMaxCharge(_maxCharge);
        peakCombiner.setResamplingFrequency(_resamplingFrequency);
    }

    @Override
	protected Collection<Feature> findPeptidesIn2DWindow(float[][] spectra,
			Scan[] scanWindow) throws InterruptedException {
    	Feature[] allPeptides = extractor.extractPeakFeatures(_scans, spectra, _mzRange);
//		Spectrum.Peak[] rawPeaks = 
//            ExtractMaxima2D.analyze(spectra, _mzRange.min, 1 / ((float) _resamplingFrequency), SmootherCreator.getThresholdSmoother(), 0.0F);
//		
		Arrays.sort(allPeptides, Spectrum.comparePeakMzAsc);        
        allPeptides = peakCombiner.createFeaturesFromPeaks(_run, allPeptides);

		java.util.List<Feature> features = new ArrayList<Feature>(); 
//		for (Spectrum.Peak rawPeak : rawPeaks) {
//			features.add(new Feature(rawPeak));
//		}
		for (Feature peptide : allPeptides) {
			features.add(peptide);
		}
		return features;
	}

	@Override
	public boolean isPeakRidgeWalkSmoothed() {
		return peakRidgeWalkSmoothed;
	}

	@Override
	public void setPeakRidgeWalkSmoothed(boolean peakRidgeWalkSmoothed) {
		this.peakRidgeWalkSmoothed = peakRidgeWalkSmoothed;
	}

}
