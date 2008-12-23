package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.FloatRange;
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
