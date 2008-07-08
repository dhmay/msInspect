package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.viewer.feature.Feature;
import org.labkey.common.tools.FloatRange;
import org.labkey.common.tools.Scan;

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
