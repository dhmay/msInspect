package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.MSRun;


/**
 * Combine peaks into peptide features
 */
public interface PeakCombiner
{
    public static int DEFAULT_MAX_CHARGE = 6;

    /**
     * What it says
     * @param run
     * @param peaksIN
     * @return
     */
    public Feature[] createFeaturesFromPeaks(MSRun run, Spectrum.Peak[] peaksIN);       
}
