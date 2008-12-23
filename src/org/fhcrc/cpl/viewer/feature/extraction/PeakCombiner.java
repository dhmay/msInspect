package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;


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
