package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * Interface for all feature-finding strategies
 */
public interface FeatureStrategy
{
    public static final int DEFAULT_RESAMPLING_FREQUENCY = 36;

    /**
     * pass all of the interesting stuff to the FeatureStrategy so that
     * the findPeptides() method doesn't need any additional parameters
     * @param run
     * @param startScan
     * @param scanCount
     * @param maxCharge
     * @param mzRange
     */
    public void init(MSRun run, int startScan,
                     int scanCount, int maxCharge,
                     FloatRange mzRange, boolean plotStatistics);

    /**
     * Find peptides
     * @return
     * @throws InterruptedException
     */
    public Feature[] findPeptides() throws InterruptedException;

    public void setStatusListener(BaseFeatureStrategy.StatusListener status);

    public void setDumpWindowSize(int dumpWindowSize);

    public int getDumpWindowSize();

    public void setPeakRidgeWalkSmoothed(boolean peakRidgeWalkSmoothed);

    public boolean isPeakRidgeWalkSmoothed();

    public void plotStatistics();

    public void setKeepStatistics(boolean keepStatistics);

    public boolean isKeepStatistics();
}
