package org.fhcrc.cpl.viewer.feature.extraction.strategy;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.viewer.feature.extraction.AccurateMassAdjuster;

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
     * @param startScanIndex
     * @param scanCount
     * @param maxCharge
     * @param mzRange
     */
    public void init(MSRun run, int startScanIndex,
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

    /**
     * Return an AccurateMassAdjuster appropriate for this feature strategy.  For different default settings
     * for different strategies.  Should return null if the strategy doesn't care about those defaults
     * @return
     */
    public AccurateMassAdjuster getAccurateMassAdjuster();

    /**
     * Get the default FeatureSelector for this FeatureStrategy.  Different default filtering is
     * appropriate for different strategies
     * @return
     */
    public FeatureSet.FeatureSelector getDefaultFeatureSelector();
}
