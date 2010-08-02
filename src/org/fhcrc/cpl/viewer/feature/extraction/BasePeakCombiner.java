package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.viewer.feature.extraction.strategy.BaseFeatureStrategy;

/**
 * superclass for peak combiners
 */
public abstract class BasePeakCombiner implements PeakCombiner
{
    protected int _maxCharge = PeakCombiner.DEFAULT_MAX_CHARGE;


    public int getMaxCharge()
    {
        return _maxCharge;
    }

    public void setMaxCharge(int maxCharge)
    {
        _maxCharge = maxCharge;
    }
}
