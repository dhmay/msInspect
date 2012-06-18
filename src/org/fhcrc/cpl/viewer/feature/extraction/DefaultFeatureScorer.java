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

package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.BaseFeatureStrategy;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * Default feature scorer.
 * If we ever implement another one, maybe some machinery can be moved
 * up into a base class*
 *
 * Has side effects on scored features: adjusts m/z, and sets kl, etc.
 *
 */
public class DefaultFeatureScorer implements FeatureScorer
{

    //maximum peaks per feature to consider
    public static final int DEFAULT_MAX_PEAKS_PER_FEATURE = 10;
    protected int _maxPeaksPerFeature = DEFAULT_MAX_PEAKS_PER_FEATURE;

    private static Logger _log = Logger.getLogger(DefaultFeatureScorer.class);


    // need to pick units to scale by, 0 means good 1 is pretty bad...
    public static final float SUMSQUARES_MZ_WEIGHT = 6; // 1/6 of a Da is bad
    public static final float SUMSQUARES_SCALED_INTENSITY_WEIGHT = 4; // 1/4 in scaled intesity is bad

    protected boolean keepStatistics = false;
    protected List<Float> unscaledMzDistances = new ArrayList<Float>();
    protected List<Float> unscaledIntensityDistances = new ArrayList<Float>();


    /**
     * score a candidate feature f, and set the list of peaks that it's comprised of
     * <p/>
     * Expecting charge=0 features and peaks
     *
     * TODO: document this better
     *
     * @param f
     * @param peaks
     */
    public float scoreFeature(Feature f, Spectrum.Peak[] peaks)
    {

        double maxResampledDistanceBetweenFeatures =
                DefaultPeakCombiner.DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS /
                        (SpectrumResampler.getResampleFrequency() - 1);

        int absCharge = Math.abs(f.charge);
        float invAbsCharge = 1.0F / absCharge;
        float mass = 0f;
        if (f.charge >= 0)
            mass = (f.mz-Spectrum.HYDROGEN_ION_MASS) * absCharge;
        else mass = f.mz * absCharge;

        Spectrum.Peak[] peptidePeaks = new Spectrum.Peak[_maxPeaksPerFeature];

        // find start position in peak list
        int p = Arrays.binarySearch(peaks, f, Spectrum.comparePeakMzAsc);
        if (p < 0) p = -(p + 1);
        p = Math.max(0, p - 1);
        int pLastFound = p;

        float mzP0 = f.mz;
        float sum = 0.0F;
        float intensityLast = 0.0F;
        float intensityHighest = 0.0F;
        boolean skippedPeak = false;
        float distSum = 0;
        float distCount = 0;
        for (int Pi = 0; Pi < peptidePeaks.length; Pi++)
        {
            float mzPi = mzP0 + Pi * invAbsCharge + (distCount > 0 ? distSum/distCount : 0);
            p = findClosestPeak(peaks, mzPi, p);
            Spectrum.Peak peakFound = peaks[p];
            float dist = Math.abs(peakFound.mz - mzPi);

//QUESTION: why 2.5?
            boolean found =
                    !peakFound.excluded &&
                            Math.abs(dist) < 2.5 * maxResampledDistanceBetweenFeatures;

            if (!found)
            {
                // miss two in a row, call it quits
                if (Pi > 0 && null == peptidePeaks[Pi-1])
                    break;
            }
            else
            {
                intensityHighest = Math.max(peakFound.intensity, intensityHighest);

                // if peak is too BIG stop:
                // CONSIDER do this at end, use fn of signal v. expected?
                // don't like hard cut-offs
                if (Pi > 2 && peakFound.intensity != intensityHighest &&
                        peakFound.intensity > intensityLast * 1.33 + f.median)
                    break;

                // fine-tune anchor MZ
                if (peakFound.intensity > intensityHighest/2)
                {
                    distSum += dist;
                    distCount++;
                }

                for (int s = pLastFound + 1; s < p; s++)
                {
                    Spectrum.Peak skipped = peaks[s];
                    if (skipped.intensity > intensityLast / 2 + f.median)
                        skippedPeak = true;
                }

                // record it
                peptidePeaks[Pi] = peakFound;
                pLastFound = p;
                intensityLast = peakFound.intensity;
            }
        }

        // Adjust MZ for feature
        float mean=0;
        float weight=0;
        for (int i=0 ; i<peptidePeaks.length ; i++)
        {
            //float mi = m0 + i * invCharge;
            if (null == peptidePeaks[i]) continue;
            mean += (peptidePeaks[i].mz - i*invAbsCharge) * peptidePeaks[i].intensity;
            weight += peptidePeaks[i].intensity;
        }
        f.mz = mean/weight;

        //
        // Calculate quality measures
        //

        //The 'signal' array will hold values for each peak, scaled to sum to 1.
        //Missing peaks are filled in with intensity 0.1 prior to scaling.
        //Maximum of 6 peaks considered.
        //This block also determines which peaks are "counted" in the Feature.peaks count.
        //Only features with intensity >= 1/50th of the most-intense peaks, and greater than twice
        //the median background noise, are counted
        int signalLength = Math.min(6, _maxPeaksPerFeature);
        float[] signal = new float[signalLength];
        int countPeaks = 0;
        for (int i = 0; i < peptidePeaks.length; i++)
        {
            Spectrum.Peak peak = peptidePeaks[i];
            if (null != peak && (peak.intensity > intensityHighest/50 &&
                    peak.intensity > 2 * f.median))
                countPeaks++;
            if (i < signal.length)
            {
                float v = null == peak ? 0 : peptidePeaks[i].intensity;
                signal[i] = Math.max(0.1F, v);
                sum += signal[i];
            }
        }
        for (int i = 0; i < signal.length; i++)
            signal[i] /= sum;

        f.kl = Spectrum.KLPoissonDistance(mass, signal);
        //TODO: Do we need another score?
        //f.score = f.kl >= 0 ? f.intensity * (float)-Math.log(f.kl + Math.exp(-1)) : 0F;
        f.setPeaks(countPeaks);
        f.skippedPeaks = skippedPeak;
        f.comprised = peptidePeaks;
        f.mzPeak0 = peptidePeaks[0].mz;

        //
        // What other metrics would be useful
        //

        float sumDist = calcSumSquaresDistance2D(f, signal, peptidePeaks);
        
        return sumDist;
    }

    /**
     * Assign a score to a particular feature, given a scaled representation of
     * its peak intensities and the peaks themselves.   This is kind of a 2-dimensional
     * sum-squares distance measure.
     *
     * Lower is better
     * @param f
     * @param scaledIntensityDistribution
     * @param peptidePeaks
     * @return
     */
    public float calcSumSquaresDistance2D(Feature f,
                                          float[] scaledIntensityDistribution,
                                          Spectrum.Peak[] peptidePeaks)
    {
        //Guess mass.  don't actually set f.mass, since this may not be the charge this feature
        //ends up having
        float mass = 0f;
        if (f.charge >= 0)
            mass = (f.mz-Spectrum.HYDROGEN_ION_MASS) * f.charge;
        else mass = f.mz * -f.charge;

        //Get a scaled poisson distribution for this mass
        float[] expectedDistribution = Spectrum.Poisson(mass);

        //Get the index of the last peak actually present. If that's beyond the end of the expected distribution,
        //move it back
        int lastPeak = 0;
        for (int i=0 ; i<scaledIntensityDistribution.length ; i++)
            if (null != peptidePeaks[i])
                lastPeak = i;
        lastPeak = Math.min(lastPeak, expectedDistribution.length-1);

        float sumDist = 0;
        for (int i=0 ; i<=lastPeak ; i++)
        {
            float dist = 1; // what should I use here....

            //only consider peaks that are actually present.
            //TODO: does this give an unfair advantage to features with missing peaks?
            if (null != peptidePeaks[i])
            {
                //compare observed peaks with theoretical peaks based on Dalton offsets from the base feature m/z
                float theoreticalPeakMz = (f.mz + i / (float) Math.abs(f.charge));
                float distMZ = Math.abs(peptidePeaks[i].mz - theoreticalPeakMz);

                //Remove the effect of resampling error on m/z accuracy
                distMZ = Math.max(0, distMZ - 1/(2*SpectrumResampler.getResampleFrequency()));

                float distIn =
                        Math.abs(scaledIntensityDistribution[i] - expectedDistribution[i]);

                if (keepStatistics)
                {
                    unscaledMzDistances.add(distMZ);
                    unscaledIntensityDistances.add(distIn);
                }

                //scale the dimensions appropriately
                distMZ *= SUMSQUARES_MZ_WEIGHT;
                distIn *= SUMSQUARES_SCALED_INTENSITY_WEIGHT;
                
                dist = distMZ * distMZ + distIn * distIn;

                //Scale the contribution of this peak by the predicted scaled intensity of the peak
                sumDist += dist * expectedDistribution[i];

//System.err.println("\t\tmass: " + mass + ", distMZ: " + distMZ + ", distIn: " + distIn + ", first int: " + scaledIntensityDistribution[0] + ", first expected int: " + expectedDistribution[0] + ", sumdist: " + sumDist);

            }
        }

        if (_log.isDebugEnabled() && Float.isNaN(sumDist))
        {
            _log.debug("scorePeakDistribution: returning NaN.  mz = " + f.mz + ", charge=" + f.charge);
        }

        return sumDist;
    }

    protected int findClosestPeak(Spectrum.Peak[] peaks, float x, int start)
    {
        float dist = Math.abs(x - peaks[start].mz);
        int i = start + 1;
        for (; i < peaks.length; i++)
        {
            float d = Math.abs(x - peaks[i].mz);
            if (d > dist)
                break;
            dist = d;
        }
        return i - 1;
    }


    public int getMaxPeaksPerFeature()
    {
        return _maxPeaksPerFeature;
    }

    public void setMaxPeaksPerFeature(int maxPeaksPerFeature)
    {
        this._maxPeaksPerFeature = maxPeaksPerFeature;
    }


    public boolean isKeepStatistics()
    {
        return keepStatistics;
    }

    public void setKeepStatistics(boolean keepStatistics)
    {
        this.keepStatistics = keepStatistics;
    }

    public void plotStatistics()
    {
        if (!keepStatistics)
            return;
        PanelWithHistogram pwhMz = new PanelWithHistogram(unscaledMzDistances,"Unscaled M/Z Dist");
        pwhMz.displayInTab();
        PanelWithHistogram pwhInt = new PanelWithHistogram(unscaledIntensityDistances,"Unscaled Int Dist");
        pwhInt.displayInTab();
//        float[] mzDistArray = new float[unscaledMzDistances.size()];
//        float[] intDistArray = new float[unscaledMzDistances.size()];
//        for (int i=0; i<unscaledMzDistances.size(); i++)
//        {
//            mzDistArray[i] = unscaledMzDistances.get(i);
//            intDistArray[i] = unscaledIntensityDistances.get(i);
//        }

        PanelWithScatterPlot pwsp = new PanelWithScatterPlot(unscaledMzDistances, unscaledIntensityDistances,
                                   "mz vs int dist");
        pwsp.setPointSize(1);
        pwsp.setAxisLabels("m/z distance","intensity distance");
        pwsp.displayInTab();
    }
}
