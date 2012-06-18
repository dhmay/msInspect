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

import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.ExtractMaxima2D;
import org.fhcrc.cpl.viewer.feature.extraction.PeakExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.SmootherCreator;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * smooth, ExtractMaxima2D.analyze()
 *   Sort by decreasing intensity
 *   Throw out peaks with < 5 scans
 *   Break up what look like double-humped elution profiles by breaking peaks
 *   at valleys in the profile
 *   For each peak:
 *           -Compute the extents
 *                  -Walking down scans, make sure that the peak's intensity
 *                   is greater than the intensities of the two m/z's around it
 *                      -Make sure that intensity is above threshold
 *                      -If you start going back up in intensity, stop
 *                  -Same for end scan
 *                  -If too short (<5 scans), toss it out
 *            -"Integrate" to determine total intensity: for each scan, add a
 *             block of intensity corresponding to the intensity on that scan
 *             multiplied by half the amount of time between the adjoining scans
 *            -Create a feature to represent the individual peak
 *   Exclude all features initially
 *   Filter features: for each feature:
 *           find all neighbors within 5 scans and 1.1 m/z units
 *           if the m/z's are too close (within 1.5 / _resamplingFrequency, for some reason), no dice
 *           make sure they're lined up, scan-wise
 *           If not enough scans combined in the two features, no dice
 *           If intensities not high enough, no dice
 *           If we got here, unexclude both
 */
public class WaveletPeakExtractor implements PeakExtractor
{
    protected static Logger _log = Logger.getLogger(WaveletPeakExtractor.class);

    protected boolean _peakRidgeWalkSmoothed = DEFAULT_PEAK_RIDGE_WALK_SMOOTHED;

    public static final float DEFAULT_THRESHOLD_OFFSET_WAVELET = 1.0f;
    public static final float DEFAULT_THRESHOLD_OFFSET_SMOOTHED = 0f;

    public static final int DEFAULT_PEAK_LENGTH_REQUIREMENT = 5;

    protected float minRidgeProportionOfMax = .05f;


    //This constant and variable control which level of the wavelet is used.  3, the default, works well for
    //a resampling frequency of 36.
    public static final int DEFAULT_WAVELET_LEVEL = 3;
    protected int waveletLevel = DEFAULT_WAVELET_LEVEL;

    //minimum number of scans
    protected int minPeakScans = DEFAULT_PEAK_LENGTH_REQUIREMENT;

    public Feature[] extractPeakFeatures(Scan[] scans, float[][] spectra,
                                         FloatRange mzRange)
            throws InterruptedException
    {
        // debug counters
        int d_rawPeaks = 0;
        int d_tooShort = 0;
        int d_filtered = 0;

        Thread currentThread = Thread.currentThread();

        int numSpectra = spectra.length;
        int spectrumHeight = spectra[0].length;


        //
        // separate background and signal components
        // we're pretty conservative about what we call background,
        // we can use background and/or median later to filter peaks
        // further.
        //
        BackgroundRemover backgroundRemover = new BackgroundRemover();
        float[][] background = backgroundRemover.removeBackground(spectra);
        float[][] median = backgroundRemover.calculateMedian(spectra);



        float[][] wavelets = extractWavelets(spectra);

        //smooth spectra
        float[][] smoothedSpectra = new float[wavelets.length][];
        for (int s = 0; s < wavelets.length; s++)
            smoothedSpectra[s] = wavelets[s].clone();
        SmootherCreator.getThresholdSmoother().smooth(smoothedSpectra);

        //
        // extract all maxima from a time-smoothed version of the wavelet transformed data
        // we have a lot of versions of the data now, but it's helpful to keep it all around
        //
        Spectrum.Peak[] rawPeaks =
                ExtractMaxima2D.analyze(smoothedSpectra, mzRange.min, SpectrumResampler.getResampleInterval(),
                        null, 0.0F);

        int numRawPeaks = rawPeaks.length;

        // ???? if no raw peaks to process, bail with empty collection to avoid out
        // of bounds errors
        if ( 0 == numRawPeaks )
            return new Feature[0];

        Arrays.sort(rawPeaks, Spectrum.comparePeakIntensityDesc);
        if (_log.isDebugEnabled())
        {
            _log.debug("raw peaks: " + rawPeaks.length);
            if (rawPeaks.length > 0)
                for (int i = 0; i < 10; i++)
                    _log.debug("  " + (i) + ": " +
                            rawPeaks[rawPeaks.length * i / 100].getIntensity());
            if (rawPeaks.length > 0)
                for (int i = 1; i < 10; i++)
                    _log.debug("  " + (i * 10) + ": " +
                            rawPeaks[rawPeaks.length * i / 10].getIntensity());
            _log.debug(" 100: " + rawPeaks[rawPeaks.length - 1].getIntensity());
        }


        //
        // We probably have a LOT of peaks at this point.  95% are extremely small.
        // The size<shortPeak check in this loop will throw out most, and I
        // don't need to be too aggressive here
        //
        // NOTE: this loop also tries to handle "double humped elution profiles"
        // by breaking peaks at "valleys" in the elution profile.  This could perhaps
        // be improved upon with a second pass that specifically handles this case and
        // makes sure that in a region the splitting is consistent.
        //

        float[][] ridge;
        float thresholdOffset;

        if (_peakRidgeWalkSmoothed)
        {
            // Use smoothed peaks for ridge-walking
            ridge = smoothedSpectra;
            thresholdOffset = DEFAULT_THRESHOLD_OFFSET_SMOOTHED;
        }
        else
        {
            // Use wavelet peaks for ridge-walking
            ridge = wavelets;
            thresholdOffset = DEFAULT_THRESHOLD_OFFSET_WAVELET;
        }

        ArrayList<Feature> waveletPeaks = new ArrayList<Feature>();
        for (int i = 0, max = rawPeaks.length; i < max; i++)
        {
            Spectrum.Peak peak = rawPeaks[i];
            int imz = Math.round((peak.mz - mzRange.min) * SpectrumResampler.getResampleFrequency());

            if (imz < 2 || imz >= spectrumHeight - 2)
                continue;
            int scan = peak.scan;
            float peakIn = spectra[peak.scan][imz];
            float peakBg = background[peak.scan][imz];
            float peakMd = median[peak.scan][imz];

            // compute extent of this feature
            // this is complicated by the fact that we may have
            // double-humped peaks we need to separate
            int scanStart = scan;
            int scanEnd = scan;
            //todo: everything in the thresholdIn calculation is questionable and performs rather badly for metabolite data.
            //It's quite easy to get multiple features that overlap, and multiple features for the same m/z that
            //don't really seem to be separated by a valley.
            float thresholdIn = thresholdOffset + 0.5F * background[scan][imz] +
                    2 * median[scan][imz];
            thresholdIn = Math.max(peakIn * minRidgeProportionOfMax, thresholdIn);
            float minIn = Float.MAX_VALUE;
            int scanMinIn = scan;

            for (; scanStart > 0; scanStart--)
            {
                float spectraIn = spectra[scanStart][imz];
                // are we still on the ridge and > 5% of max
                if (spectraIn < thresholdIn ||
                        ridge[scanStart][imz] < ridge[scanStart][imz - 2] ||
                        ridge[scanStart][imz] < ridge[scanStart][imz + 2])
                {
                    scanStart++;
                    break;
                }
                // have we crossed a 'valley'
                float smoothIn = smoothedSpectra[scanStart][imz];
                if (smoothIn < minIn)
                {
                    minIn = spectraIn;
                    scanMinIn = scanStart;
                }
                else if (smoothIn > minIn * 1.1 + thresholdIn) // made up cutoff
                {
                    scanStart = scanMinIn;
                    break;
                }
            }

            minIn = Float.MAX_VALUE;
            scanMinIn = scan;
            for (; scanEnd < numSpectra - 1; scanEnd++)
            {
                float spectraIn = spectra[scanEnd][imz];
                // are we still on the ridge and > 5% of max
                if (spectraIn < thresholdIn ||
                        ridge[scanEnd][imz] < ridge[scanEnd][imz - 2] ||
                        ridge[scanEnd][imz] < ridge[scanEnd][imz + 2])
                {
                    scanEnd--;
                    break;
                }
                // have we crossed a 'valley'
                float smoothIn = smoothedSpectra[scanEnd][imz];
                if (smoothIn < minIn)
                {
                    minIn = spectraIn;
                    scanMinIn = scanEnd;
                }
                else if (smoothIn > minIn * 1.1 + thresholdIn) // made up cutoff
                {
                    scanEnd = scanMinIn;
                    break;
                }
            }


            if (scanEnd - scanStart + 1 < minPeakScans)
            {
                d_tooShort++;
                continue;
            }

            // total intensity calculation
            // "integrate" intensity over time
            float totalIn = 0;
            for (int s = scanStart; s <= scanEnd; s++)
            {
                double in = spectra[s][imz];
                double time = 1.0;
                if (s > 0 && s + 1 < scans.length)
                    time = (scans[s + 1].getDoubleRetentionTime() -
                            scans[s - 1].getDoubleRetentionTime()) /
                            2;
                totalIn += time * in;
            }

            Feature f = new Feature(peak.scan, scanStart, scanEnd,
                    peak.mz, peakIn, 0, -1, 0);
            f.background = peakBg;
            f.median = peakMd;
            f.totalIntensity = totalIn;
            waveletPeaks.add(f);
        }

        if (currentThread.isInterrupted())
            throw new InterruptedException();

        Collections.sort(waveletPeaks, Spectrum.compareFeatureLengthDesc);

        // some stats
        if (_log.isDebugEnabled())
        {
            _log.debug("FEATURE LENGTHS " + waveletPeaks.size() + " peaks");
            if (waveletPeaks.size() > 0)
            {
                for (int i = 0; i < 10; i++)
                    _log.debug("  " + (i) + ": " +
                            rawPeaks[waveletPeaks.size() * i / 100].getIntensity());
            }
            if (waveletPeaks.size() > 0)
            {
                for (int i = 1; i < 10; i++)
                    _log.debug("  " + (i * 10) + ": " +
                            rawPeaks[waveletPeaks.size() * i / 10].getIntensity());
            }
        }

        // OK, up the bar for a feature (min length of at largest peak in the feature)
        minPeakScans++;
        _log.debug("SHORT PEAK " + minPeakScans);

        //
        // Filtering
        //
        // Sticking with the idea of persistance, eliminate peaks
        // that don't seem to have correllated peaks
        //

        Tree2D tree = new Tree2D();
        for (Feature f : waveletPeaks)
        {
            f.excluded = true;
            tree.add(f.scan, f.mz, f);
        }
        Collections.sort(waveletPeaks, Spectrum.compareFeatureLengthDesc);
//QUESTION: so only features with other related features get included... so how do we
//end up with any single-peak features?

        for (Feature f : waveletPeaks)
        {
            int lenF = f.getScanLast() - f.getScanFirst() + 1;
            ArrayList<Feature> l = (ArrayList<Feature>)
                    tree.getPoints(f.scan - 5, f.mz - 1.1F, f.scan + 5, f.mz + 1.1F);
//too much detail
//            _log.debug("Evaluating feature with " + l.size() + " neighbors: " + f);
            for (Feature neighbor : l)
            {
                int lenN = neighbor.getScanLast() - neighbor.getScanFirst() + 1;
                // we want features of different masses
//QUESTION: what is the significance of 1.5?
                if (Math.abs(neighbor.mz - f.mz) < 1.5 / SpectrumResampler.getResampleFrequency())
                {
                    _log.debug("\tneighbor too close in m/z (" + Math.abs(neighbor.mz - f.mz) + ")");
                    continue;
                }
                // are centers aligned
                int s = Math.max(f.getScanFirst(), neighbor.getScanFirst());
                int e = Math.min(f.getScanLast(), neighbor.getScanLast());
                if (f.scan < s || f.scan > e || neighbor.scan < s || neighbor.scan > e)
                {
//too much detail
//                          _log.debug("\tneighbor scans not aligned");
                    continue;
                }
                // if neither is chosen yet, add a small hurdle
                if (f.excluded && neighbor.excluded)
                {
//QUESTION: why checking the total number of scans in both?
                    if (lenF + lenN < minPeakScans)
                    {
//too much detail
//                      _log.debug("\tneighbor combined scans too short");
                        continue;
                    }
                    if (f.intensity < 1 + 0.5 * f.background + 3 * Math.max(1, f.median) &&
                            neighbor.intensity < 1 + 0.5 * f.background +
                                    3 * Math.max(1, neighbor.median))
                    {
//too much detail
//                      _log.debug("\tneighbor too low-intensity");
                        continue;
                    }
                }
                f.excluded = false;
                neighbor.excluded = false;
//too much detail
//              _log.debug("\tSuccess!");
            }
        }
        ArrayList<Feature> copy = new ArrayList<Feature>();
        for (Feature f : waveletPeaks)
        {
            if (!f.excluded)
            {
                copy.add(f);
//too much detail
//              _log.debug("Final list: adding: " + f);
            }
            else
                d_filtered++;
        }
        waveletPeaks = copy;

        Feature[] result =
                waveletPeaks.toArray(new Feature[waveletPeaks.size()]);
        return result;
    }


    /**
     *
     *         the D3 matrix is used by default to locate peaks,
     *         the D3 processing effectively sharpens the raw data, this does a great job of cleaning
     *         up messy data, e.g. wide/overlapping peaks, and stray peaks in centroided data
     *
     *         why D3?  This is based on our default resampling frequency 1Da/36 and our max expected charge
     *         which I'm presuming to be 6.  Level three represents a feature width of 2^3=8, which
     *         is near to 36/6.  A big change in resampling or desired max charge, might require a
     *         change.
     *
     * dhmay changing 20100803 to use waveletLevel rather than always using 3
     * @param spectra
     * @return
     */
    public float[][] extractWavelets(float[][] spectra)
    {
        int numSpectra = spectra.length;

        float[][] wavelets = new float[numSpectra][];
        for (int s = 0; s < spectra.length; s++)
            wavelets[s] = Spectrum.WaveletDX(spectra[s], null, waveletLevel);
        return wavelets;
    }


    /**
     * Assign the extra necessary info to a set of features representing peaks
     * so that they can be written to a file
     * @param peakFeatures
     * @param scans
     */
    public void prettyPeakFeaturesForOutput(Collection<Feature> peakFeatures,
                                            Scan[] scans)
    {
        for (Feature peak : peakFeatures)
        {
            // fix up scan
            Scan scan = scans[peak.scan];
            peak.setTime((float)scan.getDoubleRetentionTime());
            peak.scan = scan.getNum();
            peak.scanFirst = scans[peak.scanFirst].getNum();
            peak.scanLast = scans[peak.scanLast].getNum();
            peak.setPeaks(1); // otherwise these will get filtered by default
        }
    }





    public void setShortPeak(int minPeakScans)
    {
        this.minPeakScans = minPeakScans;
    }

    public int getShortPeak()
    {
        return minPeakScans;
    }

    public boolean isPeakRidgeWalkSmoothed()
    {
        return _peakRidgeWalkSmoothed;
    }

    public void setPeakRidgeWalkSmoothed(boolean peakRidgeWalkSmoothed)
    {
        _peakRidgeWalkSmoothed = peakRidgeWalkSmoothed;
    }

    public int getWaveletLevel() {
        return waveletLevel;
    }

    public void setWaveletLevel(int waveletLevel) {
        this.waveletLevel = waveletLevel;
    }
}
