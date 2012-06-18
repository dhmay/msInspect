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
package org.fhcrc.cpl.viewer.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

import java.util.*;

/**
 * User: mbellew
 * Date: Sep 7, 2004
 * Time: 12:10:06 PM
 */
public class FeatureStrategyUsingWindow extends FeatureExtractor
{
    static Logger _log = Logger.getLogger(FeatureStrategyUsingWindow.class);
    static Feature[] _emptyFeatureArray = new Feature[0];

    Scan[] _scans;
    protected double _noise = 1.0;

    static double debugMZ = 0;

    public FeatureStrategyUsingWindow(MSRun run, int scan, int count,
                                      int maxCharge, FloatRange range, double sn)
    {
        super(run, scan, count, maxCharge, range, sn);
        _scans = getScans(run, scan, count);
    }


    public Feature[] _analyze() throws InterruptedException
    {
        return analyzeScanAtATime(_scans);
    }


    protected float[] _smooth(float[] s)
    {
        return Spectrum.FFTsmooth(s, 6.0, false);
    }


    protected Collection<Feature> analyze1D(Scan scan)
    {
        Spectrum.Peak[] peaks = _pickPeaks(scan);

        return ExtractPeptideFeatures(_run, peaks);//, scan);
    }


    protected Spectrum.Peak[] _pickPeaks(Scan scan)
    {
        float[][] rawSpectrum = scan.getSpectrum();
        float[][] spectrum = Spectrum.ResampleSpectrum(rawSpectrum, _mzRange,
                SpectrumResampler.getResampleFrequency(), false);
        spectrum[1] = _smooth(spectrum[1]);

        // not sure if this is technically noise, we've already smoothed
        // and perhaps removed background,
        // but it's what we consider a small peak
        _noise = Spectrum.Noise(spectrum[1], 0, spectrum[0].length);
        _noise = Math.max(_noise, 0.1);

        //
        // Generate Peaks
        //
        Spectrum.Peak[] peaks = Spectrum.PickPeaks(spectrum, _noise);
        for (Spectrum.Peak peak : peaks)
            peak.scan = scan.getNum();
        return peaks;
    }

    /**
     * For each single-peak feature:
     * Consider all features in a small window around it
     *       If multiple have the same m/z, take the one with the lower scan, discard the others
     *       Toss out anything < 1/10 the intensity of the highest-intense peak
     *       Consider all possible charge states:
     *               For each peak, determine the mass distance from the highest-intense peak
     *               If within tolerance, create a new multi-peak feature, make sure it contains
     *               the highest-intensity peak, score it (a complicated process),
     *               add it to a list of candidates
     *       Pick the best charge state and create a feature
     *
     * @param run
     * @param peaksIN
     * @return
     */
    protected ArrayList<Feature> ExtractPeptideFeatures(MSRun run, Spectrum.Peak[] peaksIN)
    {
        _logDebug("ExtractPeptideFeatures(" + peaksIN.length + ")");

        //
        // store peaks in 2-D tree for searching
        //
        Tree2D peaks2D = new Tree2D();
        for (Spectrum.Peak p : peaksIN)
        {
            peaks2D.add(p.scan, p.mz, p);
        }

        // OUTPUT list
        ArrayList<Feature> featureList = new ArrayList<Feature>();

        Spectrum.Peak[] peaksByIntensity = peaksIN.clone();
        Arrays.sort(peaksByIntensity, Spectrum.comparePeakIntensityDesc);

        //
        // UNDONE: PROBLEM WITH VERY LARGE 1+ FEATURES
        // in addition to all the other problems, we've observed a strange
        // phenomenon with very large peaks in 1+ charge features.  There may be an "echo"
        // of the feature at about .5-.7 daltons after the first, and sometimes second peak.
        // Thie echo can be very large (as large as the second "real" peak).  Moreover, the
        // echo seems to shift or compress the previous peak downward in MZ.  If this effect
        // is large enough it can break the corret identification of the main feature.
        //
        // There's probably a name for this phenomenon, but I don't know what it is.
        //
        // TODO: ignore the "echo" feature
        // TODO: correct for the shift in the main peaks
        //

        _logDebug("ExtractPeptideFeatures: using noise = " + _noise);
        for (Spectrum.Peak peakHighest : peaksByIntensity)
        {
            if (peakHighest.excluded)
                continue;

            boolean debugPeak = false;
            if (debugMZ != 0)
            {
                // note the 100-140 check will put the feature in the middle of the detail pane
                if(Math.abs(peakHighest.mz - debugMZ) < 1 && peakHighest.scan > 100 &&
                        peakHighest.scan < 140)
                    debugPeak = true;
            }

            float mzWindowStart = peakHighest.mz - 2.1F;
            float mzWindowEnd = peakHighest.mz + 6.1F;
            float scanStart = peakHighest.scan - 9;
            float scanEnd = peakHighest.scan + 9;

            // get collection of nearby peaks
            Spectrum.Peak[] peaks;
            {
                ArrayList<Spectrum.Peak> peakList =
                        peaks2D.getPoints(scanStart, mzWindowStart, scanEnd, mzWindowEnd);
                for (int p = peakList.size()-1; p>=0 ; p--)
                {
                    if ((peakList.get(p)).excluded)
                        peakList.remove(p);
                }
                // eliminate possible duplicates and get one list of peaks sorted by mz
                Spectrum.Peak[] peaksNear = peakList.toArray(new Spectrum.Peak[peakList.size()]);
                assert peaksNear.length > 0;
                Arrays.sort(peaksNear, Spectrum.comparePeakIntensityDesc);
                Arrays.sort(peaksNear, Spectrum.comparePeakMzAsc);
                int lenPeaks = 0;
                for (int s = 1; s < peaksNear.length; s++)
                {
                    if (peaksNear[lenPeaks].mz == peaksNear[s].mz)
                    {
                        if (Math.abs(peaksNear[s].scan-peakHighest.scan) <
                                Math.abs(peaksNear[lenPeaks].scan-peakHighest.scan))
                            peaksNear[lenPeaks] = peaksNear[s];
                    }
                    else
                    {
                        peaksNear[++lenPeaks] = peaksNear[s];
                    }
                }
                ++lenPeaks;
                peaks = new Spectrum.Peak[lenPeaks];
                System.arraycopy(peaksNear, 0, peaks, 0, lenPeaks);
            }

            if (debugPeak)
            {
                _logDebug("\n\npeaks");
                for (Spectrum.Peak peak : peaks)
                {
                    _logDebug((peak == peakHighest ? " *" : "  ") + peak.toString() +
                            "\t" + (peak.mz - peakHighest.mz));
                }
                _logDebug("scored features");
            }

            //
            // generate candidate features
            //

            ArrayList<Feature> candidates = new ArrayList<Feature>();

            //dhmay -- I don't really understand this parameter, why it has this exact value
            double delta = 1.0 / (SpectrumResampler.getResampleFrequency() - 1);

            for (Spectrum.Peak p : peaks)
            {
                if (p.mz > peakHighest.mz)
                    break;
                if (p.excluded)
                    continue;
                // even at high end of mass range (6400amu) leading peak is not expected
                // to be smaller than about 1/6 highest peak, and typically much larger
                if (p.intensity < peakHighest.intensity / 10.0F)
                    continue;
                double distance = peakHighest.mz - p.mz;
                for (int charge = _maxCharge; charge >= 1; charge--)
                {
                    double r = distanceNearestFraction(distance, charge);
                    if (r < 2 * delta)
                    {
                        Feature f = newFeatureRange(peakHighest, p.mz, charge);
                        ScoreFeature(f, peaks);

                        if (f.peaks == 1)
                            continue;

                        boolean containsPeakHighest = false;
                        for (Spectrum.Peak comprisedPeak : f.comprised)
                        {
                            if (comprisedPeak == peakHighest)
                            {
                                containsPeakHighest = true;
                                break;
                            }
                        }
                        if (!containsPeakHighest)
                            continue;

                        candidates.add(f);
                    }
                }
            }

            Feature bestFeature = BestFeature(candidates);

            if (bestFeature == null)
            {
                // 0+
                bestFeature = newFeatureRange(peakHighest, peakHighest.mz, 0);
                bestFeature.comprised = new Spectrum.Peak[] {peakHighest};
            }

            // CONSIDER: use expected largest peak instead of largest actual peak
            if (peakHighest instanceof Feature)
                bestFeature.totalIntensity = ((Feature)peakHighest).totalIntensity;

            //
            // CONSIDER: back check against found features
            //   could this be part of one of those??
            //

            //
            // UNDONE: quality tests, for now we output everything for evaluation
            //
            QualityTest:
            {
                if (bestFeature.getPeaks() == 1)
                {
                    bestFeature.setCharge(0);
                    bestFeature.setKl(-1);
                }

                // Now that charge is set, compute mass
                bestFeature.updateMass();

                if (debugPeak)
                {
                    _logDebug("bestFeature");
                    _logDebug("  " + bestFeature.toString());
                }

                // Got a good one

                //
                // fix ups
                //
                if (peakHighest instanceof Feature && 0 != ((Feature)peakHighest).getTime())
                    bestFeature.setTime(((Feature)peakHighest).getTime());
                else
                {
                    try
                    {
                        // sometimes we're not handed real scan numbers so catch()
                        int index = run.getIndexForScanNum(peakHighest.scan);
                        if (index < 0)
                            index = -(index+1);
                        bestFeature.setTime((float)run.getScan(index).getDoubleRetentionTime());
                    }
                    catch (Throwable t)
                    {
                    }
                }

                featureList.add(bestFeature);
            }

            //
            // exclude all peaks explained by bestFeature (even if we don't output it)
            //
            if (debugPeak)
                _logDebug("exclude peaks");

            for (Spectrum.Peak p : bestFeature.comprised)
            {
                if (null != p)
                {
                    p.excluded = true;
                    if (debugPeak)
                        _logDebug("  " + p.toString());
                }
            }
            assert peakHighest.excluded;
        }

        return featureList;
    }


    private Feature newFeatureRange(Spectrum.Peak peakHighest, float mz, int charge)
    {
        Feature f;
        if (peakHighest instanceof Feature)
            f = new Feature((Feature)peakHighest);
        else
            f = new Feature(peakHighest);
        f.charge = charge;
        f.mz = mz;
        return f;
    }



    /**
     * I used to do this pairwise, but had problem with weird
     * non-trasitive comparison heurisitics.
     *
     * This could undoubtably be improved, maybe with a combined
     * scoring function fn(charge, peaks, kl)?  or new binary
     * feature comparision BetterFeature(a,b)
     *
     * @param candidates
     * @return
     */
    static Feature BestFeature(ArrayList<Feature> candidates)
    {
        if (candidates.size() == 0)
            return null;

        //
        // Based on very little obvservation, it seems that the distribution
        // of KL scores tends to be bimodal.  There are the terrible scores
        // for guess that are hopelessly wrong, and decent scores for features
        // with the right charge (or factor of the right charge)
        //

        //
        // Toss the bad scores. .5 is arbitrary (and maybe too big).
        // It's hard to choose dynamically as the input set is often small
        //

        double min = candidates.get(0).kl;
        double max = min;
        for (Feature candidate : candidates)
        {
            min = Math.min(candidate.kl, min);
            max= Math.max(candidate.kl, max);
        }
        for (int i=candidates.size()-1 ; i>=0 ; i--)
        {
            Feature candidate = candidates.get(i);
            if (candidate.kl > min + 0.5)
                candidates.remove(i);
        }
        if (candidates.size() == 1)
            return candidates.get(0);

        //
        // other heuristics given a set of "not awful" features
        //

        Set<Feature> prunedSet = new HashSet<Feature>(candidates);
        // NOTE we rely on the candidate list being sorted by MZ ASC,CHARGE DESC,
        // see ExtractPeptideFeatures
        for (int i=0 ; i<candidates.size()-1 ; i++)
        {
            Feature a = candidates.get(i);
            for (int j=i+1 ; j<candidates.size() ; j++)
            {
                Feature b = candidates.get(j);
                assert a.mzPeak0 < b.mzPeak0  || a.charge > b.charge;   // see note above

                // if a "contains" b, remove b
                // CONSIDER: kl difference limiter (relying on filter above for now)
                if (a.charge % b.charge == 0 && a.ContainsPeak(b.comprised[0]))
                {
                    // if the test above is true, this test should be fairly redundant
                    if (a.peaks > b.peaks || a.peaks >= 6)
                        prunedSet.remove(b);
                }
            }
        }
        if (prunedSet.size() == 1)
            return prunedSet.iterator().next();

        //
        // we've still got more than one candidate!
        // probably can't explain all relatively big peaks with just one feature...
        //

        // Here's where we need a good secondary score function
        // How do we weight peaks, kl, sum(intensity)
        // in practice it's usually an easy choice (the peaks penalty only rarely changes result)

        Feature[] finalists = prunedSet.toArray(_emptyFeatureArray);
        Arrays.sort(finalists, new Comparator<Feature>()
        {
            double _score(Feature f)
            {
                return f.peaks * 0.1 - f.kl;
            }
            public int compare(Feature o1, Feature o2)
            {
                double score1 = _score(o1);
                double score2 = _score(o2);
                return score1 < score2 ? 1 : score1 > score2 ? -1 : 0; // sort DESC by score
            }
        });

        for (int i=0 ; i<finalists.length-1 ; i++)
            finalists[i].next = finalists[i+1];
        return finalists[0];
    }


    static int _closest(Spectrum.Peak[] peaks, float x, int start)
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


    /**
     * score a candidate feature f
     * <p/>
     * Expecting charge=0 features and peaks
     *
     * TODO: document this better
     *
     * @param f
     * @param peaks
     */
    public void ScoreFeature(Feature f, Spectrum.Peak[] peaks)
    {
        int charge = f.charge;
        float invCharge = 1.0F / charge;
        float mass = (f.mz-Spectrum.HYDROGEN_ION_MASS) * charge;
        float[] expected = Spectrum.Poisson(mass);

        Spectrum.Peak[] peptidePeaks = new Spectrum.Peak[10];

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

        //dhmay -- I don't really understand this parameter, why it has this exact value
        double delta = 1.0 / (SpectrumResampler.getResampleFrequency() - 1);

        for (int Pi = 0; Pi < peptidePeaks.length; Pi++)
        {
            float mzPi = mzP0 + Pi * invCharge + (distCount > 0 ? distSum/distCount : 0);
            p = _closest(peaks, mzPi, p);
            Spectrum.Peak peakFound = peaks[p];
            float dist = Math.abs(peakFound.mz - mzPi);

            boolean found = !peakFound.excluded && Math.abs(dist) < 2.5 * delta;


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
            mean += (peptidePeaks[i].mz - i*invCharge) * peptidePeaks[i].intensity;
            weight += peptidePeaks[i].intensity;
        }
        f.mz = mean/weight;

        //
        // Calculate quality measures
        //

        float[] signal = new float[6];
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

        // some sort of sumSquares distance...
        // NOT TESTED YET, EXPERIMENTAL
        int lastPeak = 0;
        for (int i=0 ; i<signal.length ; i++)
            if (null != peptidePeaks[i])
                lastPeak = i;
        float sumDist = 0;
        for (int i=0 ; i<=lastPeak ; i++)
        {
            float dist = 1; // what should I use here....
            if (null != peptidePeaks[i])
            {
                float distMZ = Math.abs((f.mz + i*invCharge) - peptidePeaks[i].mz);
                distMZ = Math.max(0, distMZ-1/(2* SpectrumResampler.getResampleFrequency()));
                float distIn = Math.abs(signal[i] - expected[i]);
                // need to pick units to scale by, 0 means good 1 is pretty bad...
                distMZ *= 6; // 1/6 of a Da is bad
                distIn *= 4; // 1/4 in scaled intesity is bad
                dist = distMZ * distMZ + distIn * distIn;
                sumDist += dist * expected[i];
            }
        }
        f.dist = sumDist;
    }


    static double distanceNearestFraction(double x, double denom)
    {
        // scale so we can do distance to nearest integer
        double t = x * denom;
        double delta = Math.abs(t - Math.round(t));
        return delta / denom;
    }

    protected static void _logDebug(String s)
    {
        _log.debug(s);
        //System.err.println(s);
    }
}
