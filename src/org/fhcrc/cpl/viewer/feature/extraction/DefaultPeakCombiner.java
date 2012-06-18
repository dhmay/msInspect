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
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Default class for creating Features from Peaks
 *
 * 07/29/20008 dhmay enhanced to handle negatively-charged ions
 */
public class DefaultPeakCombiner extends BasePeakCombiner
{
    private static Logger _log = Logger.getLogger(DefaultPeakCombiner.class);

    public static final double DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS = 1.0;

    //default is positively charged ions
    protected boolean negativeChargeMode = false;

    protected double _maxAbsDistanceBetweenPeaks =
            DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS;


    protected FeatureScorer _featureScorer;

    protected boolean keepStatistics = false;

    //track the number of candidates evaluated for each feature
    public List<Float> numCandidatesList = new ArrayList<Float>();
//    public List<Float> klScoresList = new ArrayList<Float>();

    public DefaultPeakCombiner()
    {
        _featureScorer = new DefaultFeatureScorer();
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
    public Feature[] createFeaturesFromPeaks(MSRun run, Spectrum.Peak[] peaksIN)
    {
        _log.debug("ExtractPeptideFeatures(" + peaksIN.length + ")");

        //TODO: really should set some parameters in FeatureScorer here (maxDist and resamplingFreq), but not sure if we want to surface those in interface

        //
        // store peaks in 2-D tree for searching
        //
        Tree2D peaks2D = new Tree2D();
        for (Spectrum.Peak p : peaksIN)
        {
            peaks2D.add(p.scan, p.mz, p);
        }

        // OUTPUT list
        List<Feature> featureList = new ArrayList<Feature>();

        Spectrum.Peak[] peaksByIntensity = peaksIN.clone();
        Arrays.sort(peaksByIntensity, Spectrum.comparePeakIntensityDesc);

        //
        // UNDONE: PROBLEM WITH VERY LARGE 1+ FEATURES
        // in addition to all the other problems, we've observed a strange
        // phenomenon with very large peaks in 1+ charge features.  There may be an "echo"
        // of the feature at about .5-.7 daltons after the first, and sometimes second peak.
        // Thie echo can be very large (as large as the second "real" peak).  Moreover, the
        // echo seems to shift or compress the previous peak downward in MZ.  If this effect
        // is large enough it can break the correct identification of the main feature.
        //
        // There's probably a name for this phenomenon, but I don't know what it is.
        //
        // TODO: ignore the "echo" feature
        // TODO: correct for the shift in the main peaks
        //


        for (Spectrum.Peak peakHighest : peaksByIntensity)
        {
            if (peakHighest.excluded)
                continue;

            boolean debugPeak = false;
//            if (debugMZ != 0)
//            {
//                // note the 100-140 check will put the feature in the middle of the detail pane
//                if(Math.abs(peakHighest.mz - debugMZ) < 1 && peakHighest.scan > 100 &&
//                        peakHighest.scan < 140)
//                    debugPeak = true;
//            }



            float mzWindowStart = peakHighest.mz - 2.1F;
            float mzWindowEnd = peakHighest.mz + 6.1F;
            float scanStart = peakHighest.scan - 9;
            float scanEnd = peakHighest.scan + 9;

            // get collection of nearby peaks
            Spectrum.Peak[] peaks;

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


            if (debugPeak)
            {
                _log.debug("\n\npeaks");
                for (Spectrum.Peak peak : peaks)
                {
                    _log.debug((peak == peakHighest ? " *" : "  ") + peak.toString() +
                            "\t" + (peak.mz - peakHighest.mz));
                }
                _log.debug("scored features");
            }

            //
            // generate candidate features
            //

            ArrayList<Feature> candidates = new ArrayList<Feature>();

            double maxResampledDistanceBetweenFeatures =
                    _maxAbsDistanceBetweenPeaks / (SpectrumResampler.getResampleFrequency() - 1);
            for (Spectrum.Peak p : peaks)
            {
            	// we are looking for the mono-isotopic peak 
                if (p.mz > peakHighest.mz)
                    break;
                if (p.excluded)
                    continue;
                // even at high end of mass range (6400amu) leading peak is not expected
                // to be smaller than about 1/6 highest peak, and typically much larger
                if (p.intensity < peakHighest.intensity / 10.0F)
                    continue;
                double distance = peakHighest.mz - p.mz;

                for (int absCharge = Math.abs(_maxCharge); absCharge >= 1; absCharge--)
                {
                    double r = distanceNearestFraction(distance, absCharge);
                    if (r < 2 * maxResampledDistanceBetweenFeatures)
                    {
                        Feature f = newFeatureRange(peakHighest, p.mz,
                                negativeChargeMode ? -absCharge : absCharge);
                        f.dist = _featureScorer.scoreFeature(f, peaks);
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
//too much detail
//                       _log.debug("\tCandidate: " + f);
                    }
                }

            }



            Feature bestFeature = determineBestFeature(candidates);



            if (bestFeature == null)
            {
                // 0+
                bestFeature = newFeatureRange(peakHighest, peakHighest.mz, 0);
                bestFeature.comprised = new Spectrum.Peak[] {peakHighest};
//too much detail
//              _log.debug("No passing candidate, inventing: " + bestFeature);
            }
            else
            {
//System.err.println("DINGDINGDING!!!!!");

//too much detail
//              _log.debug("Winner: " + bestFeature);
            }
            // CONSIDER: use expected largest peak instead of largest actual peak
            if (peakHighest instanceof Feature)
                bestFeature.totalIntensity = ((Feature)peakHighest).totalIntensity;

            if (keepStatistics)
            {
                numCandidatesList.add((float) candidates.size());
//                klScoresList.add(bestFeature.kl);
            }

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
                    _log.debug("bestFeature");
                    _log.debug("  " + bestFeature.toString());
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
                _log.debug("exclude peaks");

            for (Spectrum.Peak p : bestFeature.comprised)
            {
                if (null != p)
                {
                    p.excluded = true;
                    if (debugPeak)
                        _log.debug("  " + p.toString());
                }
            }
            assert peakHighest.excluded;
        }
        return featureList.toArray(new Feature[featureList.size()]);
    }


    static double distanceNearestFraction(double x, double denom)
    {
        // scale so we can do distance to nearest integer
        double t = x * denom;
        double delta = Math.abs(t - Math.round(t));
        return delta / denom;
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
    static Feature determineBestFeature(ArrayList<Feature> candidates)
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
                assert a.mzPeak0 < b.mzPeak0  || Math.abs(a.charge) > Math.abs(b.charge);   // see note above

                // if a "contains" b, remove b
                // CONSIDER: kl difference limiter (relying on filter above for now)
                if (Math.abs(a.charge) % Math.abs(b.charge) == 0 && a.ContainsPeak(b.comprised[0]))
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

        Feature[] finalists = prunedSet.toArray(new Feature[prunedSet.size()]);
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


    /**
     * Create a new feature representing this peak at a given charge
     * @param peakHighest
     * @param mz
     * @param charge
     * @return
     */
    protected Feature newFeatureRange(Spectrum.Peak peakHighest, float mz, int charge)
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

    public double getMaxAbsDistanceBetweenPeaks()
    {
        return _maxAbsDistanceBetweenPeaks;
    }

    public FeatureScorer getFeatureScorer()
    {
        return _featureScorer;
    }

    public void setFeatureScorer(FeatureScorer featureScorer)
    {
        _featureScorer = featureScorer;
    }


    public boolean isKeepStatistics()
    {
        return keepStatistics;
    }

    public void setKeepStatistics(boolean keepStatistics)
    {
        this.keepStatistics = keepStatistics;
        if (_featureScorer != null)
            ((DefaultFeatureScorer) _featureScorer).setKeepStatistics(keepStatistics);
    }

    public void plotStatistics()
    {
        if (keepStatistics)
        {
            PanelWithHistogram pwhNumCandidates =
                    new PanelWithHistogram(numCandidatesList, "#Candidate Features");
            pwhNumCandidates.displayInTab();

//            PanelWithScatterPlot pwsp =
//                    new PanelWithScatterPlot(numCandidatesList, klScoresList, "KL vs. #Candidates");
//            pwsp.setPointSize(2);
//            pwsp.displayInTab();

            ((DefaultFeatureScorer) _featureScorer).plotStatistics();
        }

    }

    public boolean isNegativeChargeMode()
    {
        return negativeChargeMode;
    }

    public void setNegativeChargeMode(boolean negativeChargeMode)
    {
        this.negativeChargeMode = negativeChargeMode;
    }
}
