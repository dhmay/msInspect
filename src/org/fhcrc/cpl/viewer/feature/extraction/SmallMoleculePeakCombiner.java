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
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.viewer.quant.PeakOverlapCorrection;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * class for creating Features from Peaks, specific to small-molecule analysis, e.g., metabolites.
 *
 * Differs from FeatureStrategyPeakClusters in that:
 * -hardcoded maximum charge of 2, maximum m/z range 3 Thompsons
 * -most-intense peak is assumed to be monoisotope
 * -"k/l" calculation is actually a simple comparison of the log-ratio of the first two peaks with a theoretical
 * mass-logratio relationship derived from analyzing the HMDB databases 
 * -Determination of best feature is purely a decision between charge states.  Based on number of
 * peaks (more is better), then KL
 * -hard filter on candidate features if any peak intensity is more than MAX_ANY_FIRST_PEAK_PROPORTION of first
 * -Single peak features allowed, assume charge 1, given KL -1
 * -No missing peaks allowed.  comprised array has no null values
 *
 */
public class SmallMoleculePeakCombiner extends BasePeakCombiner
{
    private static Logger _log = Logger.getLogger(SmallMoleculePeakCombiner.class);

    public static final double DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS = 1.0;

    //default is positively charged ions
    protected boolean negativeChargeMode = false;

    protected double _maxAbsDistanceBetweenPeaks =
            DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS;

    //This should be higher if considering chlorine-containing ions
    public static final float MAX_ANY_FIRST_PEAK_PROPORTION = 0.5f;

    //maximum distance between peaks to be considered tied together
    protected double _maxResampledDistanceBetweenFeatures =
            DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS / (SpectrumResampler.getResampleFrequency() - 1);

    protected boolean keepStatistics = false;

    //5 peaks possible with chlorine boosting the 3rd & beyond, but there are so few chlorine-containing peaks
    //it's too damaging to consider them
    protected int maxPeaks = 4;

    //maximum distance in scans from the scan range of the base peak to the apex scan of any other peak
    protected int maxScanRangeDist = 1;

    //maximum distance from apex scan of first peak to apex scan of any other peak
    protected int maxApexScanDist = 6;

    //coefficients for the equation relating metabolite mass to the log-ratio of peak 1 intensity / peak 2 intensity
    //Empirically derived
    public static final double[] MASS_2PEAK_LOGRATIO_COEFFS = new double[]
            {
                    4.260321140289307,
                    -0.015644630417227745,
                    3.08442504319828E-5,
                    -2.9751459962312765E-8,
                    1.0629329881550742E-11
            };
    //Standard deviation of error around regression above, from HMDB database.  Empirically derived
    public static final double MASS_2PEAK_LOGRATIO_ERROR_SD = 0.40649617;

//    protected List<Integer> allCandidatePeakCounts = new ArrayList<Integer>();

    //track the number of candidates evaluated for each feature
    public List<Float> numCandidatesList = new ArrayList<Float>();      
//    public List<Float> klScoresList = new ArrayList<Float>();

    public SmallMoleculePeakCombiner()
    {
        _maxCharge = 2;
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
        _maxCharge = 2;

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

//List<Integer> scanDeviations = new ArrayList<Integer>();
        for (Spectrum.Peak peakHighest : peaksByIntensity)
        {
            if (peakHighest.excluded)
                continue;

            //dhmay I don't know why we get 0-intensity peaks, but we do. Below noise?  At any rate, for small
            //molecules, we don't want 'em
            if (peakHighest.intensity == 0)
            {
                peakHighest.excluded = true;
                continue;
            }

            boolean debugPeak = false;

            float mzWindowStart = peakHighest.mz - 1.01F;
            float mzWindowEnd = peakHighest.mz + (maxPeaks-1) + .1F;
            float scanStart = ((Feature) peakHighest).scanFirst - maxApexScanDist;
            float scanEnd = peakHighest.scan + maxApexScanDist;

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
//lensPeaks.add(lenPeaks);
            peaks = new Spectrum.Peak[lenPeaks];
            System.arraycopy(peaksNear, 0, peaks, 0, lenPeaks);
            if (debugPeak)
            {
                _log.debug("\n\npeaks=" + peaks.length);
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


            for (int absCharge = Math.abs(_maxCharge); absCharge >= 1; absCharge--)
            {
                Feature f = newFeatureRange(peakHighest, peakHighest.mz,
                        negativeChargeMode ? -absCharge : absCharge);
                fleshOutFeature(f, peaks);

                //if we've only got one peak, forget it -- just use the charge-1 version, no matter how many peaks it has
                if (absCharge > 1 && f.peaks == 1)
                    continue;
//if (absCharge > 1)
//            System.err.println("charge" + f.charge + ", peaks=" + f.peaks + ", kl=" + f.kl);

//for (int i=0; i<f.comprised.length; i++) if (f.comprised[i] == null) System.err.println("NULL!!! " + i + " / " + (f.comprised.length-1));



//dhmay commenting out the requirement that we have 2+ peaks.  BIG CHANGE!!!
//                if (f.peaks == 1)
//                    continue;

                if (f.comprised[0] == null)// || f.comprised[1] == null)
                    continue;

                //add no candidates with any peak with intensity more than MAX_ANY_FIRST_PEAK_PROPORTION
                //times the first peak intensity.
                //Also, add no candidates with any missing peaks
                boolean shouldAdd = true;

                for (int i=1; i<f.comprised.length; i++)
                {
//                    if (hasMissing && f.comprised[i] != null)
//                        hasPresentAfterMissing = true;

                    if (Math.round((f.comprised[i].getMz() - f.getMz()) * (float) absCharge) != i)
                    {
//if (absCharge > 1) System.err.println("fail, ch=" + absCharge + ", " + f.comprised[i].getMz() + ", " + f.getMz() + ", " + i + ", " + (f.comprised[i].getMz() - f.getMz()));
                        shouldAdd = false;
                        break;
                    }

                    if (f.comprised[i].intensity > MAX_ANY_FIRST_PEAK_PROPORTION * f.comprised[0].intensity)
                    {
                        shouldAdd = false;
                        break;
                    }
                }
                if (shouldAdd)
                {
//if (absCharge > 1) System.err.println("\tAdding");
                    candidates.add(f);
//allCandidatePeakCounts.add(f.peaks);
                }

            }


            Feature bestFeature = determineBestFeature(candidates);
//if (bestFeature != null && bestFeature.comprised.length > 1) System.err.println("MULTI!");

            if (bestFeature == null)
            {
                // 0+
                bestFeature = newFeatureRange(peakHighest, peakHighest.mz, 1);
                bestFeature.comprised = new Spectrum.Peak[] {peakHighest};
            }
            // CONSIDER: use expected largest peak instead of largest actual peak
            if (peakHighest instanceof Feature)
                bestFeature.totalIntensity = ((Feature)peakHighest).totalIntensity;

//if (bestFeature.peaks > 1)
//{
//    int devScans = 0;
//    int peak2Scan = bestFeature.comprised[1].scan;
//    if (peak2Scan > bestFeature.getScanLast())
//        devScans = peak2Scan - bestFeature.getScanLast();
//    if (peak2Scan < bestFeature.getScanFirst())
//        devScans = bestFeature.getScanFirst() - peak2Scan;
//    scanDeviations.add(devScans);
//}

            if (keepStatistics)
            {
                numCandidatesList.add((float) candidates.size());
            }

//boolean wasCharge2 = bestFeature.charge == 2;
//if (bestFeature.charge == 2)
//{
//System.err.println("charge2, kl=" + bestFeature.kl);
//}

            QualityTest:
            {
                if (bestFeature.getPeaks() == 1 )
                {
                    bestFeature.setCharge(negativeChargeMode ? -1 : 1);
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
//if ("1chlorine".equals(bestFeature.getDescription()))
// System.err.println("chlor mass=" + bestFeature.mass);

                annotateFeatureDescription(bestFeature);

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
//new PanelWithHistogram(scanDeviations, "scan deviations").displayInTab();

//new PanelWithHistogram(allCandidatePeakCounts, "peaks").displayInTab();
//new PanelWithHistogram(lensPeaks, "peaks").displayInTab();

//int numMulti = 0;
//        List<Float> atts = new ArrayList<Float>();
//for (Feature feature : featureList)
//{
//    if (feature.comprised.length > 1) {numMulti++; atts.add(feature.kl); System.err.println(feature.kl);}
//}
//System.err.println("Multis: " + numMulti);
//new PanelWithHistogram(atts, "atts").displayInTab();

        
        return featureList.toArray(new Feature[featureList.size()]);
    }

    protected void annotateFeatureDescription(Feature feature)
    {
        Spectrum.Peak[] peaks = feature.comprised;
//        if (peaks.length > 2) System.err.println(peaks[2].intensity/peaks[1].intensity);
        if (peaks.length > 2 && peaks[2].intensity > peaks[1].intensity * 1.2)
        {
//            if (peaks[2] > peaks[1]) System.err.println("Better! mass=" + mass + ", peaks: " + peaks[0] + ", " +peaks[1] + ", " + peaks[2]);

            feature.setDescription("1Cl");
        }
    }


    public void fleshOutFeature(Feature f, Spectrum.Peak[] peaks)
    {
        int absCharge = Math.abs(f.charge);
        float invAbsCharge = 1.0F / absCharge;
        float mass = 0f;
        if (f.charge >= 0)
            mass = (f.mz-Spectrum.HYDROGEN_ION_MASS) * absCharge;
        else mass = f.mz * absCharge;
        Spectrum.Peak[] peptidePeaks = new Spectrum.Peak[maxPeaks];

        // find start position in peak list
        int p = Arrays.binarySearch(peaks, f, Spectrum.comparePeakMzAsc);
        if (p < 0) p = -(p + 1);
        p = Math.max(0, p - 1);
        int pLastFound = p;

        int lastPeptidePeakIndexFound = 0;
        float mzP0 = f.mz;
        float sum = 0.0F;
        float intensityLast = 0.0F;
        float intensityHighest = 0.0F;
        boolean skippedPeak = false;
        float distSum = 0;
        float distCount = 0;
        float lastPeakPairScanOffset = 0;
        for (int Pi = 0; Pi < peptidePeaks.length; Pi++)
        {
            float mzPi = mzP0 + Pi * invAbsCharge + (distCount > 0 ? distSum/distCount : 0);
            p = findClosestPeak(peaks, mzPi, p);
            Spectrum.Peak peakFound = peaks[p];
            float dist = Math.abs(peakFound.mz - mzPi);

//QUESTION: why 2.5?
            boolean found =
                    !peakFound.excluded &&
                            Math.abs(dist) < 2.5 * _maxResampledDistanceBetweenFeatures;

            //dhmay adding check to make sure this new peak is within 1 scan of orig peak boundaries
            //This is specific to low-mass stuff: first peak will dominate, should elute as long
            //or longer than second
            if (Pi > 0)
            {
                int scanDist = 0;
                if (peakFound.scan > f.scanLast)
                    scanDist = peakFound.scan - f.scanLast;
                if (peakFound.scan < f.scanFirst)
                    scanDist = f.scanFirst - peakFound.scan;
                if (Math.abs(scanDist) > maxScanRangeDist)
                    found = false;
            }

            //dhmay adding check for peak offset.  From the third peak onward, do not add a peak if
            //the absolute difference in scans between it and base peak is greater than the previous difference
            //by 4 or more /or/ it's greater by 2 or more and double the previous difference
            if (Pi > 1)
            {
                float currentPeakPairScanOffset =
                        Math.abs(peakFound.scan - f.scan);
                if (currentPeakPairScanOffset >= lastPeakPairScanOffset + 4 ||
                     ((currentPeakPairScanOffset >= lastPeakPairScanOffset * 2)) &&
                    (currentPeakPairScanOffset >= lastPeakPairScanOffset + 2))
                    found = false;
            }


            if (!found)
            {
                //dhmay changing this for SmallMolecule.  NO missing peaks allowed
                // miss one, call it quits
                if (Pi > 0)
                    break;
            }
            else
            {
                //if this peak intensity is really really small, don't add it, and stop right here
                if (Pi > 0 && (peakFound.intensity < intensityHighest/200f ||
                    peakFound.intensity < 2 * f.median))
                        break;

                // if peak is too BIG stop:
                // CONSIDER do this at end, use fn of signal v. expected?
                // don't like hard cut-offs
                //dhmay: this is fine for small molecules, only checks if peak is second, or fourth & higher
                if ((Pi == 1 || Pi >= 3)  &&
                        peakFound.intensity > intensityLast * 1.33 + f.median)
                    break;

                intensityHighest = Math.max(peakFound.intensity, intensityHighest);


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
                lastPeptidePeakIndexFound = Pi;
                intensityLast = peakFound.intensity;

                if (Pi > 0)
                {
                    lastPeakPairScanOffset = Math.abs(peptidePeaks[Pi].scan - f.scan);
                }
            }
        }

        if (lastPeptidePeakIndexFound < maxPeaks)
        {
            Spectrum.Peak[] peaksNew = new Spectrum.Peak[lastPeptidePeakIndexFound+1];
            System.arraycopy(peptidePeaks, 0, peaksNew, 0,lastPeptidePeakIndexFound+1 );
            peptidePeaks = peaksNew;
//System.err.println(peaksNew.length + " of " + maxPeaks);
        }

//if (peptidePeaks.length == 4) System.err.println(peptidePeaks[0].mz + ", " + peptidePeaks[1].mz + ", " +peptidePeaks[2].mz + ", " +peptidePeaks[3].mz );        

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
        int signalLength = 2;
        float[] signal = new float[Math.min(signalLength, peptidePeaks.length)];
        for (int i = 0; i < signal.length; i++)
        {
            signal[i] = peptidePeaks[i].intensity;
            sum += signal[i];
        }
        for (int i = 0; i < signal.length; i++)
            signal[i] /= sum;
//if (Float.isNaN(signal[0])) System.err.println("***" + peptidePeaks[0].intensity + ", " + sum);
        //dhmay implementing check for chlorine
        f.kl = calcKlSmallMolecule(mass, signal);
        //TODO: Do we need another score?
        //f.score = f.kl >= 0 ? f.intensity * (float)-Math.log(f.kl + Math.exp(-1)) : 0F;
        f.setPeaks(peptidePeaks.length);
        //skippedPeaks seems to mean that a peak in the array was skipped, not that there's a peak missing
        f.skippedPeaks = skippedPeak;
        f.comprised = peptidePeaks;
        f.mzPeak0 = peptidePeaks[0].mz;
    }

    /**
     * Calculate a "KL" score based on averages from the HMDB database
     *
     * Side effect: set description of suspected chlorine-containing features
     * @param mass
     * @param peaks
     * @return
     */
    protected float calcKlSmallMolecule(float mass, float[] peaks)
    {
        if (peaks.length == 1)
            return -1;
        double idealLogPeakRatio = RegressionUtilities.mapValueUsingCoefficients(MASS_2PEAK_LOGRATIO_COEFFS, mass);
        float actualLogPeakRatio = (float) Math.log(peaks[0] / peaks[1]);
        float kl = (float) Math.abs(actualLogPeakRatio - idealLogPeakRatio) / (float) MASS_2PEAK_LOGRATIO_ERROR_SD;

        return kl;
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
    static Feature determineBestFeature(List<Feature> candidates)
    {
        if (candidates.size() == 0)
            return null;

        if (candidates.size() == 1)
            return candidates.get(0);

        //If there are any features with multiple peaks, don't consider single-peak features
        List<Feature> candidatesWithMultiplePeaks = new ArrayList<Feature>();
        for (Feature feature : candidates)
        {
            if (feature.peaks > 1)
                candidatesWithMultiplePeaks.add(feature);
        }
        if (!candidatesWithMultiplePeaks.isEmpty())
        {
             candidates = candidatesWithMultiplePeaks;
        }
        
        //sort by peaks descending, up to 3, then by kl descending.  Which basically means sort by KL
        Collections.sort(candidates, new Comparator<Feature>()
        {
            public int compare(Feature o1, Feature o2)
            {
                int klCompare = _compareDesc(Math.min(o1.peaks,3), Math.min(o2.peaks,3));
                if (klCompare == 0)
                    return _compareAsc(o1.kl, o2.kl);
                return klCompare;
            }
        });
        return candidates.get(0);
    }

    static int _compareAsc(float a, float b)
    {
        return a == b ? 0 : a < b ? -1 : 1;
    }


    static int _compareDesc(float a, float b)
    {
        return a == b ? 0 : a < b ? 1 : -1;
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

    /**
     * Side effects: sets resampling freq in SpectrumResampler and sets max resampled distance between features
     * @param resamplingFrequency
     */
    public void setResamplingFrequency(int resamplingFrequency)
    {
        SpectrumResampler.setResampleFrequency(resamplingFrequency);
        _maxResampledDistanceBetweenFeatures =
                _maxAbsDistanceBetweenPeaks / (resamplingFrequency - 1);
    }


    public double getMaxAbsDistanceBetweenPeaks()
    {
        return _maxAbsDistanceBetweenPeaks;
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
        if (keepStatistics)
        {
            PanelWithHistogram pwhNumCandidates =
                    new PanelWithHistogram(numCandidatesList, "#Candidate Features");
            pwhNumCandidates.displayInTab();

//            PanelWithScatterPlot pwsp =
//                    new PanelWithScatterPlot(numCandidatesList, klScoresList, "KL vs. #Candidates");
//            pwsp.setPointSize(2);
//            pwsp.displayInTab();

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
