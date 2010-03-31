package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.viewer.quant.PeakOverlapCorrection;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * class for creating Features from Peaks, specific to small-molecule analysis, e.g., metabolites.
 *
 * Differs from FeatureStrategyPeakClusters in that:
 * -hardcoded maximum charge of 2, maximum m/z range 3 Thompsons
 * -most-intense peak is assumed to be monoisotope
 * -k/l calculation is based on first two peaks only
 * -Determination of best feature is purely a decision between charge states.  Based on k/l, then number of
 * peaks (more is better), which effectively means just k/l
 * -hard filter on candidate features if any peak intensity is more than MAX_ANY_FIRST_PEAK_PROPORTION of first
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

    public static final float MAX_ANY_FIRST_PEAK_PROPORTION = 0.8f;

    //maximum distance between features to be considered tied together
    protected double _maxResampledDistanceBetweenFeatures =
            DEFAULT_MAX_ABS_DISTANCE_BETWEEN_PEAKS / (_resamplingFrequency - 1);

    protected boolean keepStatistics = false;

    protected int maxPeaks = 5;

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

//List<Integer> lensPeaks = new ArrayList<Integer>();
        for (Spectrum.Peak peakHighest : peaksByIntensity)
        {
            if (peakHighest.excluded)
                continue;

            boolean debugPeak = false;

            float mzWindowStart = peakHighest.mz - 1.01F;
            //5 peaks possible with chlorine boosting the 3rd & beyond
            float mzWindowEnd = peakHighest.mz + (maxPeaks-1) + .1F;
            float scanStart = peakHighest.scan - 12;
            float scanEnd = peakHighest.scan + 12;

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

                if (f.peaks == 1)
                    continue;

                //add no candidates with any peak with intensity more than MAX_ANY_FIRST_PEAK_PROPORTION
                //times the first peak intensity
                boolean shouldAdd = true;
                for (int i=1; i<f.comprised.length; i++)
                {
                    if (f.comprised[i] != null &&
                        f.comprised[i].intensity > MAX_ANY_FIRST_PEAK_PROPORTION * f.comprised[0].intensity)
                    {
                        shouldAdd = false;
                        break;
                    }
                }
                if (shouldAdd)
                {
                    candidates.add(f);
//allCandidatePeakCounts.add(f.peaks);
                }
            }


            Feature bestFeature = determineBestFeature(candidates);


            if (bestFeature == null)
            {
                // 0+
                bestFeature = newFeatureRange(peakHighest, peakHighest.mz, 0);
                bestFeature.comprised = new Spectrum.Peak[] {peakHighest};
            }
            // CONSIDER: use expected largest peak instead of largest actual peak
            if (peakHighest instanceof Feature)
                bestFeature.totalIntensity = ((Feature)peakHighest).totalIntensity;

            if (keepStatistics)
            {
                numCandidatesList.add((float) candidates.size());
            }


            QualityTest:
            {
                if (bestFeature.getPeaks() == 1 )
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
if ("1chlorine".equals(bestFeature.getDescription()))
 System.err.println("chlor mass=" + bestFeature.mass);
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
//new PanelWithHistogram(allCandidatePeakCounts, "peaks").displayInTab();
//new PanelWithHistogram(lensPeaks, "peaks").displayInTab();
        return featureList.toArray(new Feature[featureList.size()]);
    }


    public void fleshOutFeature(Feature f, Spectrum.Peak[] peaks)
    {
        int absCharge = Math.abs(f.charge);
        float invAbsCharge = 1.0F / absCharge;
        float mass = 0f;
        if (f.charge >= 0)
            mass = (f.mz-Spectrum.HYDROGEN_ION_MASS) * absCharge;
        else mass = f.mz * absCharge;

        Spectrum.Peak[] peptidePeaks = new Spectrum.Peak[3];

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
                            Math.abs(dist) < 2.5 * _maxResampledDistanceBetweenFeatures;

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
        int countPeaks = 0;
        int signalLength = 3;
        float[] signal = new float[signalLength];
        for (int i = 0; i < peptidePeaks.length; i++)
        {
            Spectrum.Peak peak = peptidePeaks[i];
            if (null != peak && (peak.intensity > intensityHighest/50f &&
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

        //dhmay implementing check for chlorine
        f.kl = calcKlSmallMolecule(f, mass, signal);
        //TODO: Do we need another score?
        //f.score = f.kl >= 0 ? f.intensity * (float)-Math.log(f.kl + Math.exp(-1)) : 0F;
        f.setPeaks(countPeaks);
        f.skippedPeaks = skippedPeak;
        f.comprised = peptidePeaks;
        f.mzPeak0 = peptidePeaks[0].mz;
    }

    /**
     * Calculate a KL score assuming no Chlorines, and then assuming 1.  Return whichever's less
     * @param mass
     * @param peaks
     * @return
     */
    protected float calcKlSmallMolecule(Feature feature, float mass, float[] peaks)
    {
        int numPeaks = Math.min(peaks.length, 3);

        float[] idealPeaks = new float[numPeaks];
        System.arraycopy(Spectrum.Poisson(mass), 0, idealPeaks, 0, numPeaks);
        float sum = 0;
        for (float ideal : idealPeaks)
            sum += ideal;
        for (int i=0; i<numPeaks; i++)
            idealPeaks[i] /= sum;
        float klNoChlorine = PeakOverlapCorrection.calcKLUsingTemplate(idealPeaks, peaks);

        if (numPeaks <= 2 || (peaks[1] > peaks[2]))
            return klNoChlorine;
//        return klNoChlorine;

        System.arraycopy(Spectrum.Poisson(mass), 0, idealPeaks, 0, numPeaks);
        sum = 1;
        for (int i=2; i<numPeaks; i++)
        {
            idealPeaks[i] += (idealPeaks[i-2] * 0.2423);
            sum += idealPeaks[i];
        }
        for (int i=0; i<numPeaks; i++)
            idealPeaks[i] /= sum;
        float kl1Chlorine = PeakOverlapCorrection.calcKLUsingTemplate(idealPeaks, peaks);
        if (kl1Chlorine < klNoChlorine)
        {
//            if (peaks[2] > peaks[1]) System.err.println("Better! mass=" + mass + ", peaks: " + peaks[0] + ", " +peaks[1] + ", " + peaks[2]);

            feature.setDescription("1chlorine");
        }
        return Math.min(klNoChlorine, kl1Chlorine);
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
    static Feature determineBestFeature(ArrayList<Feature> candidates)
    {
        if (candidates.size() == 0)
            return null;

        if (candidates.size() == 1)
            return candidates.get(0);
        
        //sort by peaks descending, then by kl ascending
        Collections.sort(candidates, new Comparator<Feature>()
        {
            public int compare(Feature o1, Feature o2)
            {
                int klCompare = _compareAsc(o1.kl, o2.kl);
                if (klCompare == 0)
                    return _compareDesc(o1.peaks, o2.peaks);
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

    public void setResamplingFrequency(int resamplingFrequency)
    {
        super.setResamplingFrequency(resamplingFrequency);
        _maxResampledDistanceBetweenFeatures =
                _maxAbsDistanceBetweenPeaks / (_resamplingFrequency - 1);
    }


    public double getMaxAbsDistanceBetweenPeaks()
    {
        return _maxAbsDistanceBetweenPeaks;
    }

    public void setMaxAbsDistanceBetweenPeaks(double maxAbsDistanceBetweenPeaks)
    {
        this._maxAbsDistanceBetweenPeaks = maxAbsDistanceBetweenPeaks;
        _maxResampledDistanceBetweenFeatures =
                _maxAbsDistanceBetweenPeaks / (_resamplingFrequency - 1);
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
