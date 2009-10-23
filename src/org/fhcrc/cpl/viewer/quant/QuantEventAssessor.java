/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

import java.util.*;
import java.util.List;


/**
 * This uses a big HACK, adding a dummy search score to peptide features to store the flag reason description
 *
 */
public class QuantEventAssessor
{
    public static Logger _log = Logger.getLogger(QuantEventAssessor.class);

    public static final int FLAG_REASON_OK = 0;
    public static final int FLAG_REASON_COELUTING = 1;
    public static final int FLAG_REASON_DISSIMILAR_KL = 2;
    public static final int FLAG_REASON_DISSIMILAR_MS1_RATIO = 3;
    public static final int FLAG_REASON_BIG_2PEAK_KL = 4;
    public static final int FLAG_REASON_UNEVALUATED = 5;
    public static final int FLAG_REASON_OTHER = 6;



    public static final String[] flagReasonDescriptions = new String[]
            {
                    "OK",
                    "Coeluting peptide",
                    "Dissimilar light/heavy KL",
                    "Singlepeak Ratio Different",
                    "Big 2-peak KL",
                    "Unevaluated",
                    "Other",
            };

    public static final String[] flagReasonCodes = new String[]
            {
                    "OK",
                    "CoelutingPeptide",
                    "DissimilarKL",
                    "MS1MS2RatioDiff",
                    "Big2PeakKL",
                    "Unevaluated",
                    "Other"
            };

    public static String getAssessmentCodeDesc(int assessmentCode)
    {
        return flagReasonDescriptions[assessmentCode];
    }

    public static int parseAssessmentCodeString(String curationStatusString)
    {
//System.err.println("   " + curationStatusString);        
        if (flagReasonCodes[FLAG_REASON_OK].equals(curationStatusString))
            return FLAG_REASON_OK;
        if (flagReasonCodes[FLAG_REASON_COELUTING].equals(curationStatusString))
            return FLAG_REASON_COELUTING;
        if (flagReasonCodes[FLAG_REASON_DISSIMILAR_KL].equals(curationStatusString))
            return FLAG_REASON_DISSIMILAR_KL;
        if (flagReasonCodes[FLAG_REASON_DISSIMILAR_KL].equals(curationStatusString))
            return FLAG_REASON_DISSIMILAR_KL;
        if (flagReasonCodes[FLAG_REASON_BIG_2PEAK_KL].equals(curationStatusString))
            return FLAG_REASON_BIG_2PEAK_KL;
        if (flagReasonCodes[FLAG_REASON_OTHER].equals(curationStatusString))
            return FLAG_REASON_OTHER;        
        return FLAG_REASON_UNEVALUATED;
    }

    public static final String REASON_DUMMY_SEARCH_SCORE_NAME = "dummy_flag_desc";

    //Maximum allowed proportion of highest peak intensity that the peak 1Da below monoisotope can have
    public static final float DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF = 0.4f;
    protected float peakBelowIntensityRatioCutoff = DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF;

    public static final int DEFAULT_NUM_SCANS_AROUND_EVENT = 0;
    protected int numScansAroundEventToConsider = DEFAULT_NUM_SCANS_AROUND_EVENT;

    public static final float DEFAULT_PEAK_PPM_TOLERANCE = 50;
    protected float peakPPMTolerance = DEFAULT_PEAK_PPM_TOLERANCE;

    //maximum KL score for either light or heavy peaks, only first two peaks used in calculation
    public static final float DEFAULT_MAX_2PEAK_KL = 10.0f;
    protected float max2PeakKL = DEFAULT_MAX_2PEAK_KL;

    //maximum KL score heavy peaks, using light peaks as template
//    public static final float DEFAULT_MAX_2PEAK_KL_HEAVY_VS_LIGHT = 5f;
//    protected float max2PeakKLHeavyVsLight = DEFAULT_MAX_2PEAK_KL_HEAVY_VS_LIGHT;

    protected int labelType = QuantitationUtilities.LABEL_LYCINE;

    protected boolean showCharts = false;


    public static final float SILAC_LABEL_MASS = 134.115092f;
    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;


    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;


    protected int numPeaksToUse = 5;
    protected int numRawPeaksToKeep = 2;

    //minimum proportion of total light feature intensity that the peak /below/ the heavy monoisotope can have,
    //in order for us /not/ to check for a peak below the heavy monoisotope
    protected float minSignificantPeakContributionBelowMonoisotope = 0.02f;

    protected float maxLightHeavyOverlapToIgnore = 0.02f;

    protected float maxKlDiff = 0.15f;
    protected float minKlRatio = 0.7f;

    //Maximum allowed difference between log ratio we calculate (simply, by most-intense peaks) and algorithm
    protected float maxLogRatioDiff = (float) (Math.log(1.5) - Math.log(1));

    //Ratios must be higher than minFlagRatio, OR lower than maxFlagRatio, to be flagged
    protected float minFlagRatio = 0f;
    protected float maxFlagRatio = 999f;


    //scaffolding for calculating ratios using regression based on one datapoint per-scan, like RelEx.
    //I think that method pretty much doesn't work very well.
//    List<Float> singlePeakSlopeRatios = new ArrayList<Float>();
//    List<Float> multiPeakSlopeRatios = new ArrayList<Float>();

    public QuantEventAssessor()
    {
    }

    /**
     * Assess this feature as a quantitative event
     * @param feature
     * @param run
     * @return
     */
    public QuantEventAssessment assessFeature(Feature feature, MSRun run)
    {
        return assessQuantEvent(new QuantEvent(feature, ""), run);
    }

    /**
     * Assess this quantitative event.  Return assessment.  Side effect: set event.algorithmicAssessment
     * @param event
     * @param run
     * @return
     */
    public QuantEventAssessment assessQuantEvent(QuantEvent event, MSRun run)
    {
        float ratio = event.getRatio();

        if (ratio < getMinFlagRatio() && ratio > getMaxFlagRatio())
        {
            _log.debug("Skipping ratio " + ratio);
            return new QuantEventAssessment(FLAG_REASON_UNEVALUATED, "Not evaluated: ratio near parity");
        }

        int reason = FLAG_REASON_OK;
        String reasonDesc = "";

        float lightMass = Math.max((event.getLightMz() - Spectrum.HYDROGEN_ION_MASS) * event.getCharge(), 0.0f);
        float heavyMass = Math.max((event.getHeavyMz() - Spectrum.HYDROGEN_ION_MASS) * event.getCharge(), 0.0f);
        QuantPeakSetSummary lightPeaksSummary = calcPeakIntensities(event.getFirstHeavyQuantScan(),
                event.getLastHeavyQuantScan(), lightMass, event.getLightMz(), event.getCharge(),
                run, numPeaksToUse);
        QuantPeakSetSummary heavyPeaksSummary = calcPeakIntensities(event.getFirstHeavyQuantScan(),
                event.getLastHeavyQuantScan(), heavyMass, event.getHeavyMz(), event.getCharge(),
                run, numPeaksToUse);

        float intensityBelowLight = lightPeaksSummary.sumIntensityPeakBelow;
        float intensityBelowHeavy = heavyPeaksSummary.sumIntensityPeakBelow;

        List<Float> lightPeakIntensities = lightPeaksSummary.peakSumIntensities;
        List<Float> heavyPeakIntensities = heavyPeaksSummary.peakSumIntensities;

        _log.debug("**light, " + lightPeakIntensities.get(0) + ", " + lightPeakIntensities.get(1) + ", " +
                lightPeakIntensities.get(2) + ", " + lightPeakIntensities.get(3));
        _log.debug("**heavy, " + heavyPeakIntensities.get(0) + ", " + heavyPeakIntensities.get(1) + ", " +
                heavyPeakIntensities.get(2) + ", " + heavyPeakIntensities.get(3));
        _log.debug("**light KL2: " + calcKL(lightMass, lightPeakIntensities, 2) + ", heavy KL2: " +
                    calcKL(heavyMass, heavyPeakIntensities, 2));


        int numPeaksSeparation = PeakOverlapCorrection.calcNumPeaksSeparation(lightPeaksSummary.monoisotopicMass,
                heavyPeaksSummary.monoisotopicMass);
                
        _log.debug("Light mass: " + lightMass + ", Heavy mass " + heavyMass + ", num peaks separation: " + numPeaksSeparation);

        double[] lightIsotopicDistribution = PeakOverlapCorrection.getIsotopicDistribution(
                lightPeaksSummary.monoisotopicMass, numPeaksSeparation+1);
        boolean lightHasSignificantPeakBelowHeavy = lightIsotopicDistribution[numPeaksSeparation-1] >
                minSignificantPeakContributionBelowMonoisotope;
//for (int i=0; i<numPeaksSeparation; i++) System.err.println("Light peak " + i + ", " + lightIsotopicDistribution[i]);        
        _log.debug("Light intrusion on heavy? " + lightHasSignificantPeakBelowHeavy + ", last light contribution: " +
                lightIsotopicDistribution[numPeaksSeparation-1] + ", cap: "  +minSignificantPeakContributionBelowMonoisotope);
        boolean lightHasSignificantHeavyOverlap = lightIsotopicDistribution[numPeaksSeparation] >
                maxLightHeavyOverlapToIgnore;

        //See if the intensity of the peak 1Da below either monoisotope is high enough to worry about.
        if (reason == FLAG_REASON_OK)
        {
            float belowIntensityRatioLight = intensityBelowLight / lightPeakIntensities.get(0);
            //Only evaluate heavy if if light/heavy peaks not overlapping
            float belowIntensityRatioHeavy = lightHasSignificantPeakBelowHeavy ? 0 :
                    intensityBelowHeavy / heavyPeakIntensities.get(0);

            _log.debug("BELOW: light=" + intensityBelowLight + ", ratio="+belowIntensityRatioLight +
                    ", heavy=" +intensityBelowHeavy + ", ratio=" +  belowIntensityRatioHeavy);
            if (Math.max(belowIntensityRatioLight, belowIntensityRatioHeavy) > peakBelowIntensityRatioCutoff)
            {
                reason = FLAG_REASON_COELUTING;
                reasonDesc = "Coeluting intensity ratio light=" + Rounder.round(belowIntensityRatioLight,3) +
                        " heavy=" + Rounder.round(belowIntensityRatioHeavy,3);
            }
        }

        float algRatio = event.getRatio();
        float singlePeakRatio = -1f;

        //Calculate a ratio from the theoretically most-intense peak, compare with algorithm ratio
        if (reason == FLAG_REASON_OK)
        {
            int highestPeakIndex = Spectrum.calcMaxIdealPeakIndex(lightPeaksSummary.monoisotopicMass);

            float lightAreaHighestPeak = lightPeakIntensities.get(highestPeakIndex);
            float heavyAreaHighestPeak = heavyPeakIntensities.get(highestPeakIndex);

            singlePeakRatio = lightAreaHighestPeak / heavyAreaHighestPeak;

            int lightIndexOfHeavyHighestPeak =
                    numPeaksSeparation + highestPeakIndex;
            float lightPeakIntrusionHeavyHighest = (float) PeakOverlapCorrection.getIsotopicDistribution(
                            lightPeaksSummary.monoisotopicMass, lightIndexOfHeavyHighestPeak+1)[lightIndexOfHeavyHighestPeak];
            if (lightPeakIntrusionHeavyHighest >= maxLightHeavyOverlapToIgnore)
            {
                singlePeakRatio = (float) PeakOverlapCorrection.correctRatioForOverlap(
                    lightPeaksSummary.monoisotopicMass, heavyPeaksSummary.monoisotopicMass,
                    singlePeakRatio, highestPeakIndex, highestPeakIndex);
            }
            float logRatioDiff = (float) Math.abs(Math.log(singlePeakRatio) - Math.log(algRatio));

            _log.debug("MS1 Ratio: " + singlePeakRatio + ", alg ratio: " + algRatio + ", log diff: " +
                    Rounder.round(logRatioDiff,3) + ", to beat: " + maxLogRatioDiff);

            if (logRatioDiff > maxLogRatioDiff)
            {
                reason = FLAG_REASON_DISSIMILAR_MS1_RATIO;
                reasonDesc = "Singlepeak=" + Rounder.round(singlePeakRatio, 3) + " algorithm=" +
                        Rounder.round(algRatio, 3);
            }
        }

        //check the KL of both "features".  If peaks overlapping, calculate a special ideal distribution for heavy
        if (reason == FLAG_REASON_OK)
        {
            float lightKl2Peaks = calcKL(lightMass, lightPeakIntensities, 2);
            float heavyKl2Peaks = 0f;
            if (lightHasSignificantHeavyOverlap)
            {
                //if there's significant overlap, calculate a new template peak distribution for heavy,
                //based on the ratio that the algorithm gives us
                float[] heavyIdealDist = Spectrum.Poisson(heavyPeaksSummary.monoisotopicMass);
                float[] lightIdealDist = Spectrum.Poisson(lightPeaksSummary.monoisotopicMass);

                int numPeaksOverlap = lightIdealDist.length - numPeaksSeparation;
                //calculate a new sum for the ideal heavy peaks
                float newHeavySum = 1;
                for (int i=0; i<numPeaksOverlap; i++)
                {
                    float thisPeakExtra = (lightIdealDist[i+numPeaksSeparation] * algRatio);
                    heavyIdealDist[i] += thisPeakExtra;
                    newHeavySum += thisPeakExtra;
                }
                //make the new heavy distribution sum to 1
                for (int i=0; i<heavyIdealDist.length; i++)
                {
                    heavyIdealDist[i] /= (newHeavySum);
                }
                heavyKl2Peaks = calcKL(heavyIdealDist, heavyPeakIntensities, 2);
            }
            else
                heavyKl2Peaks = calcKL(heavyPeaksSummary.monoisotopicMass, heavyPeakIntensities, 2);
            float klDiff = Math.abs(lightKl2Peaks - heavyKl2Peaks);
            float klRatio = lightKl2Peaks / heavyKl2Peaks;
            _log.debug("Light KL: " + lightKl2Peaks + ", Heavy KL: " + heavyKl2Peaks + ": diff=" + klDiff + ", ratio=" + klRatio);
            if (klDiff > maxKlDiff && klRatio < minKlRatio)
            {
                reason = FLAG_REASON_DISSIMILAR_KL;
                reasonDesc = "light KL=" + Rounder.round(lightKl2Peaks,3) + " heavy=" + Rounder.round(heavyKl2Peaks,3) +
                        " diff=" + Rounder.round(klDiff,3) + " ratio=" + Rounder.round(klRatio,3);
            }
            else if (lightKl2Peaks > max2PeakKL || heavyKl2Peaks > max2PeakKL)
            {
                reason = FLAG_REASON_BIG_2PEAK_KL;
                reasonDesc = "";
                if (lightKl2Peaks > max2PeakKL)
                    reasonDesc = reasonDesc + "light KL=" + lightKl2Peaks + " ";
                if (heavyKl2Peaks > max2PeakKL)
                    reasonDesc = reasonDesc + "heavy KL=" + heavyKl2Peaks;                      
            }
        }

        _log.debug("REASON: " + flagReasonDescriptions[reason]);

        QuantEventAssessment result = new QuantEventAssessment(reason, reasonDesc);
        result.setSinglePeakRatio(singlePeakRatio);
        event.setAlgorithmicAssessment(result);
        return result;
    }

    /**
     * Calculate a KL score for a list of peak intensities.  Have to normalize intensity first by dividing by peak sum
     * @param mass
     * @param peakIntensities
     * @param numPeaksToUse
     * @return
     */
    protected float calcKL(float mass, List<Float> peakIntensities, int numPeaksToUse)
    {
        return calcKL(Spectrum.Poisson(mass), peakIntensities, numPeaksToUse);
    }


    protected float calcKL(float[] idealPeaks, List<Float> peakIntensities, int numPeaksToUse)
    {             
        float[] peakIntensitiesXPeaks = new float[numPeaksToUse];
        float[] idealPeaksXPeaks = new float[numPeaksToUse];

        float sumRealPeaks = 0;
        float sumIdealPeaks = 0;
        for (int i=0; i<numPeaksToUse; i++)
        {
            if (i < peakIntensities.size())
                peakIntensitiesXPeaks[i] = peakIntensities.get(i);
            peakIntensitiesXPeaks[i] = Math.max(0.1f, peakIntensitiesXPeaks[i]);
            sumRealPeaks += peakIntensitiesXPeaks[i];

            sumIdealPeaks += idealPeaks[i];
            idealPeaksXPeaks[i] = idealPeaks[i];
        }
		for (int i = 0; i < numPeaksToUse; i++)
        {
            peakIntensitiesXPeaks[i] /= sumRealPeaks;
            idealPeaksXPeaks[i] /= sumIdealPeaks;
//System.err.println(peakIntensities6Peaks[i]);
        }
        return PeakOverlapCorrection.calcKLUsingTemplate(idealPeaksXPeaks, peakIntensitiesXPeaks);
    }


    /**
     * Return the peak quant summary, with intensities of all peaks.
     * Use raw intensities, NOT resampled
     * @param firstScan
     * @param lastScan
     * @param mass
     * @param mz
     * @param charge
     * @param run
     * @param numPeaks
     * @return
     */
    protected QuantPeakSetSummary calcPeakIntensities(int firstScan, int lastScan, float mass, float mz,
                                                      float charge, MSRun run,
                                                      int numPeaks)
    {
        //Define the scan range, making sure not to go out of bounds
        //assumes heavy and light scan extents same
        int firstScanIndex = Math.max(Math.abs(
                run.getIndexForScanNum(firstScan)) -
                                       numScansAroundEventToConsider, 0);
        int lastScanIndex = Math.min(Math.abs(
                run.getIndexForScanNum(lastScan)) +
                                       numScansAroundEventToConsider, run.getScanCount()-1);
        lastScanIndex = Math.max(firstScanIndex, lastScanIndex);

        Scan[] scans = FeatureFinder.getScans(run, firstScanIndex,
                lastScanIndex - firstScanIndex + 1);

        float mzTol = (MassUtilities.calculateAbsoluteDeltaMass(
                mass, peakPPMTolerance, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / charge;
        QuantPeakSetSummary result = new QuantPeakSetSummary();
        result.monoisotopicMass = mass;
        result.scanRetentionTimes = new double[scans.length];
        for (int i=0; i<scans.length; i++)
        {
            result.scanRetentionTimes[i] = i;
        }
        result.peakSumIntensities = new ArrayList<Float>();
        for (int i=0; i<numPeaks; i++) result.peakSumIntensities.add(0f);

        for (int scanIndex = 0; scanIndex < scans.length; scanIndex++)
        {
            float[][] spectrum = scans[scanIndex].getSpectrum();

            for (int peakIndex=-1; peakIndex<numPeaks; peakIndex++)
            {
                float peakMzCenter = mz +
                        ((peakIndex * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH) /
                                charge);

                int startIndex = Arrays.binarySearch(spectrum[0], peakMzCenter - mzTol);
                startIndex = startIndex < 0 ? -(startIndex+1) : startIndex;
                int endIndex = Arrays.binarySearch(spectrum[0], peakMzCenter + mzTol);
                endIndex = endIndex < 0 ? -(endIndex+1) : endIndex;

                float intensityMax = 0;
                for (int i=startIndex; i<=endIndex; i++)
                    intensityMax = Math.max(intensityMax, spectrum[1][i]);

                if (peakIndex == -1)
                    result.sumIntensityPeakBelow += intensityMax;
                else
                    result.peakSumIntensities.set(peakIndex, result.peakSumIntensities.get(peakIndex) + intensityMax);

            }
        }
        return result;
    }


    public float getPeakPPMTolerance()
    {
        return peakPPMTolerance;
    }

    public void setPeakPPMTolerance(float peakPPMTolerance)
    {
        this.peakPPMTolerance = peakPPMTolerance;
    }

    public static Logger getLogger()
    {
        return _log;
    }

    public float getMinFlagRatio()
    {
        return minFlagRatio;
    }

    public void setMinFlagRatio(float minFlagRatio)
    {
        this.minFlagRatio = minFlagRatio;
    }

    public float getMaxFlagRatio()
    {
        return maxFlagRatio;
    }

    public void setMaxFlagRatio(float maxFlagRatio)
    {
        this.maxFlagRatio = maxFlagRatio;
    }

    public boolean isShowCharts()
    {
        return showCharts;
    }

    public void setShowCharts(boolean showCharts)
    {
        this.showCharts = showCharts;
    }

    public int getLabelType()
    {
        return labelType;
    }

    public void setLabelType(int labelType)
    {
        this.labelType = labelType;
    }

    /**
     * Data structure to pass around summary information about one set of peaks (light or heavy)
     */
    public static class QuantPeakSetSummary
    {
        protected float monoisotopicMass;
        protected float sumIntensityPeakBelow;
        protected List<Float> peakSumIntensities;
        protected double[] scanRetentionTimes;

//        protected List<double[]> rawIntensities;

        public QuantPeakSetSummary()
        {
            peakSumIntensities = new ArrayList<Float>();
        }

//        double[] combineRawIntensitiesAllPeaks()
//        {
//            int fullLength = 0;
//            for (double[] array : rawIntensities)
//                fullLength += array.length;
//            double[] result = new double[fullLength];
//            int position=0;
//            for (double[] array : rawIntensities)
//            {
//                System.arraycopy(array, 0, result, position, array.length);
//                position += array.length;
//            }
//            return result;
//        }
    }

    /**
     * Class for communicating assessment results
     */
    public static class QuantEventAssessment
    {
        protected int status = -1;
        protected String explanation;
        protected float singlePeakRatio;

        public String toString()
        {
            return "Assessment.  Status: " +  getAssessmentCodeDesc(status) + ", single-peak: " +
                    singlePeakRatio + ", explanation: " + explanation;
        }

        public QuantEventAssessment(int status, String explanation)
        {
            this.status = status;
            this.explanation = explanation;
        }

        public int getStatus()
        {
            return status;
        }

        public void setStatus(int status)
        {
            this.status = status;
        }

        public String getExplanation()
        {
            return explanation;
        }

        public void setExplanation(String explanation)
        {
            this.explanation = explanation;
        }

        public float getSinglePeakRatio()
        {
            return singlePeakRatio;
        }

        public void setSinglePeakRatio(float singlePeakRatio)
        {
            this.singlePeakRatio = singlePeakRatio;
        }
    }    



}
