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
package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

import java.util.*;
import java.util.List;


/**
 * This class contains the code for deciding whether an isotopic label quantitation event is good or not.
 * It's aimed at running-time efficiency: there are several tests for different kinds of badness, and they're run
 * roughly in order of increasing time investment.  Therefore, all that's retained is a good-bad call, and if bad, an 
 * indication of /one/ of the (potentially several) things that is bad.
 *
 * This uses a big HACK, adding a dummy search score to peptide features to store the flag reason description.
 */
public class QuantEventAssessor
{
    public static Logger _log = Logger.getLogger(QuantEventAssessor.class);

    public static final int FLAG_REASON_OK = 0;
    public static final int FLAG_REASON_COELUTING = 1;
    public static final int FLAG_REASON_DISSIMILAR_KL = 2;
    public static final int FLAG_REASON_DISSIMILAR_MS1_RATIO = 3;
    public static final int FLAG_REASON_BIG_2PEAK_KL = 4;
    public static final int FLAG_REASON_MISSING_PEAKS = 5;
    public static final int FLAG_REASON_UNEVALUATED = 6;
    public static final int FLAG_REASON_OTHER = 7;

    public static final String[] flagReasonDescriptions = new String[]
            {
                    "OK",
                    "Coeluting peptide",
                    "Dissimilar light/heavy KL",
                    "Singlepeak Ratio Different",
                    "Big 2-peak KL",
                    "Missing peaks",
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
                    "MissingPeaks",
                    "Unevaluated",
                    "Other",
            };

    //boundaries for the one-peak ratio.  If it's outside these boundaries, snap it to them.  This is artificial,
    //of course, but it's necessary to prevent ratios in the thousands (or more!) from drowning out all other
    //ratios for a protein
    public static final float MIN_RATIO_ONE_PEAK = 0.04f;
    public static final float MAX_RATIO_ONE_PEAK = 25f;

    //Default high and low "extreme" ratios.
    public static final float DEFAULT_EXTREME_RATIO_HIGH = 2f;
    public static final float DEFAULT_EXTREME_RATIO_LOW = 1f / DEFAULT_EXTREME_RATIO_HIGH;

    //If an event has an Extreme ratio, and a problem with it could only
    //mean that the correct value is more Extreme, then we call the event OK.
    protected float extremeRatioHigh = DEFAULT_EXTREME_RATIO_HIGH;
    protected float extremeRatioLow = DEFAULT_EXTREME_RATIO_LOW;


    //Should we perform all the checks, for bad events (for capturing all info), or stop at the first one (for performance)
    protected boolean shouldPerformAllChecks = true;

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
        if (flagReasonCodes[FLAG_REASON_DISSIMILAR_MS1_RATIO].equals(curationStatusString))
            return FLAG_REASON_DISSIMILAR_MS1_RATIO;
        if (flagReasonCodes[FLAG_REASON_BIG_2PEAK_KL].equals(curationStatusString))
            return FLAG_REASON_BIG_2PEAK_KL;
        if (flagReasonCodes[FLAG_REASON_OTHER].equals(curationStatusString))
            return FLAG_REASON_OTHER;
        if (flagReasonCodes[FLAG_REASON_MISSING_PEAKS].equals(curationStatusString))
            return FLAG_REASON_MISSING_PEAKS;
        return FLAG_REASON_UNEVALUATED;
    }

    public static final String REASON_DUMMY_SEARCH_SCORE_NAME = "dummy_flag_desc";

    //Maximum allowed proportion of highest peak intensity that the peak 1Da below monoisotope can have
    //dhmay changing from 0.3 to 0.35, 20091029.  Was too sensitive
    public static final float DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF = 0.35f;
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

//    private int labelType = QuantitationUtilities.LABEL_LYCINE;

    protected boolean showCharts = false;


//    public static final float SILAC_LABEL_MASS = 134.115092f;
//    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;
//
//
//    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
//    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
//    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
//            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;


    protected int numPeaksToUse = 5;
    protected int numRawPeaksToKeep = 2;

    //for coeluting peptide check, the maximum proportion (in the theoretical distribution) of the below-monoisotope
    //peak within the whole isotopic distribution.  Want to make sure that this is distinguishable from a 1-peak wonder
    float maxCoelutingPeptideContributionOfBelowMonoPeak = 0.5f;
    //Maximum KL score of a coeluting peptide in another charge, to consider it an adequate explanation for a
    //below-monoisotope peak.
    //This value has a big effect on false positives
    float maxKLCoelutingOtherChargePeptide = 1f;

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

    //Should events be flagged if there are peaks missing and they wouldn't otherwise be flagged
    protected boolean shouldFlagIfMissingPeaks = false;
    //Should events be flagged if KL scores for light and heavy are different?
    protected boolean shouldFlagDifferentKL = true;


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
        _log.debug("**light KL2: " + scaleAndCalcKL(lightMass, lightPeakIntensities, 2) + ", heavy KL2: " +
                    scaleAndCalcKL(heavyMass, heavyPeakIntensities, 2));


        int numPeaksSeparation = PeakOverlapCorrection.calcNumPeaksSeparation(lightPeaksSummary.monoisotopicMass,
                heavyPeaksSummary.monoisotopicMass);
                
        _log.debug("Light mass: " + lightMass + ", Heavy mass " + heavyMass + ", num peaks separation: " + numPeaksSeparation);

        int highestPeakIndex = Spectrum.calcMaxIdealPeakIndex(lightPeaksSummary.monoisotopicMass);
        float lightAreaHighestPeak = lightPeakIntensities.get(highestPeakIndex);
        float heavyAreaHighestPeak = heavyPeakIntensities.get(highestPeakIndex);

        float singlePeakRatio = lightAreaHighestPeak / heavyAreaHighestPeak;
        //dhmay adding 20091027.  Set boundaries on the single-peak ratio, so we don't create such an extreme ratio
        //that it dominates the geometric mean calculation
        singlePeakRatio = Math.max(MIN_RATIO_ONE_PEAK, Math.min(singlePeakRatio, MAX_RATIO_ONE_PEAK));
        float algRatio = event.getRatio();

        double[] lightIsotopicDistribution = PeakOverlapCorrection.getIsotopicDistribution(
                lightPeaksSummary.monoisotopicMass, numPeaksSeparation+1);
        double[] heavyIsotopicDistribution = PeakOverlapCorrection.getIsotopicDistribution(
                heavyPeaksSummary.monoisotopicMass, numPeaksSeparation+1);
        //Determine whether there's a significant light peak below the heavy mono carefully: use whichever light-heavy
        //ratio is higher, algorithm or single-peak
        float lastLightPeakContributionToHeavyMono =
                (float) lightIsotopicDistribution[numPeaksSeparation-1] * Math.max(singlePeakRatio, algRatio) / (float) heavyIsotopicDistribution[0];
        boolean lightHasSignificantPeakBelowHeavy =  lastLightPeakContributionToHeavyMono >
                minSignificantPeakContributionBelowMonoisotope;
//for (int i=0; i<numPeaksSeparation; i++) System.err.println("Light peak " + i + ", " + lightIsotopicDistribution[i]);


        _log.debug("Light intrusion on heavy? " + lightHasSignificantPeakBelowHeavy + ", last light contribution: " +
                lastLightPeakContributionToHeavyMono + ", cap: " + minSignificantPeakContributionBelowMonoisotope);
        //Determine whether there's significant light-heavy distribution overlap using algorithm ratio
        boolean lightHasSignificantHeavyOverlap =
                lightIsotopicDistribution[numPeaksSeparation] * algRatio / heavyIsotopicDistribution[0] >
                        maxLightHeavyOverlapToIgnore;

        QuantEventAssessment result = new QuantEventAssessment(FLAG_REASON_OK, "");


        //dhmay adding check for 0 light or heavy area in BOTH of the first two peaks of either light or heavy
        if ((result.isGood() || shouldPerformAllChecks) && shouldFlagIfMissingPeaks)
        {
            if (lightPeakIntensities.get(0) == 0 && lightPeakIntensities.get(1) == 0)
            {
                result.setFlag(FLAG_REASON_MISSING_PEAKS, true);
                result.setCheckExplanation(FLAG_REASON_MISSING_PEAKS, "Missing first two light peaks");
            }
            else if (heavyPeakIntensities.get(0) == 0 && heavyPeakIntensities.get(1) == 0)
            {
                result.setFlag(FLAG_REASON_MISSING_PEAKS, true);
                result.setCheckExplanation(FLAG_REASON_MISSING_PEAKS, "Missing first two heavy peaks");
            }
        }

        //Calculate a ratio from the theoretically most-intense peak, compare with algorithm ratio
        if (result.isGood() ||shouldPerformAllChecks)
        {
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
                if ((algRatio >= extremeRatioHigh && singlePeakRatio >= extremeRatioHigh) ||
                        (algRatio <= extremeRatioLow && singlePeakRatio <= extremeRatioLow))
                {
                    _log.debug("Giving MS1 ratio diff a pass because of extreme ratio: algorithm: " + algRatio +
                            ", MS1: " + singlePeakRatio);
                }
                else
                {
                    //dhmay added check for extreme ratios about which we don't care, 20091029
                result.setFlag(FLAG_REASON_DISSIMILAR_MS1_RATIO, true);
                result.setCheckExplanation(FLAG_REASON_DISSIMILAR_MS1_RATIO, "Singlepeak=" +
                        Rounder.round(singlePeakRatio, 3) + " algorithm=" +
                            Rounder.round(algRatio, 3));
                }
            }
        }
//if (reason == FLAG_REASON_DISSIMILAR_MS1_RATIO)  System.err.println("****Scan " + event.getScan() + ", reasonDesc=" + reasonDesc);      

        float lightKl2Peaks = 0f;
        float heavyKl2Peaks = 0f;

        //check the KL of both "features".  If peaks overlapping, calculate a special ideal distribution for heavy
        if (result.isGood() ||shouldPerformAllChecks)
        {
            lightKl2Peaks = scaleAndCalcKL(lightMass, lightPeakIntensities, 2);
            heavyKl2Peaks = 0f;
            if (lightHasSignificantHeavyOverlap)
            {
                //if there's significant overlap, calculate a new template peak distribution for heavy,
                //based on the ratio that the algorithm gives us
                float[] heavyIdealDistNoOverwrite = Spectrum.Poisson(heavyPeaksSummary.monoisotopicMass);
                float[] heavyIdealDist = new float[heavyIdealDistNoOverwrite.length];
                System.arraycopy(heavyIdealDistNoOverwrite, 0, heavyIdealDist, 0, heavyIdealDist.length);

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
                heavyKl2Peaks = scaleAndCalcKL(heavyIdealDist, heavyPeakIntensities, 2);
            }
            else
                heavyKl2Peaks = scaleAndCalcKL(heavyPeaksSummary.monoisotopicMass, heavyPeakIntensities, 2);
            float klDiff = Math.abs(lightKl2Peaks - heavyKl2Peaks);
            float klRatio = lightKl2Peaks / heavyKl2Peaks;
            _log.debug("Light KL: " + lightKl2Peaks + ", Heavy KL: " + heavyKl2Peaks + ": diff=" + klDiff + ", ratio=" + klRatio);
            if (klDiff > maxKlDiff && klRatio < minKlRatio && shouldFlagDifferentKL)
            {
                //The KL scores are significantly different.  But if it's a very extreme ratio and the worse score is on
                //the /less/ intense one, we don't care (it doesn't matter if there's something coeluting)
                if ((algRatio < (extremeRatioLow*2) && lightKl2Peaks > heavyKl2Peaks) ||
                    (algRatio > (extremeRatioHigh*2) && heavyKl2Peaks > lightKl2Peaks))
                {
                    result.setCheckExplanation(FLAG_REASON_DISSIMILAR_KL, "KL scores different, but worse one is low-abundance");
                }
                else
                {
                    //OK, no excuses, KL scores are just different.  Flag.
                    result.setFlag(FLAG_REASON_DISSIMILAR_KL, true);

                    result.setCheckExplanation(FLAG_REASON_DISSIMILAR_KL, "light KL=" + Rounder.round(lightKl2Peaks,3) + " heavy=" + Rounder.round(heavyKl2Peaks,3) +
                            " diff=" + Rounder.round(klDiff,3) + " ratio=" + Rounder.round(klRatio,3));
                }
            }
            if ((result.isGood() || shouldPerformAllChecks) && (lightKl2Peaks > max2PeakKL || heavyKl2Peaks > max2PeakKL))

            {
                String explanation = "";
                result.setFlag(FLAG_REASON_BIG_2PEAK_KL, true);
                if (lightKl2Peaks > max2PeakKL)
                    explanation = explanation + "light KL=" + lightKl2Peaks + " ";
                if (heavyKl2Peaks > max2PeakKL)
                    explanation = explanation + "heavy KL=" + heavyKl2Peaks;
                result.setCheckExplanation(FLAG_REASON_BIG_2PEAK_KL, explanation);

            }
        }

        //dhmay moving this check down to be last, since it potentially involves feature-finding
        //See if the intensity of the peak 1Da below either monoisotope is high enough to worry about.
        if (result.isGood() || shouldPerformAllChecks)
        {
            Pair<Boolean, String> statusAndDesc = checkCoelutingPeptides(intensityBelowLight, lightPeakIntensities,
                    intensityBelowHeavy, heavyPeakIntensities,
                    lightMass, event.getLightMz(), heavyMass, event.getHeavyMz(), event.getCharge(), run, 
                    event.getFirstHeavyQuantScan(), event.getLastHeavyQuantScan(), lightHasSignificantPeakBelowHeavy,
                    algRatio, event.getScan(), lightKl2Peaks, heavyKl2Peaks);
            result.setCheckExplanation(FLAG_REASON_COELUTING, statusAndDesc.second);
            result.setFlag(FLAG_REASON_COELUTING, statusAndDesc.first);
        }

        _log.debug(result.toString());

        result.setSinglePeakRatio(singlePeakRatio);
        for (int i=0; i<flagReasonCodes.length; i++)
        {
            if (i != FLAG_REASON_OK && result.getFlagBitmap()[i])
            {
                result.setBad();
                break;
            }
        }
        event.setAlgorithmicAssessment(result);
        
        return result;
    }


    /**
     * Check for coeluting peptides that might explain a high below-monoisotope peak.
     * Pulling this out as a separate method because it got kind of complicated.
     * @param intensityBelowLight
     * @param lightPeakIntensities
     * @param intensityBelowHeavy
     * @param heavyPeakIntensities
     * @param lightMass
     * @param lightMz
     * @param heavyMass
     * @param heavyMz
     * @param charge
     * @param run
     * @param firstScan
     * @param lastScan
     * @param lightHasSignificantPeakBelowHeavy
     * @param algRatio
     * @param scan
     * @param lightKl2Peaks
     * @param heavyKl2Peaks
     * @return
     */
    protected Pair<Boolean, String> checkCoelutingPeptides(float intensityBelowLight, List<Float> lightPeakIntensities,
                                                           float intensityBelowHeavy, List<Float> heavyPeakIntensities,
                                                           float lightMass, float lightMz, float heavyMass, float heavyMz,
                                                           int charge, MSRun run,
                                                           int firstScan, int lastScan, boolean lightHasSignificantPeakBelowHeavy,
                                                           float algRatio, int scan, float lightKl2Peaks, float heavyKl2Peaks)
    {
        boolean resultIsBad = false;
        String resultDesc = "";
        float belowIntensityRatioLight = intensityBelowLight / lightPeakIntensities.get(0);
        //Only evaluate heavy if if light/heavy peaks not overlapping
        float belowIntensityRatioHeavy = lightHasSignificantPeakBelowHeavy ? 0 :
                intensityBelowHeavy / heavyPeakIntensities.get(0);

        _log.debug("BELOW: light=" + intensityBelowLight + ", ratio="+belowIntensityRatioLight +
                ", heavy=" +intensityBelowHeavy + ", ratio=" +  belowIntensityRatioHeavy);
        if (Math.max(belowIntensityRatioLight, belowIntensityRatioHeavy) > peakBelowIntensityRatioCutoff)
        {
            //extra check to ignore coeluting peptides that could only make an extreme
            //ratio more extreme
            boolean correctedRatioCouldOnlyGetSmaller = belowIntensityRatioHeavy < peakBelowIntensityRatioCutoff;
            boolean correctedRatioCouldOnlyGetBigger = belowIntensityRatioLight < peakBelowIntensityRatioCutoff;
            if (!(correctedRatioCouldOnlyGetSmaller && algRatio <= extremeRatioLow)    &&
                    !(correctedRatioCouldOnlyGetBigger && algRatio >= extremeRatioHigh))
            {

//_log.setLevel(Level.DEBUG);
                _log.debug("\n**Coelute, SCAN " + scan + ", charge=" + charge + ", lightmz: " + lightMz + ", heavymz: " + heavyMz + ", light ratio: " + belowIntensityRatioLight + ", heavy ratio: " + belowIntensityRatioHeavy);

                //OK, we've got a potentially problematic event due to coeluting peptide.
                //Now we check to see if that peak below the monoisotope is actually a problem.
                boolean lightCheckedAndOK = false;
                if (belowIntensityRatioLight >= peakBelowIntensityRatioCutoff)
                {
                    _log.debug("\tCheck light");
                    if ((lightKl2Peaks < 1f) && (scalePeaksSum1(lightPeakIntensities).get(0) >= 0.1) &&
                            checkPeakBelowExplainedByOtherChargeFeature(lightMz,
                                    charge, run, firstScan, lastScan))
                    {
                        resultDesc = "Coeluting peak below light, but explained by other feature";
                        lightCheckedAndOK = true;
                    }
                    else
                    {
                        resultIsBad = true;
                        resultDesc = "Coeluting intensity ratio LIGHT=" + Rounder.round(belowIntensityRatioLight,3) +
                                        " heavy=" + Rounder.round(belowIntensityRatioHeavy,3);
                    }
                }

                if ((resultIsBad || shouldPerformAllChecks) && belowIntensityRatioHeavy >= peakBelowIntensityRatioCutoff)
                {
                    _log.debug("\tCheck heavy");

                    if ((heavyKl2Peaks < 1f) && (scalePeaksSum1(heavyPeakIntensities).get(0) >= 0.1) &&
                            checkPeakBelowExplainedByOtherChargeFeature(heavyMz,
                                    charge, run, firstScan, lastScan))
                    {
                        if (lightCheckedAndOK)
                            resultDesc = "Coeluting peak below light and heavy, but explained by other features";
                        else
                            resultDesc = "Coeluting peak below heavy, but explained by other feature";
                    }
                    else
                    {
                        resultIsBad = true;
                        resultDesc = "Coeluting intensity ratio light=" + Rounder.round(belowIntensityRatioLight,3) +
                                        " HEAVY=" + Rounder.round(belowIntensityRatioHeavy,3);
                    }
                }
//_log.setLevel(Level.INFO);
            }
            else
                _log.debug("Giving coeluting peptide a pass because of extreme ratio: " + algRatio +
                        ", could only get smaller? " + correctedRatioCouldOnlyGetSmaller + ", bigger? " + correctedRatioCouldOnlyGetBigger);
        }
        return new Pair<Boolean, String>(resultIsBad, resultDesc);
    }

    /**
     * Scale a set of peak intensities so they sum to 1
     * @param peakIntensities
     * @return
     */
    protected List<Float> scalePeaksSum1(List<Float> peakIntensities)
    {
        List<Float> result = new ArrayList<Float>();
        float peakSum = 0f;
        for (float peakIntensity : peakIntensities)
            peakSum += peakIntensity;
        for (float peakIntensity : peakIntensities)
            result.add(peakIntensity / peakSum);
        return result;
    }

    /**
     * Check for a peptide from another charge that might explain a big below-monoisotope peak.  If the peptide's
     * from another charge, then the peaks might not overlap our peptide's peaks.   Right now, only handling charge
     * 2 vs. 3.
     *
     * Require at least 2 peaks in the clear and that the peak below the monoisotope be at most proportion
     * maxCoelutingPeptideContributionOfBelowMonoPeak of the total distribution (to ignore one-peak wonders, which are
     * indistinguishable from an actual coeluting peptide in the same charge), and that the KL of the other-charge
     * feature be low
     * @param mz
     * @param charge
     * @param run
     * @param firstScan
     * @param lastScan
     * @return
     */
    protected boolean checkPeakBelowExplainedByOtherChargeFeature(float mz,
                                                                  int charge, MSRun run, int firstScan, int lastScan)
    {
        float peakSeparationMass = (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH;

        List<Float> monoMassesToCheckPeaks = new ArrayList<Float>();
        List<Integer> numbersOfPeaksToCheck = new ArrayList<Integer>();
        float mzPeakBelow = mz - (peakSeparationMass / (float) charge);
        float[] idealPeaksThisMass = null;
        float massPeakBelow = 0f;

        int otherPeptideCharge = 0;
        switch (charge)
        {
            case 2:
                otherPeptideCharge = 3;
                massPeakBelow = MassUtilities.calcMassForMzAndCharge(mzPeakBelow, otherPeptideCharge);
                idealPeaksThisMass = Spectrum.Poisson(massPeakBelow);
                if (idealPeaksThisMass[1] < maxCoelutingPeptideContributionOfBelowMonoPeak)
                {
                    _log.debug("\tcharge 2, type 1");
                    monoMassesToCheckPeaks.add(massPeakBelow-peakSeparationMass);
                    numbersOfPeaksToCheck.add(3);
                }
                if (idealPeaksThisMass[0] < maxCoelutingPeptideContributionOfBelowMonoPeak)
                {
                    _log.debug("\tcharge 2, type 2");
                    monoMassesToCheckPeaks.add(massPeakBelow);
                    numbersOfPeaksToCheck.add(3);
                }
                break;
            case 3:
                otherPeptideCharge = 2;
                massPeakBelow = MassUtilities.calcMassForMzAndCharge(mzPeakBelow, otherPeptideCharge);
                idealPeaksThisMass = Spectrum.Poisson(massPeakBelow);                                
                if (idealPeaksThisMass[1] < maxCoelutingPeptideContributionOfBelowMonoPeak)
                {
                    _log.debug("\tcharge 3, type 1");
                    monoMassesToCheckPeaks.add(massPeakBelow-peakSeparationMass);
                    numbersOfPeaksToCheck.add(3);
                }
                if (idealPeaksThisMass[0] < maxCoelutingPeptideContributionOfBelowMonoPeak)
                {
                    _log.debug("\tcharge 3, type 2");
                    monoMassesToCheckPeaks.add(massPeakBelow);
                    numbersOfPeaksToCheck.add(2);
                }
                break;
            default:
                return false; //todo: handle charge 1, 4
        }

        //todo: make this twice as efficient by pulling out intensities for multiple features at same time
        for (int monoIndex=0; monoIndex<monoMassesToCheckPeaks.size(); monoIndex++)
        {
            float monoMassToCheck = monoMassesToCheckPeaks.get(monoIndex);
            int numPeaksToCheck = numbersOfPeaksToCheck.get(monoIndex);
            _log.debug("\tchecking " + monoMassToCheck);
            float mzTolerance = (MassUtilities.calculateAbsoluteDeltaMass(
                    monoMassToCheck, peakPPMTolerance, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / otherPeptideCharge;
            float[] mzValues = new float[numPeaksToCheck];
            for (int i=0; i<mzValues.length; i++)
            {
                mzValues[i] = MassUtilities.calcMzForMassAndCharge(monoMassToCheck + (i * peakSeparationMass), otherPeptideCharge);
            }
            float[] otherPeptideIntensities = extractPeakSumIntensities(firstScan, lastScan, mzValues, mzTolerance, run);
            if (_log.isDebugEnabled())
                for (int i=0; i<mzValues.length; i++) System.err.println("\t\tmz: " + mzValues[i] + ", int: " + otherPeptideIntensities[i]);
            float otherPeptideKL = scaleAndCalcKL(monoMassToCheck, otherPeptideIntensities);
            _log.debug("\tKL: " + otherPeptideKL);
            if (otherPeptideKL < maxKLCoelutingOtherChargePeptide)
            {
                _log.debug("GOOD!");
                return true;
            }
            else _log.debug("BAD!");
        }

        return false;
    }


    /**
     * NEVER FULLY IMPLEMENTED.  This approach doesn't look promising, because feature-finding works really badly
     * when there are multiple features jumbled together
     * @param intensityBelowLight
     * @param lightPeakIntensities
     * @param intensityBelowHeavy
     * @param heavyPeakIntensities
     * @param lightMass
     * @param lightMz
     * @param heavyMass
     * @param heavyMz
     * @param charge
     * @param run
     * @param numPeaksToUse
     * @param firstScan
     * @param lastScan
     * @return
     */
//    protected int checkCoelutingFeatures(float intensityBelowLight, List<Float> lightPeakIntensities,
//                            float intensityBelowHeavy, List<Float> heavyPeakIntensities,
//                            float lightMass, float lightMz, float heavyMass, float heavyMz,
//                            float charge, MSRun run, int numPeaksToUse,
//                            int firstScan, int lastScan)
//    {
//        //m/z tolerance around each peak, for this charge and mass
//        float peakMassTolerance = MassUtilities.calculateAbsoluteDeltaMass(heavyMass, peakPPMTolerance,
//                FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
//
//        Class featureStrategyClass = FeatureExtractor.getDefaultClass();
//        float featureFindingMzRangeLow = lightMz - (2.0f / charge) - 0.1f;
//        float featureFindingMzRangeHigh = heavyMz +  ((numPeaksToUse + 3) / charge) + 0.1f;
//        FloatRange mzRangeExtract = new FloatRange(featureFindingMzRangeLow, featureFindingMzRangeHigh);
//        int firstScanIndexWithSlop = Math.max(run.getIndexForScanNum(firstScan)-2, 0);
//        int lastScanIndexWithSlop = Math.min(run.getScanCount()-1, run.getIndexForScanNum(lastScan) + 2);
//
//System.err.println(" Range: scan: " + run.getScanNumForIndex(firstScanIndexWithSlop) + ", " + run.getScanNumForIndex(lastScanIndexWithSlop) + ".  mz: " +
//                 featureFindingMzRangeLow + "," + featureFindingMzRangeHigh);
//        FeatureFinder featureFinder =
//                new FeatureFinder(run, firstScanIndexWithSlop, lastScanIndexWithSlop - firstScanIndexWithSlop + 1,
//                        PeakCombiner.DEFAULT_MAX_CHARGE,
//                        mzRangeExtract,
//                        featureStrategyClass, false);
//        Feature[] allRegionFeatures = null;
//        try
//        {
//            allRegionFeatures = featureFinder.findPeptides().getFeatures();
//        }
//        catch (InterruptedException e)
//        {}
//
//        //todo: this is inefficient, hitting all the features in turn.  Should store these in a tree
//        Arrays.sort(allRegionFeatures, new Feature.MassAscComparator());
//        Feature lightFeature = null;
//        Feature heavyFeature = null;
//        boolean pastLightMass = false;
//        boolean pastHeavyMass = false;
//        //assign the most-intense features with m/z of the light and heavy monoisotopes
//        for (Feature feature : allRegionFeatures)
//        {
//            //first, check to see if we're beyond the m/z range we care about for feature monoisotopic masses
//            float featureMass = feature.getMass();
//            if (!pastLightMass)
//                pastLightMass = featureMass > lightMass + peakMassTolerance;
//            pastHeavyMass = featureMass > heavyMass + peakMassTolerance;
//            if (pastHeavyMass)
//                break;
//System.err.println("Mass: " + feature.getMass() + ", mz: " +  feature.getMz() + ", charge: " + feature.getCharge() + ", peaks: " + feature.getPeaks() + ", scan: " + feature.getScan());
//            if (!pastLightMass && Math.abs(featureMass - lightMass) < peakMassTolerance)
//            {
//                if (lightFeature == null || feature.getIntensity() > lightFeature.getIntensity())
//                    lightFeature = feature;
//            }
//            if (pastLightMass && Math.abs(featureMass - heavyMass) < peakMassTolerance)
//            {
//                if (heavyFeature == null || feature.getIntensity() > heavyFeature.getIntensity())
//                    heavyFeature = feature;
//            }
//        }
//System.err.println("FEATURES: " + allRegionFeatures.length + ", LIGHTMASS: " + lightMass + ", HEAVYMASS: " + heavyMass +  ", LIGHT: " + (lightFeature != null) + ", HEAVY: " + (heavyFeature != null));
//
//        return FLAG_REASON_COELUTING;
//    }


    //Utility methods for calculating KL on various things.  Move elsewhere?

    /**
     * Calculate a KL score for a list of peak intensities.  Have to normalize intensity first by dividing by peak sum
     * @param mass
     * @param peakIntensities
     * @return
     */
    protected float scaleAndCalcKL(float mass, List<Float> peakIntensities)
    {
        return scaleAndCalcKL(Spectrum.Poisson(mass), peakIntensities, peakIntensities.size());
    }

    protected float scaleAndCalcKL(float mass, float[] peakIntensities)
    {
        return scaleAndCalcKL(Spectrum.Poisson(mass), peakIntensities, peakIntensities.length);
    }

    /**
     * Calculate a KL score for a list of peak intensities.  Have to normalize intensity first by dividing by peak sum
     * @param mass
     * @param peakIntensities
     * @param numPeaksToUse
     * @return
     */
    protected float scaleAndCalcKL(float mass, List<Float> peakIntensities, int numPeaksToUse)
    {
        return scaleAndCalcKL(Spectrum.Poisson(mass), peakIntensities, numPeaksToUse);
    }


    protected float scaleAndCalcKL(float[] idealPeaks, List<Float> peakIntensities, int numPeaksToUse)
    {
        float[] peakIntensitiesFloat = new float[peakIntensities.size()];
        for (int i=0; i<peakIntensities.size(); i++)
            peakIntensitiesFloat[i] = peakIntensities.get(i);
        return scaleAndCalcKL(idealPeaks, peakIntensitiesFloat, numPeaksToUse);
    }

    protected float scaleAndCalcKL(float[] idealPeaks, float[] peakIntensities)
    {
        return scaleAndCalcKL(idealPeaks, peakIntensities, numPeaksToUse);
    }

    protected float scaleAndCalcKL(float[] idealPeaks, float[] peakIntensities, int numPeaksToUse)
    {
        float[] peakIntensitiesXPeaks = peakIntensities;
        float[] idealPeaksXPeaks = idealPeaks;

        if (numPeaksToUse < idealPeaks.length || (idealPeaks.length != peakIntensities.length))
        {
            peakIntensitiesXPeaks = new float[numPeaksToUse];
            idealPeaksXPeaks = new float[numPeaksToUse];


            float sumRealPeaks = 0;
            float sumIdealPeaks = 0;
            for (int i=0; i<numPeaksToUse; i++)
            {
                if (i < peakIntensities.length)
                    peakIntensitiesXPeaks[i] = peakIntensities[i];
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
    public QuantPeakSetSummary calcPeakIntensities(int firstScan, int lastScan, float mass, float mz,
                                                      float charge, MSRun run,
                                                      int numPeaks)
    {
        float mzTol = (MassUtilities.calculateAbsoluteDeltaMass(
                mass, peakPPMTolerance, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / charge;
        QuantPeakSetSummary result = new QuantPeakSetSummary();
        result.monoisotopicMass = mass;
        result.peakSumIntensities = new ArrayList<Float>();

        float[] mzValues = new float[numPeaks+1];
        for (int peakIndex=-1; peakIndex<numPeaks; peakIndex++)
        {
            mzValues[peakIndex+1] = mz +
                    ((peakIndex * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH) /
                            charge);
        }

        float[] peakIntensities = extractPeakSumIntensities(firstScan, lastScan, mzValues, mzTol, run);
        result.sumIntensityPeakBelow = peakIntensities[0];
        for (int i=1; i<numPeaks+1; i++)
            result.peakSumIntensities.add(peakIntensities[i]);

        return result;
    }

    /**
     * Extract the maximum intensities of the peaks at the defined mz values, within mzTolerance, summed over
     * all scans in the range firstScan..lastScan
     * @param firstScan
     * @param lastScan
     * @param mzValues
     * @param mzTolerance
     * @param run
     * @return
     */
    protected float[] extractPeakSumIntensities(int firstScan, int lastScan, float[] mzValues, float mzTolerance, MSRun run)
    {
        //Define the scan index range, making sure not to go out of bounds
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

        int numPeaks = mzValues.length;
        float[] result = new float[numPeaks];
        for (int scanIndex = 0; scanIndex < scans.length; scanIndex++)
        {
            float[][] spectrum = scans[scanIndex].getSpectrum();

            for (int i=0; i<numPeaks; i++)
            {
                float peakMzCenter = mzValues[i];

                int startIndex = Arrays.binarySearch(spectrum[0], peakMzCenter - mzTolerance);
                startIndex = startIndex < 0 ? -(startIndex+1) : startIndex;
                int endIndex = Arrays.binarySearch(spectrum[0], peakMzCenter + mzTolerance);
                endIndex = endIndex < 0 ? -(endIndex+1) : endIndex;
                //dhmay 20091222, in a hurry -- range check
                if (endIndex >= spectrum[1].length)
                    endIndex = spectrum[1].length-1;

                float maxIntensityThisPeak = 0;
                for (int j=startIndex; j<=endIndex; j++)
                    maxIntensityThisPeak = Math.max(maxIntensityThisPeak, spectrum[1][j]);
                result[i] += maxIntensityThisPeak;
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

//    public int getLabelType()
//    {
//        return labelType;
//    }

//    public void setLabelType(int labelType)
//    {
//        this.labelType = labelType;
//    }

    /**
     * Data structure to pass around summary information about one set of peaks (light or heavy)
     */
    public static class QuantPeakSetSummary
    {
        protected float monoisotopicMass;
        protected float sumIntensityPeakBelow;
        protected List<Float> peakSumIntensities;

        public QuantPeakSetSummary()
        {
            peakSumIntensities = new ArrayList<Float>();
        }

        public String toString()
        {
            StringBuffer resultBuf = new StringBuffer("QuantPeakSetSummary, monoMass = " + monoisotopicMass +
                    ", peakBelow=" + sumIntensityPeakBelow + ", peak intensities: ");
            for (int i=0; i<peakSumIntensities.size(); i++)
                resultBuf.append(" " + (i+1) + "=" + peakSumIntensities.get(i));
            return resultBuf.toString();
        }

        public List<Float> getPeakSumIntensities()
        {
            return peakSumIntensities;
        }
    }

    /**
     * Class for communicating assessment results
     */
    public static class QuantEventAssessment
    {
        protected float singlePeakRatio;

        protected boolean[] flagBitmap = new boolean[flagReasonCodes.length];
        protected String[] checkExplanations = new String[flagReasonCodes.length];

        public String toString()
        {
            //todo: explanations?
            StringBuffer result = new StringBuffer("Assessment. Flags:");
            for (int i=0; i<flagBitmap.length; i++)
                result.append(" " + flagReasonCodes[i] + "=" + flagBitmap[i]);
            result.append(", single-peak: " + singlePeakRatio);
            return result.toString();
        }

        public QuantEventAssessment(int status, String explanation)
        {
            setStatus(status);
            setExplanation(explanation);
        }

        public QuantEventAssessment(boolean[] flagBitmap, String explanation)
        {
            this.flagBitmap = flagBitmap;
            setCheckExplanation(getStatus(), explanation);
        }

        public int getStatus()
        {
            for (int i=0; i<flagBitmap.length; i++)
                if (flagBitmap[i])
                    return i;
            return -1;
        }

        public void setBad()
        {
            setFlag(FLAG_REASON_OK, false);
        }

        public void setFlag(int status, boolean value)
        {
            flagBitmap[status] = value;
        }

        public void setStatus(int status)
        {
            flagBitmap = new boolean[flagReasonCodes.length];
            flagBitmap[status] = true;
        }

        public void setCheckExplanation(int status, String description)
        {
            checkExplanations[status] = description;
        }

        public String getCheckExplanation(int status)
        {
            return checkExplanations[status];
        }

        public String getExplanation()
        {
            return getCheckExplanation(getStatus());
        }

        public void setExplanation(String explanation)
        {
            setCheckExplanation(getStatus(), explanation);
        }

        public float getSinglePeakRatio()
        {
            return singlePeakRatio;
        }

        public void setSinglePeakRatio(float singlePeakRatio)
        {
            this.singlePeakRatio = singlePeakRatio;
        }

        public boolean isGood()
        {
            return flagBitmap[FLAG_REASON_OK];
        }

        public boolean[] getFlagBitmap()
        {
            return flagBitmap;
        }

        public void setFlagBitmap(boolean[] flagBitmap)
        {
            this.flagBitmap = flagBitmap;
        }
    }


    public boolean isShouldFlagIfMissingPeaks()
    {
        return shouldFlagIfMissingPeaks;
    }

    public void setShouldFlagIfMissingPeaks(boolean shouldFlagIfMissingPeaks)
    {
        this.shouldFlagIfMissingPeaks = shouldFlagIfMissingPeaks;
    }

    public boolean isShouldFlagDifferentKL()
    {
        return shouldFlagDifferentKL;
    }

    public void setShouldFlagDifferentKL(boolean shouldFlagDifferentKL)
    {
        this.shouldFlagDifferentKL = shouldFlagDifferentKL;
    }

    public boolean isShouldPerformAllChecks()
    {
        return shouldPerformAllChecks;
    }

    public void setShouldPerformAllChecks(boolean shouldPerformAllChecks)
    {
        this.shouldPerformAllChecks = shouldPerformAllChecks;
    }

    public int getNumPeaksToUse() {
        return numPeaksToUse;
    }

    public void setNumPeaksToUse(int numPeaksToUse) {
        this.numPeaksToUse = numPeaksToUse;
    }
}
