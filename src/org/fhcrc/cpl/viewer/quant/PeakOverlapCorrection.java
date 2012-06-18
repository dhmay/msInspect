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

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;

import java.util.List;

/**
 *  Methods for correcting isotopically labeled light and heavy areas to account for peak overlap.
 *
 * NOTE: This class is not currently used by Q3, which has its own implementation. This implementation should be
 * numerically equivalent to the correction in Q3, which is done using matrices. 
 *
 * TODO: make Q3 use this class for peak overlap correction
 */
public class PeakOverlapCorrection
{
    private static Logger _log = Logger.getLogger(PeakOverlapCorrection.class);

    //dhmay adding 20100312
    public static final float INFINITE_RATIO_VALUE = 999f;

    //store the first 25 elements of the factorial sequence, so we don't have to calculate it every time
    protected static final double[] FACTORIAL_SEQUENCE_25 = new double[] {
            1.0, 1.0, 2.0, 6.0, 24.0,
            120.0, 720.0, 5040.0, 40320.0, 362880.0,
            3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
            1.307674e+12, 2.092279e+13, 3.556874e+14, 6.402374e+15, 1.216451e+17,
            2.432902e+18, 5.109094e+19, 1.124001e+21, 2.585202e+22, 6.204484e+23};

    /**
     * Get isotope intensities based on a simple poisson model. Sum of all intensities approaches 1 as len -> infinity
     *              p(x) = lambda^x exp(-lambda)/x!
     *
     * Equivalent R code:
     *
     * iso <- function(mass, len)
     * {
     *   lambda=mass/1800;
     *   result=c();
     *   for (i in 0:len-1)
     *   {
     *     result = c(result, lambda^i * exp(-lambda)/factorial(i));
     *   }
     *   result;
     * }
     * @param mass
     * @param len
     * @return
     */
    public static double[] getIsotopicDistribution(double mass, int len)
    {
        double lambda = mass / 1800;
        double lambdaExp = Math.exp(-lambda);

        double[] result = new double[len];

        for (int i = 0; i < len; i++)
            if (i < FACTORIAL_SEQUENCE_25.length)
                result[i] = Math.pow(lambda, i) / FACTORIAL_SEQUENCE_25[i] * lambdaExp;
            else
                result[i] = 0; // Not really zero, but too small to matter

        return result;
    }

    /**
     * This is trivial, but it's nice to have it to use elsewhere.  Technically it needs a slight correction
     * for mass defect, but in practice, not really unless the peak difference is gargantuan
     * @param lightMass
     * @param heavyMass
     * @return
     */
    public static int calcNumPeaksSeparation(double lightMass, double heavyMass)
    {
        return (int) Math.round((heavyMass - lightMass));
    }

    /**
     * Account for any incursion of light onto heavy area, and return adjusted areas.  Based on raw peak area
     * calculations that begin with the firstpeak (for monoisotope, use 0) and continue through the last 
     * free-and-clear light peak that's <= highestPeakMax
     * @param lightMass
     * @param heavyMass
     * @param rawLightArea
     * @param rawHeavyArea
     * @param firstPeak the peak to start with (0-based).  Normally (e.g., for Q3) this will be the monoisotope
     * @param highestPeakMax MAXIMUM highest isotopic peak used in calculation of areas (0-based).  If greater than
     * the number of Daltons of separation between light and heavy - 1, we assume only safe peaks (no overlap) used.
     * @return a pair of (light, heavy) adjusted areas.  Adjusted areas represent the inferred sum of
     * intensity across all isotopic peaks (even those that overlap), so may be very different from raw
     * areas, especially for high masses
     */
    public static double[] correctLightHeavyAreas(double lightMass, double heavyMass,
                                                  double rawLightArea, double rawHeavyArea,
                                                  int firstPeak, int highestPeakMax)
    {
        //number of Daltons separating light and heavy.  Same as number of safe peaks for light
        int numDaltonsSeparation = calcNumPeaksSeparation(lightMass, heavyMass);

        //assumption about the number of highest light and heavy peak actually used -- maxPeaks or the number of
        //light peaks that don't intrude on heavy, whichever's less
        int maxPeakUsed = Math.min(highestPeakMax, numDaltonsSeparation-1);

        //Get the full isotopic distribution of the light peaks
        //For some reason, Q3 uses the heavy mass, rather than the light, for this distribution calculation.
        //It really doesn't matter that much, so preserving that behavior
        double[] fullLightIsotopicDistribution =
                getIsotopicDistribution(heavyMass, numDaltonsSeparation + maxPeakUsed + 1);

        //head gets sum of first numPeaksUsed peak distribution values.
        //tail gets the sum of the first numPeaksUsed values starting with the heavy monoisotope
        double head = 0;
        for (int i = firstPeak; i <= maxPeakUsed; i++)
                head += fullLightIsotopicDistribution[i];
        double tail = 0;
        int lastUsedHeavyPeakIndex = numDaltonsSeparation + maxPeakUsed;
//System.err.println(numDaltonsSeparation + ", " + maxPeakUsed + ", " + lastUsedHeavyPeakIndex + ", " + firstPeak);        
        for (int i = numDaltonsSeparation + firstPeak; i <= lastUsedHeavyPeakIndex; i++)
        {
                tail += fullLightIsotopicDistribution[i];
//System.err.println("  " + fullLightIsotopicDistribution[i]);
        }
//if (tail > 0) System.err.println(tail);        

        //Corrected light area represents inferred sum across all peaks
        double correctedLightArea = rawLightArea / head;
        //Corrected heavy area first subtracts the inferred contribution of light peaks, then expands to represent
        //inferred sum over all peaks
        //dhmay changing 20100312.  The corrected heavy area can be negative.  Fix it to zero if so.  Note: be sure
        //to check for this when calculating ratio
        double correctedHeavyArea = Math.max(0, (rawHeavyArea - (tail * correctedLightArea)) / head);

        _log.debug("Raw: " + rawLightArea + ", " + rawHeavyArea + ". Sep: " + numDaltonsSeparation + ". head=" + head +
                ", tail=" + tail + 
                ". Corrected: " + correctedLightArea + ", " + correctedHeavyArea);

        return new double[] { correctedLightArea, correctedHeavyArea };
    }

    /**
     * Cover method for starting with monoisotope
     * @param lightMass
     * @param heavyMass
     * @param ratio
     * @param maxPeaks
     * @return
     */
    public static double correctRatioForOverlap(double lightMass, double heavyMass, float ratio, int maxPeaks)
    {
        return correctRatioForOverlap(lightMass, heavyMass, ratio, 0, maxPeaks);
    }

    /**
     * Correct the ratio by calling correctLightHeavyAreas, dividing light by heavy.  Absolute intensities don't
     * matter in this calculation, just the ratio
     * @param lightMass
     * @param heavyMass
     * @param ratio
     * @param firstPeak
     * @param maxPeaks
     * @return
     */
    public static double correctRatioForOverlap(double lightMass, double heavyMass, float ratio, int firstPeak,
                                                int maxPeaks)
    {
        double[] lightHeavy = correctLightHeavyAreas(lightMass, heavyMass, 1, 1f / ratio, firstPeak, maxPeaks);
        //dhmay adding 20090312
        if (lightHeavy[1] == 0)
            return INFINITE_RATIO_VALUE;
        return lightHeavy[0] / lightHeavy[1];
    }


    /**
     * Calculate a KL score, using the passed-in 6-peak intensity template as the ideal.
     *
     * Cribbed from Spectrum.KLPoissonDistance and made more flexible
     * @param idealPeaks must sum to 1
     * @param peakIntensities must sum to 1
     * @return
     */
    public static float calcKLUsingTemplate(float[] idealPeaks, float[] peakIntensities)
    {
        assert idealPeaks.length == peakIntensities.length;

        double diff = 0.0;
        double sumP = 0.0;
        double sumQ = 0.0;
        int pqMinLength = Math.min(idealPeaks.length, peakIntensities.length);
//System.err.println("KL, peaks: " + pqMinLength);

        for (int k = 0; k < pqMinLength ; k++)
        {
            diff += idealPeaks[k] * Math.log((double)idealPeaks[k] / peakIntensities[k]);
            sumP += idealPeaks[k];
            sumQ += peakIntensities[k];
//System.err.println(" Peak " + k + ", p=" + idealPeaks[k] + ", q=" + peakIntensities[k] + ", diff: " +
//  (idealPeaks[k] * Math.log((double)idealPeaks[k] / peakIntensities[k])));
        }
        double kl = diff / Spectrum.LN2;
        kl = Math.max(kl, 0.0);
//System.err.println("KL: " + kl);        
        assert Math.abs(sumP-1.0) < 0.001 && Math.abs(sumQ-1.0) < 0.001;

        return (float) kl;
    }
}
