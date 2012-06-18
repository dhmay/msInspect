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

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.apache.log4j.Logger;

import java.util.Arrays;

/**
 * For adjusting the masses of features found in a resampled space, to take
 * advantage of the higher mass accuracy to be found in the unresampled space
 *
 * Optionally recalculate both total and maximum intensity for each feature.  This is done in a completely
 * different way than the original calculation (which occurs in resampled space), and so the two types of
 * values should not be compared.  Calculation is done here for massive performance benefit over recalculating
 * afterward, and because recalculating intensity based on accurate-mass peaks only makes sense if we
 * can determine mass accurately.  Currently this is only done for profile-mode data.
 *
 * If intensity based on accurate mass is specified, and no intensity values are found in the narrow
 * mass range specified, then accurate mass is invalidated, as well.
 */
public class AccurateMassAdjuster
{
    private static Logger _log = Logger.getLogger(AccurateMassAdjuster.class);


    public static final int DEFAULT_SCAN_WINDOW_SIZE = 1;

    protected int _scanWindowSize = DEFAULT_SCAN_WINDOW_SIZE;

    //modes of calculating accurate mass in profile mode.  CENTER uses the mean m/z value from all scans, and
    //each scan's value is the weighted mean of all the values in range.  MAX uses the maximum m/z value from
    //all scans, and each scan's value is the m/z value at maximum intensity in range.
    //CENTER is more appropriate for peptides, where multiple peptides in range are unlikely.
    //MAX is more appropriate for metabolites, where multiple metabolites in range are very likely.
    public static final int PROFILE_MASS_MODE_CENTER = 0;
    public static final int PROFILE_MASS_MODE_MAX = 1;

    public static final int DEFAULT_PROFILE_MASS_MODE = PROFILE_MASS_MODE_CENTER;

    //When calculating the mass of profile-mode data, use a weighted average, or maximum?
    protected int profileMassMode = DEFAULT_PROFILE_MASS_MODE;

    //proportion of the resampling bucket size (e.g., 1/36Da for a resmpling frequency of 36) to search
    //for the highest peak
    public static final double DEFAULT_RESAMPLING_SIZE_PROPORTION = .66666667;
    protected double resamplingSizeProportion = DEFAULT_RESAMPLING_SIZE_PROPORTION;

    protected boolean shouldAdjustComprisedMasses = false;

    //dhmay adding 20100511
    protected boolean shouldRecalculateIntensities = false;
    
    //todo: This should be adjustable.
    //Right now it's only used by FeatureStrategySmallMolecule, so tight tolerance
    //is appropriate.  What tolerance, exactly?  Need to find a good way to figure it out
    protected float intensityRecalcPPMError = 3f;

    protected int maxUsedPeaks = 2;

    //dhmay adding for debugging
    protected float debugMass = 303.29f;
    protected float debugDeltaMass = 0.05f;
    protected boolean massDebugMode = false;

    public void adjustAllMasses(MSRun run, Feature[] features)
            throws InterruptedException
    {
        _log.debug("adjustAllMasses, " + toString());
        Thread currentThread = Thread.currentThread();
        Arrays.sort(features, Spectrum.comparePeakScanAsc);
        for (Feature feature : features)
        {
            float mz = calculateAccurateMz(run, feature);
            if (mz > 0)
            {
                feature.setMz(mz);
                feature.setAccurateMZ(true);
                feature.updateMass();
            }
            if (currentThread.isInterrupted())
                throw new InterruptedException();
        }
    }

    /**
     * Switch on whether the run is in centroid mode and do the appropriate thing
     * @param run
     * @param f
     * @return
     */
    public float calculateAccurateMz(MSRun run, Feature f)
    {

        if (run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
            return calculateAccurateMzCentroid(run, f);
        else
            return calculateAccurateMzProfile(run, f);
    }


    public float calculateAccurateMzCentroid(MSRun run, Feature f)
    {
        // NOTE about isotopes
        // some heavy isotopes are greater than +1Da (C=1.0033), and some are less (N=0.9971)
        // I use 1.0013 as working number (C13 is more common than N15)
        final double ISOTOPE_FACTOR = 1.0013;

        // CONSIDER: does it make sense to do any sort of averaging across a few scans?
        Scan scan = run.getScan(run.getIndexForScanNum(f.scan));
        float[][] s = scan.getSpectrum();
        double delta = .66 / SpectrumResampler.getResampleFrequency();

        int numPeaksToUse = Math.min(maxUsedPeaks, f.comprised.length);
        double mzP0 = 0;
        double[] peakSumMZs = new double[numPeaksToUse];
        double[] peakSumINs = new double[numPeaksToUse];

        for (int i = 0; i < 2 && i < f.comprised.length; i++)
        {
            double mzPi;
            if (f.comprised[i] == null)
                mzPi = f.mz + i * ISOTOPE_FACTOR / f.charge;
            else
                mzPi = f.comprised[i].mz;

            // find biggest close peak
            // could use more complicated "distance" measure
            // but in centroided data there is usually only one 'real' peak
            // in a space this small with maybe some tiny side peaks
            int p = Arrays.binarySearch(s[0], (float)(mzPi - delta));
            if (p < 0) p = -1 * (p + 1);
            double mzBiggest = 0;
            double inBiggest = 0;
            for (; p < s[0].length; p++)
            {
                float mz = s[0][p];
                if (mz > mzPi + delta)
                    break;
                float in = s[1][p];
                if (in > inBiggest)
                {
                    mzBiggest = mz;
                    inBiggest = in;
                }
            }
            if (mzBiggest == 0) // UNDONE: use wider window???
            {
                //System.out.println("missed " + i + "\t" + f.toString());
                return 0;
            }
            if (f.charge == 0)
            {
                return (float)mzBiggest;
            }
            if (i == 0)
                mzP0 = mzBiggest;

            double peakMz = inBiggest * (mzBiggest - i * ISOTOPE_FACTOR / f.charge);

            // weighted average of peaks
            peakSumMZs[i] = peakMz;
            peakSumINs[i] += inBiggest;
        }

        double avgMZ = BasicStatistics.sum(peakSumMZs) / BasicStatistics.sum(peakSumINs);
        for (int i=0; i<maxUsedPeaks; i++)
        {
            if (f.comprised[i] != null)
                f.comprised[i].setMz((float)(peakSumMZs[i] / peakSumINs[i]));

        }

        // if we got this all right we'd expect ABS(newMZ-mzP0) to be less
        // than the resolution of the machine
        // I'm going to assume that the centroided means FT (very accurate)
        if (Math.abs(avgMZ - mzP0) > mzP0 * 5.0 / 1000000.0)
        {
            return 0;
        }


        // NOTE: could return avgMZ here, however
        //  a) it's very close to mzP0 (see test)
        //  b) I suspect people would rather see the value of peak from the source data
        return (float)mzP0;
    }

    /**
     *
     * adjustment of profile-mode mass; average results from some number of adjacent scans.
     * Optionally adjust masses of comprised peaks, too.   Note: if comprised peak mass adjustment fails,
     * there's no report made of the failure.
     *
     * @param run
     * @param f
     * @return
     */
    public float calculateAccurateMzProfile(MSRun run, Feature f)
    {
        if (Math.abs(debugMass - f.mass) < debugDeltaMass)
        {
            massDebugMode = true;
            _log.debug("MASS: " + f.mass);
        }
        if (_scanWindowSize <= 0)
            return 0.f;

        int scanIndex = run.getIndexForScanNum(f.scan);

        int lowScanIndex  = (int) Math.max(scanIndex - (_scanWindowSize - 1)/2.0 + .5, 0);
        int highScanIndex = (int) Math.min(scanIndex + (_scanWindowSize - 1)/2.0 + .5, run.getScanCount() - 1);

        Pair<Float, Pair<Float, Float>> accMzAndIntensities = calculateAccurateMzProfileForMultiScanPeak(
                run, f.mz, f.charge, lowScanIndex, highScanIndex);
        float monoisotopicAdjustedMz = accMzAndIntensities.first;
        if (shouldRecalculateIntensities && monoisotopicAdjustedMz > 0)
        {
            f.intensity = accMzAndIntensities.second.first;
            f.totalIntensity = accMzAndIntensities.second.second;
        }

        if (shouldAdjustComprisedMasses)
        {
            for (int i=0; i<f.comprised.length; i++)
            {
                if (f.comprised[i] == null)
                    continue;

                if (i==0)
                {
                    if (monoisotopicAdjustedMz > 0)
                        f.comprised[i].setMz(monoisotopicAdjustedMz);
                }
                else
                {
                    Pair<Float, Pair<Float, Float>> accMzAndIntensitiesComprisedPeak =
                            calculateAccurateMzProfileForMultiScanPeak(run, f.comprised[i].mz,
                            f.charge, lowScanIndex, highScanIndex);
                    float accPeakMz = accMzAndIntensitiesComprisedPeak.first;

                    if (accPeakMz > 0)
                    {
                        f.comprised[i].setMz(accPeakMz);
                        if (shouldRecalculateIntensities)
                        {
                             f.comprised[i].intensity = accMzAndIntensitiesComprisedPeak.second.first;
                             //nothing to do with total intensity for comprised peak
                        }
                    }
                }
            }
        }

        massDebugMode = false;

        return monoisotopicAdjustedMz;
    }

    /**
     * Calculate accurate MZ for a multi-scan peak.
     *
     * Optionally adjust intensity of each peak based on accurate mass; if we do that, and can't find
     * intensity values near accurate mass, invalidate accurate mass.
     * @param run
     * @param mz
     * @param charge
     * @param lowScanIndex
     * @param highScanIndex
     * @return
     */
    protected Pair<Float, Pair<Float, Float>> calculateAccurateMzProfileForMultiScanPeak(
            MSRun run, float mz, int charge, int lowScanIndex, int highScanIndex)
    {
        //for PROFILE_MASS_MODE_CENTER
        float sumMz = 0.f;
        int n = 0;

//        float[] scanMzValues = new float[highScanIndex - lowScanIndex + 1];
//List<Float> nonzeroScanMzValues = new ArrayList<Float>();

        //for PROFILE_MASS_MODE_MAX
        float mzAtMaxInt = 0;
        float maxInt = 0;
        for (int s = lowScanIndex; s <= highScanIndex; s++)
        {
            Scan scan = run.getScan(s);
            if (massDebugMode)
                _log.debug("Scan " + s);
            Pair<Float, Float> scanMzAndInt = calculateAccurateMzProfileForScan(scan, mz);
            float scanMz = scanMzAndInt.first;
            float scanInt = scanMzAndInt.second;

            if (scanMz > 0.f)
            {
                sumMz += scanMz;
                n++;
//                scanMzValues[s-lowScanIndex] = scanMz;
//nonzeroScanMzValues.add(scanMz);
                if (scanInt > maxInt)
                {
                    maxInt = scanInt;
                    mzAtMaxInt = scanMz;
                }
            }
        }
//float massSpreadPPM = MassUtilities.convertDaToPPM((float) (BasicStatistics.max(nonzeroScanMzValues) - BasicStatistics.min(nonzeroScanMzValues)), mz);
//if (massSpreadPPM > 5) System.err.println("BIG MASS SPREAD: " + massSpreadPPM);

        float accMz = 0.f;
        if (n > 0) //n always > 0 if profileMassMode == PROFILE_MASS_MODE_MAX
        {
            switch (profileMassMode)
            {
                case PROFILE_MASS_MODE_CENTER:
                    accMz = sumMz / n;
                    break;
                case PROFILE_MASS_MODE_MAX:
                    accMz = mzAtMaxInt;
                    break;
            }
        }

        Pair<Float, Pair<Float, Float>> result = new Pair<Float, Pair<Float, Float>>(accMz, null);

        //if we found an accurate mz, then recalculate intensity if we're supposed to.
        //This could result in accMz being invalidated
        if (shouldRecalculateIntensities && (accMz > 0))
        {
            Scan[] scans = run.getScans();

            //if charge 0, assume charge 1
            int absChargeMin1 = Math.max(Math.abs(charge), 1);
            float deltaMz = MassUtilities.convertPPMToDa(intensityRecalcPPMError,
                    Feature.convertMzToMass(mz, absChargeMin1));
            float lowAccMz = accMz - deltaMz;
            float highAccMz = accMz + deltaMz;

            float totalIntensitySum = 0;
            float maxScanIntensity = 0;
            for (int i=lowScanIndex; i<=highScanIndex; i++)
            {
                Scan scan = run.getScan(i);
                float[][] spectrum = scan.getSpectrum();

                int p = Arrays.binarySearch(spectrum[0], lowAccMz);
                if (p < 0)
                    p = -1 * (p + 1);
                float intensityThisScan = 0;
                for (; p<spectrum[0].length; p++)
                {
                    if (spectrum[0][p] > highAccMz)
                        break;
                    intensityThisScan += spectrum[1][p];
                }
                maxScanIntensity = Math.max(intensityThisScan, maxScanIntensity);

                //"integrate" scan intensity over time; multiply this scan's intensity by half the difference
                //in time between the scan before and after
                if (intensityThisScan > 0)
                {
                    double time = 1.0;
                    if (i > 0 && i + 1 < run.getScanCount())
                        time = (scans[i + 1].getDoubleRetentionTime() -
                                scans[i - 1].getDoubleRetentionTime()) /
                                2;
                    totalIntensitySum += time * intensityThisScan;
                }
            }

            //if no intensities in range, something went terribly wrong, i.e., the "accurate" mass value
            // wasn't so accurate after all.  Don't correct accurate mass or assign intensity
            if (maxScanIntensity == 0)
            {
                result.first = 0f;
//System.err.println("No acc intensity!  accMz = " + accMz + ", low=" + lowAccMz + ", high=" + highAccMz + ", orig=" + mz);
//for (float val : scanMzValues)
//    System.err.println("\t" + val);
//
//            for (int i=lowScanIndex; i<=highScanIndex; i++)
//            {
//                Scan scan = run.getScan(i);
////System.err.println("\t\tscan " + i);
//                float[][] spectrum = scan.getSpectrum();
//
//                int p = Arrays.binarySearch(spectrum[0], lowAccMz);
//                if (p < 0)
//                    p = -1 * (p + 1);
//                float intensityThisScan = 0;
//                for (; p<spectrum[0].length; p++)
//                {
//                    if (spectrum[0][p] > highAccMz)
//                        break;
////System.err.println("\t\t\tmz=" + spectrum[0][p] + ", int=" + spectrum[1][p]);
//
//                    intensityThisScan += spectrum[1][p];
//                }
//            }
//

            }
            else
                result.second = new Pair<Float, Float>(maxScanIntensity, totalIntensitySum);
        }

        return result;
    }

    /**
     * adjustment of profile-mode mass using either maximum or center-of-mass, depending on profileMassMode
     * @param scan
     * @param mz
     * @return pair containing acc mass and intensity.  intensity is either sum or max, depending on profileMassMode
     */
    protected Pair<Float, Float> calculateAccurateMzProfileForScan(Scan scan, float mz)
    {
        float[][] s = scan.getSpectrum();
        double delta = resamplingSizeProportion / SpectrumResampler.getResampleFrequency();

        double lowMz = mz - delta;
        double highMz = mz + delta;

        int p = Arrays.binarySearch(s[0], (float) lowMz);
        if (p < 0)
            p = -1 * (p + 1);

        double sumMz = 0;
        double sumInt = 0;

        double mzAtMaxInt = 0;
        double maxInt = 0;

        for (; p < s[0].length; p++)
        {
            if (s[0][p] > highMz)
            {
                if (massDebugMode) _log.debug("outofrange: " + s[0][p]);
                break;
            }

            switch (profileMassMode)
            {
                case PROFILE_MASS_MODE_CENTER:
                    sumMz += s[0][p] * s[1][p]; // mz weighted by intensity
                    sumInt += s[1][p];
                    break;
                case PROFILE_MASS_MODE_MAX:
                    if (s[1][p] > maxInt)
                    {
                        mzAtMaxInt = s[0][p];
                        maxInt = s[1][p];
                    }
                    break;
            }
        }

        float result = 0;
        switch (profileMassMode)
        {
            case PROFILE_MASS_MODE_CENTER:
                if (massDebugMode)
                    _log.debug("\tmzWeightedMean=" + sumMz/sumInt + ", int=" + sumInt);
                // Couldn't figure out a decent match
                if (sumInt <= 0.0)
                    result = 0.f;
                else
                    result = (float) (sumMz/sumInt);
                return new Pair<Float, Float>(result, (float) sumInt);
            case PROFILE_MASS_MODE_MAX:
                if (massDebugMode)
                    _log.debug("\tmzMax=" + mzAtMaxInt + ", int=" + maxInt);
                result = (float) mzAtMaxInt;
                return new Pair<Float, Float>(result, (float) maxInt);
        }
        //never gets here
        return null;
    }

    public String toString()
    {
        return "AccurateMassAdjuster: adjustAllMasses, profileMassMode=" + profileMassMode + ", scanWindow: " + _scanWindowSize +
                ", resamplingProporation: " + resamplingSizeProportion;
    }

    public int getScanWindowSize()
    {
        return _scanWindowSize;
    }

    public void setScanWindowSize(int scanWindowSize)
    {
        _scanWindowSize = scanWindowSize;
    }

    public double getResamplingSizeProportion()
    {
        return resamplingSizeProportion;
    }

    public void setResamplingSizeProportion(double resamplingSizeProportion)
    {
        this.resamplingSizeProportion = resamplingSizeProportion;
    }

    public int getProfileMassMode()
    {
        return profileMassMode;
    }

    public void setProfileMassMode(int profileMassMode)
    {
        this.profileMassMode = profileMassMode;
    }


    public static Logger get_log()
    {
        return _log;
    }

    public static void set_log(Logger _log)
    {
        AccurateMassAdjuster._log = _log;
    }

    public boolean isShouldAdjustComprisedMasses() {
        return shouldAdjustComprisedMasses;
    }

    public void setShouldAdjustComprisedMasses(boolean shouldAdjustComprisedMasses) {
        this.shouldAdjustComprisedMasses = shouldAdjustComprisedMasses;
    }

    public boolean isShouldRecalculateIntensities() {
        return shouldRecalculateIntensities;
    }

    public void setShouldRecalculateIntensities(boolean shouldRecalculateIntensities) {
        this.shouldRecalculateIntensities = shouldRecalculateIntensities;
    }

    public float getIntensityRecalcPPMError() {
        return intensityRecalcPPMError;
    }

    public void setIntensityRecalcPPMError(float intensityRecalcPPMError) {
        this.intensityRecalcPPMError = intensityRecalcPPMError;
    }
}
