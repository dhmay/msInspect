package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.apache.log4j.Logger;

import java.util.Arrays;

/**
 * For adjusting the masses of features found in a resampled space, to take
 * advantage of the higher mass accuracy to be found in the unresampled space
 */
public class AccurateMassAdjuster
{
    private static Logger _log = Logger.getLogger(AccurateMassAdjuster.class);


    public static final int DEFAULT_RESAMPLING_FREQUENCY = 36;
    public static final int DEFAULT_SCAN_WINDOW_SIZE = 1;

    protected int _resamplingFrequency = DEFAULT_RESAMPLING_FREQUENCY;

    protected int _scanWindowSize = DEFAULT_SCAN_WINDOW_SIZE;

    public static final int PROFILE_MASS_MODE_CENTER = 0;
    public static final int PROFILE_MASS_MODE_MAX = 1;

    public static final int DEFAULT_PROFILE_MASS_MODE = PROFILE_MASS_MODE_CENTER;

    //When calculating the mass of profile-mode data, use a weighted average, or maximum?
    protected int profileMassMode = DEFAULT_PROFILE_MASS_MODE;

    //proportion of the resampling bucket size (e.g., 1/36Da for a resmpling frequency of 36) to search
    //for the highest peak
    public static final double DEFAULT_RESAMPLING_SIZE_PROPORTION = .66666667;
    protected double resamplingSizeProportion = DEFAULT_RESAMPLING_SIZE_PROPORTION;


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
            float mz = calculateAccurateMass(run, feature);
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
    public float calculateAccurateMass(MSRun run, Feature f)
    {

        if (run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
            return calculateAccurateMassCentroid(run, f);
        else
            return calculateAccurateMassProfile(run, f);
    }


    public float calculateAccurateMassCentroid(MSRun run, Feature f)
    {
        // NOTE about isotopes
        // some heavy isotopes are greater than +1Da (C=1.0033), and some are less (N=0.9971)
        // I use 1.0013 as working number (C13 is more common than N15)
        final double ISOTOPE_FACTOR = 1.0013;

        // CONSIDER: does it make sense to do any sort of averaging across a few scans?
        Scan scan = run.getScan(run.getIndexForScanNum(f.scan));
        float[][] s = scan.getSpectrum();
        double delta = .66 / _resamplingFrequency;



        double mzP0 = 0;
        double sumMZ = 0;
        double sumIN = 0;
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
            sumMZ += peakMz;
            sumIN += inBiggest;
        }
        double avgMZ = sumMZ / sumIN;
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
     * adjustment of profile-mode mass; average results from some number of adjacent scans
     */
    public float calculateAccurateMassProfile(MSRun run, Feature f)
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

        float sumMz = 0.f;
        int n = 0;

        for (int s = lowScanIndex; s <= highScanIndex; s++)
        {
            Scan scan = run.getScan(s);
            if (massDebugMode)
                _log.debug("Scan " + s);
            float maxMz = calculateAccurateMassProfileForScan(scan, f);


            if (maxMz > 0.f)
            {
                sumMz += maxMz;
                n++;
            }
        }
        massDebugMode = false;
        return n > 0 ? sumMz / n : 0.f;
    }

    /**
     * adjustment of profile-mode mass using either maximum or center-of-mass, depending on profileMassMode
     * @param scan
     * @param f
     * @return
     */
    protected float calculateAccurateMassProfileForScan(Scan scan, Feature f)
    {
        float[][] s = scan.getSpectrum();
        double delta = resamplingSizeProportion / _resamplingFrequency;

        double lowMz = f.mz - delta;
        double highMz = f.mz + delta;

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
                break;
            case PROFILE_MASS_MODE_MAX:
                if (massDebugMode)
                    _log.debug("\tmzMax=" + mzAtMaxInt + ", int=" + maxInt);
                result = (float) mzAtMaxInt;
                break;
        }

        return result;
    }

    public String toString()
    {
        return "AccurateMassAdjuster: adjustAllMasses, profileMassMode=" + profileMassMode + ", scanWindow: " + _scanWindowSize +
                ", resamplingProporation: " + resamplingSizeProportion;
    }


    public int getResamplingFrequency()
    {
        return _resamplingFrequency;
    }

    public void setResamplingFrequency(int resamplingFrequency)
    {
        _resamplingFrequency = resamplingFrequency;
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
}
