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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.apache.commons.lang.math.IntRange;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;

import java.util.*;

public class Feature extends Spectrum.Peak implements Cloneable
{
    static Logger _log = Logger.getLogger(Feature.class);

    public float time = 0;
    public float mass = 0;
    public int charge = 0;
    protected boolean accurateMZ = false;  // do I consider this highly accurate (e.g. from centroided input data)

    // Measures filled in by ScoreFeature()
    public float kl = -1.0F;
    public float dist = -1.0F;
    public int peaks = 1;
    public boolean skippedPeaks = false;

    public int chargeStates = 1;
    protected String description;
    public int scanFirst;
    public int scanLast;
    public float totalIntensity = 0F;
    public int scanCount = -1;
    public static float MZ_RANGE = .04f;

    // internal, used by ScoreFeature and ExtractPeptideFeatures
    public float mzPeak0;
    public Spectrum.Peak[] comprised;
    public Feature next; // alternate features (for algorihtm development)

    // one-dimensional window around each feature in m/z space
//    public String mzWindow = null;
    public int intensityLeadingPeaks;
    public int intensityTrailingPeaks;
    public float[] intensityWindow;

    //Map of name->object properties for use in situations where you need to store
    //information about a feature that doesn't _always_ need to be stored about a
    //feature.  Initialized to null
    protected Map<String, Object> _propertyMap = null;


    //Aminoacid modifications, usually loaded from pepxml.
    //This remains null unless there are actually modifications.  If there are,
    //then all elements are null except the actual modifications.  Index+1 = position of mod
    protected ModifiedAminoAcid[] mModifiedAminoAcids = null;

    public Feature()
    {
    }


    public Feature(Spectrum.Peak p)
    {
        this.scan = p.scan;
        this.scanFirst = p.scan;
        this.scanLast = p.scan;
        this.scanCount = 1;
        this.mz = p.mz;
        this.intensity = p.intensity;
        this.totalIntensity = p.intensity;
        this.kl = 0; //Should this be the default????
        this.median = p.median;
        this.background = p.background;
    }


    public Feature(Feature feature)
    {
        this.setTime(feature.getTime());
        this.scan = feature.scan;
        this.mz = feature.mz;
        this.intensity = feature.intensity;
        this.charge = feature.charge;
        this.kl = feature.kl;
        this.dist = feature.dist;
        this.background = feature.background;
        this.median = feature.median;
        this.peaks = feature.peaks;
        this.skippedPeaks = feature.skippedPeaks;
        this.description = feature.description;
        this.chargeStates = feature.chargeStates;
        this.totalIntensity = feature.totalIntensity;
        this.scanFirst = feature.scanFirst;
        this.scanLast = feature.scanLast;
        this.scanCount = feature.scanCount;
        this.accurateMZ = feature.accurateMZ;
        this.mass = feature.mass;

        //dhmay adding 20100412
        //Medium-depth copy: Keeps the existing comprised peak objects, but creates a new array to hold them
        if (feature.comprised != null)
        {
            comprised = new Spectrum.Peak[feature.comprised.length];
            System.arraycopy(feature.comprised, 0, comprised, 0, feature.comprised.length);
        }

        //pull all extra data from the other feature
        for (String propertyName : feature.getPropertyMap().keySet())
        {
            setProperty(propertyName,
                        feature.getPropertyMap().get(propertyName));
        }

    }


    public Feature(int scan, float mz, float intensity)
    {
        this.scan = scan;
        this.mz = mz;
        this.intensity = intensity;
        this.totalIntensity = intensity;
    }


    public Feature(int scan, int scanFirst, int scanLast, float mz, float intensity, int charge, float kl, float totalIntensity)
    {
        this.scan = scan;
        this.mz = mz;
        this.intensity = intensity;
        this.charge = charge;
        this.kl = kl;
        this.scanFirst = scanFirst;
        this.scanLast = scanLast;
        this.totalIntensity = totalIntensity;
        //NOTE: This isn't really scan count if scanFirst/Last are not sequential
        this.setScanCount(scanLast - scanFirst + 1);
        updateMass();
    }

    /**
     * Query the properties set on this Feature to determine what extra
     * info types it has.  This isn't just a getter, it takes a small
     * amount of time
     * @return
     */
    public FeatureExtraInformationDef[] determineExtraInformationTypes()
    {
        Set<FeatureExtraInformationDef> extraInfoSet =
                new HashSet<FeatureExtraInformationDef>();
        for (String propertyName : getPropertyMap().keySet())
        {
            FeatureExtraInformationDef extraInfoDef =
                FeatureExtraInformationDef.getInfoTypeForColumn(propertyName);
            if (extraInfoDef != null)
                extraInfoSet.add(extraInfoDef); 
        }
        return extraInfoSet.toArray(new FeatureExtraInformationDef[0]);
    }

    public void afterPopulate()
    {
        if (charge > 0)
        {
            if (mass == 0 && mz > 0)
                updateMass();
            else if (mz == 0 && mass > 0)
                updateMz();
        }
    }


    /**
     * Make the mass agree with the m/z.
     *
     * Note: this method USED to subtract out the mass of the light label, for isotopically
     * labeled peptides.  That no longer is the case.
     *
     * dhmay updating 7/29/2008 to handle negatively charged ions
     */
    public void updateMass()
    {
        mass = convertMzToMass(mz, charge);
    }

    public static float convertMzToMass(float mz, int charge)
    {
        if (charge > 0)
        {
            float m = (mz - Spectrum.HYDROGEN_ION_MASS) * charge;
            m = Math.max(0.0F, m);
            return m;
        }
        else if (charge == 0)
        {
            return 0;
        }
        else //charge < 0, no proton mass to account for
        {
            //dhmay changing the way negative ion mass is calculated.  Assuming M-H ions
            float m = (mz + Spectrum.HYDROGEN_ION_MASS) * -charge;
//            float m = mz * -charge;
            m = Math.max(0.0F, m);
            return m;
        }
    }


    /**
     * Make the m/z agree with the mass.
     *
     * Note: this method USED to subtract out the mass of the light label, for isotopically
     * labeled peptides.  That no longer is the case.
     *
     * dhmay updating 7/29/2008 to handle negatively charged ions      
     */
    public void updateMz()
    {
        mz = convertMassToMz(mass, charge);
    }

    public static float convertMassToMz(float mass, int charge)
    {
        if (charge > 0)
        {
            float m = (mass / charge) + Spectrum.HYDROGEN_ION_MASS;
            m = Math.max(0.0F, m);
            return m;
        }
        else if (charge == 0)
        {
            return 0;
        }
        else //charge < 0, no proton mass to account for
        {

            float m = (mass / -charge);
            m = Math.max(0.0F, m);
            return m;
        }
    }


    public boolean isFeatureInRange(Feature feature)
    {
        return isFeatureInRange(feature, 3, MZ_RANGE);
    }


    public boolean isFeatureInRange(Feature feature, int maxScanGap, float mzGap)
    {
        return feature.charge == charge && feature.getScanLast() >= getScanFirst() - maxScanGap && feature.getScanFirst() <= getScanLast() + maxScanGap && Math.abs(this.mz - feature.mz) <= mzGap;
    }


    public void addFeatureToRange(Feature feature)
    {
        this.setTotalIntensity(this.getTotalIntensity() + feature.getTotalIntensity());
        this.setScanCount(this.getScanCount() + feature.getScanCount());
        if (feature.intensity > this.intensity)
        {
            this.kl = feature.kl;
            this.intensity = feature.intensity;
            this.mz = feature.mz;
            this.mass = feature.mass;
            this.scan = feature.scan;
            this.time = feature.time;
            this.dist = feature.dist;
        }
        if (feature.getScanFirst() < getScanFirst())
            setScanFirst(feature.getScanFirst());
        if (feature.getScanLast() > getScanLast())
            setScanLast(feature.getScanLast());
    }


    private final static String BASE_FEATURE_HEADER =
            "scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist";
    private final static String DESCRIPTION_FEATURE_HEADER =
            "\tdescription";

    public static String getFeatureHeader(FeatureExtraInformationDef[] extraInfoToDisplay)
    {
        StringBuffer resultBuf = new StringBuffer(BASE_FEATURE_HEADER);

        if (extraInfoToDisplay != null && extraInfoToDisplay.length > 0)
        {
            for (FeatureExtraInformationDef extraInfo : extraInfoToDisplay)
            {
                for (String columnName : extraInfo.getColumnNames())
                    resultBuf.append("\t" + columnName);
            }
        }
        //description comes last.  This is for historical reasons... actually it probably
        //makes more sense to put the extra fields last, now, but that would be a change
        //to the format, and there may be scripts out there that depend on the ordering
        resultBuf.append(DESCRIPTION_FEATURE_HEADER);
        return resultBuf.toString();
    }

    /**
     * no-arg toString() doesn't display any extra data
     * @return
     */
    public String toString()
    {
        return toString(null);
    }

    /**
     * In addition to basic fields, displays the extra info it's told to
     * @param extraInfoToDisplay
     * @return
     */
    public String toString(FeatureExtraInformationDef[] extraInfoToDisplay)
    {
        String txt = null;
        if (null != getDescription())
        {
            txt = getDescription().replaceAll("[\n\r]+", "\\\\");
            txt = getDescription().replace('"', '\'');
        }
        StringBuffer out = new StringBuffer();
        out.append(scan + "\t" + getTime() + "\t" + mz + "\t" + (isAccurateMZ() ? "true" : "false") + "\t" + mass);
        out.append("\t" + intensity + "\t" + charge + "\t" + chargeStates);
        out.append("\t" + kl + "\t" + background + "\t" + getMedian() + "\t" + getPeaks());
        out.append("\t" + getScanFirst() + "\t" + getScanLast() + "\t" + getScanCount() + "\t" + getTotalIntensity());
        out.append("\t" + getSumSquaresDist());

        if (extraInfoToDisplay != null)
        {
            for (FeatureExtraInformationDef extraInfo : extraInfoToDisplay)
            {
                for (String propertyName : extraInfo.getColumnNames())
                {
                    String propertyValue =
                            extraInfo.convertToString(propertyName, getProperty(propertyName));
                    if (propertyValue == null)
                        propertyValue = "";
                    out.append("\t" + propertyValue);
                }
            }
        }
        //description
        out.append("\t" + (null == txt ? "" : "\"" + txt + "\""));

        return out.toString();
    }


    public Object clone()
    {
        return new Feature(this);
    }

    public int getScanCount()
    {
        if (scanCount == -1)
            return scanFirst - scanLast + 1;
        return scanCount;
    }

    public void setScanCount(int scanCount)
    {
        this.scanCount = scanCount;
    }

    public int getScanFirst()
    {
        return scanFirst;
    }

    public void setScanFirst(int scanFirst)
    {
        this.scanFirst = scanFirst;
    }

    public int getScanLast()
    {
        return scanLast;
    }

    public void setScanLast(int scanLast)
    {
        this.scanLast = scanLast;
    }

    public float getTotalIntensity()
    {
        return totalIntensity;
    }

    public void setTotalIntensity(float totalIntensity)
    {
        this.totalIntensity = totalIntensity;
    }


    /**
     * Needs to be defined when redefining equals. Seems likely to be unique...
     *
     * There are definitely pieces of code that depend on this method NOT hashing
     * features with the same scan, mz and intensity but different charge/KL to the same
     * value (FeatureStrategyPeakClusters).
     *
     * There may or may not be code out there that depends on this method hashing things
     * that agree on all these values, but have other values that are different, hashing
     * to the same value.
     *
     * TODO: make this agree with equals()
     */
    public int hashCode()
    {
        return scan ^ Float.floatToIntBits(mz) ^ (Float.floatToIntBits(intensity) << charge)
                ^ Float.floatToIntBits(kl);
    }

    /**
     * This equals() implementation is based only on the _known_ fields of Feature.
     *
     * So if two fields agree in all those fields but disagree in their extra info,
     * they will be considered equal.
     * @param o
     * @return
     */
    public boolean equals(Object o)
    {
        if (null == o || !(o instanceof Feature))
            return false;

        Feature feature = (Feature) o;
        return
                this.scan == feature.scan &&
                        this.mz == feature.mz &&
                        this.intensity == feature.intensity &&
                        this.charge == feature.charge &&
                        this.kl == feature.kl &&
                        this.dist == feature.dist &&
                        this.background == feature.background &&
                        this.median == feature.median &&
                        this.peaks == feature.peaks &&
                        this.skippedPeaks == feature.skippedPeaks &&
                        this.description == feature.description &&
                        this.chargeStates == feature.chargeStates &&
                        this.totalIntensity == feature.totalIntensity &&
                        this.scanFirst == feature.scanFirst &&
                        this.scanLast == feature.scanLast &&
                        this.scanCount == feature.scanCount &&
                        this.accurateMZ == feature.accurateMZ &&
                        this.mass == feature.mass;
    }


    //getters and setters for regular fields

    public void setMz(float mz)
    {
        this.mz = mz;
    }


    public float getMass()
    {
        return mass;
    }


    public void setMass(float mass)
    {
        // we need setMass as well as setMZ, for when we might load externally generated feature files
        this.mass = mass;
    }


    public int getCharge()
    {
        return charge;
    }

    public void setCharge(int charge)
    {
        this.charge = charge;
    }

    public float getKl()
    {
        return kl;
    }

    public void setKl(float kl)
    {
        this.kl = kl;
    }

    public int getPeaks()
    {
        return peaks;
    }

    public void setPeaks(int peaks)
    {
        this.peaks = peaks;
    }

    public String getDescription()
    {
        return description;
    }

    public void setDescription(String description)
    {
        this.description = description;
    }

    public int getChargeStates()
    {
        return chargeStates;
    }

    public void setChargeStates(int chargeStates)
    {
        this.chargeStates = chargeStates;
    }

    public float getTime()
    {
        return time;
    }

    public void setTime(float time)
    {
        this.time = time;
    }


    public float getSumSquaresDist()
    {
        return dist;
    }

    public void setSumSquaresDist(float dist)
    {
        this.dist = dist;
    }

    public boolean ContainsPeak(Spectrum.Peak p)
    {
        Spectrum.Peak[] peaks = comprised;
        if (null == peaks) return false;
        for (int i = 0; i < peaks.length; i++)
        {
            Spectrum.Peak peak = peaks[i];
            if (peak == p)
                return true;
        }
        return false;
    }

    public boolean isAccurateMZ()
    {
        return accurateMZ;
    }

    public void setAccurateMZ(boolean accurateMZ)
    {
        this.accurateMZ = accurateMZ;
    }

    /**
     * Compare features based on "quality".
     * "Quality" is currently defined as the number of peaks, with KL scores
     * breaking ties
     */
    public static class PeaksKLComparatorDesc implements Comparator<Feature>
    {
        protected static PeaksKLComparatorDesc singletonInstance = null;

        public static PeaksKLComparatorDesc getSingletonInstance()
        {
            if (singletonInstance == null)
                singletonInstance = new PeaksKLComparatorDesc();
            return singletonInstance;
        }

        public int compare(Feature o1, Feature o2)
        {
            int comp = (o1.getPeaks() == o2.getPeaks() ? 0 : o1.getPeaks() < o2.getPeaks() ? 1 : -1);
            if (comp != 0)
                return comp;
            comp = (o1.getKl() == o2.getKl() ? 0 : o1.getKl() < o2.getKl() ? 1 : -1);
            if (comp != 0)
                return comp;
            //if peaks and KL equals, return the "first one better" to break ties
            return -1;
        }
    }

    public static class IntensityDescComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            return _compareDesc(o1.intensity, o2.intensity);
        }
    }

    public static class MassAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            if (o1.mass == o2.mass)
                return 0;
            return _compareAsc(o1.mass, o2.mass);
        }
    }

    public static class MzDistAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            if (o1.mz == o2.mz)
                return _compareAsc(o1.kl, o2.kl);
            return _compareAsc(o1.mz, o2.mz);
        }
    }

    public static class MzAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            if (o1.mz == o2.mz)
                return 0;
            return _compareAsc(o1.mz, o2.mz);
        }
    }

    public static class MzScanAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            float mzDiff = o1.mz - o2.mz;

            if (0 == mzDiff)
                return _compareAsc(o1.scan, o2.scan);

            return mzDiff > 0 ? 1 : -1;
        }
    }

    public static class ScanAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            return _compareAsc(o1.scan, o2.scan);
        }
    }

    /**
     * This one is useful for associating equivalent features in multiple files with each other,
     * through equivalent ordering
     */
    public static class ScanChargeMzAscComparator implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            float scanDiff = o1.scan - o2.scan;

            if (0 == scanDiff)
            {
                float chargeDiff = o1.charge - o2.charge;
                if (0 == chargeDiff)
                    return _compareAsc(o1.mz, o2.mz);
                return chargeDiff > 0 ? 1 : -1;
            }

            return scanDiff > 0 ? 1 : -1;
        }
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
     * A helper method to get a FloatRange of a certain size around a point
     *
     * @param mz
     * @param mzWindow
     * @return a range of size mzWindow*2 around mz
     */
    protected static FloatRange getMZRange(float mz, int mzWindow)
    {
        return new FloatRange(mz - mzWindow, mz + mzWindow);
    }

    /**
     * dhmay adding 01/05/2006
     * Start with an array of scans "MSScan[] scans",
     * "_mzRange" set to a FloatRange indicating the
     * range of mz values to consider, and
     * a float "mz" set to the particular monoistopic mz
     * of the feature under consideration.
     *
     * @param run
     * @param mzWindowSize
     * @param scanRange
     * @param resamplingFrequency
     * @return calculated intensity
     */
    public float calculateFeatureIntensityInRange(MSRun run, int mzWindowSize,
                                                  IntRange scanRange,
                                                  int resamplingFrequency)
    {
        //if there are no scans in the range, there's no intensity
        if (scanRange == null)
            return 0;

        // This checks to see if the data have been "centroided". To save space and
        // make some downstream processing easier.
        boolean centroided = run.getHeaderInfo().getDataProcessing().getCentroided() == 1;

        int numScans = scanRange.getMaximumInteger() - scanRange.getMinimumInteger() + 1;
        MSRun.MSScan[] scans = run.getPartialScanArray(scanRange.getMinimumInteger(),
                scanRange.getMaximumInteger());

        assert(scans.length == numScans);

        float[][] spectra = new float[numScans][];
        FloatRange mzRange = getMZRange(mz, mzWindowSize);

        for (int i = 0; i < numScans; i++)
        {
            if (scans[i] == null)
                spectra[i] = new float[0];
            else
            {
                float[][] raw = scans[i].getSpectrum();
                // If the data has been centroided, we might try to make it a little easier
                // to resample and process.
// However, cleanSpectrum takes an ENORMOUS amount of time (more than everything else combined)
// AND it seems to throw things way off.  I do not like it, no sir.
//                if (centroided)
//                    raw = FeatureStrategyCentroided.cleanSpectrum(raw);

                spectra[i] = Spectrum.Resample(raw, mzRange, resamplingFrequency);
            }
        }

        // Spectrum.RemoveBackground() removes the "background" noise from spectra as a side-effect
        // Note: RemoveBackground requires that every entry in spectra be not null;
        float[][] background = Spectrum.RemoveBackground(spectra);

        // Convert mz to an index into our resample spectra array
        int imz = Math.round((mz - mzRange.min) * resamplingFrequency);

        // total intensity calculation
        // "integrate" intensity over time
        float totalIn = 0;
        for (int s = 0; s < spectra.length; s++)
        {
            //Note: this case is for missing or bad scans.  Not likely.
            if (spectra[s] == null || spectra[s].length < imz + 1)
            {
                continue;
            }

            double in = spectra[s][imz];
            double time = 1.0;

            //Note:  for first and last scans, time=1, since no scans on either side
            if (s > 0 && s + 1 < numScans && scans[s - 1] != null && scans[s + 1] != null)
                time = (scans[s + 1].getDoubleRetentionTime() -
                        scans[s - 1].getDoubleRetentionTime()) / 2;
            totalIn += time * in;
        }
        return totalIn;
    }

    /**
     * Helper method to find the overlap between the scan ranges of two features
     *
     * @param feature1
     * @param feature2
     * @return the overlap between the scan ranges of the two features.  If no overlap,
     *         returns null
     */
    public static IntRange findOverlappingScanRange(Feature feature1, Feature feature2)
    {
        int commonStartScan = Math.max(feature1.getScanFirst(), feature2.getScanFirst());
        int commonEndScan = Math.min(feature1.getScanLast(), feature2.getScanLast());

        if (commonEndScan < commonStartScan)
            return null;
        return new IntRange(commonStartScan, commonEndScan);
    }


    /**
     * Check to see if a given property is on the allowed list
     * This is no longer used -- allowing any property to be set
     *
     * @param propertyName
     * @return
     */
/*
    protected boolean isAllowedProperty(String propertyName)
    {
        for (int i = 0; i < _allowedProperties.length; i++)
        {
            if (_allowedProperties[i].equals(propertyName))
                return true;
        }
        return false;
    }
*/


    public boolean hasProperty(String propertyName)
    {
        return getPropertyMap().containsKey(propertyName);
    }


    public Map<String, Object> getPropertyMap()
    {
        //initialize the property map if it's not already created
        if (_propertyMap == null)
            _propertyMap = new HashMap<String, Object>();
        return _propertyMap;
    }

    /**
     * Set a property of this feature, initializing the property map if necessary.
     * If the value is of type String,
     *
     * @param propertyName
     * @param propertyValue
     */
    public void setProperty(String propertyName, Object propertyValue)
    {
        getPropertyMap().put(propertyName, propertyValue);
    }

    /**
     * Unset a property and return the old value if there was one
     *
     * @param propertyName
     */
    public Object unsetProperty(String propertyName)
    {
        Object oldValue = getProperty(propertyName);
        getPropertyMap().remove(propertyName);
        return oldValue;
    }


    /**
     * Indicates whether a given property is set for this feature
     * @param propertyName
     * @return
     */
    public boolean isPropertySet(String propertyName)
    {
        return getPropertyMap().containsKey(propertyName);
    }

    /**
     * Get a property of this feature.  If not set, return null
     *
     * @param propertyName
     * @return
     */
    public Object getProperty(String propertyName)
    {
        return getProperty(propertyName, null);
    }

    public Object get(Object propertyName)
    {
        return getProperty((String) propertyName);
    }

    /**
     * Get a property of this feature.  If not set, return defaultValue
     *
     * @param propertyName
     * @return
     */
    public Object getProperty(String propertyName, Object defaultValue)
    {
        if (_propertyMap == null)
            return defaultValue;
        Object result = _propertyMap.get(propertyName);
        if (result == null)
            result = defaultValue;
        return result;
    }

    public int getIntProperty(String propertyName, int defaultValue)
    {
        return (Integer) getProperty(propertyName, defaultValue);
    }

    public double getDoubleProperty(String propertyName, double defaultValue)
    {
        Object propertyVal = getProperty(propertyName, defaultValue);
        if (propertyVal != null && propertyVal instanceof Float)
            propertyVal = ((Float) propertyVal).doubleValue();
        return (Double) propertyVal;
    }

    public float getFloatProperty(String propertyName, float defaultValue)
    {
        Object propertyVal = getProperty(propertyName, defaultValue);
        if (propertyVal != null && propertyVal instanceof Double)
            propertyVal = ((Double) propertyVal).floatValue();
        return (Float) propertyVal;
    }

    public String getStringProperty(String propertyName, String defaultValue)
    {
        return (String) getProperty(propertyName, defaultValue);
    }


}
