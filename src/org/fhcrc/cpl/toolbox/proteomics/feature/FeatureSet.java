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

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.arguments.BooleanArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.*;
import org.apache.commons.beanutils.BeanUtils;
import org.apache.commons.lang.math.IntRange;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * User: migra
 * Date: Sep 13, 2004
 * Time: 3:05:56 PM
 */
public class FeatureSet implements Cloneable
{
    static Logger _log = Logger.getLogger(FeatureSet.class);

    public static final int FEATURE_FILE_FORMAT_MSINSPECTTSV = 0;
    public static final int FEATURE_FILE_FORMAT_APML = 1;
    public static final int FEATURE_FILE_FORMAT_HARDKLOR = 2;

    public static final int DEFAULT_FEATURE_FILE_FORMAT = FEATURE_FILE_FORMAT_MSINSPECTTSV;   



    //maintain loading status
    protected int _loadStatus = FEATURESET_LOAD_NOT_LOADED;
    protected String _loadStatusMessage = null;

    protected Feature[] _features;
    private Map<String, Object> _properties = new HashMap<String, Object>();

    private String _tag; // optional tag for this feature set
    private Color _color;
    private int _style = 0;
    private File _sourceFile;
    private boolean _displayed = true;


    //loading status codes
    public static final int FEATURESET_LOAD_NOT_LOADED = -1;
    public static final int FEATURESET_LOAD_SUCCESS = 0;
    public static final int FEATURESET_LOAD_ERROR_FILE_NOT_FOUND = 1;
    public static final int FEATURESET_LOAD_ERROR_BAD_FILE_FORMAT = 2;
    public static final int FEATURESET_LOAD_ERROR_NO_FEATURES_FOUND = 3;
    public static final int FEATURESET_LOAD_ERROR_UNKNOWN = 4;

    //different types of intensities to use in quantitation
    public static final int TOTAL_INTENSITY=0;
    public static final int MAX_INTENSITY=1;
    public static final int RECALCULATED_INTENSITY=2;
    //default to totalIntensity
    public static final int DEFAULT_INTENSITY_TYPE=TOTAL_INTENSITY;

    // How many daltons to extract from run on either side of each feature,
    // when comparing features
    static final int MZ_WINDOW_SIZE = 3;

    // Each scan contains a list of (mz, intensity) pairs. The m/z values are typically not sampled
    // regularly, so we often resample onto a regular grid before doing other processing. This value
    // sets the number of samples to provide per one m/z unit. The value 36 is typical throughout
    // msInspect (changing it could actually be a little tricky).
    static final int RESAMPLING_FREQUENCY = 36;

    //track the known extra information types for this FeatureSet
    protected List<FeatureExtraInformationDef> extraInformationTypes;




    //no-arg constructor used by feature-file handlers
    public FeatureSet()
    {

    }

    public FeatureSet(File file, Color color)
    {
        //initialize loading status values
        setLoadStatus(FEATURESET_LOAD_NOT_LOADED);
        setLoadStatusMessage("Features not yet loaded");

        //this used to throw a NullPointerException
        if (file == null)
        {
            setLoadStatus(FEATURESET_LOAD_ERROR_FILE_NOT_FOUND);
            setLoadStatusMessage("Error loading features: no file specified");
            return;
        }

        try
        {
            _sourceFile = file;
//          ApplicationContext.setMessage("Loading file " + file);
            loadFeatureFile(file);
            //check for load success before continuing
            if (getLoadStatus() == FEATURESET_LOAD_SUCCESS)
            {
                Arrays.sort(_features, new Feature.MzScanAscComparator());
                _color = color;
            }
        }
        catch (Exception e)
        {
            setLoadStatus(FEATURESET_LOAD_ERROR_UNKNOWN);
            _log.error("Feature-loading exception",e);
            setLoadStatusMessage("Unknown error loading features from file ");
        }
        if (getLoadStatus() != FEATURESET_LOAD_SUCCESS)
        {
            //all the non-success statuses require the filename appended 
            setLoadStatusMessage(getLoadStatusMessage() + file.getName());
        }
    }

    public FeatureSet(File file) throws Exception
    {
        this(file, null);
    }

    public FeatureSet(Feature[] features)
    {
        _features = features;
        inferExtraInformationTypesFromFeatures();
    }

    /**
     * For each feature, determine its set of extra information types,
     * using the property names set on the feature.
     *
     * Set the full list of information types on this FeatureSet to
     * that set.
     *
     * Return value is null so as to avoid confusion -- this does affect the set
     */
    public void inferExtraInformationTypesFromFeatures()
    {
        Set<FeatureExtraInformationDef> extraInfoSet =
                new HashSet<FeatureExtraInformationDef>();
        for (Feature feature : _features)
        {
            for (FeatureExtraInformationDef featureInfoDef :
                    feature.determineExtraInformationTypes())
            {
                extraInfoSet.add(featureInfoDef);
            }
        }
        removeAllExtraInformationTypes();
        for (FeatureExtraInformationDef infoDef : extraInfoSet)
        {
            addExtraInformationType(infoDef);
        }
    }


    public FeatureSet(Feature[] features, Color color)
    {
        this(features);
        _color = color;
    }


    public FeatureSet(Spectrum.Peak[] peaks, Color color)
    {
        _color = color;
        _features = new Feature[peaks.length];
        for (int i = 0; i < peaks.length; i++)
        {
            Spectrum.Peak p = peaks[i];
            _features[i] = new Feature(p.scan, p.scan, p.scan, p.mz, p.intensity, 1, 0.0F, p.intensity);
        }
    }

    //dmay adding 12/19/2005
    public int getLoadStatus()
    {
        return _loadStatus;
    }

    public void setLoadStatus(int loadStatus)
    {
        _loadStatus = loadStatus;
    }

    public String getLoadStatusMessage()
    {
        return _loadStatusMessage;
    }

    public void setLoadStatusMessage(String loadStatusMessage)
    {
        _loadStatusMessage = loadStatusMessage;
    }

    /**
     * Populate the Time attribute of every Feature in a FeatureSet
     * with the time corresponding to its scan.
     *
     * This is rather inefficient.  It would be much, much better if MSRun
     * kept ahold of a map of scan numbers to times.
     * //TODO: when we start getting pepXml files with retention times
        //TODO: populated, need to check and see if they're populated
        //TODO: before going back to run
     * @param run
     */
    public void populateTimesForMS2Features(MSRun run)
    {

//        Feature[] featuresCopy = new Feature[featureSet.getFeatures().length];
//        System.arraycopy(featureSet.getFeatures(), 0, featuresCopy,
//                         0, featuresCopy.length);
//        Arrays.sort(featuresCopy, new Feature.ScanAscComparator());
        for (Feature feature : getFeatures())
        {
            //first hoping that time is populated.
            //If not (if time <= 0), recalculating from scan
            //if this throws an NPE, so be it
            if (feature.getTime() <= 0)
            {
                int ms2ScanIndex =
                        run.getIndexForMS2ScanNum(feature.getScan());
                if (ms2ScanIndex < 0)
                    ms2ScanIndex = -ms2ScanIndex;
                MSRun.MSScan scan = run.getMS2Scan(ms2ScanIndex);
                feature.setTime((float)
                        scan.getDoubleRetentionTime());
            }
        }
    }

    /**
     * Populate the startTime and endTime properties of every Feature in a FeatureSet
     * with the times corresponding to its start and end scans
     * @param run
     */
    public void populateTimesForMS1Features(MSRun run)
    {
        for (Feature feature : getFeatures())
        {
            //first, set the time appropriately (if it's already set, believe that setting)
            //Then, check whether scanFirst and scanLast are the same as scan.  If so, set the
            //start and end times the same.  Otherwise, check the run to find the real ones out.
            //All this complexity is to avoid hitting the run, which is expensive.
            if (feature.getTime() <= 0)
                feature.setTime((float) run.getScan(run.getIndexForScanNum(feature.getScan())).getDoubleRetentionTime());
            if (feature.scanFirst == feature.scan)
                TimeExtraInfoDef.setStartTime(feature, feature.getTime());
            else
            {
                int fScanFirstIndex = run.getIndexForScanNum(feature.scanFirst);
                TimeExtraInfoDef.setStartTime(feature,
                        run.getScan(fScanFirstIndex).getDoubleRetentionTime());
            }
            if (feature.scanLast == feature.scan)
                TimeExtraInfoDef.setEndTime(feature, feature.getTime());
            else
            {
                int fScanLastIndex = run.getIndexForScanNum(feature.scanLast);
                TimeExtraInfoDef.setEndTime(feature,
                        run.getScan(fScanLastIndex).getDoubleRetentionTime());
            }
        }
    }



    /**
     *
     * Select a subset of a feature array.
     * Input array MUST be sorted with the Spectrum.Feature.MzScanAscComparator
     *
     * @param features
     * @param sel      What features to select. Features will not be merged
     * @return
     */
    public static Feature[] selectFeatures(Feature[] features, FeatureSelector sel)
    {
        ArrayList<Feature> list = new ArrayList<Feature>();
        for (Feature f : features)
        {
            if (f.intensity >= sel.getMinIntensity() &&
                    f.charge <= sel.getMaxCharge() && f.charge >= sel.getMinCharge() &&
                    f.mz >= sel.getMinMz() && f.mz <= sel.getMaxMz() &&
                    f.mass >= sel.getMinMass() && f.mass <= sel.getMaxMass() &&
                    f.peaks >= sel.getMinPeaks() && f.peaks <= sel.getMaxPeaks() &&    
                    f.scan >= sel.getScanFirst() && f.scan <= sel.getScanLast() &&
                    f.kl <= sel.getMaxKL() && f.scanCount >= sel.getMinScans() &&
                    f.totalIntensity >= sel.getMinTotalIntensity() &&
                    f.time >= sel.getMinTime() && f.time <= sel.getMaxTime() &&
                    MS2ExtraInfoDef.getPeptideProphet(f) >= sel.getMinPProphet() &&
                    (sel.getMaxAMTFDR() == 1 || (AmtExtraInfoDef.hasMatchFDR(f) &&
                            AmtExtraInfoDef.getMatchFDR(f) < sel.getMaxAMTFDR())) &&
                    (sel.getMaxMassDeviationPPM() == Integer.MAX_VALUE ||
                            Math.abs(MassCalibrationUtilities.calculateMassDefectDeviationPPM(f.getMass(),
                                    MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH)) <=
                                    sel.getMaxMassDeviationPPM()) &&
                    f.getSumSquaresDist() <= sel.getMaxSumSquaresDist() &&
                    (!sel.isAccurateMzOnly() || f.isAccurateMZ())
                    )
            {
                list.add(f);
            }
        }
        return list.toArray(new Feature[list.size()]);
    }

    public static Feature[] mergeFeatures(Feature[] features, FeatureSelector sel)
    {
        if (null == features)
            return null;

        ArrayList<Feature> featureRangeList;
        while (true) //Loop until we can't compress any more
        {
            int trail = 0;
            featureRangeList = new ArrayList<Feature>();
            for (int i = 0; i < features.length; i++)
            {
                Feature f = features[i];
                if (f.intensity >= sel.getMinIntensity() && f.charge <= sel.getMaxCharge() &&
                        f.mz >= sel.getMinMz() && f.mz <= sel.getMaxMz() && f.kl <= sel.getMaxKL()
                        && f.peaks >= sel.getMinPeaks()
                        && f.scan >= sel.getScanFirst() && f.scan <= sel.getScanLast()
                        //TODO: properly this should be handled in such a way that I don't have to
                        //reference this column name
                        && f.getIntProperty("peptideprophet",0) >= sel.getMinPProphet())
                {
                    boolean newRange = true;
                    //Scan all the feature ranges & merge.
                    int j;
                    for (j = trail; j < featureRangeList.size(); j++)
                    {
                        Feature fr = (Feature) featureRangeList.get(j);
                        float ppmGap = (f.mz - fr.mz);
                        if (ppmGap > sel.maxMzGap)
                        {
                            //Don't need to look at this next time since we're sorted.
                            trail = j + 1;
                            continue;
                        }

                        if (fr.isFeatureInRange(f, sel.maxScanGap, sel.maxMzGap))
                        {
                            fr.addFeatureToRange(f);
                            newRange = false;
                            break;
                        }

                        if (-ppmGap > sel.maxMzGap)
                            break;
                    }

                    if (newRange)
                        featureRangeList.add(new Feature(f));
                }

            }

            if (featureRangeList.size() == features.length)
                break; //OK, no more compressing can be done

            features = (Feature[]) featureRangeList.toArray(new Feature[featureRangeList.size()]);
            Arrays.sort(features, new Feature.MzScanAscComparator());
        }

        //Now filter by length
        ArrayList<Feature> filteredByLength = new ArrayList<Feature>();
        for (Feature fr : featureRangeList)
        {
            if (fr.getScanCount() >= sel.getMinScans())
                filteredByLength.add(fr);

        }

        //Re-sort because merging might lead to subtly out-of-order peaks
        Feature[] ranges = (Feature[]) filteredByLength.toArray(new Feature[filteredByLength.size()]);
        Arrays.sort(ranges, new Feature.MzScanAscComparator());
        return ranges;
    }

    public Feature[] getFeatures(FeatureSelector sel)
    {
        if (null == _features)
            return null;

        return selectFeatures(_features, sel);
    }

    public Feature[] getFeatures()
    {
        return _features;
    }

    public void setFeatures(Feature[] features)
    {
        _features = features;
    }

    public String getTag()
    {
        return _tag;
    }

    public void setTag(String tag)
    {
        _tag = tag;
    }

    public FeatureSet filter(FeatureSelector sel)
    {
        FeatureSet fs = (FeatureSet)this.clone();
        Feature[] features = getFeatures(sel);
        fs.setFeatures(features);
        fs.setColor(this.getColor());
        Map<String,Object> properties = new HashMap<String,Object>();
        properties.putAll(this.getProperties());
        if (null != this.getSourceFile())
            properties.put("origSourceFile", this.getSourceFile().getPath());
        properties.put("filter", sel.toString());

        fs.setProperties(properties);
        return fs;
    }

    public Feature findNearestFeature(int scan, float mz)
    {
        return findNearestFeature(scan, mz, Integer.MAX_VALUE, Float.MAX_VALUE);
    }

    public int findNearestFeatureIndex(int scan, float mz, int maxScanDistance, float maxMzDistance)
    {
        Feature feature = new Feature(scan, mz, 1);
        int index = Arrays.binarySearch(_features, feature, new Feature.MzScanAscComparator());
        double minDistance = Float.MAX_VALUE;
        int nearestFeature = -1;

        if (index >= 0)
            return index;

        int pos = -index - 1;

        //Didn't find an exact match. Search for nearest feature in 2 dimensions
        for (int i = pos - 1; i >= 0; i--)
        {
            Feature f = _features[i];
            float mzDist = Math.abs(mz - f.mz);
            if (mzDist > maxMzDistance || mzDist > minDistance)
                break;

            int scanDist = Math.abs(scan - f.scan);
            if (scanDist > maxScanDistance)
                continue;

            double distance = Math.sqrt(scanDist * scanDist + mzDist * mzDist);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestFeature = i;
            }
        }

        for (int i = pos; i < _features.length; i++)
        {
            Feature f = _features[i];
            float mzDist = Math.abs(mz - f.mz);
            if (mzDist > minDistance || mzDist > maxMzDistance) //Not going to get any closer since sorted by Mz
                return nearestFeature;

            int scanDist = Math.abs(scan - f.scan);
            if (scanDist > maxScanDistance)
                continue;

            double distance = Math.sqrt(scanDist * scanDist + mzDist * mzDist);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestFeature = i;
            }
        }

        return nearestFeature;
    }

    public Feature findNearestFeature(int scan, float mz, int maxScanDistance, float maxMzDistance)
    {
        int index = findNearestFeatureIndex(scan, mz, maxScanDistance, maxMzDistance);
        return index < 0 ? null : _features[index];
    }

    public Object clone()
    {
        try
        {
            FeatureSet fs = (FeatureSet) super.clone();
            fs.setProperties(this._properties);
            fs.setTag(this._tag);
            fs.setSourceFile(this._sourceFile); // Used for column headers
            return fs;
        }
        catch (CloneNotSupportedException x)
        {
            return null;
        }
    }

    /**
     * Deep-copies everything except the properties array.  The properties array
     * is shallow-copied.  This is done because we've done nothing to make sure that
     * the items in the properties array will clone() nicely.  Sheer laziness, really.
     * @return
     */
    public FeatureSet deepCopy()
    {
        Feature[] features = new Feature[_features.length];
        for (int i = 0; i < _features.length; i++)
            features[i] = (Feature) _features[i].clone();

        FeatureSet fs = new FeatureSet(features, _color);
        fs.setDisplayed(_displayed);
        fs.setSourceFile(new File(_sourceFile.getAbsolutePath()));

        //shallow copy of the properties array.
        HashMap propertiesCopy = new HashMap();
        for (Object propertyKey : _properties.keySet())
            propertiesCopy.put(propertyKey, _properties.get(propertyKey));

        fs.setProperties(propertiesCopy);

        return fs;
    }


    /**
     * Combine features that represent the same peptide.  Features
     * must be within a scanDiff,massDiff window to be considered same.
     * <p/>
     * This removes redunancy caused by multiple charge states.  It
     * also combines features of the same mass/charge that are close
     * together (e.g. a small feature that is split by ion competition)
     * <p/>
     * CONSIDER(mbellew) the mass (mz) tolerance could be tighter
     * for features w/ same charge
     *
     * @param scanDiff
     * @param massDiff
     * @param sumIntensities If this is true, resulting features have intensity and
     * totalIntensity values that are the sums of component features.  If false,
     * resulting features retain their original intensities
     * charge state
     * @return
     */
    public FeatureSet deconvolute(int scanDiff, double massDiff, boolean sumIntensities)
    {
        FeatureGrouper grouper = new FeatureGrouper();
        grouper.addSet(this);
        grouper.setGroupByMass(true);


        grouper.split2D(massDiff, scanDiff);
        Clusterer2D.BucketSummary[] buckets = grouper.summarize();
        Feature[] deconvoluted = new Feature[buckets.length];

        int numPeptideConflicts = 0;
        int numPreservedPeptides = 0;

        for (int i = 0; i < buckets.length; i++)
        {
            Clusterer2D.BucketSummary bucket = buckets[i];
            Feature deconvolutedFeature = null;
            if (bucket.featureCount == 1)
            {
                deconvolutedFeature = (Feature) FeatureGrouper.getFeatures(bucket)[0].clone();
                deconvolutedFeature.setChargeStates(1);
            }
            else
            {
                Feature[] bucketFeatures = FeatureGrouper.getFeatures(bucket);
                Feature best = bucketFeatures[0];
                float sumIntensity = 0.0f;
                float sumTotalIntensity = 0.0f;

                //NOTE:  this won't work with stuff of charge > 10.  Then again, what will?
                int[] chargeStateCounts = new int[10];

                String description = "";

                for (Feature f : bucketFeatures)
                {
                    if (description.length() > 0)
                        description += ", ";
                    chargeStateCounts[f.charge]++;
                    description += f.charge;
                    if (null != f.getDescription())
                        description += " (" + f.getDescription() + ")";
                    sumIntensity += f.intensity;
                    sumTotalIntensity += f.totalIntensity;

                    if (f.totalIntensity > best.totalIntensity)
                        best = f;
                }
                deconvolutedFeature = (Feature) best.clone();

                int numChargeStates = 0;
                for (int j=0; j<chargeStateCounts.length; j++)
                    if (chargeStateCounts[j] > 0)
                        numChargeStates++;

                if (sumIntensities)
                {
                    deconvolutedFeature.setIntensity(sumIntensity);
                    deconvolutedFeature.setTotalIntensity(sumTotalIntensity);
                }
                deconvolutedFeature.setChargeStates(numChargeStates);
                deconvolutedFeature.setDescription(description);


                //if there's MS2 data in this FeatureSet, then we need to make sure
                //that the collapsed feature contains the peptide and protein ID carried
                //by its components.
                //If there are conflicts, we leave the existing ID on the collapsed feature
                //alone, or if it had none, don't assign
                //TODO: somehow move this to MS2ExtraInfoDef?
                if (this.hasExtraInformationType(MS2ExtraInfoDef.getSingletonInstance()))
                {
                    Set<String> featurePeptides = new HashSet<String>();
                    Set<String> featureProteins = new HashSet<String>();

                    for (Feature f : bucketFeatures)
                    {
                        String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(f);
                        if (featurePeptide != null)
                        {
                            featurePeptides.add(featurePeptide);

                            String featureProtein = MS2ExtraInfoDef.getFirstProtein(f);
                            if (featureProtein != null)
                                featureProteins.add(featureProtein);
                        }
                    }

                    if (featurePeptides.size() == 1 &&
                            MS2ExtraInfoDef.getFirstPeptide(deconvolutedFeature) == null)
                    {
                        MS2ExtraInfoDef.setSinglePeptide(deconvolutedFeature,
                                featurePeptides.iterator().next());
                        numPreservedPeptides++;

                        if (featureProteins.size() == 1 &&
                                MS2ExtraInfoDef.getFirstProtein(deconvolutedFeature) == null)
                            MS2ExtraInfoDef.addProtein(deconvolutedFeature,
                                    featureProteins.iterator().next());
                    }
                    else
                    {
                        if (featurePeptides.size() > 1)
                            numPeptideConflicts++;
                    }
                }
            }
            deconvolutedFeature.comprised = FeatureGrouper.getFeatures(bucket);
            deconvoluted[i] = deconvolutedFeature;
        }

        //reporting on peptides preserved and conflicts
        if (this.hasExtraInformationType(MS2ExtraInfoDef.getSingletonInstance()))
        {
            _log.debug("deconvolute: peptides actively preserved: " + numPreservedPeptides);
            _log.debug("deconvolute: peptide conflicts: " + numPeptideConflicts);
        }

        FeatureSet fs = (FeatureSet) this.clone();
//Make map modifiable & set properties.
        Map props = new HashMap();
        props.putAll(this.getProperties());
        if (null != this.getSourceFile())
            props.put("origSourceFile", this.getSourceFile());
        props.put("deconvoluteScanDiff", String.valueOf(scanDiff));
        props.put("deconvoluteMassDiff", String.valueOf(massDiff));

        fs.setFeatures(deconvoluted);
        fs.setProperties(props);
        return fs;
    }


    /**
     * Relative quantitation using ICAT label defaults...
     */
    public FeatureSet icat()
    {
        return quant(AnalyzeICAT.icatLabel);
    }

    /**
     * Relative quantitiation using explicit label...
     */
    public FeatureSet quant(float light, float heavy, char residue, int maxLabelCount)
    {
        return quant(light, heavy, residue, maxLabelCount, null);
    }

    /**
     * Relative quantitiation using explicit label...
     */
    public FeatureSet quant(float light, float heavy, char residue, int maxLabelCount, MSRun run)
    {
        float delta = heavy - light;
        AnalyzeICAT.IsotopicLabel label = new AnalyzeICAT.IsotopicLabel(light, delta, residue, maxLabelCount);
        return quant(label, run);
    }

    public FeatureSet quant(AnalyzeICAT.IsotopicLabel label)
    {
        return quant(label, null);
    }

    //This comparator is used to sort pairs in increasing order of first scan.
    //This makes accessing the MS1 data faster, since data access is more localized
    static Comparator comparePairScanAsc = new Comparator()
    {
        public int compare(Object a, Object b)
        {
            Feature lightA = (Feature) ((Pair) a).first;
            Feature heavyA = (Feature) ((Pair) a).second;
            Feature lightB = (Feature) ((Pair) b).first;
            Feature heavyB = (Feature) ((Pair) b).second;
            IntRange aRange = Feature.findOverlappingScanRange(lightA, heavyA);
            IntRange bRange = Feature.findOverlappingScanRange(lightB, heavyB);
            int aScan = (aRange == null ? 0 : aRange.getMinimumInteger());
            int bScan = (bRange == null ? 0 : bRange.getMinimumInteger());
            return aScan == bScan ? 0 : aScan < bScan ? -1 : 1;
        }
    };

    //by default, use totalIntensity of each partner
    public FeatureSet quant(AnalyzeICAT.IsotopicLabel label, MSRun run)
    {
        return quant(label, TOTAL_INTENSITY, run);
    }

    /**
     * Relative quantitiation using explicit label...
     */
    public FeatureSet quant(AnalyzeICAT.IsotopicLabel label, int intensityType, MSRun run)
    {
        return quant(label, intensityType, run, AnalyzeICAT.DEFAULT_DELTA_MASS,
                AnalyzeICAT.DEFAULT_DELTA_MASS_TYPE, AnalyzeICAT.DEFAULT_DELTA_TIME);
    }

    public FeatureSet quant(AnalyzeICAT.IsotopicLabel label, int intensityType, MSRun run,
                            float massTolerance, int massToleranceType, float timeTolerance)
    {
        boolean pairsOnly = false;

        ArrayList pairs = AnalyzeICAT.analyze(getFeatures(), label, massTolerance,
                massToleranceType, timeTolerance);

        // TODO: option to output unpaired features
        Map<Feature,Feature> icatFeatures = new IdentityHashMap<Feature,Feature>(3 * pairs.size());

        ArrayList<Feature> list = new ArrayList<Feature>(_features.length);

        //
        // output paired features
        //
        // Consider May want to make ratio MAX_VALUE or MIN_VALUE when heavy is zero.
        //          +/- Inf (or Nan) may confuse downstream tools.


        // Give us the option of using either the total, recalculated, or maximum intensity
        // of each partner. Total has problems when one partner shows longer
        // elution time (or runs into a different co-eluting peptide).
        // TODO: This is just a surrogate for using the same range of scans
        // for each partner (max intensity has problems too).
        // Recalculated intensity requires access to the MS1 run
        if (intensityType == TOTAL_INTENSITY)
            Collections.sort(pairs,comparePairScanAsc);


        for (int i = 0; i < pairs.size(); i++)
        {
            Pair p = (Pair)pairs.get(i);
            Feature light = (Feature)p.first;
            Feature heavy = (Feature)p.second;
            Feature f = new Feature(light);
            f.setTotalIntensity(heavy.totalIntensity + light.totalIntensity);

            if (intensityType == TOTAL_INTENSITY)
            {
                IsotopicLabelExtraInfoDef.setHeavyIntensity(f,heavy.totalIntensity);
                IsotopicLabelExtraInfoDef.setLightIntensity(f,light.totalIntensity);
            }
            else if (intensityType == RECALCULATED_INTENSITY)
            {
                if (run == null)
                {
                    _log.error("No run specified, unable to recalculate intensities.");
                    return null;
                }
                IntRange overlappingScanRange = Feature.findOverlappingScanRange(light, heavy);
                IsotopicLabelExtraInfoDef.setLightIntensity(f,
                        light.calculateFeatureIntensityInRange(run, MZ_WINDOW_SIZE, overlappingScanRange,
                                RESAMPLING_FREQUENCY));
                IsotopicLabelExtraInfoDef.setHeavyIntensity(f,
                        heavy.calculateFeatureIntensityInRange(run, MZ_WINDOW_SIZE, overlappingScanRange,
                                RESAMPLING_FREQUENCY));
            }
            else if (intensityType == MAX_INTENSITY)
            {
                IsotopicLabelExtraInfoDef.setHeavyIntensity(f,heavy.intensity);
                IsotopicLabelExtraInfoDef.setLightIntensity(f,light.intensity);

                // smearing the total out doesn't seem to work better than just using the max
                // f.setHeavyIntensity(heavy.totalIntensity/heavy.scanCount);
                // f.setLightIntensity(light.totalIntensity/light.scanCount);
            }

            IsotopicLabelExtraInfoDef.setRatio(f,
                    IsotopicLabelExtraInfoDef.getLightIntensity(f) /
                            IsotopicLabelExtraInfoDef.getHeavyIntensity(f)); 
            IsotopicLabelExtraInfoDef.setLabelCount(f,
                    Math.round((heavy.mass - light.mass) / label.getHeavy()));
            f.setProperty("label", label);
            f.setChargeStates(Math.max(light.getChargeStates(), heavy.getChargeStates()));

//we used to do this in order to subtract out the light label.  Not doing that any more
//(mass will be consistent with m/z and charge), so no need to update mass at all            
//            f.updateMass();

            //deal with any peptide identifications, likely supplied by AMT.
            //If both light and heavy have the same ID, or any in common, or only one has an ID,
            //keep it.  If light and heavy have different IDs, toss them both out
            List<String> heavyPeptides = MS2ExtraInfoDef.getPeptideList(heavy);
            List<String> lightPeptides = MS2ExtraInfoDef.getPeptideList(heavy);
            if (heavyPeptides != null || lightPeptides != null)
            {
                if (heavyPeptides == null)
                {
                    MS2ExtraInfoDef.setPeptideList(f, lightPeptides.get(0));
                }
                else if (lightPeptides == null)
                {
                    MS2ExtraInfoDef.setPeptideList(f, heavyPeptides.get(0));
                }
                else
                {
                    //both heavy and light peptides exist.
                    if (heavyPeptides.size() == 1 && lightPeptides.size() == 1)
                    {
                        if (heavyPeptides.get(0).equals(lightPeptides.get(0)))
                            MS2ExtraInfoDef.setPeptideList(f, heavyPeptides.get(0));
                        else
                            MS2ExtraInfoDef.removeAllPeptides(f);
                    }
                    else
                    {
                        Set<String> commonPeptides = new HashSet<String>();
                        for (String heavyPeptide : heavyPeptides)
                            if (lightPeptides.contains(heavyPeptide))
                                commonPeptides.add(heavyPeptide);
                        if (commonPeptides.size() == 0)
                            MS2ExtraInfoDef.removeAllPeptides(f);
                        else
                            MS2ExtraInfoDef.setPeptideList(f, commonPeptides.iterator().next());
                    }
                }

                //now that we've figured out what peptide to assign, make sure it has the
                //right number of labeled residues.  If not, unset.
                if (MS2ExtraInfoDef.getFirstPeptide(f) != null)
                {
                    int numLabeledResidues = 0;
                    String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(f);
                    for (int j=0; j<featurePeptide.length(); j++)
                        if (featurePeptide.charAt(j) == label.getResidue())
                            numLabeledResidues++;
                    if (numLabeledResidues != IsotopicLabelExtraInfoDef.getLabelCount(f))
                    {
//                        if (numLabeledResidues > 0) System.err.println("Tossing: " + featurePeptide + ", " + numLabeledResidues + ", " + IsotopicLabelExtraInfoDef.getLabelCount(f));

                        MS2ExtraInfoDef.removeAllPeptides(f);
                    }
//                    else
//                        System.err.println("Saving: " + featurePeptide + ", " + numLabeledResidues);
                }
            }

            list.add(f);
            icatFeatures.put(light, light);
            icatFeatures.put(heavy, heavy);

        }


        //
        // output remaining features
        //

        if (!pairsOnly)
        {
            for (int i=0 ; i<_features.length ; i++)
            {
                Feature f = _features[i];
                if (icatFeatures.containsKey(f))
                    continue;
                list.add(new Feature(f));
            }
        }

        FeatureSet fs = (FeatureSet)this.clone();
        fs.addExtraInformationType(new IsotopicLabelExtraInfoDef());
        fs.getProperties().put("label", label.toString());
        fs.setFeatures(list.toArray(new Feature[0]));
        return fs;
    }


    /**
     * Return the best hit from a List of FeatureSets
     */
    public static Feature hitTest(java.util.List featureSets, int x, float y, int maxScan, float maxMz)
    {
        Feature feature = null;
        double minDistance = Double.MAX_VALUE;

        if (null == featureSets)
            return null;

        for (int i = 0; i < featureSets.size(); i++)
        {
            FeatureSet fs = (FeatureSet)featureSets.get(i);
            Feature feat = fs.findNearestFeature(x, (float)y, maxScan, maxMz);
            if (null == feat)
                continue;
            double distance = Math.sqrt(Math.pow(feat.scan - x, 2) + Math.pow(feat.mz - y, 2));
            if (distance < minDistance)
                feature = feat;
        }

        return feature;

    }

    /*
    * @param	file	the file containing features to load
    */
    public void loadFeatureFile(File file) throws Exception
    {
        //first check if the file exists
        if (file == null || !file.exists())
        {
            setLoadStatus(FEATURESET_LOAD_ERROR_FILE_NOT_FOUND);
            setLoadStatusMessage("Error loading features: unable to find file ");
            return;
        }

        FeatureSetFileHandler fileHandler = null;

        if (PepXMLFeatureFileHandler.getSingletonInstance().canHandleFile(file))
        {
            //try loading it as a pepXML file
            _log.debug("Loading as PepXML file");
            fileHandler = PepXMLFeatureFileHandler.getSingletonInstance();
        }
        else if (APMLFeatureFileHandler.getSingletonInstance().canHandleFile(file))
        {
            //try loading it as an APML file
            _log.debug("Loading as APML file");
            fileHandler = APMLFeatureFileHandler.getSingletonInstance();
        }
        else if (HardklorFeatureFileHandler.getSingletonInstance().canHandleFile(file))
        {
            //try loading it as a Hardklor file
            _log.debug("Loading as Hardklor file");
            fileHandler = HardklorFeatureFileHandler.getSingletonInstance();
        }
        else if (NativeTSVFeatureFileHandler.getSingletonInstance().canHandleFile(file))
        {
            //if not an xml file, assume tab-separated value file
            _log.debug("Loading as msInspect .tsv file");            
            fileHandler = NativeTSVFeatureFileHandler.getSingletonInstance();
        }
        else
        {
            throw new IllegalArgumentException("Unknown feature file type.  Doesn't seem to be APML, pepXML or msInspect TSV file.  Quitting.");
        }

        try
        {
            FeatureSet loadedFeatureSet = fileHandler.loadFeatureSet(file);

            //This is a bit cumbersome: load up the file in the handler, then take
            //the resulting FeatureSet and copy all the important stuff here.
            _features = loadedFeatureSet.getFeatures();
            _log.debug("Loaded " + _features.length + " features from file");
            setProperties(loadedFeatureSet._properties);
            for(FeatureExtraInformationDef infoType : loadedFeatureSet.getExtraInformationTypes())
                addExtraInformationType(infoType);
            setTag(loadedFeatureSet._tag);
            setSourceFile(file);

            //if we got here, load was successful
            setLoadStatus(FeatureSet.FEATURESET_LOAD_SUCCESS);
            setLoadStatusMessage(_features.length + " Features loaded successfully");
        }
        catch (IOException e)
        {
            //in case of Exceptions, assume no features found,
            //problem with file format
            _log.error("User attempted to load bad feature file, filename = " + file.getName() +
                    ", exception message = " + e.getMessage(),e);
            setLoadStatus(FeatureSet.FEATURESET_LOAD_ERROR_BAD_FILE_FORMAT);
            setLoadStatusMessage("Error loading features: bad file format in file ");
        }

        if (getLoadStatus() != FEATURESET_LOAD_SUCCESS)
            return;

        //"success" case
        //if no features found, report
        if (_features == null || _features.length == 0)
        {
            _log.info("User attempted to load file with no features");
            setLoadStatus(FEATURESET_LOAD_ERROR_NO_FEATURES_FOUND);
            setLoadStatusMessage("Error loading features: no features found in file ");
            return;
        }
    }

    public Color getColor()
    {
        return null == _color ? Color.RED : _color;
    }

    public void setColor(Color color)
    {
        _color = color;
    }

    public int getStyle()
    {
        return _style;
    }

    public void setStyle(int style)
    {
        _style = style;
    }

    public File getSourceFile()
    {
        return _sourceFile;
    }

    public void setSourceFile(File sourceFile)
    {
        _sourceFile = sourceFile;
    }

    public boolean isDisplayed()
    {
        return _displayed;
    }

    public void setDisplayed(boolean displayed)
    {
        _displayed = displayed;
    }

    public Map<String,Object> getProperties()
    {
        return _properties;
    }

    public Object getProperty(String propertyName)
    {
        if (null == _properties)
            return null;
        return _properties.get(propertyName);
    }

    public void setProperties(Map<String,Object> properties)
    {
        if (null == properties)
            _properties = new HashMap<String,Object>();
        else
            _properties = new HashMap<String,Object>(properties);
    }

    public void setProperty(String propertyName, Object propertyValue)
    {
        if (null == _properties)
            _properties = new HashMap<String,Object>();
        _properties.put(propertyName, propertyValue);
        if (_log.isDebugEnabled())
        {
            String className = null;
            if (propertyValue != null)
                className = propertyValue.getClass().getName();
            _log.debug("setProperty: " + propertyName + ", class " + className);
        }
    }


    public void save() throws IOException
    {
        save(getSourceFile());
    }

    public void save(File file) throws IOException
    {
        save(file, false);
    }

    public void save(File file, boolean dumpWindow) throws IOException
    {
        save(file, dumpWindow, NativeTSVFeatureFileHandler.FILE_TYPE_NAME);
    }

    public void save(PrintWriter outPW, boolean dumpWindow, String fileType) throws IOException
    {
        FeatureSetFileHandler fileHandler = null;

        if (APMLFeatureFileHandler.FILE_TYPE_NAME.equals(fileType))
        {
            fileHandler = new APMLFeatureFileHandler();
        }
        else if (HardklorFeatureFileHandler.FILE_TYPE_NAME.equals(fileType))
        {
            fileHandler = new HardklorFeatureFileHandler();
        }
        else
        {
            fileHandler = new NativeTSVFeatureFileHandler();
        }
        fileHandler.setDumpWindow(dumpWindow);
        fileHandler.saveFeatureSet(this, outPW);
    }

    public void save(File outFile, boolean dumpWindow, String fileType) throws IOException
    {
        FeatureSetFileHandler fileHandler = null;

        if (APMLFeatureFileHandler.FILE_TYPE_NAME.equals(fileType))
        {
            fileHandler = new APMLFeatureFileHandler();
        }
        else if (HardklorFeatureFileHandler.FILE_TYPE_NAME.equals(fileType))
        {
            fileHandler = new HardklorFeatureFileHandler();
        }
        else
        {
            fileHandler = new NativeTSVFeatureFileHandler();
        }
        fileHandler.setDumpWindow(dumpWindow);
        fileHandler.saveFeatureSet(this, outFile);
    }

    public void save(PrintWriter out)
    {
        save(out, false);
    }

    public void save(PrintWriter out, boolean dumpWindow)
    {
        NativeTSVFeatureFileHandler tsvFileHandler = new NativeTSVFeatureFileHandler();
        tsvFileHandler.setDumpWindow(dumpWindow);
        tsvFileHandler.saveFeatureSet(this, out);
    }


    public void savePepXml(File outFile)
            throws IOException
    {
        savePepXml(outFile, 1);
    }

    /**
     * For saving in pepXml format
     * @param outFile
     */
    public void savePepXml(File outFile, int firstSpectrumQueryIndex)
            throws IOException
    {
        PepXMLFeatureFileHandler pepXmlFileHandler = new PepXMLFeatureFileHandler();
        pepXmlFileHandler.setFirstSpectrumQueryIndex(firstSpectrumQueryIndex);
        pepXmlFileHandler.saveFeatureSet(this, outFile);
    }

    public static class FeatureSelector implements Cloneable
    {
        int minCharge = -10;
        int maxCharge = 10;
        // UNDONE: default these ranges based on getLowMz(), getHighMz()
        float maxMz = 10000;
        float minMz = 0f;
        float minIntensity = 0f;
        private float minTotalIntensity = 0f;
        int minScans = 0;
        int scanFirst = 0;
        int scanLast = Integer.MAX_VALUE;
        double maxKL = Double.MAX_VALUE;
        int minPeaks = 0;
        int maxPeaks = Integer.MAX_VALUE;        
        float minMass = 0f;
        float maxMass = Float.MAX_VALUE;
        float minTime = 0f;
        float maxTime = Float.MAX_VALUE;
        int maxMassDeviationPPM = Integer.MAX_VALUE;
        //dhmay adding 2/25/2007
        float maxSumSquaresDist = Float.MAX_VALUE;

        float minPProphet = 0;
        float maxAMTFDR = 1f;

        int maxScanGap = 3;
        float maxMzGap = .12f;

        //dhmay adding 2/25/2009
        boolean accurateMzOnly = false;

        public boolean equals(Object o)
        {
            if (null == o || !(o instanceof FeatureSelector))
                return false;

            FeatureSelector fs = (FeatureSelector) o;
            return getMinCharge() == fs.getMinCharge() && getMaxCharge() == fs.getMaxCharge() &&
                    getMaxMz() == fs.getMaxMz() && getMinMz() == fs.getMinMz() &&
                    getMinIntensity() == fs.getMinIntensity() && getMinScans() == fs.getMinScans() &&
                    maxScanGap == fs.maxScanGap && maxMzGap == fs.maxMzGap && getScanFirst() == fs.getScanFirst() &&
                    getScanLast() == fs.getScanLast() &&
                    getMaxKL() == fs.getMaxKL() && getMinPeaks() == fs.getMinPeaks()
                     && getMaxPeaks() == fs.getMaxPeaks()
                    && getMinTotalIntensity() == fs.getMinTotalIntensity()
                    && getMinMass() == fs.getMinMass() && getMaxMass() == fs.getMaxMass()
                    && getMinTime() == fs.getMinTime() && getMaxTime() == fs.getMaxTime()
                    && getMinPProphet() == fs.getMinPProphet() &&
                    getMaxMassDeviationPPM() == fs.getMaxMassDeviationPPM() &&
                    getMaxSumSquaresDist() == fs.getMaxSumSquaresDist() &&
                    isAccurateMzOnly() == fs.isAccurateMzOnly() &&
                    fs.getMaxAMTFDR() == getMaxAMTFDR();
        }

        public String toString()
        {
            StringBuffer sb = new StringBuffer();
            FeatureSelector unchanged = new FeatureSelector();
            String[] props = new String[] {"minCharge", "maxCharge", "minMz", "maxMz", "minMass", "maxMass", "minIntensity", "minTotalIntensity", "maxKL",
                    "minPeaks", "maxPeaks", "scanFirst", "scanLast", "minTime", "maxTime", "minScans", "minPProphet",
                    "maxMassDeviationPPM","maxSumSquaresDist","accurateMzOnly","maxAMTFDR"};
            try
            {
                for (String prop : props)
                {
                    Object val = BeanUtils.getSimpleProperty(this, prop);
                    Object orig = BeanUtils.getSimpleProperty(unchanged, prop);
                    if (!val.equals(orig))
                        sb.append(" --" + prop + "=" + val.toString());
                }
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage("FeatureSelector: ", x);
            }
            return sb.toString();
        }

        public boolean setFilterParam(String paramName, String paramVal)
        {
            //TODO: Should use reflection here since all names match
            if ("--minMz".equalsIgnoreCase(paramName))
                setMinMz(Float.parseFloat(paramVal));
            else if ("--maxMz".equalsIgnoreCase(paramName))
                setMaxMz(Float.parseFloat(paramVal));
            else if ("--minMass".equalsIgnoreCase(paramName))
                setMinMass(Float.parseFloat(paramVal));
            else if ("--maxMass".equalsIgnoreCase(paramName))
                setMaxMass(Float.parseFloat(paramVal));
            else if ("--minCharge".equalsIgnoreCase(paramName))
                setMinCharge(Integer.parseInt(paramVal));
            else if ("--maxCharge".equalsIgnoreCase(paramName))
                setMaxCharge(Integer.parseInt(paramVal));
            else if ("--minPeaks".equalsIgnoreCase(paramName))
                setMinPeaks(Integer.parseInt(paramVal));
            else if ("--maxPeaks".equalsIgnoreCase(paramName))
                setMaxPeaks(Integer.parseInt(paramVal));
            else if ("--minScanCount".equalsIgnoreCase(paramName))
                setMinScans(Integer.parseInt(paramVal));
            else if ("--scanFirst".equalsIgnoreCase(paramName))
                setScanFirst(Integer.parseInt(paramVal));
            else if ("--scanLast".equalsIgnoreCase(paramName))
                setScanLast(Integer.parseInt(paramVal));
            else if ("--minScans".equalsIgnoreCase(paramName))
                setMinScans(Integer.parseInt(paramVal));
            else if ("--maxKL".equalsIgnoreCase(paramName))
                setMaxKL(Double.parseDouble(paramVal));
            else if ("--minIntensity".equalsIgnoreCase(paramName))
                setMinIntensity(Float.parseFloat(paramVal));
            else if ("--minTime".equalsIgnoreCase(paramName))
                setMinTime(Float.parseFloat(paramVal));
            else if ("--maxTime".equalsIgnoreCase(paramName))
                setMaxTime(Float.parseFloat(paramVal));
            else if ("--minTotalIntensity".equalsIgnoreCase(paramName))
                setMinTotalIntensity(Float.parseFloat(paramVal));
            else if ("--minPProphet".equalsIgnoreCase(paramName))
                setMinPProphet(Float.parseFloat(paramVal));
            else if ("--maxMassDeviationPPM".equalsIgnoreCase(paramName))
                setMaxMassDeviationPPM(Integer.parseInt(paramVal));
            else if ("--maxSumSquaresDist".equalsIgnoreCase(paramName))
                setMaxSumSquaresDist(Integer.parseInt(paramVal));
            else if ("--accMzOnly".equalsIgnoreCase(paramName))
                try
                {
                    setAccurateMzOnly((Boolean) new BooleanArgumentDefinition("dummy").convertArgumentValue(paramVal));
                }
                catch (ArgumentValidationException e)
                {

                }
            else if ("--maxamtfdr".equalsIgnoreCase(paramName))
                setMaxAMTFDR(Float.parseFloat(paramVal));
            else
                return false;

            return true;
        }


        public Object clone()
        {
            try
            {
                return super.clone();
            }
            catch (Exception x)
            {
                //Impossible
                return null;
            }
        }


        public int getMaxMassDeviationPPM()
        {
            return maxMassDeviationPPM;
        }

        public void setMaxMassDeviationPPM(int maxMassDeviationPPM)
        {
            this.maxMassDeviationPPM = maxMassDeviationPPM;
        }

        public int getMinCharge()
        {
            return minCharge;
        }

        public void setMinCharge(int minCharge)
        {
            this.minCharge = minCharge;
        }

        public int getMaxCharge()
        {
            return maxCharge;
        }

        public void setMaxCharge(int maxCharge)
        {
            this.maxCharge = maxCharge;
        }

        public float getMaxMz()
        {
            return maxMz;
        }

        public void setMaxMz(float maxMz)
        {
            this.maxMz = maxMz;
        }

        public float getMinMz()
        {
            return minMz;
        }

        public void setMinMz(float minMz)
        {
            this.minMz = minMz;
        }

        public float getMinPProphet()
        {
            return minPProphet;
        }

        public void setMinPProphet(float minPProphet)
        {
            this.minPProphet = minPProphet;
        }


        public float getMinIntensity()
        {
            return minIntensity;
        }

        public void setMinIntensity(float minIntensity)
        {
            this.minIntensity = minIntensity;
        }

        public int getMinScans()
        {
            return minScans;
        }

        public void setMinScans(int minScans)
        {
            this.minScans = minScans;
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

        public double getMaxKL()
        {
            return maxKL;
        }

        public void setMaxKL(double maxKL)
        {
            this.maxKL = maxKL;
        }

        public int getMinPeaks()
        {
            return minPeaks;
        }

        public void setMinPeaks(int minPeaks)
        {
            this.minPeaks = minPeaks;
        }

        public int getMaxPeaks()
        {
            return maxPeaks;
        }

        public void setMaxPeaks(int maxPeaks)
        {
            this.maxPeaks = maxPeaks;
        }

        public float getMinTotalIntensity()
        {
            return minTotalIntensity;
        }

        public void setMinTotalIntensity(float minTotalIntensity)
        {
            this.minTotalIntensity = minTotalIntensity;
        }

        public float getMinMass()
        {
            return minMass;
        }

        public void setMinMass(float minMass)
        {
            this.minMass = minMass;
        }

        public float getMaxMass()
        {
            return maxMass;
        }

        public void setMaxMass(float maxMass)
        {
            this.maxMass = maxMass;
        }

        public float getMinTime()
        {
            return minTime;
        }

        public void setMinTime(float minTime)
        {
            this.minTime = minTime;
        }

        public float getMaxTime()
        {
            return maxTime;
        }

        public void setMaxTime(float maxTime)
        {
            this.maxTime = maxTime;
        }


        public float getMaxSumSquaresDist()
        {
            return maxSumSquaresDist;
        }

        public void setMaxSumSquaresDist(float maxSumSquaresDistance)
        {
            this.maxSumSquaresDist = maxSumSquaresDistance;
        }

        public boolean isAccurateMzOnly()
        {
            return accurateMzOnly;
        }

        public void setAccurateMzOnly(boolean accurateMzOnly)
        {
            this.accurateMzOnly = accurateMzOnly;
        }

        public float getMaxAMTFDR()
        {
            return maxAMTFDR;
        }

        public void setMaxAMTFDR(float maxAMTFDR)
        {
            this.maxAMTFDR = maxAMTFDR;
        }
    }

    public List<FeatureExtraInformationDef> getExtraInformationTypes()
    {
        if (extraInformationTypes == null)
        {
            extraInformationTypes = new ArrayList<FeatureExtraInformationDef>();
        }
        return extraInformationTypes;
    }

    public FeatureExtraInformationDef[] getExtraInformationTypesArray()
    {
        return getExtraInformationTypes().toArray(new FeatureExtraInformationDef[0]);
    }

    public void addExtraInformationType(FeatureExtraInformationDef infoType)
    {
        if (!getExtraInformationTypes().contains(infoType))
            getExtraInformationTypes().add(infoType);
    }

    public boolean hasExtraInformationType(FeatureExtraInformationDef infoType)
    {
        return getExtraInformationTypes().contains(infoType);
    }

    public void removeAllExtraInformationTypes()
    {
        extraInformationTypes = new ArrayList<FeatureExtraInformationDef>();
    }



}
