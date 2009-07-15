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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandModuleUtilities;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.Window2DFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

import java.io.File;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class FlagQuantEventsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FlagQuantEventsCLM.class);

    protected File[] featureFiles;
    protected File outFile;
    protected File outDir;
    protected File mzXmlDir;

    public static final int FLAG_REASON_NONE = 0;
    public static final int FLAG_REASON_COELUTING = 1;
    public static final int FLAG_REASON_DISSIMILAR_KL = 2;
    public static final int FLAG_REASON_DISSIMILAR_MS1_RATIO = 3;
    public static final int FLAG_REASON_DIFF_PEAKDIST = 4;


    protected int numScansSlopForCoeluting = 2;
    protected float massSlopPPM = 25f;
    protected float massSlopDa = 0.05f;

    //Maximum allowed proportion of highest peak intensity that the peak 1Da below monoisotope can have
    public static final float DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF = 0.5f;
    protected float peakBelowIntensityRatioCutoff = DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF;

    public static final int DEFAULT_NUM_SCANS_AROUND_EVENT = 2;
    protected int numScansAroundEventToConsider = DEFAULT_NUM_SCANS_AROUND_EVENT;

    public static final float DEFAULT_PEAK_PPM_TOLERANCE = 25;
    protected float peakPPMTolerance = DEFAULT_PEAK_PPM_TOLERANCE;


    //todo: move
    public static final float MOD_MASS_SLOP = 0.1f;



    protected int labelType = LABEL_LYCINE;

    protected boolean showCharts = false;

    public static final int LABEL_ACRYLAMIDE = 0;
    public static final int LABEL_LYCINE = 1;

    public static final float SILAC_LABEL_MASS = 134.115092f;
    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;


    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;


    protected String[] labelStrings = new String[]
            {
                    "acrylamide",
                    "silac"
            };

    protected String[] labelExplanations = new String[]
            {
                    "Acrylamide (3.0106Da on C)",
                    "SILAC Lycine labeling (134.115092 on K, i.e. 6Da SILAC)"
            };

    public static final String[] flagReasonDescs = new String[]
            {
                    "None",
                    "Coeluting peptide",
                    "Dissimilar KL",
                    "MS1 Ratio different from MS2",
                    "Light/Heavy have different peak distributions",
            };

    protected int maxDaltonsUpCoelute = 3;
    protected int maxDaltonsDownCoelute = 3;





//scaffolding
    private List<Float> daltonsOff = new ArrayList<Float>();


    public FlagQuantEventsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "flagquant";
        mShortDescription = "Flag questionable quantitation events";
        mHelpMessage = "Flag questionable quantitation events";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "MS2 Feature files"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "mzXML directory"),
                        new FileToWriteArgumentDefinition("out", false, "Output file"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
                       new EnumeratedValuesArgumentDefinition("label", false, labelStrings, labelExplanations,
                               "silac"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (featureFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }

        if (outFile != null)
            assertArgumentAbsent("outdir");
        if (outDir != null)
            assertArgumentAbsent("out");

        mzXmlDir = getFileArgumentValue("mzxmldir");

        labelType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label"));

        showCharts = getBooleanArgumentValue("showcharts");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File featureFile : featureFiles)
        {
            File outputFile = outFile;
            if (outputFile == null)
                outputFile = new File(outDir, featureFile.getName());

            Iterator<FeatureSet> featureSetIterator = null;
            try
            {
                featureSetIterator =
                        new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failure loading pepXML file " + featureFile.getAbsolutePath(),e);
            }

            int flaggedFeaturesThisFile = 0;
            int totalFeaturesThisFile = 0;

            List<File> tempFeatureFiles = new ArrayList<File>();
            int numSetsProcessed = 0;
            while (featureSetIterator.hasNext())
            {
                FeatureSet featureSet = featureSetIterator.next();
                totalFeaturesThisFile += featureSet.getFeatures().length;

                MSRun run = null;
                try
                {
                    File featureSetFile = featureSet.getSourceFile();
                    if (MS2ExtraInfoDef.getFeatureSetBaseName(featureSet) != null)
                        featureSetFile = new File(MS2ExtraInfoDef.getFeatureSetBaseName(featureSet) + ".pep.xml");
                    File mzXmlFile = ViewerCommandModuleUtilities.findCorrespondingMzXmlFile(
                            featureSetFile, mzXmlDir);
                    ApplicationContext.infoMessage("Loading mzXml file " + mzXmlFile.getAbsolutePath());
                    run = MSRun.load(mzXmlFile.getAbsolutePath());
                    ApplicationContext.infoMessage("Loaded.");
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Can't find or open mzXml file for file " +
                            featureFile.getAbsolutePath(), e);
                }


                ApplicationContext.infoMessage("\tProcessing fraction " + (numSetsProcessed+1) + "...");

                String baseName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
                if (baseName == null)
                {
                    baseName = featureFile.getName();
                    if (numSetsProcessed > 0 || featureSetIterator.hasNext())
                        baseName = baseName + "_" + numSetsProcessed;
                }
                //if PeptideProphet was run from a directory below the directory containing the
                //mzXML files, we may have ../ in the baseName, which causes trouble in saving
                //the temporary files
                while (baseName.contains(".." + File.separator))
                    baseName.replaceFirst(".." + File.separator, "");

                File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);

                processFeatureSet(featureSet, run);

                flaggedFeaturesThisFile += featureSet.getFeatures().length;

                try
                {
                    featureSet.savePepXml(thisFractionFile);
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failure writing file " +
                            thisFractionFile.getAbsolutePath(),e);
                }

                _log.debug("Saved fraction file as " + thisFractionFile.getAbsolutePath());
                tempFeatureFiles.add(thisFractionFile);

                numSetsProcessed++;
            }

            ApplicationContext.infoMessage("Flagged " + flaggedFeaturesThisFile + " out of " + totalFeaturesThisFile +
                       " features (" + ((float) flaggedFeaturesThisFile *100f / (float) totalFeaturesThisFile) +
                       "%) in this file");

            ApplicationContext.infoMessage("Saving output file " + outputFile.getAbsolutePath() + "...");

            try
            {
                if (numSetsProcessed == 1)
                {
                    FileReader in = new FileReader(tempFeatureFiles.get(0));
                    FileWriter out = new FileWriter(outputFile);
                    int c;

                    while ((c = in.read()) != -1)
                        out.write(c);

                    in.close();
                    out.close();
                }
                else
                {
                    ApplicationContext.infoMessage("\tCombining individual fraction files... " +
                            outputFile.getAbsolutePath() + "...");
                    new PepXMLFeatureFileHandler().combinePepXmlFiles(tempFeatureFiles, outputFile);
                }
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save output file " + outputFile.getAbsolutePath(), e);
            }
            ApplicationContext.infoMessage("Done.");
        }
    }

    /**
     * Side effect: ms2FeatureSet retains only flagged features
     * @param ms2FeatureSet
     * @param run
     * @throws CommandLineModuleExecutionException
     */
    protected void processFeatureSet(FeatureSet ms2FeatureSet, MSRun run)
            throws CommandLineModuleExecutionException
    {
        List<Feature> flaggedFeatures = new ArrayList<Feature>();

        Window2DFeatureSetMatcher fsm = new Window2DFeatureSetMatcher();
        fsm.setMatchWithinChargeOnly(true);
        fsm.setMassDiffType(FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);
        fsm.setMinMassDiff(-maxDaltonsUpCoelute - massSlopDa);
        fsm.setMaxMassDiff(maxDaltonsDownCoelute + massSlopDa);

        fsm.setElutionMode(Window2DFeatureSetMatcher.ELUTION_MODE_SCAN);
        fsm.setElutionRangeMode(Window2DFeatureSetMatcher.ELUTION_RANGE_MODE_RANGE);
        fsm.setMinElutionDiff(-numScansSlopForCoeluting);
        fsm.setMaxElutionDiff(numScansSlopForCoeluting);

        Map<Feature, Feature> origOppositeFeatureMap = new HashMap<Feature, Feature>();
        List<Feature> origAndOppositeFeatures = new ArrayList<Feature>();
        Map<Feature, Boolean> featureIsLightMap = new HashMap<Feature, Boolean>();
        for (Feature feature : ms2FeatureSet.getFeatures())
        {
            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                continue;
            Feature oppFeature = (Feature) feature.clone();
            boolean isLight = false;
            if (isHeavyLabeled(feature))
                oppFeature.setMass(feature.getMass() - calcHeavyLightMassDiff(feature));
            else if (isLightLabeled(feature))
            {
                isLight = true;
                oppFeature.setMass(feature.getMass() + calcHeavyLightMassDiff(feature));
            }
            else throw new RuntimeException("Neither light nor heavy labeled?!\n" + feature.toString());
            oppFeature.updateMz();

            featureIsLightMap.put(feature, isLight);
            origOppositeFeatureMap.put(feature, oppFeature);
            origAndOppositeFeatures.add(feature);
            origAndOppositeFeatures.add(oppFeature);
        }

        //WARNING! messing with ms2FeatureSet features here
        ms2FeatureSet.setFeatures(origAndOppositeFeatures.toArray(new Feature[origAndOppositeFeatures.size()]));
        List<Float> flaggedReasons = new ArrayList<Float>();

        for (Feature feature : origOppositeFeatureMap.keySet())
        {
            Feature oppFeature = origOppositeFeatureMap.get(feature);

            int flagReason = processFeature(feature, oppFeature,
                    featureIsLightMap.get(feature), run);
            flaggedReasons.add((float) flagReason);
            if (flagReason != FLAG_REASON_NONE)
            {
                feature.setDescription(flagReasonDescs[flagReason]);
                flaggedFeatures.add(feature);
            }
        }

        if (showCharts)
        {
            new PanelWithHistogram(flaggedReasons, "flag reason").displayInTab();
        }
        ms2FeatureSet.setFeatures(flaggedFeatures.toArray(new Feature[flaggedFeatures.size()]));
    }

    protected int processFeature(Feature feature,// List<Feature> ms1MatchedFeatures,
                                 Feature oppFeature,// List<Feature> oppMs1MatchedFeatures,
                                 boolean baseIsLight,
                                 MSRun run)
    {



//List<Integer> badScans = Arrays.asList(new Integer[] {3424, 3516,8436,2304,3638,7493,8978,5545,2829,4879});
//Level origLevel = _log.getLevel();
//if (badScans.contains(feature.getScan()))
//{
//_log.setLevel(Level.DEBUG);
//_log.debug("***** scan " + feature.getScan());
//}
        if (baseIsLight)
            _log.debug("Feature is light");
        else
            _log.debug("Feature is heavy");

        int reason = FLAG_REASON_NONE;
            //calcPeakIntensities' result starts with one peak below
            List<Float> peakIntensities = calcPeakIntensities(feature, run);
            List<Float> oppPeakIntensities = calcPeakIntensities(oppFeature, run);

            float intensityBelow = peakIntensities.remove(0);
            float intensityBelowOpp = oppPeakIntensities.remove(0);

        _log.debug("**" + peakIntensities.get(0) + ", " + peakIntensities.get(1) + ", " + peakIntensities.get(2) + ", " + peakIntensities.get(3));
        _log.debug("**" + oppPeakIntensities.get(0) + ", " + oppPeakIntensities.get(1) + ", " + oppPeakIntensities.get(2) + ", " + oppPeakIntensities.get(3));

        _log.debug("************************" + calcKL(feature.getMass(), peakIntensities) + ", " + calcKL(oppFeature.getMass(), oppPeakIntensities));

        //check the KL of both "features"
        float kl = calcKL(feature.getMass(), peakIntensities);
        float klOpp = calcKL(oppFeature.getMass(), oppPeakIntensities);
        float klDiff = Math.abs(kl - klOpp);
        float klRatio = kl / klOpp;
        if (kl > klOpp)
            klRatio = 1/klRatio;
        _log.debug("KL: " + kl + ", klOPP: " + klOpp + ": diff=" + klDiff + ", ratio=" + klRatio);
        if (klDiff > 0.25 && klRatio < 0.7)
        {
            reason = FLAG_REASON_DIFF_PEAKDIST;
            _log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
            return reason;
        }

        //See if the intensity of the peak 1Da below monoisotopic is high enough to worry about
        float maxPeakIntensity = BasicStatistics.max(peakIntensities);
        float maxPeakIntensityOpp = BasicStatistics.max(oppPeakIntensities);

        float belowIntensityRatio = intensityBelow / maxPeakIntensity;
        float belowIntensityRatioOpp = intensityBelowOpp / maxPeakIntensityOpp;

        _log.debug("BELOW_RATIO: " + belowIntensityRatio + ", OPP: " + belowIntensityRatioOpp);
        if (Math.max(belowIntensityRatio, belowIntensityRatioOpp) > 0.5)
        {
            reason = FLAG_REASON_COELUTING;
            _log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
            return reason;
        }
        _log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
        return reason;
    }

    /**
     * Calculate a KL score for a list of peak intensities
     * @param mass
     * @param peakIntensities
     * @return
     */
    protected float calcKL(float mass, List<Float> peakIntensities)
    {
        int signalLength = peakIntensities.size();
        float[] signal = new float[signalLength];
        float sum = 0;
        for (int i = 0; i < peakIntensities.size(); i++)
        {
            if (i < signal.length)
            {
                signal[i] = Math.max(0.1F, peakIntensities.get(i));
                sum += signal[i];
            }
        }
        for (int i = 0; i < signal.length; i++)
            signal[i] /= sum;

        return Spectrum.KLPoissonDistance(mass, signal);
    }


    /**
     * Return the list of peak intensities, starting with 1Da below mass.
     * First, resample spectra onto a high-res grid.
     * @param feature
     * @param run
     * @return
     */
    protected List<Float> calcPeakIntensities(Feature feature, MSRun run)
    {
        //Define the m/z range that I want to look at.  1Da down and 6Da up.
        //Todo: narrow this on the up side?

        //This better be bigger in all cases than peakPPMTolerance
        float spilloverMz = 0.2f;
        float lowMz = feature.getMz() -
                ((1 * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH + spilloverMz) /
                        feature.getCharge());
        float highMz = feature.getMz() +
                ((6 * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH + spilloverMz) /
                        feature.getCharge());
        FloatRange mzRange = new FloatRange(lowMz, highMz);

        //Define the scan range, making sure not to go out of bounds
        //assumes heavy and light scan extents same
        int firstScanIndex = Math.max(Math.abs(
                run.getIndexForScanNum(IsotopicLabelExtraInfoDef.getHeavyFirstScan(feature))) -
                                       numScansAroundEventToConsider, 0);
        int lastScanIndex = Math.min(Math.abs(
                run.getIndexForScanNum(IsotopicLabelExtraInfoDef.getHeavyLastScan(feature))) +
                                       numScansAroundEventToConsider, run.getScanCount()-1);
                                     
        Scan[] scans = FeatureFinder.getScans(run, firstScanIndex,
                lastScanIndex - firstScanIndex + 1);
        _log.debug("Scans: " + scans.length);

        SpectrumResampler spectrumResampler = new SpectrumResampler(mzRange);
        int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;
        spectrumResampler.setResamplingFrequency(resolution);
        spectrumResampler.setUseMedianSmooth(false);
        float[][] resampledSpectra = null;

        try
        {
            resampledSpectra =
                    spectrumResampler.resampleSpectra(scans);
        }
        catch (InterruptedException e)
        {
            throw new RuntimeException("Unexpectedly interrupted while resampling");
        }

        List<Float> peakIntensities = new ArrayList<Float>();
        for (int i=-1; i<4; i++)
        {
            float intensitySum = 0;
            //m/z of this peak
            float mz = feature.getMz() +
                    ((float) i * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH /
                            (float) feature.getCharge());
            //tolerance around each peak to grab.  Convert ppm to m/z, based on feature mass and charge
            float mzTol = (MassUtilities.calculateAbsoluteDeltaMass(
                    feature.getMass(), peakPPMTolerance, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / feature.getCharge();
            float minPeakMz = mz - mzTol;
            float maxPeakMz = mz + mzTol;

            //calculate the resampled spectra buckets that we need to add up
            float lowBucket = (minPeakMz - mzRange.min) * resolution;
            int lowBucketIndex = (int)Math.floor(lowBucket);
            float highBucket = (maxPeakMz - mzRange.min) * resolution;
            int highBucketIndex = (int)Math.floor(highBucket);
//System.err.println(i + ", " + mz + ", " + minPeakMz + ", " + maxPeakMz + ", lowb=" + lowBucket + ", highb=" + highBucket + ", buckets=" + resampledSpectra[0].length +  ", highmzrecalc=" + (mzRange.min + (highBucketIndex / (float) resolution)));

            for (float[] scanSpectra : resampledSpectra)
            {
                for (int bucket = lowBucketIndex; bucket <= highBucketIndex; bucket++)
                    intensitySum += scanSpectra[bucket];
            }

//_log.debug(i + ", " + mz + ", " + minPeakMz + ", " + maxPeakMz + ", lowb=" + lowBucket + ", highb=" + highBucket  + ", int=" + intensitySum);

            peakIntensities.add(intensitySum);
        }
        return peakIntensities;
    }

    /**
     * Returns false if the peptide does not contain the residue
     * todo: move
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    protected boolean checkForModAllResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return false;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                boolean foundIt = false;
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    return false;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                    {
                        foundIt = true;
                        break;
                    }
                if (!foundIt)
                    return false;
            }
        }
        return true;
    }

    /**
     * todo: move
     * @param feature
     * @return
     */
    protected float calcHeavyLightMassDiff(Feature feature)
    {
        String residue = "";
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                residue = "C";
                break;
            case LABEL_LYCINE:
                residue = "K";
                break;
        }
        int numLabels = 0;
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
        for (int i=0; i<peptide.length(); i++)
        {
            if (residue.equals("" + peptide.charAt(i)))
                numLabels++;
        }

        return numLabels * (labelType == LABEL_ACRYLAMIDE ? ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE :
                SILAC_LABEL_MASSDIFF_PERRESIDUE);
    }


    /**
     * Returns false if the peptide does not contain the residue
     * todo: move
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    protected boolean checkForModNoResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return true;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    continue;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                        return false;
            }
        }
        return true;
    }


    /**
     * todo: move
     * @param feature
     * @return
     */
    protected boolean isLightLabeled(Feature feature)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean lightLabeled = false;
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                lightLabeled =  checkForModAllResidues(peptide, mods, 'C', ACRYLAMIDE_LABEL_LIGHTMASS);
                break;
            case LABEL_LYCINE:
                lightLabeled =  checkForModNoResidues(peptide, mods, 'K', SILAC_LABEL_MASS);
                break;
        }
        return lightLabeled;
    }

    /**
     * todo: move
     * @param feature
     * @return
     */
    protected boolean isHeavyLabeled(Feature feature)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean heavyLabeled = false;
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                heavyLabeled = checkForModAllResidues(peptide, mods, 'C', ACRYLAMIDE_LABEL_HEAVYMASS);
                break;
            case LABEL_LYCINE:
                heavyLabeled = checkForModAllResidues(peptide, mods, 'K', SILAC_LABEL_MASS);
                break;
        }
        return heavyLabeled;
    }
}
