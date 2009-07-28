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
import org.apache.log4j.Level;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
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
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

import java.io.File;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;


/**
 * This uses a big HACK, adding a dummy search score to peptide features to store the flag reason description
 *
 * todo: fix KL analysis and single-peak ratio calculation for 3Da-separated acrylamide peptides
 */
public class FlagQuantEventsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
//scaffolding
//lane d
//List<Integer> badScans = Arrays.asList(new Integer[] {3424, 3516,8436,2304,3638,7493,8978,5545,2829,4879});
//lane e
//List<Integer> badScans = Arrays.asList(new Integer[] {    9603, 6106, 6177, 8008, 1062, 1130, 6923, 7176, 7917, 3139, 3141, 4457, 4673, 2694,});
//lane f
List<Integer> badScans = Arrays.asList(new Integer[] {        2431, 3652, 3680, 9026, 9650, 5552, 5558, 5621, 5586,});

//lane g
//List<Integer> badScans = Arrays.asList(new Integer[] { 1312, 1280, 1332, 3871, 9333, 10577, 5382, 1672, 4144, 6134, 6193, 5079,});
//lane h
//List<Integer> badScans = Arrays.asList(new Integer[] { 6812, 8080, 10360, 11287, 1861, 4710, 4717, 1291, 2933, 4106, });
//lane h falsepos scans
//List<Integer> badScans = Arrays.asList(new Integer[] {11002, 679, 7976, 9687, 1839, 9980, 10337, 7686, 12358, 4972, 4482,});

//lane i
//List<Integer> badScans = Arrays.asList(new Integer[] {     10105, 10282, 7822, 1758, 4587, 5548, 8247, 7551,});

    List<Float> flaggedReasons = new ArrayList<Float>();


List<Integer> flaggedScans = new ArrayList<Integer>();
List<Integer> badFlagged = new ArrayList<Integer>();

List<Float> reasonsTruePos = new ArrayList<Float>();
List<Float> reasonsFalsePos = new ArrayList<Float>();


    protected static Logger _log = Logger.getLogger(FlagQuantEventsCLM.class);

    protected File[] featureFiles;
    protected File outFile;
    protected File outDir;
    protected File mzXmlDir;

    public static final int FLAG_REASON_NONE = 0;
    public static final int FLAG_REASON_COELUTING = 1;
    public static final int FLAG_REASON_DISSIMILAR_KL = 2;
    public static final int FLAG_REASON_DISSIMILAR_MS1_RATIO = 3;

    public static final String REASON_DUMMY_SEARCH_SCORE_NAME = "dummy_flag_desc";


    protected float massSlopPPM = 25f;
    protected float massSlopDa = 0.05f;

    //Maximum allowed proportion of highest peak intensity that the peak 1Da below monoisotope can have
    public static final float DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF = 0.4f;
    protected float peakBelowIntensityRatioCutoff = DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF;

    public static final int DEFAULT_NUM_SCANS_AROUND_EVENT = 0;
    protected int numScansAroundEventToConsider = DEFAULT_NUM_SCANS_AROUND_EVENT;

    public static final float DEFAULT_PEAK_PPM_TOLERANCE = 50;
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
            };


    protected float maxKlDiff = 0.2f;
    protected float minKlRatio = 0.7f;

    //Maximum allowed difference between log ratio we calculate (simply, by most-intense peaks) and algorithm
    protected float maxLogRatioDiff = (float) (Math.log(1.5) - Math.log(1));


    protected int numPeaksForKLCalc = 3;

    //Ratios must be higher than minFlagRatio, OR lower than maxFlagRatio, to be flagged
    protected float minFlagRatio = 0f;
    protected float maxFlagRatio = 999f;

    //For calculating statistics on flagged vs non-flagged peptides
    Map<String, List<Float>> peptidePerFractionRatioMap = new HashMap<String, List<Float>>();
    Map<String, Integer> peptideNumFractionsFlaggedMap = new HashMap<String, Integer>();



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
                        new BooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                        new DecimalArgumentDefinition("minflagratio", false,
                                "Ratios must be higher than this, OR lower than maxflagratio, or both, to flag", minFlagRatio),
                        new DecimalArgumentDefinition("maxflagratio", false,
                                "Ratios must be lower than this, OR higher than maxflagratio, or both, to flag", maxFlagRatio),
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
        if (labelType == LABEL_ACRYLAMIDE)
            ApplicationContext.infoMessage("WARNING: peak distribution analysis and MS1 ratio comparison will " +
                    "be inappropriate for high-mass 3Da-separated light and heavy isotopes.");

        minFlagRatio = getFloatArgumentValue("minflagratio");
        maxFlagRatio = getFloatArgumentValue("maxflagratio");
        if (hasArgumentValue("minflagratio") || hasArgumentValue("maxflagratio"))
            ApplicationContext.infoMessage("NOTE: only ratios higher than " + minFlagRatio + " or lower than " +
                    maxFlagRatio + " (or both) will be flagged");


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
                    _log.debug("\tCombining individual fraction files... " +
                            outputFile.getAbsolutePath() + "...");
                    new PepXMLFeatureFileHandler().combinePepXmlFiles(tempFeatureFiles, outputFile);
                }
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save output file " + outputFile.getAbsolutePath(), e);
            }
            ApplicationContext.infoMessage("Done.");
System.err.println("Bad: " + badScans.size() + ", bad flagged: " + badFlagged.size() + ", sens=" + ((float)badFlagged.size() / (float)badScans.size()) + ", spec=" + ((float)badFlagged.size() / (float)flaggedFeaturesThisFile));
System.err.println("False negatives: ");
for (int scan : badScans) if (!badFlagged.contains(scan)) System.err.println(scan);
System.err.println("False positives: ");
for (int scan : flaggedScans) if (!badScans.contains(scan)) System.err.println(scan);


            if (showCharts)
            {
//                new PanelWithHistogram(reasonsTruePos, "TruePos reasons").displayInTab();
//                new PanelWithHistogram(reasonsFalsePos, "FalsePos reasons").displayInTab();
                   new PanelWithHistogram(flaggedReasons, "Flag reasons").displayInTab();
            }
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

        Map<Feature, Feature> origOppositeFeatureMap = new HashMap<Feature, Feature>();
        List<Feature> origAndOppositeFeatures = new ArrayList<Feature>();
        Map<Feature, Boolean> featureIsLightMap = new HashMap<Feature, Boolean>();

        Map<String, List<Float>> peptideRatiosMap = new HashMap<String, List<Float>>();

        for (Feature feature : ms2FeatureSet.getFeatures())
        {
            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                continue;
            float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);

            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

            List<Float> ratiosThisPeptide = peptideRatiosMap.get(peptide);
            if (ratiosThisPeptide == null)
            {
                ratiosThisPeptide = new ArrayList<Float>();
                peptideRatiosMap.put(peptide, ratiosThisPeptide);
            }
            ratiosThisPeptide.add(ratio);

            if (ratio < minFlagRatio && ratio > maxFlagRatio)
            {
                _log.debug("Skipping ratio " + ratio);
                continue;
            }
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

        for (Feature feature : origOppositeFeatureMap.keySet())
        {
            Feature oppFeature = origOppositeFeatureMap.get(feature);

            Pair<Integer, String> flagReasonAndDesc = processFeature(feature, oppFeature,
                    featureIsLightMap.get(feature), run);
            int flagReason = flagReasonAndDesc.first;
            String flagReasonDesc = flagReasonAndDesc.second;
            flaggedReasons.add((float) flagReason);
            if (flagReason != FLAG_REASON_NONE)
            {
                feature.setDescription(flagReasonDesc);
                //This is a HACK, adding a dummy search score to store the flag reason description
                MS2ExtraInfoDef.addSearchScore(feature, REASON_DUMMY_SEARCH_SCORE_NAME, flagReasonDesc);
                flaggedFeatures.add(feature);
            }
        }

        for (String peptide : peptideRatiosMap.keySet())
        {
            List<Float> allFracRatiosThisPeptide = peptidePerFractionRatioMap.get(peptide);
            if (allFracRatiosThisPeptide == null)
            {
                allFracRatiosThisPeptide = new ArrayList<Float>();
                peptidePerFractionRatioMap.put(peptide, allFracRatiosThisPeptide);
            }
            allFracRatiosThisPeptide.add((float)BasicStatistics.geometricMean(peptideRatiosMap.get(peptide)));
        }

        ms2FeatureSet.setFeatures(flaggedFeatures.toArray(new Feature[flaggedFeatures.size()]));
    }

    /**
     *
     * @param feature
     * @param oppFeature
     * @param baseIsLight
     * @param run
     * @return a pair containing the result code and a String that gives details
     */
    protected Pair<Integer, String> processFeature(Feature feature,// List<Feature> ms1MatchedFeatures,
                                 Feature oppFeature,// List<Feature> oppMs1MatchedFeatures,
                                 boolean baseIsLight,
                                 MSRun run)
    {





Level origLevel = _log.getLevel();
if (badScans.contains(feature.getScan()))
{
_log.setLevel(Level.DEBUG);
_log.debug("***** scan " + feature.getScan());
}
        if (baseIsLight)
            _log.debug("Feature is light");
        else
            _log.debug("Feature is heavy");

        int reason = FLAG_REASON_NONE;
        String reasonDesc = "";

        //calcPeakIntensities' result starts with one peak below
        List<Float> lightPeakIntensities = calcPeakIntensities(baseIsLight ? feature : oppFeature, run);
        List<Float> heavyPeakIntensities = calcPeakIntensities(baseIsLight ? oppFeature : feature, run);

        float intensityBelowLight = lightPeakIntensities.remove(0);
        float intensityBelowHeavy = heavyPeakIntensities.remove(0);

_log.debug("**light, " + lightPeakIntensities.get(0) + ", " + lightPeakIntensities.get(1) + ", " + lightPeakIntensities.get(2) + ", " + lightPeakIntensities.get(3));
 _log.debug("**heavy, " + heavyPeakIntensities.get(0) + ", " + heavyPeakIntensities.get(1) + ", " + heavyPeakIntensities.get(2) + ", " + heavyPeakIntensities.get(3));


        if (reason == FLAG_REASON_NONE)
        {
            //See if the intensity of the peak 1Da below monoisotopic is high enough to worry about
            float belowIntensityRatioLight = intensityBelowLight / lightPeakIntensities.get(0);
            float belowIntensityRatioHeavy = intensityBelowHeavy / heavyPeakIntensities.get(0);

            _log.debug("BELOW: light=" + intensityBelowLight + ", ratio="+belowIntensityRatioLight +
                    ", heavy=" +intensityBelowHeavy + ", ratio=" +  belowIntensityRatioHeavy);
            if (Math.max(belowIntensityRatioLight, belowIntensityRatioHeavy) > peakBelowIntensityRatioCutoff)
            {
                reason = FLAG_REASON_COELUTING;
                reasonDesc = "COELUTE.  intensity ratio light=" + Rounder.round(belowIntensityRatioLight,3) +
                        ", heavy=" + Rounder.round(belowIntensityRatioHeavy,3);
            }
        }


        if (reason == FLAG_REASON_NONE)
        {
            //Compare intensities of the peak that has theoretical max intensity
            int theoreticalMaxIndex = Spectrum.calcMaxIdealPeakIndex(feature.getMass());
//            Spectrum.KLPoissonDistance()
//            int maxPeakIndexLight = 0;
//            int maxPeakIndexHeavy = 0;
//            float maxPeakIntensityLight = 0;
//            float maxPeakIntensityHeavy = 0;
//            for (int i=0; i<lightPeakIntensities.size(); i++)
//            {
//                if (lightPeakIntensities.get(i) > maxPeakIntensityLight)
//                {
//                    maxPeakIntensityLight = lightPeakIntensities.get(i);
//                    maxPeakIndexLight = i;
//                }
//                if (heavyPeakIntensities.get(i) > maxPeakIntensityHeavy)
//                {
//                    maxPeakIntensityHeavy = heavyPeakIntensities.get(i);
//                    maxPeakIndexHeavy = i;
//                }
//            }
//
//            //simple ratio based on one peak: the peak that has the highest single intensity among all peaks in both features
//            int maxPeakIndexHeavyLight = maxPeakIndexLight;
//            if (maxPeakIntensityHeavy > maxPeakIntensityLight)
//                maxPeakIndexHeavyLight = maxPeakIndexHeavy;
//            float ms1Ratio = peakIntensities.get(maxPeakIndexHeavyLight) / oppPeakIntensities.get(maxPeakIndexHeavyLight);
            float ms1Ratio = lightPeakIntensities.get(theoreticalMaxIndex) / heavyPeakIntensities.get(theoreticalMaxIndex);

            float algRatio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
            float logRatioDiff = (float) Math.abs(Math.log(ms1Ratio) - Math.log(algRatio));
            _log.debug("MS1 Ratio: " + ms1Ratio + ", alg ratio: " + algRatio + ", log diff: " + logRatioDiff + ", to beat: " + maxLogRatioDiff);
            if (logRatioDiff > maxLogRatioDiff)
            {
                reason = FLAG_REASON_DISSIMILAR_MS1_RATIO;
                reasonDesc = "DIFF SINGLEPEAK RATIO. singlepeak=" + Rounder.round(ms1Ratio, 3) + ", algorithm=" +
                        Rounder.round(algRatio, 3) + ", log diff: " + Rounder.round(logRatioDiff,3);
            }
        }


        if (reason == FLAG_REASON_NONE)
        {
            //check the KL of both "features"

            float lightKl = calcKL(feature.getMass(), lightPeakIntensities);
            float heavyKl = calcKL(oppFeature.getMass(), heavyPeakIntensities);
            float klDiff = Math.abs(lightKl - heavyKl);
            float klRatio = lightKl / heavyKl;
            _log.debug("Light KL: " + lightKl + ", Heavy KL: " + heavyKl + ": diff=" + klDiff + ", ratio=" + klRatio);
            if (klDiff > maxKlDiff && klRatio < minKlRatio)
            {
                reason = FLAG_REASON_DISSIMILAR_KL;
                reasonDesc = "DIFF KL. light=" + Rounder.round(lightKl,3) + ", heavy=" + Rounder.round(heavyKl,3) +
                        ", diff=" + Rounder.round(klDiff,3) + ", ratio=" + Rounder.round(klRatio,3);
            }
        }


        _log.debug("REASON: " + flagReasonDescs[reason]);
_log.setLevel(origLevel);
if (reason != FLAG_REASON_NONE)
{
    flaggedScans.add(feature.getScan());
    if (badScans.contains(feature.getScan()))
    {
        badFlagged.add(feature.getScan());
        reasonsTruePos.add((float)reason);
    }
    else
        reasonsFalsePos.add((float)reason);
}

        return new Pair<Integer, String>(reason, reasonDesc);
    }

    /**
     * Calculate a KL score for a list of peak intensities
     * @param mass
     * @param peakIntensities
     * @return
     */
    protected float calcKL(float mass, List<Float> peakIntensities)
    {
        //use only numPeaksForKLCalc peaks
        int numPeaks = Math.min(numPeaksForKLCalc, peakIntensities.size());
        float[] signal = new float[numPeaks];
        float sum = 0f;
        for (int i=0; i<numPeaks; i++)
        {
            signal[i] = peakIntensities.get(i);
            sum += signal[i];
        }
_log.debug("Sum: " + sum);

        for (int i = 0; i < signal.length; i++)
        {
            signal[i] /= sum;
_log.debug(signal[i]);                        
        }

        //Sometimes KL comes back negative
        return Math.max(Spectrum.KLPoissonDistance(mass, signal), 0f);
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

        //This better be bigger in all cases than peakPPMTolerance.  Could calc here.
        float spilloverMz = 0.2f;
        float lowMz = feature.getMz() -
                ((1 * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH + spilloverMz) /
                        feature.getCharge());
        //Todo: narrow this on the up side?
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
        lastScanIndex = Math.max(firstScanIndex, lastScanIndex);
                                     
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
            int lowBucketIndex = (int)Math.max(0, Math.floor(lowBucket));
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
