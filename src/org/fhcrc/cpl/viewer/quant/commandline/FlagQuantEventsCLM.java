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
import org.fhcrc.cpl.viewer.quant.PeakOverlapCorrection;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
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
import java.util.List;


/**
 * This uses a big HACK, adding a dummy search score to peptide features to store the flag reason description
 *
 * todo: fix KL analysis for 3Da-separated acrylamide peptides
 */
public class FlagQuantEventsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{

    //Controls whether the badScans list is used to calculate sensitivity and specificity.
    //This is for refining flagquant
    protected boolean shouldCalcSensSpecWithHardcodedScans = false;

//    List<Integer> badScans = Arrays.asList(new Integer[] {});
//scaffolding
//lane d
//List<Integer> badScans = Arrays.asList(new Integer[] {3424, 3516,8436,2304,3638,7493,8978,5545,2829,4879});
//lane e
List<Integer> badScans = Arrays.asList(new Integer[] {    9603, 6106, 6177, 8008, 1062, 1130, 6923, 7176, 7917, 3139, 3141, 4457, 4673, 2694,});
//lane f
//List<Integer> badScans = Arrays.asList(new Integer[] {        2431, 3652, 3680, 9026, 9650, 5552, 5558, 5621, 5586,});

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

    List<Float> algRatios = new ArrayList<Float>();
    List<Float> weightedGeomMeanRatios = new ArrayList<Float>();



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

    //Maximum allowed proportion of highest peak intensity that the peak 1Da below monoisotope can have
    public static final float DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF = 0.4f;
    protected float peakBelowIntensityRatioCutoff = DEFAULT_PEAKBELOW_INTENSITY_RATIO_CUTOFF;

    //todo: set this back to 0 when done playing with scan-by-scan ratio stuff
    public static final int DEFAULT_NUM_SCANS_AROUND_EVENT = 15;
    protected int numScansAroundEventToConsider = DEFAULT_NUM_SCANS_AROUND_EVENT;

    public static final float DEFAULT_PEAK_PPM_TOLERANCE = 50;
    protected float peakPPMTolerance = DEFAULT_PEAK_PPM_TOLERANCE;


    protected int labelType = QuantitationUtilities.LABEL_LYCINE;

    protected boolean showCharts = false;


    public static final float SILAC_LABEL_MASS = 134.115092f;
    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;


    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;

    //scaffolding
    protected List<Float> lightPrevPeakRatios = new ArrayList<Float>();
    protected List<Float> heavyPrevPeakRatios = new ArrayList<Float>();


    protected int numPeaksToUse = 5;
    protected int numRawPeaksToKeep = 2;

    //minimum proportion of total light feature intensity that the peak /below/ the heavy monoisotope can have,
    //in order for us /not/ to check for a peak below the heavy monoisotope
    protected float minSignificantPeakContributionBelowMonoisotope = 0.02f;

    protected float maxLightHeavyOverlapToIgnore = 0.02f;



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
    Map<String, List<Pair<Float, Boolean>>> peptidePerFractionRatioFlaggedMap = new HashMap<String, List<Pair<Float,Boolean>>>();
    Map<String, Integer> peptideNumFractionsFlaggedMap = new HashMap<String, Integer>();
    List<Float> stDevsFromMeanRatioNoFlag = new ArrayList<Float>();
    List<Float> stDevsFromMeanRatioFlag = new ArrayList<Float>();
    List<Float> corrFlagged = new ArrayList<Float>();
    List<Float> corrUnFlagged = new ArrayList<Float>();
    //for sens, spec
    List<Float> corrBad = new ArrayList<Float>();
    List<Float> corrGood = new ArrayList<Float>();

    List<Float> singlePeakRatios = new ArrayList<Float>();
    List<Float> singlePeakRatiosNoCorrect = new ArrayList<Float>();



    //scaffolding for calculating ratios using regression based on one datapoint per-scan, like RelEx.
    //I think that method pretty much doesn't work very well.
//    List<Float> singlePeakSlopeRatios = new ArrayList<Float>();
//    List<Float> multiPeakSlopeRatios = new ArrayList<Float>();

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
                        new FileToWriteArgumentDefinition("out", false, "Output pepXML file"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
                        new EnumeratedValuesArgumentDefinition("label", false, QuantitationUtilities.ALL_LABEL_CODES,
                               QuantitationUtilities.ALL_LABEL_EXPLANATIONS, QuantitationUtilities.LABEL_LYCINE_CODE),
                        new BooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                        new DecimalArgumentDefinition("minflagratio", false,
                                "Ratios must be higher than this, or lower than maxflagratio, or both, to flag",
                                minFlagRatio),
                        new DecimalArgumentDefinition("maxflagratio", false,
                                "Ratios must be lower than this, or higher than maxflagratio, or both, to flag",
                                maxFlagRatio),
                        new DecimalArgumentDefinition("peakppm", false,
                                "Mass tolerance around each theoretical peak (ppm)", DEFAULT_PEAK_PPM_TOLERANCE),
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

        if (outFile == null)
            assertArgumentPresent("outdir");
        if (outDir == null)
            assertArgumentPresent("out");        

        mzXmlDir = getFileArgumentValue("mzxmldir");

        labelType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label"));
        if (labelType == QuantitationUtilities.LABEL_ACRYLAMIDE)
            ApplicationContext.infoMessage("WARNING: this tool does not currently know how to account for overlap of " +
                    "light and heavy peaks with 3Da separation.  Peak distribution analysis and MS1 ratio comparison " +
                    "will not be performed for 3Da-separated light and heavy isotopes.");

        minFlagRatio = getFloatArgumentValue("minflagratio");
        maxFlagRatio = getFloatArgumentValue("maxflagratio");
        if (hasArgumentValue("minflagratio") || hasArgumentValue("maxflagratio"))
            ApplicationContext.infoMessage("NOTE: only ratios higher than " + minFlagRatio + " or lower than " +
                    maxFlagRatio + " (or both) will be flagged");

        peakPPMTolerance = getFloatArgumentValue("peakppm");
        showCharts = getBooleanArgumentValue("showcharts");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File featureFile : featureFiles)
        {
            ApplicationContext.infoMessage("Handling file " + featureFile.getAbsolutePath());
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
                    baseName = baseName.replaceFirst(".." + File.separator, "");

                int numFeaturesThisFraction = featureSet.getFeatures().length;
                File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);
                processFeatureSet(featureSet, run);

                int numFlaggedFeaturesThisFraction = featureSet.getFeatures().length;

                ApplicationContext.infoMessage("Flagged " + numFlaggedFeaturesThisFraction + " out of " +
                        numFeaturesThisFraction + " features this fraction");
                flaggedFeaturesThisFile += numFlaggedFeaturesThisFraction;

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
            if (shouldCalcSensSpecWithHardcodedScans && !badScans.isEmpty())
            {
                System.err.println("Bad: " + badScans.size() + ", bad flagged: " + badFlagged.size() + ", sens=" + ((float)badFlagged.size() / (float)badScans.size()) + ", spec=" + (((float) flaggedFeaturesThisFile - (float)badFlagged.size()) / (float) flaggedFeaturesThisFile));
                System.err.println("False negatives: ");
                for (int scan : badScans) if (!badFlagged.contains(scan)) System.err.println(scan);
                System.err.println("False positives: ");
                for (int scan : flaggedScans) if (!badScans.contains(scan)) System.err.println(scan);
            }


            if (showCharts)
            {
//                new PanelWithHistogram(reasonsTruePos, "TruePos reasons").displayInTab();
//                new PanelWithHistogram(reasonsFalsePos, "FalsePos reasons").displayInTab();
                   new PanelWithHistogram(flaggedReasons, "Flag reasons").displayInTab();


                List<Float> cvs = new ArrayList<Float>();
                List<Float> numFlagged = new ArrayList<Float>();
                List<Float> distFromMeanFlagged = new ArrayList<Float>();
                List<Float> distFromMeanNotFlagged = new ArrayList<Float>();

                for (String peptide : peptidePerFractionRatioFlaggedMap.keySet())
                {
                    List<Pair<Float, Boolean>> ratioFlaggeds = peptidePerFractionRatioFlaggedMap.get(peptide);
                    List<Float> logRatios = new ArrayList<Float>();

                    for (Pair<Float, Boolean> pair : ratioFlaggeds)
                    {
                        logRatios.add((float)Math.log(pair.first));

                    }
                    cvs.add((float)BasicStatistics.coefficientOfVariation(logRatios));
                    int numFracsFlagged = 0;
                    if (peptideNumFractionsFlaggedMap.containsKey(peptide))
                        numFracsFlagged = peptideNumFractionsFlaggedMap.get(peptide);
                    numFlagged.add((float) numFracsFlagged);

                    if (logRatios.size() > 2)
                    {
                        float meanLog = (float) BasicStatistics.mean(logRatios);
                        float stddevLog = (float) BasicStatistics.standardDeviation(logRatios);
                        if (stddevLog > 0)
                        {
                            for (Pair<Float, Boolean> pair : ratioFlaggeds)
                            {
                                float scaled = ((float) Math.log(pair.first) - meanLog) / stddevLog;

                                if (pair.second)
                                {
                                    distFromMeanFlagged.add(Math.abs(scaled));
//     System.err.println("**" +Math.log(pair.first) + ", " + meanLog + ", " +  scaled);

                                }
                                else
                                {
                                    distFromMeanNotFlagged.add(Math.abs(scaled));
//       System.err.println(Math.log(pair.first) + ", " + meanLog + ", " + scaled);
                                }
                            }
                        }
                    }

                }
//                new PanelWithHistogram(distFromMeanFlagged, "dist flagged").displayInTab();
//                new PanelWithHistogram(distFromMeanNotFlagged, "dist notflagged").displayInTab();
                ApplicationContext.infoMessage("Mean dist from mean, flagged: " + BasicStatistics.mean(distFromMeanFlagged) + ", not: " + BasicStatistics.mean(distFromMeanNotFlagged));

                if (showCharts)
                {
                    new PanelWithScatterPlot(cvs, numFlagged, "CV vs num fractions flagged").displayInTab();
                    
//                    System.err.println("Mean cor flagged: " + BasicStatistics.mean(corrFlagged) + ", unflagged: " + BasicStatistics.mean(corrUnFlagged));
//                    new PanelWithHistogram(corrFlagged, "cor flagged").displayInTab();
//                    new PanelWithHistogram(corrUnFlagged, "cor unflagged").displayInTab();
                    if (shouldCalcSensSpecWithHardcodedScans && !badScans.isEmpty())
                    {
                        System.err.println("Mean cor bad: " + BasicStatistics.mean(corrBad) + ", good: " + BasicStatistics.mean(corrGood));

                        new PanelWithHistogram(corrBad, "cor bad").displayInTab();
                        new PanelWithHistogram(corrGood, "cor good").displayInTab();
                    }

                    if (!lightPrevPeakRatios.isEmpty())
                    {
                        new PanelWithHistogram(lightPrevPeakRatios, "light belowpeak ratios").displayInTab();
                        new PanelWithHistogram(heavyPrevPeakRatios, "heavy belowpeak ratios").displayInTab();
                        System.err.println("Light belowpeak median: " + BasicStatistics.median(lightPrevPeakRatios) + ", heavy: " + BasicStatistics.median(heavyPrevPeakRatios));
                        new PanelWithScatterPlot(lightPrevPeakRatios,heavyPrevPeakRatios, "light vs heavy belowpeak ratio").displayInTab();
                    }

//                    new PanelWithScatterPlot(algorithmRatios, multiPeakSlopeRatios, "alg vs multislope").displayInTab();
                    new PanelWithScatterPlot(algRatios, singlePeakRatios, "alg vs singlepeak").displayInTab();
                    new PanelWithScatterPlot(singlePeakRatiosNoCorrect, singlePeakRatios, "singlepeak orig vs correct").displayInTab();

//                    new PanelWithScatterPlot(singlePeakRatios, multiPeakSlopeRatios, "singlepeak vs multislope").displayInTab();
//                    new PanelWithScatterPlot(singlePeakSlopeRatios, multiPeakSlopeRatios, "slope vs multislope").displayInTab();
//                    new PanelWithScatterPlot(algRatios, weightedGeomMeanRatios, "alg vs weightedmean").displayInTab();
                }
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

        Map<String, Map<Integer, List<Float>>> peptideChargeRatiosMap = new HashMap<String, Map<Integer, List<Float>>>();
        Set<String> flaggedPeptides = new HashSet<String>();

        for (Feature feature : ms2FeatureSet.getFeatures())
        {

            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                continue;
            float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);

            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

            Map<Integer, List<Float>> chargeRatiosMap = peptideChargeRatiosMap.get(peptide);
            if (chargeRatiosMap == null)
            {
                chargeRatiosMap = new HashMap<Integer, List<Float>>();
                peptideChargeRatiosMap.put(peptide, chargeRatiosMap);
            }
            List<Float> ratiosThisCharge = chargeRatiosMap.get(feature.charge);
            if (ratiosThisCharge == null)
            {
                ratiosThisCharge = new ArrayList<Float>();
                chargeRatiosMap.put(feature.charge, ratiosThisCharge);
            }
            ratiosThisCharge.add(ratio);

            if (ratio < minFlagRatio && ratio > maxFlagRatio)
            {
                _log.debug("Skipping ratio " + ratio);
                continue;
            }
            Feature oppFeature = (Feature) feature.clone();
            boolean isLight = false;
            if (IsotopicLabelExtraInfoDef.isHeavyLabeled(feature, labelType))
                oppFeature.setMass(feature.getMass() - QuantitationUtilities.calcHeavyLightMassDiff(peptide, labelType));
            else if (IsotopicLabelExtraInfoDef.isLightLabeled(feature, labelType))
            {
                isLight = true;
                oppFeature.setMass(feature.getMass() + QuantitationUtilities.calcHeavyLightMassDiff(peptide, labelType));
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
        Set<String> peptidesFlaggedThisSet = new HashSet<String>();

        //sort features by scan for performance (loading scans into memory)
        List<Feature> featuresByScan = new ArrayList<Feature>(origOppositeFeatureMap.keySet());
        Collections.sort(featuresByScan, new Feature.ScanAscComparator());

        int i=0;
        int numFeatures = featuresByScan.size();
        ApplicationContext.infoMessage("Processing " + numFeatures + " features....");
        for (Feature feature : featuresByScan)
        {
            Feature oppFeature = origOppositeFeatureMap.get(feature);
            if (numFeatures > 10 && i % (int) (numFeatures/10f) == 0)
                ApplicationContext.infoMessage((10 * (i / (int) (numFeatures/10f))) + "%");
            i++;
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
                flaggedPeptides.add(MS2ExtraInfoDef.getFirstPeptide(feature));
                peptidesFlaggedThisSet.add(MS2ExtraInfoDef.getFirstPeptide(feature));
            }
        }

        for (String peptide : flaggedPeptides)
        {
            int oldNum = 0;
            if (peptideNumFractionsFlaggedMap.containsKey(peptide))
                oldNum = peptideNumFractionsFlaggedMap.get(peptide);
            peptideNumFractionsFlaggedMap.put(peptide, oldNum+1);
        }

        for (String peptide : peptideChargeRatiosMap.keySet())
        {
            Map<Integer, List<Float>> chargeRatiosMap = peptideChargeRatiosMap.get(peptide);
            for (int charge :chargeRatiosMap.keySet() )
            {
            List<Pair<Float, Boolean>> allFracRatiosFlaggedThisPeptide = peptidePerFractionRatioFlaggedMap.get(peptide);
            if (allFracRatiosFlaggedThisPeptide == null)
            {
                allFracRatiosFlaggedThisPeptide = new ArrayList<Pair<Float, Boolean>>();
                peptidePerFractionRatioFlaggedMap.put(peptide, allFracRatiosFlaggedThisPeptide);
            }
            allFracRatiosFlaggedThisPeptide.add(new Pair<Float, Boolean>(
                    (float)BasicStatistics.geometricMean(chargeRatiosMap.get(charge)),
                    peptidesFlaggedThisSet.contains(peptide)));
            }
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
        //Retain original debug level for setting back later.  Override to DEBUG level here if we're tracking sens/spec
        Level origLevel = _log.getLevel();
        if (shouldCalcSensSpecWithHardcodedScans && badScans.contains(feature.getScan()))
        {
            _log.setLevel(Level.DEBUG);
            _log.debug("***** bad scan " + feature.getScan());
        }
        if (baseIsLight)
            _log.debug("Feature is light");
        else
            _log.debug("Feature is heavy");

        int reason = FLAG_REASON_NONE;
        String reasonDesc = "";

        QuantPeakSetSummary lightPeaksSummary = calcPeakIntensities(baseIsLight ? feature : oppFeature, run,
                numPeaksToUse);
        QuantPeakSetSummary heavyPeaksSummary = calcPeakIntensities(baseIsLight ? oppFeature : feature, run,
                numPeaksToUse);

        float intensityBelowLight = lightPeaksSummary.sumIntensityPeakBelow;
        float intensityBelowHeavy = heavyPeaksSummary.sumIntensityPeakBelow;

        List<Float> lightPeakIntensities = lightPeaksSummary.peakSumIntensities;
        List<Float> heavyPeakIntensities = heavyPeaksSummary.peakSumIntensities;

        _log.debug("**light, " + lightPeakIntensities.get(0) + ", " + lightPeakIntensities.get(1) + ", " +
                lightPeakIntensities.get(2) + ", " + lightPeakIntensities.get(3));
        _log.debug("**heavy, " + heavyPeakIntensities.get(0) + ", " + heavyPeakIntensities.get(1) + ", " +
                heavyPeakIntensities.get(2) + ", " + heavyPeakIntensities.get(3));

        int numPeaksSeparation = PeakOverlapCorrection.calcNumPeaksSeparation(lightPeaksSummary.monoisotopicMass,
                heavyPeaksSummary.monoisotopicMass);
        double[] lightIsotopicDistribution = PeakOverlapCorrection.getIsotopicDistribution(
                lightPeaksSummary.monoisotopicMass, numPeaksSeparation+1);
        boolean lightHasSignificantPeakBelowHeavy = lightIsotopicDistribution[numPeaksSeparation-1] >
                minSignificantPeakContributionBelowMonoisotope;
        boolean lightHasSignificantHeavyOverlap = lightIsotopicDistribution[numPeaksSeparation] >
                maxLightHeavyOverlapToIgnore;

        //See if the intensity of the peak 1Da below either monoisotope is high enough to worry about.
        if (reason == FLAG_REASON_NONE)
        {
            float belowIntensityRatioLight = (float)Math.log(intensityBelowLight / lightPeakIntensities.get(0));
            //Only evaluate heavy if if light/heavy peaks not overlapping                        
            float belowIntensityRatioHeavy = lightHasSignificantPeakBelowHeavy ? 0 :
                    (float) Math.log(intensityBelowHeavy / heavyPeakIntensities.get(0));

            //record-keeping.  Don't include pairs where there's significant overlap
            if (!lightHasSignificantPeakBelowHeavy &&
                !Float.isNaN(belowIntensityRatioLight) && !Float.isInfinite(belowIntensityRatioLight) &&
                !Float.isNaN(belowIntensityRatioHeavy) && !Float.isInfinite(belowIntensityRatioHeavy))
            {
                lightPrevPeakRatios.add(belowIntensityRatioLight);
                heavyPrevPeakRatios.add(belowIntensityRatioHeavy);
            }

            _log.debug("BELOW: light=" + intensityBelowLight + ", ratio="+belowIntensityRatioLight +
                    ", heavy=" +intensityBelowHeavy + ", ratio=" +  belowIntensityRatioHeavy);
            if (Math.max(belowIntensityRatioLight, belowIntensityRatioHeavy) > peakBelowIntensityRatioCutoff)
            {
                reason = FLAG_REASON_COELUTING;
                reasonDesc = "COELUTE.  intensity ratio light=" + Rounder.round(belowIntensityRatioLight,3) +
                        ", heavy=" + Rounder.round(belowIntensityRatioHeavy,3);
            }
        }

        float algRatio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);

        //Calculate a ratio from the theoretically most-intense peak, compare with algorithm ratio
        if (reason == FLAG_REASON_NONE)
        {
            int highestPeakIndex = Spectrum.calcMaxIdealPeakIndex(lightPeaksSummary.monoisotopicMass);

            float lightAreaHighestPeak = lightPeakIntensities.get(highestPeakIndex);
            float heavyAreaHighestPeak = heavyPeakIntensities.get(highestPeakIndex);

            float singlePeakRatio = lightAreaHighestPeak / heavyAreaHighestPeak;
            //store ratio prior to correction
            float singlePeakRatioNoCorrect = singlePeakRatio;

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
                    logRatioDiff + ", to beat: " + maxLogRatioDiff);
if (!Float.isInfinite((float)Math.log(singlePeakRatio)) && !Float.isInfinite((float)Math.log(algRatio)) &&
        !Float.isInfinite((float) Math.log(singlePeakRatioNoCorrect)) && !Float.isNaN((float) Math.log(singlePeakRatioNoCorrect)))
{
singlePeakRatiosNoCorrect.add((float)Math.log(singlePeakRatioNoCorrect));
singlePeakRatios.add((float)Math.log(singlePeakRatio));
algRatios.add((float)Math.log(algRatio));
}


            if (logRatioDiff > maxLogRatioDiff)
            {
                reason = FLAG_REASON_DISSIMILAR_MS1_RATIO;
                reasonDesc = "DIFF SINGLEPEAK RATIO. single=" + Rounder.round(singlePeakRatio, 3) + ", algorithm=" +
                        Rounder.round(algRatio, 3) + ", log diff: " + Rounder.round(logRatioDiff,3);
            }
        }

        //check the KL of both "features".  If peaks overlapping, calculate a special ideal distribution for heavy
        if (reason == FLAG_REASON_NONE)
        {
            float lightKl = calcKL(feature.getMass(), lightPeakIntensities);
            float heavyKl = 0f;
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
                heavyKl = calcKL(heavyIdealDist, heavyPeakIntensities);
            }
            else
                heavyKl = calcKL(heavyPeaksSummary.monoisotopicMass, heavyPeakIntensities);
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

            if (shouldCalcSensSpecWithHardcodedScans)
            {
                if (badScans.contains(feature.getScan()))
                {
                    badFlagged.add(feature.getScan());
                    reasonsTruePos.add((float)reason);
                }
                else
                    reasonsFalsePos.add((float)reason);
            }
        }

//scaffolding for saving one chart for every single event to a file       
//        if (showCharts)
//        {
//            float cor = (float) BasicStatistics.correlationCoefficient(lightPeaksSummary.getRawIntensitiesHighestPeak(), heavyPeaksSummary.getRawIntensitiesHighestPeak());
//            if (reason == FLAG_REASON_NONE)
//                corrUnFlagged.add(cor);
//            else         corrFlagged.add(cor);
//            if (badScans.contains(feature.getScan()))
//                corrBad.add(cor);
//            else corrGood.add(cor);
//            File outDir = new File("/home/dhmay/temp/flaggedplots");
//            if (reason == FLAG_REASON_NONE)
//                outDir = new File("/home/dhmay/temp/unflaggedplots");
//            File corrPlotFile = new File(outDir, feature.scan + "_" + MS2ExtraInfoDef.getFirstPeptide(feature) + ".png");
//            try
//            {
//                PanelWithScatterPlot pwsp = new PanelWithScatterPlot(lightPeaksSummary.rawIntensitiesHighestPeaks.get(0),
//                        heavyPeaksSummary.rawIntensitiesHighestPeaks.get(0), "corr");
//                pwsp.addData(lightPeaksSummary.rawIntensitiesHighestPeaks.get(1),
//                        heavyPeaksSummary.rawIntensitiesHighestPeaks.get(1), "");
//                pwsp.addRegressionLine(0, true);
//                pwsp.addRegressionLine(1, true);
//                pwsp.addLine(algRatio, 0, 0, BasicStatistics.max(heavyPeaksSummary.rawIntensitiesHighestPeaks.get(0)));
//                pwsp.setSeriesColor(4, Color.darkGray);
//
//                double[] logRatios = new double[lightPeaksSummary.rawIntensitiesHighestPeaks.get(0).length];
//                double[] sums = new double[lightPeaksSummary.rawIntensitiesHighestPeaks.get(0).length];
//                double[] ratios2 = new double[lightPeaksSummary.rawIntensitiesHighestPeaks.get(0).length];
//                double[] sums2 = new double[lightPeaksSummary.rawIntensitiesHighestPeaks.get(0).length];
//                float maxSum = 0;
//                List<Float> ratiosToMean = new ArrayList<Float>();
//                List<Float> sumsForWeights = new ArrayList<Float>();
//
//                for (int i=0; i<lightPeaksSummary.rawIntensitiesHighestPeaks.get(0).length; i++)
//                {
//                    double logRatio = Math.log(lightPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i] / heavyPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i]);
//                    if (Double.isInfinite(logRatio) || Double.isNaN(logRatio) || logRatio>7)
//                        logRatio = 5;
//                    else
//                    {
//                        ratiosToMean.add((float)logRatio);
//                        sumsForWeights.add((float) (lightPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i] + heavyPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i]));
//                    }
//                    logRatios[i] =logRatio;
//                    double sum = (lightPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i] + heavyPeaksSummary.rawIntensitiesHighestPeaks.get(0)[i]);
//                    if (Double.isInfinite(sum) || Double.isNaN(sum))
//                        sum = 0;
//                    sums[i] = sum;
//
//                    double ratio2 = Math.log(lightPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i] / heavyPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i]);
//                    if (Double.isInfinite(ratio2) || Double.isNaN(ratio2) || ratio2>7)
//                        ratio2 = 5;
//                    else
//                    {
//                        ratiosToMean.add((float)ratio2);
//                        sumsForWeights.add((float) (lightPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i] + heavyPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i]));
//                    }
//                    ratios2[i] =ratio2;
//                    double sum2 = (lightPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i] + heavyPeaksSummary.rawIntensitiesHighestPeaks.get(1)[i]);
//                    if (Double.isInfinite(sum2) || Double.isNaN(sum2))
//                        sum2 = 0;
//                    sums2[i] = sum2;
//                    maxSum = (float) Math.max(maxSum, Math.max(sum, sum2));
//                }
//                PanelWithScatterPlot pwsp2 = new PanelWithScatterPlot(sums, logRatios, "ratio_v_sum");
//                pwsp2.addData(sums2, ratios2, "peak2");
//                pwsp2.addLine(0, Math.log(algRatio), 0, maxSum);
//                pwsp2.setSeriesColor(0, Color.blue); pwsp2.setSeriesColor(1, Color.green);  pwsp2.setSeriesColor(2, Color.red);
//                if (!ratiosToMean.isEmpty())
//                {
//                    double weightedMean = BasicStatistics.weightedMean(ratiosToMean, sumsForWeights);
//                    weightedGeomMeanRatios.add((float) Math.log(weightedMean));
//                    algRatios.add((float)Math.log(algRatio));
//                    pwsp2.addLine(0, weightedMean, 0, maxSum);
//                    pwsp2.setSeriesColor(3, Color.darkGray);
//                }
//
//                pwsp2.saveChartToImageFile(new File(outDir, feature.scan + "_" + MS2ExtraInfoDef.getFirstPeptide(feature) + ".ratiosum.png"));
//                pwsp.saveChartToImageFile(corrPlotFile);
//
//            }
//            catch(IOException e)
//
//            {
//                throw new RuntimeException(e);
//            }
//        }


        return new Pair<Integer, String>(reason, reasonDesc);
    }

    /**
     * Calculate a KL score for a list of peak intensities.  Have to normalize intensity first by dividing by peak sum
     * @param mass
     * @param peakIntensities
     * @return
     */
    protected float calcKL(float mass, List<Float> peakIntensities)
    {
        return calcKL(Spectrum.Poisson(mass), peakIntensities);
    }

    protected float calcKL(float[] idealPeaks, List<Float> peakIntensities)
    {
        float[] peakIntensities6Peaks = new float[6];
        float sum = 0;
        for (int i=0; i<peakIntensities6Peaks.length; i++)
        {
            if (i < peakIntensities.size())
                peakIntensities6Peaks[i] = peakIntensities.get(i);
            peakIntensities6Peaks[i] = Math.max(0.1f, peakIntensities6Peaks[i]);
            sum += peakIntensities6Peaks[i];
        }
		for (int i = 0; i < peakIntensities6Peaks.length; i++)
        {
            peakIntensities6Peaks[i] /= sum;
//System.err.println(peakIntensities6Peaks[i]);
        }
        return PeakOverlapCorrection.calcKLUsingTemplate(idealPeaks, peakIntensities6Peaks);
    }


    /**
     * Return the peak quant summary, with intensities of all peaks.
     * Use raw intensities, NOT resampled
     * @param feature
     * @param run
     * @param numPeaks
     * @return
     */
    protected QuantPeakSetSummary calcPeakIntensities(Feature feature, MSRun run,
                                                      int numPeaks)
    {
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

        float mzTol = (MassUtilities.calculateAbsoluteDeltaMass(
                feature.getMass(), peakPPMTolerance, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / feature.getCharge();
        QuantPeakSetSummary result = new QuantPeakSetSummary();
        result.monoisotopicMass = feature.getMass();
        result.scanRetentionTimes = new double[scans.length];
        for (int i=0; i<scans.length; i++)
        {
            result.scanRetentionTimes[i] = i;
        }
//        result.rawIntensities = new ArrayList<double[]>();
//        for (int i = 0; i<numPeaksRaw; i++)
//            result.rawIntensities.add(new double[scans.length]);
        result.peakSumIntensities = new ArrayList<Float>();
        for (int i=0; i<numPeaks; i++) result.peakSumIntensities.add(0f);

        for (int scanIndex = 0; scanIndex < scans.length; scanIndex++)
        {
            float[][] spectrum = scans[scanIndex].getSpectrum();

            for (int peakIndex=-1; peakIndex<numPeaks; peakIndex++)
            {
                float peakMzCenter = feature.getMz() +
                        ((peakIndex * (float) MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH) /
                                feature.getCharge());

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

//                if (peakIndex < numPeaksRaw)
//                    result.rawIntensities.get(peakIndex)[scanIndex] = intensityMax;
            }
        }
        return result;
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

        protected List<double[]> rawIntensities;

        public QuantPeakSetSummary()
        {
            peakSumIntensities = new ArrayList<Float>();
        }

        public QuantPeakSetSummary(float peakBelow, List<Float> intensitySums, List<double[]> rawIntensities,
                                   float monoisotopicMass)
        {
            sumIntensityPeakBelow = peakBelow;
            peakSumIntensities = intensitySums;
            this.rawIntensities = rawIntensities;
            this.monoisotopicMass = monoisotopicMass;
        }

        double[] combineRawIntensitiesAllPeaks()
        {
            int fullLength = 0;
            for (double[] array : rawIntensities)
                fullLength += array.length;
            double[] result = new double[fullLength];
            int position=0;
            for (double[] array : rawIntensities)
            {
                System.arraycopy(array, 0, result, position, array.length);
                position += array.length;
            }
            return result;
        }
    }


}
