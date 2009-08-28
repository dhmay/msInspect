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
import org.fhcrc.cpl.viewer.quant.QuantEventAssessor;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.ApplicationContext;
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
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;

import java.io.File;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;
import java.util.List;


/**
 * This uses a big HACK, adding a dummy search score to peptide features to store the flag reason description
 *
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

    QuantEventAssessor quantEventAssessor = new QuantEventAssessor();

    protected static Logger _log = Logger.getLogger(FlagQuantEventsCLM.class);

    protected File[] featureFiles;
    protected File outFile;
    protected File outDir;
    protected File mzXmlDir;


    public static final String REASON_DUMMY_SEARCH_SCORE_NAME = "dummy_flag_desc";


    //scaffolding
    protected List<Float> lightPrevPeakRatios = new ArrayList<Float>();
    protected List<Float> heavyPrevPeakRatios = new ArrayList<Float>();


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
                        new BooleanArgumentDefinition("showcharts", false, "Show charts?", false),
                        new DecimalArgumentDefinition("minflagratio", false,
                                "Ratios must be higher than this, or lower than maxflagratio, or both, to flag",
                                0f),
                        new DecimalArgumentDefinition("maxflagratio", false,
                                "Ratios must be lower than this, or higher than maxflagratio, or both, to flag",
                                999f),
                        new DecimalArgumentDefinition("peakppm", false,
                                "Mass tolerance around each theoretical peak (ppm)",
                                QuantEventAssessor.DEFAULT_PEAK_PPM_TOLERANCE),
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

        quantEventAssessor.setLabelType(((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label")));

        quantEventAssessor.setPeakPPMTolerance(getFloatArgumentValue("peakppm"));
        quantEventAssessor.setShowCharts(getBooleanArgumentValue("showcharts"));
        quantEventAssessor.setMinFlagRatio(getFloatArgumentValue("minflagratio"));
        quantEventAssessor.setMaxFlagRatio(getFloatArgumentValue("maxflagratio"));
        if (hasArgumentValue("minflagratio") || hasArgumentValue("maxflagratio"))
            ApplicationContext.infoMessage("NOTE: only ratios higher than " + quantEventAssessor.getMinFlagRatio() +
                    " or lower than " + quantEventAssessor.getMaxFlagRatio() + " (or both) will be flagged");

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


            if (quantEventAssessor.isShowCharts())
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

                if (quantEventAssessor.isShowCharts())
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

        Map<String, Map<Integer, List<Float>>> peptideChargeRatiosMap = new HashMap<String, Map<Integer, List<Float>>>();
        Set<String> flaggedPeptides = new HashSet<String>();

        Arrays.sort(ms2FeatureSet.getFeatures(), new Feature.ScanChargeMzAscComparator());

        int numFeatures = ms2FeatureSet.getFeatures().length;
        ApplicationContext.infoMessage("Processing " + numFeatures + " features....");
        int i=0;
        for (Feature feature : ms2FeatureSet.getFeatures())
        {
            if (numFeatures > 10 && i % (int) (numFeatures/10f) == 0)
                ApplicationContext.infoMessage((10 * (i / (int) (numFeatures/10f))) + "%");
            i++;

            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                continue;

            QuantEventAssessor.QuantEventAssessment assessment = quantEventAssessor.assessFeature(feature, run);
            int flagReason = assessment.getStatus();
            String flagReasonDesc = assessment.getExplanation();
            if (flagReason != QuantEventAssessor.FLAG_REASON_NONE &&
                flagReason != QuantEventAssessor.FLAG_REASON_UNEVALUATED)
            {
                flaggedReasons.add((float) flagReason);                
                feature.setDescription(flagReasonDesc);
                //This is a HACK, adding a dummy search score to store the flag reason description
                MS2ExtraInfoDef.addSearchScore(feature, REASON_DUMMY_SEARCH_SCORE_NAME, flagReasonDesc);
                flaggedFeatures.add(feature);
                flaggedPeptides.add(MS2ExtraInfoDef.getFirstPeptide(feature));

                if (assessment.getSinglePeakRatio() != -1)
                {
                    algRatios.add((float)IsotopicLabelExtraInfoDef.getRatio(feature));
                    singlePeakRatios.add(assessment.getSinglePeakRatio());
                }
            }
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
                    flaggedPeptides.contains(peptide)));
            }
        }

        ms2FeatureSet.setFeatures(flaggedFeatures.toArray(new Feature[flaggedFeatures.size()]));
    }
}
