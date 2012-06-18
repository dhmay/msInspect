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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandModuleUtilities;
import org.fhcrc.cpl.viewer.quant.QuantEventAssessor;
import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.fhcrc.cpl.viewer.quant.turk.TurkUtilities;
import org.fhcrc.cpl.viewer.quant.gui.QuantitationVisualizer;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;

import java.io.*;
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

    protected static final String FEATURE_PROPERTY_QUANTASSESSMENT = "FLAGQUANT_QUANTASSESSMENT";

    public static final String REASON_DUMMY_SEARCH_SCORE_NAME = "SEARCHSCORE_FLAGQUANT_REASON";


//This is scaffolding, for comparing against some known bad events.  Probably no longer effective
List<Integer> badScans = Arrays.asList(new Integer[] {    9603, 6106, 6177, 8008, 1062, 1130, 6923, 7176, 7917, 3139, 3141, 4457, 4673, 2694,});


    List<Float> flaggedReasons = new ArrayList<Float>();
    List<Integer> flaggedScans = new ArrayList<Integer>();
    List<Integer> badFlagged = new ArrayList<Integer>();

    List<Float> reasonsTruePos = new ArrayList<Float>();
    List<Float> reasonsFalsePos = new ArrayList<Float>();

    List<Float> logAlgRatiosWithSinglePeakRatios = new ArrayList<Float>();
    List<Float> weightedGeomMeanRatios = new ArrayList<Float>();

    QuantEventAssessor quantEventAssessor = new QuantEventAssessor();

    protected static Logger _log = Logger.getLogger(FlagQuantEventsCLM.class);

    protected File[] featureFiles;
    protected File outFlaggedFile;
    protected File outNoFlaggedFile;
    protected File outNoFlaggedDir;

    protected File outBadEventTurkFile;
    protected PrintWriter outBadEventTurkPW;
    protected List<QuantEvent> badQuantEvents = new ArrayList<QuantEvent>();
    protected File outBadEventQurateFile;
    protected File outTurkImageDir;

    protected File mzXmlDir;




    //scaffolding
    protected List<Float> lightPrevPeakRatios = new ArrayList<Float>();
    protected List<Float> heavyPrevPeakRatios = new ArrayList<Float>();


    //For calculating statistics on flagged vs non-flagged peptides
    Map<String, List<Pair<Float, Boolean>>> peptidePerFractionRatioFlaggedMap = new HashMap<String, List<Pair<Float,Boolean>>>();

    //all log-ratios and indications of whether flagged or not, for chart
    List<Pair<Float,Boolean>> allLogRatiosFlaggeds = new ArrayList<Pair<Float,Boolean>>();

    Map<String, Integer> peptideNumFractionsFlaggedMap = new HashMap<String, Integer>();
    List<Float> stDevsFromMeanRatioNoFlag = new ArrayList<Float>();
    List<Float> stDevsFromMeanRatioFlag = new ArrayList<Float>();
    List<Float> corrFlagged = new ArrayList<Float>();
    List<Float> corrUnFlagged = new ArrayList<Float>();
    //for sens, spec
    List<Float> corrBad = new ArrayList<Float>();
    List<Float> corrGood = new ArrayList<Float>();

    List<Float> logSinglePeakRatios = new ArrayList<Float>();

    protected boolean shouldWriteBadTurk;
    protected int badTurkId = 0;
    protected QuantitationVisualizer quantVisualizer = null;

    protected Map<String, List<Integer>> fractionScanListMapEventsToProcess = new HashMap<String, List<Integer>>();
    protected boolean onlyProcessSpecifiedEvents = false;

    //scaffolding for calculating ratios using regression based on one datapoint per-scan, like RelEx.
    //I think that method pretty much doesn't work very well.
//    List<Float> singlePeakSlopeRatios = new ArrayList<Float>();
//    List<Float> multiPeakSlopeRatios = new ArrayList<Float>();

    protected Set<String> peptidesWithGoodQuantEvents = new HashSet<String>();
    protected Set<String> proteinsWithGoodQuantEvents = new HashSet<String>();



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
                        createUnnamedSeriesFileArgumentDefinition(true, "pepXML file"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "mzXML directory"),
                        new FileToWriteArgumentDefinition("outflagged", false,
                                "Output pepXML file containing flagged events only"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory for flagged features"),
//                        new EnumeratedValuesArgumentDefinition("label", false, QuantitationUtilities.ALL_LABEL_CODES,
//                               QuantitationUtilities.ALL_LABEL_EXPLANATIONS, QuantitationUtilities.LABEL_LYCINE_CODE),
                        new BooleanArgumentDefinition("showcharts", false, "Show charts?", false),
                        new DecimalArgumentDefinition("minflagratio", false,
                                "Ratios must be higher than this, or lower than maxflagratio, or both, to flag",
                                0f),
                        new DecimalArgumentDefinition("maxflagratio", false,
                                "Ratios must be lower than this, or higher than minflagratio, or both, to flag",
                                999f),
                        new DecimalArgumentDefinition("peakppm", false,
                                "Mass tolerance around each theoretical peak (ppm)",
                                QuantEventAssessor.DEFAULT_PEAK_PPM_TOLERANCE),
                        new FileToWriteArgumentDefinition("outnoflagged", false,
                                "Output pepXML file containing all input features (with and without ratios) " +
                                        "EXCEPT flagged features"),
                        new DirectoryToWriteArgumentDefinition("outbadturkdir", false,
                                "Output directory for Mechanical Turk files for bad events"),
                        new DirectoryToWriteArgumentDefinition("outnoflaggeddir", false,
                                "Output directory for 'outnoflagged' files (for multiple input files)"),
                        //todo: remove this special-purpose arg
                        new FileToReadArgumentDefinition("quratetoprocess", false, "Qurate file.  If this file is present, only process events from this file"),

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        outFlaggedFile = getFileArgumentValue("outflagged");
        File outDir = getFileArgumentValue("outdir");

        if (outFlaggedFile != null)
            assertArgumentAbsent("outdir", "out");

        if (featureFiles.length > 1 &&
                (hasArgumentValue("outdir") || hasArgumentValue("outnoflagged") || hasArgumentValue("outflagged")))
            throw new ArgumentValidationException("arguments ourdir, outflagged and outnoflagged can only be " +
                    "specified with a single input file");

        if (outDir != null)
        {
            outFlaggedFile = new File(outDir, featureFiles[0].getName());
        }

        outNoFlaggedFile = getFileArgumentValue("outnoflagged");
        outNoFlaggedDir = getFileArgumentValue("outnoflaggeddir");
        if (outNoFlaggedDir != null && featureFiles[0].getAbsolutePath().contains(outNoFlaggedDir.getAbsolutePath()))
            ApplicationContext.infoMessage("Warning: is outnoflaggeddir same as directory of input files?  This would cause problems");

        if (outNoFlaggedFile != null && outNoFlaggedDir != null)
             throw new ArgumentValidationException("only outnoflagged or outnoflaggeddir may be specified, not both");


        mzXmlDir = getFileArgumentValue("mzxmldir");

//        quantEventAssessor.setLabelType(((EnumeratedValuesArgumentDefinition)
//                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label")));

        quantEventAssessor.setPeakPPMTolerance(getFloatArgumentValue("peakppm"));
        quantEventAssessor.setShowCharts(getBooleanArgumentValue("showcharts"));
        quantEventAssessor.setMinFlagRatio(getFloatArgumentValue("minflagratio"));
        quantEventAssessor.setMaxFlagRatio(getFloatArgumentValue("maxflagratio"));
        if (hasArgumentValue("minflagratio") || hasArgumentValue("maxflagratio"))
            ApplicationContext.infoMessage("NOTE: only ratios higher than " + quantEventAssessor.getMinFlagRatio() +
                    " or lower than " + quantEventAssessor.getMaxFlagRatio() + " (or both) will be flagged");

        if (hasArgumentValue("outbadturkdir"))
        {
            shouldWriteBadTurk = true;
            outTurkImageDir = getFileArgumentValue("outbadturkdir");
            outBadEventTurkFile = new File(outTurkImageDir, "bad_event_turk.csv");
            outBadEventQurateFile = new File(outTurkImageDir, "bad_event_qurate.tsv");
            try
            {
                outBadEventTurkPW = new PrintWriter(outBadEventTurkFile);
                outBadEventTurkPW.println(TurkUtilities.createTurkHITFileHeaderLine());             
            }
            catch (IOException e)
            {
                throw new ArgumentValidationException("Failed to open file for turk HITs", e);
            }
        }

        //This is pretty special-purpose.  Ignore events unless they appear in a qurate .tsv file
        //todo: remove this when no longer needed
        if (hasArgumentValue("quratetoprocess"))
        {
            onlyProcessSpecifiedEvents = true;
            fractionScanListMapEventsToProcess = new HashMap<String, List<Integer>>();
            try
            {
                List<QuantEvent> quantEvents = QuantEvent.loadQuantEvents(getFileArgumentValue("quratetoprocess"));
                for (QuantEvent quantEvent : quantEvents)
                {
                    List<Integer> scansThisFraction = fractionScanListMapEventsToProcess.get(quantEvent.getFraction());
                    if (scansThisFraction == null)
                    {
                        scansThisFraction = new ArrayList<Integer>();
                        fractionScanListMapEventsToProcess.put(quantEvent.getFraction(), scansThisFraction);
                    }
                    scansThisFraction.add(quantEvent.getScan());
                }
                ApplicationContext.infoMessage("Only processing " + quantEvents.size() + " specified events");
                
            }
            catch (IOException e)
            {
                throw new ArgumentValidationException("Failed to load Qurate file",e);
            }
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        quantVisualizer = new QuantitationVisualizer();
        for (File featureFile : featureFiles)
            handleFile(featureFile);
    }

    protected void handleFile(File featureFile) throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Handling file " + featureFile.getAbsolutePath());

        List<FeatureSet> featureSets = null;
//        Iterator<FeatureSet> featureSetIterator = null;
        try
        {
            featureSets = PepXMLFeatureFileHandler.getSingletonInstance().loadAllFeatureSets(featureFile);

//            featureSetIterator =
//                    new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failure loading pepXML file " + featureFile.getAbsolutePath(),e);
        }

        int flaggedFeaturesThisFile = 0;
        int totalFeaturesThisFile = 0;

        List<File> tempFlaggedFeatureFiles = new ArrayList<File>();

        //only gets populated if we're saving unflagged features
        List<File> tempUnFlaggedFeatureFiles = new ArrayList<File>();

        int numSetsProcessed = 0;
        for (FeatureSet featureSet : featureSets)
        {
//            FeatureSet featureSet = featureSetIterator.next();
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
                if (featureSets.size() > 1)
                    baseName = baseName + "_" + numSetsProcessed;
            }
            //if PeptideProphet was run from a directory below the directory containing the
            //mzXML files, we may have ../ in the baseName, which causes trouble in saving
            //the temporary files
            while (baseName.contains(".." + File.separator))
                baseName = baseName.replaceFirst(".." + File.separator, "");

            int numFeaturesThisFraction = featureSet.getFeatures().length;

//            Pair<FeatureSet, FeatureSet> flaggedUnflaggedSets = processFeatureSet(featureSet, run);
            int numFlaggedFeaturesThisFraction = processFeatureSet(featureSet, run);
            run = null;
            System.gc();

//            FeatureSet flaggedFeatureSet = flaggedUnflaggedSets.first;
//            FeatureSet unflaggedFeatureSet = flaggedUnflaggedSets.second;


//            int numFlaggedFeaturesThisFraction = flaggedFeatureSet.getFeatures().length;

            ApplicationContext.infoMessage("Flagged " + numFlaggedFeaturesThisFraction + " out of " +
                    numFeaturesThisFraction + " features this fraction");
            flaggedFeaturesThisFile += numFlaggedFeaturesThisFraction;


//            if (outFlaggedFile != null)
//            {
//                File thisFractionFlaggedFeatureFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);
//                try
//                {
//                    flaggedFeatureSet.savePepXml(thisFractionFlaggedFeatureFile);
//                    _log.debug("Saved fraction flagged file as " + thisFractionFlaggedFeatureFile.getAbsolutePath());
//                    tempFlaggedFeatureFiles.add(thisFractionFlaggedFeatureFile);
//
//                }
//                catch (IOException e)
//                {
//                    throw new CommandLineModuleExecutionException("Failure writing file " +
//                            thisFractionFlaggedFeatureFile.getAbsolutePath(),e);
//                }
//            }
//
//            if (outNoFlaggedFile != null)
//            {
//                File thisFractionUnFlaggedFeatureFile =
//                        TempFileManager.createTempFile(baseName + ".unflagged.pep.xml", this);
//                try
//                {
//                    unflaggedFeatureSet.savePepXml(thisFractionUnFlaggedFeatureFile);
//                    _log.debug("Saved fraction unflagged file as " + thisFractionUnFlaggedFeatureFile.getAbsolutePath());
//                    tempUnFlaggedFeatureFiles.add(thisFractionUnFlaggedFeatureFile);
//                }
//                catch (IOException e)
//                {
//                    throw new CommandLineModuleExecutionException("Failure writing file " +
//                            thisFractionUnFlaggedFeatureFile.getAbsolutePath(),e);
//                }
//            }

            numSetsProcessed++;
        }

        ApplicationContext.infoMessage("Flagged " + flaggedFeaturesThisFile + " out of " + totalFeaturesThisFile +
                " features (" + ((float) flaggedFeaturesThisFile *100f / (float) totalFeaturesThisFile) +
                "%) in this file");

        int numFlaggedWithoutOtherPeptideSupport = 0;
        Set<String> peptidesFlaggedWithoutOtherSupport = new HashSet<String>();
        Set<String> allQuantPeptides = new HashSet<String>();
        int numFlaggedWithoutOtherProteinSupport = 0;
        Set<String> proteinsFlaggedWithoutOtherSupport = new HashSet<String>();
        Set<String> allQuantProteins = new HashSet<String>();
        for (FeatureSet featureSet : featureSets)
        {
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                    allQuantPeptides.add(peptide);
                    allQuantProteins.addAll(MS2ExtraInfoDef.getProteinList(feature));
                    QuantEventAssessor.QuantEventAssessment assessment = (QuantEventAssessor.QuantEventAssessment)feature.getProperty(FEATURE_PROPERTY_QUANTASSESSMENT);
                    if (assessment.getStatus() != QuantEventAssessor.FLAG_REASON_OK)
                    {
                        allQuantPeptides.add(peptide);
                        if (!peptidesWithGoodQuantEvents.contains(peptide))
                        {
                            numFlaggedWithoutOtherPeptideSupport++;
                            peptidesFlaggedWithoutOtherSupport.add(peptide);
                        }
                        boolean hasUnsupportedProtein = false;
                        for (String protein : MS2ExtraInfoDef.getProteinList(feature))
                        {
                            if (!proteinsWithGoodQuantEvents.contains(protein))
                            {
                                hasUnsupportedProtein = true;
                                proteinsFlaggedWithoutOtherSupport.add(protein);
                            }
                        }
                        if (hasUnsupportedProtein) numFlaggedWithoutOtherProteinSupport++;
                    }
                }
            }
        }
        ApplicationContext.infoMessage("Overall, " + numFlaggedWithoutOtherPeptideSupport + " out of " + flaggedFeaturesThisFile +
               " flagged features (" + (100f * numFlaggedWithoutOtherPeptideSupport/flaggedFeaturesThisFile) + "%) have no unflagged ratio peptide support");
        ApplicationContext.infoMessage(peptidesFlaggedWithoutOtherSupport.size() + " out of " + allQuantPeptides.size() + " quant peptides would lose ratios if flagged removed");
        ApplicationContext.infoMessage("Overall, " + numFlaggedWithoutOtherProteinSupport + " out of " + flaggedFeaturesThisFile +
               " flagged features (" + (100f * numFlaggedWithoutOtherProteinSupport/flaggedFeaturesThisFile) + "%) have no unflagged ratio protein support");
        ApplicationContext.infoMessage(proteinsFlaggedWithoutOtherSupport.size() + " out of " + allQuantProteins.size() + " quant proteins would lose ratios if flagged removed");



        //combine flagged feature files
        if (outFlaggedFile != null)
        {
            ApplicationContext.infoMessage("Saving flagged features...");            
            List<FeatureSet> flaggedFeatureSets = new ArrayList<FeatureSet>();
            for (FeatureSet featureSet : featureSets)
            {
                FeatureSet flaggedFeatureSet = (FeatureSet) featureSet.clone();
                List<Feature> keptFeatures = new ArrayList<Feature>();
                for (Feature feature : flaggedFeatureSet.getFeatures())
                {
                    boolean keep = false;
                    if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        QuantEventAssessor.QuantEventAssessment assessment = (QuantEventAssessor.QuantEventAssessment)feature.getProperty(FEATURE_PROPERTY_QUANTASSESSMENT);
                        if (assessment.getStatus() != QuantEventAssessor.FLAG_REASON_OK)
                            keep = true;
                    }
                    if (keep)
                        keptFeatures.add(feature);
                }
                flaggedFeatureSet.setFeatures(keptFeatures.toArray(new Feature[keptFeatures.size()]));
                flaggedFeatureSets.add(flaggedFeatureSet);
            }
            try
            {
                PepXMLFeatureFileHandler.getSingletonInstance().saveFeatureSets(flaggedFeatureSets, outFlaggedFile);
                            ApplicationContext.infoMessage("Saved unflagged features to file " +
                    outFlaggedFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save output file " + outFlaggedFile.getAbsolutePath(), e);
            }
        }

        File outNoFlaggedFileThisFile = outNoFlaggedFile;
        if (outNoFlaggedFileThisFile == null && outNoFlaggedDir != null)
        {
            outNoFlaggedFileThisFile = new File(outNoFlaggedDir, featureFile.getName());
        }

        //combine unflagged feature files
        if (outNoFlaggedFileThisFile != null)
        {
            ApplicationContext.infoMessage("Saving unflagged features to file " + outNoFlaggedFileThisFile + " ...");
            List<FeatureSet> unflaggedFeatureSets = new ArrayList<FeatureSet>();
            for (FeatureSet featureSet : featureSets)
            {
                FeatureSet unflaggedFeatureSet = (FeatureSet) featureSet.clone();
                List<Feature> keptFeatures = new ArrayList<Feature>();
                for (Feature feature : unflaggedFeatureSet.getFeatures())
                {
                    boolean keep = true;

                    if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        QuantEventAssessor.QuantEventAssessment assessment = (QuantEventAssessor.QuantEventAssessment)feature.getProperty(FEATURE_PROPERTY_QUANTASSESSMENT);
                        if (assessment.getStatus() != QuantEventAssessor.FLAG_REASON_OK)
                            keep = false;
                    }
                    if (keep)
                        keptFeatures.add(feature);
                }
                unflaggedFeatureSet.setFeatures(keptFeatures.toArray(new Feature[keptFeatures.size()]));
                unflaggedFeatureSets.add(unflaggedFeatureSet);

            }
            try
            {
                PepXMLFeatureFileHandler.getSingletonInstance().saveFeatureSets(unflaggedFeatureSets, outNoFlaggedFileThisFile);
            ApplicationContext.infoMessage("Saved unflagged features to file " +
                    outNoFlaggedFileThisFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save output file " + outNoFlaggedFileThisFile.getAbsolutePath(), e);
            }

        //combine flagged feature files
//        if (outFlaggedFile != null)
//        {
//            ApplicationContext.infoMessage("Saving flagged features...");
//            try
//            {
//                if (numSetsProcessed == 1)
//                {
//                    FileReader in = new FileReader(tempFlaggedFeatureFiles.get(0));
//                    FileWriter out = new FileWriter(outFlaggedFile);
//                    int c;
//
//                    while ((c = in.read()) != -1)
//                        out.write(c);
//
//                    in.close();
//                    out.close();
//                }
//                else
//                {
//                    _log.debug("\tCombining individual fraction files... ");
//                    new PepXMLFeatureFileHandler().combinePepXmlFiles(tempFlaggedFeatureFiles, outFlaggedFile);
//                    ApplicationContext.infoMessage("Saved flagged features to file " +
//                            outFlaggedFile.getAbsolutePath());
//                }
//            }
//            catch (IOException e)
//            {
//                throw new CommandLineModuleExecutionException("Failed to save output file " + outFlaggedFile.getAbsolutePath(), e);
//            }
//        }

        //combine unflagged feature files
//        if (outNoFlaggedFile != null)
//        {
//            ApplicationContext.infoMessage("Saving unflagged features...");
//            try
//            {
//                if (numSetsProcessed == 1)
//                {
//                    FileReader in = new FileReader(tempUnFlaggedFeatureFiles.get(0));
//                    FileWriter out = new FileWriter(outNoFlaggedFile);
//                    int c;
//
//                    while ((c = in.read()) != -1)
//                        out.write(c);
//
//                    in.close();
//                    out.close();
//                }
//                else
//                {
//                    _log.debug("\tCombining individual fraction files... ");
//                    new PepXMLFeatureFileHandler().combinePepXmlFiles(tempUnFlaggedFeatureFiles, outNoFlaggedFile);
//                    ApplicationContext.infoMessage("Saved unflagged features to file " +
//                            outNoFlaggedFile.getAbsolutePath());
//                }
//            }
//            catch (IOException e)
//            {
//                throw new CommandLineModuleExecutionException("Failed to save unflagged output file " +
//                        outNoFlaggedFile.getAbsolutePath(), e);
//            }
        }
        ApplicationContext.infoMessage("Done.");

        
        if (shouldCalcSensSpecWithHardcodedScans && !badScans.isEmpty())
        {
            System.err.println("Bad: " + badScans.size() + ", bad flagged: " + badFlagged.size() +
                    ", sens=" + ((float)badFlagged.size() / (float)badScans.size()) + ", spec=" + (((float) flaggedFeaturesThisFile - (float)badFlagged.size()) / (float) flaggedFeaturesThisFile));
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

            new PanelWithScatterPlot(cvs, numFlagged, "CV vs num fractions flagged").displayInTab();

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

            new PanelWithScatterPlot(logAlgRatiosWithSinglePeakRatios,
                    logSinglePeakRatios, "alg vs singlepeak").displayInTab();



            Collections.sort(allLogRatiosFlaggeds, new Comparator<Pair<Float,Boolean>>()
            {
                public int compare(Pair<Float,Boolean> o1, Pair<Float,Boolean> o2)
                {
                    return o1.first > o2.first ? 1 : (o1.first == o2.first ? 0 : -1);
                }
            });

            int rollingWindowSize = 30;
            List<Boolean> windowAssessments = new ArrayList<Boolean>(rollingWindowSize);
            List<Float> windowLogRatios = new ArrayList<Float>(rollingWindowSize);

            List<Float> allRatiosWithWindows = new ArrayList<Float>();
            List<Float> allProportionsBad = new ArrayList<Float>();
            List<Float> allBadLogRatios = new ArrayList<Float>();
            List<Float> allGoodLogRatios = new ArrayList<Float>();

            for (Pair<Float,Boolean> ratioAndAssessment : allLogRatiosFlaggeds)
            {
                if (ratioAndAssessment.second)
                    allBadLogRatios.add(ratioAndAssessment.first);
                else
                    allGoodLogRatios.add(ratioAndAssessment.first);

                windowAssessments.add(ratioAndAssessment.second);
                windowLogRatios.add(ratioAndAssessment.first);
                if (windowAssessments.size() > rollingWindowSize)
                {
                    windowAssessments.remove(0);
                    windowLogRatios.remove(0);
                }
                if (windowAssessments.size() == rollingWindowSize)
                {
                    int numBad = 0;
                    for (boolean assessment : windowAssessments)
                        if (assessment)
                        {
                            numBad++;
                        }
                    float medianLogRatio = (float) BasicStatistics.median(windowLogRatios);
                    if (!Float.isInfinite(medianLogRatio))
                    {
                        allRatiosWithWindows.add((float) BasicStatistics.median(windowLogRatios));
                        allProportionsBad.add((float) numBad / (float) rollingWindowSize);
                    }
                }
            }
            PanelWithLineChart pwlcProportionBad = new PanelWithLineChart();
            pwlcProportionBad.setName("LogRatio vs. Proportion Bad");
            pwlcProportionBad.addData(allRatiosWithWindows, allProportionsBad, "LogRatio vs. Proportion Bad");
            pwlcProportionBad.displayInTab();

            new PanelWithHistogram(allBadLogRatios, "BadLogRatios").displayInTab();
            new PanelWithHistogram(allGoodLogRatios, "GoodLogRatios").displayInTab();

        }


        if (shouldWriteBadTurk)
        {
            outBadEventTurkPW.close();
            try
            {
                QuantEvent.saveQuantEventsToTSV(badQuantEvents, outBadEventQurateFile, true, false);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save Qurate file for bad events", e);
            }
        }
    }

    /**
     * @param ms2FeatureSet
     * @param run
     * @throws CommandLineModuleExecutionException
     * @return a pair of FeatureSets, one containing flagged and one containing unflagged features
     */
    protected int processFeatureSet(FeatureSet ms2FeatureSet, MSRun run)
            throws CommandLineModuleExecutionException
    {
//        List<Feature> flaggedFeatures = new ArrayList<Feature>();
//        List<Feature> unflaggedFeatures = new ArrayList<Feature>();


        Map<String, Map<Integer, List<Float>>> peptideChargeRatiosMap = new HashMap<String, Map<Integer, List<Float>>>();
        Set<String> flaggedPeptides = new HashSet<String>();

        Arrays.sort(ms2FeatureSet.getFeatures(), new Feature.ScanChargeMzAscComparator());

        int numFeatures = ms2FeatureSet.getFeatures().length;
        ApplicationContext.infoMessage("Processing " + numFeatures + " features....");
        int i=0;

        String featureSetBaseName = MS2ExtraInfoDef.getFeatureSetBaseName(ms2FeatureSet);
        if (onlyProcessSpecifiedEvents && !fractionScanListMapEventsToProcess.containsKey(featureSetBaseName))
        {
            ApplicationContext.infoMessage("WARNING!!! Fraction " + featureSetBaseName + " not found in qurate file.  Skipping.");
            FeatureSet dummy = (FeatureSet)ms2FeatureSet.clone();
            dummy.setFeatures(new Feature[0]);
            return 0;
        }
        int numFlaggedFeatures = 0;
        for (Feature feature : ms2FeatureSet.getFeatures())
        {
            //todo: remove this when done with it (for turk)
            if (onlyProcessSpecifiedEvents && !fractionScanListMapEventsToProcess.get(featureSetBaseName).contains(feature.getScan()))
                continue;


            if (numFeatures > 10 && i % (int) (numFeatures/10f) == 0)
                ApplicationContext.infoMessage((10 * (i / (int) (numFeatures/10f))) + "%");
            i++;

            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
//                unflaggedFeatures.add(feature);
                continue;
            }

            QuantEvent quantEvent = new QuantEvent(feature, MS2ExtraInfoDef.getFeatureSetBaseName(ms2FeatureSet));
            QuantEventAssessor.QuantEventAssessment assessment = quantEventAssessor.assessQuantEvent(quantEvent, run);
            feature.setProperty(FEATURE_PROPERTY_QUANTASSESSMENT, assessment);
            int flagReason = assessment.getStatus();
            String flagReasonDesc = assessment.getExplanation();

            float algRatio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
            allLogRatiosFlaggeds.add(new Pair<Float, Boolean>((float)Math.log(algRatio),
                    flagReason != QuantEventAssessor.FLAG_REASON_OK));
            if (flagReason == QuantEventAssessor.FLAG_REASON_OK)
            {
                peptidesWithGoodQuantEvents.add(MS2ExtraInfoDef.getFirstPeptide(feature));
                proteinsWithGoodQuantEvents.addAll(MS2ExtraInfoDef.getProteinList(feature));
//                unflaggedFeatures.add(feature);
            }
            else
            {                        
                badQuantEvents.add(quantEvent);
                flaggedReasons.add((float) flagReason);                
                feature.setDescription(flagReasonDesc);
                numFlaggedFeatures++;
//                flaggedFeatures.add(feature);
                flaggedPeptides.add(MS2ExtraInfoDef.getFirstPeptide(feature));

                if (assessment.getSinglePeakRatio() >0 && algRatio > 0)
                {
                    logAlgRatiosWithSinglePeakRatios.add((float)Math.log(algRatio));
                    logSinglePeakRatios.add((float) Math.log(assessment.getSinglePeakRatio()));
                }
                if (shouldWriteBadTurk)
                {
                    PanelWithSpectrumChart spectrumPanel = quantVisualizer.createPanelWithSpectrumChart(run, quantEvent);
                    try
                    {
                        String turkFileLine = quantVisualizer.saveTurkImage(quantEvent, spectrumPanel,
                                outTurkImageDir, badTurkId++);
                        outBadEventTurkPW.println(turkFileLine);
                    }
                    catch (IOException e)
                    {
                        throw new CommandLineModuleExecutionException("Failed to create Turk image", e);
                    }
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

//        FeatureSet flaggedFeatureSet = (FeatureSet) ms2FeatureSet.clone();
//        flaggedFeatureSet.setFeatures(flaggedFeatures.toArray(new Feature[flaggedFeatures.size()]));
//
//        FeatureSet unflaggedFeatureSet = (FeatureSet) ms2FeatureSet.clone();
//        unflaggedFeatureSet.setFeatures(unflaggedFeatures.toArray(new Feature[unflaggedFeatures.size()]));

        return numFlaggedFeatures;
    }
}
