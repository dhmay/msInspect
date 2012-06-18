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
package org.fhcrc.cpl.viewer.amt.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.toolbox.gui.chart.MultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.IOException;


/**
 * Command linemodule for matching AMT databases to MS1 data from various sources
 */
public class AmtDatabaseMatcherCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log =
            Logger.getLogger(AmtDatabaseMatcherCLM.class);

    boolean dummyMatch = false;

    protected File outFile;
    protected File outDir;
    protected File inDir;
    protected File ms2Dir;
    protected File mzXmlDir;
    protected File amtDBFile;

    protected AmtDatabase amtDB;

    protected AmtDatabaseMatcher amtDatabaseMatcher;



    protected FeatureSet ms1FeatureSet;
    protected FeatureSet embeddedMs2FeatureSet;
    protected Feature[] amtDBBaseFeatures;

    protected float looseDeltaMass = 20; 
    protected int deltaMassType = AmtDatabaseMatcher.DEFAULT_2D_MATCH_DELTA_MASS_TYPE;
    protected float looseDeltaElution = .15f;
    protected int rTimeoutMilliseconds = AmtMatchProbabilityAssigner.MAX_R_PROB_ASSIGNMENT_MILLIS;


    protected boolean useAmtDBMods = false;
    protected MS2Modification[] ms2ModificationsForMatching;

    protected boolean calibrateMassesUsingMatches = true;

    //this is for a mode that calculates the sensitivity of the model.
    protected boolean calcFDR = false;

    public static final MS2Modification[] defaultMS2ModificationsForMatching =
            {
                    ModificationListArgumentDefinition.safeCreateModificationFromString("C57.021"),
                    ModificationListArgumentDefinition.safeCreateModificationFromString("M16V"),
            };

    //for embedded ms2
    public static final double DEFAULT_MIN_EMBEDDED_MS2_PEPTIDE_PROPHET = .9;
    protected double minMS2PeptideProphet = DEFAULT_MIN_EMBEDDED_MS2_PEPTIDE_PROPHET;
    protected MSRun run = null;


    protected boolean removeFractions = false;

//    protected double amtFeatureMassAdjustment = 0;

    //only if we're removing runs from the database
    protected int minRunPeptideMatchPercent = 0;
    protected int maxEntriesInMassMatchedDatabase = Integer.MAX_VALUE;
    protected int maxRunsInMassMatchedDatabase = Integer.MAX_VALUE;

    protected boolean restrictDatabaseRuns = false;

    protected Fractionation2DUtilities.FractionatedAMTDatabaseStructure amtDatabaseStructure;

    //mode stuff
    protected int mode=-1;

    protected final static String[] modeStrings =
            {
                    "singlems1",
                    "ms1dir",
            };

    protected final static String[] modeExplanations =
            {
                    "Match against a single MS1 feature file",
                    "Match against a directory of MS1 files, one by one",
            };
    protected static final int MODE_MATCH_AGAINST_MS1 = 0;
    protected static final int MODE_MATCH_AGAINST_MS1_DIR = 1;


    protected int nonlinearMappingPolynomialDegree = AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;


    protected boolean showCharts = false;
    protected File saveChartsDir = null;

    protected boolean nonlinearAlignment = true;

    protected int minRunsToKeep = 1;
    protected int maxRunsToKeep = 0;

    protected Set<String> thisRunMatchedPeptides = new HashSet<String>();
    protected Set<String> allMatchedPeptides = new HashSet<String>();

    protected boolean populateMs2Times = true;


    public AmtDatabaseMatcherCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matchamt";
        mShortDescription = "Match between AMT databases and MS1 feature files";
        mHelpMessage = "This tool matches MS1 features to an AMT database using mass and Normalized Retention Time.  " +
                "It assigns peptide IDs to the features matched, with confidence values from 0 to 1.\n" +
                "If embedded MS2 features are available for the MS1 features, they can be used to establish " +
                "the Time->Hydrophobicity mapping -- either via the MS2 times themselves or via the times of MS1 " +
                "features that match to them.\n" +
                "Otherwise, mapping is done using a two-step process that starts with strict mass-only matching " +
                "and then uses quantile regression to find the best mapping.";

        CommandLineArgumentDefinition[] basicArgDefs =
                {
                        new EnumeratedValuesArgumentDefinition("mode",true,modeStrings,
                                modeExplanations),
                        createUnnamedFileArgumentDefinition(
                                true, "AMT database for matching"),
                        new FileToWriteArgumentDefinition("out",false,
                                "Output filepath for matching results (for 'singlems1' mode)"),
                        new DirectoryToWriteArgumentDefinition("outdir",false,
                                "Output directory for all matching result files (for 'ms1dir' mode)"),
                        new FeatureFileArgumentDefinition("ms1",false,
                                "Input MS1 feature file for matching (for 'singlems1' mode)"),
                        new DirectoryToReadArgumentDefinition("ms1dir",false,
                                "Directory full of ms1 feature files for matching (for 'ms1dir' mode)"),
                        new FeatureFileArgumentDefinition("embeddedms2",false,
                                "Embedded MS2 feature file from the same run that the MS1 features from matching " +
                                "are from.  This will be used to help develop the mapping between Retention Time and " +
                                "Retention Time and Normalized Retention Time (the scale of the AMT database) for " +
                                "this run. (for 'singlems1' mode)"),
                        new DirectoryToReadArgumentDefinition("ms2dir",false,
                                "Embedded MS2 feature file (e.g., pepxml) directory. For 'ms1dir' mode. For each MS1 " +
                                "file, the corresponding MS2 feature file will be located in this directory by " +
                                "filename and used to aid in the mapping between Retention Time and Normalized " +
                                "Retention Time for the run."),
                        new DirectoryToReadArgumentDefinition("mzxmldir",false,
                                "Directory containing mzXML files from the same runs as the MS1 feature files. " +
                                        "For each embedded MS2 feature file, the corresponding mzXML file will be " +
                                        "located by filename and used to populate MS2 scan times, for use in " +
                                        "developing the RT->NRT map.  This argument is not necessary if the MS2 " +
                                        "feature files already contain retention times (see the 'populatems2times' " +
                                        "command). (for 'ms1dir' mode)"),
                        new DecimalArgumentDefinition("minpprophet",false,
                                "Minimum PeptideProphet score to use from the embedded MS2 feature file(s) in " +
                                        "building the RT->NRT map",
                                minMS2PeptideProphet),
                        new FileToReadArgumentDefinition("mzxml",false,
                                "mzXML file used to populate MS2 scan times, for use in " +
                                        "developing the RT->NRT map.  This argument is not necessary if the MS2 " +
                                        "feature file already contains retention times (see the 'populatems2times' " +
                                        "command). (for 'singlems1' mode)"),
                        new DecimalArgumentDefinition("loosedeltaelution", false,
                                "A loose deltaElution value.  The EM probability model will be " +
                                        "fit using all AMT matches in the loose match.  The appropriate value for this " +
                                        "parameter may depend somewhat on data, so you should adjust it accordingly " +
                                        "if you find true matches being excluded or all of the true matches clustering " +
                                        "near the center of the RT distribution.",
                                looseDeltaElution),
                        new DeltaMassArgumentDefinition("loosedeltamass", false,
                                "A loose mass tolerance for the initial AMT match.  The EM probability model will be " +
                                        "fit using all AMT matches in the loose match.  The appropriate value for this " +
                                        "parameter is largely data-independent.",
                                new DeltaMassArgumentDefinition.DeltaMassWithType(looseDeltaMass, deltaMassType)),                        
                        new DeltaMassArgumentDefinition("massmatchdeltamass", false,
                                "A mass tolerance value to be used when developing the RT->NRT map based on " +
                                        "mass-only matching.  This value is not used if the 'embeddedms2' or " +
                                        "'ms2dir' argument is specified.",
                                new DeltaMassArgumentDefinition.DeltaMassWithType(
                                        AmtDatabaseMatcher.DEFAULT_MASS_MATCH_DELTA_MASS,
                                        AmtDatabaseMatcher.DEFAULT_MASS_MATCH_DELTA_MASS_TYPE)),
                        new ModificationListArgumentDefinition("modifications", false,
                                "A list of modifications to use in matching.  This list should contain all of the " +
                                        "expected modifications in your MS1 data.  During matching, AMT feature " +
                                        "masses will be adjusted to account for any static mods specified here, " +
                                        "and multiple features will be created to account for any variable mods.  " +
                                        "The default value is appropriate for non-isotopically-labeled MS1 data with " +
                                        "iodoacetylated Cysteine and possibly oxidized Methionine.", 
                                defaultMS2ModificationsForMatching),
                        new BooleanArgumentDefinition("showcharts", false,
                                "Show useful charts created when matching?  Not recommended when matching large " +
                                        "numbers of files",
                                showCharts),
                        new DirectoryToReadArgumentDefinition("savechartsdir", false,
                                "Directory to save charts to.  This can be used with or without 'showcharts'"),
                        new BooleanArgumentDefinition("calibratematches",false,
                                "Calibrate MS1 feature masses using AMT matches?  This is a good idea if you " +
                                        "suspect there might be a miscalibration in the MS1 data.  Even a small " +
                                        "miscalibration can have a significant effect on matching by breaking " +
                                        "the assumption that mass match error is distributed normally.", 
                                calibrateMassesUsingMatches),
                        new DecimalArgumentDefinition("minmatchprob", false,
                                "Minimum AMT match probability to keep main output",
                                AmtMatchProbabilityAssigner.DEFAULT_MIN_MATCH_PROBABILITY),
                        new DecimalArgumentDefinition("maxmatchfdr", false,
                                "Maximum AMT match FDR to keep in output",
                                AmtMatchProbabilityAssigner.DEFAULT_MAX_MATCH_FDR),
                        new BooleanArgumentDefinition("usems1foralignment", false,
                                "Use MS1 times, rather than MS2 times, for alignment?  This is done by matching " +
                                        "MS1 and MS2 in a tight window",
                                AmtDatabaseMatcher.DEFAULT_USE_MS1_TIMES_FOR_ALIGNMENT),


                };

        CommandLineArgumentDefinition[] advancedArgDefs =
                {
                        new IntegerArgumentDefinition("minfractionstokeep", false,
                                "Minimum number of fractions to keep (for \"removefractions\")",
                                minRunsToKeep),
                        new IntegerArgumentDefinition("maxfractionstokeep", false,
                                "Maximum number of fractions to keep (for \"removefractions\").  Default is " +
                                        "actually the number of runs in the database",
                                maxRunsToKeep),
                        new BooleanArgumentDefinition("calcFDR", false,
                                "Calculate FDR for all results.  We do this by making half of " +
                                "the database into decoy features, then repeating with the other half.  This is for " +
                                        "purposes of evaluating the EM model _only_!  The probabilities calculated " +
                                        "this way are not as accurate as the probabilities calculated without " +
                                        "the decoy hits, because the decoy hits add to the background complexity " +
                                        "of the null distribution.",
                                calcFDR),
                        new BooleanArgumentDefinition("removefractions", false,
                                "Remove fractions from the database that are unlike the MS1 matching fraction.  " +
                                        "This is only recommended for very dense, extensively fractionated " +
                                        "databases built from hundreds of runs.",
                                removeFractions),
                        new StringArgumentDefinition("amtdbstructure", false,
                                "For multi-fraction AMT databases from one or more experiments.  Defines the " +
                                        "arrangement of runs within the AMT database.  This is only used to produce " +
                                        "fancy heatmap charts of matches.  Of the format " +
                                        "'#rows,#cols,row|col,#experiments', e.g., '12,11,col,2' for a " +
                                        "database with two experiments, each with 12 rows and 11 columns, whose runs " +
                                        "are in order by column within the database"),
                        new IntegerArgumentDefinition("mappingpolynomialdegree", false,
                                "The degree of the polynomial to fit when mapping time to hydrophobicity " +
                                        "nonlinearly. If you notice the mapping overfitting, you may reduce " +
                                        "this value.  If you think the mapping is not capturing all of the nonlinear " +
                                        "quirks of the data, try increasing it.",
                                nonlinearMappingPolynomialDegree),
                        new BooleanArgumentDefinition("dummymatch", false,
                                "Do a dummy match against a mass-shifted database, rather than a real match.  " +
                                        "This is only used for visualizing the false match density.",
                                dummyMatch),
                        new DecimalArgumentDefinition("maxregressionleverage", false,
                                "Maximum leverage /numerator/ (denominator is N) for features included in the modal " +
                                        "regression to map RT to Hydrophobicity.  If you have spurious features at " +
                                        "the beginning or end of the gradient, you may want to reduce this value.",
                                AmtDatabaseMatcher.DEFAULT_MAX_LEVERAGE_NUMERATOR),
                        new DecimalArgumentDefinition("maxregressionstudres", false,
                                "Maximum studentized residual for features included in the modal regression to map " +
                                        "RT to Hydrophobicity.  You may want to decrease this value if many " +
                                        "spurious matches are throwing off the regression, or increase it if " +
                                        "legitimate features are being excluded.",
                                AmtDatabaseMatcher.DEFAULT_MAX_STUDENTIZED_RESIDUAL),
                        new DecimalArgumentDefinition("maxsecondbestprob", false,
                                "Maximum probability of the second-best AMT match, in order to keep best match",
                                AmtMatchProbabilityAssigner.DEFAULT_MAX_SECONDBEST_PROBABILITY),
                        new DecimalArgumentDefinition("minsecondbestprobdiff", false,
                                "Minimum difference between best and secodn-best probability, " +
                                        "in order to keep best match",
                                AmtMatchProbabilityAssigner.DEFAULT_MIN_SECONDBEST_PROBABILITY_DIFFERENCE),
                        new IntegerArgumentDefinition("minemiterations", false,
                                "Minimum number of iterations for the EM algorithm deciding probability values",
                                AmtMatchProbabilityAssigner.DEFAULT_MIN_EM_ITERATIONS),
                        new IntegerArgumentDefinition("maxemiterations", false,
                                "Maximum number of iterations for the EM algorithm deciding probability values",
                                AmtMatchProbabilityAssigner.DEFAULT_MAX_EM_ITERATIONS),
                        new IntegerArgumentDefinition("maxrmillis", false,
                                "Maximum number of milliseconds to wait for R to assign match probabilities",
                                rTimeoutMilliseconds),
                        new DecimalArgumentDefinition("deltamassms1ms2ppm", false,
                                "Mass tolerance for MS1 feature match with MS2 in order to retrieve MS1 feature times",
                                AmtDatabaseBuilder.DEFAULT_MS1_MS2_MASS_TOLERANCE_PPM),
                        new DecimalArgumentDefinition("deltatimems1ms2", false,
                                "Time tolerance (in seconds) for MS1 feature match with MS2 in order to retrieve MS1 feature " +
                                        "times", AmtDatabaseBuilder.DEFAULT_MS1_MS2_TIME_TOLERANCE_SECONDS),
                };

        addArgumentDefinitions(basicArgDefs);
        addArgumentDefinitions(advancedArgDefs, true);


    }



    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));



        dummyMatch = getBooleanArgumentValue("dummymatch");
        calcFDR = getBooleanArgumentValue("calcfdr");
        if (calcFDR && mode ==MODE_MATCH_AGAINST_MS1_DIR)
             throw new ArgumentValidationException("calcsensitivity is only allowed in single-file mode");

        removeFractions = getBooleanArgumentValue("removefractions");

        minMS2PeptideProphet = getDoubleArgumentValue("minpprophet");
        ms1FeatureSet = getFeatureSetArgumentValue("ms1");

        populateMs2Times =  hasArgumentValue("mzxml") || hasArgumentValue("mzxmldir");
        mzXmlDir = getFileArgumentValue("mzxmldir");

        calibrateMassesUsingMatches =  getBooleanArgumentValue("calibratematches");



//        amtFeatureMassAdjustment = getDoubleArgumentValue("amtfeaturemassadjustment");

        embeddedMs2FeatureSet = getFeatureSetArgumentValue("embeddedms2");


        nonlinearMappingPolynomialDegree = getIntegerArgumentValue("mappingpolynomialdegree");


        if (hasArgumentValue("amtdbstructure"))
        {
            amtDatabaseStructure =
                    Fractionation2DUtilities.digestAmtStructureArgument(getStringArgumentValue("amtdbstructure"));
            ApplicationContext.infoMessage("AMT database structure: " + amtDatabaseStructure);
        }        

        if (hasArgumentValue("modifications"))
        {
            MS2Modification[] specifiedMods =
                    getModificationListArgumentValue("modifications");
            if (specifiedMods != null)
            {
                for (MS2Modification specifiedMod : specifiedMods)
                {
                    _log.debug("Including user-specified modification: " + specifiedMod);
                }
                ms2ModificationsForMatching = specifiedMods;
            }
        }
        else
        {
            ms2ModificationsForMatching = defaultMS2ModificationsForMatching;
        }

        _log.debug("Total modifications for matching: " + ms2ModificationsForMatching.length);



        amtDBFile = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");


        ms2Dir = getFileArgumentValue("ms2dir");
        inDir = getFileArgumentValue("ms1dir");




        

        switch (mode)
        {
            case MODE_MATCH_AGAINST_MS1:
                assertArgumentPresent("ms1");
                assertArgumentAbsent("ms1dir");
                if (hasArgumentValue("ms2dir"))
                {
                    assertArgumentAbsent("embeddedms2");
                    try
                    {
                        embeddedMs2FeatureSet = findCorrespondingFeatureSet(ms1FeatureSet.getSourceFile(), ms2Dir);
                    }
                    catch (IOException e)
                    {
                        throw new ArgumentValidationException(e);
                    }

                }
                else if (hasArgumentValue("embeddedms2"))
                    assertArgumentAbsent("ms2dir");
//                else assertArgumentPresent("embeddedms2","mode");
                break;
            case MODE_MATCH_AGAINST_MS1_DIR:
                assertArgumentPresent("ms1dir");
                assertArgumentAbsent("ms1");
                assertArgumentAbsent("embeddedms2");
                assertArgumentAbsent("mzxml");
                break;
        }

        File mzXmlFile = null;
        if (embeddedMs2FeatureSet != null)
        {
            //need times from run
            if (populateMs2Times)
            {
                if (hasArgumentValue("mzxml"))
                {
                    mzXmlFile = getFileArgumentValue("mzxml");
                }
                else
                {
                    try
                    {
                        mzXmlFile =findCorrespondingMzXmlFile(embeddedMs2FeatureSet.getSourceFile());
                    }
                    catch (IOException e)
                    {
                        throw new ArgumentValidationException(e.getMessage());
                    }
                }

                try
                {
                    run = MSRun.load(mzXmlFile.getAbsolutePath());
                    embeddedMs2FeatureSet.populateTimesForMS2Features(run);

                }
                catch (Exception e)
                {
                    throw new ArgumentValidationException("failed to load run from file " + mzXmlFile.getAbsolutePath());
                }
            }
        }

        try
        {
            amtDB = new AmtXmlReader(amtDBFile).getDatabase();
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }

        minRunsToKeep = getIntegerArgumentValue("minfractionstokeep");
        maxRunsToKeep = getIntegerArgumentValue("maxfractionstokeep");
        if (!hasArgumentValue("maxfractionstokeep"))
            maxRunsToKeep = amtDB.numRuns();

//        adjustAmtMassesForMassDefect = getBooleanArgumentValue("correctamtmasses");




        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassWithType =
                getDeltaMassArgumentValue("loosedeltamass");
        looseDeltaMass = deltaMassWithType.getDeltaMass();
        deltaMassType = deltaMassWithType.getDeltaMassType();
        looseDeltaElution = getFloatArgumentValue("loosedeltaelution");



        showCharts = getBooleanArgumentValue("showCharts");
        saveChartsDir = getFileArgumentValue("savechartsdir");

        if (saveChartsDir != null)
        {
            ApplicationContext.infoMessage("Charts will be saved to directory " + saveChartsDir);
            if (!showCharts)
                MultiChartDisplayPanel.setHiddenMode(true);
            showCharts = true;
        }

        DeltaMassArgumentDefinition.DeltaMassWithType massMatchDeltaMassWithType =
            getDeltaMassArgumentValue("massmatchdeltamass");

        amtDatabaseMatcher = new AmtDatabaseMatcher();
        amtDatabaseMatcher.setMassMatchDeltaMass(
                massMatchDeltaMassWithType.getDeltaMass());
        amtDatabaseMatcher.setMassMatchDeltaMassType(
                massMatchDeltaMassWithType.getDeltaMassType());
        amtDatabaseMatcher.setRealMatchDeltaMass(looseDeltaMass);
        amtDatabaseMatcher.setRealMatchDeltaMassType(deltaMassType);
        amtDatabaseMatcher.setRealMatchDeltaElution(looseDeltaElution);
        amtDatabaseMatcher.setBuildCharts(showCharts);
        amtDatabaseMatcher.setMaxRegressionLeverageNumerator(
                getDoubleArgumentValue("maxregressionleverage"));
        amtDatabaseMatcher.setMaxRegressionStudentizedResidual(
                getDoubleArgumentValue("maxregressionstudres"));

        amtDatabaseMatcher.setMinMatchProbabilityToKeep(getFloatArgumentValue("minmatchprob"));
        amtDatabaseMatcher.setMaxMatchFDRToKeep(getFloatArgumentValue("maxmatchfdr"));

        amtDatabaseMatcher.setMaxSecondBestProbability(getFloatArgumentValue("maxsecondbestprob"));
        amtDatabaseMatcher.setMinSecondBestProbabilityDifference(getFloatArgumentValue("minsecondbestprobdiff"));
        amtDatabaseMatcher.setMinEMIterations(getIntegerArgumentValue("minemiterations"));
        amtDatabaseMatcher.setMaxEMIterations(getIntegerArgumentValue("maxemiterations"));
        amtDatabaseMatcher.setUseMs1TimesForAlignment(getBooleanArgumentValue("usems1foralignment"));
        amtDatabaseMatcher.setMs1Ms2MassTolerancePPM(getFloatArgumentValue("deltamassms1ms2ppm"));
        amtDatabaseMatcher.setMs1Ms2TimeToleranceSeconds(getFloatArgumentValue("deltatimems1ms2"));
        amtDatabaseMatcher.setNonlinearMappingPolynomialDegree(nonlinearMappingPolynomialDegree);
        amtDatabaseMatcher.setMaxRProbAssignmentMillis(rTimeoutMilliseconds);

        amtDatabaseMatcher.setDecoyMatch(dummyMatch);

        if (amtDatabaseStructure != null)
            amtDatabaseMatcher.defineAMTDBStructure(amtDatabaseStructure);


        AmtDatabaseFeatureSetGenerator featureGenerator =
                new AmtDatabaseFeatureSetGenerator(amtDB);
        amtDBBaseFeatures = featureGenerator.createFeaturesForModifications(
                ms2ModificationsForMatching);
//        if (dummyMatch)
//            for (Feature feature : amtDBBaseFeatures)
//                feature.setMass(feature.getMass() + 10);
    }


    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            ApplicationContext.infoMessage("AMT Database: " + amtDB.toString());

            switch(mode)
            {
                case  MODE_MATCH_AGAINST_MS1:
                    if (embeddedMs2FeatureSet != null)
                        ApplicationContext.setMessage("Using MS2 FeatureSet with " +
                                embeddedMs2FeatureSet.getFeatures().length + " features");
                    File outputFile = outFile;
                    if (outputFile == null && outDir != null)
                    {
                        outputFile = new File(outDir.getAbsolutePath() + File.separatorChar +
                                ms1FeatureSet.getSourceFile().getName().substring(0,
                                        ms1FeatureSet.getSourceFile().getName().length() - 3) +
                                "matched.tsv");
                    }
                    try
                    {
                        if (calcFDR)
                        {
                            amtDatabaseMatcher.matchWithFDRCalc(
                                    amtDB, amtDBBaseFeatures,
                                    ms1FeatureSet,
                                    embeddedMs2FeatureSet,
                                    (float) minMS2PeptideProphet,
                                    ms2ModificationsForMatching,
                                    outputFile,
                                    removeFractions, minRunsToKeep, maxRunsToKeep,
                                    calibrateMassesUsingMatches,
                                    showCharts);
                        }
                        else
                        {
                            amtDatabaseMatcher.matchAgainstMs1(
                                    amtDB, amtDBBaseFeatures,
                                    ms1FeatureSet,
                                    embeddedMs2FeatureSet,
                                    (float) minMS2PeptideProphet,
                                    ms2ModificationsForMatching,
                                    outputFile,
                                    removeFractions, minRunsToKeep, maxRunsToKeep,
                                    calibrateMassesUsingMatches,
                                    showCharts);
                        }
                        ApplicationContext.setMessage("Quantile Correlation Scores: Mass=" +
                                amtDatabaseMatcher.getProbabilityAssigner().getQuantileCorrX() + ", NRT=" +
                                amtDatabaseMatcher.getProbabilityAssigner().getQuantileCorrY());
                        ApplicationContext.setMessage("Quantile Regression Beta Scores: Mass=" +
                                amtDatabaseMatcher.getProbabilityAssigner().getQuantileBetaX() + ", NRT=" +
                                amtDatabaseMatcher.getProbabilityAssigner().getQuantileBetaY());
                        int numIterations = amtDatabaseMatcher.getProbabilityAssigner().getNumIterations();
                        if (amtDatabaseMatcher.getProbabilityAssigner().isConverged())
                            ApplicationContext.setMessage("EM algorithm converged after " + numIterations +
                                    " iterations.");
                        else
                            ApplicationContext.setMessage("EM algorithm did not converge in " + numIterations +
                                    " iterations!");
                        ApplicationContext.setMessage("Done.");
//                        if (showCharts)
//                            amtDatabaseMatcher.createMassTimeErrorPlots(amtDatabaseMatcher.featureMatchingResult);

                        if (saveChartsDir != null)
                        {
                            MultiChartDisplayPanel.getSingletonInstance().saveAllChartsToFiles(saveChartsDir);
                            ApplicationContext.infoMessage("Saved all charts to directory " +
                                    saveChartsDir.getAbsolutePath());
                        }
                    }
                    catch(RuntimeException e)
                    {
                        throw new CommandLineModuleExecutionException(e);
                    }
                    break;
                case  MODE_MATCH_AGAINST_MS1_DIR:
                    //sequentially match against each MS1 file
                    int totalNumPeptidesMatched = 0;

                    Map<File,Exception> filesThatFailed = new HashMap<File,Exception>();
                    int numAttemptedFiles = 0;

                    int numFilesTried=0;
                    List<Float> ksScoreXList = new ArrayList<Float>();
                    List<Float> ksScoreYList = new ArrayList<Float>();

                    List<Float> quantCorrXList = new ArrayList<Float>();
                    List<Float> quantCorrYList = new ArrayList<Float>();

                    List<Float> quantBetaXList = new ArrayList<Float>();
                    List<Float> quantBetaYList = new ArrayList<Float>();

                    List<Float> emIterationsList = new ArrayList<Float>();


                    for (File ms1File : inDir.listFiles())
                    {
                        if (showCharts)
                        {
                            MultiChartDisplayPanel.destroySingleton();
                        }
                        if (!ms1File.getName().endsWith(".tsv"))
                        {
                            ApplicationContext.infoMessage("Skipping non-tsv file " + ms1File.getName());
                            continue;
                        }
                        ApplicationContext.infoMessage("Processing file " + ms1File.getName());
                        numFilesTried++;
                        numAttemptedFiles++;
                        FeatureSet thisMs1FeatureSet = new FeatureSet(ms1File);

                        FeatureSet currentEmbeddedMs2FeatureSet  = null;
                        if (ms2Dir != null)
                        {
                            currentEmbeddedMs2FeatureSet = findCorrespondingFeatureSet(ms1File, ms2Dir);


                            if (populateMs2Times)
                            {
                                File mzXmlFile =
                                        findCorrespondingMzXmlFile(currentEmbeddedMs2FeatureSet.getSourceFile());
                                run = MSRun.load(mzXmlFile.getAbsolutePath());
                                currentEmbeddedMs2FeatureSet.populateTimesForMS2Features(run);
                            }

                        }

                        File currentOutputFile = null;
                        if (outDir != null)
                            currentOutputFile = new File(outDir.getAbsolutePath() + File.separatorChar +
                                    ms1File.getName().substring(0, ms1File.getName().length() - 3) +
                                    "matched.tsv");
                        try
                        {
                            amtDatabaseMatcher.matchAgainstMs1(
                                    amtDB, amtDBBaseFeatures,
                                    thisMs1FeatureSet,
                                    currentEmbeddedMs2FeatureSet,
                                    (float) minMS2PeptideProphet,
                                    ms2ModificationsForMatching,
                                    currentOutputFile,
                                    removeFractions, minRunsToKeep, maxRunsToKeep,
                                    calibrateMassesUsingMatches,
                                    showCharts);
                            ksScoreXList.add(amtDatabaseMatcher.getProbabilityAssigner().getKsScoreX());
                            ksScoreYList.add(amtDatabaseMatcher.getProbabilityAssigner().getKsScoreY());
                            quantCorrXList.add(amtDatabaseMatcher.getProbabilityAssigner().getQuantileCorrX());
                            quantCorrYList.add(amtDatabaseMatcher.getProbabilityAssigner().getQuantileCorrY());
                            quantBetaXList.add(amtDatabaseMatcher.getProbabilityAssigner().getQuantileBetaX());
                            quantBetaYList.add(amtDatabaseMatcher.getProbabilityAssigner().getQuantileBetaY());

                            int numIterations = amtDatabaseMatcher.getProbabilityAssigner().getNumIterations();
                            emIterationsList.add((float) numIterations);

                            ApplicationContext.infoMessage("EM converged after " + numIterations + " iterations.");
                            ApplicationContext.infoMessage("Completed matching file " + ms1File.getAbsolutePath());
                        }
                        catch (Exception e)
                        {
                            filesThatFailed.put(thisMs1FeatureSet.getSourceFile(), e);
                            ApplicationContext.infoMessage("\nWARNING: File " + thisMs1FeatureSet.getSourceFile().getName() +
                                    " Failed matching\n");
                            e.printStackTrace(System.err);
                        }
                        totalNumPeptidesMatched += thisRunMatchedPeptides.size();
                        ApplicationContext.infoMessage("Running count: " + (numFilesTried - filesThatFailed.size()) +
                                " out of " + numFilesTried + " files so far succeeded");

                        if (saveChartsDir != null)
                        {
                            MultiChartDisplayPanel.getSingletonInstance().saveAllChartsToFiles(saveChartsDir);
                            ApplicationContext.infoMessage("Saved all charts to directory " +
                                    saveChartsDir.getAbsolutePath());
                        }

                    }

                    if (filesThatFailed.size() > 0)
                    {
                        ApplicationContext.infoMessage("\n" + filesThatFailed.size() +
                                " files (out of " + numAttemptedFiles + ") failed matching:");
                        for (File file : filesThatFailed.keySet())
                        {
                            ApplicationContext.infoMessage("\t" + file.getName() + "\t" +
                                    filesThatFailed.get(file).getMessage());
                        }
                        ApplicationContext.infoMessage("");
                    }

                    ApplicationContext.infoMessage("KS score distribution: ");
                    float minKsX = Float.MAX_VALUE;
                    float maxKsX = Float.MIN_VALUE;
                    for (Float ks : ksScoreXList)
                    {
                        if (ks < minKsX)
                            minKsX = ks;
                        if (ks > maxKsX)
                            maxKsX = ks;
                    }
                    ApplicationContext.infoMessage("\tX: min: " + minKsX + ", max: " + maxKsX + ", median: " +
                            BasicStatistics.median(ksScoreXList) + ", mean: " + BasicStatistics.mean(ksScoreXList));

                    float minKsY = Float.MAX_VALUE;
                    float maxKsY = Float.MIN_VALUE;
                    for (Float ks : ksScoreYList)
                    {
                        if (ks < minKsY)
                            minKsY = ks;
                        if (ks > maxKsY)
                            maxKsY = ks;
                    }
                    List<Float> allKsScoreList = new ArrayList<Float>();
                    allKsScoreList.addAll(ksScoreYList);
                    allKsScoreList.addAll(ksScoreXList);
                    ApplicationContext.infoMessage("\tY: min: " + minKsY + ", max: " + maxKsY + ", median: " +
                            BasicStatistics.median(ksScoreYList) + ", mean: " + BasicStatistics.mean(ksScoreYList));
                    ApplicationContext.infoMessage("\tOverall median: " + BasicStatistics.median(allKsScoreList) +
                        ", mean: " + BasicStatistics.mean(allKsScoreList));


                    ApplicationContext.infoMessage("Quantile Correlation distribution: ");
                    float minQuantCorrX = Float.MAX_VALUE;
                    float maxQuantCorrX = Float.MIN_VALUE;
                    for (Float quantCorr : quantCorrXList)
                    {
                        if (quantCorr < minQuantCorrX)
                            minQuantCorrX = quantCorr;
                        if (quantCorr > maxQuantCorrX)
                            maxQuantCorrX = quantCorr;
                    }
                    ApplicationContext.infoMessage("\tX: min: " + minQuantCorrX + ", max: " + maxQuantCorrX + ", median: " +
                            BasicStatistics.median(quantCorrXList) + ", mean: " + BasicStatistics.mean(quantCorrXList));

                    float minQuantCorrY = Float.MAX_VALUE;
                    float maxQuantCorrY = Float.MIN_VALUE;
                    for (Float quantCorr : quantCorrYList)
                    {
                        if (quantCorr < minQuantCorrY)
                            minQuantCorrY = quantCorr;
                        if (quantCorr > maxQuantCorrY)
                            maxQuantCorrY = quantCorr;
                    }
                    ApplicationContext.infoMessage("\tY: min: " + minQuantCorrY + ", max: " + maxQuantCorrY + ", median: " +
                            BasicStatistics.median(quantCorrYList) + ", mean: " + BasicStatistics.mean(quantCorrYList));

                    List<Float> allQuantCorrScoreList = new ArrayList<Float>();
                    allQuantCorrScoreList.addAll(quantCorrXList);
                    allQuantCorrScoreList.addAll(quantCorrYList);
                    ApplicationContext.infoMessage("\tOverall median: " + BasicStatistics.median(allQuantCorrScoreList) +
                        ", mean: " + BasicStatistics.mean(allQuantCorrScoreList));


                    ApplicationContext.infoMessage("Quantile Beta distribution: ");
                    float minQuantBetaX = Float.MAX_VALUE;
                    float maxQuantBetaX = Float.MIN_VALUE;
                    for (Float quantBeta : quantBetaXList)
                    {
                        if (quantBeta < minQuantBetaX)
                            minQuantBetaX = quantBeta;
                        if (quantBeta > maxQuantBetaX)
                            maxQuantBetaX = quantBeta;
                    }
                    ApplicationContext.infoMessage("\tX: min: " + minQuantBetaX + ", max: " + maxQuantBetaX + ", median: " +
                            BasicStatistics.median(quantBetaXList) + ", mean: " + BasicStatistics.mean(quantBetaXList));

                    float minQuantBetaY = Float.MAX_VALUE;
                    float maxQuantBetaY = Float.MIN_VALUE;
                    for (Float quantBeta : quantBetaYList)
                    {
                        if (quantBeta < minQuantBetaY)
                            minQuantBetaY = quantBeta;
                        if (quantBeta > maxQuantBetaY)
                            maxQuantBetaY = quantBeta;
                    }
                    ApplicationContext.infoMessage("\tY: min: " + minQuantBetaY + ", max: " + maxQuantBetaY + ", median: " +
                            BasicStatistics.median(quantBetaYList) + ", mean: " + BasicStatistics.mean(quantBetaYList));

                    List<Float> allQuantBetaScoreList = new ArrayList<Float>();
                    allQuantBetaScoreList.addAll(quantBetaXList);
                    allQuantBetaScoreList.addAll(quantBetaYList);
                    ApplicationContext.infoMessage("\tOverall median: " + BasicStatistics.median(allQuantBetaScoreList) +
                        ", mean: " + BasicStatistics.mean(allQuantBetaScoreList));                    

                    float minIterations = Float.MAX_VALUE;
                    float maxIterations = Float.MIN_VALUE;
                    for (Float numIterations : emIterationsList)
                    {
                        if (numIterations < minIterations)
                            minIterations = numIterations;
                        if (numIterations > maxIterations)
                            maxIterations = numIterations;
                    }
                    ApplicationContext.infoMessage("Iterations: min: " + minIterations + ", max: " + maxIterations +
                            ", median: " +
                            BasicStatistics.median(emIterationsList) + ", mean: " + BasicStatistics.mean(emIterationsList));
                    int numMinIterations = 0;
                    for (float numIt : emIterationsList)
                        if (numIt == AmtMatchProbabilityAssigner.DEFAULT_MIN_EM_ITERATIONS) numMinIterations++;
                    if (!emIterationsList.isEmpty())
                    	ApplicationContext.infoMessage("Number of runs converged after minimum iterations: " +
                            numMinIterations + "(" + (numMinIterations * 100 / emIterationsList.size()) + "%)");

                    ApplicationContext.infoMessage("Done.");
//PanelWithHistogram pspX = new PanelWithHistogram(quantCorrXList, "NRT corr scores");
//pspX.displayInTab();
//PanelWithHistogram pspY = new PanelWithHistogram(quantCorrYList, "Mass corr scores");
//pspY.displayInTab();
                    break;
            }


        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

                }

    protected FeatureSet findCorrespondingFeatureSet(File ms1File, File directory)
            throws IOException
    {
        String pepXmlFilename =
                (ms1File.getName().substring(0,
                        ms1File.getName().indexOf(".")) + ".pep.xml");
        //".mod.pep.xml" is the default extension output by postprocesspepxml
        String modPepXmlFilename =
                (ms1File.getName().substring(0,
                        ms1File.getName().indexOf(".")) + ".mod.pep.xml");
        String tsvFilename =
                (ms1File.getName().substring(0,
                        ms1File.getName().indexOf(".")) + ".tsv");
        String filteredTsvFilename =
                (ms1File.getName().substring(0,
                        ms1File.getName().indexOf(".")) + ".filtered.tsv");

        File resultFile = null;
        for (String potentialMs2Filename : directory.list())
        {
            if (potentialMs2Filename.equalsIgnoreCase(pepXmlFilename) ||
                    potentialMs2Filename.equalsIgnoreCase(tsvFilename) ||
                    potentialMs2Filename.equalsIgnoreCase(filteredTsvFilename) ||
                    potentialMs2Filename.equalsIgnoreCase(modPepXmlFilename))
            {
                resultFile = new File(directory.getAbsolutePath() + File.separatorChar +
                        potentialMs2Filename);
                break;
            }
        }
        if (resultFile == null)
            throw new IOException("No corresponding feature file for file " + ms1File.getAbsolutePath());

        FeatureSet result = null;
        try
        {
            result = new FeatureSet(resultFile);
            ApplicationContext.setMessage("Located feature file " +
                    resultFile.getAbsolutePath() +
                    " with " + result.getFeatures().length + " features");
        }
        catch (Exception e)
        {
            throw new IOException("Failed to load feature file " + resultFile.getAbsolutePath() +
                    ": " + e.getMessage());
        }

        return result;
    }

    protected File findCorrespondingMzXmlFile(File ms2File)
            throws IOException
    {
        String ms2Filename = ms2File.getName();
        String mzXmlFileName =
                (ms2Filename.substring(0, ms2Filename.indexOf(".")) +
                        ".mzXML");
        File result = null;
        for (File potentialMzXmlFile : mzXmlDir.listFiles())
        {
            String potentialMzXmlFilename = potentialMzXmlFile.getName();
            if (potentialMzXmlFilename.equalsIgnoreCase(mzXmlFileName))
            {
                result = potentialMzXmlFile;
                ApplicationContext.setMessage("Located mzXML file " +
                        potentialMzXmlFile.getAbsolutePath());
            }
        }
        if (result == null)
            throw new IOException("No corresponding mzXML file for MS2 file " + ms2File.getAbsolutePath());
        return result;
    }


            //unnecessary?
    protected void annotateFractionConcordanceForMatches(
            AmtDatabase amtDatabase,
            Feature[] ms2Features,
            FeatureSetMatcher.FeatureMatchingResult matchingResult)
    {
        float[] percentMatchedPerRun = new float[amtDatabase.numRuns()];
        AmtRunEntry[] runEntries = amtDatabase.getRuns();

        final Map<AmtRunEntry, Integer> runPeptideMatchNumberMap =
                new HashMap<AmtRunEntry, Integer>();
        final Map<AmtRunEntry, Float> runPeptideMatchPercentMap =
                new HashMap<AmtRunEntry, Float>();

        Set<String> ms2Peptides = AmtDatabaseMatcher.createPeptideSetFromFeatures(ms2Features);

        Map<String, Float> matchingPeptideMaxRunConcordancePercentMap =
                new HashMap<String, Float>();

        Set<Feature> matchedAmtFeatures = matchingResult.getSlaveSetFeatures();
        Set<String> matchedAmtPeptides = new HashSet<String>();
        for (Feature matchedAmtFeature : matchedAmtFeatures)
            matchedAmtPeptides.add(MS2ExtraInfoDef.getFirstPeptide(matchedAmtFeature));

        for (int i=0; i<runEntries.length; i++)
        {
            if (i % (Math.max(runEntries.length/10,1)) == 0)
            {
                ApplicationContext.setMessage(Rounder.round(((double) i * 100.0 / (double) runEntries.length),0) +
                        "% done evaluating runs");
            }
            AmtRunEntry runEntry = runEntries[i];

            Set<String> matchedPeptidesInThisRun = new HashSet<String>();

            int numPeptidesThisRun = 0;
            int numPeptidesInCommon = 0;
            for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
            {
                AmtPeptideEntry.AmtPeptideObservation obs =
                        peptideEntry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    numPeptidesThisRun++;
                    String peptideSequence = peptideEntry.getPeptideSequence();                                       
                    if (ms2Peptides.contains(peptideSequence))
                        numPeptidesInCommon++;
                    if (matchedAmtPeptides.contains(peptideSequence))
                        matchedPeptidesInThisRun.add(peptideSequence);
                }
            }

            float peptideMatchesPercent =
                    (((float) numPeptidesInCommon) /
                            ((float) numPeptidesThisRun) * 100);
//            System.err.println("Matches this run: " + numMassMatchedFeatures + " (" + massMatchesPercent + "%)");
            percentMatchedPerRun[i] = peptideMatchesPercent;

            runPeptideMatchPercentMap.put(runEntry, peptideMatchesPercent);
            runPeptideMatchNumberMap.put(runEntry, numPeptidesInCommon);

            for (String peptideSequence : matchedPeptidesInThisRun)
            {
                Float oldMaxConcordancePercent =
                        matchingPeptideMaxRunConcordancePercentMap.get(peptideSequence);
                if (oldMaxConcordancePercent == null)
                    oldMaxConcordancePercent = 0f;
                matchingPeptideMaxRunConcordancePercentMap.put(peptideSequence, Math.max(oldMaxConcordancePercent, peptideMatchesPercent));
            }
        }

        if (showCharts && amtDatabaseStructure != null)
        {
            double[] peptideMatchesPerRun = new double[amtDatabase.numRuns()];
            for (int i=0; i<peptideMatchesPerRun.length; i++)
                peptideMatchesPerRun[i] = runPeptideMatchNumberMap.get(amtDatabase.getRunBySequence(i+1));

            Fractionation2DUtilities.showHeatMapChart(
                    amtDatabaseStructure,
                    peptideMatchesPerRun, "AMT Database Run Peptides In Common with This Run",
                    true);

            double[] peptidePercentMatchesPerRun = new double[amtDatabase.numRuns()];
            for (int i=0; i<peptidePercentMatchesPerRun.length; i++)
                peptidePercentMatchesPerRun[i] = runPeptideMatchPercentMap.get(amtDatabase.getRunBySequence(i+1));

            Fractionation2DUtilities.showHeatMapChart(
                    amtDatabaseStructure,
                    peptidePercentMatchesPerRun, "AMT Database Run Peptide % In Common with This Run",
                    true);
        }

        if (showCharts)
        {
            List<Double> concordancePercentValues = new ArrayList<Double>();
            List<Double> deltaMassValues = new ArrayList<Double>();
            List<Double> deltaHValues = new ArrayList<Double>();
            List<Double> deltaMassTimesDeltaHValues = new ArrayList<Double>();


            for (Feature masterSetFeature : matchingResult.getMasterSetFeatures())
            {
                List<Feature> matchesThisFeature = matchingResult.get(masterSetFeature);
                if (matchesThisFeature.size() == 1)
                {
                    Feature amtMatchedFeature = matchesThisFeature.get(0);
                    concordancePercentValues.add((double)matchingPeptideMaxRunConcordancePercentMap.get(MS2ExtraInfoDef.getFirstPeptide(amtMatchedFeature)));

                    double deltaMass =  (double)(masterSetFeature.getMass() - amtMatchedFeature.getMass()) *
                                         (1000000 / masterSetFeature.getMass());
                    double deltaH = (AmtExtraInfoDef.getObservedHydrophobicity(masterSetFeature) -
                            (AmtExtraInfoDef.getObservedHydrophobicity(amtMatchedFeature)));
                    deltaMassValues.add(deltaMass);
                    deltaHValues.add(deltaH);
                    deltaMassTimesDeltaHValues.add(deltaMass * deltaH);
                }
            }
            ScatterPlotDialog spdMass =
                    new ScatterPlotDialog(deltaMassValues, concordancePercentValues, "Delta Mass vs. Concordance");
            spdMass.setVisible(true);
            ScatterPlotDialog spdH = new ScatterPlotDialog(deltaHValues, concordancePercentValues, "Delta H vs. Concordance");
            spdH.setVisible(true);

                        ScatterPlotDialog spdMassH = new ScatterPlotDialog(deltaMassTimesDeltaHValues, concordancePercentValues, "Delta Mass*H vs. Concordance");
            spdMassH.setVisible(true);

            float[] numMatchesWithPercentConcordance = new float[100];
            float[] percentValues = new float[100];
            for (int i=0; i<100; i++)
            {
                percentValues[i] = i;
                int numMatchesThisPercent = 0;
                for (Double percent : concordancePercentValues)
                {
                    if (percent < i)
                        numMatchesThisPercent++;
                }
                numMatchesWithPercentConcordance[i] = numMatchesThisPercent;
            }
            PanelWithLineChart pwlc = new PanelWithLineChart(percentValues, numMatchesWithPercentConcordance, "num matches with percent, vs. percent");
            ChartDialog cd = new ChartDialog(pwlc);
            cd.setVisible(true);
        }


    }


    public FeatureSet getMs1FeatureSet()
    {
        return ms1FeatureSet;
    }

    public void setMs1FeatureSet(FeatureSet ms1FeatureSet)
    {
        this.ms1FeatureSet = ms1FeatureSet;
    }

    public AmtDatabaseMatcher getAmtDatabaseMatcher()
    {
        return amtDatabaseMatcher;
    }

}
