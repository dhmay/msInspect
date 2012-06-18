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

import org.fhcrc.cpl.viewer.commandline.modules.FeatureSelectionParamsCommandLineModule;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FastaFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.IOException;


/**
 * Command line module for creating AMT databases, either from pepXml (or other feature)
 * files, or from other AMT databases, or from random peptides pulled out of a
 * FASTA file
 */
public class AmtDatabaseCreatorCommandLineModule extends
        FeatureSelectionParamsCommandLineModule
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AmtDatabaseCreatorCommandLineModule.class);


    protected File mzXmlDir = null;
    protected File ms2FeaturesDir = null;
    protected File ms1FeaturesDir = null;

    File pepXmlOrMS2FeatureFile = null;
    FeatureSet ms1FeatureSet = null;
    File mzXmlFile = null;
    File[] inputFiles = null;
    MSRun run = null;

    protected File outFile = null;

    protected int scanOrTimeMode = ProteomicsRegressionUtilities.DEFAULT_REGRESSION_MODE;

    //TODO: get rid of this?
    protected boolean robustRegression = false;

    protected double maxStudentizedResidualForRegression;
    protected double maxStudentizedResidualForInclusion;

    protected int mode=-1;

    protected static final int CREATE_AMTXML_FROM_DIRECTORIES_MODE=0;
    protected static final int CREATE_AMTXML_FROM_MS2_FEATURES_MODE=1;
    protected static final int CREATE_AMTXML_FROM_MULTIPLE_AMT_XMLS_MODE=2;
    protected static final int CREATE_AMTXML_FROM_RANDOM_PEPTIDES_MODE=3;
    protected static final int CREATE_AMTXML_FROM_PROTEIN_TRYPTIC_PEPTIDES_MODE=4;

    //min and max peptide length for tryptic digest
    protected int minPeptideLength = 6;
    protected int maxPeptideLength = 25;


    protected final static String[] modeStrings = {"directories",
             "ms2featurefile","amtxmls","randompeptides","proteintryptic"};

    protected final static String[] modeExplanations =
            {
                    "supply directories for MS2 and mzXML files",
                    "create an AMT database from a single MS2 and mzXML file",
                    "combine multiple existing AMT databases",
                    "create a database of random peptides from a FASTA file",
                    "create a database of tryptic peptides from a list of proteins in a FASTA file",
            };

    protected boolean align = false;

    //for building random peptide databases
    protected int numRandomPeptides = 50000;
    protected Protein[] proteinsFromFasta;

    protected List<Protein> proteinsToDigest;
    protected int maxMissedCleavages = 0;

    protected boolean populateMs2Times = true;

    protected AmtDatabaseBuilder databaseBuilder;




    public AmtDatabaseCreatorCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        super.init();
        mCommandName = "createamt";
        mShortDescription = "Create an AMT database.";
        mHelpMessage = "Create an AMT database to store peptide observations from several, or many, LC-MS/MS runs. " +
                "Many of the arguments for this command have to do with filtering the peptides used in the database.";

        CommandLineArgumentDefinition[] childBasicArgDefs =
        {
                new EnumeratedValuesArgumentDefinition("mode",true,
                       modeStrings, modeExplanations),
                new FileToWriteArgumentDefinition("out", false, "output file"),
                new DirectoryToReadArgumentDefinition("mzxmldir", false,
                        "Directory of mzXML files (for 'directories' mode), only necessary if retention times are " +
                                "not populated in the MS2 feature file"),
                new DirectoryToReadArgumentDefinition("ms2dir", false,
                        "Directory of MS2 feature files (for 'directories' mode)"),
                new FileToReadArgumentDefinition("ms2features", false,
                        "Input MS2 feature file (for 'ms2featurefile' mode)"),
                new FileToReadArgumentDefinition("mzxml", false,
                        "Input mzXml file (for 'ms2featurefile' mode), only necessary if retention times are not " +
                                "populated in the MS2 feature file"),
                createUnnamedSeriesFileArgumentDefinition(false,
                        "Input file (for 'ms2features' mode)"),
                new FastaFileArgumentDefinition("fasta", false,
                        "FASTA file to pull random peptides from ('randompeptides' mode only"),
                new DirectoryToReadArgumentDefinition("ms1dir", false,
                        "Directory of MS1 feature files (for 'directories' mode)"),
                new FileToReadArgumentDefinition("ms1features", false,
                        "Input MS1 feature file (for 'ms1featurefile' mode)"),
        };
        addArgumentDefinitions(childBasicArgDefs);


        CommandLineArgumentDefinition[] childAdvancedArgDefs =
        {
                new EnumeratedValuesArgumentDefinition("scanortimemode",false,
                        "Use scans or times from features (default 'time')",
                        new String[]{"scan","time"}),
                new BooleanArgumentDefinition("align", false,
                        "use nonlinear alignment when mapping time to hydrophobicity.  This is not necessarily " +
                                "recommended, as the manageamt command has a mode ('alignallruns') for nonlinearly " +
                                "aligning all runs to a single scale that is much more effective.", align),
                new IntegerArgumentDefinition("numpeptides", false,
                        "Number of random peptides to use in database creation ('randompeptides' mode only)",
                        numRandomPeptides),
                new DecimalArgumentDefinition("maxsrforregression", false,
                        "maximum studentized residual for use in regression calculation for transforming RT to NRT",
                        AmtDatabaseBuilder.DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_REGRESSION),
                new DecimalArgumentDefinition("maxsrforinclusion", false,
                        "maximum studentized residual for inclusion in database.  Any observation with a higher " +
                                "studentized residual, based on the RT->NRT regression, will be excluded",
                        AmtDatabaseBuilder.DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_INCLUSION),
                new DecimalArgumentDefinition("deltamassppm", false,
                        "Mass tolerance for MS1 feature match with MS2 in order to retrieve MS1 feature times",
                        AmtDatabaseBuilder.DEFAULT_MS1_MS2_MASS_TOLERANCE_PPM),
                new DecimalArgumentDefinition("deltatime", false,
                        "Time tolerance (in seconds) for MS1 feature match with MS2 in order to retrieve MS1 feature " +
                                "times", AmtDatabaseBuilder.DEFAULT_MS1_MS2_TIME_TOLERANCE_SECONDS),
                new BooleanArgumentDefinition("ignoreunknownmods", false,
                        "Ignore modifications on individual peptide IDs that are not declared at the top of the " +
                        "pepXML file? (otherwise, fail when encountering these modifications)",
                        AmtDatabaseBuilder.DEFAULT_IGNORE_UNKNOWN_MODIFICATIONS),
                new IntegerArgumentDefinition("maxmissedcleavages", false,
                        "Maximum number of missed cleavages, for generating peptides from a list of proteins",
                        maxMissedCleavages),
                new StringArgumentDefinition("proteins", false,
                        "Comma-separated list of protein identifiers"),
        };
        addArgumentDefinitions(childAdvancedArgDefs, true);       
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        //feature-filtering parameters
        super.assignArgumentValues();

        boolean usingMs1Times = false;

        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        maxStudentizedResidualForRegression =
                getDoubleArgumentValue("maxsrforregression");
        maxStudentizedResidualForInclusion =
                getDoubleArgumentValue("maxsrforinclusion");

        if (hasArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
        {
            inputFiles = getUnnamedSeriesFileArgumentValues();
        }

        populateMs2Times = hasArgumentValue("mzxml") || hasArgumentValue("mzxmldir");



        switch (mode)
        {
            case CREATE_AMTXML_FROM_DIRECTORIES_MODE:
            {
                assertArgumentPresent("ms2dir");
                assertArgumentAbsent("ms1features");
                assertArgumentAbsent("ms2features");
                if (populateMs2Times)
                    assertArgumentPresent("mzxmldir");
                else
                    assertArgumentAbsent("mzxmldir");
                mzXmlDir = getFileArgumentValue("mzxmldir");
                ms2FeaturesDir = getFileArgumentValue("ms2dir");
                ms1FeaturesDir = getFileArgumentValue("ms1dir");

                break;
            }
            case CREATE_AMTXML_FROM_MS2_FEATURES_MODE:
            {
                assertArgumentPresent("ms2features");
                assertArgumentAbsent("ms2dir");

                pepXmlOrMS2FeatureFile = getFileArgumentValue("ms2features");
                if (hasArgumentValue("ms1features"))
                {
                    try
                    {
                        ms1FeatureSet = new FeatureSet(getFileArgumentValue("ms1features"));
                    }
                    catch (Exception e)
                    {
                        throw new ArgumentValidationException("Failed to load MS1 features from file " +
                            getFileArgumentValue("ms1features").getAbsolutePath(),e);
                    }
                }

                ms1FeaturesDir = getFileArgumentValue("ms1dir");
                if (ms1FeatureSet != null && ms1FeaturesDir != null)
                    throw new ArgumentValidationException("Can't specify both ms1dir and ms1features");


                mzXmlFile = getFileArgumentValue("mzxml");
                if (mzXmlFile != null)
                {
                    try
                    {
                        run = MSRun.load(mzXmlFile.getAbsolutePath());
                    }
                    catch (IOException e)
                    {
                        throw new ArgumentValidationException("Failed to load mzXml file", e);
                    }
                }

                break;
            }
            case CREATE_AMTXML_FROM_MULTIPLE_AMT_XMLS_MODE:
            {
                assertArgumentPresent(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
                break;
            }
            case CREATE_AMTXML_FROM_RANDOM_PEPTIDES_MODE:
            {
                assertArgumentPresent("fasta");
                numRandomPeptides = getIntegerArgumentValue("numpeptides");
                proteinsFromFasta = getFastaFileArgumentValue("fasta");
                break;
            }
            case CREATE_AMTXML_FROM_PROTEIN_TRYPTIC_PEPTIDES_MODE:
            {
                assertArgumentPresent("fasta");                
                assertArgumentPresent("proteins");
                maxMissedCleavages = getIntegerArgumentValue("maxmissedcleavages");
                String proteinsString = getStringArgumentValue("proteins");
                String[] proteinNamesToDigest = proteinsString.split(",");

                //todo: this could be a lot more efficient.  Whatever.
                proteinsFromFasta = getFastaFileArgumentValue("fasta");
                proteinsToDigest = new ArrayList<Protein>();
                Set<String> proteinsFound = new HashSet<String>();
                for (Protein protein : proteinsFromFasta)
                {
                    for (String proteinName : proteinNamesToDigest)
                    {
                        if (protein.getLookup().equalsIgnoreCase(proteinName))
                        {
                            proteinsToDigest.add(protein);
                            proteinsFound.add(protein.getLookup());
                        }
                    }
                }

                for (String protein : proteinNamesToDigest)
                {
                    if (!proteinsFound.contains(protein))
                        ApplicationContext.infoMessage("WARNING! protein " + protein + " not found in FASTA");
                }
            }
        }

        if (ms1FeatureSet != null || ms1FeaturesDir != null)
        {
            ApplicationContext.infoMessage("Using MS1 features, not MS2 features, for retention times");
            usingMs1Times = true;
        }
        else
            ApplicationContext.infoMessage("WARNING: Using MS2 features (no MS1 features), for retention times.  " +
                    "If MS1 features are available, they are generally preferred.");

        databaseBuilder = new AmtDatabaseBuilder();
        databaseBuilder.setIgnoreUnknownModifications(getBooleanArgumentValue("ignoreunknownmods"));

        if (usingMs1Times)
        {
            databaseBuilder.setMs1Ms2MassTolerancePPM(getFloatArgumentValue("deltamassppm"));
            databaseBuilder.setMs1Ms2TimeToleranceSeconds(getFloatArgumentValue("deltatime"));
        }

        if (hasArgumentValue("scanortimemode"))
            if ("scan".equalsIgnoreCase(getStringArgumentValue("scanortimemode")))
                scanOrTimeMode = ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN;

        outFile = getFileArgumentValue("out");

        align = getBooleanArgumentValue("align");
        
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        AmtDatabase amtDB = null;


        switch (mode)
        {
            case CREATE_AMTXML_FROM_DIRECTORIES_MODE:
            {
                try
                {
                    amtDB = databaseBuilder.createAmtDBFromDirectories(
                                    ms2FeaturesDir, ms1FeaturesDir, mzXmlDir,
                                    scanOrTimeMode, featureSelector,
                                    robustRegression,
                                    maxStudentizedResidualForRegression,
                                    maxStudentizedResidualForInclusion,
                                    align);
                }
                catch (Exception e)
                {
                    _log.error("Failed to build AMT database: ", e);
                }
                break;
            }            
            case CREATE_AMTXML_FROM_MS2_FEATURES_MODE:
            {
                try
                {
                    try
                    {
                        amtDB = databaseBuilder.createAmtDBFromPepXml(
                            pepXmlOrMS2FeatureFile, ms1FeaturesDir,  run, scanOrTimeMode, featureSelector,
                            robustRegression,
                            maxStudentizedResidualForRegression,
                            maxStudentizedResidualForInclusion);
                        if (ms1FeatureSet != null)
                            ApplicationContext.infoMessage("WARNING!!! You supplied an MS1 feature file.  " +
                                    "It was NOT USED.  Please supply a DIRECTORY of MS1 files if using pepXML " +
                                    "to build the database");
                    }
                    catch (Exception e)
                    {
                        if (ms1FeaturesDir != null)
                              throw new CommandLineModuleExecutionException("ERROR.  You supplied a non-pepXML MS2 " +
                                      "file as well as a directory of MS1 files.  Please supply a single MS1 " +
                                      "file instead.");

                        //Failed.  Try again, treating the "pepxml" file as a feature file
                        FeatureSet featureSet = new FeatureSet(pepXmlOrMS2FeatureFile);
                        Pair<Integer,Integer> numFeaturesDummy = new Pair<Integer,Integer>(0,0);
                        if (populateMs2Times)
                            featureSet.populateTimesForMS2Features(run);
                        
                        amtDB = databaseBuilder.createAmtDatabaseForRun(featureSet, ms1FeatureSet,
                            scanOrTimeMode, robustRegression, numFeaturesDummy,
                            maxStudentizedResidualForRegression,
                            maxStudentizedResidualForInclusion);
                    }
                }
                catch (Exception e)
                {
                    _log.error("Failed to build AMT database: ", e);
                }
                break;
            }
            case CREATE_AMTXML_FROM_MULTIPLE_AMT_XMLS_MODE:
            {
                amtDB = new AmtDatabase();
                try
                {
                    for (File amtXmlFile : inputFiles)
                    {
                        ApplicationContext.infoMessage("Processing file " + amtXmlFile.getAbsolutePath() + "...");
                        AmtXmlReader amtXmlReader = new AmtXmlReader(amtXmlFile);
                        AmtDatabase additionalAmtDatabase = amtXmlReader.getDatabase();

                        amtDB.addObservationsFromAnotherDatabase(additionalAmtDatabase);
                        ApplicationContext.infoMessage("\tadded database from file " +
                                amtXmlFile.getName() + ": ");
                        ApplicationContext.infoMessage(additionalAmtDatabase.toString());
                    }
                    ApplicationContext.infoMessage("Combined database:");
                    ApplicationContext.infoMessage(amtDB.toString());

                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
                break;
            }
            case CREATE_AMTXML_FROM_RANDOM_PEPTIDES_MODE:
            {
                Set<String> randomPeptideStrings = new HashSet<String>();
                Set<Peptide> randomPeptides = new HashSet<Peptide>();
                List<Peptide>[] digestedProteins = new List[proteinsFromFasta.length];

                //the logic for choosing random peptides is right here.  Should ideally
                //be moved down into AmtDatabaseBuilder
                while (randomPeptides.size() < numRandomPeptides)
                {
                    int proteinIndex =
                            (int) Math.round(Math.random() * proteinsFromFasta.length);
                    //this skews things slightly, but, whatever
                    if (proteinIndex == proteinsFromFasta.length) proteinIndex--;

                    if (digestedProteins[proteinIndex] == null)
                    {
                        List<Peptide> peptideList =
                                ProteinMatcher.generatePeptidesFromProtein(
                                        proteinsFromFasta[proteinIndex], 0, 5);
                        digestedProteins[proteinIndex] = peptideList;
                    }

                    if (digestedProteins[proteinIndex].size() == 0)
                        continue;

                    int tries=0;
                    boolean chosen=false;
                    while (!chosen && tries < 20)
                    {
                        int peptideIndex =
                                (int) Math.round(Math.random() *
                                digestedProteins[proteinIndex].size());
                        if (peptideIndex == digestedProteins[proteinIndex].size())
                            peptideIndex--;
                        Peptide randomPeptide =
                                digestedProteins[proteinIndex].get(peptideIndex);
                        if (!randomPeptideStrings.contains(
                                new String(randomPeptide.getChars())))
                        {
                            chosen=true;
                            randomPeptides.add(randomPeptide);
                            randomPeptideStrings.add(new String(randomPeptide.getChars()));
                        }
                        tries++;
                    }
                }

                //now we've got numRandomPeptides peptides chosen
                amtDB = new AmtDatabase();
                AmtRunEntry dummyRun = new AmtRunEntry(new double[]{0,1},new MS2Modification[0]);
                amtDB.addRunEntry(dummyRun);
                for (Peptide randomPeptide : randomPeptides)
                    amtDB.addObservation(new String(randomPeptide.getChars()),
                            null, .95,
                            AmtUtilities.calculateNormalizedHydrophobicity(randomPeptide),
                            dummyRun, 1, 3000);
            }
            case CREATE_AMTXML_FROM_PROTEIN_TRYPTIC_PEPTIDES_MODE:
            {
                amtDB = new AmtDatabase();
                AmtRunEntry dummyRun = new AmtRunEntry(new double[]{0,1},new MS2Modification[0]);
                amtDB.addRunEntry(dummyRun);
                PeptideGenerator peptideGenerator = new PeptideGenerator();
                peptideGenerator.setMaxMissedCleavages(maxMissedCleavages);
                peptideGenerator.setMinResidues(minPeptideLength);
                peptideGenerator.setMaxResidues(maxPeptideLength);

                for (Protein protein : proteinsToDigest)
                {
                    Peptide[] peptides = peptideGenerator.digestProtein(protein);
                    for (Peptide peptide : peptides)
                    {
                        amtDB.addObservation(new String(peptide.getChars()),
                                null, .95,
                                AmtUtilities.calculateNormalizedHydrophobicity(peptide),
                                dummyRun, 1, 3000);
                    }
                }
            }
        }

        if (outFile != null && amtDB != null)
            writeAmtDatabase(amtDB, outFile);
    }



    /**
     *
     * @param amtDatabase
     * @param outAmtXmlFile
     */
    protected static void writeAmtDatabase(AmtDatabase amtDatabase,
                                           File outAmtXmlFile)
    {
        try
        {
            AmtXmlWriter amtXmlWriter = new AmtXmlWriter(amtDatabase);
            amtXmlWriter.write(outAmtXmlFile);
            ApplicationContext.infoMessage("Wrote " +
                    amtDatabase.numEntries() + " entries to amtxml file " +
                    outAmtXmlFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            ApplicationContext.infoMessage("Error writing amt file " +
                    outAmtXmlFile.getAbsolutePath());
        }
    }


}
