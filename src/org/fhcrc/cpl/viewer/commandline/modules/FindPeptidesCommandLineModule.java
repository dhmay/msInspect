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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFindingBroker;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategy;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.FloatRange;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;


/**
 * Command Module for feature finding.
 *
 * Not much of the actual work of feature finding is done directly in this class.  This class serves mainly to
 * parse user arguments and control the flow of work.  This design philosophy is used throughout Command Modules
 * in the msInspect platform
 */
public class FindPeptidesCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    //used for multiple levels of log messages
    protected static Logger _log = Logger.getLogger(FindPeptidesCommandLineModule.class);

    //Set of input mzXML files to process
    protected File[] mzXmlFiles = null;
    //Single output file (for only one mzXML file)
    protected File outFile = null;
    //Output directory -- if this is provided, an output filename will be chosen corresponding to each input
    //filename.  Required if more than one input file
    protected File outDir = null;
    //Scan to begin with
    protected int startScan = 1;
    //Number of scans to process
    protected int scanCount = Integer.MAX_VALUE;
    //Minimum and maximum m/z values for features
    protected float minMz = 0;
    protected float maxMz = Float.MAX_VALUE;
    //Size of the window of raw spectra around each feature to dump to the file
    protected int dumpWindowSize = 0;
    //If accurate mass adjustment is performed, the number of scans around each feature to check for
    //the most-intense peak in order to adjust the mass
    protected int accurateMassAdjustmentScans =
            FeatureFinder.DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS;
    //The feature strategy to use
    protected Class<? extends FeatureStrategy> featureStrategyClass =
            FeatureFinder.DEFAULT_FEATURE_FINDING_CLASS;
    //Should we perform ridge-walking in smoothed space?
    protected boolean peakRidgeWalkSmoothed = false;
    //Should we plot statistics about the found features?
    protected boolean plotStatistics = false;
    //Should we NOT calculate accurate mass in the unresampled spectra?
    protected boolean noAccurateMass = false;

    //There was a shift in the way feature strategies were implemented, and it was necessary to leave
    //the "old school" strategies in place.  This code assumes the strategy is "new school" unless told to look
    //among the "old school" strategies.
    boolean oldSchoolStrategy = false;


    /**
     * Just calls initializer
     */
    public FindPeptidesCommandLineModule()
    {
        init();
    }

    /**
     * Define the command name, usage message, help message, and short description, and define all arguments
     */
    protected void init()
    {
        //command name
        mCommandName = "findpeptides";        

        //This takes the place of the auto-generated usage message based on argument definitions.  This is
        //defined explicitly because we don't want this message to change at all from version to version --
        //legacy support.  Most modules will not define this variable.
        mUsageMessage =
                "--findPeptides [--dumpWindow=windowSize] [--out=outfilename] [--outdir=outdirpath] [--start=startScan] [--count=scanCount]\n" +
                "               [--minMz=minMzVal] [--maxMz=maxMzVal] [--strategy=className] [--noAccurateMass] [--accurateMassScans=<int>] [--walkSmoothed] mzxmlfile\n";

        //A longer help message
        mHelpMessage =
                "The findpeptides command finds peptide features in an mzXML file, based on the criteria supplied";

        //A short (single-sentence) description of this command
        mShortDescription = "Find features in an mzXML file based on various critera";

        //Define all of the basic arguments.  The UnnamedSeries argument allows multiple values to be specified without
        //names.  These will all be interpreted as files and validated accordingly
        CommandLineArgumentDefinition[] basicArgDefs =
            {
                    createUnnamedSeriesFileArgumentDefinition(true, "Input mzXML file(s)"),
                    createFileToWriteArgumentDefinition("out", false, "Output File"),
                    createDecimalArgumentDefinition("minmz", false,
                            "Minimum M/Z Value (default: the minimum m/z value in the file)"),
                    createDecimalArgumentDefinition("maxmz", false,
                            "Maximum M/Z Value (default: the maximum m/z value in the file)"),
                    createIntegerArgumentDefinition("start", false,
                            "Minimum scan number", startScan),
                    createIntegerArgumentDefinition("count", false,
                            "Number of scans to search", scanCount),
                    createStringArgumentDefinition("strategy", false,
                            "Class name of a feature-finding strategy implementation"),
                    createDirectoryToReadArgumentDefinition("outdir", false,
                            "Output Directory (for finding features in multiple files)"),
            };
        //add the basic arguments
        addArgumentDefinitions(basicArgDefs);
        //add advanced arguments.  Treated just like basic arguments, but marked as 'advanced' in the UI.  Most
        //users should be happy with the default values for these
        CommandLineArgumentDefinition[] advancedArgDefs =
            {
                    createIntegerArgumentDefinition("dumpwindow", false,
                            "Number of scans around each feature to dump to the file",
                            dumpWindowSize),
                    createBooleanArgumentDefinition("noaccuratemass", false,
                            "Do NOT attempt mass-accuracy adjustment after default peak finding strategy " +
                                    "(by default, adjustment is done)",
                            noAccurateMass),
                    createIntegerArgumentDefinition("accuratemassscans", false,
                            "When attempting to improve mass-accuracy, consider a neighborhood of <int> scans",
                            accurateMassAdjustmentScans),
                    createBooleanArgumentDefinition("plotstats", false,
                             "Plot statistics related to feature-finding",
                             plotStatistics),
                    createBooleanArgumentDefinition("walksmoothed", false,
                             "When calculating feature extents, use smoothed rather than wavelet-transformed spectra)",
                             peakRidgeWalkSmoothed),
            };
        //add the advanced arguments
        addArgumentDefinitions(advancedArgDefs, true);
    }

    /**
     * After the lower-level argument validation has taken place, assign the values to variables, and throw an
     * exception if there's any trouble
     * @throws ArgumentValidationException
     */
    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        //Assign the input and output file argument values
        mzXmlFiles = getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        //Implement the rules about when the "outdir" and "out" args are appropriate.  If there's more than one
        //input file, we want "outdir".  Either is acceptable if there's only one file.  In no case may both
        //be specified.
        if (mzXmlFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }
        else
        {
            if (hasArgumentValue("out"))
                assertArgumentAbsent("outdir");
            else
                assertArgumentPresent("outdir");
        }

        //Assign the variables that represent the m/z and scan boundaries, if present
        if (hasArgumentValue("minmz"))
            minMz = (float) getDoubleArgumentValue("minmz");
        if (hasArgumentValue("maxmz"))
            maxMz = (float) getDoubleArgumentValue("maxmz");
        startScan = getIntegerArgumentValue("start");
        if (hasArgumentValue("count"))
            scanCount = getIntegerArgumentValue("count");

        //parameters related to accurate mass adjustment in the unresampled space
        //The "noaccuratemass" argument is not stored separately, just used to adjust accurateMassAdjustmentScans
        accurateMassAdjustmentScans = getIntegerArgumentValue("accuratemassscans");
        if (hasArgumentValue("noaccuratemass") && getBooleanArgumentValue("noaccuratemass"))
            accurateMassAdjustmentScans = 0;

        dumpWindowSize = getIntegerArgumentValue("dumpwindow");
        peakRidgeWalkSmoothed = getBooleanArgumentValue("walksmoothed");
        //Should we plot statistics on the feature-finding?
        plotStatistics = getBooleanArgumentValue("plotstats");

        //Determine the strategy to be used, and try to instantiate it.  Throw an exception if we fail
        String strategy = getStringArgumentValue("strategy");       

        if (strategy != null)
        {
            try
            {
                featureStrategyClass = FeatureFindingBroker.getFeatureStrategyClass(strategy);
            }
            catch (ClassNotFoundException e)
            {
                throw new ArgumentValidationException("Could not instantiate Feature strategy with name " +
                        strategy, e);
            }

        }
    }


    /**
     * Controller method to oversee the work, which is done in findFeaturesInFile
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        //If there's only one output file, there must be only one input file.  Handle it.
        if (outFile != null)
            findFeaturesInFile(mzXmlFiles[0], outFile);
        else
        {
            //For each input file, create an appropriate output file and call findFeaturesInFile with the pair.
            //Let the user know what's going on.
            for (File mzXmlFile : mzXmlFiles)
            {
                ApplicationContext.setMessage("Processing mzXml file " +
                        mzXmlFile.getAbsolutePath() + " ...");                
                String mzXmlFileName = mzXmlFile.getName();
                int mzXmlFileNameLength = mzXmlFileName.length();
                String outputFileName;
                if (mzXmlFileName.toLowerCase().endsWith(".mzxml"))
                    outputFileName =
                            mzXmlFileName.substring(0, mzXmlFileNameLength -
                                    ".mzxml".length());
                else if (mzXmlFileName.endsWith(".xml"))
                    outputFileName =
                            mzXmlFileName.substring(0, mzXmlFileNameLength -
                                    ".xml".length());
                else
                    outputFileName = mzXmlFileName;
                outputFileName = outputFileName + ".peptides.tsv";
                File outputFile = new File(outDir, outputFileName);
                findFeaturesInFile(mzXmlFile, outputFile);
                ApplicationContext.setMessage("Saved feature file " +
                        outputFile.getAbsolutePath());
            }
        }

    }

    /**
     * Workhorse method that operates on a single input file and output file
     * @param mzXmlFile
     * @param outFile
     * @throws CommandLineModuleExecutionException
     */
    protected void findFeaturesInFile(File mzXmlFile, File outFile)
            throws CommandLineModuleExecutionException
    {
        //Create the PrintWriter that we'll use
        PrintWriter outPW = null;
        try
        {
            //superclass convenience method to create a PrintWriter for a given file
            outPW = getPrintWriter(outFile);
        }
        catch (FileNotFoundException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        //Try to load the MSRun object from the mzXML file.  An MSRun represents all the scans in the run.
        //A soft-reference cache is maintained so not all the scans are in memory at once.
        MSRun run;
        try
        {
            run = MSRun.load(mzXmlFile.getAbsolutePath());
            if (run == null)
                throw new CommandLineModuleExecutionException("Error opening run from file " +
                        mzXmlFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        //Load the m/z range from the file, and, if thee user has specified a range, define the range
        //as the intersection of both ranges.
        FloatRange range = FeatureExtractor.getMzExtractionRange(run);

        float thisRunMinMz = range.min;
        //careful protection from stupid values of minmz and maxmz
        if (hasArgumentValue("minmz"))
        {
            if (minMz > range.max)
                thisRunMinMz = range.max;
            else if (minMz < range.min)
                thisRunMinMz = range.min;
            else
                thisRunMinMz = minMz;

            if (thisRunMinMz != minMz)
                ApplicationContext.infoMessage("Specified minMz value is outside the range for this run.  Adjusting to " +
                        thisRunMinMz);
        }

        float thisRunMaxMz = range.max;
        if (hasArgumentValue("maxmz"))
        {
            if (maxMz > range.max)
            {
                thisRunMaxMz = range.max;
            }
            else if (maxMz < range.min)
                thisRunMaxMz = range.min;
            else
                thisRunMaxMz = maxMz;

            if (thisRunMaxMz != maxMz)
                ApplicationContext.infoMessage("Specified maxMz value is outside the range for this run.  Adjusting to " +
                        thisRunMaxMz);
        }

        if (thisRunMinMz >= thisRunMaxMz)
        {
            throw new CommandLineModuleExecutionException("Empty m/z range specified, so feature finding is not possible.  Quitting");
        }

        //load the scan count from the run, overriding with the user's if provided
        int thisRunScanCount = run.getScanCount();
        if (hasArgumentValue("count") && scanCount < Integer.MAX_VALUE)
            thisRunScanCount = scanCount;

        outPW.flush();

        try
        {
            //Delegate to FeatureFindingBroker for the actual finding of features            
            FeatureSet featureSet = FeatureFindingBroker.findPeptides(
                    run, startScan, thisRunScanCount,
                    PeakCombiner.DEFAULT_MAX_CHARGE,
                    new FloatRange(thisRunMinMz, thisRunMaxMz),
                    dumpWindowSize, accurateMassAdjustmentScans, featureStrategyClass,
                    (null != outFile), peakRidgeWalkSmoothed, plotStatistics);
            //Save the found features to the specified file
            featureSet.save(outPW, dumpWindowSize > 0);
        }
        catch (Exception e)
        {
            //if there are any problems, throw the appropriate type of exception
            ApplicationContext.infoMessage("Error while finding features");
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            //close the PrintWriter if, in fact, it's open
            if (null != outPW)
                outPW.close();
        }
    }
}

