/* 
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFindingBroker;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategy;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.FloatRange;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;


/**
 * Command linemodule for feature finding
 */
public class FindPeptidesCommandLineModule extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FindPeptidesCommandLineModule.class);

    protected File[] mzXmlFiles = null;
    protected File outFile = null;
    protected File outDir = null;       
    protected PrintWriter outPW = null;
    protected int startScan = 1;
    protected int scanCount = Integer.MAX_VALUE;
    protected float minMz = 0;
    protected float maxMz = Float.MAX_VALUE;
    protected int dumpWindowSize = 0;
    protected int accurateMassAdjustmentScans =
            FeatureFinder.DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS;
    protected Class<? extends FeatureStrategy> featureStrategyClass =
            FeatureFinder.DEFAULT_FEATURE_FINDING_CLASS;
    protected boolean peakRidgeWalkSmoothed = false;
    protected boolean plotStatistics = false;

    protected boolean noAccurateMass = false;

    boolean oldSchoolStrategy = false;

    public FindPeptidesCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "findpeptides";        

        mUsageMessage =
                "--findPeptides [--dumpWindow=windowSize] [--out=outfilename] [--outdir=outdirpath] [--start=startScan] [--count=scanCount]\n" +
                "               [--minMz=minMzVal] [--maxMz=maxMzVal] [--strategy=className] [--noAccurateMass] [--accurateMassScans=<int>] [--walkSmoothed] mzxmlfile\n";

        mHelpMessage =
                "The findpeptides command finds peptide features in an mzXML file, based on the criteria supplied";

        mShortDescription = "Find features in an mzXML file based on various critera";

        CommandLineArgumentDefinition[] argDefs =
            {
                    createUnnamedSeriesArgumentDefinition(
                            ArgumentDefinitionFactory.FILE_TO_READ,
                            true, "Input mzXML file(s)"),
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
                    createDirectoryToReadArgumentDefinition("outdir", false,
                            "Output Directory (for finding features in multiple files)"),
                    createBooleanArgumentDefinition("walksmoothed", false,
                             "When calculating feature extents, use smoothed rather than wavelet-transformed spectra)",
                             peakRidgeWalkSmoothed),
                    createBooleanArgumentDefinition("plotstats", false,
                             "Plot statistics related to feature-finding",
                             plotStatistics),

            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");


        mzXmlFiles = getUnnamedSeriesFileArgumentValues();

        if (hasArgumentValue("minmz"))
            minMz = (float) getDoubleArgumentValue("minmz");
        if (hasArgumentValue("maxmz"))
            maxMz = (float) getDoubleArgumentValue("maxmz");
        if (hasArgumentValue("count"))
            scanCount = getIntegerArgumentValue("count");
        accurateMassAdjustmentScans = getIntegerArgumentValue("accuratemassscans");
        if (hasArgumentValue("noaccuratemass") && getBooleanArgumentValue("noaccuratemass"))
            accurateMassAdjustmentScans = 0;

        startScan = getIntegerArgumentValue("start");

        dumpWindowSize = getIntegerArgumentValue("dumpwindow");
        peakRidgeWalkSmoothed = getBooleanArgumentValue("walksmoothed");
        plotStatistics = getBooleanArgumentValue("plotstats");


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
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (outFile != null)
            findFeaturesInFile(mzXmlFiles[0], outFile);
        else
        {
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

    protected void findFeaturesInFile(File mzXmlFile, File outFile)
            throws CommandLineModuleExecutionException
    {
        try
        {
            outPW = getPrintWriter(outFile);
        }
        catch (FileNotFoundException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

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

        if (thisRunMinMz == thisRunMaxMz)
        {
            throw new CommandLineModuleExecutionException("Empty m/z range specified, so feature finding is not possible.  Quitting");
        }


        int thisRunScanCount = run.getScanCount();
        if (hasArgumentValue("count") && scanCount < Integer.MAX_VALUE)
            thisRunScanCount = scanCount;

        outPW.flush(); // test printwriter

        try
        {
            FeatureSet featureSet = FeatureFindingBroker.findPeptides(
                    run, startScan, thisRunScanCount,
                    PeakCombiner.DEFAULT_MAX_CHARGE,
                    new FloatRange(thisRunMinMz, thisRunMaxMz),
                    dumpWindowSize, accurateMassAdjustmentScans, featureStrategyClass,
                    (null != outFile), peakRidgeWalkSmoothed, plotStatistics);

            featureSet.save(outPW, dumpWindowSize > 0);

            outPW.close();
        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Error while finding features");
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (null != outPW)
                outPW.close();
        }
    }
}

