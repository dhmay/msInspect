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
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for feature filtering
 */
public class FilterFeaturesCommandLineModule extends FeatureSelectionParamsCommandLineModule
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FilterFeaturesCommandLineModule.class);

    protected File[] inFeatureFiles = null;
    protected File outFile = null;
    protected File outDir = null;

    protected double minSearchScore = Float.MIN_VALUE;
    protected double maxSearchScore = Float.MAX_VALUE;

    protected String searchScoreName = null;


    public static final int OUT_FORMAT_MSINSPECT=0;
    public static final int OUT_FORMAT_PEPXML=1;

    protected int outFormat = OUT_FORMAT_MSINSPECT;

    protected final static String[] outFormatStrings = {"msinspect",
            "pepxml"};
    protected final static String[] outFormatExplanations = {"msInspect .tsv format",
            "pepXML format"};

    public FilterFeaturesCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        super.init();
        mCommandName = "filter";
        mUsageMessage =
                "--filter [--out=filename] [--outdir=directory] [--minMz=float]\n" +
                        "[--maxMz=float] [--minMass=float] [--maxMass=float]\n" +
                        "[--minPeaks=int] [--minCharge=int] [--maxCharge=int] [--maxKL=float]\n" +
                        "[--minIntensity=float] [--minTotalIntensity=float]\n" +
                        "[--minTime=float] [--maxTime=float]\n" +
                        "[--scanFirst=int] [scanLast=int] [--minScans=int]\n" +
                        "[--minpprophet=float] [--maxmassdeviationppm] [--outputpepxml]\n" +
                        "[--maxsumsquaresdist=float] [--minsearchscore=float] [--maxsearchscore=float] " +
                        "[--searchscorename=string] [--outformat=msinspect|pepxml]\n" +
                        "featurefile [featurefile] [featurefile] ...";

        mHelpMessage =
                "The filter command allows you to one or more feature files based on various criteria.\n" +
                "The out parameter specifies a feature file for output.  You may filter features\n" +
                "based on ma, mass, peaks, charge, kl, intensity, time, or scans.\n" +
                "If you wish to filter more than one file at once, you will need to specify\n" +
                "the outdir parameter, to indicate the output directory.  The output filenames\n" +
                "will be based on the input filenames, ending in '.filtered.tsv'";

        mShortDescription = "Filter feature sets based on various critera";

        CommandLineArgumentDefinition[] argDefs =
            {
                    createUnnamedSeriesArgumentDefinition(
                            ArgumentDefinitionFactory.FILE_TO_READ,
                            true, "Input feature file(s)"),
                    createFileToWriteArgumentDefinition("out", false, "Output File"),
                    createDirectoryToReadArgumentDefinition("outdir", false,
                            "Output Directory (for filtering multiple files)"),
                    createDecimalArgumentDefinition("minsearchscore", false,
                            "Minimum search score", minSearchScore),
                    createDecimalArgumentDefinition("maxsearchscore", false,
                            "Maximum search score", maxSearchScore),
                    createStringArgumentDefinition("searchscorename", false,
                            "Search score name (for minsearchscore or maxsearchscore)"),
                    createEnumeratedArgumentDefinition("outformat",false,outFormatStrings,
                            outFormatExplanations)

            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        super.assignArgumentValues();
        inFeatureFiles = getUnnamedSeriesFileArgumentValues();

        minSearchScore = getDoubleArgumentValue("minsearchscore");
        maxSearchScore = getDoubleArgumentValue("maxsearchscore");
        searchScoreName = getStringArgumentValue("searchscorename");

        if (hasArgumentValue("outformat"))
        {
            outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(getStringArgumentValue("outformat"));
        }


        if (hasArgumentValue("searchscorename"))
        {
            if (!hasArgumentValue("minsearchscore") && !hasArgumentValue("maxsearchscore"))
                throw new ArgumentValidationException("If a searchscorename is specified, a minsearchscore or " +
                        "maxsearchscore (or both) must also be specified");
        }
        if (hasArgumentValue("maxsearchscore") || hasArgumentValue("maxsearchscore"))
        {
            if (!hasArgumentValue("searchscorename"))
                throw new ArgumentValidationException("If a minsearchscore or " +
                        "maxsearchscore (or both) is specified, searchscorename must also be specified");
        }

        if (hasArgumentValue("out"))
        {
            assertArgumentAbsent("outdir");
        }
        else
        {
            assertArgumentPresent("outdir");
        }

        if (hasArgumentValue("outdir"))
        {
            assertArgumentAbsent("out");
        }

        outDir = getFileArgumentValue("outdir");
        outFile = getFileArgumentValue("out");

        if (inFeatureFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        FeatureSet currentFeatureSet = null;
        if (inFeatureFiles.length == 1)
        {
            try
            {
                currentFeatureSet = new FeatureSet(inFeatureFiles[0]);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Error opening feature file " +
                        inFeatureFiles[0], e);
            }
            if (outFile == null)
                outFile = new File(outDir, createFilteredFeatureFileFilename(inFeatureFiles[0].getName(), outFormat));
            filterFeatureFile(currentFeatureSet, outFile);
        }
        else
        {
            for (File inFeatureFile : inFeatureFiles)
            {
                File inputFile = inFeatureFile;
                try
                {
                    currentFeatureSet = new FeatureSet(inFeatureFile);
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Error opening feature file " +
                            inFeatureFile, e);
                }



                String inputFileName = inputFile.getName();

                File outputFile = new File(outDir, createFilteredFeatureFileFilename(inputFileName, outFormat));
                filterFeatureFile(currentFeatureSet, outputFile);
                ApplicationContext.setMessage("Saved filtered feature file " +
                        outputFile);
            }
        }
        ApplicationContext.setMessage("Done saving filtered features.");
    }

    /**
     * Decide the filename of the output file, if not specified by user
     * @param inputFileName
     * @return
     */
    public static String createFilteredFeatureFileFilename(String inputFileName, int outFormat)
    {
        int inputFileNameLength = inputFileName.length();
        String outputFileName;

        if (inputFileName.endsWith(".filtered.tsv"))
            outputFileName = inputFileName;
        else
        {
            if (inputFileName.endsWith(".peptides.tsv"))
                outputFileName =
                        inputFileName.substring(0, inputFileNameLength -
                                ".peptides.tsv".length());

            else if (inputFileName.endsWith(".filtered.tsv"))
                outputFileName = inputFileName;
            else if (inputFileName.endsWith(".features.tsv"))
                outputFileName =
                        inputFileName.substring(0, inputFileNameLength -
                                ".features.tsv".length());
            else if (inputFileName.endsWith(".tsv"))
                outputFileName =
                        inputFileName.substring(0, inputFileNameLength -
                                ".tsv".length());
            else if (inputFileName.endsWith(".pep.xml"))
                outputFileName =
                        inputFileName.substring(0, inputFileNameLength -
                                ".pep.xml".length());
            else if (inputFileName.endsWith(".xml"))
                outputFileName =
                        inputFileName.substring(0, inputFileNameLength -
                                ".xml".length());
            else
                outputFileName = inputFileName;

            switch (outFormat)
            {
                case OUT_FORMAT_PEPXML:
                    outputFileName = outputFileName + ".filtered.pep.xml";
                    break;
                case OUT_FORMAT_MSINSPECT:
                    outputFileName = outputFileName + ".filtered.tsv";
                    break;
            }
        }
        return outputFileName;
    }

    protected void filterFeatureFile(FeatureSet inFeatureSet, File outputFile)
        throws CommandLineModuleExecutionException
    {
        PrintWriter out;

        try
        {
            out = getPrintWriter(outputFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Error opening file " + outFile + " for writing",e);
        }
        //dhmay adding failure handling
        if (inFeatureSet.getLoadStatus() == FeatureSet.FEATURESET_LOAD_SUCCESS)
        {
            ApplicationContext.setMessage("Filtering");
            inFeatureSet = inFeatureSet.filter(featureSelector);

            if (searchScoreName != null)
            {
                List<Feature> passingFeatures = new ArrayList<Feature>();
                for (Feature feature : inFeatureSet.getFeatures())
                {
                    float searchScore = Float.parseFloat(MS2ExtraInfoDef.getSearchScore(feature, searchScoreName));
                    if (searchScore <= maxSearchScore && searchScore >= minSearchScore)
                    {
                         passingFeatures.add(feature);
                    }
                }
                inFeatureSet.setFeatures(passingFeatures.toArray(new Feature[passingFeatures.size()]));
            }

            switch (outFormat)
            {
                case OUT_FORMAT_MSINSPECT:
                    inFeatureSet.save(out);
                    break;
                case OUT_FORMAT_PEPXML:
                    try
                    {
                        inFeatureSet.savePepXml(outputFile);
                    }
                    catch (IOException e)
                    {
                        throw new CommandLineModuleExecutionException(e);
                    }
                    break;

            }
        }
        else
        {
            throw new CommandLineModuleExecutionException("Failure loading features.  Cause: " + inFeatureSet.getLoadStatusMessage(),false);
        }

        out.flush();
        if (null != out)
            out.close();
    }

    /*
    * dhmay adding to centralize PrintWriter creation
    * Utility method: given an output filename, checks if null.  If not, returns a PrintWriter
    * on that file. If so, returns a PrintWriter on System.err
    *
    * @param outFileName the output filename specified on the command line
    * @return an appropriate printwriter
    */
    protected static PrintWriter getPrintWriter(File outFile) throws java.io.FileNotFoundException
    {
      PrintWriter pw = null;
      if (null != outFile)
      {
        try
        {
          pw = new PrintWriter(new FileOutputStream(outFile));
        }
        catch (java.io.FileNotFoundException e)
        {
          System.err.println("Error creating PrintWriter from file " + outFile.getAbsolutePath() + ", file not found");
          throw e;
        }
      }
      else
        pw = new PrintWriter(System.out);

      return pw;
    }
}