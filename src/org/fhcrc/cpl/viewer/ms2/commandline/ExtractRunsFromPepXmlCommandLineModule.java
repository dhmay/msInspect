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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.FeatureSelectionParamsCommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;
import org.apache.log4j.Logger;

import java.io.File;


public class ExtractRunsFromPepXmlCommandLineModule extends FeatureSelectionParamsCommandLineModule
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ExtractRunsFromPepXmlCommandLineModule.class);

    protected File inPepXmlFile = null;
    protected File outDirectory = null;
    protected String sourceFileName = null;

    protected static final int FILE_FORMAT_TSV = 0;
    protected static final int FILE_FORMAT_PEPXML = 1;

    protected int outFormat = 1;

    protected boolean populateTimes = false;
    protected File mzXmlDir = null;


    protected static final String[] formatStrings = {
            "tsv",
            "pepxml"
        };

    public ExtractRunsFromPepXmlCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        super.init();
        mCommandName = "extractrunsfrompepxml";
        mShortDescription = "Extract individual runs from a multi-fraction PepXML file";
        mHelpMessage = "Extracts individual runs from a PepXML file containing multiple runs as fractions.  Saves the results as individual PepXML files.";

        CommandLineArgumentDefinition[] argDefs =
            {
                    new FileToReadArgumentDefinition(
                            CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,
                            true, "Input pepXml file"),
                    new DirectoryToWriteArgumentDefinition("outdir", true, "Output Directory"),
                    new StringArgumentDefinition("sourcefilename", false, "Source File Name (without .xml)"),
                    new EnumeratedValuesArgumentDefinition("outformat", false, "Output format", formatStrings, "pepxml"),
                    new BooleanArgumentDefinition("populatetimes", false, "Populate times using mzXML file", populateTimes),
                    new DirectoryToReadArgumentDefinition("mzxmldir", false, "Directory to search for mzXML files (for populating times)")
            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        super.assignArgumentValues();
        inPepXmlFile = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);

        outDirectory = getFileArgumentValue("outdir");
        if (!outDirectory.isDirectory())
        {
            throw new ArgumentValidationException(outDirectory.getAbsolutePath() + " is not a directory.");
        }

        sourceFileName = getStringArgumentValue("sourcefilename");

        if (sourceFileName != null && sourceFileName.endsWith(".xml"))
        {
            System.err.println("removing .xml from source file name, since extension is not stored in combined pepxml file");
            sourceFileName = sourceFileName.substring(0, sourceFileName.lastIndexOf(".xml"));
        }

        String outFormatString = getStringArgumentValue("outformat");
        for (int i=0; i<formatStrings.length; i++)
        {
            String formatString = formatStrings[i];
            if (formatString.equalsIgnoreCase(outFormatString))
            {
                outFormat = i;
                break;
            }
        }

        populateTimes = getBooleanArgumentValue("populatetimes");
        if (populateTimes)
        {
            assertArgumentPresent("mzxmldir","populatetimes");
            mzXmlDir = getFileArgumentValue("mzxmldir");
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (sourceFileName == null)
        {
            ApplicationContext.infoMessage("Extracting all fractions...");
        }
        try
        {
            PepXmlLoader pepXmlLoader = new PepXmlLoader(inPepXmlFile, _log);
            PepXmlLoader.FractionIterator fractionIterator = pepXmlLoader.getFractionIterator();
            boolean foundFraction = false;
            PepXMLFeatureFileHandler featureFileHandler =
               new PepXMLFeatureFileHandler();
            int numFeaturesExtractedSoFar = 0;
            while (fractionIterator.hasNext())
            {
                PepXmlLoader.PepXmlFraction fraction =
                        (PepXmlLoader.PepXmlFraction) fractionIterator.next();
                if (sourceFileName == null ||
                    sourceFileName.equals(fraction.getDataBasename()))
                {
                    foundFraction = true;
                    FeatureSet outFeatureSet =
                       featureFileHandler.createFeatureSetFromPepXMLFraction(fraction, pepXmlLoader);
                    MS2ExtraInfoDef.setFeatureSetModifications(outFeatureSet,
                            fraction.getModifications().toArray(new MS2Modification[0]));
                    MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(outFeatureSet, fraction.getDatabaseLocalPath());

                    outFeatureSet = outFeatureSet.filter(featureSelector);

                    if (populateTimes)
                    {
                        String mzXmlFileName =
                                (fraction.getDataBasename() + ".mzXML");
                        File mzXmlFile = null;

                        boolean foundIt = false;
                        for (String potentialMzXmlFilename : mzXmlDir.list())
                        {
                            if (potentialMzXmlFilename.equalsIgnoreCase(mzXmlFileName))
                            {
                                mzXmlFile = new File(mzXmlDir.getAbsolutePath() + File.separatorChar +
                                        potentialMzXmlFilename);
                                MSRun run = MSRun.load(mzXmlFile.getAbsolutePath());
                                ApplicationContext.setMessage("Located mzXML file " + potentialMzXmlFilename);

                                outFeatureSet.populateTimesForMS2Features(run);
                                foundIt = true;
                                break;
                            }
                        }
                        if (!foundIt)
                            throw new CommandLineModuleExecutionException("Couldn't find source mzXML file for fraction " +fraction.getDataBasename() );
                    }




                    String outSuffix = null;
                    switch(outFormat)
                    {
                        case FILE_FORMAT_PEPXML:
                            outSuffix = ".pep.xml";
                            break;
                        case FILE_FORMAT_TSV:
                            outSuffix = ".pep.tsv";
                            break;
                    }
                    String fractionBaseName = fraction.getDataBasename();
                    if (fractionBaseName.contains("" + File.separatorChar) && fractionBaseName.lastIndexOf(File.separatorChar) < fractionBaseName.length()-1)
                        fractionBaseName =  fractionBaseName.substring(fractionBaseName.lastIndexOf(File.separatorChar) + 1);
                    File outFile =
                            new File(outDirectory.getAbsolutePath() +
                                     File.separatorChar +
                                     fractionBaseName + outSuffix);
                    switch(outFormat)
                    {
                        case FILE_FORMAT_PEPXML:
                            outFeatureSet.savePepXml(outFile, numFeaturesExtractedSoFar+1);
                            break;
                        case FILE_FORMAT_TSV:
                            outFeatureSet.save(outFile);
                            break;
                    }
                    ApplicationContext.infoMessage("Wrote fractional output to " + outFile.getAbsolutePath());
                    if (sourceFileName != null)
                        break;
                    numFeaturesExtractedSoFar += outFeatureSet.getFeatures().length;
                }
            }

            if (!foundFraction)
                ApplicationContext.infoMessage("Could not find fraction with source file " + sourceFileName + ", quitting");
            ApplicationContext.infoMessage("Done.");
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Error!", e);
        }
    }

}
