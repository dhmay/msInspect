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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.MzXmlWriter;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Command linemodule for Saving pieces of mzXML files
 */
public class SaveMzxmlWindowCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SaveMzxmlWindowCommandLineModule.class);

    protected File outFile = null;
    protected File outDir = null;
    protected int startScan = -1, endScan = -1;
    protected float minMz = -1, maxMz = -1;
    protected boolean excludeMS1Scans = false;
    protected File[] inFiles = null;

    public SaveMzxmlWindowCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "savemzxmlwindow";

        mHelpMessage =
                "This command saves out a partial mzXML file to a new file, using the specified\n" +
                "start and end scan numbers and minimum and maximum MZ values.  Any values not \n" +
                "specified default to the start/end values from the original file";

        mShortDescription = "Save a partial window from an mzXml file, limited by scan and MZ";


        CommandLineArgumentDefinition[] argDefs =
            {
                    createUnnamedSeriesFileArgumentDefinition(true, "Input mzXML file(s)"),
                    new FileToWriteArgumentDefinition("out", false, "Output File"),
                    new FileToWriteArgumentDefinition("outdir", false, "Output Directory (for multiple inputs"),
                    new DecimalArgumentDefinition("minmz", false, "Minimum M/Z value"),
                    new DecimalArgumentDefinition("maxmz", false, "Maximum M/Z value"),
                    new IntegerArgumentDefinition("minscan", false, "Minimum scan number"),
                    new IntegerArgumentDefinition("maxscan", false, "Maximum scan number"),
                    new BooleanArgumentDefinition("excludems1", false, "Exclude MS1 scans from output",
                            excludeMS1Scans),
            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFiles = this.getUnnamedSeriesFileArgumentValues();


        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (hasArgumentValue("outdir"))
            assertArgumentAbsent("out");
        else
            assertArgumentPresent("out");

        if (hasArgumentValue("minscan"))
            startScan = getIntegerArgumentValue("minscan");
        if (hasArgumentValue("maxscan"))
            endScan = getIntegerArgumentValue("maxscan");
        if (hasArgumentValue("minmz"))
            minMz = getFloatArgumentValue("minmz");
        if (hasArgumentValue("maxmz"))
            maxMz = getFloatArgumentValue("maxmz");




        //Tell the user about the window
        ApplicationContext.infoMessage("\n" + TextProvider.getText("SCAN_RANGE_START_END",
                                                            "MIN_SCAN",""+startScan,
                                                            "MAX_SCAN",""+endScan));
        ApplicationContext.infoMessage(TextProvider.getText("MZ_RANGE_START_END",
                                                            "MIN_MZ",""+minMz,
                                                            "MAX_MZ",""+maxMz));

        excludeMS1Scans = getBooleanArgumentValue("excludems1");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File file : inFiles)
        {
            File outputFile = outFile;
            if (outputFile == null)
                outputFile = CommandLineModuleUtilities.createOutputFile(file, "filtered.mzXML", outDir);
            ApplicationContext.infoMessage("Handling input file: " + file.getName() + ".  Output file: " + outputFile.getAbsolutePath());
            handleFile(file, outputFile);
        }
    }

    public void handleFile(File inputFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        MSRun run = null;

        try
        {
            run = MSRun.load(inputFile.getAbsolutePath());
            if (run == null)
                throw new CommandLineModuleExecutionException(TextProvider.getText("ERROR_LOADING_FILE"));
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(TextProvider.getText("ERROR_LOADING_FILE"));
        }

        int thisRunStartScan = startScan, thisRunEndScan = endScan;
        float thisRunMinMz = minMz, thisRunMaxMz = maxMz;

        //fix up window boundaries if they're funky or unset
        if (thisRunStartScan < 0)
            thisRunStartScan = 0;
        if (thisRunEndScan < 0)
            thisRunEndScan = run.getScan(run.getScanCount()-1).getNum();
        if (thisRunMinMz < 0)
            thisRunMinMz = run.getMzRange().min;
        if (thisRunMaxMz < 0)
            thisRunMaxMz = run.getMzRange().max;


        MzXmlWriter mzXmlWriter = new MzXmlWriter(run);
        mzXmlWriter.setExcludeMS1Scans(excludeMS1Scans);
        try
        {
            mzXmlWriter.writeSubregion(outputFile ,thisRunStartScan, thisRunEndScan,
                                       thisRunMinMz, thisRunMaxMz);
            ApplicationContext.infoMessage("Done.  Wrote file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
