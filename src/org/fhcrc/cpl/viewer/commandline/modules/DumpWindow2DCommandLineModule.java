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
import org.fhcrc.cpl.viewer.DumpWindow2D;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;


/**
     * Dump a 2D window (scans in rows, m/z in columns) around each feature.
     * Only for use by D0/D3 quantitation code; this will be going away.
 */
public class DumpWindow2DCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(DumpWindow2DCommandLineModule.class);

    protected File mzXmlFile;
    protected File featureFile;
    protected File outFile;

    int leadingScans=4;
    int trailingScans=4;
    int leadingPeaks=1;
    int trailingPeaks=3;


    public DumpWindow2DCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "dumpwindow2d";

        mHelpMessage =
                "The dumpwindow2d command creates a feature file that contains a 2D window of spectra around each feature";

        mShortDescription = "Create a feature file that contains a 2D window of spectra around each feature";

        CommandLineArgumentDefinition[] argDefs =
               {
                    createUnnamedFileArgumentDefinition(true, "mzXML File"),
                    new FileToReadArgumentDefinition("features", true, "Feature file"),
                    new FileToWriteArgumentDefinition("out", false, "output file"),
                    new IntegerArgumentDefinition("leadingscans", false,
                            "Leading Scans", leadingScans),
                    new IntegerArgumentDefinition("trailingscans", false,
                            "Trailing Scans", trailingScans),
                    new IntegerArgumentDefinition("leadingpeaks", false,
                            "Leading Peaks", leadingPeaks),
                    new IntegerArgumentDefinition("trailingpeaks", false,
                            "Leading Peaks", leadingPeaks),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzXmlFile = getFileArgumentValue(
                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        featureFile = getFileArgumentValue("features");
        outFile = getFileArgumentValue("out");
        if (null == outFile)
        {
            String outFilePath = mzXmlFile.getPath();
            if ( outFilePath.endsWith(".mzXML") )
                outFilePath = outFilePath.substring(0, outFilePath.length() - ".mzXML".length());
            outFilePath += ".window2D.tsv";
            outFile = new File(outFilePath);
        }

        leadingScans = getIntegerArgumentValue("leadingscans");
        trailingScans = getIntegerArgumentValue("trailingscans");
        leadingPeaks = getIntegerArgumentValue("leadingpeaks");
        trailingPeaks = getIntegerArgumentValue("trailingpeaks");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            DumpWindow2D.dump(mzXmlFile, featureFile, leadingScans, trailingScans,
                              leadingPeaks, trailingPeaks, new File(outFile.getAbsolutePath()));
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

    }

}
