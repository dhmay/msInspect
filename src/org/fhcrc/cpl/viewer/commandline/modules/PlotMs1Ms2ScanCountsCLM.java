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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.DirectoryToReadArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToWriteArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;


import java.io.File;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PlotMs1Ms2ScanCountsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PlotMs1Ms2ScanCountsCLM.class);

    protected File[] mzXmlFiles;
    protected File outFile;

    public PlotMs1Ms2ScanCountsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "plotms1ms2scancounts";
        mShortDescription = "Plots scan counts";
        mHelpMessage = "argh";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(false, null),
                        new DirectoryToReadArgumentDefinition("mzxmldir", false, "directory of mzXml Files"),
                        new FileToWriteArgumentDefinition("out",false, null),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzXmlFiles = getUnnamedSeriesFileArgumentValues();
        if (mzXmlFiles == null)
        {
            assertArgumentPresent("mzxmldir");
            File mzXmlDir = getFileArgumentValue("mzxmldir");
            mzXmlFiles = mzXmlDir.listFiles();
        }
        else
            assertArgumentAbsent("mzxmldir");
        List<File> filesToProcessList = new ArrayList<File>();
        for (File file : mzXmlFiles)
            if (file.getName().toLowerCase().endsWith("mzxml"))
                filesToProcessList.add(file);
        mzXmlFiles = filesToProcessList.toArray(new File[filesToProcessList.size()]);
        outFile = getFileArgumentValue("out");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PanelWithLineChart pwlc = new PanelWithLineChart();

        double[] ms1ScanCounts = new double[mzXmlFiles.length];
        double[] ms2ScanCounts = new double[mzXmlFiles.length];
        double[] xValues = new double[mzXmlFiles.length];

        try
        {
            for (int i=0; i<mzXmlFiles.length; i++)
            {
                xValues[i] = i+1;

                MSRun run = MSRun.load(mzXmlFiles[i].getAbsolutePath());
                ms1ScanCounts[i] = run.getScanCount();
                ms2ScanCounts[i] = run.getMS2Scans().length;

                _log.debug(mzXmlFiles[i].getAbsolutePath() + ": ms1=" + ms1ScanCounts[i] +
                        ", ms2=" + ms2ScanCounts[i]);
            }

            pwlc.addData(xValues, ms1ScanCounts, "MS1 Scans");
            pwlc.addData(xValues, ms2ScanCounts, "MS2 Scans");
            pwlc.displayInTab();
            if (outFile != null)
                pwlc.saveChartToImageFile(outFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
