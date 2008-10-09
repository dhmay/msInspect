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
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.gui.util.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.awt.*;


/**
 * test
 */
public class ViewSpectrumChartCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ViewSpectrumChartCLM.class);

    protected MSRun run;
    protected int minScan;
    protected int maxScan;
    protected int minMz;
    protected int maxMz;
    protected int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;

    protected int scanLine1 = 0;
    protected int scanLine2 = 0;

    protected boolean staticMode = false;
    protected boolean showScans = true;

    protected boolean showCharts = true;
    protected File outFile = null;
    protected File outScansFile = null;



    protected int dialogWidth = 1000;
    protected int dialogHeight = 1000;
    protected int scansFileImageHeight = 100;
    protected int maxScansImageHeight = 4000;


    public ViewSpectrumChartCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "viewspectrum";

        mHelpMessage ="viewspectrum";
        mShortDescription = "PeptideSpectrumChart";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, "mzXML file"),
                        createIntegerArgumentDefinition("minscan", true, "minimum scan number"),
                        createIntegerArgumentDefinition("maxscan", true, "maximum scan number"),
                        createIntegerArgumentDefinition("minmz", true, "minimum m/z"),
                        createIntegerArgumentDefinition("maxmz", true, "maximum m/z"),
                        createIntegerArgumentDefinition("resolution", false, "resolution (number of breaks per Thompson",
                                resolution),
                        createIntegerArgumentDefinition("scanline1", false, "Line for marking a particular scan",
                                scanLine1),
                        createIntegerArgumentDefinition("scanline2", false, "Another line for marking a particular scan",
                                scanLine2),

                        createBooleanArgumentDefinition("static", false, "Only show static image?  For remote invocation",
                                staticMode),

                        createIntegerArgumentDefinition("height", false, "Window height (also used for spectrum file)",
                                dialogHeight),
                        createIntegerArgumentDefinition("width", false, "Window width", dialogWidth),
                        createIntegerArgumentDefinition("scansfileimageheight", false,
                                "Height of EACH per-scan image, in the output file",
                                scansFileImageHeight),
                        createBooleanArgumentDefinition("showscans", false, "Show individual scan spectra?",
                                showScans),
                        createBooleanArgumentDefinition("showcharts", false, "Show charts at all?",
                                showCharts),
                        createFileToWriteArgumentDefinition("out", false, "Output image file for heatmap"),
                        createFileToWriteArgumentDefinition("outscans", false,
                                "Output image file for single-scan line plots"),
                        createIntegerArgumentDefinition("maxscansimageheight", false,
                                "Maximum overall height for the all-scans line plot image (overrides scansfileimageheight)",
                                maxScansImageHeight)


                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        minScan = getIntegerArgumentValue("minscan");
        maxScan = getIntegerArgumentValue("maxscan");
        minMz = getIntegerArgumentValue("minmz");
        maxMz = getIntegerArgumentValue("maxmz");
        resolution = getIntegerArgumentValue("resolution");
        try
        {
            run = MSRun.load(getUnnamedFileArgumentValue().getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load run",e);
        }

        scanLine1 = getIntegerArgumentValue("scanline1");
        scanLine2 = getIntegerArgumentValue("scanline2");
        staticMode = getBooleanArgumentValue("static");
        showScans = getBooleanArgumentValue("showscans");


        dialogHeight = getIntegerArgumentValue("height");
        dialogWidth = getIntegerArgumentValue("width");
        scansFileImageHeight = getIntegerArgumentValue("scansfileimageheight");
        maxScansImageHeight = getIntegerArgumentValue("maxscansimageheight");

        showCharts = getBooleanArgumentValue("showcharts");
        outFile = getFileArgumentValue("out");
        outScansFile = getFileArgumentValue("outscans");


        if (outScansFile != null && !showScans)
            throw new ArgumentValidationException("Can't write scans to file if we don't show them.  Quitting");
        if (!showCharts && outFile == null && outScansFile == null)
            throw new ArgumentValidationException("Not showing charts, no file to write... nothing to do!  Quitting.");
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PanelWithSpectrumChart spectrumPanel =
                new PanelWithSpectrumChart(run, minScan, maxScan, minMz, maxMz, resolution, scanLine1, scanLine2,
                        showScans);
        spectrumPanel.setName("Spectrum");
        spectrumPanel.setVisible(true);

        DropdownMultiChartDisplayPanel multiChartPanelForScans = new DropdownMultiChartDisplayPanel();
        multiChartPanelForScans.setDisplaySlideshowButton(true);
       

        int spacerHeight = 50;
        int topChartHeight = (dialogHeight - spacerHeight) * 7 / 10 - 10;
        int bottomChartHeight = (dialogHeight - spacerHeight) - topChartHeight;
        int chartWidth = dialogWidth-20;

        if (!showScans)
            topChartHeight = dialogHeight - 20;

        spectrumPanel.setSize(new Dimension(chartWidth, topChartHeight));
        spectrumPanel.setMinimumSize(new Dimension(chartWidth, topChartHeight));

        PanelWithChart spectrumChartToDisplay = spectrumPanel;

        if (showScans)
        {
            Map<Integer, PanelWithLineChart> scanChartMap = spectrumPanel.getScanLineChartMap();
            List<Integer> allScans = new ArrayList<Integer>(scanChartMap.keySet());
            Collections.sort(allScans);

            for (int scan : allScans)
                multiChartPanelForScans.addChartPanel(scanChartMap.get(scan));
        }

        if (staticMode)
        {
            File chartFile = TempFileManager.createTempFile("spectrum_chart.png", this);
            try
            {
                spectrumPanel.saveChartToImageFile(chartFile);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failure saving chart file " + chartFile.getAbsolutePath());
            }

            PanelWithBlindImageChart imageChart = new PanelWithBlindImageChart(chartFile, "spectrum");
            imageChart.setSize(new Dimension(chartWidth, topChartHeight));
            imageChart.setMinimumSize(new Dimension(chartWidth, topChartHeight));

            spectrumChartToDisplay = imageChart;
        }

        ApplicationContext.infoMessage("Saving charts to files...");
        if (outFile != null)
        {
            try
            {
                spectrumChartToDisplay.setSize(dialogWidth, dialogHeight);
                spectrumChartToDisplay.saveChartToImageFile(outFile);
                ApplicationContext.infoMessage("Wrote spectrum to image " + outFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to save image file",e);
            }
        }

        if (outScansFile != null)
        {
            try
            {
                spectrumPanel.savePerScanSpectraImage(dialogWidth, scansFileImageHeight,
                        maxScansImageHeight, outScansFile);
                ApplicationContext.infoMessage("Wrote scans to image " + outScansFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to write image file " +
                        outScansFile.getAbsolutePath(),e);
            }
        }

        if (showCharts)
        {
            JDialog dialog = new JDialog();
            dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

            dialog.setLayout(new GridBagLayout());

            dialog.setTitle("Scans " + minScan + "-" + maxScan);
            GridBagConstraints fullRowGBC = new GridBagConstraints();
            fullRowGBC.gridwidth = GridBagConstraints.REMAINDER;
            dialog.setSize(dialogWidth,dialogHeight);
            dialog.add(spectrumChartToDisplay, fullRowGBC);

            if (showScans)
            {
                multiChartPanelForScans.setSize(dialogWidth, bottomChartHeight);
                multiChartPanelForScans.setMinimumSize(new Dimension(dialogWidth, bottomChartHeight));
                dialog.add(multiChartPanelForScans, fullRowGBC);
            }
            dialog.setVisible(true);
        }





//        spectrumPanel.displayDialog("asdf");
    }

}
