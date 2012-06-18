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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureMassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBoxAndWhiskerChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;


import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PlotMassCalibrationCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PlotMassCalibrationCLM.class);

    protected FeatureSet[] featureSets;

    protected File outFile;
    protected File outBoxWhiskersPlotFile;

    protected double theoreticalMassWavelength =
            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH;

    protected boolean displayPPMLine = false;
    protected double ppmForLine = -1;
    protected boolean showCharts=true;
    protected boolean useMz = false;

    protected PanelWithScatterPlot scatterPlotChartPanel = null;
    protected PanelWithBoxAndWhiskerChart boxWhiskersChartPanel = null;


    public PlotMassCalibrationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "plotmasscalibration";
        mShortDescription = "Visualizes the mass calibration of a feature file as a scatterplot";
        mHelpMessage = "Visualizes the mass calibration of a feature file as a scatterplot, based on work by Wolski et al., 2006.  X axis is feature mass.  Y axis is the distance from the nearest theoretical mass cluster.\n" +
                       " Optionally display lines indicating a given ppm deviation from theoretical mass"; 

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFeatureFileArgumentDefinition(false, null),
                        new DirectoryToReadArgumentDefinition("indir",false, "Directory containing input files, all of which will be used"),
                        new FileToWriteArgumentDefinition("out",false, null),
                        new DecimalArgumentDefinition("theoreticalwavelength",false,
                                "Theoretical mass cluster wavelength",
                                theoreticalMassWavelength),
                        new DecimalArgumentDefinition("ppmline", false,
                                "PPM cutoff to display on plot (default none)", ppmForLine),
                        new FileToWriteArgumentDefinition("outboxwhiskersplot",false, "File to save box-and-whiskers plot to"),
                        new BooleanArgumentDefinition("showcharts",false,"Show charts?", showCharts),
                        new BooleanArgumentDefinition("usemz",false,"Plot m/z instead of mass", useMz),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        Object[] featureSetObjects = getUnnamedSeriesArgumentValues();
        if (featureSetObjects != null)
        {
            assertArgumentAbsent("indir");
            featureSets = new FeatureSet[featureSetObjects.length];
            for (int i=0; i<featureSetObjects.length; i++)
                featureSets[i] = (FeatureSet) featureSetObjects[i];
        }
        else
        {
            assertArgumentPresent("indir");
            File[] featureFiles = getFileArgumentValue("indir").listFiles();
            List<FeatureSet> featureSetsList = new ArrayList<FeatureSet>();
            for (File possibleFeatureFile : featureFiles)
            {
                if (possibleFeatureFile.isDirectory())
                    continue;
                try
                {
                    featureSetsList.add(new FeatureSet(possibleFeatureFile));
                }
                catch (Exception e)
                {
                    throw new ArgumentValidationException("Error loading feature file " + possibleFeatureFile.getName(), e);
                }
            }
            featureSets = featureSetsList.toArray(new FeatureSet[featureSetsList.size()]);
        }
        
        outFile = getFileArgumentValue("out");
        outBoxWhiskersPlotFile = getFileArgumentValue("outboxwhiskersplot");
        showCharts = getBooleanArgumentValue("showcharts");
        useMz = getBooleanArgumentValue("usemz");

        theoreticalMassWavelength = getDoubleArgumentValue("theoreticalwavelength");

        ppmForLine = getDoubleArgumentValue("ppmline");
        if (ppmForLine >= 0)
            displayPPMLine = true;

       
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<Feature[]> featureArrayList = new ArrayList<Feature[]>(featureSets.length);
        for (FeatureSet featureSet : featureSets)
            featureArrayList.add(featureSet.getFeatures());

        scatterPlotChartPanel = FeatureMassCalibrationUtilities.plotMassDefectDeviation(
                featureArrayList, theoreticalMassWavelength, displayPPMLine, ppmForLine, useMz);
        scatterPlotChartPanel.setName("By mass");

        if (showCharts)
            scatterPlotChartPanel.displayInTab();
        if (outFile != null)
        {
            try
            {
                scatterPlotChartPanel.saveChartToImageFile(outFile);
                System.err.println("Saved output image");
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        PanelWithScatterPlot pwspScan = new PanelWithScatterPlot();


        boxWhiskersChartPanel =
                new PanelWithBoxAndWhiskerChart("Mass cluster deviation (Da)");

        for (int i=0; i<featureArrayList.size(); i++)
        {
            Feature[] features = featureArrayList.get(i);
            double[][] scatterPlotData = new double[2][features.length];

            for (int j=0; j<features.length; j++)
            {
                int scan = features[j].getScan();

                float mass = features[j].getMass();
                if (useMz)
                    mass = features[j].getMz();
                double massOffset = (mass -  MassCalibrationUtilities.THEORETICAL_MASS_DEFECT_INTERCEPT)
                        % theoreticalMassWavelength;


                if (massOffset >= 0.5)
                    massOffset -= 1;
                scatterPlotData[0][j] = (float) (scan);
                scatterPlotData[1][j] = (float) massOffset;
            }
            pwspScan.addData(scatterPlotData[0], scatterPlotData[1], "Set " + (i+1) + " (" + featureSets[i].getSourceFile().getName() + ")");
            boxWhiskersChartPanel.addData(scatterPlotData[1], "Set " + (i+1) + " (" + featureSets[i].getSourceFile().getName() + ")");
        }

        pwspScan.setAxisLabels("Scan","Mass Cluster Deviation (Da)");
        pwspScan.setName("By scan");
        if (showCharts)
            pwspScan.displayInTab();

        if (showCharts && featureSets.length > 1)
            boxWhiskersChartPanel.displayInTab();
        if (outBoxWhiskersPlotFile != null)
            try
            {
                boxWhiskersChartPanel.saveChartToImageFile(outBoxWhiskersPlotFile);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save box-and-whiskers plot to file",e);
            }
    }


    public PanelWithBoxAndWhiskerChart getBoxWhiskersChartPanel() {
        return boxWhiskersChartPanel;
    }

    public PanelWithScatterPlot getScatterPlotChartPanel() {
        return scatterPlotChartPanel;
    }
}
