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

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;


import java.io.File;
import java.io.IOException;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PlotDeltaMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PlotDeltaMassesCLM.class);

    protected FeatureSet[] featureSets;

    protected File outFile;

    protected boolean showRegressionLine = false;
    protected double ppmForLine = -1;


    protected boolean usePPM = false;
    protected boolean fractional = true;


    public PlotDeltaMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "plotdeltamasses";
        mShortDescription = "Plot deltaMass values against feature masses";
        mHelpMessage = "Plot deltaMass values against feature masses, and against time.  Also histograms delta " +
                       "masses. This tool can be used to " +
                       "assess mass calibration as well as precursor mass accuracy";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFeatureFileArgumentDefinition(true, null),
                        new FileToWriteArgumentDefinition("out",false, "Output file for chart"),
                        new DecimalArgumentDefinition("ppmline", false,
                                "PPM cutoff to display on plot (default none)", ppmForLine),
                        new BooleanArgumentDefinition("useppm", false,
                                "Use PPM instead of Daltons", usePPM),
                        new BooleanArgumentDefinition("showregressionline", false,
                                "Show a simple regression line through the data", showRegressionLine),
                        new BooleanArgumentDefinition("fractional", false,
                                "Fractional masses?", fractional),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        Object[] featureSetObjects = getUnnamedSeriesArgumentValues();
        featureSets = new FeatureSet[featureSetObjects.length];
        for (int i=0; i<featureSetObjects.length; i++)
            featureSets[i] = (FeatureSet) featureSetObjects[i];
        outFile = getFileArgumentValue("out");

        ppmForLine = getDoubleArgumentValue("ppmline");
        showRegressionLine = getBooleanArgumentValue("showregressionline");
        usePPM = getBooleanArgumentValue("useppm");
        fractional = getBooleanArgumentValue("fractional");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
        PanelWithScatterPlot pwspXTime = new PanelWithScatterPlot();


        for (FeatureSet featureSet : featureSets)
        {
            double[] featureCalculatedMasses =
                    new double[featureSet.getFeatures().length];
            double[] featureScanNumbers =
                    new double[featureSet.getFeatures().length];
            double[] featureDeltaMasses =
                    new double[featureSet.getFeatures().length];
            double[] featurePeptideProphets =
                    new double[featureSet.getFeatures().length];

            for (int i=0; i<featureCalculatedMasses.length; i++)
            {
                Feature feature = featureSet.getFeatures()[i];
                featureCalculatedMasses[i] = feature.getMass();
                featureScanNumbers[i] = feature.getScan();
                float deltaMass = MS2ExtraInfoDef.getDeltaMass(feature);

                featureDeltaMasses[i] = deltaMass;

                if (fractional)
                {

                    featureDeltaMasses[i] = ((deltaMass + .5) %1) - .5;
                    if (featureDeltaMasses[i] < -.5)
                        featureDeltaMasses[i] += 1;
                }

                if (usePPM)
                {
                    featureDeltaMasses[i] =
                            featureDeltaMasses[i] * 1000000.0 / feature.getMass();
                }

                featurePeptideProphets[i] = MS2ExtraInfoDef.getPeptideProphet(feature);

            }
            ApplicationContext.infoMessage("Mean deltamass: " + BasicStatistics.mean(featureDeltaMasses));

            pwsp.addData(featureCalculatedMasses, featureDeltaMasses,
                    featureSet.getSourceFile().getName());
            pwspXTime.addData(featureScanNumbers, featureDeltaMasses,
                    featureSet.getSourceFile().getName());

            double sumAbsoluteDeltaMasses = 0;
            double[] absoluteDeltaMasses = new double[featureDeltaMasses.length];
            for (int i=0; i<featureDeltaMasses.length; i++)
            {
                absoluteDeltaMasses[i] = Math.abs(featureDeltaMasses[i]);
                sumAbsoluteDeltaMasses += absoluteDeltaMasses[i];

            }
            ApplicationContext.infoMessage("Mean absolute deltamass: " +
                    (sumAbsoluteDeltaMasses / featureDeltaMasses.length));

            if (!fractional)
            {
                int[] countsEachNumDaltonsOffPlusFive = new int[10];
                for (double deltaMass : featureDeltaMasses)
                {
                    int index = (int) Math.round(deltaMass) + 5;
                    if (index < 0 || index > 9)
                        _log.info("out-of-range deltamass " + deltaMass);
                    countsEachNumDaltonsOffPlusFive[index]++;
                }
                int[] numsDaltonsOff = new int[10];
                ApplicationContext.setMessage("Features each # of Daltons off:");
                for (int i=0; i< 10; i++)
                {
                    numsDaltonsOff[i] = (i-5);
                    ApplicationContext.setMessage("\t " + (i-5) + ": " + countsEachNumDaltonsOffPlusFive[i]);
                }
                PanelWithLineChart pwlc = new PanelWithLineChart(numsDaltonsOff, countsEachNumDaltonsOffPlusFive,
                        "Counts each # of Daltons Off");
                pwlc.setName("Counts each # of Daltons Off");
                pwlc.displayInTab();
            }

            PanelWithHistogram histPanelAbs = new PanelWithHistogram(absoluteDeltaMasses, "Absolute deltaMass", 200);
            histPanelAbs.setName(usePPM ? "Absolute deltaMass (ppm)" : "Absolute deltaMass (Da)");
            histPanelAbs.displayInTab();

            PanelWithHistogram histPanel = new PanelWithHistogram(featureDeltaMasses, "deltaMass", 200);
            histPanel.setName(usePPM ? "deltaMass (ppm)" : "deltaMass (Da)");
            histPanel.displayInTab();

            PanelWithScatterPlot pwspPepProph = new PanelWithScatterPlot(featurePeptideProphets, featureDeltaMasses,
                    "deltaMass vs. pProphet");
            pwspPepProph.setAxisLabels("peptideProphet", usePPM ? "deltaMass (ppm)" : "deltaMass (Da)");
            pwspPepProph.displayInTab();

        }
        if (showRegressionLine)
            pwsp.addRegressionLine();
        String yAxisTitle = "Fractional Delta Mass";
        if (usePPM)
            yAxisTitle = "Fractional Delta Mass (PPM)";
        pwsp.setAxisLabels("Peptide Mass", yAxisTitle);
        pwsp.setName("deltaMass vs. Mass");
        pwsp.displayInTab();

        pwspXTime.setAxisLabels("Scan Number", yAxisTitle);
        pwspXTime.setName("deltaMass vs. Time");
        pwspXTime.displayInTab();

        if (outFile != null)
        {
            try
            {
                pwsp.saveChartToImageFile(outFile);
                System.err.println("Saved output image");
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }


    }

}
