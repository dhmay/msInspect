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
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.viewer.amt.AmtUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.util.*;


/**
 * Show a plot that lets the user determine the mass accuracy of the instrumentation
 */
public class MassAccuracyCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MassAccuracyCommandLineModule.class);

    protected FeatureSet ms1Features = null;
    protected FeatureSet ms2Features = null;
    protected MSRun run = null;
    protected double minPeptideProphet = 0;
    float deltaMass = FeatureSetMatcher.DEFAULT_DELTA_MASS_PPM;
    int deltaMassType = FeatureSetMatcher.DEFAULT_DELTA_MASS_TYPE;
    double deltaScan = FeatureSetMatcher.DEFAULT_DELTA_SCAN;
    double deltaTime = FeatureSetMatcher.DEFAULT_DELTA_TIME;
    double deltaHydrophobicity = FeatureSetMatcher.DEFAULT_DELTA_HYDROPHOBICITY;



    public MassAccuracyCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "massaccuracy";
        mShortDescription = "Determine MS1 mass accuracy by matching to MS2";
        mUsageMessage = CommandLineModule.MODULE_USAGE_AUTOMATIC;
        mHelpMessage =
                "calculate the mass accuracy of ms2 features";

        CommandLineArgumentDefinition[] argDefs =
            {
                    new FeatureFileArgumentDefinition("ms1Features", true, "MS1 features"),
                    new FeatureFileArgumentDefinition("ms2Features", true, "MS2 features"),
                    new FileToReadArgumentDefinition("mzxml", true, "mzXML file"),
                    new DecimalArgumentDefinition("minpprophet", false, "Minimum PeptideProphet score"),
                    new DecimalArgumentDefinition("deltatime", false, "Maximum time between matched features"),
                    new DeltaMassArgumentDefinition("deltamass", false,
                            "Maximum mass difference between matched features (in units of da [Daltons] or ppm [parts per million]",
                            new DeltaMassArgumentDefinition.DeltaMassWithType(deltaMass, deltaMassType))
            };
        addArgumentDefinitions(argDefs);
    }



    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        ms1Features = getFeatureSetArgumentValue("ms1features");
        ms2Features = getFeatureSetArgumentValue("ms2features");
        try
        {
            run = MSRun.load((getFileArgumentValue("mzxml")).getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }
        if (hasArgumentValue("minpprophet"))
            minPeptideProphet = getDoubleArgumentValue("minpprophet");

        FeatureSet.FeatureSelector peptideProphetFeatureSelector = new FeatureSet.FeatureSelector();
        peptideProphetFeatureSelector.setMinPProphet((float) minPeptideProphet);
        ms2Features.filter(peptideProphetFeatureSelector);
            ms2Features.populateTimesForMS2Features(run);
            ms1Features.populateTimesForMS1Features(run);
            Map<String,Double> regressionLineMap =
                    AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                            ms2Features.getFeatures(),
                            ProteomicsRegressionUtilities.REGRESSION_MODE_TIME,
                            false);
            AmtUtilities.recordHydrophobicities(ms2Features,
                                                regressionLineMap,
                                                ProteomicsRegressionUtilities.REGRESSION_MODE_TIME);
            AmtUtilities.recordHydrophobicities(ms1Features,
                                                regressionLineMap,
                                                ProteomicsRegressionUtilities.REGRESSION_MODE_TIME);

        if (hasArgumentValue("deltatime"))
            deltaTime = getDoubleArgumentValue("deltatime");
        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassWithType =
                getDeltaMassArgumentValue("deltamass");
        deltaMass = deltaMassWithType.getDeltaMass();
        deltaMassType = deltaMassWithType.getDeltaMassType();
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        int NUM_TICKS = 50;

        float[][] plotData = new float[3][NUM_TICKS];
    deltaMassType = FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE;
//This block is for plotting curves of mass accuracy over changes in parameters
        for (int i=0; i<NUM_TICKS; i++)
        {
            deltaMass = (float) .04 * i;
            ClusteringFeatureSetMatcher featureSetMatcher =
                    new ClusteringFeatureSetMatcher();
            featureSetMatcher.init(deltaMass, deltaMassType,
                    (float) deltaHydrophobicity);
            FeatureSetMatcher.FeatureMatchingResult featureMatchingResult =
                    featureSetMatcher.matchFeatures(ms1Features, ms2Features);

            double[] massDiffRatios = new double[featureMatchingResult.size()];
            int j=0;
            for (Feature ms1Feature : featureMatchingResult.getMasterSetFeatures())
            {
                Feature ms2Feature = featureMatchingResult.getBestMatch(ms1Feature);
                massDiffRatios[j++] =
                        1000000 * ((Math.abs(ms1Feature.getMass() - ms2Feature.getMass())) / ms2Feature.getMass());
            }

            plotData[0][i] = deltaMass;

            plotData[1][i] = (float) BasicStatistics.median(massDiffRatios);
            plotData[2][i] = (float) BasicStatistics.mean(massDiffRatios);
        }

        XYSeries xySeries = new XYSeries("median");
        XYSeries xySeries1 = new XYSeries("mean");
        for (int i=0; i<plotData[0].length; i++)
        {
            xySeries.add(plotData[0][i],plotData[1][i]);
            xySeries1.add(plotData[0][i],plotData[2][i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(xySeries);
        dataset.addSeries(xySeries1);
        XYPlot plot = new XYPlot(dataset, new NumberAxis() , new NumberAxis(),
                new XYLineAndShapeRenderer());

        ChartDialog cd0 = new ChartDialog(plot);
        cd0.setTitle("Mean and median deviation, for various mass match tolerances");
        cd0.setVisible(true);

    }

}
