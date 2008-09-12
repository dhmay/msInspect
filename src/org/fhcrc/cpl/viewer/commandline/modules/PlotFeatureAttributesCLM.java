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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ViewerArgumentDefinitionFactory;
import org.fhcrc.cpl.toolbox.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.apache.log4j.Logger;


import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PlotFeatureAttributesCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PlotFeatureAttributesCLM.class);

    protected FeatureSet featureSet;

    protected File outFile;

    protected int breaks=100;

    protected static final int PEPTIDE_PROPHET = 0;
    protected static final int FRACTIONAL_DELTA_MASS = 1;
    protected static final int INTENSITY = 2;
    protected static final int TIME = 3;
    protected static final int LOGRATIO = 4;
    protected static final int LIGHTAREA = 5;
    protected static final int CHARGE= 6;
    protected static final int SUMSQUARES_DISTANCE= 7;
    protected static final int KL= 8;
    protected static final int PEAKS=9;
    protected static final int SEARCHSCORE=10;
    protected static final int FVAL=11;
    protected static final int MASS=12;
    protected static final int MZ=13;



    protected static final int MODE_HISTOGRAM = 0;
    protected static final int MODE_SCATTER = 1;

    //Bit of a hack, of course.
    protected static final float FEATURE_NO_ATTRIBUTE = -9999.93285f;



    protected final static String[] modeStrings =
            {
                    "histogram",
                    "scatter"
            };

    protected final static String[] modeExplanations =
            {
                    "Histogram",
                    "Scatterplot"
            };




    protected final static String[] attrTypeStrings =
            {
                    "pprophet",
                    "fracdeltamass",
                    "intensity",
                    "time",
                    "ratio",
                    "lightarea",
                    "charge",
                    "sumsquaresdist",
                    "kl",
                    "peaks",
                    "searchscore",
                    "fval",
                    "mass",
                    "mz"
            };

    protected final static String[] attrTypeExplanations =
            {
                    "PeptideProphet probabilities",
                    "Fractional Delta Mass",
                    "(Maximum) Intensity",
                    "Retention Time",
                    "Ratio",
                    "Light Area (quantitated peptides)",
                    "Charge",
                    "Sum-of-squares distance",
                    "K/L Score",
                    "Number of peaks",
                    "Search score (must provide score name)",
                    "fval",
                    "mass",
                    "mz"                    
            };

    protected int mode=-1;

    protected int xAttrType=-1;
    protected int yAttrType=-1;

    protected boolean logMode = false;


    protected String searchScoreName = null;




    public PlotFeatureAttributesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "plotfeatureattributes";
        mShortDescription = "Histograms attributes of featuresets";
        mHelpMessage = "asdfasdf";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createEnumeratedArgumentDefinition("attribute",true,attrTypeStrings,
                                attrTypeExplanations),
                        createEnumeratedArgumentDefinition("attribute2",false,attrTypeStrings,
                                attrTypeExplanations),
                        createBooleanArgumentDefinition("logmode",false, "Log mode", logMode),
                        createStringArgumentDefinition("searchscore",false,
                                "Search score name (for searchscore mode)"),
                        createUnnamedArgumentDefinition(
                                ViewerArgumentDefinitionFactory.FEATURE_FILE,true,
                                "Feature file"),
                        createFileToWriteArgumentDefinition("out",false, null),
                        createIntegerArgumentDefinition("breaks",false, "Number of breaks", breaks)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        if (hasArgumentValue("attribute2"))
            mode = MODE_SCATTER;
        else
            mode = MODE_HISTOGRAM;

        xAttrType = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("attribute")).getIndexForArgumentValue(getStringArgumentValue("attribute"));

        if (hasArgumentValue("attribute2"))
            yAttrType =  ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("attribute2")).getIndexForArgumentValue(getStringArgumentValue("attribute2"));

        searchScoreName=getStringArgumentValue("searchscore");

        if (xAttrType == SEARCHSCORE || yAttrType == SEARCHSCORE)
            assertArgumentPresent("searchscore");
        else
            assertArgumentAbsent("searchscore");

        featureSet = getFeatureSetArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        outFile = getFileArgumentValue("out");
        breaks = getIntegerArgumentValue("breaks");

        logMode = getBooleanArgumentValue("logmode");
    }

    protected float getAttributeFromFeature(Feature feature, int attrType)
            throws CommandLineModuleExecutionException
    {
        float result = FEATURE_NO_ATTRIBUTE;

        try
        {
        switch(attrType)
        {
            case PEPTIDE_PROPHET:
                if (MS2ExtraInfoDef.hasPeptideProphet(feature))
                    result = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                break;
            case FRACTIONAL_DELTA_MASS:
                try
                {
                    double deltaMass =
                            MS2ExtraInfoDef.getDeltaMass(feature);
                    if (deltaMass < 0)
                        deltaMass = -deltaMass;
                    while (deltaMass > 0.5)
                        deltaMass -= 1;
                    if (deltaMass < 0)
                        deltaMass = -deltaMass;
                    result = (float)deltaMass;
                }
                catch (Exception e){}
                break;
            case INTENSITY:
                result = feature.getIntensity();
                break;
            case TIME:
                result = feature.getTime();
                break;
            case LOGRATIO:
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    result =  (float)
                            Math.max(.001,IsotopicLabelExtraInfoDef.getRatio(feature));
                break;
            case LIGHTAREA:
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    result =  (float)
                            Math.max(.001,IsotopicLabelExtraInfoDef.getLightIntensity(feature));
                break;
            case CHARGE:
                result = (float) feature.getCharge();
                break;
            case SUMSQUARES_DISTANCE:
                result =  feature.getSumSquaresDist();
                break;
            case KL:
                result =  feature.getKl();
                break;
            case PEAKS:
                result = (float) feature.getPeaks();
                break;
            case SEARCHSCORE:
                result = Float.parseFloat(MS2ExtraInfoDef.getSearchScore(feature, searchScoreName));
                break;
            case FVAL:
                result = (float) MS2ExtraInfoDef.getFval(feature);
                break;
            case MASS:
                result = feature.getMass();
            case MZ:
                result = feature.getMz();                
        }

        if (logMode)
            result = (float) Math.log(result);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to load attribute " + attrType,e);
        }
        return result;
    }

    protected String getAttributeTitle(int attrType)
    {
        String title = "";
        switch(attrType)
        {
            case PEPTIDE_PROPHET:
                title = "PeptideProphet Probabilities";
                break;
            case FRACTIONAL_DELTA_MASS:
                title = "Fractional Delta Masses";
                break;
            case INTENSITY:
                title = "Intensities";
                break;
            case TIME:
                title = "Retention times";
                break;
            case LOGRATIO:
                title = "Ratios";
                break;
            case LIGHTAREA:
                title= "Light Areas";
                break;
            case CHARGE:
                title= "Charge";
                break;
            case SUMSQUARES_DISTANCE:
                title= "Sum-of-squares Distance";
                break;
            case KL:
                title= "K/L Score";
                break;
            case PEAKS:
                title= "Peaks";
                break;
            case SEARCHSCORE:
                title=searchScoreName;
                break;
            case FVAL:
                title="fval";
                break;
        }

        if (logMode)
            title = title + " (log)";
        return title;
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (logMode)
            ApplicationContext.infoMessage("Plotting in log mode");

        PanelWithChart panelWithChart = null;
        switch (mode)
        {
            case MODE_HISTOGRAM:
                panelWithChart = buildHistogram();
                break;
            case MODE_SCATTER:
                panelWithChart = buildScatter();
                break;
        }

        panelWithChart.displayDialog(panelWithChart.getName());

        if (outFile != null)
        {
            try
            {
                panelWithChart.saveChartToImageFile(outFile);
                System.err.println("Saved output image");
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
    }

    protected PanelWithChart buildHistogram()
            throws CommandLineModuleExecutionException
    {

        String title = getAttributeTitle(xAttrType);



        List<Float> featureAttributes = new ArrayList<Float>();
        for (int i=0; i<featureSet.getFeatures().length; i++)
        {
            Feature feature = featureSet.getFeatures()[i];

            float xAttr = getAttributeFromFeature(feature, xAttrType);
            if (xAttr != FEATURE_NO_ATTRIBUTE)
                featureAttributes.add(xAttr);
        }
        ApplicationContext.infoMessage("# datapoints: " + featureAttributes.size());
        ApplicationContext.infoMessage("Mean: " + BasicStatistics.mean(featureAttributes));
        ApplicationContext.infoMessage("Median: " + BasicStatistics.median(featureAttributes));


        PanelWithHistogram panelWithHistogram =
                new PanelWithHistogram(featureAttributes, title, breaks);

        return panelWithHistogram;
    }

    protected PanelWithChart buildScatter()
            throws CommandLineModuleExecutionException
    {

        String title = getAttributeTitle(xAttrType) + " (x) vs. " + getAttributeTitle(yAttrType) + " (y)";

        List<Float> featureXAttributes = new ArrayList<Float>();
        List<Float> featureYAttributes = new ArrayList<Float>();

        for (int i=0; i<featureSet.getFeatures().length; i++)
        {
            Feature feature = featureSet.getFeatures()[i];
            float xAttr = getAttributeFromFeature(feature, xAttrType);
            float yAttr = getAttributeFromFeature(feature, yAttrType);

            if (xAttr != FEATURE_NO_ATTRIBUTE && yAttr != FEATURE_NO_ATTRIBUTE)
            {
                featureXAttributes.add(getAttributeFromFeature(feature, xAttrType));
                featureYAttributes.add(getAttributeFromFeature(feature, yAttrType));
            }
        }
        ApplicationContext.infoMessage("# datapoints: " + featureXAttributes.size());
        ApplicationContext.infoMessage("X Mean: " + BasicStatistics.mean(featureXAttributes));
        ApplicationContext.infoMessage("X Median: " + BasicStatistics.median(featureXAttributes));
        ApplicationContext.infoMessage("Y Mean: " + BasicStatistics.mean(featureYAttributes));
        ApplicationContext.infoMessage("Y Median: " + BasicStatistics.median(featureYAttributes));


        PanelWithScatterPlot panelWithScatterPlot =
                new PanelWithScatterPlot(featureXAttributes, featureYAttributes, title);

        return panelWithScatterPlot;
    }

}
