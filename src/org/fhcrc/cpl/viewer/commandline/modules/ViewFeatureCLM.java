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
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.viewer.gui.FeatureViewerFrame;
import org.fhcrc.cpl.viewer.gui.FeatureViewer;
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
public class ViewFeatureCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ViewFeatureCLM.class);

    protected MSRun run;
    protected int numScansBeforeAfter = 10;
    protected float mzPaddingBelowAbove = .5f;
    protected int peakPaddingAbove = 2;
    protected int peakPaddingBelow = 2;
    protected int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;


    protected float deltaMassPPM = 40;


    protected boolean showScans = true;
    protected boolean showCharts = true;
    protected File outFile = null;

    protected Feature feature;

    protected int maxPeaks = 10;


    public ViewFeatureCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "viewfeature";

        mHelpMessage ="viewfeature";
        mShortDescription = "viewfeature";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, "feature file"),
                        new FileToReadArgumentDefinition("mzxml", true, "mzXML file"),
                        new DecimalArgumentDefinition("mass", false, "feature mass"),
                        new DecimalArgumentDefinition("deltamassppm", false, "delta mass (ppm)", deltaMassPPM),
                        new IntegerArgumentDefinition("scan", false, "feature scan"),

                        new IntegerArgumentDefinition("scanspadding", false, "scan padding", numScansBeforeAfter),
                        new IntegerArgumentDefinition("peakpadbelow", false, "peak padding below", peakPaddingBelow),
                        new IntegerArgumentDefinition("peakpadabove", false, "peak padding above", peakPaddingAbove),
                        new DecimalArgumentDefinition("mzpadding", false, "m/z paddinb below and above", mzPaddingBelowAbove),
                        new IntegerArgumentDefinition("resolution", false, "resolution (number of breaks per Thompson",
                                resolution),
                        new BooleanArgumentDefinition("showscans", false, "Show individual scan spectra?",
                                showScans),
                        new BooleanArgumentDefinition("showcharts", false, "Show charts onscreen at all?",
                                showCharts),
                        new FileToWriteArgumentDefinition("out", false, "Output image file for heatmap"),
                        new IntegerArgumentDefinition("maxpeaks", false, "max peaks to show", maxPeaks),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        numScansBeforeAfter = getIntegerArgumentValue("scanspadding");
        peakPaddingBelow = getIntegerArgumentValue("peakpadbelow");
        peakPaddingAbove = getIntegerArgumentValue("peakpadabove");

        mzPaddingBelowAbove = getFloatArgumentValue("mzpadding");
        resolution = getIntegerArgumentValue("resolution");


        showScans = getBooleanArgumentValue("showscans");
        showCharts = getBooleanArgumentValue("showcharts");
        outFile = getFileArgumentValue("out");

        if (!hasArgumentValue("mass") && !hasArgumentValue("scan"))
            throw new ArgumentValidationException("Must provide mass or scan");
        Float mass = getFloatArgumentValue("mass");
        Integer scan = getIntegerArgumentValue("scan");

        deltaMassPPM = getFloatArgumentValue("deltamassppm");

        maxPeaks = getIntegerArgumentValue("maxpeaks");

        float deltaMassDa = 0f;
        float minMass = 0;
        float maxMass = 0;
        if (mass != null)
        {
            deltaMassDa = MassUtilities.convertPPMToDa(deltaMassPPM, mass);
            minMass =  mass-deltaMassDa;
            maxMass =  mass+deltaMassDa;
        }

        System.err.println("delta mass in Da: " + deltaMassDa + ".  Min=" + minMass + ", Max=" + maxMass);

        try
        {
            FeatureSet featureSet = new FeatureSet(getUnnamedFileArgumentValue());
            List<Feature> matchingFeatures = new ArrayList<Feature>();
            for (Feature feature : featureSet.getFeatures())
            {
                if ((mass == null || (feature.getMass() >= minMass && feature.getMass() <= maxMass)) &&
                   (scan == null || scan == feature.getScan()))
                {
                    ApplicationContext.infoMessage("Matched: " + feature);
                    matchingFeatures.add(feature);
                }
            }
            if (matchingFeatures.size() > 1 || matchingFeatures.isEmpty())
            {
                ApplicationContext.infoMessage(matchingFeatures.size() + " features matched, quitting.");
                return;
            }
            feature = matchingFeatures.get(0);
            feature.setPeaks(Math.min(maxPeaks,feature.peaks));
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Failed to load feature file");
        }

        try
        {
            run = MSRun.load(getFileArgumentValue("mzxml").getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load run",e);
        }
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {        
        if (feature == null)
            return;
//        JDialog dialog = new JDialog();
//        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
//        dialog.setSize(new Dimension(FeatureViewer.DEFAULT_WIDTH+10, FeatureViewer.DEFAULT_HEIGHT+10));
//        dialog.setTitle("Feature");

        FeatureViewerFrame featureViewerFrame = new FeatureViewerFrame(run);
        featureViewerFrame.getFeatureViewer().setNumScansBeforeAfter(numScansBeforeAfter);
        featureViewerFrame.getFeatureViewer().setMzPaddingBelowAbove(mzPaddingBelowAbove );
        featureViewerFrame.getFeatureViewer().setPeakPaddingAbove( peakPaddingAbove);
        featureViewerFrame.getFeatureViewer().setPeakPaddingBelow( peakPaddingBelow );
        featureViewerFrame.getFeatureViewer().setResolution( resolution);
        featureViewerFrame.getFeatureViewer().setShowIndividualScans( false);

        featureViewerFrame.displayFeature(feature);

        ApplicationContext.infoMessage("Displaying feature " + feature);
//        dialog.add(featureViewerFrame);
try
{
        featureViewerFrame.setVisible(true);
}catch (Exception e) { ApplicationContext.errorMessage("Failed to set visible!",e);}

//        spectrumPanel.displayDialog("asdf");
    }

}
