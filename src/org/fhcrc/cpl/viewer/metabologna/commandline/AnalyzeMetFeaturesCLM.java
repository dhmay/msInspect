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
package org.fhcrc.cpl.viewer.metabologna.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 */
public class AnalyzeMetFeaturesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AnalyzeMetFeaturesCLM.class);

    protected FeatureSet featureSet;
    protected File outFile;

    protected boolean shouldCollapseByMax = false;

    public AnalyzeMetFeaturesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "analyzemetfeatures";
        mShortDescription = "Analyze metabolite features";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFeatureFileArgumentDefinition(true, "Input file"),
                        new FileToWriteArgumentDefinition("out", false, "out"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureSet = this.getUnnamedFeatureSetArgumentValue();
        outFile = getFileArgumentValue("out");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        FeatureSet.FeatureSelector fselManyPeaks = new FeatureSet.FeatureSelector();
        fselManyPeaks.setMinPeaks(4);
        Feature[] features4PlusPeaks = featureSet.getFeatures(fselManyPeaks);


        List<Integer> distsWithinPair = new ArrayList<Integer>();
        List<Integer> distsBetweenPairs = new ArrayList<Integer>();
        List<Feature> badFeatures = new ArrayList<Feature>();
        for (Feature feature : features4PlusPeaks)
        {
            Spectrum.Peak[] peaks = feature.comprised;
//            Feature[] peakFeatures = (Feature[]) peaks;
            int scanDist12 = Math.abs(peaks[0].scan - peaks[1].scan);
            int scanDist23 = Math.abs(peaks[2].scan - peaks[1].scan);
            int scanDist34 = Math.abs(peaks[2].scan - peaks[3].scan);

            int scanDistMax12_34 = Math.max(scanDist12, scanDist34);

            distsWithinPair.add(scanDistMax12_34);
            distsBetweenPairs.add(scanDist23);

            if (scanDist23 > 1 && ((scanDist23 >= scanDistMax12_34 + 2) && scanDist23 >= scanDistMax12_34 * 2))
            {
                badFeatures.add(feature);
            System.err.println(peaks[0].scan + " " + peaks[1].scan + " " +peaks[2].scan + " " +peaks[3].scan + " " );
            }
        }
        new PanelWithScatterPlot(distsWithinPair, distsBetweenPairs, "1-2 diff vs 1-2,3-4 diff").displayInTab();
        new PanelWithHistogram(distsBetweenPairs, "ratio").displayInTab();

        featureSet.setFeatures(badFeatures.toArray(new Feature[badFeatures.size()]));
//        new PanelWithScatterPlot(peaks12Diffs, peaks32Over21DiffPPM, "1-2 diff vs 1-2,3-4 diff").displayInTab();
//        new PanelWithHistogram(peaks32Over21DiffPPM, "ratio").displayInTab();
        if (outFile != null)
        {
            try {featureSet.save(outFile, false,"APML");} catch (IOException e) {throw new CommandLineModuleExecutionException(e);}
        }
    }
}


