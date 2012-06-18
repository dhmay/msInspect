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

import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class CountMissedCleavagesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CountMissedCleavagesCLM.class);

    protected FeatureSet[] featureSets;


    public CountMissedCleavagesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "countmissedcleavages";
        mShortDescription = "Count missed cleavages in ID'd peptides";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFeatureFileArgumentDefinition(true, null),

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

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (FeatureSet featureSet : featureSets) {
            Map<Integer, Set<String>> numCleavagesPeptideSetMap = new HashMap<Integer, Set<String>>();
            PeptideGenerator pg = new PeptideGenerator();
            pg.setMaxMissedCleavages(0);
            Set<String> allPeptides = new HashSet<String>();
            for (Feature feature : featureSet.getFeatures()) {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                allPeptides.add(peptide);

                Protein protein = new Protein("dummy", peptide.getBytes());
                Peptide[] digestedPeps = pg.digestProtein(protein);
                int numCleavages = digestedPeps.length-1;
                if (numCleavages < 0) ApplicationContext.infoMessage("Negative cleavages?! " + peptide);
                Set<String> peptidesThisNum = numCleavagesPeptideSetMap.get(numCleavages);
                if (peptidesThisNum == null) {
                    peptidesThisNum = new HashSet<String>();
                    numCleavagesPeptideSetMap.put(numCleavages, peptidesThisNum);
                }
                peptidesThisNum.add(peptide);
            }
            ApplicationContext.infoMessage(featureSet.getSourceFile().getName() + ": cleavages:");
            for (int i=0; i<30; i++) {
                if (numCleavagesPeptideSetMap.containsKey(i)) {
                    int count = numCleavagesPeptideSetMap.get(i).size();
                    ApplicationContext.infoMessage("\t" + i + ": " + count +
                       "   (" + Rounder.round(100.0 * (float) (count) / (float) (allPeptides.size()), 1) + "%)");
                }
            }
        }


    }

}
