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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class SummarizeProtXmlCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SummarizeProtXmlCLM.class);

    File protXmlFile;

    public SummarizeProtXmlCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "summarizeprotxml";
        mShortDescription = "Summarize the contents of a protXML file";
        mHelpMessage = "Summarize the contents of a protXML file";
        CommandLineArgumentDefinition[] argDefs =
                {
                    createUnnamedFileArgumentDefinition(true, "protxml file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFile = getUnnamedFileArgumentValue();
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            PanelWithChart sensSpecChart =
                    ProteinUtilities.generateSensSpecChart(protXmlFile);
            sensSpecChart.setName("Sens/Spec");
            sensSpecChart.displayInTab();

            ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

            List<Float> groupProbabilityList = new ArrayList<Float>();

            List<Float> proteinProbabilityList = new ArrayList<Float>();
            List<Float> proteinSpectralCountList = new ArrayList<Float>();

            Iterator<ProteinGroup> iterator = protXmlReader.iterator();
            while (iterator.hasNext())
            {
                ProteinGroup pg = iterator.next();
                groupProbabilityList.add(pg.getProbability());
                
                for (ProtXmlReader.Protein protein : pg.getProteins())
                {
                    proteinProbabilityList.add(protein.getProbability());
                    proteinSpectralCountList.add((float) protein.getTotalNumberPeptides());
                }
            }

            PanelWithHistogram pwh = new PanelWithHistogram(groupProbabilityList, "Probabilities of Proteins");
            pwh.displayInTab();

            int numAbovePoint5 = 0;
            int numAbovePoint75 = 0;
            int numAbovePoint9 = 0;
            int numAbovePoint95 = 0;

            for (int i=0; i<groupProbabilityList.size(); i++)
            {
                float probability = groupProbabilityList.get(i);
                if (probability >= .5f)
                {
                    numAbovePoint5++;
                    if (probability >= .75f)
                    {
                        numAbovePoint75++;
                        if (probability >= .9f)
                        {
                            numAbovePoint9++;
                            if (probability >= .95f)
                            {
                                numAbovePoint95++;
                            }
                        }
                    }
                }
            }

            List<Float> spectralCountsAbovePoint9 = new ArrayList<Float>();

            for (int i=0; i<proteinProbabilityList.size(); i++)
            {
                float probability = proteinProbabilityList.get(i);
                if (probability >= 0.9f)
                    spectralCountsAbovePoint9.add(proteinSpectralCountList.get(i));

            }

            ApplicationContext.infoMessage("Total protein groups: " + groupProbabilityList.size()) ;
            ApplicationContext.infoMessage("\tProb >= .5: " + numAbovePoint5);
            ApplicationContext.infoMessage("\tProb >= .75: " + numAbovePoint75);
            ApplicationContext.infoMessage("\tProb >= .9: " + numAbovePoint9);
            ApplicationContext.infoMessage("\tProb >= .95: " + numAbovePoint95);

            PanelWithHistogram pwhSpecCount = new PanelWithHistogram(spectralCountsAbovePoint9, "Spectral Counts prob .9");
            pwhSpecCount.displayInTab();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
