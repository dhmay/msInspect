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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
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

    protected File protXmlFiles[];
    protected boolean showCharts = false;
    protected float minProteinProphet = 0;

    public SummarizeProtXmlCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "summarizeprotxml";
        mShortDescription = "Summarize the contents of one or more protXML files";
        mHelpMessage = "Summarize the contents of one or more protXML files";
        CommandLineArgumentDefinition[] argDefs =
                {
                    createUnnamedSeriesFileArgumentDefinition(true, "protxml file"),
                        createBooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        createDecimalArgumentDefinition("minpprophet", false, "Min proteinprophet for MA plot",
                                minProteinProphet),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();
        showCharts = getBooleanArgumentValue("showcharts");
        minProteinProphet = getFloatArgumentValue("minpprophet");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("File\tGroups\tPoint5\tPoint75\tPoint9\tPoint95");
        for (File protXmlFile : protXmlFiles)
            summarizeFile(protXmlFile);
    }

    protected void summarizeFile(File protXmlFile) throws CommandLineModuleExecutionException
    {
        try
        {
            if (showCharts)
            {
                PanelWithChart sensSpecChart =
                        ProteinUtilities.generateSensSpecChart(protXmlFile);
                sensSpecChart.setName("Sens/Spec");
                sensSpecChart.displayInTab();
            }

            ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

            List<Float> groupProbabilityList = new ArrayList<Float>();

            List<Float> proteinProbabilityList = new ArrayList<Float>();
            List<Float> proteinSpectralCountList = new ArrayList<Float>();
            List<Float> proteinRatioForMAList = new ArrayList<Float>();
            List<Float> proteinSpectralCountForMAList = new ArrayList<Float>();


            Iterator<ProteinGroup> iterator = protXmlReader.iterator();
            while (iterator.hasNext())
            {
                ProteinGroup pg = iterator.next();

                groupProbabilityList.add(pg.getProbability());

                for (ProtXmlReader.Protein protein : pg.getProteins())
                {

                    proteinProbabilityList.add(protein.getProbability());
                    proteinSpectralCountList.add((float) protein.getTotalNumberPeptides());
                    if (protein.getProbability() >= minProteinProphet && protein.getQuantitationRatio() != null)
                    {
                        proteinRatioForMAList.add(protein.getQuantitationRatio().getRatioMean());
                        proteinSpectralCountForMAList.add((float) protein.getTotalNumberPeptides());
                    }

                }
            }

            if (showCharts)
            {
                PanelWithHistogram pwh = new PanelWithHistogram(groupProbabilityList, "Probabilities of Proteins");
                pwh.displayInTab();

                if (proteinRatioForMAList.size() > 0)
                {
                    List<Float> logRatios = new ArrayList<Float>();
                    for (Float ratio : proteinRatioForMAList)
                        logRatios.add((float) Math.log(ratio));
                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(proteinSpectralCountForMAList,
                            logRatios, "MAPlot");
                    pwsp.displayInTab();
                }
            }

            if (proteinRatioForMAList.size() > 0)
            {
                int numWithin25Percent=0;
                int numWithin50Percent=0;
                int numWithinTwofold=0;
                int numWithinThreefold=0;
                for (Float ratio : proteinRatioForMAList)
                {
                    Float ratioBigOverSmall = (float) Math.exp(Math.abs(Math.log(ratio)));
                    if (ratioBigOverSmall < 3)
                    {
                        numWithinThreefold++;
                        if (ratioBigOverSmall < 2)
                        {
                            numWithinTwofold++;
                            if (ratioBigOverSmall < 1.5)
                            {
                                numWithin50Percent++;
                                if (ratioBigOverSmall < 1.25)
                                    numWithin25Percent++;
                            }
                        }
                    }
                }
                int numRatios = proteinRatioForMAList.size();
                int percentWithin25Percent=100 * numWithin25Percent / numRatios;
                int percentWithin50Percent=100 * numWithin50Percent / numRatios;
                int percentWithinTwofold=100 * numWithinTwofold / numRatios;
                int percentWithinThreefold=100 * numWithinThreefold / numRatios;
                ApplicationContext.infoMessage("Ratios within 3fold: " + percentWithinThreefold + "%, 2fold: " +
                        percentWithinTwofold + "%, 50%: " + percentWithin50Percent + "%, 25%: " + percentWithin25Percent + "%");

            }

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

            ApplicationContext.infoMessage(protXmlFile.getName() + "\t" +  groupProbabilityList.size() + "\t" + numAbovePoint5 + "\t" +
                                        numAbovePoint75 + "\t" + numAbovePoint9 + "\t" + numAbovePoint95);

            if (showCharts)
            {
                PanelWithHistogram pwhSpecCount = new PanelWithHistogram(spectralCountsAbovePoint9, "Spectral Counts prob .9");
                pwhSpecCount.displayInTab();
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
