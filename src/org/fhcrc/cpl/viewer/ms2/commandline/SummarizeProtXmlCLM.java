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

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
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
public class SummarizeProtXmlCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SummarizeProtXmlCLM.class);

    protected File protXmlFiles[];
    protected boolean showCharts = false;
    protected boolean listProteins = false;
    protected float minProteinProphet = 0;
    protected File protGeneFile;
    protected Map<String,List<String>> protGenesMap = null;

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
                        createUnnamedSeriesFileArgumentDefinition(true, "protxml file(s)"),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new DecimalArgumentDefinition("minpprophet", false, "Min proteinprophet for MA plot",
                                minProteinProphet),
                        new FileToReadArgumentDefinition("protgenefile", false,
                                "File associating gene symbols with protein accession numbers"),
                        new BooleanArgumentDefinition("listproteins", false, "List proteins to stderr?",
                                listProteins),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();
        showCharts = getBooleanArgumentValue("showcharts");
        listProteins = getBooleanArgumentValue("listproteins");
        minProteinProphet = getFloatArgumentValue("minpprophet");
        protGeneFile = getFileArgumentValue("protgenefile");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("File\tGroups\tPoint1\tPoint5\tPoint75\tPoint9\tPoint95");
        if (protGeneFile != null)
        {
            try
            {
                ApplicationContext.infoMessage("Loading IPI-gene map...");
                protGenesMap = QAUtilities.loadIpiGeneListMap(protGeneFile);
                ApplicationContext.infoMessage("Done");

            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to load protein-gene map file",e);
            }
        }
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

            List<Float> proteinRatioVarianceList = new ArrayList<Float>();
            List<Float> proteinWithRatioProbabilityList = new ArrayList<Float>();


            List<Float> proteinSpectralCountForMAList = new ArrayList<Float>();
            List<Float> proteinRatioNumPeptidesList = new ArrayList<Float>();



            Iterator<ProteinGroup> iterator = protXmlReader.iterator();

            Set<Integer> groupsWithProbabilityGreaterThanPoint1 = new HashSet<Integer>();
            Set<Integer> groupsWith2PlusPeptides = new HashSet<Integer>();
            Set<Integer> groupsWithZeroOrOneGenes = new HashSet<Integer>();

            while (iterator.hasNext())
            {
                ProteinGroup pg = iterator.next();
                Set<String> genesThisGroup = new HashSet<String>();

                Set<String> allPeptidesThisGroup = new HashSet<String>();

                groupProbabilityList.add(pg.getProbability());

                for (ProtXmlReader.Protein protein : pg.getProteins())
                {
                    proteinProbabilityList.add(protein.getProbability());
                    proteinSpectralCountList.add((float) protein.getTotalNumberPeptides());
                    for (ProtXmlReader.Peptide peptide : protein.getPeptides())
                    {
                        allPeptidesThisGroup.add(peptide.getPeptideSequence());
                    }
                    if (protein.getProbability() >= minProteinProphet && protein.getQuantitationRatio() != null)
                    {
                        proteinRatioForMAList.add(protein.getQuantitationRatio().getRatioMean());
                        proteinSpectralCountForMAList.add((float) protein.getTotalNumberPeptides());
                        proteinRatioVarianceList.add(protein.getQuantitationRatio().getRatioStandardDev());
                        proteinRatioNumPeptidesList.add((float) protein.getQuantitationRatio().getRatioNumberPeptides());

                        proteinWithRatioProbabilityList.add(pg.getProbability());
                    }
                    List<String> genesThisProtein = null;
                    if (protGenesMap != null && protGenesMap.containsKey(protein.getProteinName()))
                    {
                        genesThisProtein = protGenesMap.get(protein.getProteinName());
                        genesThisGroup.addAll(genesThisProtein);
                    }
                    if (listProteins)
                    {
                        String genesString = "";
                        if (genesThisProtein != null)
                        {
                            StringBuffer genesStringBuf = new StringBuffer();
                            for (int i=0; i<genesThisProtein.size(); i++)
                            {
                                if (i>0) genesStringBuf.append(",");
                                genesStringBuf.append(genesThisProtein.get(i));
                            }
                            genesString = genesStringBuf.toString();
                        }
                        System.err.println(protein.getProteinName() + "\t" + genesString + "\t" +
                                protein.getProbability() + "\t" +
                                protein.getPeptides().size());
                    }
                }
                if (allPeptidesThisGroup.size() > 1)
                    groupsWith2PlusPeptides.add(pg.getGroupNumber());
                if (pg.getGroupProbability() >= 0.1)
                    groupsWithProbabilityGreaterThanPoint1.add(pg.getGroupNumber());
                if (protGenesMap != null && genesThisGroup.size() < 2)
                    groupsWithZeroOrOneGenes.add(pg.getGroupNumber());
            }
            ApplicationContext.infoMessage("Protein Groups with probability >=0.1: " +
                    groupsWithProbabilityGreaterThanPoint1.size());
            ApplicationContext.infoMessage("Protein Groups with probability >=0.1 and 2 or more peptides: " +
                    setIntersection(groupsWith2PlusPeptides, groupsWithProbabilityGreaterThanPoint1).size());
            if (protGenesMap != null)
            {
                ApplicationContext.infoMessage("Protein Groups with probability >=0.1 and mapping to zero or one genes: " +
                        setIntersection(groupsWithZeroOrOneGenes, groupsWithProbabilityGreaterThanPoint1).size());
                ApplicationContext.infoMessage("Protein Groups with probability >=0.1 and mapping to zero or one genes " +
                        " and 2 or more peptides: " +
                        setIntersection(groupsWith2PlusPeptides, setIntersection(groupsWithZeroOrOneGenes, groupsWithProbabilityGreaterThanPoint1)).size());                
            }

            if (showCharts)
            {
                PanelWithHistogram pwh = new PanelWithHistogram(groupProbabilityList, "Probabilities of Proteins");
                pwh.displayInTab();

                if (proteinRatioForMAList.size() > 0)
                {
                    List<Float> logRatios = new ArrayList<Float>();
                    List<Float> logSpectralCounts = new ArrayList<Float>();
                    List<Float> logVariance = new ArrayList<Float>();


                    for (Float ratio : proteinRatioForMAList)
                        logRatios.add((float) Math.log(Math.max(ratio,0.000001)));
                    for (Float count : proteinSpectralCountForMAList)
                        logSpectralCounts.add((float) Math.log(Math.max(1, count)));
                    for (Float var : proteinRatioVarianceList)
                        logVariance.add((float) Math.log(Math.max(var,0.000001)));
                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(logSpectralCounts,
                            logRatios, "MAPlot");
                    pwsp.displayInTab();

                    PanelWithHistogram pwhRatios = new PanelWithHistogram(logRatios, "Log Ratios");
                    pwhRatios.displayInTab();

                    PanelWithScatterPlot pwsp2 = new PanelWithScatterPlot(logRatios, logVariance,
                            "LogRatio vs LogVariance");
                    pwsp2.setAxisLabels("Log Ratio","Log Ratio Variance");
                    pwsp2.displayInTab();

                    PanelWithScatterPlot pwsp3 = new PanelWithScatterPlot(proteinWithRatioProbabilityList, logRatios, 
                            "Probability vs LogRatio");
                    pwsp3.setAxisLabels("Probability","Log Ratio");
                    pwsp3.displayInTab();

                    PanelWithScatterPlot pwsp4 = new PanelWithScatterPlot(proteinRatioNumPeptidesList, proteinRatioVarianceList, 
                            "Quant Events vs Variance");
                    pwsp4.setAxisLabels("Quant Events","Variance");
                    pwsp4.displayInTab();

                    PanelWithScatterPlot pwsp5 = new PanelWithScatterPlot(proteinRatioNumPeptidesList, logRatios, 
                            "Quant Events vs LogRatio");
                    pwsp5.setAxisLabels("Quant Events","Log Ratio");
                    pwsp5.displayInTab();

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
                ApplicationContext.infoMessage("Quantitated proteins: " + numRatios);
                ApplicationContext.infoMessage("Ratios within 3fold: " + percentWithinThreefold + "%, 2fold: " +
                        percentWithinTwofold + "%, 50%: " + percentWithin50Percent + "%, 25%: " + percentWithin25Percent + "%");

            }

            int numAbovePoint1 = 0;
            int numAbovePoint5 = 0;
            int numAbovePoint75 = 0;
            int numAbovePoint9 = 0;
            int numAbovePoint95 = 0;

            for (int i=0; i<groupProbabilityList.size(); i++)
            {
                float probability = groupProbabilityList.get(i);
                if (probability >= .5f)
                {
                    numAbovePoint1++;
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
            }

            List<Float> spectralCountsAbovePoint9 = new ArrayList<Float>();

            for (int i=0; i<proteinProbabilityList.size(); i++)
            {
                float probability = proteinProbabilityList.get(i);
                if (probability >= 0.9f)
                    spectralCountsAbovePoint9.add(proteinSpectralCountList.get(i));

            }

            ApplicationContext.infoMessage(protXmlFile.getName() + "\t" +  groupProbabilityList.size()+ "\t" + numAbovePoint1 + "\t" + numAbovePoint5 + "\t" +
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

    public static <T> Set<T> setIntersection(Set<? extends T> a, Set<? extends T> b)
    {
        Set result = new HashSet<T>(a);
        result.retainAll(b);
        return result;
    }


}
