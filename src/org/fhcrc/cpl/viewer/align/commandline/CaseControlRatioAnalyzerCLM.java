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
package org.fhcrc.cpl.viewer.align.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBoxAndWhiskerChart;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule for feature finding
 */
public class CaseControlRatioAnalyzerCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CaseControlRatioAnalyzerCLM.class);

    protected FeatureSet pepXmlFeatureSet;
    protected File protXmlFile;
    protected File pepArrayFile;
    List<String> caseRunNameList;
    List<String> controlRunNameList;

    protected String proteinToAnalyze;

    protected boolean showCharts = false;

    Object[] arrayRows;


    public CaseControlRatioAnalyzerCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "analyzecasecontrolratio";

        mShortDescription = "analyzecasecontrolratio";
        mHelpMessage = "analyzecasecontrolratio";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("protxml",true,"protxml"),
                        new FeatureFileArgumentDefinition("pepxml",true,"pepxml"),
                        new FileToReadArgumentDefinition("peparray",true,"peptide array"),

                        new StringArgumentDefinition("protein", true, "Protein to analyze"),

                        new FileToReadArgumentDefinition("caserunlistfile",true,"File containing the names of runs in the case group, one per line"),
                            new FileToReadArgumentDefinition("controlrunlistfile",true,"File containing the names of runs in the control group, one per line"),

                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        proteinToAnalyze = getStringArgumentValue("protein");

        showCharts = getBooleanArgumentValue("showcharts");

        protXmlFile = getFileArgumentValue("protXml");
        pepArrayFile = getFileArgumentValue("peparray");

        pepXmlFeatureSet = getFeatureSetArgumentValue("pepxml");



        File caseRunListFile = getFileArgumentValue("caserunlistfile");
        File controlRunListFile = getFileArgumentValue("controlrunlistfile");



        try
        {
            BufferedReader br = new BufferedReader(new FileReader(caseRunListFile));
            caseRunNameList = new ArrayList<String>();
            while (true)
            {
                String line = br.readLine();
                if (null == line)
                    break;
                if (line.length() == 0 || line.charAt(0) == '#')
                    continue;
                caseRunNameList.add(line);
            }

            br = new BufferedReader(new FileReader(controlRunListFile));
            controlRunNameList = new ArrayList<String>();
            while (true)
            {
                String line = br.readLine();
                if (null == line)
                    break;
                if (line.length() == 0 || line.charAt(0) == '#')
                    continue;
                controlRunNameList.add(line);
            }

            ApplicationContext.infoMessage("Case runs:");
            for (String caseRun : caseRunNameList)
                ApplicationContext.infoMessage("\t" + caseRun);
            ApplicationContext.infoMessage("Control runs:");
            for (String controlRun : controlRunNameList)
                ApplicationContext.infoMessage("\t" + controlRun);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PeptideArrayAnalyzer pepArrayAnalyzer  = null;
        Map<String, List<Map<String, Object>>> peptideRowsMap = new HashMap<String, List<Map<String, Object>>>();
        try
        {
            pepArrayAnalyzer = new PeptideArrayAnalyzer(pepArrayFile);
            arrayRows = pepArrayAnalyzer.getRowMaps();

                for (Object rowMapObj : arrayRows)
                {
                    Map<String, Object> rowMap = (Map<String, Object>) rowMapObj;
                    Set<String> rowPeptides = pepArrayAnalyzer.getAllRowPeptides(rowMap);
                    if (rowPeptides.size() == 1)
                    {
                        String peptide = rowPeptides.iterator().next();
                        List<Map<String, Object>> rowsThisPeptide = peptideRowsMap.get(peptide);
                        if (rowsThisPeptide == null)
                        {
                            rowsThisPeptide = new ArrayList<Map<String, Object>>();
                            peptideRowsMap.put(peptide, rowsThisPeptide);
                        }
                        rowsThisPeptide.add(rowMap);
                    }
                }
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to open peptide array file");
        }
        try
        {
            ProtXmlReader.Protein protXmlProtein = ProteinUtilities.loadFirstProteinOccurrence(protXmlFile, proteinToAnalyze);
            float proteinRatio = protXmlProtein.getQuantitationRatio().getRatioMean();
            ApplicationContext.infoMessage("Protein: " + protXmlProtein.getProteinName() + ", prob: " +
                    protXmlProtein.getProbability() + ", ratio: " + proteinRatio);
            Map<String, Set<String>> proteinPeptideMap = ProteinUtilities.loadProteinPeptideSequenceMapFromProtXML(protXmlFile, 0f);
            Set<String> proteinPeptides =proteinPeptideMap.get(proteinToAnalyze);

            Map<String, List<Pair<List<Float>, List<Float>>>> peptideCaseControlIntensitiesMap =
                    new HashMap<String, List<Pair<List<Float>, List<Float>>>>();
            Map<String, List<Float>> peptideLogRatiosMap = new HashMap<String, List<Float>>();

            for (Feature feature : pepXmlFeatureSet.getFeatures())
            {
                if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                    continue;
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (!proteinPeptides.contains(peptide))
                    continue;
                float peptideRatio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                List<Float> logRatiosThisPeptide = peptideLogRatiosMap.get(peptide);
                if (logRatiosThisPeptide == null)
                {
                    logRatiosThisPeptide = new ArrayList<Float>();
                    peptideLogRatiosMap.put(peptide, logRatiosThisPeptide);

                }
                logRatiosThisPeptide.add((float) Math.log(peptideRatio));


            }


            PanelWithScatterPlot pwspPeptideLogRatios = new PanelWithScatterPlot();
            pwspPeptideLogRatios.setName("log-ratios by peptide");
            pwspPeptideLogRatios.setPointSize(4);
            pwspPeptideLogRatios.getChart().getXYPlot().setDomainCrosshairValue(Math.log(proteinRatio), true);
            pwspPeptideLogRatios.getChart().getXYPlot().setDomainCrosshairVisible(true);
            int peptideIndex = 0;
            float maxAbsValue = -5000;
            for (String peptide : peptideLogRatiosMap.keySet())
            {
                List<Float> logRatiosThisPeptide = peptideLogRatiosMap.get(peptide);
                float meanRatioThisPeptide = (float) Math.exp(BasicStatistics.mean(logRatiosThisPeptide));
                ApplicationContext.infoMessage("Peptide " + peptide + ": geo mean ratio: " + meanRatioThisPeptide + ", observations: " + logRatiosThisPeptide.size());

                List<Float> yValues = new ArrayList<Float>();
                for (float blah : logRatiosThisPeptide)
                    yValues.add((float) peptideIndex);

                pwspPeptideLogRatios.addData(logRatiosThisPeptide, yValues, peptide);

                for (float value : logRatiosThisPeptide)
                    maxAbsValue = Math.max(maxAbsValue, Math.abs(value));

                peptideIndex++;
            }

            float maxChartAxis = (float) Math.max(Math.log(10), maxAbsValue);
            float minChartAxis = -maxChartAxis;
            pwspPeptideLogRatios.getChart().getXYPlot().getDomainAxis().setLowerBound(minChartAxis);
            pwspPeptideLogRatios.getChart().getXYPlot().getDomainAxis().setUpperBound(maxChartAxis);
            pwspPeptideLogRatios.displayInTab();




            //todo: cut this out when moving to ProtXmlProteinInfoCLM
            for (String peptide : peptideLogRatiosMap.keySet())
            {
                List<Pair<List<Float>, List<Float>>> caseControlIntensitiesThisPeptide = new ArrayList<Pair<List<Float>, List<Float>>>();
                peptideCaseControlIntensitiesMap.put(peptide, caseControlIntensitiesThisPeptide);
                for (Map<String, Object> rowMap : peptideRowsMap.get(peptide))
                {
                    List<Float> caseIntensities = new ArrayList<Float>();
                    for (String run : caseRunNameList)
                    {
                        if (rowMap.containsKey("intensity_" + run))
                            caseIntensities.add(Float.parseFloat(rowMap.get("intensity_" + run).toString()));
                    }
                    List<Float> controlIntensities = new ArrayList<Float>();
                    for (String run : controlRunNameList)
                    {
                        if (rowMap.containsKey("intensity_" + run))
                            controlIntensities.add(Float.parseFloat(rowMap.get("intensity_" + run).toString()));
                    }

                    if (caseIntensities.size() > 1 && controlIntensities.size() > 1)
                    {

                        Pair<List<Float>, List<Float>> caseControlIntensities = new Pair<List<Float>, List<Float>>(caseIntensities, controlIntensities);
                        caseControlIntensitiesThisPeptide.add(caseControlIntensities);

                        PanelWithBoxAndWhiskerChart intensityBoxPlots = new PanelWithBoxAndWhiskerChart("intensities");
                        intensityBoxPlots.addData(caseIntensities, "case");
                        intensityBoxPlots.addData(controlIntensities, "control");

                        intensityBoxPlots.displayInTab();
                    }
                }
            }
            //todo: end cut out part
            
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
