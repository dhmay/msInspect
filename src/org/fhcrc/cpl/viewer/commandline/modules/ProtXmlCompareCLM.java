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

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.apache.log4j.Logger;


import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class ProtXmlCompareCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProtXmlCompareCLM.class);

    File[] protXmlFiles;

    protected float minProteinProphet = 0.0f;

    protected boolean listUnique2Proteins = false;

    public ProtXmlCompareCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "protxmlcompare";
        mShortDescription = "protxmlcompare";
        mHelpMessage = "protxmlcompare";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "protxml files to compare"),
                        createDecimalArgumentDefinition("minpprophet", false,
                                "Minimum ProteinProphet score.", minProteinProphet),
                        createBooleanArgumentDefinition("listunique2proteins", false,
                                "List the proteins unique to the second file", listUnique2Proteins)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();

        minProteinProphet = (float) getDoubleArgumentValue("minpprophet");
        if (minProteinProphet > 0)
            ApplicationContext.infoMessage("Only considering proteins with ProteinProphet > " + minProteinProphet);
        listUnique2Proteins = getBooleanArgumentValue("listunique2proteins");
    }


    protected Map<String, Pair<Float, Float>> loadProteinScoreCoverageMap(File protXmlFile)
            throws XMLStreamException, FileNotFoundException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each protein to its highest score
        Map<String, Pair<Float, Float>> proteinScorePercentMap = new HashMap<String, Pair<Float, Float>>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        while (iterator.hasNext())
        {
            ProteinGroup group = iterator.next();   
            float groupScore = group.getProbability();

            if (groupScore < minProteinProphet)
                continue;

            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {
                Set<String> proteinNames = new HashSet<String>();
                //add first protein
                proteinNames.add(protXmlReaderProtein.getProteinName());
                //add indistinguishable proteins
                proteinNames.addAll(protXmlReaderProtein.getIndistinguishableProteinNames());

                for (String proteinName : proteinNames)
                {
                    float scoreThisProteinThisTime = Math.max(protXmlReaderProtein.getProbability(), groupScore);
                    Float coverageThisProteinThisTime = protXmlReaderProtein.getPercentCoverage();
                    if (!proteinName.equals(protXmlReaderProtein.getProteinName()))
                        coverageThisProteinThisTime = null;

                    Pair<Float, Float> scoreAndPercentCoverage = proteinScorePercentMap.get(proteinName);
                    if (scoreAndPercentCoverage == null)
                    {
                        proteinScorePercentMap.put(proteinName,
                                new Pair<Float, Float>(scoreThisProteinThisTime, coverageThisProteinThisTime));
                    }
                    else
                    {
                        scoreAndPercentCoverage.first = Math.max(scoreThisProteinThisTime, scoreAndPercentCoverage.first);
                        scoreAndPercentCoverage.second = scoreAndPercentCoverage == null ? coverageThisProteinThisTime :
                                Math.max(coverageThisProteinThisTime == null ? -1 : coverageThisProteinThisTime, scoreAndPercentCoverage.second);
                    }
//if (proteinName.equals("IPI00021891"))
//{
//    System.err.println("IPI00021891");
//    for (ProtXmlReader.Peptide peptide : protXmlReaderProtein.getPeptides())
//    {
//        System.err.println("\t" + peptide.getPeptideSequence() + ", " + peptide.isContributingEvidence() + ", " + peptide.isNondegenerateEvidence());
//    }
//}
                }
            }
        }
        return proteinScorePercentMap;
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        float goodProbability = 0.9f;
        try
        {
            try
            {
//                PanelWithChart sensSpecChart1 = ProteinUtilities.generateSensSpecChart(protXmlFiles[0]);
                PanelWithChart sensSpecChart =
                        ProteinUtilities.generateSensSpecChart(protXmlFiles[0],protXmlFiles[1]);

                sensSpecChart.setName("Sens/Spec");
//                sensSpecChart2.setName("Sens/Spec 2");

                sensSpecChart.displayInTab();
//                sensSpecChart2.displayInTab();
            }
            catch (Exception e)
            {
                ApplicationContext.setMessage("Failed to generate sens/spec charts.  Error: " + e.getMessage() + ", type: " + e.getClass().getName());
            }




            Map<String, Pair<Float,Float>> proteinScorePercentMap1 =
                    loadProteinScoreCoverageMap(protXmlFiles[0]);
            List<Float> percentCoveragesGoodProb1 = new ArrayList<Float>();
            for (Pair<Float,Float> scorePercent : proteinScorePercentMap1.values())
            {
                if (scorePercent.first > goodProbability && scorePercent.second != null)
                    percentCoveragesGoodProb1.add(scorePercent.second);
            }
            ApplicationContext.infoMessage("Percent coverages (good proteins) 1: 5th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb1, 5))+ "%, 10th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb1, 10)) + "%, median: " +
                    (100 * BasicStatistics.median(percentCoveragesGoodProb1)) + "%");
            ApplicationContext.infoMessage("Min % coverage 1: " + (100 * BasicStatistics.min(percentCoveragesGoodProb1)) +
                    "%, max: " + (100 * BasicStatistics.max(percentCoveragesGoodProb1)));

            Map<String, Pair<Float,Float>> proteinScorePercentMap2 =
                    loadProteinScoreCoverageMap(protXmlFiles[1]);
            List<Float> percentCoveragesGoodProb2 = new ArrayList<Float>();
            for (Pair<Float,Float> scorePercent : proteinScorePercentMap2.values())
            {
                if (scorePercent.first > goodProbability && scorePercent.second != null)
                    percentCoveragesGoodProb2.add(scorePercent.second);
            }
            ApplicationContext.infoMessage("Percent coverages (good proteins) 2: 5th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb2, 5))+ "%, 10th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb2, 10)) + "%, median: " +
                    (100 * BasicStatistics.median(percentCoveragesGoodProb2)) + "%");
            ApplicationContext.infoMessage("Min % coverage 2: " + (100 * BasicStatistics.min(percentCoveragesGoodProb2)) +
                    "%, max: " + (100 * BasicStatistics.max(percentCoveragesGoodProb2)));




            Set<String> commonProteins = new HashSet<String>();

            for (String protein : proteinScorePercentMap1.keySet())
            {
                if (proteinScorePercentMap2.containsKey(protein))
                {
                    commonProteins.add(protein);
                }
            }


            List<Float> scoresInCommon1 = new ArrayList<Float>();
            List<Float> scoresInCommon2 = new ArrayList<Float>();

            List<Float> percentsInCommon1 = new ArrayList<Float>();
            List<Float> percentsInCommon2 = new ArrayList<Float>();
            for (String protein : commonProteins)
            {
                scoresInCommon1.add(proteinScorePercentMap1.get(protein).first);
                scoresInCommon2.add(proteinScorePercentMap2.get(protein).first);

                if (proteinScorePercentMap1.get(protein).second != null &&
                    proteinScorePercentMap2.get(protein).second != null)
//          && proteinScorePercentMap1.get(protein).first >= goodProbability && proteinScorePercentMap2.get(protein).first >= goodProbability)
                {
                    percentsInCommon1.add(proteinScorePercentMap1.get(protein).second);
                    percentsInCommon2.add(proteinScorePercentMap2.get(protein).second);
//if (proteinScorePercentMap2.get(protein).second < proteinScorePercentMap1.get(protein).second)
//    System.err.println("Decreasing % coverage: " + protein);
                }
            }

            ApplicationContext.infoMessage("Proteins in file 1: " + proteinScorePercentMap1.size() +
                    ", file 2: " + proteinScorePercentMap2.size() + ", common: " + commonProteins.size());
//for (String protIn1 : proteinScoreMap1.keySet()) if (!commonProteins.contains(protIn1)) System.err.println("**" + protIn1);
            PanelWithScatterPlot psp = new PanelWithScatterPlot(scoresInCommon1, scoresInCommon2, "protein scores");
            psp.setAxisLabels("Protein Probability, File 1","Protein Probability, File 2");
            psp.setPointSize(2);
            psp.displayInTab();

            List<Float> percentIncreases = new ArrayList<Float>();
            for (int i=0; i<percentsInCommon1.size(); i++)
            {
                float percent1 = percentsInCommon1.get(i);
                float percent2 = percentsInCommon2.get(i);

                if (percent2 > percent1)
                {
                    percentIncreases.add((percent2 - percent1) * 100/ percent1);
                }
            }
            if (!percentIncreases.isEmpty())
            {
                ApplicationContext.infoMessage("Coverage-increasing proteins: " + percentIncreases.size() +
                        " (" + (percentIncreases.size() * 100 / percentsInCommon1.size()) + "%)");
                ApplicationContext.infoMessage("Median % coverage increase: " + BasicStatistics.median(percentIncreases) +
                        "%, 5th percentile: " + BasicStatistics.percentile(percentIncreases,5) +
                        "%, 95th percentile: " + BasicStatistics.percentile(percentIncreases,95) + "%");
                
            }

            List<Float> foldChangesInCommon2 = new ArrayList<Float>();
            for (int i=0; i<percentsInCommon2.size(); i++)
                foldChangesInCommon2.add((percentsInCommon2.get(i) - percentsInCommon1.get(i))/ percentsInCommon1.get(i));
            PanelWithScatterPlot pspCoverage2 = new PanelWithScatterPlot(percentsInCommon1, foldChangesInCommon2, "prot coverage fold change");
            pspCoverage2.setAxisLabels("Protein % Coverage, File 1","Coverage fold change, File 2");
            pspCoverage2.setPointSize(2);
            pspCoverage2.displayInTab();

            List<Float> foldIncreases = new ArrayList<Float>();
            for (float foldChange : foldChangesInCommon2)
            {
                if (foldChange > 0)
                    foldIncreases.add(foldChange);
            }
            if (!foldIncreases.isEmpty())
            {
                ApplicationContext.infoMessage("Median coverage FOLD increase: " + BasicStatistics.median(foldIncreases) +
                        ", 5th percentile: " + BasicStatistics.percentile(foldIncreases,5) +
                        ", 95th percentile: " + BasicStatistics.percentile(foldIncreases,95));

            }


            List<Float> scoresIn2ButNot1 = new ArrayList<Float>();
            List<Float> coveragesIn2ButNot1 = new ArrayList<Float>();
            Set<String> proteinsIn2ButNot1 = new HashSet<String>();

            for (String protein : proteinScorePercentMap2.keySet())
            {
                if (!commonProteins.contains(protein))
                {
                    proteinsIn2ButNot1.add(protein);
                    scoresIn2ButNot1.add(proteinScorePercentMap2.get(protein).first);
                    if (proteinScorePercentMap2.get(protein).second != null)
                        coveragesIn2ButNot1.add(proteinScorePercentMap2.get(protein).second);
                }
            }
            if (scoresIn2ButNot1.size() > 0)
            {
                PanelWithHistogram pwh = new PanelWithHistogram(scoresIn2ButNot1, "scores in 2 but not 1");
                pwh.displayInTab();

                if (!coveragesIn2ButNot1.isEmpty())
                {
                    PanelWithHistogram pwh2 = new PanelWithHistogram(coveragesIn2ButNot1, "coverage in 2 but not 1");
                    pwh2.displayInTab();
                }
            }


            int scoresAbovePoint751 = 0;
            int scoresAbovePoint752 = 0;

            List<Float> scores1 = new ArrayList<Float>();
            List<Float> scores2 = new ArrayList<Float>();


            Set<String> proteinsAbovePoint951 = new HashSet<String>();
            Set<String> proteinsAbovePoint952 = new HashSet<String>();



            for (String protein  : proteinScorePercentMap1.keySet())
            {
                float score = proteinScorePercentMap1.get(protein).first;
                if (score >= .75f)
                    scoresAbovePoint751++;
                if (score >= goodProbability)
                    proteinsAbovePoint951.add(protein);
                scores1.add(score);
            }

            for (String protein  : proteinScorePercentMap2.keySet())
            {
                float score = proteinScorePercentMap2.get(protein).first;
                if (score >= .75f)
                    scoresAbovePoint752++;
                if (score >= goodProbability)
                    proteinsAbovePoint952.add(protein);
                scores2.add(score);
            }

            Set<String> commonProteinsAbovePoint95 = new HashSet<String>();
            for (String protein : proteinsAbovePoint951)
            {
                if (proteinsAbovePoint952.contains(protein))
                    commonProteinsAbovePoint95.add(protein);
            }

            ApplicationContext.infoMessage("Scores above " + goodProbability + ": Set 1 = " + proteinsAbovePoint951.size() +
                                           ", Set 2 = " + proteinsAbovePoint952.size() + ", common=" +
                                           commonProteinsAbovePoint95.size() + ", unique to 2: " +
                                           ( proteinsAbovePoint952.size() - commonProteinsAbovePoint95.size()));
            ApplicationContext.infoMessage("Scores above .75: Set 1 = " + scoresAbovePoint751 +
                                           ", Set 2 = " + scoresAbovePoint752);

            if (listUnique2Proteins)
            {
                ApplicationContext.infoMessage("Proteins unique to second file:");
                if (proteinsIn2ButNot1.isEmpty())
                    ApplicationContext.infoMessage("\tNone.");
                for (String prot2 : proteinsIn2ButNot1)
                    ApplicationContext.infoMessage("\t" + prot2);


            }

            PanelWithHistogram pwh2 = new PanelWithHistogram(scores1, "Scores in 1");
            pwh2.displayInTab();
            PanelWithHistogram pwh3 = new PanelWithHistogram(scores2, "Scores in 2");
            pwh3.displayInTab();


            List<ProtXmlReader.Protein> proteinList1 =
                    ProteinUtilities.loadProtXmlProteinsFromProtXML(protXmlFiles[0]);
            List<ProtXmlReader.Protein> proteinList2 =
                    ProteinUtilities.loadProtXmlProteinsFromProtXML(protXmlFiles[1]);
            List<Float> allLogRatios1 = new ArrayList<Float>();
            Set<String> point95RatioProteins1 = new HashSet<String>();



            List<ProteinGroup> proteinGroups1 = ProteinUtilities.loadProteinGroupsFromProtXML(protXmlFiles[0]);
            List<ProteinGroup> proteinGroups2 = ProteinUtilities.loadProteinGroupsFromProtXML(protXmlFiles[1]);

            List<ProteinGroup> commonProtGroupsByIPI = new ArrayList<ProteinGroup>();
            int numGroupsAbovePoint951 = 0;
            int commonGroupsAbovePoint95 = 0;

            int numQuantGroupsAbovePoint951 = 0;

            List<Float> numPeptidesPerGroup1 = new ArrayList<Float>();

            for (ProteinGroup proteinGroup1 : proteinGroups1)
            {
                boolean quantified = false;
                boolean abovePoint95 = proteinGroup1.getProbability() >= goodProbability;
                if (abovePoint95)
                    numGroupsAbovePoint951++;
                for (ProtXmlReader.Protein protein : proteinGroup1.getProteins())
                {
                    int numPeptides =  protein.getUniquePeptidesCount();
                    if (numPeptides < 50)
                        numPeptidesPerGroup1.add((float) numPeptides);

                    if (protein.getQuantitationRatio() != null)
                        quantified=true;
                }
                if (quantified && abovePoint95)
                    numQuantGroupsAbovePoint951++;
            }





            if (true)
            {
                List<Float> ratios1 = new ArrayList<Float>();
                List<Float> ratios2 = new ArrayList<Float>();

                List<Double> numRatioPeptides1 = new ArrayList<Double>();

                List<Double> numPeptides1 = new ArrayList<Double>();

//                Set<ProtXmlReader.Protein> proteinsWithRatios2 = new HashSet<ProtXmlReader.Protein>();
                Map<String, ProtXmlReader.Protein> ratioIdProteinMap1 =
                        new HashMap<String, ProtXmlReader.Protein>();

                for (ProtXmlReader.Protein protein1 : proteinList1)
                {
                    if (protein1.getProbability() < minProteinProphet)
                        continue;
                    ProtXmlReader.QuantitationRatio ratio1 = protein1.getQuantitationRatio();
                    numPeptides1.add((double) protein1.getPeptides().size());
                    if (ratio1 != null)
                    {
                        float logRatio1 = (float) Math.log(ratio1.getRatioMean());
                        if (logRatio1 < -10f) logRatio1 = -10f;
                        if (logRatio1 > 10f) logRatio1 = 10f;
                        allLogRatios1.add(logRatio1);
                        if (protein1.getProbability() > goodProbability)
                            point95RatioProteins1.add(protein1.getProteinName());
                        numRatioPeptides1.add((double) ratio1.getRatioNumberPeptides());
                        ratioIdProteinMap1.put(protein1.getProteinName(), protein1);
                    }

                }

                List<Float> allLogRatios2 = new ArrayList<Float>();
                Set<String> point95RatioProteins2 = new HashSet<String>();

                List<Double> numRatioPeptides2 = new ArrayList<Double>();
                List<Double> numPeptides2 = new ArrayList<Double>();
                Map<String, ProtXmlReader.Protein> ratioIdProteinMap2 =
                        new HashMap<String, ProtXmlReader.Protein>();
                for (ProtXmlReader.Protein protein2 : proteinList2)
                {
                    if (protein2.getProbability() < minProteinProphet)
                        continue;
                    ProtXmlReader.QuantitationRatio ratio2 = protein2.getQuantitationRatio();
                    numPeptides2.add((double) protein2.getPeptides().size());
                    if (ratio2 != null)
                    {
                        float logRatio2 = (float) Math.log(ratio2.getRatioMean());
                        if (logRatio2 < -10f) logRatio2 = -10f;
                        if (logRatio2 > 10f) logRatio2 = 10f;
                        allLogRatios2.add(logRatio2);
                        if (protein2.getProbability() > goodProbability)
                            point95RatioProteins2.add(protein2.getProteinName());
                        numRatioPeptides2.add((double) ratio2.getRatioNumberPeptides());
                        ratioIdProteinMap2.put(protein2.getProteinName(), protein2);
                    }

                }

                List<Float> scoresWithRatios1 = new ArrayList<Float>();
                List<Float> scoresWithRatios2 = new ArrayList<Float>();

                for (String proteinName : ratioIdProteinMap1.keySet())
                {

                    ProtXmlReader.Protein protein1 =  ratioIdProteinMap1.get(proteinName);
//                    if (protein1.getProbability() < 0.9) continue;
                    double ratio1 = protein1.getQuantitationRatio().getRatioMean();
                    scoresWithRatios1.add(protein1.getProbability());


                    ProtXmlReader.Protein protein2 = ratioIdProteinMap2.get(proteinName);
                    if (protein2 != null)
                    {
//                    if (protein2.getProbability() < 0.9) continue;

                        double ratio2 = protein2.getQuantitationRatio().getRatioMean();
                        ratios1.add((float)ratio1);
                        ratios2.add((float)ratio2);

                        scoresWithRatios2.add(protein2.getProbability());

                    }
                    else scoresWithRatios2.add(0.0f);
                }

                for (String proteinName : ratioIdProteinMap2.keySet())
                {

                    if (!ratioIdProteinMap1.containsKey(proteinName))
                    {
                        ProtXmlReader.Protein protein2 = ratioIdProteinMap2.get(proteinName);
                        scoresWithRatios2.add(protein2.getProbability());
                        scoresWithRatios1.add(0.0f);
                    }
                }


                PanelWithScatterPlot probsWithRatiosSP = new PanelWithScatterPlot(scoresWithRatios1, scoresWithRatios2, "Probs w/ratios");
                probsWithRatiosSP.displayInTab();

                Set<String> commonConfidentRatioProteinNames = new HashSet<String>();
                for (String proteinName : point95RatioProteins1)
                    if (point95RatioProteins2.contains(proteinName))
                         commonConfidentRatioProteinNames.add(proteinName);


                ApplicationContext.setMessage("Ratios in file 1: " + allLogRatios1.size() +
                        ", file 2: " + allLogRatios2.size());
                ApplicationContext.setMessage("Ratios with prob>" + goodProbability + ", file 1: " + point95RatioProteins1.size() +
                        ", file 2: " + point95RatioProteins2.size() + ", common=" +
                        commonConfidentRatioProteinNames.size() + ", unique to 2: " +
                        (point95RatioProteins2.size() - commonConfidentRatioProteinNames.size()));
                ApplicationContext.setMessage("Mean peptides, file 1: " + BasicStatistics.mean(numPeptides1) +
                         ", file2: " + BasicStatistics.mean(numPeptides2));                
                ApplicationContext.setMessage("Mean ratio peptides, file 1: " + BasicStatistics.mean(numRatioPeptides1) +
                         ", file2: " + BasicStatistics.mean(numRatioPeptides2));

                if (ratios1.size() >= 1)
                {

                    PanelWithScatterPlot ratioSP = new PanelWithScatterPlot(ratios1, ratios2, "Quantitation ratios");
                    ratioSP.displayInTab();
                    System.err.println("Ratios to compare: " + ratios1.size());

                    double[] logRatios1 = new double[ratios1.size()];
                    double[] logRatios2 = new double[ratios2.size()];

                    for (int i=0; i<ratios1.size(); i++)
                    {
                        logRatios1[i] = Math.log(Math.max(ratios1.get(i),.001));
                        logRatios2[i] = Math.log(Math.max(ratios2.get(i),.001));
                    }

                    PanelWithScatterPlot logRatioSP = new PanelWithScatterPlot(logRatios1, logRatios2, "log Quantitation ratios");
                    logRatioSP.displayInTab();
                }
                else
                    ApplicationContext.setMessage("No ratios to compare on plot");

                if (allLogRatios1.size() >=1)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(allLogRatios1, "Log Ratios 1");
                    pwh.displayInTab();
                }
               if (allLogRatios2.size() >=1)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(allLogRatios2, "Log Ratios 2");
                    pwh.displayInTab();
                }
            }





            int numGroupsAbovePoint952 = 0;
            int numQuantGroupsAbovePoint952 = 0;
            int commonQuantGroupsAbovePoint95 = 0;


            List<Float> numPeptidesPerGroup2 = new ArrayList<Float>();

            for (ProteinGroup proteinGroup2 : proteinGroups2)
            {
                boolean abovePoint95 = false;
                if (proteinGroup2.getProbability() >= goodProbability)
                    abovePoint95 = true;
                if (abovePoint95)
                    numGroupsAbovePoint952++;
                Collection<String> proteinNames1 = proteinScorePercentMap1.keySet();
                boolean inCommon = false;
                boolean inCommonAbovePoint95 = false;
                boolean inCommonQuantAbovePoint95 = false;


                boolean quantified = false;
                for (ProtXmlReader.Protein protein : proteinGroup2.getProteins())
                {
                    int numPeptides =  protein.getUniquePeptidesCount();
                    if (numPeptides < 50)
                        numPeptidesPerGroup2.add((float) numPeptides);
                    
                    if (proteinNames1.contains(protein.getProteinName()))
                        inCommon = true;

                    if (abovePoint95 && proteinsAbovePoint951.contains(protein.getProteinName()))
                        inCommonAbovePoint95 = true;
                    if (protein.getQuantitationRatio() != null)
                        quantified=true;
                    if (abovePoint95 && quantified && point95RatioProteins1.contains(protein.getProteinName()))
                        inCommonQuantAbovePoint95 = true;
                }
                if (inCommon) commonProtGroupsByIPI.add(proteinGroup2);
                if (inCommonAbovePoint95) commonGroupsAbovePoint95++;

                if (quantified && abovePoint95)
                {
                    numQuantGroupsAbovePoint952++;
                    if (inCommonQuantAbovePoint95)
                        commonQuantGroupsAbovePoint95++;
                }
            }


            ApplicationContext.infoMessage("ProteinGroups: Set 1 = " + proteinGroups1.size() +
                    ", Set 2 = " + proteinGroups2.size() + ", common (by IPI) = " + commonProtGroupsByIPI.size() +
                    ", Unique to 2: " + (proteinGroups2.size() - commonProtGroupsByIPI.size()));
            ApplicationContext.infoMessage("ProteinGroups above " + goodProbability + ": Set 1 = " + numGroupsAbovePoint951 +
                    ", Set 2 = " + numGroupsAbovePoint952 + ", common (by IPI) = " + commonGroupsAbovePoint95 +
                    ", Unique to 2: " + (numGroupsAbovePoint952 - commonGroupsAbovePoint95));
            ApplicationContext.infoMessage("Quantified ProteinGroups above " + goodProbability + ": Set 1 = " + numQuantGroupsAbovePoint951 +
                    ", Set 2 = " + numQuantGroupsAbovePoint952 + ", common (by IPI) = " + commonQuantGroupsAbovePoint95 +
                    ", Unique to 2: " + (numQuantGroupsAbovePoint952 - commonQuantGroupsAbovePoint95));

            PanelWithHistogram peptidesPerGroupHist1 = new PanelWithHistogram(numPeptidesPerGroup1, "#Peptides1");
            peptidesPerGroupHist1.displayInTab();

            PanelWithHistogram peptidesPerGroupHist2 = new PanelWithHistogram(numPeptidesPerGroup2, "#Peptides2");
            peptidesPerGroupHist2.displayInTab();


//            Map<String, Set<String>> peptideGoodProtMap1 =
//                ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFiles[0], goodProbability);
//            Map<String, Set<String>> peptideGoodProtMap2 =
//                ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFiles[0], goodProbability);

//            int numPeptides2NewProts = 0;
//            int numPeptides2ExistingProts = 0;
//            Set<String> peptides1 = new HashSet<String>();
//            for (String peptide : peptideGoodProtMap2.keySet())
//            {
//                if (peptideGoodProtMap1.containsKey(peptide))
//                    numPeptides2ExistingProts++;
//                else
//                    numPeptides2NewProts++;
//            }
//            ApplicationContext.infoMessage("");

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
