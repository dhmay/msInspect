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

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
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
public class ProtXmlCompareCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProtXmlCompareCLM.class);

    File[] protXmlFiles;

    protected float minProteinProphet = 0.0f;

    protected boolean listUnique2Proteins = false;

    protected boolean showCharts = true;

    protected File outSpecCountFile;
    protected File outProteinRatioFile;


    public ProtXmlCompareCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "protxmlcompare";
        mShortDescription = "Compare two protXML files";
        mHelpMessage = "Compare two protXML files -- compare probabilities for each protein, " +
                "show protein overlap, etc";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "protxml files to compare"),
                        new DecimalArgumentDefinition("minpprophet", false,
                                "Minimum ProteinProphet score.", minProteinProphet),
                        new BooleanArgumentDefinition("listunique2proteins", false,
                                "List the proteins unique to the second file", listUnique2Proteins),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new FileToWriteArgumentDefinition("outspeccountfile", false, "output file for spectral count comparison"),
                        new FileToWriteArgumentDefinition("outproteinratiofile", false, "output file for comparing protein ratios"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();

        minProteinProphet = getFloatArgumentValue("minpprophet");
        if (minProteinProphet > 0)
            ApplicationContext.infoMessage("Only considering proteins with ProteinProphet > " + minProteinProphet);
        listUnique2Proteins = getBooleanArgumentValue("listunique2proteins");
        showCharts = getBooleanArgumentValue("showcharts");
        outSpecCountFile = getFileArgumentValue("outspeccountfile");
        outProteinRatioFile = getFileArgumentValue("outproteinratiofile");

    }

    protected static class ProteinInfo
    {
        protected float score;
        protected float coverage;
        protected int spectralCount;
        protected int numUniquePeptides;

        public ProteinInfo(float score, float coverage, int spectralCount, int numUniquePeptides)
        {
            this.score = score;
            this.coverage = coverage;
            this.spectralCount = spectralCount;
            this.numUniquePeptides = numUniquePeptides;
        }

        public float getScore()
        {
            return score;
        }

        public float getCoverage()
        {
            return coverage;
        }

        public int getSpectralCount()
        {
            return spectralCount;
        }

        public void setScore(float score)
        {
            this.score = score;
        }

        public void setCoverage(float coverage)
        {
            this.coverage = coverage;
        }

        public void setSpectralCount(int spectralCount)
        {
            this.spectralCount = spectralCount;
        }

        public int getNumUniquePeptides()
        {
            return numUniquePeptides;
        }

        public void setNumUniquePeptides(int numUniquePeptides)
        {
            this.numUniquePeptides = numUniquePeptides;
        }
    }


    protected Map<String, ProteinInfo> loadProteinInfoMap(File protXmlFile)
            throws XMLStreamException, FileNotFoundException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each protein to its highest score
        Map<String, ProteinInfo> proteinInfoMap = new HashMap<String, ProteinInfo>();

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
                    int spectralCountThisProteinThisTime = protXmlReaderProtein.getTotalNumberPeptides();
                    int numUniquePeptidesThisProteinThisTime = protXmlReaderProtein.getUniquePeptidesCount();
                    if (!proteinName.equals(protXmlReaderProtein.getProteinName()))
                        coverageThisProteinThisTime = null;
                    ProteinInfo proteinInfo = proteinInfoMap.get(proteinName);
                    if (proteinInfo == null)
                    {
                        proteinInfoMap.put(proteinName,
                                new ProteinInfo(scoreThisProteinThisTime,
                                        coverageThisProteinThisTime == null ? 0 : coverageThisProteinThisTime,
                                        spectralCountThisProteinThisTime, numUniquePeptidesThisProteinThisTime));
                    }
                    else
                    {
                        proteinInfo.setScore(Math.max(scoreThisProteinThisTime, proteinInfo.getScore()));
                        proteinInfo.setCoverage(Math.max(coverageThisProteinThisTime == null ? -1 : coverageThisProteinThisTime, proteinInfo.getCoverage()));
                        proteinInfo.setSpectralCount(Math.max(spectralCountThisProteinThisTime, proteinInfo.getSpectralCount()));
                        proteinInfo.setNumUniquePeptides((Math.max(numUniquePeptidesThisProteinThisTime, proteinInfo.getNumUniquePeptides())));
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
        return proteinInfoMap;
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        float goodProbability = 0.1f;
        try
        {
            if (showCharts)
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
            }




            Map<String, ProteinInfo> proteinInfoMap1 =
                    loadProteinInfoMap(protXmlFiles[0]);

            List<Float> percentCoveragesGoodProb1 = new ArrayList<Float>();
            List<Float> spectralCountsGoodProb1 = new ArrayList<Float>();
            Map<String, Integer> spectralCountMap1 = new HashMap<String, Integer>();

            for (String protein : proteinInfoMap1.keySet())
            {
                ProteinInfo proteinInfo = proteinInfoMap1.get(protein);
                if (proteinInfo.getScore() > goodProbability && proteinInfo.getCoverage() != 0 && proteinInfo.getSpectralCount() != 0)
                {
                    percentCoveragesGoodProb1.add(proteinInfo.getCoverage());
                    spectralCountsGoodProb1.add((float)proteinInfo.getSpectralCount());
                    spectralCountMap1.put(protein, proteinInfo.getSpectralCount());
                }
            }
            ApplicationContext.infoMessage("Percent coverages (good proteins) 1: 5th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb1, 5))+ "%, 10th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb1, 10)) + "%, median: " +
                    (100 * BasicStatistics.median(percentCoveragesGoodProb1)) + "%");
            ApplicationContext.infoMessage("Min % coverage 1: " + (100 * BasicStatistics.min(percentCoveragesGoodProb1)) +
                    "%, max: " + (100 * BasicStatistics.max(percentCoveragesGoodProb1)));

            Map<String, ProteinInfo> proteinInfoMap2 =
                    loadProteinInfoMap(protXmlFiles[1]);
            List<Float> percentCoveragesGoodProb2 = new ArrayList<Float>();
            List<Float> spectralCountsGoodProb2 = new ArrayList<Float>();
            Map<String, Integer> spectralCountMap2 = new HashMap<String, Integer>();

            for (String protein : proteinInfoMap2.keySet())
            {
                ProteinInfo proteinInfo = proteinInfoMap2.get(protein);
                if (proteinInfo.getScore() > goodProbability && proteinInfo.getCoverage() != 0  && proteinInfo.getSpectralCount() != 0)
                {
                    percentCoveragesGoodProb2.add(proteinInfo.getCoverage());
                    spectralCountsGoodProb2.add((float)proteinInfo.getSpectralCount());
                    spectralCountMap2.put(protein, proteinInfo.getSpectralCount());

                }
            }
            ApplicationContext.infoMessage("Percent coverages (good proteins) 2: 5th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb2, 5))+ "%, 10th percentile: " +
                    (100 * BasicStatistics.percentile(percentCoveragesGoodProb2, 10)) + "%, median: " +
                    (100 * BasicStatistics.median(percentCoveragesGoodProb2)) + "%");
            ApplicationContext.infoMessage("Min % coverage 2: " + (100 * BasicStatistics.min(percentCoveragesGoodProb2)) +
                    "%, max: " + (100 * BasicStatistics.max(percentCoveragesGoodProb2)));




            Set<String> commonProteins = new HashSet<String>();

            for (String protein : proteinInfoMap1.keySet())
            {
                if (proteinInfoMap2.containsKey(protein))
                {
                    commonProteins.add(protein);
                }
            }

            if (outSpecCountFile != null)
            {
                PrintWriter pw = null;
                try
                {
                    pw = new PrintWriter(outSpecCountFile);
                    pw.println("protein\tspec1\tspec2\tratio");
                    for (String protein : spectralCountMap1.keySet())
                    {
                        if (spectralCountMap2.containsKey(protein))
                        {
                         pw.println(protein + "\t" + spectralCountMap1.get(protein) + "\t" + spectralCountMap2.get(protein) + "\t" +
                             (((float) spectralCountMap1.get(protein)) / ((float)spectralCountMap2.get(protein))) );
                        pw.flush();
                        }
                    }
                    pw.close();
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failed to write spectral count file",e);
                }
            }


            List<Float> scoresInCommon1 = new ArrayList<Float>();
            List<Float> scoresInCommon2 = new ArrayList<Float>();

            List<Float> percentsInCommon1 = new ArrayList<Float>();
            List<Float> percentsInCommon2 = new ArrayList<Float>();

            List<Float> spectralCountsInCommon1 = new ArrayList<Float>();
            List<Float> spectralCountsInCommon2 = new ArrayList<Float>();
            List<Float> spectralCountFoldChanges2 = new ArrayList<Float>();


            for (String protein : commonProteins)
            {
                scoresInCommon1.add(proteinInfoMap1.get(protein).getScore());
                scoresInCommon2.add(proteinInfoMap2.get(protein).getScore());

                if (proteinInfoMap1.get(protein).getCoverage() != 0 &&
                    proteinInfoMap2.get(protein).getCoverage() != 0)
//          && proteinScorePercentMap1.get(protein).first >= goodProbability && proteinScorePercentMap2.get(protein).first >= goodProbability)
                {
                    percentsInCommon1.add(proteinInfoMap1.get(protein).getCoverage());
                    percentsInCommon2.add(proteinInfoMap2.get(protein).getCoverage());
//if (proteinScorePercentMap2.get(protein).second < proteinScorePercentMap1.get(protein).second)
//    System.err.println("Decreasing % coverage: " + protein);
                }

                if (proteinInfoMap1.get(protein).getSpectralCount() != 0 &&
                    proteinInfoMap2.get(protein).getSpectralCount() != 0)
                {
                    spectralCountsInCommon1.add((float)proteinInfoMap1.get(protein).getSpectralCount());
                    spectralCountsInCommon2.add((float)proteinInfoMap2.get(protein).getSpectralCount());
                    spectralCountFoldChanges2.add(((float)proteinInfoMap2.get(protein).getSpectralCount() - (float)proteinInfoMap1.get(protein).getSpectralCount())/ (float)proteinInfoMap1.get(protein).getSpectralCount());
                }
            }

            ApplicationContext.infoMessage("Proteins in file 1: " + proteinInfoMap1.size() +
                    ", file 2: " + proteinInfoMap2.size() + ", common: " + commonProteins.size());
//for (String protIn1 : proteinScoreMap1.keySet()) if (!commonProteins.contains(protIn1)) System.err.println("**" + protIn1);

            ApplicationContext.infoMessage("Covariance by spectral count: " + BasicStatistics.covariance(scoresInCommon1, scoresInCommon2));

            if (showCharts)
            {
                PanelWithScatterPlot psp = new PanelWithScatterPlot(scoresInCommon1, scoresInCommon2, "protein scores");
                psp.setAxisLabels("Protein Probability, File 1","Protein Probability, File 2");
                psp.setPointSize(2);
                psp.displayInTab();


//System.err.println("Spectral counts in common: " + spectralCountsInCommon1.size());

                PanelWithScatterPlot pspSpecCount = new PanelWithScatterPlot(spectralCountsInCommon1, spectralCountsInCommon2, "spectral counts");
                pspSpecCount.setAxisLabels("Spectral Count, File 1","Spectral Count file 2");
                pspSpecCount.setPointSize(2);
                pspSpecCount.displayInTab();

                List<Float> spectralCountLogRatios = new ArrayList<Float>(spectralCountsInCommon1.size());
                List<Float> spectralCountSums = new ArrayList<Float>(spectralCountsInCommon1.size());
                for (int i=0; i<spectralCountsInCommon1.size(); i++)
                {
                    spectralCountLogRatios.add((float)Math.log((spectralCountsInCommon1.get(i) + (Math.random()/5))/ (spectralCountsInCommon2.get(i) + (Math.random()/5))));
                    spectralCountSums.add((float) Math.log(spectralCountsInCommon1.get(i) + spectralCountsInCommon2.get(i) + (Math.random()/2)));
                }
                PanelWithScatterPlot pspSpecCountRatioMA = new PanelWithScatterPlot(spectralCountSums, spectralCountLogRatios, "spectral count MA");
                pspSpecCountRatioMA.setAxisLabels("Spectral Count Sum","Log Spectral Count Ratio");
                pspSpecCountRatioMA.setPointSize(3);
                pspSpecCountRatioMA.displayInTab();
            }

            List<Float> spectralCountPercentilesInCommon1 = new ArrayList<Float>();
            List<Float> spectralCountPercentilesInCommon2 = new ArrayList<Float>();
            List<Float> percentileCutoffs1 = new ArrayList<Float>();
            List<Float> percentileCutoffs2 = new ArrayList<Float>();
            for (int i=1; i<101; i++)
            {
                percentileCutoffs1.add((float)BasicStatistics.percentile(spectralCountsGoodProb1, i));
                percentileCutoffs2.add((float)BasicStatistics.percentile(spectralCountsGoodProb2, i));
            }



            for (int i=0; i<spectralCountsInCommon1.size(); i++)
            {
                spectralCountPercentilesInCommon1.add(
                        (float) Math.abs(Collections.binarySearch(percentileCutoffs1, spectralCountsInCommon1.get(i))) + 1);
                spectralCountPercentilesInCommon2.add(
                        (float) Math.abs(Collections.binarySearch(percentileCutoffs2, spectralCountsInCommon2.get(i))) + 1);
//System.err.println("***" + spectralCountsInCommon1.get(i) + "-> " + spectralCountPercentilesInCommon1.get(i));
            }
            if (showCharts)
            {
                PanelWithScatterPlot pspSpecCountPercentile =
                        new PanelWithScatterPlot(spectralCountPercentilesInCommon1, spectralCountPercentilesInCommon2, "spectral count percentiles");
                pspSpecCountPercentile.setAxisLabels("Spectral Count %ile, File 1","Spectral Count %ile file 2");
                pspSpecCountPercentile.setPointSize(2);
                pspSpecCountPercentile.displayInTab();

//for (String protIn1 : proteinScoreMap1.keySet()) if (!commonProteins.contains(protIn1)) System.err.println("**" + protIn1);
                PanelWithScatterPlot pspSpecCountFoldChange = new PanelWithScatterPlot(spectralCountsInCommon1, spectralCountFoldChanges2, "spectral count fold");
                pspSpecCountFoldChange.setAxisLabels("Spectral Count, File 1","Spectral Count fold change, file 2");
                pspSpecCountFoldChange.setPointSize(2);
                pspSpecCountFoldChange.displayInTab();

                new PanelWithHistogram(spectralCountPercentilesInCommon1, "Speccount Percentile 1").displayInTab();
                new PanelWithHistogram(spectralCountPercentilesInCommon2, "Speccount Percentile 2").displayInTab();
            }

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
            if (showCharts)
            {
                PanelWithScatterPlot pspCoverage2 = new PanelWithScatterPlot(percentsInCommon1, foldChangesInCommon2, "prot coverage fold change");
                pspCoverage2.setAxisLabels("Protein % Coverage, File 1","Coverage fold change, File 2");
                pspCoverage2.setPointSize(2);
                pspCoverage2.displayInTab();
            }

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

            for (String protein : proteinInfoMap2.keySet())
            {
                if (!commonProteins.contains(protein))
                {
                    proteinsIn2ButNot1.add(protein);
                    scoresIn2ButNot1.add(proteinInfoMap2.get(protein).getScore());
                    if (proteinInfoMap2.get(protein).getCoverage() != 0)
                        coveragesIn2ButNot1.add(proteinInfoMap2.get(protein).getCoverage());
                }
            }
            if (scoresIn2ButNot1.size() > 0 && showCharts)
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


            Set<String> proteinsAboveGoodProb1 = new HashSet<String>();
            Set<String> proteinsAboveGoodProb2 = new HashSet<String>();
            Set<String> proteinsProb01 = new HashSet<String>();
            Set<String> proteinsProb02 = new HashSet<String>();
            Set<String> proteinsBelowPoint11 = new HashSet<String>();
            Set<String> proteinsBelowPoint12 = new HashSet<String>();

            Set<String> proteinsFallingBelowGoodProbIn2 = new HashSet<String>();


            for (String protein  : proteinInfoMap1.keySet())
            {
                float score = proteinInfoMap1.get(protein).getScore();
                if (score >= .75f)
                    scoresAbovePoint751++;
                if (score >= goodProbability)
                    proteinsAboveGoodProb1.add(protein);
                if (score == 0f)
                    proteinsProb01.add(protein);
                if (score < 0.2f)
                    proteinsBelowPoint11.add(protein);
                scores1.add(score);
            }

            for (String protein  : proteinInfoMap2.keySet())
            {
                float score = proteinInfoMap2.get(protein).getScore();
                if (score >= .75f)
                    scoresAbovePoint752++;

                if (score >= goodProbability)
                    proteinsAboveGoodProb2.add(protein);
                else if (proteinsAboveGoodProb1.contains(protein))
                    proteinsFallingBelowGoodProbIn2.add(protein);

                if (score == 0f)
                    proteinsProb02.add(protein);
                if (score < 0.2f)
                    proteinsBelowPoint12.add(protein);
                scores2.add(score);
            }


            Set<String> commonProteinsAbovePoint95 = new HashSet<String>();
            for (String protein : proteinsAboveGoodProb1)
            {
                if (proteinsAboveGoodProb2.contains(protein))
                    commonProteinsAbovePoint95.add(protein);
            }
            Set<String> proteinsAboveGoodProbUniqueTo2 = new HashSet<String>();
            for (String protein : proteinsAboveGoodProb2)
            {
                if (!commonProteinsAbovePoint95.contains(protein))
                    proteinsAboveGoodProbUniqueTo2.add(protein);
            }
            ApplicationContext.infoMessage("Scores above " + goodProbability + ": Set 1 = " + proteinsAboveGoodProb1.size() +
                                           ", Set 2 = " + proteinsAboveGoodProb2.size() + ", common=" +
                                           commonProteinsAbovePoint95.size() + ", unique to 2: " +
                                           ( proteinsAboveGoodProb2.size() - commonProteinsAbovePoint95.size()));
            int numCommonProteinsProb0 = 0;
            for (String protein01 : proteinsProb01)
                if (proteinsProb02.contains(protein01))
                    numCommonProteinsProb0++;
            ApplicationContext.infoMessage("Proteins with score 0: Set 1: " + proteinsProb01.size() + ", Set 2: " +
                    proteinsProb02.size() + ", common: " + numCommonProteinsProb0);

            int numCommonProteinsProbBelowPoint11 = 0;
            for (String proteinBelowPoint11 : proteinsBelowPoint11)
                if (proteinsBelowPoint12.contains(proteinBelowPoint11))
                    numCommonProteinsProbBelowPoint11++;
            ApplicationContext.infoMessage("Proteins with score below 0.1: Set 1: " + proteinsBelowPoint11.size() +
                    ", Set 2: " + proteinsBelowPoint12.size() + ", common: " + numCommonProteinsProbBelowPoint11);

            if (_log.isDebugEnabled() && !proteinsAboveGoodProbUniqueTo2.isEmpty())
            {
                System.err.println("High-confidence proteins unique to 2:");
                for (String prot : proteinsAboveGoodProbUniqueTo2)
                    System.err.println("\t" + prot);
            }
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

            if (showCharts)
            {
                PanelWithHistogram pwh2 = new PanelWithHistogram(scores1, "Scores in 1");
                pwh2.displayInTab();
                PanelWithHistogram pwh3 = new PanelWithHistogram(scores2, "Scores in 2");
                pwh3.displayInTab();
            }


            List<ProtXmlReader.Protein> proteinList1 =
                    ProteinUtilities.loadProtXmlProteinsFromProtXML(protXmlFiles[0]);
            List<ProtXmlReader.Protein> proteinList2 =
                    ProteinUtilities.loadProtXmlProteinsFromProtXML(protXmlFiles[1]);
            List<Float> allLogRatios1 = new ArrayList<Float>();
            Set<String> point95RatioProteins1 = new HashSet<String>();

/*
            List<Float> coveragesWithProbFallingBelowCutoffIn2 = new ArrayList<Float>();
            for (String proteinFallingBelowGoodProbIn2 : proteinsFallingBelowGoodProbIn2)
                coveragesWithProbFallingBelowCutoffIn2.add((float) proteinInfoMap2.get(proteinFallingBelowGoodProbIn2).getNumUniquePeptides() - (float) proteinInfoMap1.get(proteinFallingBelowGoodProbIn2).getNumUniquePeptides());
            PanelWithHistogram fallingProbCoverageHistogram =
                    new PanelWithHistogram(coveragesWithProbFallingBelowCutoffIn2, "File 1 numpeps, prob falling file 2");
            fallingProbCoverageHistogram.displayInTab();

            List<Float> coveragesWithHighProbIn1 = new ArrayList<Float>();
            for (String proteinAboveGoodProb1 : proteinsAboveGoodProb1)
            if (proteinInfoMap2.containsKey(proteinAboveGoodProb1))
                    coveragesWithHighProbIn1.add((float) proteinInfoMap2.get(proteinAboveGoodProb1).getNumUniquePeptides() - (float) proteinInfoMap1.get(proteinAboveGoodProb1).getNumUniquePeptides());
            PanelWithHistogram highProb1CoverageHistogram =
                    new PanelWithHistogram(coveragesWithHighProbIn1, "File 1 numpeps, good prob");
            highProb1CoverageHistogram.displayInTab();
*/            


            List<ProteinGroup> proteinGroups1 = ProteinUtilities.loadProteinGroupsFromProtXML(protXmlFiles[0]);
            List<ProteinGroup> proteinGroups2 = ProteinUtilities.loadProteinGroupsFromProtXML(protXmlFiles[1]);

            //todo: actually calculate this
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
                List<String> namesWithRatios = new ArrayList<String>();


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
                        namesWithRatios.add(proteinName);

                        scoresWithRatios2.add(protein2.getProbability());
//if (Math.abs(ratio2 - ratio1) > 0.0001)
//    System.err.println(proteinName + ", " + ratio1 + ", " + ratio2 + ", " + (ratio2-ratio1));
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

                if (showCharts)
                    new PanelWithScatterPlot(scoresWithRatios1, scoresWithRatios2, "Probs w/ratios").displayInTab();

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



                    if (showCharts)
                        new PanelWithScatterPlot(ratios1, ratios2, "Quantitation ratios").displayInTab();
                    System.err.println("Ratios to compare: " + ratios1.size());

                    double[] logRatios1 = new double[ratios1.size()];
                    double[] logRatios2 = new double[ratios2.size()];

                    for (int i=0; i<ratios1.size(); i++)
                    {
                        logRatios1[i] = Math.log(Math.max(ratios1.get(i),.001));
                        logRatios2[i] = Math.log(Math.max(ratios2.get(i),.001));
                    }

                    //Write a file with names and ratios for proteins seen with ratios in both files
                    if (outProteinRatioFile != null)
                    {
                        try
                        {
                            PrintWriter pw = new PrintWriter(outProteinRatioFile);
                            pw.println("protein\tratio1\tratio2");
                            for (int i=0; i<ratios1.size(); i++)
                            {
                                pw.println(namesWithRatios.get(i) + "\t" + ratios1.get(i) + "\t" + ratios2.get(i));
                                pw.flush();
                            }
                            pw.close();
                            ApplicationContext.infoMessage("Wrote protein ratio file " + outProteinRatioFile.getAbsolutePath());
                        }
                        catch (IOException e)
                        {
                            throw new CommandLineModuleExecutionException("Failed to write protein ratio file",e);
                        }

                    }

                    if (showCharts)
                        new PanelWithScatterPlot(logRatios1, logRatios2, "log ratios").displayInTab();

                    ApplicationContext.infoMessage("Log ratio correlation coefficient: " +
                            BasicStatistics.correlationCoefficient(logRatios1, logRatios2));
                }
                else
                    ApplicationContext.setMessage("No ratios to compare on plot");

                if (allLogRatios1.size() >=1 && showCharts)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(allLogRatios1, "Log Ratios 1");
                    pwh.displayInTab();
                }
               if (allLogRatios2.size() >=1 && showCharts)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(allLogRatios2, "Log Ratios 2");
                    pwh.displayInTab();
                }
            }





            int numGroupsAbovePoint952 = 0;
            int numQuantGroupsAbovePoint952 = 0;
            int commonQuantGroupsAbovePoint95 = 0;


            List<Float> numPeptidesPerGroup2 = new ArrayList<Float>();
            List<Integer> uniqueProteinGroups1 = new ArrayList<Integer>();
            for (ProteinGroup proteinGroup2 : proteinGroups2)
            {
                boolean abovePoint95 = false;
                if (proteinGroup2.getProbability() >= goodProbability)
                    abovePoint95 = true;
                if (abovePoint95)
                    numGroupsAbovePoint952++;
                Collection<String> proteinNames1 = proteinInfoMap1.keySet();
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

                    if (abovePoint95 && proteinsAboveGoodProb1.contains(protein.getProteinName()))
                        inCommonAbovePoint95 = true;
                    if (protein.getQuantitationRatio() != null)
                        quantified=true;
                    if (abovePoint95 && quantified && point95RatioProteins1.contains(protein.getProteinName()))
                        inCommonQuantAbovePoint95 = true;
                }
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
                    " (" + Rounder.round((100f * (float) commonGroupsAbovePoint95 / (float) (numGroupsAbovePoint951 +
                    numGroupsAbovePoint952 - commonGroupsAbovePoint95)),2) + "%)" +                     
                    ", Unique to 2: " + (numGroupsAbovePoint952 - commonGroupsAbovePoint95));
            ApplicationContext.infoMessage("Quantified ProteinGroups above " + goodProbability + ": Set 1 = " + numQuantGroupsAbovePoint951 +
                    ", Set 2 = " + numQuantGroupsAbovePoint952 + ", common (by IPI) = " + commonQuantGroupsAbovePoint95 +
                    ", Unique to 2: " + (numQuantGroupsAbovePoint952 - commonQuantGroupsAbovePoint95));

            if (showCharts)
            {
                PanelWithHistogram peptidesPerGroupHist1 = new PanelWithHistogram(numPeptidesPerGroup1, "#Peptides1");
                peptidesPerGroupHist1.displayInTab();

                PanelWithHistogram peptidesPerGroupHist2 = new PanelWithHistogram(numPeptidesPerGroup2, "#Peptides2");
                peptidesPerGroupHist2.displayInTab();
            }


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
