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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command line moeule
 */
public class PepArrayMultiSampleAnalyzerCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PepArrayMultiSampleAnalyzerCLM.class);

    protected File file;
    protected File outFile;
    protected File outDir;

    protected int maxMissedCleavages = 0;


    protected int minRunsPerGroup = 2;

    protected boolean showCharts = false;

    Object[] arrayRows;

    //can be populated from a FASTA database or a protXML search results file
    protected Map<String, Collection<String>> peptideProteinMap = null;

    protected PeptideArrayAnalyzer arrayAnalyzer = null;



    protected boolean collapseRowsByPeptide = true;

    public PepArrayMultiSampleAnalyzerCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "peparraymultisample";

        mShortDescription = "Compares two sets of samples.  Performs t-test for each peptide";
        mHelpMessage = mShortDescription;

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, null),
                        new FileToWriteArgumentDefinition("out",false, "output file"),
                        new FileToWriteArgumentDefinition("outdir",false, "output directory"),

                        new FileToReadArgumentDefinition("caserunlistfile",true,"File containing the names of runs in the case group, one per line"),
                        new FileToReadArgumentDefinition("controlrunlistfile",true,"File containing the names of runs in the control group, one per line"),


                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new IntegerArgumentDefinition("minrunspergroup", false,
                                "Minimum number of runs in each group in which a feature must be located to be counted",2),
                        new FileToReadArgumentDefinition("fasta", false, "FASTA database (must be provided if protxml not provided)"),
                        new FileToReadArgumentDefinition("protxml", false, "protXML search results file (must be provided if fasta not provided)"),
                        new BooleanArgumentDefinition("collapsebypeptide", false, "Collapse rows by peptide, using median value?", collapseRowsByPeptide),
                        new IntegerArgumentDefinition("missedcleavages", false,"maximum missed cleavages for peptides (if fasta provided)", maxMissedCleavages),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        file = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        minRunsPerGroup = getIntegerArgumentValue("minrunspergroup");

        showCharts = getBooleanArgumentValue("showcharts");

        File fastaFile = getFileArgumentValue("fasta");
        maxMissedCleavages = getIntegerArgumentValue("missedcleavages");
        File protXmlFile = getFileArgumentValue("protxml");
        if (fastaFile != null && protXmlFile != null)
        {
            assertArgumentAbsent("fasta","protxml");
        }
        if (fastaFile == null && protXmlFile == null)
        {
            assertArgumentPresent("protxml","fasta");
        }
        try
        {
            if (fastaFile != null)
            {
                //really stupid to have to do this, but for some reason map<string,set<string>>
                // doesn't cast to map<string,collection<string>> .  Same with list->collection
                peptideProteinMap = new HashMap<String,Collection<String>>(ProteinUtilities.loadTrypticPeptideProteinMapFromFasta(fastaFile, maxMissedCleavages));
                ApplicationContext.infoMessage("Loaded peptide-protein map from FASTA.  Warning: protein summary " +
                        "may include peptides that ProteinProphet would not assign");
            }
            else
            {
                peptideProteinMap = new HashMap<String,Collection<String>>(ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFile, 0.1));
                ApplicationContext.infoMessage("Loaded peptide-protein map from protXML file");
            }
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("failed to load peptide-protein mapping",e);
        }


        collapseRowsByPeptide = getBooleanArgumentValue("collapsebypeptide");

        try
        {
            arrayAnalyzer = new PeptideArrayAnalyzer(getUnnamedFileArgumentValue());
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("failed to load peptide array",e);
        }


        try
        {
            arrayAnalyzer.loadCaseControlRunListsFromFiles(getFileArgumentValue("caserunlistfile"),
                    getFileArgumentValue("controlrunlistfile"));
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load case or control run list file", e);
        }

        ApplicationContext.infoMessage("Case runs:");
        for (String caseRun : arrayAnalyzer.getCaseRunNames())
            ApplicationContext.setMessage("\t" + caseRun);
        ApplicationContext.infoMessage("Control runs:");
        for (String controlRun : arrayAnalyzer.getControlRunNames())
            ApplicationContext.setMessage("\t" + controlRun);
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            float[][] resultMatrix = arrayAnalyzer.runTwoGroupTTestAllRows(2);

            Map<String, List<Float>> peptideTScoresMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideRatiosMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideCaseMeanIntensitiesMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideControlMeanIntensitiesMap = new HashMap<String, List<Float>>();


            List<Float> noNaP = new ArrayList<Float>();
            List<Float> noNaQ = new ArrayList<Float>();
            List<Float> noNaT = new ArrayList<Float>();

            for (int i=0; i<resultMatrix.length; i++)
            {
                if (!Float.isNaN(resultMatrix[i][0]))
                {
                    noNaP.add(resultMatrix[i][0]);
                    noNaT.add(resultMatrix[i][1]);
                    noNaQ.add(resultMatrix[i][2]);
                }
            }

            if (showCharts)
            {
                new PanelWithHistogram(noNaP, "pep p-values").displayInTab();
                new PanelWithHistogram(noNaQ, "pep q-values").displayInTab();
                new PanelWithHistogram(noNaT, "pep t-scores").displayInTab();
            }

            Map<String, List<Integer>> peptideRowIndexesMap = arrayAnalyzer.loadPeptideRowIndexesMap();
            for (String peptide : peptideRowIndexesMap.keySet())
            {
                if (!peptideProteinMap.containsKey(peptide))
                {
                    _log.debug("Skipping unknown peptide " + peptide);
                    continue;
                }

                List<Float> tScoresThisPeptide = new ArrayList<Float>();
                peptideTScoresMap.put(peptide, tScoresThisPeptide);

                List<Float> ratiosThisPeptide = new ArrayList<Float>();
                peptideRatiosMap.put(peptide, ratiosThisPeptide);

                List<Float> caseMeanIntensitiesThisPeptide = new ArrayList<Float>();
                peptideCaseMeanIntensitiesMap.put(peptide, caseMeanIntensitiesThisPeptide);

                List<Float> controlMeanIntensitiesThisPeptide = new ArrayList<Float>();
                peptideControlMeanIntensitiesMap.put(peptide, controlMeanIntensitiesThisPeptide);

                for (int i : peptideRowIndexesMap.get(peptide))
                {
                    if (!Float.isNaN(resultMatrix[i][1]))
                    {
                        tScoresThisPeptide.add(resultMatrix[i][1]);
                        ratiosThisPeptide.add(
                                arrayAnalyzer.calcRowCaseControlRatioOfGeomMeans(
                                        arrayAnalyzer.getRowMap(i)).floatValue());
                        caseMeanIntensitiesThisPeptide.add((float)
                                BasicStatistics.geometricMean(arrayAnalyzer.getPresentCaseIntensities(
                                        arrayAnalyzer.getRowMap(i))));
                        controlMeanIntensitiesThisPeptide.add((float)
                                BasicStatistics.geometricMean(arrayAnalyzer.getPresentControlIntensities(
                                        arrayAnalyzer.getRowMap(i))));
                    }
                }
            }

            Set<String> observedProteinsSet = new HashSet<String>();
            for (String peptide : peptideTScoresMap.keySet())
            {
                if (!peptideTScoresMap.get(peptide).isEmpty())
                {
                    if (collapseRowsByPeptide)
                    {
                    //This is a bit hacky.  Collapsing to a single peptide.  Rather than explicitly
                    //change the data model, I'm just holding a single value per peptide, the median of all the
                    //non-NaN values we had
                        List<Float> newList = new ArrayList<Float>();
                        newList.add((float) BasicStatistics.median(peptideTScoresMap.get(peptide)));
                        peptideTScoresMap.put(peptide, newList);

                        List<Float> newListR = new ArrayList<Float>();
                        newListR.add((float) BasicStatistics.geometricMean(peptideRatiosMap.get(peptide)));
                        peptideRatiosMap.put(peptide, newListR);

                        List<Float> newListCaseI = new ArrayList<Float>();
                        newListCaseI.add((float) BasicStatistics.geometricMean(peptideCaseMeanIntensitiesMap.get(peptide)));
                        peptideCaseMeanIntensitiesMap.put(peptide, newListCaseI);

                        List<Float> newListControlI = new ArrayList<Float>();
                        newListControlI.add((float) BasicStatistics.geometricMean(peptideControlMeanIntensitiesMap.get(peptide)));
                        peptideControlMeanIntensitiesMap.put(peptide, newListControlI);
                    }
                    observedProteinsSet.addAll(peptideProteinMap.get(peptide));
                }
            }
            List<String> observedProteins = new ArrayList<String>(observedProteinsSet);

            writeGSEAFile( observedProteins,  peptideTScoresMap, peptideRatiosMap,
                    peptideCaseMeanIntensitiesMap, peptideControlMeanIntensitiesMap); 

        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Error running t test", e);
        }
    }

    protected void writeGSEAFile(List<String> observedProteins,
                                 Map<String, List<Float>> peptideTScoresMap,
                                 Map<String, List<Float>> peptideRatiosMap,
                                 Map<String, List<Float>> peptideCaseMeanIntensitiesMap,
                                 Map<String, List<Float>> peptideControlMeanIntensitiesMap)
            throws CommandLineModuleExecutionException
    {
            PrintWriter outPW = null;
            try
            {
                outPW = new PrintWriter(outFile);
            }
            catch (FileNotFoundException e)
            {
                throw new CommandLineModuleExecutionException("Failed to open output file " + outFile.getAbsolutePath());
            }
            StringBuffer headerLineBuf = new StringBuffer("peptide\tratio\ttscore\tintmeancase\tintmeancontrol");
            for (String protein : observedProteins)
            {
                headerLineBuf.append("\t" + protein);
            }
            outPW.println(headerLineBuf);
            outPW.flush();

            for (String peptide : peptideTScoresMap.keySet())
            {
                StringBuffer proteinsPartBuf = new StringBuffer();
                Collection<String> proteinsThisPeptide = peptideProteinMap.get(peptide);
                for (String protein : observedProteins)
                {
                    if (proteinsThisPeptide.contains(protein))
                        proteinsPartBuf.append("\t1");
                    else
                        proteinsPartBuf.append("\t0");
                }
                for (int i=0; i<peptideTScoresMap.get(peptide).size(); i++)
                {
                    float tscore = peptideTScoresMap.get(peptide).get(i);
                    float ratio = peptideRatiosMap.get(peptide).get(i);
                    float caseMean = peptideCaseMeanIntensitiesMap.get(peptide).get(i);
                    float controlMean = peptideControlMeanIntensitiesMap.get(peptide).get(i);

                    if (!Float.isNaN(tscore))
                        outPW.println(peptide + "\t" + ratio + "\t" + tscore + "\t" + caseMean + "\t" + controlMean + proteinsPartBuf.toString());
                }
            }
            outPW.close();
    }

}
