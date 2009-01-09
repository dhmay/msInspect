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
package org.fhcrc.cpl.viewer.align;

import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule for feature finding
 */
public class PeptideArrayAnalyzer
{
    protected static Logger _log =
            Logger.getLogger(PeptideArrayAnalyzer.class);

    protected Map<String,Object>[] rowMaps = null;
    protected List<String> runNames = null;

    protected String[] caseRunNames;
    protected String[] controlRunNames;

    protected boolean showCharts = false;

    public static final int CONSENSUS_INTENSITY_MODE_MEAN = 0;
    public static final int CONSENSUS_INTENSITY_MODE_FIRST = 1;


    public PeptideArrayAnalyzer(File arrayFile)
            throws IOException
    {
        parseArray(arrayFile);
    }

    public void parseArray(File arrayFile) throws IOException
    {
        TabLoader tabLoader = new TabLoader(arrayFile);
        tabLoader.setReturnElementClass(HashMap.class);
        Object[] rows = tabLoader.load();

        runNames = new ArrayList<String>();

        for (TabLoader.ColumnDescriptor column : tabLoader.getColumns())
        {
            _log.debug("loading column " + column.name);
            if (column.name.startsWith("intensity_"))
            {
                runNames.add(column.name.substring("intensity_".length()));
                _log.debug("adding run " + runNames.get(runNames.size() - 1));
            }
        }

        rowMaps = (Map<String, Object>[]) rows;
    }



    public int countConflictRows()
    {
        int conflictRows = 0;

        for (Map<String,Object> rowMap : rowMaps)
        {
            String rowPeptide = null;
            for (String runName : runNames)
            {
                String runPeptide = (String) rowMap.get("peptide_" + runName);
                if (runPeptide != null)
                {
                    if (rowPeptide == null)
                    {
                        rowPeptide = runPeptide;
                    }
                    else
                    {
                        if (!runPeptide.equals(rowPeptide))
                        {
                            conflictRows++;
                            break;
                        }
                    }
                }
            }
        }
        return conflictRows;
    }

    /**
     * Return the sets of peptides matched in each run
     * @return
     */
    public Map<String,Set<String>> getRunMatchedPeptides()
    {
        Map<String,Set<String>> runOriginalPeptides =
                new HashMap<String,Set<String>>();

        for (String runName : runNames)
        {
            runOriginalPeptides.put(runName, new HashSet<String>());
        }
        for (Map<String,Object> rowMap : rowMaps)
        {
            for (String runName : runNames)
            {
                String runPeptide = (String) rowMap.get("peptide_" + runName);
                if (runPeptide != null)
                {
                    runOriginalPeptides.get(runName).add(runPeptide);
                }
            }
        }
        return runOriginalPeptides;
    }

    /**
     * identify the peptides matched in each of these runs, then fill in matches
     * across a single row: if some runs in a row match the same peptide ID, but
     * others have no ID, fill in that ID for the ones with no match
     * @param includeOriginalPeptides
     * @return
     */
    public Map<String,Set<String>> fillInMatchedPeptides(boolean includeOriginalPeptides)
    {
        Map<String,Set<String>> runNewPeptides = new HashMap<String,Set<String>>();
        Map<String,Set<String>> runOriginalPeptides = getRunMatchedPeptides();

        for (String runName : runNames)
        {
            runNewPeptides.put(runName, new HashSet<String>());
        }

        int conflictRows = 0;
        for (Map<String,Object> rowMap : rowMaps)
        {
            boolean agreement = true;
            String rowPeptide = null;
            for (String runName : runNames)
            {
                String runPeptide = (String) rowMap.get("peptide_" + runName);
                if (runPeptide != null)
                {
                    if (rowPeptide == null)
                    {
                        rowPeptide = runPeptide;
                    }
                    else
                    {
                        if (!runPeptide.equals(rowPeptide))
                        {
                            conflictRows++;
                            agreement = false;
                            break;
                        }
                    }
                }
            }
            if (rowPeptide != null && agreement)
            {
                for (String runName : runNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + runName);
                    if (runIntensity == null)
                        continue;

 //                   double intensity = Double.parseDouble(runIntensity.toString());
 //                       includeOriginalPeptides

                    if (
                        !runOriginalPeptides.get(runName).contains(rowPeptide))
                    {
                        runNewPeptides.get(runName).add(rowPeptide);
                    }
                }
            }


        }
//        for (String runName : runNames)
//            System.err.println("\t" + runName + ": " + runNewPeptides.get(runName).size());        
//            System.err.println("Conflict rows: " + conflictRows);
//        System.err.println("NEW peptide counts per run:");

        if (!includeOriginalPeptides)
            return runNewPeptides;

        //just change name
        Map<String,Set<String>> allPeptides = runNewPeptides;
        for (String runName : runNames)
        {
            Set<String> runPeptides = allPeptides.get(runName);
            for (String originalPeptide : runOriginalPeptides.get(runName))
                runPeptides.add(originalPeptide);
        }

        return allPeptides;
    }

    /**
     * Create feature files for features that occur in all runs, with no conflicts
     */
    public FeatureSet createConsensusFeatureSet(File detailsFile,
                                           int minConsensusRuns,
                                           int intensityMode)
            throws IOException
    {
        TabLoader detailsTabLoader = new TabLoader(detailsFile);
        detailsTabLoader.setReturnElementClass(HashMap.class);
        Object[] detailRows = detailsTabLoader.load();
        Map<String, Object>[] detailsRowMaps = (Map<String, Object>[]) detailRows;

        Map<String, List<Feature>> runFeatureLists =
                new HashMap<String, List<Feature>>(runNames.size());
        for (String runName : runNames)
            runFeatureLists.put(runName, new ArrayList<Feature>());

        int currentDetailsRowId = -1;
        int currentDetailsArrayIndex = -1;
        Map<String, Object> currentDetailsRow = null;

        List<Feature> resultFeatureList = new ArrayList<Feature>();
        for (int i=0; i<rowMaps.length; i++)
        {
            int rowId = (Integer) rowMaps[i].get("id");
            Map<String, Object> rowMap = rowMaps[i];
            //minimum consensus runs achieved
            if ((Integer) rowMap.get("setCount") >= minConsensusRuns)
            {
                Map<String,List<Feature>> thisRowRunFeatureMap =
                        new HashMap<String,List<Feature>>();
                while (currentDetailsRowId < rowId)
                {
                    currentDetailsArrayIndex++;
                    currentDetailsRow = detailsRowMaps[currentDetailsArrayIndex];
                    currentDetailsRowId = (Integer) currentDetailsRow.get("id");
                }
                while (currentDetailsRowId == rowId && currentDetailsArrayIndex < detailsRowMaps.length)
                {
                    String fileName = (String) currentDetailsRow.get("file");
                    String runName =
                         fileName.substring(0, fileName.indexOf("."));
                    Feature feature = new Feature(
                            (Integer) currentDetailsRow.get("scan"),
                            (Integer) currentDetailsRow.get("scanFirst"),
                            (Integer) currentDetailsRow.get("scanLast"),
                            ((Double) currentDetailsRow.get("mz")).floatValue(),
                            ((Double) currentDetailsRow.get("intensity")).floatValue(),
                            (Integer) currentDetailsRow.get("charge"),
                            ((Double) currentDetailsRow.get("kl")).floatValue(),
                            ((Double) currentDetailsRow.get("totalIntensity")).floatValue()
                            );
                    feature.setTime(((Double) currentDetailsRow.get("time")).floatValue());
if (feature.getMass() > 194 && feature.getMass() < 195) System.err.println(feature);
                    if (!thisRowRunFeatureMap.containsKey(runName))
                            thisRowRunFeatureMap.put(runName, new ArrayList<Feature>());
                    //todo: peptides, proteins
                    thisRowRunFeatureMap.get(runName).add(feature);

                    currentDetailsArrayIndex++;
                    if (currentDetailsArrayIndex < detailsRowMaps.length)
                    {
                        currentDetailsRow = detailsRowMaps[currentDetailsArrayIndex];
                        currentDetailsRowId = (Integer) currentDetailsRow.get("id");
                    }
                }

                List<Double> thisFeatureIntensities = new ArrayList<Double>();
                Feature firstFeatureOccurrence = null;
                for (String runName : thisRowRunFeatureMap.keySet())
                {
                    List<Feature> featuresThisRowThisRun = thisRowRunFeatureMap.get(runName);
                    //todo: NOT restricting to those runs that only have one feature.  OK?
                    if (featuresThisRowThisRun.size() >= 1)
                    {
                        Feature feature = featuresThisRowThisRun.get(0);
                        runFeatureLists.get(runName).add(feature);
                        if (firstFeatureOccurrence == null)
                            firstFeatureOccurrence = feature;
                        thisFeatureIntensities.add((double) feature.getIntensity());
                    }
                }
                Feature consensusFeature = firstFeatureOccurrence;
                if (consensusFeature == null) continue;

                switch (intensityMode)
                {
                    case CONSENSUS_INTENSITY_MODE_MEAN:
                        double[] featureIntensitiesArray = new double[thisFeatureIntensities.size()];
                        for (int j=0; j<thisFeatureIntensities.size(); j++)
                            featureIntensitiesArray[j] = thisFeatureIntensities.get(j);
                        consensusFeature.setIntensity((float) BasicStatistics.mean(featureIntensitiesArray));
                        break;
                    case CONSENSUS_INTENSITY_MODE_FIRST:
                        consensusFeature.setIntensity(firstFeatureOccurrence.getIntensity());
                        break;
                }
                resultFeatureList.add(consensusFeature);

            }
        }

        //this was for creating separate featuresets for each run
//        for (String runName : runFeatureLists.keySet())
//        {
//            FeatureSet runFeatureSet =
//                    new FeatureSet(runFeatureLists.get(runName).toArray(new Feature[runFeatureLists.size()]));
//            File outFile = new File(outDir, runName + ".peptides.tsv");
//            runFeatureSet.save(outFile);
//        }
        return new FeatureSet(resultFeatureList.toArray(new Feature[resultFeatureList.size()]));
    }

    /**
     * Create feature files based on all features that match nonconflicting peptides,
     * write them to outdir
     */
    public void createFeatureFilesAllMatched(File outDir)
    {
        Map<String,Set<String>> runNewPeptideMap = fillInMatchedPeptides(false);

        for (String runName : runNames)
        {
            List<Feature> runFeatureList = new ArrayList<Feature>();
            for (String peptide : runNewPeptideMap.get(runName))
            {
                        Feature feature = new Feature();
                        feature.setScan(1);
                        feature.setMass(1);
                        feature.setMz(1);
                        feature.setPeaks(1);
                        feature.setIntensity((float) 100);
                        MS2ExtraInfoDef.addPeptide(feature, peptide);
                runFeatureList.add(feature);
            }


            FeatureSet runFeatureSet = new FeatureSet(runFeatureList.toArray(new Feature[0]));
            File outFeatureFile = new File(outDir,runName + ".fromarray.tsv");
            try
            {
                runFeatureSet.save(outFeatureFile);
                System.err.println("Saved output feature file " + outFeatureFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                System.err.println("Error saving feature file");
            }
        }
    }


    /**
     * TODO: make this more generic, support arrays with many runs
     */
    public void analyzeMs1()
    {


        int num1Only = 0;
        int num2Only = 0;
        int numMatched = 0;
        int[] numRowsWithX = new int[runNames.size() + 1];
        List<Pair<Double,Double>> intensityPairs = new ArrayList<Pair<Double, Double>>();
        for (Object rowObj : rowMaps)
        {
            HashMap rowMap = (HashMap) rowObj;

            int numRunsWithFeature = 0;
            for (String runName : runNames)
                if (rowMap.get("intensity_" + runName) != null)
                    numRunsWithFeature++;
            numRowsWithX[numRunsWithFeature]++;

            Object intensityObject1 = rowMap.get("intensity_" + runNames.get(0));
            if (intensityObject1 != null)
            {
                Double intensity1;
                try
                {
                    intensity1 = (Double) intensityObject1;
                }
                catch (ClassCastException cce)
                {
                    intensity1 = Double.parseDouble((String) intensityObject1);
                }
                Object intensityObject2 = rowMap.get("intensity_" + runNames.get(1));
                if (intensityObject2 != null)
                {
                    Double intensity2;

                    try
                    {
                        intensity2 = (Double) intensityObject2;
                    }
                    catch (ClassCastException cce)
                    {
                        intensity2 = Double.parseDouble((String) intensityObject2);
                    }
                    numMatched++;
                    intensityPairs.add(new Pair<Double,Double>(intensity1, intensity2));
                }
                else
                    num1Only++;
            }
            else
                if (rowMap.get("intensity_" + runNames.get(1)) != null)
                    num2Only++;
        }

        System.err.println("\tNumber of rows with X features:");
        for (int i=0; i<numRowsWithX.length; i++)
            System.err.println("\t" + i + ":\t" + numRowsWithX[i]);


        System.err.println("\tTwo-run stuff:");
        ApplicationContext.infoMessage("MS1 2-run SUMMARY:");
        ApplicationContext.infoMessage("Number matched: " + numMatched);
        ApplicationContext.infoMessage("run " + runNames.get(0) + " only: " + num1Only);
        ApplicationContext.infoMessage("run " + runNames.get(0) + " total: " + (num1Only + numMatched));
        ApplicationContext.infoMessage("run " + runNames.get(1) + " only: " + num2Only);
        ApplicationContext.infoMessage("run " + runNames.get(1) + " total: " + (num2Only + numMatched));

        ApplicationContext.infoMessage("");

        double[] intensities1 = new double[intensityPairs.size()];
        double[] intensities2 = new double[intensityPairs.size()];
        for (int i=0; i<intensityPairs.size(); i++)
        {
            intensities1[i] = intensityPairs.get(i).first;
            intensities2[i] = intensityPairs.get(i).second;
        }
        ScatterPlotDialog spd =
                new ScatterPlotDialog(intensities1, intensities2, "Matched intensities");
        spd.setVisible(true);
        ApplicationContext.infoMessage("Matched intensity correlation coeff: " +
                BasicStatistics.correlationCoefficient(intensities1, intensities2));
    }


    public void analyzeMs2()
    {
        int numWithPeptideID = 0;
        int numWithConflict = 0;
        int[] peptideIDBreakdown = new int[runNames.size()+1];

        Map<String,Integer> peptideRowOccurrenceMap =
                new HashMap<String,Integer>();

        for (Map<String,Object> rowMap : rowMaps)
        {
            Set<String> peptidesThisRow = new HashSet<String>();

            int numRunsWithIDThisRow=0;

            String thisRowPeptide = null;
            boolean thisRowHasPeptide = false;
            boolean thisRowHasConflict = false;

            for (String runName : runNames)
            {
                String peptide = (String) rowMap.get("peptide_" + runName);
                if (peptide != null)
                {
                    peptidesThisRow.add(peptide);
                    numRunsWithIDThisRow++;
                    if (thisRowPeptide == null)
                    {
                        thisRowPeptide = peptide;
                        thisRowHasPeptide = true;
                    }
                    else
                    {
                        if (!thisRowPeptide.equals(peptide))
                            thisRowHasConflict = true;
                    }
                }
            }
            if (thisRowHasPeptide)
            {
                for (String peptideThisRow : peptidesThisRow)
                {
                    if (peptideRowOccurrenceMap.containsKey(peptideThisRow))
                    {
                        peptideRowOccurrenceMap.put(peptideThisRow,
                                peptideRowOccurrenceMap.get(peptideThisRow) + 1);
                    }
                    else
                    {
                        peptideRowOccurrenceMap.put(peptideThisRow,1);
                    }
                }
                numWithPeptideID++;
                if (thisRowHasConflict)
                    numWithConflict++;
                else
                    peptideIDBreakdown[numRunsWithIDThisRow] ++;
            }

        }

        int maxOccurrences = 0;
        int[] peptideOccurrenceArray = new int[1000];
        for (int occurrenceNum : peptideRowOccurrenceMap.values())
        {
            peptideOccurrenceArray[occurrenceNum]++;
            if (occurrenceNum > maxOccurrences)
                maxOccurrences++;
        }

        ApplicationContext.infoMessage("Number of peptides occurring on X different rows:");
        for (int i=1; i<=maxOccurrences; i++)
        {
            ApplicationContext.infoMessage(i + ": " + peptideOccurrenceArray[i]);
        }
        ApplicationContext.infoMessage("");

        ApplicationContext.infoMessage("Total rows: " + rowMaps.length);
        ApplicationContext.infoMessage("Rows with a peptide ID " + numWithPeptideID);
        ApplicationContext.infoMessage("Rows with a peptide ID conflict " + numWithConflict);
        ApplicationContext.infoMessage("Rows with X agreeing peptides and no conflict:");
        for (int i=1; i<peptideIDBreakdown.length; i++)
        {
            ApplicationContext.infoMessage("  " + i + ": " + peptideIDBreakdown[i]);
        }




    }



    public double[] compareIntensities(boolean add1)
    {
        double ratioSum = 0;

        List<Pair<Double,Double>> caseControlIntensityPairs =
                new ArrayList<Pair<Double,Double>>();

        for (Map<String,Object> rowMap : rowMaps)
        {
            double intensitySumCase = 0;
            for (String caseRunName : caseRunNames)
            {
                Object runIntensity = rowMap.get("intensity_" + caseRunName);
                if (runIntensity == null)
                    continue;
                intensitySumCase += Double.parseDouble(runIntensity.toString());
            }

            if (intensitySumCase == 0 && !add1)
                continue;

            double intensitySumControl = 0;
            for (String controlRunName : controlRunNames)
            {
                Object runIntensity = rowMap.get("intensity_" + controlRunName);
                if (runIntensity == null)
                    continue;
                intensitySumControl += Double.parseDouble(runIntensity.toString());
            }

            if (intensitySumControl == 0 && !add1)
                continue;

            if (add1)
            {
                if (intensitySumCase == 0 && intensitySumControl == 0)
                    continue;
                intensitySumCase++;
                intensitySumControl++;
            }

            caseControlIntensityPairs.add(new Pair<Double,Double>(intensitySumCase, intensitySumControl));
            ratioSum += intensitySumCase / intensitySumControl;
        }

        int numCommonPeptides = caseControlIntensityPairs.size();

        double[] caseIntensities = new double[numCommonPeptides];
        double[] controlIntensities = new double[numCommonPeptides];
        double[] caseControlIntensityRatios = new double[numCommonPeptides];

        for (int i=0; i<numCommonPeptides; i++)
        {
            caseIntensities[i] = caseControlIntensityPairs.get(i).first;
            controlIntensities[i] = caseControlIntensityPairs.get(i).second;
            caseControlIntensityRatios[i] = caseIntensities[i] / controlIntensities[i];
        }

        double[] logCaseIntensities = new double[numCommonPeptides];
        double[] logControlIntensities = new double[numCommonPeptides];


//double[] martyStat = new double[numCommonPeptides];
        for (int i=0; i<numCommonPeptides; i++)
        {
            logCaseIntensities[i] = Math.log(caseControlIntensityPairs.get(i).first);
            logControlIntensities[i] = Math.log(caseControlIntensityPairs.get(i).second);

//martyStat[i] = Math.abs((logCaseIntensities[i] - logControlIntensities[i]) / ((logCaseIntensities[i] + logControlIntensities[i]) / 2));
        }


        double meanRatio = ratioSum / caseIntensities.length;


        ApplicationContext.infoMessage("Peptides in common: " + caseIntensities.length);
        ScatterPlotDialog spd2 = new ScatterPlotDialog(logControlIntensities, logCaseIntensities, "X is log control, y is log case");
        spd2.setVisible(true);
//try {spd2.saveChartToImageFile(new File("/home/dhmay/temp/intensityplot.png"));} catch (Exception e) {ApplicationContext.errorMessage("error:",e);}

        if (showCharts)
        {
            ApplicationContext.infoMessage("Average case-control intensity ratio: " + meanRatio);
            ScatterPlotDialog spd = new ScatterPlotDialog(controlIntensities, caseIntensities, "X is control, y is case");            
            spd.setVisible(true);
            spd2.setVisible(true);

        }
//System.err.println("Mean Marty Stat: " + BasicStatistics.mean(martyStat));
//System.err.println("Median Marty Stat: " + BasicStatistics.median(martyStat));

        ApplicationContext.infoMessage("correlation coefficient of log intensities:" +  BasicStatistics.correlationCoefficient(logCaseIntensities, logControlIntensities));

        return caseControlIntensityRatios;
    }


    /**
     * Compare mean intensities of certain columns against each other in rows
     * in which the same peptide is identified in enough runs.
     *
     * A lot of this is hardcoded for a particular purpose right now.
     */
    public Map<String, Double> compareIntensitiesSamePeptide(double minSignificantRatio)
    {
        int minPeptideSupport = 1;
        int minFeatureSupport = 1;

        Set<String> peptidesHigherInCase =
                new HashSet<String>();
        Set<String> peptidesHigherInControl =
                new HashSet<String>();
        int rowsHigherInControl = 0;
        int rowsHigherInCase = 0;

        Map<String, Double> result = new HashMap<String,Double>();

        _log.debug("compareIntensitiesSamePeptide 1");

        try
        {
            List<Double> intensitiesCase = new ArrayList<Double>();
            List<Double> intensitiesControl = new ArrayList<Double>();
            int numPeptidesInAgreement = 0;


            for (Map<String,Object> rowMap : rowMaps)
            {
                int featureSupportCase = 0;
                int peptideSupportCase = 0;
                String peptide = null;
                double intensitySumCase = 0;

                for (String caseRunName : caseRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + caseRunName);
                    if (runIntensity == null)
                        continue;
                    intensitySumCase += Double.parseDouble(runIntensity.toString());
                    featureSupportCase++;
                    String thisPeptide = (String) rowMap.get("peptide_" + caseRunName);
                    if (thisPeptide == null)
                        continue;
                    if (peptide == null)
                    {
                        peptide = thisPeptide;
                        peptideSupportCase++;
                    }
                    else if (peptide.equals(thisPeptide))
                        peptideSupportCase++;
                    else
                        peptideSupportCase =0;
                        break;
                }
                if (featureSupportCase < minFeatureSupport)
                    continue;

                double intensityMeanCase = intensitySumCase / peptideSupportCase;
                double intensitySumControl = 0;

                int featureSupportControl = 0;
                int peptideSupportControl = 0;
                for (String controlRunName : controlRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + controlRunName);
                    if (runIntensity == null)
                        continue;
                    intensitySumControl += Double.parseDouble(runIntensity.toString());
                    featureSupportControl++;
                    String thisPeptide = (String) rowMap.get("peptide_" + controlRunName);
                    if (peptide == null || thisPeptide == null)
                        continue;

                    if (peptide.equals(thisPeptide))
                        peptideSupportControl++;
                    else
                        peptideSupportControl = 0;
                        break;
                }

                if (featureSupportControl < minFeatureSupport)
                    continue;

                boolean peptideAgreement = false;
                if (peptideSupportControl >= minPeptideSupport && peptideSupportCase >= minPeptideSupport)
                {
                    numPeptidesInAgreement++;
                    peptideAgreement=true;
                }

                double intensityMeanControl = intensitySumControl / peptideSupportControl;

                double caseControlRatio = intensityMeanCase / intensityMeanControl;
                if (caseControlRatio > minSignificantRatio)
                {
                    rowsHigherInCase++;
                    if (peptideAgreement)
                        peptidesHigherInCase.add(peptide);
                }
                else if (1 / caseControlRatio > minSignificantRatio)
                {
                    rowsHigherInControl++;
                    if (peptideAgreement)
                        peptidesHigherInControl.add(peptide);
                }

                if (peptideAgreement)
                {
//                    _log.debug("Agreeing peptide " + peptide);

                    if (peptide.contains("["))
                    {
                        _log.debug("Stripping brackets");
                        peptide = peptide.substring(peptide.indexOf("[") +1, peptide.indexOf("]"));
                    }
                    if (peptide.contains(";"))
                    {
                        _log.debug("Multiple peptides, taking first one arbitrarily");

                        peptide = peptide.substring(0, peptide.indexOf(";"));
                    }

                    if (result.containsKey(peptide))
                    {
                        _log.debug("Augmenting peptide " + peptide + " with second row");
                        double existingIntensity = result.get(peptide);
                        result.put(peptide, existingIntensity + (intensityMeanCase / intensityMeanControl));
                    }
                    else
                        result.put(peptide, intensityMeanCase / intensityMeanControl);
                }

                intensitiesCase.add(intensityMeanCase);
                intensitiesControl.add(intensityMeanControl);
            }

//            File upregulatedPeptidesHtmlFile =
//                    TempFileManager.createTempFile("peptides_upregulated_in_cases.html",this);
//            FileOutputStream fos = new FileOutputStream(upregulatedPeptidesHtmlFile);
//            StringBuffer htmlToWrite = new StringBuffer();
//            htmlToWrite.append("<html>");
//            htmlToWrite.append("<head><title>" + "Peptides Upregulated in Cases" + "</head></title>");
//            htmlToWrite.append("</body>");
////            fos.write(htmlToWrite.toString().getBytes());
////            htmlToWrite = new StringBuffer();
//            htmlToWrite.append("<table>");
//            for (String pep : peptidesHigherInCase)
//                htmlToWrite.append("<tr><td>" + pep + "</td></tr>");
//            htmlToWrite.append("</table>");
////            fos.write(htmlToWrite.toString().getBytes());
//
////            htmlToWrite = new StringBuffer();
//            htmlToWrite.append("</body></html>");
////            fos.write(htmlToWrite.toString().getBytes());
////            fos.flush();


//            System.err.println("Peptides higher in case:");
//            for (String pep : peptidesHigherInCase)
//                System.err.println(pep);

            System.err.println("# Peptides in agreement: " + numPeptidesInAgreement);
            System.err.println("# Peptides higher in case: " + peptidesHigherInCase.size() + " out of " + rowsHigherInCase + " array rows");
            System.err.println("# Peptides higher in control: " + peptidesHigherInControl.size() + " out of " + rowsHigherInControl + " array rows");
            System.err.println("# Peptides higher in one or the other: " + (peptidesHigherInCase.size() + peptidesHigherInControl.size() +
                               " out of " + (rowsHigherInCase + rowsHigherInControl) + " array rows"));


            double[] intensitiesCasearray = new double[intensitiesCase.size()];
            double[] intensitiesControlarray = new double[intensitiesControl.size()];
int numInsideTwofold=0;
            for (int i=0; i<intensitiesCasearray.length; i++)
            {
                intensitiesCasearray[i] = intensitiesCase.get(i);
                intensitiesControlarray[i] = intensitiesControl.get(i);
double intensitiesRatio = intensitiesCasearray[i] / intensitiesControlarray[i];
if (intensitiesRatio > 0.5 && intensitiesRatio < 2.0)
    numInsideTwofold++;
            }

            if (showCharts)
            {
                ScatterPlotDialog spd = new ScatterPlotDialog();

                spd.addData(intensitiesControlarray, intensitiesCasearray, "intensities");
                spd.setVisible(true);
            }
            ApplicationContext.infoMessage("Same peptide intensity summary:");
System.err.println("Within twofold: " + numInsideTwofold + " out of " + intensitiesCasearray.length);


        }
        catch (Exception e)
        {
            System.err.println("Exception: " + e.getMessage());
            e.printStackTrace(System.err);
        }
        return result;
    }

    public Map<String, Object>[] getRowMaps()
    {
        return rowMaps;
    }

    public List<String> getRunNames()
    {
        return runNames;
    }

    public String[] getCaseRunNames()
    {
        return caseRunNames;
    }

    public String[] getControlRunNames()
    {
        return controlRunNames;
    }

    public void setControlRunNames(String[] controlRunNames)
    {
        this.controlRunNames = controlRunNames;
    }

    public void setCaseRunNames(String[] caseRunNames)
    {
        this.caseRunNames = caseRunNames;
    }

    public boolean shouldShowCharts()
    {
        return showCharts;
    }

    public void setShowCharts(boolean showCharts)
    {
        this.showCharts = showCharts;
    }
}
