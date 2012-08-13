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
package org.fhcrc.cpl.viewer.align;

import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureAsMap;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
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

    protected List<String> columnNames = null;
    protected Map<String,Object>[] rowMaps = null;
    protected List<String> runNames = null;

    protected String[] caseRunNames;
    protected String[] controlRunNames;

    //Indexes of case and control runs.
    double[] caseRunIndexes1Based;
    double[] controlRunIndexes1Based;

    protected boolean showCharts = false;


    public static final int CONSENSUS_INTENSITY_MODE_MEAN = 0;
    public static final int CONSENSUS_INTENSITY_MODE_FIRST = 1;
    public static final int CONSENSUS_INTENSITY_MODE_MEDIAN = 2;
    public static final int CONSENSUS_INTENSITY_MODE_SUM = 3;
    public static final int CONSENSUS_INTENSITY_MODE_MAX = 4;


    public static final String[] INTENSITY_MODE_STRINGS =
            {
                    "mean",
                    "first",
                    "median",
                    "sum",
                    "max"
            };

    public static final String[] INTENSITY_MODE_EXPLANATIONS = 
            {
                    "mean",
                    "first",
                    "median",
                    "sum",
                    "max"
            };

    protected  Map<Integer, List<Map<String, Object>>> idDetailsRowsMap = null;


    public static final float MAX_Q_VALUE = 0.05f;
    protected float maxQValue = MAX_Q_VALUE;
    protected File outLowQValueArrayFile;
    protected File outLowQValueAgreeingPeptidePepXMLFile;
    protected File detailsFile = null;

    //need a reference for this in order to do things with the array in R
    protected File arrayFile = null;

    protected boolean arrayHasMs2 = false;

    public PeptideArrayAnalyzer(File arrayFile)
            throws IOException
    {
        parseArray(arrayFile);
        this.arrayFile = arrayFile;
        //default details file location.  This can be overridden
        detailsFile = new File(BucketedPeptideArray.calcDetailsFilepath(arrayFile.getAbsolutePath()));        
    }

    public static List<String> loadArrayRunNames(File arrayFile) throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(arrayFile));
        String line = null;
        while ((line = br.readLine()).startsWith("#"))
            continue;
        String[] chunks = line.split("\t");
        List<String> result = new ArrayList<String>();

        for (String chunk : chunks)
            if (chunk.startsWith("intensity_"))
                result.add(chunk.substring("intensity_".length()));
        br.close();
        return result;
    }

    public void parseArray(File arrayFile) throws IOException
    {
        TabLoader tabLoader = new TabLoader(arrayFile);
        tabLoader.setReturnElementClass(HashMap.class);
        Object[] rows = tabLoader.load();

        runNames = new ArrayList<String>();
        columnNames  = new ArrayList<String>();

        for (TabLoader.ColumnDescriptor column : tabLoader.getColumns())
        {
            _log.debug("loading column " + column.name);
            columnNames.add(column.name);
            if (column.name.startsWith("intensity_"))
            {
                runNames.add(column.name.substring("intensity_".length()));
                _log.debug("adding run " + runNames.get(runNames.size() - 1));
            }
            if (column.name.startsWith("peptide"))
                arrayHasMs2 = true;
        }

        rowMaps = (Map<String, Object>[]) rows;
    }

    /**
     *
     * @param rowMap
     * @return an empty set if none
     */
    public Set<String> getAllRowPeptides(Map<String, Object> rowMap)
    {
        Set<String> result = new HashSet<String>();
        for (String runName :runNames)
        {
            String peptide = this.getRunPeptide(rowMap, runName);
            if (peptide != null)
                result.add(peptide);
        }
        return result;
    }

    /**
     * Build a map from peptides to the indices of rows that contain them.  If a row has multiple peptides,
     * they both get an entry.
     * @return
     */
    public Map<String, List<Integer>> loadPeptideRowIndexesMap()
    {
        Map<String, List<Integer>> result = new HashMap<String, List<Integer>>();

        for (int i=0; i<getRowCount(); i++)
        {
            Map<String, Object> rowMap = rowMaps[i];
            for (String peptide : getAllRowPeptides(rowMap))
            {
                List<Integer> rowsThisPeptide = result.get(peptide);
                if (rowsThisPeptide == null)
                {
                    rowsThisPeptide = new ArrayList<Integer>();
                    result.put(peptide, rowsThisPeptide);
                }
                rowsThisPeptide.add(i);
            }
        }
        return result;
    }

    public int getRowCount()
    {
        return rowMaps.length;
    }

    /**
     * Requires exactly one peptide in the row, and at least one case and one control intensity
     * @return
     */
    public List<Pair<String, Pair<List<Double>, List<Double>>>> loadPeptidesCaseControlIntensitiesPerRow()
    {
        List<Pair<String, Pair<List<Double>, List<Double>>>> result = new ArrayList<Pair<String, Pair<List<Double>, List<Double>>>>();
        for (Map<String, Object> rowMap : rowMaps)
        {
            Set<String> peptides = getAllRowPeptides(rowMap);
            if (!(peptides.size() == 1))
                continue;
            List<Double> caseIntensities = new ArrayList<Double>();
            List<Double> controlIntensities = new ArrayList<Double>();

            for (String caseRunName : caseRunNames)
            {
                Double intensity = getRunIntensity(rowMap,caseRunName);
                if (intensity != null) caseIntensities.add(intensity);
            }
            for (String controlRunName : controlRunNames)
            {
                Double intensity = getRunIntensity(rowMap,controlRunName);
                if (intensity != null) controlIntensities.add(intensity);
            }
            if (caseIntensities.isEmpty() || controlIntensities.isEmpty())
                continue;

            String peptide = (String)peptides.toArray()[0];
            Pair<String, Pair<List<Double>, List<Double>>> resultRow = new Pair<String, Pair<List<Double>, List<Double>>>(
                    peptide, new Pair<List<Double>, List<Double>>(caseIntensities, controlIntensities));
            result.add(resultRow);
        }
        return result;
    }

    public void loadDetailsMap() throws IOException
    {
        _log.debug("Loading details, file=" + detailsFile.getAbsolutePath());        
        TabLoader tabLoader = new TabLoader(detailsFile);
        tabLoader.setReturnElementClass(HashMap.class);
        Object[] rows = tabLoader.load();

        runNames = new ArrayList<String>();
        columnNames  = new ArrayList<String>();

        for (TabLoader.ColumnDescriptor column : tabLoader.getColumns())
        {
            _log.debug("loading column " + column.name);
            columnNames.add(column.name);
            if (column.name.startsWith("intensity_"))
            {
                runNames.add(column.name.substring("intensity_".length()));
                _log.debug("adding run " + runNames.get(runNames.size() - 1));
            }
        }

        Map<String, Object>[] detailsRowMaps = (Map<String, Object>[]) rows;
        idDetailsRowsMap = new HashMap<Integer, List<Map<String, Object>>>();
        for (Map<String, Object> row : detailsRowMaps)
        {
            int id = Integer.parseInt(row.get("id").toString());
            List<Map<String, Object>> rowMaps = idDetailsRowsMap.get(id);
            if (rowMaps == null)
            {
                rowMaps = new ArrayList<Map<String, Object>>();
                idDetailsRowsMap.put(id, rowMaps);
            }
            rowMaps.add(row);
        }
    }

    public Double getRunIntensity(Map<String, Object> rowMap, String runName)
    {
        Double result = null;
        Object intensityObj = rowMap.get("intensity_" + runName);
        if (intensityObj != null)
            result = Double.parseDouble(intensityObj.toString());
        return result;
    }

    /**
     * Return the "column" of intensity values for a specified run.  This is Double objects rather than
     * doubles so there can be nulls
     * @param runName
     * @return
     */
    public Double[] getRunIntensities(String runName)
    {
        Double[] result = new Double[rowMaps.length];
        for (int i=0; i<result.length; i++)
            result[i] = getRunIntensity(rowMaps[i], runName);
        return result;
    }

    public String getRunPeptide(Map<String, Object> rowMap, String runName)
    {
        return (String) rowMap.get("peptide_" + runName);
    }


    public int countConflictRows()
    {
        int conflictRows = 0;

        for (Map<String,Object> rowMap : rowMaps)
        {
            String rowPeptide = null;
            for (String runName : runNames)
            {
                String runPeptide = getRunPeptide(rowMap, runName);
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

    public Map<String, List<Feature>> createFileFeatureMapForRow(int rowId)
    {
        Map<String, List<Feature>> result = new HashMap<String, List<Feature>>();
        List<Map<String, Object>> detailsRows = idDetailsRowsMap.get(rowId);
        if (detailsRows != null)
        {
            for (Map<String, Object> detailsRow : detailsRows)
            {
                try {
                    Feature feature = createFeatureFromDetailsFileRow(detailsRow);
                    String fileString = (String) detailsRow.get("file");
                    List<Feature> featuresThisFile = result.get(fileString);
                    if (featuresThisFile == null)
                    {
                        featuresThisFile = new ArrayList<Feature>();
                        result.put(fileString, featuresThisFile);
                    }
                    featuresThisFile.add(feature);
                }  catch (NullPointerException e) {
                    ApplicationContext.errorMessage("Failed to create feature for detailsRow with id " + rowId +
                            ". Null? " + (detailsRow == null),e);
                    throw e;

                }
            }
        }
        return result;
    }

    public Feature createFeatureFromDetailsFileRow(Map<String, Object> currentDetailsRow)
    {
        int scan=-1, scanFirst=-1, scanLast=-1,charge=-1;
        float mz=-1, intensity=-1, kl=-1, totalIntensity=-1;


        try {
            scan =  (Integer) currentDetailsRow.get("scan");
            scanFirst = (Integer) currentDetailsRow.get("scanFirst");
            scanLast =        (Integer) currentDetailsRow.get("scanLast");
            mz =        ((Double) currentDetailsRow.get("mz")).floatValue();
            intensity =        ((Double) currentDetailsRow.get("intensity")).floatValue();
            charge =         (Integer) currentDetailsRow.get("charge");
            kl =        ((Double) currentDetailsRow.get("kl")).floatValue();
            totalIntensity =        ((Double) currentDetailsRow.get("totalIntensity")).floatValue();
            Feature feature = new Feature(
                    scan, scanFirst, scanLast, mz, intensity, charge, kl, totalIntensity);
            feature.setTime(((Double) currentDetailsRow.get("time")).floatValue());
            if (currentDetailsRow.get("peptide") != null)
                MS2ExtraInfoDef.setPeptideList(feature, (String) currentDetailsRow.get("peptide"));
//dmtest
            if (currentDetailsRow.get("peptideprophet") != null)
                MS2ExtraInfoDef.setPeptideProphet(feature, Double.valueOf(currentDetailsRow.get("peptideprophet").toString()));
            if (currentDetailsRow.get("protein") != null)
                MS2ExtraInfoDef.setProteinList(feature, (String) currentDetailsRow.get("protein"));
            if (currentDetailsRow.get("modifiedaminoacids") != null)
                MS2ExtraInfoDef.setModifiedAminoAcids(feature,
                        MS2ExtraInfoDef.parsePositionModifiedAminoAcidListMapString(
                                (String)currentDetailsRow.get("modifiedaminoacids")));
            return feature;

        } catch (NullPointerException e) {
            ApplicationContext.infoMessage("ERROR: NPE while trying to create feature. Attributes: " +
                    "scan=" + scan + ", scanFirst=" + scanFirst + ", scanLast=" + scanLast + ", mz=" + mz +
                    "intensity=" + intensity + ", charge=" + charge + ", kl=" + kl + ", totalIntensity=" + totalIntensity + "\ncols:");
            for (String key : currentDetailsRow.keySet()) ApplicationContext.infoMessage("\t" + key);
            throw e;
        }

    }

    /**
     * Loads a map from row ID to a map from run name to features
     * @return
     * @throws IOException
     */
    public  Map<Integer, Map<String, List<Feature>>> loadDetailsRowMapsById()
            throws IOException
    {
        TabLoader detailsTabLoader = new TabLoader(detailsFile);
        detailsTabLoader.setReturnElementClass(HashMap.class);
        Object[] detailRows = detailsTabLoader.load();
        Map<String, Object>[] detailsRowMaps = (Map<String, Object>[]) detailRows;
        Map<Integer, Map<String, List<Feature>>> result = new HashMap<Integer, Map<String, List<Feature>>>();
        for (Map<String, Object> detailsRow : detailsRowMaps)
        {
            int id = Integer.parseInt(detailsRow.get("id").toString());
            Map<String, List<Feature>> runFeatureMapThisId = result.get(id);
            if (runFeatureMapThisId == null)
            {
                runFeatureMapThisId = new HashMap<String, List<Feature>>();
                result.put(id, runFeatureMapThisId);
            }
            String fileName = (String) detailsRow.get("file");
            String runName =
                    fileName.substring(0, fileName.indexOf("."));
            Feature feature = createFeatureFromDetailsFileRow(detailsRow);

            List<Feature> featuresThisIdThisRun = runFeatureMapThisId.get(runName);
            if (featuresThisIdThisRun == null)
            {
                featuresThisIdThisRun = new ArrayList<Feature>();
                runFeatureMapThisId.put(runName, featuresThisIdThisRun);
            }
            featuresThisIdThisRun.add(feature);
        }
        return result;
    }


    /**
     * Create feature files for features that occur in all runs, with no conflicts
     *
     * todo: make this use loadDetailsRowMapsById.  It's less efficient, but it will simplify things

     *
     * @param detailsFile
     * @param minConsensusRuns
     * @param intensityMode
     * @param requireSamePeptide
     * @param shouldAllowNegScanAndTime if this is false, and a feature has scan<1 or time<0, peg to those values
     * @return
     * @throws IOException
     */
    public FeatureSet createConsensusFeatureSet(File detailsFile,
                                           int minConsensusRuns,
                                           int intensityMode, boolean requireSamePeptide,
                                           boolean shouldAllowNegScanAndTime)
            throws IOException
    {

        this.detailsFile = detailsFile;
        TabLoader detailsTabLoader = new TabLoader(detailsFile);
        detailsTabLoader.setReturnElementClass(HashMap.class);
        Object[] detailRows = detailsTabLoader.load();
        Map<String, Object>[] detailsRowMaps = (Map<String, Object>[]) detailRows;

        boolean hasPeptides = false;
        for (TabLoader.ColumnDescriptor colDes : detailsTabLoader.getColumns())
            if (colDes.name.equals("peptide"))
            {
                hasPeptides = true;
                break;
            }
        if (requireSamePeptide && !hasPeptides)
            throw new IOException("Array details file " + detailsFile.getAbsolutePath() + " does not have peptide information");

        Map<String, List<Feature>> runFeatureLists =
                new HashMap<String, List<Feature>>(runNames.size());
        for (String runName : runNames)
            runFeatureLists.put(runName, new ArrayList<Feature>());

        int currentDetailsRowId = -1;
        int currentDetailsArrayIndex = -1;
        Map<String, Object> currentDetailsRow = null;

        List<Feature> resultFeatureList = new ArrayList<Feature>();

        //for reporting intensity correlation, iff there are 2 runs
        List<Float> intensities1 = new ArrayList<Float>();
        List<Float> intensities2 = new ArrayList<Float>();


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
                    Feature feature = createFeatureFromDetailsFileRow(currentDetailsRow);


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

                //for 2-run intensity correlation
                if (thisRowRunFeatureMap.size() == 2)
                {
                    List<String> runs = new ArrayList<String>(thisRowRunFeatureMap.keySet());
                    if (thisRowRunFeatureMap.get(runs.get(0)).size() == 1 &&
                        thisRowRunFeatureMap.get(runs.get(1)).size() == 1)
                    {
                        intensities1.add(thisRowRunFeatureMap.get(runs.get(0)).get(0).getIntensity());
                        intensities2.add(thisRowRunFeatureMap.get(runs.get(1)).get(0).getIntensity());
                    }
                }




                List<Feature> allFeaturesThisRow = new ArrayList<Feature>();

                    //todo: NOT restricting to those runs that only have one feature.  OK?

                for (String runName : thisRowRunFeatureMap.keySet())
                {
                    allFeaturesThisRow.addAll(thisRowRunFeatureMap.get(runName));
                    runFeatureLists.get(runName).addAll(thisRowRunFeatureMap.get(runName));
                }

                List<Double> thisFeatureIntensities = new ArrayList<Double>();
                Set<String> thisFeaturePeptides = new HashSet<String>();

                Collections.sort(allFeaturesThisRow, new Feature.IntensityDescComparator());

                //begin with the most-intense feature
                Feature mostIntenseFeature = allFeaturesThisRow.get(0);

                //update first and last scans based on everything we're glomming together. But leave max scan,
                //time and intensity alone
                Feature firstFeature = mostIntenseFeature;
                for (Feature feature : allFeaturesThisRow) {
                    if (feature.scan < firstFeature.scan)
                        firstFeature = feature;
                    mostIntenseFeature.scanFirst = Math.min(feature.scanFirst, mostIntenseFeature.scanFirst);
                    mostIntenseFeature.scanLast = Math.min(feature.scanLast, mostIntenseFeature.scanLast);

                    thisFeatureIntensities.add((double) feature.getIntensity());
                    if (requireSamePeptide && (MS2ExtraInfoDef.getPeptideList(feature) != null))
                        thisFeaturePeptides.addAll(MS2ExtraInfoDef.getPeptideList(feature));
                }
                Feature consensusFeature = mostIntenseFeature;
                if (consensusFeature == null) continue;

                //set description to row id, unless it already has a non-empty value
                if (consensusFeature.getDescription() == null || consensusFeature.getDescription().length() == 0)
                    consensusFeature.setDescription("" + rowId);
                if (!shouldAllowNegScanAndTime && consensusFeature.getScan() < 1)
                    consensusFeature.setScan(1);
                if (!shouldAllowNegScanAndTime && consensusFeature.getTime() < 0)
                    consensusFeature.setTime(0);


                if (requireSamePeptide && (thisFeaturePeptides.size() > 1))
                {
                    _log.debug("SKIPPING: features have " + thisFeaturePeptides.size() + " unique peptides");
                    continue;
                }
                switch (intensityMode)
                {
                    case CONSENSUS_INTENSITY_MODE_MEAN:
                        double[] featureIntensitiesArray = new double[thisFeatureIntensities.size()];
                        for (int j=0; j<thisFeatureIntensities.size(); j++)
                            featureIntensitiesArray[j] = thisFeatureIntensities.get(j);
                        consensusFeature.setIntensity((float) BasicStatistics.mean(featureIntensitiesArray));
                        break;
                    case CONSENSUS_INTENSITY_MODE_FIRST:
                        consensusFeature.setIntensity(firstFeature.getIntensity());
                        break;
                    case CONSENSUS_INTENSITY_MODE_MEDIAN:
                        double[] featureIntensitiesArray2 = new double[thisFeatureIntensities.size()];
                        for (int j=0; j<thisFeatureIntensities.size(); j++)
                            featureIntensitiesArray2[j] = thisFeatureIntensities.get(j);
                        consensusFeature.setIntensity((float) BasicStatistics.median(featureIntensitiesArray2));
                        break;
                    case CONSENSUS_INTENSITY_MODE_SUM:
                        consensusFeature.setIntensity((float) BasicStatistics.sum(thisFeatureIntensities));
                        break;
                    case CONSENSUS_INTENSITY_MODE_MAX:
                        consensusFeature.setIntensity((float) BasicStatistics.max(thisFeatureIntensities));
                        break;
                }
                resultFeatureList.add(consensusFeature);

            }
        }

        if (! intensities1.isEmpty() && !intensities2.isEmpty())
        {
            List<Float> cvs = new ArrayList<Float>();
            for (int i=0; i<intensities1.size(); i++)
            {
                cvs.add((float) BasicStatistics.coefficientOfVariation(new double[] {intensities1.get(i), intensities2.get(i)}));
            }
            ApplicationContext.infoMessage("2-run feature intensity correlation: " +
                    BasicStatistics.correlationCoefficient(intensities1, intensities2) + ", Mean CV: " +
                    BasicStatistics.mean(cvs) + ", Median CV: " + BasicStatistics.median(cvs));
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
        int[] numRowsWithXRuns = new int[runNames.size() + 1];
        int mostFeaturesInARow = 0;
        int[] numRowsWithXFeatures = new int[10 * runNames.size()];
        List<Integer> numRunsPerRow = new ArrayList<Integer>();
        List<Pair<Double,Double>> intensityPairs = new ArrayList<Pair<Double, Double>>();

        Map<String,List<Float>> boxPlotValueMap = new HashMap<String,List<Float>>();
        for (String runName : runNames)
            boxPlotValueMap.put(runName, new ArrayList<Float>());
        for (Object rowObj : rowMaps)
        {
            HashMap rowMap = (HashMap) rowObj;

            int numRunsWithFeature = 0;
            for (String runName : runNames)
            {
                if (rowMap.get("intensity_" + runName) != null)
                {
                    numRunsWithFeature++;
                    boxPlotValueMap.get(runName).add((float)Math.log(Double.parseDouble(rowMap.get("intensity_" + runName).toString())));
                }
            }
            numRowsWithXRuns[numRunsWithFeature]++;
            int featureCount = (Integer) rowMap.get("featureCount");
            mostFeaturesInARow = Math.max(mostFeaturesInARow, featureCount);
            numRowsWithXFeatures[featureCount]++;
            numRunsPerRow.add(numRunsWithFeature);

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
        int numTotalRows = rowMaps.length;
        System.err.println("Number of rows with features in X runs:");
        for (int i=1; i<numRowsWithXRuns.length; i++)
            System.err.println("\t" + i + ":\t" + numRowsWithXRuns[i] + " (" + ( numRowsWithXRuns[i] * 100f/ (float) numTotalRows) + "%)");
        System.err.println("Number of rows with X features:");        
        for (int i=1; i<=mostFeaturesInARow; i++)
            System.err.println("\t" + i + ":\t" + numRowsWithXFeatures[i] + " (" + ( numRowsWithXFeatures[i] * 100f/ (float) numTotalRows) + "%)");
        new PanelWithHistogram(numRunsPerRow, "# runs per row").displayInTab();

        if (runNames.size() == 2)
        {
            System.err.println("\tTwo-run array analysis:");
            ApplicationContext.infoMessage("MS1 2-run SUMMARY:");
            ApplicationContext.infoMessage("Number matched: " + numMatched);
            ApplicationContext.infoMessage("run " + runNames.get(0) + " only: " + num1Only);
            ApplicationContext.infoMessage("run " + runNames.get(0) + " total: " + (num1Only + numMatched));
            ApplicationContext.infoMessage("run " + runNames.get(1) + " only: " + num2Only);
            ApplicationContext.infoMessage("run " + runNames.get(1) + " total: " + (num2Only + numMatched));

            ApplicationContext.infoMessage("");

            double[] intensities1 = new double[intensityPairs.size()];
            double[] intensities2 = new double[intensityPairs.size()];
            double[] logIntensities1 = new double[intensityPairs.size()];
            double[] logIntensities2 = new double[intensityPairs.size()];
            for (int i=0; i<intensityPairs.size(); i++)
            {
                intensities1[i] = intensityPairs.get(i).first;
                intensities2[i] = intensityPairs.get(i).second;
                logIntensities1[i] = Math.log(intensities1[i]);
                logIntensities2[i] = Math.log(intensities2[i]);
            }
            PanelWithScatterPlot spd =
                    new PanelWithScatterPlot(intensities1, intensities2, "Matched intensities");
            spd.displayInTab();
            PanelWithScatterPlot spd2 =
                    new PanelWithScatterPlot(logIntensities1, logIntensities2, "Matched log-intensities");
            spd2.displayInTab();
            ApplicationContext.infoMessage("Matched intensity correlation coeff: " +
                    BasicStatistics.correlationCoefficient(intensities1, intensities2));
        }

        PanelWithBoxAndWhiskerChart boxPlot = new PanelWithBoxAndWhiskerChart("Log Intensity Boxplot");
        boxPlot.displayInTab();
        {
            //Counts for each run.  Not needed if two runs, since given above
            if (runNames.size() > 2)
                ApplicationContext.infoMessage("Features per run:");
            for (String runName : runNames)
            {
                List<Float> featureSetData = boxPlotValueMap.get(runName);
                boxPlot.addData(featureSetData, runName);
                    ApplicationContext.infoMessage("\t" + featureSetData.size());
            }
        }
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

    /**
     * Combine multiple features into a single feature.  The highest-intensity feature is used, and
     * it gets intensity and totalIntensities that are sums of its component features' intensities and
     * totalIntensities.  For peptide and protein assignments, if no conflicts, assign to the combined
     * feature.  If conflicts, don't.
     * @param features
     * @return
     */
    protected Feature combineFeaturesSumIntensities(List<Feature> features)
    {
        float intensitySum = 0;
        float intensityMax = -1;
        float totalIntensitySum = 0;
        Feature maxIntensityFeature = null;
        for (Feature feature : features)
        {
            intensitySum += feature.getIntensity();
            if (feature.getIntensity() > intensityMax)
            {
                intensityMax = feature.getIntensity();
                maxIntensityFeature = feature;
                intensitySum += feature.getIntensity();
                totalIntensitySum += feature.getTotalIntensity();
            }
        }
        maxIntensityFeature.setIntensity(intensitySum);
        maxIntensityFeature.setTotalIntensity(totalIntensitySum);


        Set<String> featurePeptides = new HashSet<String>();
        Set<String> featureProteins = new HashSet<String>();

        for (Feature f : features)
        {
            String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(f);
            if (featurePeptide != null)
            {
                featurePeptides.add(featurePeptide);

                String featureProtein = MS2ExtraInfoDef.getFirstProtein(f);
                if (featureProtein != null)
                    featureProteins.add(featureProtein);
            }
        }

        if (featurePeptides.size() == 1 &&
                MS2ExtraInfoDef.getFirstPeptide(maxIntensityFeature) == null)
        {
            MS2ExtraInfoDef.setSinglePeptide(maxIntensityFeature,
                    featurePeptides.iterator().next());

            if (featureProteins.size() == 1 &&
                    MS2ExtraInfoDef.getFirstProtein(maxIntensityFeature) == null)
                MS2ExtraInfoDef.addProtein(maxIntensityFeature,
                        featureProteins.iterator().next());
        }

        return maxIntensityFeature;
    }

    public static class FeatureAsMapWithIdAndFile extends FeatureAsMap
    {
        protected int id;
        protected File file;

        public FeatureAsMapWithIdAndFile()
        {
            
        }

        public int getId()
        {
            return id;
        }

        public void setId(int id)
        {
            this.id = id;
        }

        public File getFile()
        {
            return file;
        }

        public void setFile(File file)
        {
            this.file = file;
        }
    }

    protected List<Map<String, List<Feature>>> loadRunFeaturesFromDetailsArray()
    {
        try
        {
            TabLoader loader = new TabLoader(detailsFile, FeatureAsMapWithIdAndFile.class);

            Iterator it = loader.iterator();
            List<Map<String, List<Feature>>> result = new ArrayList<Map<String, List<Feature>>>();
            Map<String, List<Feature>> curMap = null;
            int curId = 0;
            while (it.hasNext())
            {
                FeatureAsMapWithIdAndFile feature = (FeatureAsMapWithIdAndFile) it.next();
                feature.afterPopulate();
                if (feature.getId() > curId)
                {
                    curMap = new HashMap<String, List<Feature>>();
                    result.add(curMap);
                    curId = feature.getId();
                }
                //todo: WARNING!  This code here means that we can't have "." in run names.  Bad.
                String runName = feature.getFile().getName();
                if (runName.contains("."))
                    runName = runName.substring(0, runName.indexOf("."));
                List<Feature> featuresThisRun = curMap.get(runName);
                if (featuresThisRun == null)
                {
                    featuresThisRun = new ArrayList<Feature>();
                    curMap.put(runName, featuresThisRun);
                }
                featuresThisRun.add(feature);
            }
            return result;
        }
        catch (IOException e)
        {
            throw new RuntimeException("Failure loading details array file " +
                detailsFile.getAbsolutePath(),e);
        }
    }


    /**
     *
     * @param minRunsPerGroup
     * @param showCharts
     */
    public void calcTScoresQValues(int minRunsPerGroup, boolean showCharts, boolean samePeptideOnly)
    {


        List<Map<String, List<Feature>>> idRunFeaturesMaps = loadRunFeaturesFromDetailsArray();
        List<Map<String, Feature>> runFeatureMapsIndexedByR = new ArrayList<Map<String, Feature>>();
        List<double[]> allCaseIntensities = new ArrayList<double[]>();
        List<double[]> allControlIntensities = new ArrayList<double[]>();
        List<Float> numMinFeaturesPerGroup = new ArrayList<Float>();

        ApplicationContext.infoMessage("Rows with peptides in array: " + idRunFeaturesMaps.size());
        for (Map<String, List<Feature>> runFeaturesMap : idRunFeaturesMaps)
        {
            double[] thisRowCase = new double[caseRunNames.length];
            double[] thisRowControl = new double[controlRunNames.length];
            Map<String, Feature> runSummaryFeatureMap = new HashMap<String, Feature>();

            int numCaseFeatures = 0;
            int numControlFeatures = 0;

            int runIndex = 0;
            for (String runName : caseRunNames)
            {
                thisRowCase[runIndex] = Double.NaN;
                List<Feature> featuresThisRun = runFeaturesMap.get(runName);
                if (featuresThisRun != null)
                {
                    Feature representativeFeature = combineFeaturesSumIntensities(featuresThisRun);
                    thisRowCase[runIndex] = representativeFeature.getIntensity();
                    runSummaryFeatureMap.put(runName, representativeFeature);
                    numCaseFeatures++;
                }
                runIndex++;
            }
            if (numCaseFeatures < minRunsPerGroup)
                continue;

            runIndex = 0;
            for (String runName : controlRunNames)
            {
                thisRowControl[runIndex] = Double.NaN;
                List<Feature> featuresThisRun = runFeaturesMap.get(runName);
                if (featuresThisRun != null)
                {
                    Feature representativeFeature = combineFeaturesSumIntensities(featuresThisRun);
                    thisRowControl[runIndex] = representativeFeature.getIntensity();
                    runSummaryFeatureMap.put(runName, representativeFeature);
                    numControlFeatures++;
                }
                runIndex++;
            }
            if (numControlFeatures < minRunsPerGroup)
                continue;
            allCaseIntensities.add(thisRowCase);
            allControlIntensities.add(thisRowControl);
            numMinFeaturesPerGroup.add((float) Math.min(numCaseFeatures, numControlFeatures));
            runFeatureMapsIndexedByR.add(runSummaryFeatureMap);
        }



        ApplicationContext.infoMessage("Running t-test on " + allCaseIntensities.size() + " out of " +
                allCaseIntensities.size() + " rows...");
        double[][] caseIntensitiesArrayForTTest = new double[allCaseIntensities.size()][allCaseIntensities.get(0).length];
        for (int i=0; i<allCaseIntensities.size(); i++)
            caseIntensitiesArrayForTTest[i] = allCaseIntensities.get(i);
        double[][] controlIntensitiesArrayForTTest = new double[allControlIntensities.size()][allControlIntensities.get(0).length];
        for (int i=0; i<allControlIntensities.size(); i++)
        {
            controlIntensitiesArrayForTTest[i] = allControlIntensities.get(i);
        }

    /**
     * Perform a 2-sided t-test on each row of a couple matrices, one at a time.  Test statistics are necessary in
     * order to determine direction of difference
        */
        Map<String, double[][]> matrixVarMap = new HashMap<String, double[][]>();
        matrixVarMap.put("case", caseIntensitiesArrayForTTest);
        matrixVarMap.put("control", controlIntensitiesArrayForTTest);

        String rResultString = RInterface.evaluateRExpression("pvalues<-c(); tstats<-c();" +
                "for (i in nrow(case):1) { ttestresult= t.test(case[i,],control[i,]); " +
                "pvalues=c(ttestresult$p.value,pvalues); tstats=c(ttestresult$statistic, tstats); }; " +
                "list(resultt=as.vector(tstats), resultp=pvalues, resultq=qvalue(pvalues)$qvalues);\n",
                null, matrixVarMap, new String[] {"qvalue"}, 150000);
        Map<String, String> varStrings =
                RInterface.extractVariableStringsFromListOutput(rResultString);
        List<Float> tScores =   RInterface.parseNumericList(varStrings.get("resultt"));
        List<Float> pValues =   RInterface.parseNumericList(varStrings.get("resultp"));
        List<Float> qValues =   RInterface.parseNumericList(varStrings.get("resultq"));

        ApplicationContext.infoMessage("Done running t-test.");
        if (showCharts)
        {
            PanelWithHistogram pwhTScores = new PanelWithHistogram(tScores,  "t-scores");
            pwhTScores.displayInTab();
            PanelWithHistogram pwhPValues = new PanelWithHistogram(pValues,  "p-values");
            pwhPValues.displayInTab();
            PanelWithHistogram pwhQValues = new PanelWithHistogram(qValues,  "q-values");
            pwhQValues.displayInTab();

            PanelWithBoxAndWhiskerChart pwbawc = new PanelWithBoxAndWhiskerChart("q-values by min per group");
            for (int minRuns=minRunsPerGroup; minRuns<caseRunNames.length; minRuns++)
            {
                //inefficient, but so what?
                List<Float> qvaluesThisCount = new ArrayList<Float>();
                for (int i=0; i<numMinFeaturesPerGroup.size(); i++)
                {
                    if (numMinFeaturesPerGroup.get(i) == minRuns)
                        qvaluesThisCount.add(qValues.get(i));
                }
                if (!qvaluesThisCount.isEmpty())
                    pwbawc.addData(qvaluesThisCount, "" + minRuns);
System.err.println("Min runs represented per group: " + minRuns + ", " + qvaluesThisCount.size() + " rows");
            }
            pwbawc.displayInTab();
        }
/*
        List<Pair<Boolean, Set<String>>> upIndicationsAndPeptides = new ArrayList<Pair<Boolean, Set<String>>>();
        List<Float> charges = new ArrayList<Float>();
        if (outLowQValueArrayFile != null)
        {
            String[] caseControlRunNames = new String[caseRunNames.length + controlRunNames.length];
            System.arraycopy(caseRunNames, 0, caseControlRunNames, 0, caseRunNames.length);
            System.arraycopy(controlRunNames, 0, caseControlRunNames, caseRunNames.length, controlRunNames.length);

            try
            {
                PrintWriter outPW = new PrintWriter(outLowQValueArrayFile);
PrintWriter outPWAll = new PrintWriter(new File(outLowQValueArrayFile.getParentFile(),
                        outLowQValueArrayFile.getName() + ".allq.tsv"));

                StringBuffer headerLineBuf = new StringBuffer("tscore\tqvalue");
                for (String runName : caseControlRunNames)
                {
                    headerLineBuf.append("\t" + runName + "_intensity" + "\t" + runName + "_peptide");
                }
                outPW.println(headerLineBuf);
                outPW.flush();
outPWAll.println(headerLineBuf); outPWAll.flush();

                int numLowQValueRows = 0;
                List<Float> numPeptidesPerRow = new ArrayList<Float>();

                for (int i=0; i<allControlIntensities.size(); i++)
                {
                    float qValue = qValues.get(i);
                    float tScore = tScores.get(i);
                    StringBuffer lineBuf = new StringBuffer(tScore + "\t" + qValue);


                    Set<String> peptidesThisRow = new HashSet<String>();
                    numLowQValueRows++;
                    outPW.print(qValue + "\t" + tScore);
                    int charge = -1;
                    for (String runName : caseControlRunNames)
                    {
                        String featureIntensity = "";
                        String featurePeptide = "";
                        if (runFeatureMap.containsKey(runName))
                        {
                            Feature feature = runFeatureMap.get(runName);
                            charge = feature.getCharge();
                            featureIntensity = "" + feature.getIntensity();
                            featurePeptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                            if (featurePeptide != null)
                            {
                                peptidesThisRow.add(featurePeptide);
                            }
                        }
                        lineBuf.append("\t" + featureIntensity + "\t" + featurePeptide);
                    }
outPWAll.println(lineBuf);outPWAll.flush();
                    if (qValue < maxQValue)
                    {
                        outPW.println(lineBuf);
                        outPW.flush();
                        upIndicationsAndPeptides.add(new Pair<Boolean, Set<String>>(tScore > 0, peptidesThisRow));
                        numPeptidesPerRow.add((float) peptidesThisRow.size());
                        charges.add((float) charge);
                    }
                }
                ApplicationContext.infoMessage("Wrote " + numLowQValueRows + " rows with q-value < " + maxQValue + " to file " +
                        outLowQValueArrayFile.getAbsolutePath());
                outPW.close();
outPWAll.close();
                if (numPeptidesPerRow.size() > 0)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(numPeptidesPerRow, "Peptides per low-q row");
                    pwh.displayInTab();
                }
            }
            catch (FileNotFoundException e)
            {
                throw new RuntimeException("ERROR!! Failed to write low-q-value array file " +
                        outLowQValueArrayFile.getAbsolutePath());
            }
        }

        if (!charges.isEmpty())
        {
            PanelWithHistogram pwh2 = new PanelWithHistogram(charges, "low-q charges");
            pwh2.displayInTab();
        }

        if (outLowQValueAgreeingPeptidePepXMLFile != null)
        {
            try
            {
                List<Feature> agreeingFeatures = new ArrayList<Feature>();
                List<Float> dummyRatios = new ArrayList<Float>();
                for (int i=0; i<upIndicationsAndPeptides.size(); i++)
                {
                    Pair<Boolean, Set<String>> upIndicationAndPeptides = upIndicationsAndPeptides.get(i);
                    int charge = charges.get(i).intValue();
                    boolean up = upIndicationAndPeptides.first;
                    Set<String> allPeptidesThisRow = upIndicationAndPeptides.second;
                    if (allPeptides.size() == 1)
                    {
                        String peptide = allPeptidesThisRow.iterator().next();
                        Protein dummyProtein = new Protein("dummy", peptide.getBytes());

                        Feature dummyFeature = new Feature(agreeingFeatures.size()+1, (float) dummyProtein.getMass(), 200);
                        dummyFeature.setCharge(charge);
                        MS2ExtraInfoDef.addPeptide(dummyFeature, peptide);
                        MS2ExtraInfoDef.setPeptideProphet(dummyFeature, 0.95f);
                        float dummyRatio = 20.0f;
                        if (!up)
                            dummyRatio = 1.0f / dummyRatio;
                        dummyRatios.add(dummyRatio);
                        IsotopicLabelExtraInfoDef.setRatio(dummyFeature, dummyRatio);
                        IsotopicLabelExtraInfoDef.setLabel(dummyFeature, new AnalyzeICAT.IsotopicLabel(57.0201f, 60.0201f, 'C', 1));
                        IsotopicLabelExtraInfoDef.setHeavyIntensity(dummyFeature, 10000f);
                        IsotopicLabelExtraInfoDef.setLightIntensity(dummyFeature, dummyRatio * 10000f);
                        agreeingFeatures.add(dummyFeature);
                    }

                }

                if (dummyRatios.size() > 0)
                {
                    PanelWithHistogram pwhDummyRatios = new PanelWithHistogram(dummyRatios, "Dummy ratios");
                    pwhDummyRatios.displayInTab();
                }
                
                FeatureSet outFeatureSet = new FeatureSet(agreeingFeatures.toArray(new Feature[agreeingFeatures.size()]));

                outFeatureSet.savePepXml(outLowQValueAgreeingPeptidePepXMLFile);
                ApplicationContext.infoMessage("Wrote " + agreeingFeatures.size() + " peptide-agreeing features out of " +
                        upIndicationsAndPeptides.size() + " features with q-value < " + maxQValue +
                        " to file " +
                        outLowQValueAgreeingPeptidePepXMLFile.getAbsolutePath() +
                        "\nDummy 'ratio' is 'case' to 'control': high 'ratio' means higher in case.");
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage(e.getMessage(), e);
                throw new RuntimeException("ERROR!! Failed to write low-q-value array file " +
                        outLowQValueArrayFile.getAbsolutePath());
            }
        }
*/
    }



    public double[] compareIntensities(boolean add1, int minRunsPerGroup)
    {
        double ratioSum = 0;

        double[] result = null;

        List<Pair<Double,Double>> caseControlIntensityPairs =
                new ArrayList<Pair<Double,Double>>();
        List<List<Double>> pairsPlotDataRows = new ArrayList<List<Double>>();

        List<Float> cvs = new ArrayList<Float>();
        List<Float> means = new ArrayList<Float>();
        List<Float> meansWithIds = new ArrayList<Float>();
        List<Float> numSamplesList = new ArrayList<Float>();
        List<String> ids = new ArrayList<String>();

        List<Float> deviationsFromMeanOverMean = new ArrayList<Float>();
        List<Float> deviationsFromMeanOverMeanWithIds = new ArrayList<Float>();
        List<Float> logDeviationsFromMeanOverMeanWithIds = new ArrayList<Float>();
        List<Float> logDeviationsFromMeanOverMean = new ArrayList<Float>();





        List<double[]> allCaseIntensities = new ArrayList<double[]>();
        List<double[]> allControlIntensities = new ArrayList<double[]>();
        List<Integer> origArrayIndices = new ArrayList<Integer>();

        for (int rowIndex=0; rowIndex<rowMaps.length; rowIndex++)
        {
            Map<String,Object> rowMap = rowMaps[rowIndex];
            Set<String> peptideIdsThisRow = new HashSet<String>();

            if (caseRunNames != null && controlRunNames != null)
            {
                List<Double> caseIntensities = new ArrayList<Double>();
                List<Double> controlIntensities = new ArrayList<Double>();
                int numCaseRunValues = 0;
                int numControlRunValues = 0;
                for (String caseRunName : caseRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + caseRunName);
                    if (runIntensity == null)
                        caseIntensities.add(Double.NaN);
                    else
                    {
                        caseIntensities.add(Double.parseDouble(runIntensity.toString()));
                        numCaseRunValues++;
                    }
                }

                for (String controlRunName : controlRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + controlRunName);
                    if (runIntensity == null)
                        controlIntensities.add(Double.NaN);
                    else
                    {
                        controlIntensities.add(Double.parseDouble(runIntensity.toString()));
                        numControlRunValues++;
                    }

                }


                if (caseIntensities.size() == 0 && !add1)
                    continue;

                if (numCaseRunValues < minRunsPerGroup || numControlRunValues < minRunsPerGroup)
                    continue;

                if (caseIntensities.size() == 0 && !add1)
                    continue;

                double meanIntensityCase = BasicStatistics.mean(caseIntensities);
                double meanIntensityControl = BasicStatistics.mean(controlIntensities);
                caseControlIntensityPairs.add(new Pair<Double,Double>(meanIntensityCase, meanIntensityControl));

                //these don't really make any sense
                float cv = (float) BasicStatistics.coefficientOfVariation(new double[] { meanIntensityCase,
                        meanIntensityControl });
                cvs.add(cv);

                ratioSum += meanIntensityCase / meanIntensityControl;

                double[] caseIntensitiesDouble = new double[caseIntensities.size()];
                for (int i=0; i<caseIntensities.size(); i++)
                     caseIntensitiesDouble[i] = caseIntensities.get(i);
                double[] controlIntensitiesDouble = new double[controlIntensities.size()];
                for (int i=0; i<controlIntensities.size(); i++)
                     controlIntensitiesDouble[i] = controlIntensities.get(i);
                allCaseIntensities.add(caseIntensitiesDouble);
                allControlIntensities.add(controlIntensitiesDouble);
                origArrayIndices.add(rowIndex);



            }
            else
            {
                int numSamplesThisRow = 0;
                List<Double> intensities = new ArrayList<Double>();
                for (String key : rowMap.keySet())
                {
                    if (key.startsWith("intensity_"))
                    {
                        String intensityString =
                            rowMap.get(key).toString();
                        if (intensityString != null)
                            intensities.add(Double.parseDouble(intensityString));
                        numSamplesThisRow++;
                    }
                    else if (key.startsWith("peptide_"))
                    {
                        String peptideString =
                            rowMap.get(key).toString();
                        if (peptideString != null)
                        {
                            peptideIdsThisRow.add(peptideString);
                        }
                    }
                }
                boolean hasId = !peptideIdsThisRow.isEmpty();
                if (intensities.size() >= 2)
                {
                    cvs.add((float) BasicStatistics.coefficientOfVariation(intensities));
                    float intensityMean = (float) BasicStatistics.mean(intensities);
                    for (double intensity : intensities)
                    {
                        means.add(intensityMean);
                        deviationsFromMeanOverMean.add((float) (intensity - intensityMean) / intensityMean);
                        logDeviationsFromMeanOverMean.add((float) ((Math.log(intensity) - Math.log(intensityMean)) / Math.log(intensityMean)));

                        if (hasId)
                        {
                            meansWithIds.add(intensityMean);
                            deviationsFromMeanOverMeanWithIds.add((float) (intensity - intensityMean) / intensityMean);
                            logDeviationsFromMeanOverMeanWithIds.add((float) ((Math.log(intensity) - Math.log(intensityMean)) / Math.log(intensityMean)));

                        }
                    }
                    numSamplesList.add((float)numSamplesThisRow);
                    if (peptideIdsThisRow.size() == 1)
                        ids.add(peptideIdsThisRow.iterator().next());
                    else ids.add("");
                }

            }



            List<Double> pairsPlotDataRow = new ArrayList<Double>();
            for (String runName : runNames)
            {
                Object runIntensity = rowMap.get("intensity_" + runName);

                pairsPlotDataRow.add((runIntensity == null) ? Double.NaN : Double.parseDouble(runIntensity.toString()));
            }
            pairsPlotDataRows.add(pairsPlotDataRow);

            
        }
//System.err.println("id\tnumsamples\tcv");
//for (int i=0; i<ids.size(); i++)
//{
//    if (ids.get(i).length() > 0)
//        System.err.println(ids.get(i) + "\t" + numSamplesList.get(i) + "\t" + cvs.get(i));
//}


        ApplicationContext.infoMessage("Coeffs. of Variation: mean: " + BasicStatistics.mean(cvs) + ", median: " + BasicStatistics.median(cvs));

        if (showCharts) {
            PanelWithHistogram pwh = new PanelWithHistogram(cvs, "CVs");
            pwh.displayInTab();
        }



        if (caseRunNames != null && controlRunNames != null && caseRunNames.length > 1 && controlRunNames.length > 1)
        {
            calcTScoresQValues(minRunsPerGroup, true, false);


            int numCommonPeptides = caseControlIntensityPairs.size();

            double[] caseMeanIntensities = new double[numCommonPeptides];
            double[] controlMeanIntensities = new double[numCommonPeptides];
            double[] caseControlIntensityRatios = new double[numCommonPeptides];

            for (int i=0; i<numCommonPeptides; i++)
            {
                caseMeanIntensities[i] = caseControlIntensityPairs.get(i).first;
                controlMeanIntensities[i] = caseControlIntensityPairs.get(i).second;
                caseControlIntensityRatios[i] = caseMeanIntensities[i] / controlMeanIntensities[i];
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


            double meanRatio = ratioSum / caseMeanIntensities.length;


            ApplicationContext.infoMessage("Peptides in common: " + caseMeanIntensities.length);
            PanelWithScatterPlot spd2 = new PanelWithScatterPlot(logControlIntensities, logCaseIntensities, "X is log control, y is log case");
            spd2.displayInTab();
//try {spd2.saveChartToImageFile(new File("/home/dhmay/temp/intensityplot.png"));} catch (Exception e) {ApplicationContext.errorMessage("error:",e);}
            if (showCharts)
            {
                ApplicationContext.infoMessage("Average case-control intensity ratio: " + meanRatio);
                ScatterPlotDialog spd = new ScatterPlotDialog(controlMeanIntensities, caseMeanIntensities, "X is control, y is case");
                spd.setVisible(true);
                spd2.setVisible(true);

                List<Float> ratios = new ArrayList<Float>();
                List<Float> logRatios = new ArrayList<Float>();
                for (int i=0; i<caseMeanIntensities.length; i++)
                {
                    ratios.add((float) (caseMeanIntensities[i] / controlMeanIntensities[i]));
                    logRatios.add((float) Math.log(ratios.get(i)));
                }
                if (showCharts) {
                    new PanelWithHistogram(ratios, "ratios").displayInTab();
                    new PanelWithHistogram(logRatios, "logratios").displayInTab();
                }

            }
            ApplicationContext.infoMessage("correlation coefficient of intensities:" +  BasicStatistics.correlationCoefficient(caseMeanIntensities, controlMeanIntensities));
            ApplicationContext.infoMessage("correlation coefficient of log intensities:" +  BasicStatistics.correlationCoefficient(logCaseIntensities, logControlIntensities));
            result = caseControlIntensityRatios;
        }
        else
        {
            List<Float> logMeans = new ArrayList<Float>();
            for (float mean : means) logMeans.add((float) Math.log(mean));
            List<Float> logMeansWithIds = new ArrayList<Float>();
            for (float mean : meansWithIds) logMeansWithIds.add((float) Math.log(mean));

            if (showCharts) {

                PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
                pwsp.setName("'MA' plot");
                if (!logMeansWithIds.isEmpty())
                    pwsp.addData(logMeansWithIds, logDeviationsFromMeanOverMeanWithIds, "Identified features");
                pwsp.addData(logMeans, logDeviationsFromMeanOverMean, "All features");


                pwsp.setAxisLabels("Mean Intensity (log)", "Deviation / Mean");
                pwsp.displayInTab();
            }

        }

        if (pairsPlotDataRows.get(0).size() > 2) {
            PanelWithRPairsPlot pairsPlot = new PanelWithRPairsPlot();
            for (int i=0; i<pairsPlotDataRows.size(); i++)
            {
                List<Double> pairsPlotDataRow = pairsPlotDataRows.get(i);
                for (int j=0; j<pairsPlotDataRow.size(); j++)
                    if (!Double.isNaN(pairsPlotDataRow.get(j)))
                        pairsPlotDataRow.set(j, Math.log(pairsPlotDataRow.get(j)));
            }
            pairsPlot.setName("Run Pairs Log Int");
            pairsPlot.setChartHeight(900);
            pairsPlot.setChartWidth(900);
            ApplicationContext.infoMessage("Building run pairs plot...");
            pairsPlot.plot(pairsPlotDataRows);
            ApplicationContext.infoMessage("\tDone.");
            pairsPlot.displayInTab();
        }

        if (pairsPlotDataRows.get(0).size() == 2)
        {
            List<Float> ratios = new ArrayList<Float>();
            List<Float> logRatios = new ArrayList<Float>();

            List<Float> logIntensities1 = new ArrayList<Float>();
            List<Float> logIntensities2 = new ArrayList<Float>();

            for (int i=0; i<pairsPlotDataRows.size(); i++)
            {
                ratios.add((float) (pairsPlotDataRows.get(i).get(0) / pairsPlotDataRows.get(i).get(1)));
                logRatios.add((float) Math.log(ratios.get(i)));

                float logIntensity1 = (float) Math.log((pairsPlotDataRows.get(i).get(0)));
                float logIntensity2 = (float) Math.log((pairsPlotDataRows.get(i).get(1)));

                if (!Float.isInfinite(logIntensity1) && !Float.isInfinite(logIntensity2) &&
                        !Float.isNaN(logIntensity1) && !Float.isNaN(logIntensity2)) {
                    logIntensities1.add(logIntensity1);
                    logIntensities2.add(logIntensity2);
                }
            }
            if (showCharts) {
                new PanelWithHistogram(ratios, "ratios").displayInTab();
                new PanelWithHistogram(logRatios, "logratios").displayInTab();
                if (!logIntensities1.isEmpty()) {
                    new PanelWithScatterPlot(logIntensities1, logIntensities2, "logintensities", "logint 1", "logint2").displayInTab();
                    System.err.println("Log intensity correlation: " + BasicStatistics.correlationCoefficient(logIntensities1, logIntensities2));
                }
            }

        }

//System.err.println("Mean Marty Stat: " + BasicStatistics.mean(martyStat));
//System.err.println("Median Marty Stat: " + BasicStatistics.median(martyStat));


        return  result;
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

                int runIndex = 0;
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


    public double[] loadRunIndexes1Based(String[] runNamesToLoad) throws IllegalArgumentException
    {
        List<String> runNameList = new ArrayList<String>(runNamesToLoad.length);
        for (String runName : runNamesToLoad)
            runNameList.add(runName);
        return loadRunIndexes1Based(runNameList);
    }

    /**
     * Load the indexes of the runs with the specified names.  Side effect: outputs messages with missing run names
     * @param runNameList
     * @return
     * @throws IllegalArgumentException  if runs not found
     */
    public double[] loadRunIndexes1Based(List<String> runNameList) throws IllegalArgumentException
    {
        double[] runIndexes1Based = new double[runNameList.size()];

        Set<String> foundRuns = new HashSet<String>();

        for (int i=0; i<runNames.size(); i++)
        {
            if (runNameList.contains(runNames.get(i)))
            {
                runIndexes1Based[runNameList.indexOf(runNames.get(i))] = ((double) (i+1));
                foundRuns.add(runNames.get(i));
            }
        }
        if (foundRuns.size() != runNameList.size())
        {
            if (foundRuns.size() < runNameList.size())
            {
                ApplicationContext.setMessage("Missing Run: ");
                for (String run : runNameList)
                {
                    if (!foundRuns.contains(run))
                        ApplicationContext.infoMessage("\t" + run);
                }
            }
            throw new IllegalArgumentException("Not all runs found in array file.");
        }
        return runIndexes1Based;
    }

    /**
     * Loads case and control run name lists from two different files
     * @param caseRunListFile
     * @param controlRunListFile
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void loadCaseControlRunListsFromFiles(File caseRunListFile, File controlRunListFile)
            throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(caseRunListFile));
        List<String> caseRunNameList = new ArrayList<String>();
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
        List<String> controlRunNameList = new ArrayList<String>();
        while (true)
        {
            String line = br.readLine();
            if (null == line)
                break;
            if (line.length() == 0 || line.charAt(0) == '#')
                continue;
            controlRunNameList.add(line);
        }

        setCaseRunNames(caseRunNameList.toArray(new String[caseRunNameList.size()]));
        setControlRunNames(controlRunNameList.toArray(new String[controlRunNameList.size()]));
    }


    /**
     * Runs a 2-group t-test between the case and control runs on every row in the array with at least
     * minPerGroup run values present in case and control runs.  Calculates t-scores, p-values and q-values
     * @param minPerGroup 
     * @return
     * @throws IOException
     */
    public float[][] runTwoGroupTTestAllRows(int minPerGroup)
            throws IOException
    {
        Map<String, Object> scalarVariableValues = new HashMap<String, Object>();
        scalarVariableValues.put("file", "'" + arrayFile.getAbsolutePath() + "'");
        scalarVariableValues.put("minPerGroup", minPerGroup);

        Map<String, double[]> vectorVariableValues = new HashMap<String, double[]>();
        vectorVariableValues.put("cases", caseRunIndexes1Based);
        vectorVariableValues.put("controls", controlRunIndexes1Based);

        String[] dependentPackageNames = new String[] {"qvalue"};
        StringBuffer rCommandBuf =
                new StringBuffer(RInterface.readResourceFile("/org/fhcrc/cpl/viewer/align/ttest_array.R"));
        rCommandBuf.append("options(max.print = 9999999);\n");
        rCommandBuf.append("ttest_array(file, minPerGroup, cases, controls);\n");

        String rResult = RInterface.evaluateRExpression(rCommandBuf.toString(), scalarVariableValues,
                vectorVariableValues, null, dependentPackageNames, 60000);
        float[][] resultMatrix = RInterface.processRMatrixResponse(rResult);
        return resultMatrix;
    }

    public Map<String, Object>[] getRowMaps()
    {
        return rowMaps;
    }

    /**
     * If the array has been messed with, then id will not be equivalent to line number
     * @param id
     * @return
     */
    public Map<String, Object> getRowMapWithID(int id)
    {
        for (Map<String, Object> rowMap : rowMaps)
            if ((Integer) rowMap.get("id") == id)
                return rowMap;
        return null;
    }

    public Map<String, Object> getRowMap(int i)
    {
        return rowMaps[i];
    }

    public List<String> getRunNames()
    {
        return runNames;
    }

    public int getRunCount()
    {
        return runNames.size();
    }

    public String[] getCaseRunNames()
    {
        return caseRunNames;
    }

    public String[] getControlRunNames()
    {
        return controlRunNames;
    }

    /**
     *
     * @param rowMap
     * @return
     */
    public List<Double> getPresentCaseIntensities(Map<String, Object> rowMap)
    {
        List<Double> result = new ArrayList<Double>();
        for (String runName : caseRunNames)
        {
            Double intensity = getRunIntensity(rowMap, runName);
            if (intensity != null)
                result.add(intensity);
        }
        return result;
    }

    /**
     *
     * @param rowMap
     * @return
     */
    public List<Double> getPresentControlIntensities(Map<String, Object> rowMap)
    {
        List<Double> result = new ArrayList<Double>();
        for (String runName : controlRunNames)
        {
            Double intensity = getRunIntensity(rowMap, runName);
            if (intensity != null)
                result.add(intensity);
        }
        return result;
    }

    /**
     * Return the ratio of the geometric mean case intensity to the mean control intensity,
     * or null if either totally missing
     * @param rowMap
     * @return
     */
    public Double calcRowCaseControlRatioOfGeomMeans(Map<String, Object> rowMap)
    {
        List<Double> caseIntensities = getPresentCaseIntensities(rowMap);
        List<Double> controlIntensities = getPresentControlIntensities(rowMap);

        if (caseIntensities.isEmpty() || controlIntensities.isEmpty())
            return null;
        return BasicStatistics.geometricMean(caseIntensities) / BasicStatistics.geometricMean(controlIntensities);
    }


    /**
     * Side effect: populates indexes of case runs. Must not be run before parsing array
     * @param controlRunNames
     * @throws IllegalArgumentException
     */
    public void setControlRunNames(String[] controlRunNames)  throws IllegalArgumentException
    {
        if (controlRunNames != null)
            controlRunIndexes1Based = loadRunIndexes1Based(controlRunNames);
        this.controlRunNames = controlRunNames;
    }

    /**
     * Side effect: populates indexes of control runs. Must not be run before parsing array
     * @param caseRunNames
     */
    public void setCaseRunNames(String[] caseRunNames)  throws IllegalArgumentException
    {
        if (caseRunNames != null)
            caseRunIndexes1Based = loadRunIndexes1Based(caseRunNames);
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

    public File getOutLowQValueArrayFile()
    {
        return outLowQValueArrayFile;
    }

    public void setOutLowQValueArrayFile(File outLowQValueArrayFile)
    {
        this.outLowQValueArrayFile = outLowQValueArrayFile;
    }

    public float getMaxQValue()
    {
        return maxQValue;
    }

    public void setMaxQValue(float maxQValue)
    {
        this.maxQValue = maxQValue;
    }

    public File getOutLowQValueAgreeingPeptidePepXMLFile()
    {
        return outLowQValueAgreeingPeptidePepXMLFile;
    }

    public void setOutLowQValueAgreeingPeptidePepXMLFile(File outLowQValueAgreeingPeptidePepXMLFile)
    {
        this.outLowQValueAgreeingPeptidePepXMLFile = outLowQValueAgreeingPeptidePepXMLFile;
    }

    public File getDetailsFile()
    {
        return detailsFile;
    }

    public void setDetailsFile(File detailsFile)
    {
        this.detailsFile = detailsFile;
    }

    public boolean doesArrayHaveMs2()
    {
        return arrayHasMs2;
    }

    public void setArrayHasMs2(boolean arrayHasMs2)
    {
        this.arrayHasMs2 = arrayHasMs2;
    }
}
