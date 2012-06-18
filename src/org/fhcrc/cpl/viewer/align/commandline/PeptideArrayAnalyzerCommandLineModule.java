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
import org.fhcrc.cpl.viewer.align.BucketedPeptideArray;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule for feature finding
 */
public class PeptideArrayAnalyzerCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptideArrayAnalyzerCommandLineModule.class);

    protected File file;
    protected File outFile;
    protected File outDir;
    protected File detailsFile;
    boolean allowHalfMatched=false;
    String[] caseRunNames;
    String[] controlRunNames;
    protected double minSignificantRatio = 3;
    protected boolean shouldRequireSamePeptide = false;




    protected boolean caseControlRunsSpecified = false;

    protected int minRunsPerGroup = 1;

    protected int minRunsForConsensusFeature = 2;

    PeptideArrayAnalyzer peptideArrayAnalyzer = null;

    protected boolean showCharts = false;

    Object[] arrayRows;



    protected int minPeptideSupport = 1;
    protected int minFeatureSupport = 1;

    protected static final int MODE_GET_INFO = 0;
    protected static final int MODE_CREATE_FEATURE_FILES_ALL_MATCHED = 1;
    protected static final int MODE_COUNT_MULTIPLY_OBSERVED_PEPTIDEES = 2;
    protected static final int MODE_COMPARE_INTENSITIES_SAME_PEPTIDE = 3;
    protected static final int MODE_CREATE_CONSENSUS_FEATURE_FILE = 4;
    protected static final int MODE_COMPARE_ALL_INTENSITIES = 5;
    protected static final int MODE_COMPARE_ALL_INTENSITIES_ADD_1 = 6;
    protected static final int MODE_COMPARE_INTENSITIES_SAME_PEPTIDE_ADD_1 = 7;
    protected static final int MODE_COMPARE_NON_PEPTIDE_INTENSITIES = 8;
    protected static final int MODE_CREATE_PEPTIDE_RATIO_FILE = 9;
    protected static final int MODE_CALC_CV = 10;




    protected final static String[] modeStrings =
            {
                    "getinfo",
                    "createfeaturefilesallmatched",
                    "countmultiplyobservedpeptides",
                    "comparepeptideintensities",
                    "createconsensusfeaturefile",
                    "compareallintensities",
                    "compareallintensitiesadd1",
                    "comparepeptideintensitiesadd1",
                    "comparenonpeptideintensities"   ,
                    "createpeptideratiofile",
                    "calccv",
            };

    public static final String[] modeExplanations =
            {
                    "Get basic peptide array information",
                    "Create feature files, one per run, containing only those features matched across all runs",
                    "Count the distinct peptides observed in multiple runs",
                    "Compare peptide intensities matched across runs",
                    "Create a 'consensus' feature file, containing features (with details from the first run) that align across at least minconsensusfeatureruns runs",
                    "Compare intensity in the case runs vs. intensity in the control runs",
                    "Compare intensity in the case runs vs. intensity in the control runs, adding 1, so that logs can be used",
                    "Compare intensity in the case runs vs. the control runs for features mapped to the same peptide",
                    "Compare intensity in the case runs vs. the control runs for features mapped to the same peptide, adding 1, so that logs can be used",
                    "Compare intensity in the case runs vs. the control runs for features that do NOT have a peptide assignment",
                    "Create a file with two columns: peptide, and ratio from one array row for the geometric mean peptide intensity of case to control.  Multiple rows per peptide possible",
                    "Calculate CVs for every row.  Give summaries separately for each number of runs in which feature occurs",
            };

    protected int mode=-1;

    public PeptideArrayAnalyzerCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "analyzepeparray";

        mShortDescription = "Tools for analyzing peptide arrays.";
        mHelpMessage = "Tools for analyzing peptide arrays, for comparing MS2 results in one set of runs to "
                + "another and for summarizing overlap of MS1 features found in different runs.";

        CommandLineArgumentDefinition[] argDefs =
               {
                    new EnumeratedValuesArgumentDefinition("mode",true,modeStrings, modeExplanations),
                    createUnnamedFileArgumentDefinition(true, null),
                    new FileToWriteArgumentDefinition("out",false, "output file"),
                    new FileToWriteArgumentDefinition("outdir",false, "output directory"),

                    new BooleanArgumentDefinition("allowhalfmatched",false,"When determining whether peptide matches are made, should it be considered a match when one run has an ID and another run has an intensity but no ID?",
                            allowHalfMatched),
                    new FileToReadArgumentDefinition("caserunlistfile",false,"File containing the names of runs in the case group, one per line"),
                    new FileToReadArgumentDefinition("controlrunlistfile",false,"File containing the names of runs in the control group, one per line"),
                    new StringArgumentDefinition("caserun", false, "Case run"),
                       new StringArgumentDefinition("controlrun", false, "Case run"),

                    new DecimalArgumentDefinition("minsignificantratio",false,"Minimum ratio of intensities considered interesting",
                            minSignificantRatio),
                    new IntegerArgumentDefinition("minconsensusfeatureruns",false,
                            "Minimum number of runs required for a feature to be included in the consensus feature set",
                            minRunsForConsensusFeature),
                    new IntegerArgumentDefinition("minpeptidesupport",false,
                            "Minimum number of runs for which the same peptide was identified",
                            minPeptideSupport),
                    new IntegerArgumentDefinition("minfeaturesupport",false,
                            "Minimum number of runs for which a non-peptide-conflicting feature was identified",
                            minFeatureSupport),
                       new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                       new IntegerArgumentDefinition("minrunspergroup", false,
                               "Minimum number of runs in each group in which a feature must be located to be counted",1),
                       new DecimalArgumentDefinition("maxqvalue", false,
                               "Maximum q-value to keep (differential peptide intensities)", PeptideArrayAnalyzer.MAX_Q_VALUE),
                       new FileToWriteArgumentDefinition("outlowqvaluearrayfile", false,
                               "Output peptide array with low q-values only"),
                       new FileToWriteArgumentDefinition("outqvaluepepxmlfile", false,
                               "Output pepXML file with low q-values only, agreeing peptides"),
                       new BooleanArgumentDefinition("samepeptide", false, "Require all features in row to have same peptide? (for consensus feature set)", shouldRequireSamePeptide)

               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));
        ApplicationContext.infoMessage("Mode: " + modeExplanations[mode]);
        file = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        minRunsPerGroup = getIntegerArgumentValue("minrunspergroup");
        shouldRequireSamePeptide = getBooleanArgumentValue("samepeptide");

        showCharts = getBooleanArgumentValue("showcharts");

        try
        {
            peptideArrayAnalyzer = new PeptideArrayAnalyzer(file);
            peptideArrayAnalyzer.setShowCharts(showCharts);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }

        allowHalfMatched = getBooleanArgumentValue("allowhalfmatched");
        File caseRunListFile = getFileArgumentValue("caserunlistfile");
        File controlRunListFile = getFileArgumentValue("controlrunlistfile");

        caseControlRunsSpecified = hasArgumentValue("caserun") || hasArgumentValue("controlrunlistfile");

        if (hasArgumentValue("caserun"))
        {
            String caseRun = getStringArgumentValue("caserun");
            String controlRun = getStringArgumentValue("controlrun");
            assertArgumentAbsent("caserunlistfile");
            assertArgumentAbsent("controlrunlistfile");
            assertArgumentPresent("controlrun");

            caseRunNames = new String[] {caseRun};
            controlRunNames = new String[] {controlRun};

        }


        if (caseRunListFile != null)
        {
            assertArgumentAbsent("caserun");
            assertArgumentAbsent("controlrun");
            assertArgumentPresent("controlrunlistfile");


            try
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
                caseRunNames = caseRunNameList.toArray(new String[0]);

                if (controlRunListFile != null)
                {
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
                    controlRunNames = controlRunNameList.toArray(new String[0]);
                }
                else
                {
                    List<String> controlRunNameList = new ArrayList<String>();
                    for (String runName : peptideArrayAnalyzer.getRunNames())
                        if (!caseRunNameList.contains(runName))
                            controlRunNameList.add(runName);
                    controlRunNames = controlRunNameList.toArray(new String[controlRunNameList.size()]);
                }

                ApplicationContext.infoMessage("Case runs:");
                for (String caseRun : caseRunNames)
                    ApplicationContext.infoMessage("\t" + caseRun);
                ApplicationContext.infoMessage("Control runs:");
                for (String controlRun : controlRunNames)
                    ApplicationContext.infoMessage("\t" + controlRun);
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }
        else
        {
            List<String> runNames = peptideArrayAnalyzer.getRunNames();
            if (runNames.size() == 2)
            {
                ApplicationContext.setMessage("No case/control run names specified.  Assuming run 1 is control, run 2 is case");
                caseRunNames = new String[1];
                caseRunNames[0] = runNames.get(0);                
                controlRunNames = new String[1];
                controlRunNames[0] = runNames.get(1);
                ApplicationContext.setMessage("Control run: " + controlRunNames[0] + ", Case run: " + caseRunNames[0]);
            }
        }

        peptideArrayAnalyzer.setCaseRunNames(caseRunNames);
        peptideArrayAnalyzer.setControlRunNames(controlRunNames);
        peptideArrayAnalyzer.setMaxQValue(getFloatArgumentValue("maxqvalue"));
        peptideArrayAnalyzer.setOutLowQValueArrayFile(getFileArgumentValue("outlowqvaluearrayfile"));
        peptideArrayAnalyzer.setOutLowQValueAgreeingPeptidePepXMLFile(getFileArgumentValue("outqvaluepepxmlfile"));

        detailsFile = new File(BucketedPeptideArray.calcDetailsFilepath(file.getAbsolutePath()));
        if (detailsFile.exists())
        {
            peptideArrayAnalyzer.setDetailsFile(detailsFile);
            ApplicationContext.infoMessage("Located details file " + detailsFile.getAbsolutePath());
        }
        else
            ApplicationContext.infoMessage("No details file found, looked for: " + detailsFile.getAbsolutePath());

        minSignificantRatio = getDoubleArgumentValue("minsignificantratio");

        minRunsForConsensusFeature = getIntegerArgumentValue("minconsensusfeatureruns");

        minPeptideSupport = getIntegerArgumentValue("minpeptidesupport");
        minFeatureSupport = getIntegerArgumentValue("minfeaturesupport");

        try
        {
            switch (mode)
            {
                case MODE_CREATE_FEATURE_FILES_ALL_MATCHED:
                    assertArgumentPresent("outdir");
                    break;
                case MODE_CREATE_CONSENSUS_FEATURE_FILE:
                    assertArgumentPresent("out");
                    assertArgumentAbsent("outdir");

                    String arrayFilePrefix =
                            file.getName().substring(0, file.getName().indexOf("."));
                    for (String fileName : file.getAbsoluteFile().getParentFile().list())
                    {
                        if (fileName.contains(arrayFilePrefix) &&
                            fileName.endsWith(".details.tsv"))
                        {
                            detailsFile = new File(file.getParentFile(), fileName);                 
                            break;
                        }
                    }
                    if (detailsFile == null)
                        throw new ArgumentValidationException("Unable to find details file for array " +
                                file.getName());
                    break;
                case MODE_GET_INFO:
                    if (caseRunNames != null)
                    {
                        List<String> controlRunNameList =
                                new ArrayList<String>();
                        for (String runName : peptideArrayAnalyzer.getRunNames())
                        {
                            boolean isCase= false;
                            for (String caseRunName : caseRunNames)
                            {
                                if (caseRunName.equals(runName))
                                {
                                    isCase = true;
                                    break;
                                }
                            }
                            if (!isCase)
                                controlRunNameList.add(runName);
                        }
                        controlRunNames = controlRunNameList.toArray(new String[0]);
                    }
                    else
                    {
                        if (peptideArrayAnalyzer.getRunNames().size() == 2)
                        {
                            controlRunNames = new String[] { peptideArrayAnalyzer.getRunNames().get(0) };
                            caseRunNames = new String[] { peptideArrayAnalyzer.getRunNames().get(1) };
                        }
                    }
                    if (caseRunNames != null)
                    {
                        System.err.println("Case runs:");
                        for (String caseRunName : caseRunNames)
                            System.err.println("   " + caseRunName);
                        System.err.println("Control runs:");
                        for (String controlRunName : controlRunNames)
                            System.err.println("   " + controlRunName);
                    }
                    break;
            }
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
        try
        {
            TabLoader tabLoader = new TabLoader(file);
            tabLoader.setReturnElementClass(HashMap.class);
            Object[] rows = tabLoader.load();
            ApplicationContext.setMessage("Array rows: " + rows.length);
            List<String> runNames = new ArrayList<String>();

            for (TabLoader.ColumnDescriptor  column : tabLoader.getColumns())
            {
                _log.debug("loading column " + column.name);
                if (column.name.startsWith("intensity_"))
                {
                    runNames.add(column.name.substring("intensity_".length()));
                    _log.debug("adding run " + runNames.get(runNames.size()-1));
                }
            }
            switch (mode)
            {
                case MODE_GET_INFO:
                    getInfo();
                    break;
                case MODE_CREATE_FEATURE_FILES_ALL_MATCHED:
                    peptideArrayAnalyzer.createFeatureFilesAllMatched(outDir);
                    break;
                case MODE_CREATE_CONSENSUS_FEATURE_FILE:
                    //todo: parameterize the last arg
                    FeatureSet consensusFeatureSet =
                        peptideArrayAnalyzer.createConsensusFeatureSet( detailsFile,
                            minRunsForConsensusFeature,
                                PeptideArrayAnalyzer.CONSENSUS_INTENSITY_MODE_MEAN,shouldRequireSamePeptide, true);
                    ApplicationContext.infoMessage("Created consensus feature set with " +
                        consensusFeatureSet.getFeatures().length + " features");
                    consensusFeatureSet.save(outFile);
                    ApplicationContext.infoMessage("Saved consensus feature set to file " +
                        outFile.getAbsolutePath());
                    break;
                case MODE_COUNT_MULTIPLY_OBSERVED_PEPTIDEES:
                    Map<String, Set<String>> runMatchedPeptides =
                            peptideArrayAnalyzer.getRunMatchedPeptides();

                    Set<String> allMatchedPeptides = new HashSet<String>();
                    Map<String, Integer> peptideMatchCountMap =
                            new HashMap<String, Integer>();
                    //This is horribly inefficient, but correct
                    for (Set<String> runPeptides : runMatchedPeptides.values())
                    {
                        for (String peptideInRun : runPeptides)
                        {
                            if (peptideMatchCountMap.keySet().contains(peptideInRun))
                            {
                                peptideMatchCountMap.put(peptideInRun,
                                        peptideMatchCountMap.get(peptideInRun) + 1);
                            }
                            else
                                peptideMatchCountMap.put(peptideInRun,1);
                        }
                    }
                    Set<String> multiplyMatchedPeptides = new HashSet<String>();
                    Set<String> peptidesMatchedAtLeastThrice = new HashSet<String>();

                    for (String peptide : peptideMatchCountMap.keySet())
                    {
                        allMatchedPeptides.add(peptide);
                        int count = peptideMatchCountMap.get(peptide);
                        if (count>2)
                            multiplyMatchedPeptides.add(peptide);
                        if (count>3)
                            peptidesMatchedAtLeastThrice.add(peptide);
                    }
                    System.err.println("Multiply matched peptides: " + multiplyMatchedPeptides.size());
                    System.err.println("Peptides matched at least 3 times: " + peptidesMatchedAtLeastThrice.size());

                    System.err.println("(out of " + allMatchedPeptides.size() + " total matched)");
                    break;
                case MODE_COMPARE_INTENSITIES_SAME_PEPTIDE:
                    compareIntensitiesSamePeptide(rows, caseRunNames, controlRunNames, false);
                    break;
                case MODE_COMPARE_INTENSITIES_SAME_PEPTIDE_ADD_1:
                    compareIntensitiesSamePeptide(rows, caseRunNames, controlRunNames, true);
                    break;                    
                case MODE_COMPARE_ALL_INTENSITIES:
                    peptideArrayAnalyzer.compareIntensities(false, minRunsPerGroup);
                    break;
                case MODE_COMPARE_ALL_INTENSITIES_ADD_1:
                    peptideArrayAnalyzer.compareIntensities(true, minRunsPerGroup);
                    break;
                case MODE_COMPARE_NON_PEPTIDE_INTENSITIES:
                    compareNonPeptideIntensities();
                    break;
                case MODE_CREATE_PEPTIDE_RATIO_FILE:
                    createPeptideRatioFile();
                    break;
                case MODE_CALC_CV:
                    analyzeCVs();
                    break;
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

    protected void analyzeCVs()
    {
        List<Integer> numbersOfRuns = new ArrayList<Integer>();
        List<Double> cvs = new ArrayList<Double>();
        List<Double> logGeoMeanIntensities = new ArrayList<Double>();

        Map<String, Object>[] rowMaps = peptideArrayAnalyzer.getRowMaps();
        for (Map<String, Object> rowMap : rowMaps)
        {
            List<Double> intensitiesThisRow = new ArrayList<Double>();
            for (String run : peptideArrayAnalyzer.getRunNames())
            {
                Double intensity = peptideArrayAnalyzer.getRunIntensity(rowMap, run);
                if (intensity != null)
                    intensitiesThisRow.add(intensity);
            }
            if (intensitiesThisRow.size() > 1)
            {
                numbersOfRuns.add(intensitiesThisRow.size());
                cvs.add(BasicStatistics.coefficientOfVariation(intensitiesThisRow));
                logGeoMeanIntensities.add(Math.log(BasicStatistics.geometricMean(intensitiesThisRow)));
            }
        }

        for (int i=2; i<=peptideArrayAnalyzer.getRunNames().size(); i++)
        {
            List<Double> cvsThisCount = new ArrayList<Double>();
            List<Double> logGeoMeanIntensitiesThisCount = new ArrayList<Double>();
            for (int j=0; j<cvs.size(); j++)
            {
                if (numbersOfRuns.get(j) == i)
                {
                    cvsThisCount.add(cvs.get(j));
                    logGeoMeanIntensitiesThisCount.add(logGeoMeanIntensities.get(j));
                }

            }
            ApplicationContext.infoMessage("Run count " + i + ": " + cvsThisCount.size() + " rows");
            if (cvsThisCount.size() > 1)
            {
                ApplicationContext.infoMessage("mean CV: " + Rounder.round(BasicStatistics.mean(cvsThisCount),3)  + ", median: " + Rounder.round(BasicStatistics.median(cvsThisCount),3));
                new PanelWithScatterPlot(logGeoMeanIntensitiesThisCount, cvsThisCount, i + "runs").displayInTab();
            }
        }
        ApplicationContext.infoMessage("Overall: " + numbersOfRuns.size() + " rows");
        ApplicationContext.infoMessage("mean CV: " + Rounder.round(BasicStatistics.mean(cvs),3) + ", median: " +
                Rounder.round(BasicStatistics.median(cvs),3));
        new PanelWithScatterPlot(logGeoMeanIntensities, cvs, "overall").displayInTab();


    }

    protected void getInfo()
    {
        if (peptideArrayAnalyzer.doesArrayHaveMs2())
        {
            peptideArrayAnalyzer.analyzeMs2();
            System.err.println("Conflict rows: " + peptideArrayAnalyzer.countConflictRows());
        }
        peptideArrayAnalyzer.analyzeMs1();
        if (peptideArrayAnalyzer.doesArrayHaveMs2())
        {
            if (caseControlRunsSpecified)
                peptideArrayAnalyzer.compareIntensitiesSamePeptide(minSignificantRatio);
        }
    }

    protected void createPeptideRatioFile() throws CommandLineModuleExecutionException
    {
        List<Pair<String, Pair<List<Double>, List<Double>>>> rowPeptidesIntensities = peptideArrayAnalyzer.loadPeptidesCaseControlIntensitiesPerRow();
        PrintWriter outPW = null;
        try
        {
            outPW = new PrintWriter(outFile);
            outPW.println("peptide\tratio");
            for (Pair<String, Pair<List<Double>, List<Double>>> peptideAndIntensityLists : rowPeptidesIntensities)
            {
                String peptide = peptideAndIntensityLists.first;
                List<Double> caseIntensities = peptideAndIntensityLists.second.first;
                List<Double> controlIntensities = peptideAndIntensityLists.second.second;

                double ratio = BasicStatistics.geometricMean(caseIntensities) / BasicStatistics.geometricMean(controlIntensities);
                outPW.println(peptide + "\t" + ratio);
                outPW.flush();
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed writing " + outFile.getAbsolutePath());
        }
        finally
        {
            if (outPW  != null)
                outPW.close();
        }
    }




    /**
     * Compare mean intensities of certain columns against each other in rows
     * in which the same peptide is identified in enough runs.
     *
     * A lot of this is hardcoded for a particular purpose right now.
     * @param rows
     * @throws CommandLineModuleExecutionException
     */

    protected Pair<double[], double[]> compareIntensitiesSamePeptide(Object[] rows,
                                                 String[] caseRunNames,
                                                 String[] controlRunNames,
                                                 boolean add1)
            throws CommandLineModuleExecutionException
    {
        if (caseRunNames == null || controlRunNames == null)
            throw new CommandLineModuleExecutionException("Error: You must define case and control runs");

        Set<String> peptidesHigherInCase =
                new HashSet<String>();
        Set<String> peptidesHigherInControl =
                new HashSet<String>();
        int rowsHigherInControl = 0;
        int rowsHigherInCase = 0;


        double[] logIntensitiesCase = null;
        double[] logIntensitiesControl = null;

            double[] intensitiesCasearray = null;
            double[] intensitiesControlarray = null;

        if (minRunsPerGroup > 1)
        {
            peptideArrayAnalyzer.calcTScoresQValues(minRunsPerGroup, showCharts, true);
        }
        try
        {
            List<Double> intensitiesCase = new ArrayList<Double>();
            List<Double> intensitiesControl = new ArrayList<Double>();

            List<Feature> featuresCase = new ArrayList<Feature>();
            List<Feature> featuresControl = new ArrayList<Feature>();
            int numPeptidesInAgreement = 0;
            int numPeptideConflicts=0;
            for (Object rowObj : rows)
            {
                HashMap rowMap = (HashMap) rowObj;

                int featureSupportCase = 0;
                int peptideSupportCase = 0;
                String peptide = null;
                double intensitySumCase = 0;

                for (String caseRunName : caseRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + caseRunName);
                    if (runIntensity == null && !add1)
                        continue;
                    try
                    {
                        intensitySumCase += Double.parseDouble(runIntensity.toString());
                    }
                    catch (Exception e){}
                    featureSupportCase++;
                    String thisPeptide = (String) rowMap.get("peptide_" + caseRunName);
                    if (thisPeptide == null && !add1)
                        continue;


                    if (add1)
                    {
                        if (peptide == null)
                        {
                            peptide = thisPeptide;
                        }                        
                        if (peptide != null && thisPeptide != null &&
                            !peptide.equals(thisPeptide))
                        {
                                numPeptideConflicts++;
                                peptideSupportCase = 0;
                                break;                            
                        }
                    }
                    else
                    {
                        if (peptide == null)
                        {
                            peptide = thisPeptide;
                            peptideSupportCase++;
                        }
                        else
                        {
                            try
                            {
                                if (peptide.equals(thisPeptide))
                                    peptideSupportCase++;
                                else
                                {
                                    numPeptideConflicts++;
                                    peptideSupportCase =0;
                                    break;
                                }
                            }
                            catch (Exception e) {}
                        }
                    }
                }

                double intensityMeanCase = 0;
                if (featureSupportCase >= minFeatureSupport)
                {
                        intensityMeanCase = intensitySumCase / featureSupportCase;
                }

                double intensitySumControl = 0;
                int featureSupportControl = 0;
                int peptideSupportControl = 0;
                boolean peptideConflict = false;
                for (String controlRunName : controlRunNames)
                {
                    Object runIntensity = rowMap.get("intensity_" + controlRunName);
                    if (runIntensity == null && !add1)
                        continue;
                    try
                    {
                        intensitySumControl += Double.parseDouble(runIntensity.toString());
                    }
                    catch (Exception e){}

                    featureSupportControl++;
                    String thisPeptide = (String) rowMap.get("peptide_" + controlRunName);
                    if (thisPeptide == null && !add1)
                        continue;

                    if (add1)
                    {
                        if (peptide == null)
                        {
                            peptide = thisPeptide;
                        }
                        if (peptide != null && thisPeptide != null &&
                            !peptide.equals(thisPeptide))
                        {
                                numPeptideConflicts++;
                                peptideSupportControl = 0;
                                peptideConflict = true;
                                break;
                        }
                    }
                    else
                    {
                        if (peptide == null)
                        {
                            peptide = thisPeptide;
                            peptideSupportControl++;
                        }
                        else
                        {
                            if (peptide.equals(thisPeptide))
                                peptideSupportControl++;
                            else
                            {
                                numPeptideConflicts++;
                                peptideSupportControl = 0;
                                peptideConflict = true;
                                break;
                            }
                        }
                    }


                }

                double intensityMeanControl = 0;
                if (featureSupportControl >= minFeatureSupport)
                {
                    intensityMeanControl = intensitySumControl / featureSupportControl;
                }

//if (peptide != null) System.err.println(peptide);

                boolean peptideAgreement = false;
                if (peptideSupportControl + peptideSupportCase >= minPeptideSupport &&
                        featureSupportCase >= minFeatureSupport &&
                        featureSupportControl >= minFeatureSupport)
                {
                    numPeptidesInAgreement++;
                    peptideAgreement=true;
                }
                if (add1)
                {
                    intensityMeanControl++;
                    intensityMeanCase++;                     
                }

                if (peptide != null && (peptideAgreement || add1))
                {
//if (!peptideAgreement) System.err.println("no-agree peptide, peptide = " + peptide + ", add1= " + add1);
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
                        peptidesHigherInControl.add(peptide);
                    }

                    intensitiesCase.add(intensityMeanCase);
                    intensitiesControl.add(intensityMeanControl);
                }

                if (!peptideConflict && (peptide != null))
                {
                    if (add1 || peptideAgreement)
                    {
//if (!peptideAgreement) System.err.println("Adding one we wouldn't otherwise, add1 is " + add1 + ", mode is " + mode);
                        Feature controlFeature = new Feature(1,1000,(float) intensityMeanControl);
                        MS2ExtraInfoDef.addPeptide(controlFeature, peptide);
                        featuresControl.add(controlFeature);

                        Feature caseFeature = new Feature(1,1000,(float) intensityMeanCase);
                        MS2ExtraInfoDef.addPeptide(caseFeature, peptide);
                        featuresCase.add(caseFeature);
                    }
                }


            }


            ApplicationContext.infoMessage("# Peptides in agreement: " + numPeptidesInAgreement);
            ApplicationContext.infoMessage("# Peptide conflicts: " + numPeptideConflicts);

            ApplicationContext.infoMessage("# Peptides higher in case: " + peptidesHigherInCase.size() + " peptides in " + rowsHigherInCase + " array rows");
            ApplicationContext.infoMessage("# Peptides higher in control: " + peptidesHigherInControl.size() + " peptides in " + rowsHigherInControl + " array rows");
            ApplicationContext.infoMessage("# Peptides higher in one or the other: " + (peptidesHigherInCase.size() + peptidesHigherInControl.size() +
                               " peptides in " + (rowsHigherInCase + rowsHigherInControl) + " array rows"));


            intensitiesCasearray = new double[featuresCase.size()];
            intensitiesControlarray = new double[featuresCase.size()];

            logIntensitiesCase = new double[featuresCase.size()];
            logIntensitiesControl = new double[featuresCase.size()];

            int numInsideTwofold=0;
            List<Float> cvs = new ArrayList<Float>();
            for (int i=0; i<intensitiesCasearray.length; i++)
            {
                intensitiesCasearray[i] = featuresCase.get(i).getIntensity();
                intensitiesControlarray[i] = featuresControl.get(i).getIntensity();

                logIntensitiesCase[i] = Math.log(intensitiesCasearray[i]);
                logIntensitiesControl[i] = Math.log(intensitiesControlarray[i]);                

                double intensitiesRatio = intensitiesCasearray[i] / intensitiesControlarray[i];

//intensitiesRatio = intensitiesRatio * 5.0/4.0;

                if (intensitiesRatio > 0.5 && intensitiesRatio < 2.0)
                    numInsideTwofold++;


                float cv = (float) BasicStatistics.coefficientOfVariation(new double[] { intensitiesCasearray[i], intensitiesControlarray[i] });
                cvs.add(cv);
            }

            ApplicationContext.infoMessage("Correlation coefficient: " +
                    BasicStatistics.correlationCoefficient(intensitiesCasearray, intensitiesControlarray));

            ApplicationContext.infoMessage("Coeffs. of Variation: mean: " + BasicStatistics.mean(cvs) + ", median: " + BasicStatistics.median(cvs));

            if (showCharts)
            {
                System.err.println("dots on plot: " + intensitiesControlarray.length);
                ScatterPlotDialog spd = new ScatterPlotDialog();
                spd.addData(intensitiesControlarray, intensitiesCasearray, "intensities: x is control, y is case");
                spd.setVisible(true);
                ScatterPlotDialog spd2 = new ScatterPlotDialog();
                spd2.addData(logIntensitiesControl, logIntensitiesCase, "intensities: x is log control, y is log case");

                List<Float> ratios = new ArrayList<Float>();
                List<Float> logRatios = new ArrayList<Float>();
                for (int i=0; i<intensitiesControlarray.length; i++)
                {
                    ratios.add((float) (intensitiesCasearray[i] / intensitiesControlarray[i]));
                    logRatios.add((float) Math.log(ratios.get(i)));
                }
                new PanelWithHistogram(ratios, "ratios").displayInTab();
                new PanelWithHistogram(logRatios, "logratios").displayInTab();




                PanelWithHistogram pwh = new PanelWithHistogram(cvs, "CVs");
                pwh.displayInTab();

////this is special-purpose
//double[] lineX = new double[1300];
//double[] lineY = new double[1300];
//for (int c=1; c<lineX.length; c++)
//{
//    lineX[c] = Math.log(c);
//    lineY[c] = Math.log(4.0 * c / 5.0);
//}
//spd.addData(lineX, lineY,"4:5 line");


                spd2.setVisible(true);
            }
            System.err.println("Same peptide intensity summary:");
System.err.println("Within twofold: " + numInsideTwofold + " out of " + intensitiesCasearray.length);

            if (outFile != null)
            {
                ApplicationContext.infoMessage("Writing ratios to file " + outFile.getAbsolutePath());
                PrintWriter outPW = null;
                try
                {

                    outPW = new PrintWriter(outFile);
                    outPW.println("peptide\tratio");
                    for (int i=0; i<featuresCase.size(); i++)
                    {
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(featuresCase.get(i));
                        outPW.println(peptide + "\t" + (intensitiesCasearray[i] / intensitiesControlarray[i]));
                    }
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
                finally
                {
                    if (outPW != null)
                        outPW.close();
                }
                ApplicationContext.infoMessage("Done writing ratios");
            }

//            if (outDir != null)
//            {
//                FeatureSet caseFeatureSet =
//                        new FeatureSet(featuresCase.toArray(new Feature[featuresCase.size()]));
//                FeatureSet controlFeatureSet =
//                        new FeatureSet(featuresControl.toArray(new Feature[featuresControl.size()]));
//
//                File caseFile = new File(outDir, "case.peptides.tsv");
//                caseFeatureSet.save(caseFile);
//                ApplicationContext.setMessage("Wrote case features to " + caseFile.getAbsolutePath());
//
//                File controlFile = new File(outDir, "control.peptides.tsv");
//                controlFeatureSet.save(controlFile);
//                ApplicationContext.setMessage("Wrote control features to " + controlFile.getAbsolutePath());
//            }

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        return new Pair<double[], double[]>(intensitiesCasearray, intensitiesControlarray);
    }


    public void compareNonPeptideIntensities()
    {
        Map<String,Object>[] rowMaps = peptideArrayAnalyzer.getRowMaps();
        double ratioSum = 0;

        List<Pair<Double,Double>> caseControlIntensityPairs =
                new ArrayList<Pair<Double,Double>>();

        for (Map<String,Object> rowMap : rowMaps)
        {
            Object peptideObject = null;
            double intensitySumCase = 0;
            for (String caseRunName : caseRunNames)
            {
                Object runIntensity = rowMap.get("intensity_" + caseRunName);
                if (runIntensity == null)
                    continue;
                intensitySumCase += Double.parseDouble(runIntensity.toString());
                if (peptideObject == null)
                {
                    peptideObject = rowMap.get("peptide_" + caseRunName);
                }
            }

            if (peptideObject != null)
            {
//System.err.println("Tossing peptide " + peptideObject.toString());
                continue;
            }

            if (intensitySumCase == 0)
                continue;

            double intensitySumControl = 0;
            for (String controlRunName : controlRunNames)
            {
                Object runIntensity = rowMap.get("intensity_" + controlRunName);
                if (runIntensity == null)
                    continue;
                intensitySumControl += Double.parseDouble(runIntensity.toString());
                if (peptideObject == null)
                {
                    peptideObject = rowMap.get("peptide_" + controlRunName);
                }
            }

            if (intensitySumControl == 0)
                continue;
            if (peptideObject != null)
            {
                continue;
            }

            caseControlIntensityPairs.add(new Pair<Double,Double>(intensitySumCase, intensitySumControl));
            ratioSum += intensitySumCase / intensitySumControl;
        }

        double[] caseIntensities = new double[caseControlIntensityPairs.size()];
        double[] controlIntensities = new double[caseControlIntensityPairs.size()];

        for (int i=0; i<caseControlIntensityPairs.size(); i++)
        {
            caseIntensities[i] = caseControlIntensityPairs.get(i).first;
            controlIntensities[i] = caseControlIntensityPairs.get(i).second;
        }

        double[] logCaseIntensities = new double[caseControlIntensityPairs.size()];
        double[] logControlIntensities = new double[caseControlIntensityPairs.size()];

        for (int i=0; i<caseControlIntensityPairs.size(); i++)
        {
            logCaseIntensities[i] = Math.log(caseControlIntensityPairs.get(i).first);
            logControlIntensities[i] = Math.log(caseControlIntensityPairs.get(i).second);
        }


        double meanRatio = ratioSum / caseIntensities.length;
        ApplicationContext.infoMessage("Average case-control intensity ratio: " + meanRatio);
        ScatterPlotDialog spd = new ScatterPlotDialog(controlIntensities, caseIntensities, "X is control, y is case");
        spd.setVisible(true);

        ScatterPlotDialog spd2 = new ScatterPlotDialog(logControlIntensities, logCaseIntensities, "X is log control, y is log case");
        spd2.setVisible(true);
    }

}
