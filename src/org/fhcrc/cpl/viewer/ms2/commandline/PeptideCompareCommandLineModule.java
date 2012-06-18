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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.normalize.Normalizer;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.BrowserController;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.util.*;
import java.util.List;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.awt.*;


/**
 */
public class PeptideCompareCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule                                                     
{
    protected static Logger _log = Logger.getLogger(PeptideCompareCommandLineModule.class);

    protected FeatureSet[] featureSets=null;

    protected List<File> featureFiles1=null;
    protected List<File> featureFiles2=null;


    protected double minPeptideProphet=0;
    protected boolean normalize = true;
    protected boolean includeUnmatched = true;

    protected boolean outFormatPepXML = false;

    public static final float INFINITE_RATIO = 20f;
    public static final float ZERO_RATIO = .05f;

    protected boolean boundPeptideRatios=true;

    protected boolean shouldListAllCommon=false;

    protected List<Float> allSet1Values=new ArrayList<Float>();
    protected List<Float> allSet2Values=new ArrayList<Float>();
    protected List<Float> set1OnlyValues=new ArrayList<Float>();
    protected List<Float> set2OnlyValues=new ArrayList<Float>();

    protected List<Float> allSet1LogValues=new ArrayList<Float>();
    protected List<Float> allSet2LogValues=new ArrayList<Float>();

    protected boolean summaryOnly = false;

    protected boolean withinCharge = false;

    protected String xAxisLabel = null;
    protected String yAxisLabel = null;





    protected final static String[] modeStrings = {
            "showoverlap",
            "createfeaturefileforpeptideunion",
            "plottimes",
            "plotintensities",
            "plottotalintensities",
            "plotpprophet",
            "plotfval",
            "plotkscoreorxcorr",
            "calcratios",
            "plotratios",
            "plotlightarea",
            "idcluster",
            "plotspectralcounts",
            "plotpairwiseratios",
            "plotexpect",

//                                                   "plotpairedspectralcounts",
//                                                   "plotgroupedspectralcounts"
    };

    protected static final String[] modeExplanations =
            {
                    "Show overlap between peptide identifications in two or more files",
                    "Create a feature file containing all the features in the first file whose peptides are found in all other files, as well",
                    "Plot time in one set vs. time in the other set, by peptide",
                    "Plot intensity in one set vs. intensity in the other set, by peptide",
                    "Plot total intensity in one set vs. total intensity in the other set, by peptide",
                    "Plot PeptideProphet probability vs. PeptideProphet probability in the other set, by peptide",
                    "Plot PeptideProphet fval vs. PeptideProphet fval in the other set, by peptide",                    
                    "Plot kscore or xcorr score in one set (which is available) against kscore or xcorr score in the other set, by peptide",
                    "Calculate ratios (run 1 : run 2) for each peptide",
                    "Plot peptide ratios against each other",
                    "Plot light areas (for peptides with ratios)",
                    "Cluster runs by peptide identifications",
                    "Plot spectral counts",
                    "Plot pairwise ratios (for multiple sets)",
                    "Plot expect value, by peptide",
            };

    protected static final int MODE_SHOW_OVERLAP = 0;
    protected static final int MODE_CREATE_FEATURES_FOR_UNION = 1;
    protected static final int MODE_PLOT_TIMES = 2;
    protected static final int MODE_PLOT_INTENSITIES = 3;
    protected static final int MODE_PLOT_TOTAL_INTENSITIES = 4;
    protected static final int MODE_PLOT_PEPTIDEPROPHET = 5;
    protected static final int MODE_PLOT_FVAL= 6;
    protected static final int MODE_PLOT_KSCORE_OR_XCORR= 7;
    protected static final int MODE_CALCULATE_RATIOS=8;
    protected static final int MODE_PLOT_RATIOS=9;
    protected static final int MODE_PLOT_LIGHT_AREAS=10;
    protected static final int MODE_ID_CLUSTER=11;
    protected static final int MODE_PLOT_SPECTRAL_COUNTS=12;
    protected static final int MODE_PLOT_PAIRWISE_RATIOS=13;
    protected static final int MODE_PLOT_EXPECT=14;







    protected boolean showCharts = false;

//    protected static final int MODE_PLOT_PAIRED_SPECTRAL_COUNTS = 4;
//    protected static final int MODE_PLOT_GROUPED_SPECTRAL_COUNTS = 5;


    protected int mode=-1;

    protected File outFile = null;

    protected File fastaFile = null;


    public PeptideCompareCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "peptidecompare";
        mShortDescription = "Compare peptide IDs between multiple runs";
        mHelpMessage = "Perform a comparison of peptides identified between multiple runs (.tsv or .pep.xml files).  You can analyze the " +
                "overlap between peptide IDs, plot various characteristics of the identified peptides against each " +
                "other (for 2 runs only), etc.  The 'minpprophet' argument determines the minimum PeptideProphet value " +
                "for all peptides to be compared.";
        mUsageMessage = CommandLineModule.MODULE_USAGE_AUTOMATIC;

        addArgumentDefinition(new EnumeratedValuesArgumentDefinition("mode",true,modeStrings,
                modeExplanations));
        addArgumentDefinition(createUnnamedSeriesFeatureFileArgumentDefinition(
                false, "input files"));
        addArgumentDefinition(new DecimalArgumentDefinition("minpprophet",false,"Minimum PeptideProphet score",
                minPeptideProphet));
        addArgumentDefinition(new FileToWriteArgumentDefinition("out",false,"Output file"));
        addArgumentDefinition(new BooleanArgumentDefinition("showcharts", false, "show charts",showCharts));
        addArgumentDefinition(new BooleanArgumentDefinition("normalize", false, "normalize intensities (for ratio mode)", normalize));
        addArgumentDefinition(new BooleanArgumentDefinition("includeunmatched", false, "include unmatched peptides (for ratio mode)", includeUnmatched));
        addArgumentDefinition(new BooleanArgumentDefinition("outformatpepxml", false, "output to pepxml", outFormatPepXML));
        addArgumentDefinition(new FileToReadArgumentDefinition("fasta", false, "fasta file"));
        addArgumentDefinition(new DirectoryToReadArgumentDefinition("featuresdir1", false,
                "First directory full of feature files"));
        addArgumentDefinition(new DirectoryToReadArgumentDefinition("featuresdir2", false,
                "Second directory full of feature files"));
        addArgumentDefinition(new BooleanArgumentDefinition("summaryonly", false,
                "Summary-level info only? For directories", summaryOnly));
        addArgumentDefinition(new BooleanArgumentDefinition("withincharge", false, "Only compare peptides within charge states", withinCharge));
        addArgumentDefinition(new StringArgumentDefinition("xaxislabel",false, "label for X axis"));
                addArgumentDefinition(new StringArgumentDefinition("yaxislabel",false, "label for Y axis"));
        addArgumentDefinition(new BooleanArgumentDefinition("listallcommon", false, "List all peptides common to all runs?",
                shouldListAllCommon));

    }



    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        _log.debug("Mode: " + modeStrings[mode]);
        outFile = getFileArgumentValue("out");


        shouldListAllCommon = getBooleanArgumentValue("listallcommon");

        minPeptideProphet = getDoubleArgumentValue("minpprophet");
        xAxisLabel = getStringArgumentValue("xaxislabel");
        yAxisLabel = getStringArgumentValue("yaxislabel");

        if (xAxisLabel != null)
            xAxisLabel = xAxisLabel.replace("_", " ");
        if (yAxisLabel != null)
            yAxisLabel = yAxisLabel.replace("_", " ");


        if (this.hasUnnamedSeriesArgumentValue())
        {
            Object[] unnamedSeriesObjects = getUnnamedSeriesArgumentValues();
            featureSets = new FeatureSet[unnamedSeriesObjects.length];
            for (int i=0; i<unnamedSeriesObjects.length; i++)
                featureSets[i] = (FeatureSet) unnamedSeriesObjects[i];
        }
        else
        {
            featureSets = null;
            assertArgumentPresent("featuresdir1");
            assertArgumentPresent("featuresdir2");

            File featureFilesDir1 = getFileArgumentValue("featuresdir1");
            File featureFilesDir2 = getFileArgumentValue("featuresdir2");

            featureFiles1 = new ArrayList<File>();
            featureFiles2 = new ArrayList<File>();


            for (File file1 : featureFilesDir1.listFiles())
            {
                if (file1.isDirectory())
                    continue;
                File file2 = null;
                try
                {
                    file2 = CommandLineModuleUtilities.findFileLikeFile(file1, featureFilesDir2, "");
                }
                catch (FileNotFoundException e)
                {
                    ApplicationContext.infoMessage("Warning: unable to find match for file " + file1.getName());
                    continue;
                }
                featureFiles1.add(file1);
                featureFiles2.add(file2);
            }
            if (featureFiles1.size() == 0)
                throw new ArgumentValidationException("No pairs of files identified to process");

        }


        switch (mode)
        {
            case MODE_CREATE_FEATURES_FOR_UNION:
                assertArgumentPresent("out");
                break;
            case MODE_CALCULATE_RATIOS:
                if (featureSets.length != 2)
                    throw new ArgumentValidationException("You must provide exactly two files for this mode");
                break;
            case MODE_ID_CLUSTER:
                    if (featureSets.length < 3)
                        throw new ArgumentValidationException("Need at least 3 feature sets to cluster");
                break;
        }
        showCharts = getBooleanArgumentValue("showcharts");
        includeUnmatched = getBooleanArgumentValue("includeunmatched");
        normalize = getBooleanArgumentValue("normalize");
        outFormatPepXML = getBooleanArgumentValue("outFormatPepXML");
        fastaFile = getFileArgumentValue("fasta");

        summaryOnly = getBooleanArgumentValue("summaryonly");

        withinCharge = getBooleanArgumentValue("withincharge");
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (mode == MODE_PLOT_TIMES||
                mode == MODE_PLOT_INTENSITIES||
                mode == MODE_PLOT_TOTAL_INTENSITIES||
                mode == MODE_PLOT_PEPTIDEPROPHET||
                mode == MODE_PLOT_FVAL||
                mode == MODE_PLOT_KSCORE_OR_XCORR||
                mode == MODE_PLOT_RATIOS || mode == MODE_PLOT_LIGHT_AREAS ||
                mode == MODE_ID_CLUSTER || mode == MODE_PLOT_SPECTRAL_COUNTS ||
                mode == MODE_PLOT_PAIRWISE_RATIOS |
                mode == MODE_PLOT_EXPECT)
            showCharts=true;
        if (featureSets != null)
        {
            handleFiles(featureSets);
        }
        else
        {
            if (summaryOnly)
                showCharts = false;

            boolean oldShowCharts = showCharts;
            showCharts = true;

            Set<String> peptidesDir1 = new HashSet<String>();
            Set<String> peptidesDir2 = new HashSet<String>();

            Set<String> quantPeptidesDir1 = new HashSet<String>();
            Set<String> quantPeptidesDir2 = new HashSet<String>();

            FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
                    sel.setMinPProphet((float)minPeptideProphet);

            ApplicationContext.setMessage("Counting peptides across ALL files");
            for (File featureFile : featureFiles1)
            {
                try
                {
                    FeatureSet featureSet = new FeatureSet(featureFile);
                    featureSet = featureSet.filter(sel);
                    
                    for (Feature feature : featureSet.getFeatures())
                    {
                        List<String> peptideList = MS2ExtraInfoDef.getPeptideList(feature);
                        if (peptideList != null)
                        {
                            for (String peptide : peptideList)
                            {
                                peptidesDir1.add(peptide);
                                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                                    quantPeptidesDir1.add(peptide);
                            }
                        }
                    }
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
            }



            for (File featureFile : featureFiles2)
            {
                try
                {
                    FeatureSet featureSet = new FeatureSet(featureFile);
                    featureSet = featureSet.filter(sel);
                    for (Feature feature : featureSet.getFeatures())
                    {
                        List<String> peptideList = MS2ExtraInfoDef.getPeptideList(feature);
                        if (peptideList != null)
                        {
                            for (String peptide : peptideList)
                            {
                                peptidesDir2.add(peptide);
                                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                                    quantPeptidesDir2.add(peptide);
                            }
                        }
                    }
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
            }

            Set<String> commonPeptides = new HashSet<String>();

            for (String peptide1 : peptidesDir1)
                if (peptidesDir2.contains(peptide1))
                    commonPeptides.add(peptide1);

            Set<String> commonQuantPeptides = new HashSet<String>();

            for (String peptide1 : quantPeptidesDir1)
                if (quantPeptidesDir2.contains(peptide1))
                    commonQuantPeptides.add(peptide1);





            showCharts = oldShowCharts;

            if (!summaryOnly)
            {
                Set<String> pepsIn2sButNot1s = new HashSet<String>();
                Set<String> quantPepsIn2sButNot1s = new HashSet<String>();

                for (int i=0; i<featureFiles1.size(); i++)
                {
                    try
                    {
                        FeatureSet featureSet1 = new FeatureSet(featureFiles1.get(i));
                        FeatureSet featureSet2 = new FeatureSet(featureFiles2.get(i));

                        ApplicationContext.infoMessage("Processing pair (" + featureFiles1.get(i).getName() +
                                ", " + featureFiles2.get(i));

                        Pair<Set<String>,Set<String>> bothPepsIn2ButNot1 =
                                handleFiles(new FeatureSet[]{ featureSet1, featureSet2 });
                        if (bothPepsIn2ButNot1 != null)
                        {
                            Set<String> pepsIn2ButNot1 = bothPepsIn2ButNot1.first;
                            Set<String> quantPepsIn2ButNot1 = bothPepsIn2ButNot1.second;


                            if (pepsIn2ButNot1 != null)
                                pepsIn2sButNot1s.addAll(pepsIn2ButNot1);
                            if (quantPepsIn2ButNot1 != null)
                                quantPepsIn2sButNot1s.addAll(quantPepsIn2ButNot1);
                        }



                    }
                    catch (Exception e)
                    {
                        ApplicationContext.infoMessage("WARNING: failed to load features from either file" +
                                featureFiles1.get(i).getAbsolutePath() + " or file " + featureFiles2.get(i).getAbsolutePath() +
                                ".  Message: " + e.getMessage());
                        e.printStackTrace(System.err);
                    }
                }
                if (!pepsIn2sButNot1s.isEmpty())
                {
                    pepsIn2sButNot1s.retainAll(commonPeptides);
                    ApplicationContext.infoMessage("Total peptides ever seen in a run from dir 2 and not the corresponding dir 1 run, but in some dir 1 run: " +
                            pepsIn2sButNot1s.size());
                }
                if (!quantPepsIn2sButNot1s.isEmpty())
                {
                    quantPepsIn2sButNot1s.retainAll(commonQuantPeptides);

                    ApplicationContext.infoMessage("Total QUANTITATED peptides ever seen in a run from dir 2 and not the corresponding dir 1 run, but in some dir 1 run: " +
                            quantPepsIn2sButNot1s.size());
                }

            }
            if (showCharts &&
                    (mode == MODE_PLOT_TIMES||
                            mode == MODE_PLOT_INTENSITIES||
                            mode == MODE_PLOT_TOTAL_INTENSITIES||
                            mode == MODE_PLOT_PEPTIDEPROPHET||
                            mode == MODE_PLOT_FVAL||
                            mode == MODE_PLOT_KSCORE_OR_XCORR||
                            mode == MODE_PLOT_RATIOS |
                    mode == MODE_PLOT_EXPECT))
            {
                PanelWithScatterPlot allPSP =
                        new PanelWithScatterPlot(allSet1Values, allSet2Values, "All Common Peptide values");
                allPSP.setPointSize(1);

                int numBothHigh=0;
                int num1High2Low=0;
                int num2High1Low=0;
                
                for (int i=0; i<allSet1Values.size(); i++)
                {
                    float val1 = allSet1Values.get(i);
                    float val2 = allSet2Values.get(i);

                    if (val1 >= 0.9f && val2 >= 0.9f)
                        numBothHigh++;
                    else if (val1 >= 0.9f && val2 < 0.9f)
                        num1High2Low++;
                    else if (val2 >= 0.9f && val1 < 0.9f)
                        num2High1Low++;                    
                }

                int num1High2Gone=0;
                for (Float val :set1OnlyValues)
                    if (val >= 0.9f)
                        num1High2Gone++;
                int num2High1Gone=0;
                for (Float val :set2OnlyValues)
                    if (val >= 0.9f)
                        num2High1Gone++;

                int numTotal=numBothHigh+num1High2Low+num2High1Low+num1High2Gone+num2High1Gone;

                ApplicationContext.infoMessage("Total points (on and off chart): " + numTotal);
                ApplicationContext.infoMessage("High in Both: " + numBothHigh + " (" + (numBothHigh*100/numTotal) + "%)");
                ApplicationContext.infoMessage("High in 1, low in 2: " + num1High2Low + " (" + (num1High2Low*100/numTotal) + "%)");
                ApplicationContext.infoMessage("High in 2, low in 1: " + num2High1Low + " (" + (num2High1Low*100/numTotal) + "%)");
                ApplicationContext.infoMessage("Only in 1, and high: " + num1High2Gone + " (" + (num1High2Gone*100/numTotal) + "%)");
                ApplicationContext.infoMessage("Only in 2, and high: " + num2High1Gone + " (" + (num2High1Gone*100/numTotal) + "%)");


                ApplicationContext.infoMessage("Correlation Coefficient::: " +
                    BasicStatistics.correlationCoefficient(allSet1Values, allSet2Values));

                String xAxisLabelToUse = "Sets 1";
                if (xAxisLabel != null)
                    xAxisLabelToUse = xAxisLabel;
                String yAxisLabelToUse = "Sets 2";
                if (yAxisLabel != null)
                    yAxisLabelToUse = yAxisLabel;
                allPSP.setAxisLabels(xAxisLabelToUse, yAxisLabelToUse);
                allPSP.displayInTab();                

                if (mode == MODE_PLOT_RATIOS)
                {
                    List<Float> allReasonableSet1Values = new ArrayList<Float>();
                    List<Float> allReasonableSet2Values = new ArrayList<Float>();
                    List<Float> allReasonablePercentChanges = new ArrayList<Float>();
                    int numUnder5PercentChange=0;
                    for (int i=0; i<allSet1Values.size(); i++)
                    {
                        float val1 = allSet1Values.get(i);
                        float val2 = allSet2Values.get(i);

                        if (val1 > 0.0000001 && val1 < 100 && val2 > 0.0000001 && val2 < 100)
                        {
                            allReasonableSet1Values.add(val1);
                            allReasonableSet2Values.add(val2);
                            float percentChange = Math.abs(100 * (val2-val1)/val1);
                            allReasonablePercentChanges.add(percentChange);
                            if (percentChange < 5)
                                numUnder5PercentChange++;
                        }
                    }
                    ApplicationContext.infoMessage("Correlation Coefficient of 'reasonable' values: " +
                            BasicStatistics.correlationCoefficient(allReasonableSet1Values, allReasonableSet2Values));

                    ApplicationContext.infoMessage("95th percentile of % changes: " +
                            BasicStatistics.percentile(allReasonablePercentChanges, 95) + "%");
                    ApplicationContext.infoMessage("% of % changes under 5%: " + (100 * numUnder5PercentChange / allReasonablePercentChanges.size()) + "%");

                    PanelWithScatterPlot reasonablePSP =
                            new PanelWithScatterPlot(allReasonableSet1Values, allReasonableSet2Values, "Reasonable Common Peptide values");
                    reasonablePSP.setPointSize(1);
                    reasonablePSP.displayInTab();


                    List<Float> allReasonableSet1LogValues = new ArrayList<Float>();
                    List<Float> allReasonableSet2LogValues = new ArrayList<Float>();
                    for (int i=0; i<allReasonableSet1Values.size(); i++)
                    {
                        float log1 = (float)Math.log(allReasonableSet1Values.get(i));
                        float log2 = (float)Math.log(allReasonableSet2Values.get(i));

                        if (!Float.isNaN(log1) && !Float.isNaN(log2) && !Float.isInfinite(log1) && !Float.isInfinite(log2))
                        {
                            allReasonableSet1LogValues.add(log1);
                            allReasonableSet2LogValues.add(log2);
                        }
                    }                                                   
                    ApplicationContext.infoMessage("Correlation Coefficient of LOG 'reasonable' values: " +
                            BasicStatistics.correlationCoefficient(allReasonableSet1LogValues, allReasonableSet2LogValues));
                    //fishy
                   allSet1LogValues = allReasonableSet1LogValues;
                   allSet1LogValues = allReasonableSet2LogValues;
                }


                ApplicationContext.infoMessage("Values on all-values plot: " + allSet1Values.size());

                if (mode == MODE_PLOT_INTENSITIES || mode == MODE_PLOT_TOTAL_INTENSITIES ||
                        mode == MODE_PLOT_RATIOS)
                {
                    PanelWithScatterPlot allLogPSP =
                        new PanelWithScatterPlot(allSet1LogValues, allSet2LogValues, "All Common Peptide values (log)");
                    allLogPSP.setPointSize(1);

                    xAxisLabelToUse = "Sets 1 Log";
                    if (xAxisLabel != null)
                        xAxisLabelToUse = xAxisLabel;
                    yAxisLabelToUse = "Sets 2 Log";
                    if (yAxisLabel != null)
                    yAxisLabelToUse = yAxisLabel;
                    allLogPSP.setAxisLabels(xAxisLabelToUse, yAxisLabelToUse);
                    allLogPSP.displayInTab();
                    ApplicationContext.infoMessage("Values on all-logvalues plot: " + allSet1LogValues.size());

                }
            }            

            System.err.println("Summary (" + featureFiles1.size() + " files per dir) Peptide Numbers:");
            System.err.println("\tDirectory 1: " + peptidesDir1.size() + ".  Unique: " + (peptidesDir1.size() - commonPeptides.size()));
            System.err.println("\tDirectory 2: " + peptidesDir2.size() + ".  Unique: " + (peptidesDir2.size() - commonPeptides.size()));
            System.err.println("\tCommon: " + commonPeptides.size());

            System.err.println("Summary (" + featureFiles1.size() + " files per dir) Quantified Peptide Numbers:");
            System.err.println("\tDirectory 1: " + quantPeptidesDir1.size() + ".  Unique: " + (quantPeptidesDir1.size() - commonQuantPeptides.size()));
            System.err.println("\tDirectory 2: " + quantPeptidesDir2.size() + ".  Unique: " + (quantPeptidesDir2.size() - commonQuantPeptides.size()));
            System.err.println("\tCommon: " + commonQuantPeptides.size());

            

        }
    }

    /**
     * returns set of peptides in 2 but not 1
     * @param featureSetsToHandle
     * @return
     * @throws CommandLineModuleExecutionException
     */
    public Pair<Set<String>,Set<String>> handleFiles(FeatureSet[] featureSetsToHandle)
            throws CommandLineModuleExecutionException
    {
        Pair<Set<String>,Set<String>> result = null;
        switch (mode)
        {
            case MODE_SHOW_OVERLAP:
                result = showOverlap(featureSetsToHandle);
                break;
            case MODE_CREATE_FEATURES_FOR_UNION:
                createFeatureFileForPeptideUnion(featureSetsToHandle);
                break;
            case MODE_PLOT_TIMES:
            case MODE_PLOT_INTENSITIES:
            case MODE_PLOT_TOTAL_INTENSITIES:                
            case MODE_PLOT_PEPTIDEPROPHET:
            case MODE_PLOT_EXPECT:
            case MODE_PLOT_FVAL:
            case MODE_PLOT_KSCORE_OR_XCORR:
            case MODE_PLOT_RATIOS:
            case MODE_PLOT_LIGHT_AREAS:
                plotSomeAttribute(featureSetsToHandle);
                break;
            case MODE_CALCULATE_RATIOS:
                calculateRatios(featureSetsToHandle);
                break;
//            case MODE_PLOT_PAIRED_SPECTRAL_COUNTS:
//                plotPairedSpectralCounts();
//                break;
//            case MODE_PLOT_GROUPED_SPECTRAL_COUNTS:
//                plotGroupedSpectralCounts();
//                break;
            case MODE_ID_CLUSTER:
                double[][] idOverlaps = calcOverlaps(featureSetsToHandle);
                double[][] dissimilarities = new double[idOverlaps.length][idOverlaps[0].length];
                for (int i=0; i<idOverlaps.length; i++)
                {
                    for (int j=0; j<idOverlaps[0].length; j++)
                    {
                        double overlap = idOverlaps[i][j];
                        double uniqueCount = idOverlaps[i][i] + idOverlaps[j][j] - overlap;
                        dissimilarities[i][j] = (uniqueCount - overlap) / uniqueCount;
//    System.err.println(i + ", " + j + ": " + overlap + ", " + uniqueCount + ", " + dissimilarities[i][j]);

                    }
                }
                Map<String, double[][]> matrixMap = new HashMap<String, double[][]>();
                matrixMap.put("dissimilarities", dissimilarities);
String headerLine = "";
for (int i=0; i<dissimilarities[0].length; i++)
{
    if (i>0) headerLine = headerLine + "\t";
    String baseName = featureSets[i].getSourceFile().getName();
    if (baseName.contains("."))
        baseName = baseName.substring(0, baseName.indexOf("."));
    headerLine = headerLine + baseName;
}
System.err.println("overlap");
                System.err.println(headerLine);

for (int i=0; i<idOverlaps.length; i++)
{
    String line = "";
    for (int j=0; j<idOverlaps[i].length; j++)
    {
        if (j > 0)  line = line + "\t";
        line = line + idOverlaps[i][j];
    }
    System.err.println(line);
}
                System.err.println();
System.err.println("differences");
                System.err.println(headerLine);

for (int i=0; i<dissimilarities.length; i++)
{
    String line = "";
    for (int j=0; j<dissimilarities[i].length; j++)
    {
        if (j > 0)  line = line + "\t";
        line = line + dissimilarities[i][j];
    }
    System.err.println(line);
}

                File imageFile = TempFileManager.createTempFile("clust.png", this);
                String rExpression = "png('" + imageFile.getAbsolutePath()+ "',1000,1000); " +
                        "plot(hclust(as.dist(dissimilarities))); " +
                        "dev.off();";
//System.err.println(rExpression);

                RInterface.evaluateRExpression(rExpression, null, matrixMap, null);
                PanelWithBlindImageChart pwbic = new PanelWithBlindImageChart(imageFile, "Cluster");
                pwbic.displayInTab();
                TempFileManager.deleteTempFiles(this);
                break;
            case MODE_PLOT_SPECTRAL_COUNTS:
                if (featureSetsToHandle.length > 2)
                    return result;
                Map<String, Integer>[] peptideSpectralCountsMaps = new HashMap[2];
                for (int i=0; i<2; i++)
                {
                    Map<String, Integer> peptideCountMap = new HashMap<String, Integer>();
                    for (Feature feature : featureSetsToHandle[i].getFeatures())
                    {
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                        if (peptide == null)
                            continue;
                        if (!peptideCountMap.containsKey(peptide))
                            peptideCountMap.put(peptide, 0);
                        peptideCountMap.put(peptide, peptideCountMap.get(peptide) + 1);
                    }
                    peptideSpectralCountsMaps[i] = peptideCountMap;
                }
                List<Float> counts1 = new ArrayList<Float>();
                List<Float> counts2 = new ArrayList<Float>();

                for (String peptide : peptideSpectralCountsMaps[0].keySet())
                {
                    
                    if (peptideSpectralCountsMaps[1].containsKey(peptide))
                    {
                        counts1.add(peptideSpectralCountsMaps[0].get(peptide).floatValue() + (float) (Math.random() * 0.5 - 0.25));
                        counts2.add(peptideSpectralCountsMaps[1].get(peptide).floatValue()+ (float)(Math.random() * 0.5 - 0.25));

                    }
                }
                PanelWithScatterPlot pwsp = new PanelWithScatterPlot(counts1, counts2, "Spectral Counts");
                pwsp.setAxisLabels("Set 1", "Set 2");
                pwsp.displayInTab();
                break;
            case MODE_PLOT_PAIRWISE_RATIOS:
                List<Map<String, Float>> peptideRatioMaps = new ArrayList<Map<String, Float>>();
                Set<String> allPeptides = new HashSet<String>();
                for (FeatureSet featureSet : featureSets)
                {
                    Map<String, List<Double>> peptideRatiosMap = new HashMap<String, List<Double>>();
                    for (Feature feature : featureSet.getFeatures())
                    {
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                        double ratio = IsotopicLabelExtraInfoDef.getRatio(feature);
                        if (peptide == null || ratio == IsotopicLabelExtraInfoDef.NO_RATIO_FOR_FEATURE)
                            continue;
                        List<Double> ratios = peptideRatiosMap.get(peptide);
                        if (ratios == null)
                        {
                            ratios = new ArrayList<Double>();
                            peptideRatiosMap.put(peptide, ratios);
                        }
                        ratios.add(ratio);
                    }
                    Map<String, Float> peptideMeanRatioMap = new HashMap<String, Float>();
                    for (String peptide : peptideRatiosMap.keySet())
                    {
                        allPeptides.add(peptide);
                        peptideMeanRatioMap.put(peptide, (float) BasicStatistics.geometricMean(peptideRatiosMap.get(peptide)));
                    }
                    peptideRatioMaps.add(peptideMeanRatioMap);
                }
                List<List<Double>> pairwisePlotLogRatioLists = new ArrayList<List<Double>>();
                PanelWithRPairsPlot pairsPlot = new PanelWithRPairsPlot();
                pairsPlot.setName("Pairwise Ratios");
                for (Map<String, Float> peptideMeanRatioMap : peptideRatioMaps)
                {
                    List<Double> meanRatios = new ArrayList<Double>();
                    for (String peptide : allPeptides)
                    {
                        if (peptideMeanRatioMap.containsKey(peptide))
                            meanRatios.add(Math.log(peptideMeanRatioMap.get(peptide)));
                        else
                            meanRatios.add(Double.NaN);
                    }
                    pairwisePlotLogRatioLists.add(meanRatios);
                }
                pairsPlot.setChartWidth(1000);
                pairsPlot.setChartHeight(1000);                

                pairsPlot.plot(pairwisePlotLogRatioLists, true);
                pairsPlot.displayInTab();

                break;
        }

        return result;
    }

    /**
     * Build a peptide intensity map for a single charge state
     * @param features
     * @param charge
     * @return
     */
    protected Map<String, Float> createPeptideIntensityMapForChargeState(
            Feature[] features,
            int charge)
    {
        Map<String, Float> peptideIntensityMap =
                new HashMap<String,Float>();
        for (Feature feature : features)
        {
//if (feature.getIntensity() == 0f) System.err.println("0");            
            if (feature.charge != charge)
                continue;
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null || MS2ExtraInfoDef.getPeptideList(feature).size() > 1)
                continue;
            Float existingIntensity = peptideIntensityMap.get(peptide);
            if (existingIntensity == null)
                existingIntensity = 0f;
            peptideIntensityMap.put(peptide, existingIntensity + feature.getIntensity());
//if (peptideIntensityMap.get(peptide) == 0f) System.err.println("GGAAAAAAHHHH");
        }
        return peptideIntensityMap;
    }

    /**
     * Heavy-duty stuff.  Ratio calculation.  Very important.  Should move elsewhere
     * @throws CommandLineModuleExecutionException
     */
    protected void calculateRatios(FeatureSet[] featureSetsToHandle)
            throws CommandLineModuleExecutionException
    {
        FeatureSet featureSet1 = featureSetsToHandle[0];
        FeatureSet featureSet2 = featureSetsToHandle[1];

        List<Float> allIntensities1 = new ArrayList<Float>();
        List<Float> allIntensities2 = new ArrayList<Float>();

        Map<String, Map<Integer,Float>> peptideIntensitiesAllCharges1 =
                new HashMap<String, Map<Integer,Float>>();
        Map<String, Map<Integer,Float>> peptideIntensitiesAllCharges2 =
                new HashMap<String, Map<Integer,Float>>();

        Set<String> allPeptides = new HashSet<String>();
        Set<String> commonPeptides = new HashSet<String>();
        Set<String> peptides1 = new HashSet<String>();


        //calculate ratios for each charge state.  Keep track of log intensity pairs
        //hardcoded max charge state of 9.  Boo hoo.
        for (int charge=1; charge<10; charge++)
        {
            Map<String, Float> peptideIntensityMap1ThisCharge =
                      createPeptideIntensityMapForChargeState(featureSet1.getFeatures(), charge);
            Map<String, Float> peptideIntensityMap2ThisCharge =
                      createPeptideIntensityMapForChargeState(featureSet2.getFeatures(), charge);
            for (String peptide1 : peptideIntensityMap1ThisCharge.keySet())
            {
                allPeptides.add(peptide1);
                peptides1.add(peptide1);
                float intensity1 = peptideIntensityMap1ThisCharge.get(peptide1);
                Map<Integer,Float> intensitiesThisPeptide1 = peptideIntensitiesAllCharges1.get(peptide1);
                if (intensitiesThisPeptide1 == null)
                {
                    intensitiesThisPeptide1 = new HashMap<Integer,Float>();
                    peptideIntensitiesAllCharges1.put(peptide1, intensitiesThisPeptide1);
                }
                Map<Integer,Float> intensitiesThisPeptide2 = peptideIntensitiesAllCharges2.get(peptide1);
                if (intensitiesThisPeptide2 == null)
                {
                    intensitiesThisPeptide2 = new HashMap<Integer,Float>();
                    peptideIntensitiesAllCharges2.put(peptide1, intensitiesThisPeptide2);
                }

                Float intensity2 = peptideIntensityMap2ThisCharge.get(peptide1);
                if (intensity2 == null)
                    intensity2 = 0f;
                intensitiesThisPeptide1.put(charge, intensity1);
                intensitiesThisPeptide2.put(charge, intensity2);

                allIntensities1.add(intensity1);
                allIntensities2.add(intensity2);
            }
            for (String peptide2 : peptideIntensityMap2ThisCharge.keySet())
            {
                allPeptides.add(peptide2);
                if (peptides1.contains(peptide2))
                    commonPeptides.add(peptide2);
                float intensity2 = peptideIntensityMap2ThisCharge.get(peptide2);

                if (!peptideIntensityMap1ThisCharge.containsKey(peptide2))
                {
                    Map<Integer,Float> intensitiesThisPeptide1 = peptideIntensitiesAllCharges1.get(peptide2);
                    if (intensitiesThisPeptide1 == null)
                    {
                        intensitiesThisPeptide1 = new HashMap<Integer,Float>();
                        peptideIntensitiesAllCharges1.put(peptide2, intensitiesThisPeptide1);
                    }
                    Map<Integer,Float> intensitiesThisPeptide2 = peptideIntensitiesAllCharges2.get(peptide2);
                    if (intensitiesThisPeptide2 == null)
                    {
                        intensitiesThisPeptide2 = new HashMap<Integer,Float>();
                        peptideIntensitiesAllCharges2.put(peptide2, intensitiesThisPeptide2);
                    }

                    intensitiesThisPeptide1.put(charge, 0f);
                    intensitiesThisPeptide2.put(charge, intensity2);

                    allIntensities1.add(0f);
                    allIntensities2.add(intensity2);
                }
            }
        }

        ApplicationContext.infoMessage("Peptides in common: " + commonPeptides.size());


        if (normalize)
        {
            List<float[]> rowsForNormalization = new ArrayList<float[]>();

            for (String peptide : allPeptides)
            {
                Map<Integer,Float> thisPeptideIntensitiesAllCharges1 =
                        peptideIntensitiesAllCharges1.get(peptide);
                Map<Integer,Float> thisPeptideIntensitiesAllCharges2 =
                        peptideIntensitiesAllCharges2.get(peptide);

                for (int i=0; i<10; i++)
                {
                    if (thisPeptideIntensitiesAllCharges1.containsKey(i) &&
                        thisPeptideIntensitiesAllCharges2.containsKey(i))
                    {
                        float[] thisRow = new float[2];

                        thisRow[0] = thisPeptideIntensitiesAllCharges1.get(i);
                        thisRow[1] = thisPeptideIntensitiesAllCharges2.get(i);

                        rowsForNormalization.add(thisRow);
                    }
                }
            }

            Normalizer.normalize(rowsForNormalization);

            int i=0;
            for (String peptide : allPeptides)
            {
                Map<Integer,Float> thisPeptideIntensitiesAllCharges1 =
                        peptideIntensitiesAllCharges1.get(peptide);
                Map<Integer,Float> thisPeptideIntensitiesAllCharges2 =
                        peptideIntensitiesAllCharges2.get(peptide);
                int numAssigned=0;
                for (int j=0; j<10; j++)
                {
                    if (thisPeptideIntensitiesAllCharges1.containsKey(j) &&
                        thisPeptideIntensitiesAllCharges2.containsKey(j))
                    {
//System.err.println("Changing " + thisPeptideIntensitiesAllCharges1.get(j) + " to " + rowsForNormalization.get(i+j)[0]);
                        thisPeptideIntensitiesAllCharges1.put(j, rowsForNormalization.get(i)[0]);
                        thisPeptideIntensitiesAllCharges2.put(j, rowsForNormalization.get(i)[1]);
                        i++;
                    }
                }
            }

            ApplicationContext.infoMessage("Done with  normalization.");

//            if (showCharts)
//            {
//                float[] run1LogIntensities = new float[peptidesToNormalizeList.size()];
//                float[] run2LogIntensities = new float[peptidesToNormalizeList.size()];
//
//                int i = 0;
//                for (String peptide : peptidesToNormalizeList)
//                {
//                    run1LogIntensities[i] = (float) Math.log(Math.max(peptideIntensity1Map.get(peptide), 0.0001));
//                    run2LogIntensities[i] = (float) Math.log(Math.max(peptideIntensity2Map.get(peptide), 0.0001));
//                    i++;
//                }
//
//                ScatterPlotDialog spdLog =
//                        new ScatterPlotDialog(run1LogIntensities, run2LogIntensities, "After Normalization: Log intensities");
//                spdLog.setAxisLabels("Run 1", "Run 2");
//                spdLog.setVisible(true);
//            }
        }




        Map<String, Float> finalPeptideRatioMap = new HashMap<String, Float>();
        Map<String, Float> finalPeptideIntensity1Map = new HashMap<String, Float>();
        Map<String, Float> finalPeptideIntensity2Map = new HashMap<String, Float>();



        List<Integer> numChargeStatesList = new ArrayList<Integer>();
        int numActualRatios = 0;
        for (String peptide : allPeptides)
        {
            Map<Integer,Float> thisPeptideIntensitiesAllCharges1 =
                    peptideIntensitiesAllCharges1.get(peptide);
            Map<Integer,Float> thisPeptideIntensitiesAllCharges2 =
                    peptideIntensitiesAllCharges2.get(peptide);

            List<Float> actualRatios = new ArrayList<Float>();
            int numMaxes = 0;
            int numMins = 0;

            float intensity1IfActual=0;
            float intensity2IfActual = 0;
            float intensity1IfMax=0;
            float intensity2IfMin = 0;
            for (int i=0; i<10; i++)
            {
                if (thisPeptideIntensitiesAllCharges1.containsKey(i) &&
                    thisPeptideIntensitiesAllCharges2.containsKey(i))
                {
                    float intensity1 = thisPeptideIntensitiesAllCharges1.get(i);
                    float intensity2 = thisPeptideIntensitiesAllCharges2.get(i);
//System.err.println(intensity1 + " vs. " + intensity2);
                    if (intensity1 == 0f)
                    {
                        numMins++;
                        intensity2IfMin = intensity2;
                    }
                    else if (intensity2 == 0f)
                    {
                        numMaxes++;
                        intensity1IfMax = intensity1;
                    }
                    else
                    {
                        actualRatios.add(intensity1 / intensity2);
                        //not ideal.
                        intensity1IfActual = intensity1;
                        intensity2IfActual = intensity2;
                    }
                }
            }
            float ratio = -1f;
            if (actualRatios.size() > 0)
            {
                numActualRatios++;
                ratio = (float)BasicStatistics.mean(actualRatios);
                numChargeStatesList.add(actualRatios.size());
                        finalPeptideIntensity1Map.put(peptide, intensity1IfActual);
                        finalPeptideIntensity2Map.put(peptide, intensity2IfActual);
            }
            else
            {
                if (numMaxes > 0)
                {
                    if (numMins == 0)
                    {
                        finalPeptideIntensity1Map.put(peptide, intensity1IfMax);
                        finalPeptideIntensity2Map.put(peptide, 0f);
                        ratio = INFINITE_RATIO;
                    }
                }
                else
                {
                    if (numMins > 0)
                    {
                        finalPeptideIntensity2Map.put(peptide, intensity2IfMin);
                        finalPeptideIntensity1Map.put(peptide, 0f);
                        ratio = ZERO_RATIO;
                    }
                }
            }


            if (ratio > -1f)
            {
                if (boundPeptideRatios)
                {
                    if (ratio <= .05f)
                    {
                        ratio = .05f;
                        finalPeptideIntensity1Map.put(peptide, .05f * finalPeptideIntensity2Map.get(peptide));
                    }
                    else if (ratio >= 20f)
                    {

                        ratio = 20f;
                        finalPeptideIntensity2Map.put(peptide, .05f * finalPeptideIntensity1Map.get(peptide));
                    }
                }
                finalPeptideRatioMap.put(peptide, ratio);
            }
        }

        ApplicationContext.infoMessage("Comparable peptides: " + finalPeptideRatioMap.size());
        ApplicationContext.infoMessage("Non-infinite ratios: " + numActualRatios);


        if (showCharts)
        {
            float[] numChargeStatesArray = new float[numChargeStatesList.size()];
            for (int i=0; i<numChargeStatesList.size(); i++)
                numChargeStatesArray[i] = numChargeStatesList.get(i);
            PanelWithHistogram pwhChargeStates = new PanelWithHistogram(numChargeStatesArray);
            ChartDialog chargeStateDialog = new ChartDialog(pwhChargeStates);
            chargeStateDialog.setTitle("# charge states per peptide");
            chargeStateDialog.setVisible(true);

            float[] logIntensities1Bounded = new float[allIntensities1.size()];
            float[] logIntensities2Bounded = new float[allIntensities2.size()];
            for (int i=0; i<allIntensities1.size(); i++)
            {
                logIntensities1Bounded[i] = (float) Math.log(Math.max(allIntensities1.get(i), 0.001));
                logIntensities2Bounded[i] = (float) Math.log(Math.max(allIntensities2.get(i), 0.001));
            }

            ScatterPlotDialog spdLog =
                    new ScatterPlotDialog(logIntensities1Bounded, logIntensities2Bounded, "Log intensities");
                String xAxisLabelToUse = "Run 1";
                if (xAxisLabel != null)
                    xAxisLabelToUse = xAxisLabel;
                String yAxisLabelToUse = "Run 2";
                if (yAxisLabel != null)
                    yAxisLabelToUse = yAxisLabel;
            spdLog.setAxisLabels(xAxisLabelToUse, yAxisLabelToUse);
            spdLog.setVisible(true);
        }



        PrintWriter outPW = null;
        boolean tabOutput = (outFile != null && !outFormatPepXML);
        boolean pepXmlOutput = (outFile != null && outFormatPepXML);
        if (tabOutput)
        {
            try
            {
                outPW = new PrintWriter(outFile);
                outPW.println("peptide\tratio\tintensity1\tintensity2");
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        List<Feature> ratioFeatures = new ArrayList<Feature>();
        for (String peptide : finalPeptideRatioMap.keySet())
        {

            float intensity1 = finalPeptideIntensity1Map.get(peptide);
            float intensity2 = finalPeptideIntensity2Map.get(peptide);
            float ratio = finalPeptideRatioMap.get(peptide);

            if (tabOutput)
                outPW.println(peptide + "\t" + ratio + "\t" + intensity1 + "\t" + intensity2);


            Feature ratioFeature = new Feature();
            MS2ExtraInfoDef.addPeptide(ratioFeature, peptide);
            MS2ExtraInfoDef.setPeptideProphet(ratioFeature, .95);
            IsotopicLabelExtraInfoDef.setRatio(ratioFeature, ratio);
            IsotopicLabelExtraInfoDef.setLightIntensity(ratioFeature, intensity1);
            IsotopicLabelExtraInfoDef.setHeavyIntensity(ratioFeature, intensity2);
            float featureMass =
                    (float) PeptideGenerator.computeMass(peptide.getBytes(), 0, peptide.length(), PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
            ratioFeature.setMass(featureMass);
            ratioFeature.setScan(1);
            ratioFeature.setScanFirst(1);
            ratioFeature.setScanLast(1);
            ratioFeature.setCharge(2);
            
            ratioFeatures.add(ratioFeature);
        }

        if (pepXmlOutput)
        {
            FeatureSet outRatioFeatureSet = new FeatureSet(ratioFeatures.toArray(new Feature[0]));
            if (fastaFile != null)
                MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(outRatioFeatureSet, fastaFile.getAbsolutePath());
                            try
                {
                    outRatioFeatureSet.savePepXml(outFile);
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }

        }

        if (showCharts)
        {
            double[] logRatios = new double[finalPeptideRatioMap.size()];
            int i=0;
            for(String peptide : finalPeptideRatioMap.keySet())
            {
                logRatios[i++] = Math.log(Math.max(finalPeptideRatioMap.get(peptide),.0001));
            }
            PanelWithHistogram pwh = new PanelWithHistogram(logRatios, "Log Ratios");
            ChartDialog hd = new ChartDialog(pwh);
            hd.setVisible(true);
        }

        if (outPW != null)
        {
            outPW.flush();
            outPW.close();
        }


    }


    protected void plotSomeAttribute(FeatureSet[] featureSetsToHandle)
    {
        Set<String> commonPeptides = new HashSet<String>();

        List<Map<String,Map<Integer, List<Feature>>>> peptidesInSetsChargeMap =
                new ArrayList<Map<String,Map<Integer, List<Feature>>>>();
        List<Map<String,List<Feature>>> peptidesInSetsMap =
                new ArrayList<Map<String,List<Feature>>>();


        List<Set<String>> peptidesInSets = new ArrayList<Set<String>>();

        List<Set<String>> comparablePeptidesInSets = new ArrayList<Set<String>>();


        int fSetNum=0;
        for (FeatureSet featureSet : featureSetsToHandle)
        {
            fSetNum++;
System.err.println("Set: " + featureSet.getSourceFile());
System.err.println("Features: " + featureSet.getFeatures().length);
            Map<String,Map<Integer,List<Feature>>> thisSetPeptideChargeFeatureMap =
                    new HashMap<String,Map<Integer,List<Feature>>>();
            Set<String> thisSetPeptides = new HashSet<String>();
            Set<String> thisSetComparablePeptides = new HashSet<String>();

            Map<String,List<Feature>> thisSetPeptideFeatureMap =
                    new HashMap<String,List<Feature>>();
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null)
                {
                    thisSetPeptides.add(peptide);
                    boolean isComparable = false;
                    switch(mode)
                    {
                        case MODE_PLOT_INTENSITIES:
                        case MODE_PLOT_TOTAL_INTENSITIES:
                        case MODE_PLOT_TIMES:
                        case MODE_PLOT_PEPTIDEPROPHET:
                        case MODE_PLOT_EXPECT:
                            //not actually correct
                            isComparable = true;
                        case MODE_PLOT_FVAL:
                            //not actually correct
                            isComparable = true;
                            break;
                        case MODE_PLOT_KSCORE_OR_XCORR:
                            //not actually correct
                            break;
                        case MODE_PLOT_RATIOS:
                            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                                isComparable= true;
                            break;
                        case MODE_PLOT_LIGHT_AREAS:
                            if (IsotopicLabelExtraInfoDef.getLightIntensity(feature) > 0)
                                isComparable = true;
                            break;

                    }


                    if (isComparable)
                    {
                        thisSetComparablePeptides.add(peptide);

                        if (withinCharge)
                        {
                            if (feature.getCharge() > 0)
                            {
                                Map<Integer,List<Feature>> chargeFeatureListMap = thisSetPeptideChargeFeatureMap.get(peptide);
                                if (chargeFeatureListMap == null)
                                {
                                    chargeFeatureListMap = new HashMap<Integer, List<Feature>>();
                                    thisSetPeptideChargeFeatureMap.put(peptide, chargeFeatureListMap);
                                }
                                List<Feature> featureListThisCharge = chargeFeatureListMap.get(feature.getCharge());
                                if (featureListThisCharge == null)
                                {
                                    featureListThisCharge = new ArrayList<Feature>();
                                    chargeFeatureListMap.put(feature.getCharge(), featureListThisCharge);
                                }
                                featureListThisCharge.add(feature);
                                thisSetPeptides.add(peptide);
                            }
                        }
                        else
                        {
                            List<Feature> featureList = thisSetPeptideFeatureMap.get(peptide);
                            if (featureList == null)
                            {
                                featureList = new ArrayList<Feature>();
                                thisSetPeptideFeatureMap.put(peptide,featureList);
                            }
                            featureList.add(feature);
                        }

                    }
                }
            }
            peptidesInSets.add(thisSetPeptides);
            comparablePeptidesInSets.add(thisSetComparablePeptides);

            if (withinCharge)
                peptidesInSetsChargeMap.add(thisSetPeptideChargeFeatureMap);
            else
                peptidesInSetsMap.add(thisSetPeptideFeatureMap);

System.err.println("peptides: " + thisSetPeptides.size() + ", comparable: " + thisSetComparablePeptides.size());
        }


        commonPeptides.addAll(peptidesInSets.get(0));
        for (int i=1; i<peptidesInSets.size(); i++)
        {
            Set<String> setPeptides = peptidesInSets.get(i);
            Set<String> commonPeptidesTemp = new HashSet<String>();
            for (String peptide : commonPeptides)
            {
                if (setPeptides.contains(peptide))
                    commonPeptidesTemp.add(peptide);
            }
            commonPeptides = commonPeptidesTemp;
        }

        List<Float> set1Values = new ArrayList<Float>();
        List<Float> set2Values = new ArrayList<Float>();
        int i=0;

        Set<String> commonComparablePeptides = new HashSet<String>();
        for (String pep1 : comparablePeptidesInSets.get(0))
            if (comparablePeptidesInSets.get(1).contains(pep1))
                commonComparablePeptides.add(pep1);
            else
                if (!withinCharge)
                {
                    List<Feature> features = peptidesInSetsMap.get(0).get(pep1);
                    float maxVal = 0f;
                    for (Feature feature : features)
                    {
                        float val =  (float)getRelevantValueFromFeature(feature);
                        if (val > maxVal)
                             maxVal = val;
                    }
                    set1OnlyValues.add(maxVal);
                }
        for (String pep2 : comparablePeptidesInSets.get(1))
            if (!commonComparablePeptides.contains(pep2))
                if (!withinCharge)
                {
                    List<Feature> features = peptidesInSetsMap.get(1).get(pep2);
                    float maxVal = 0f;
                    for (Feature feature : features)
                    {
                        float val =  (float)getRelevantValueFromFeature(feature);
                        if (val > maxVal)
                             maxVal = val;
                    }
                    set2OnlyValues.add(maxVal);
                }


        boolean set1HasKscore = false;
        boolean set2HasKscore = false;

        List<List<Float>> uncomparableValuesInSets = new ArrayList<List<Float>>();

        int setIndex = 0;
        for (FeatureSet featureSet : featureSetsToHandle)
        {
            List<Float> uncomparableValues = new ArrayList<Float>();

            Set<String> thisSetComparablePeptides = comparablePeptidesInSets.get(setIndex++);
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null && !commonComparablePeptides.contains(peptide) &&
                        thisSetComparablePeptides.contains(peptide))
                     uncomparableValues.add((float)getRelevantValueFromFeature(feature));

            }

            uncomparableValuesInSets.add(uncomparableValues);

        }

        int numEqual=0;
        for (String peptide : commonComparablePeptides)
        {
            if (withinCharge)
            {
                for (int charge = 0; charge<7; charge++)
                {
                    List<Feature> features1 = peptidesInSetsChargeMap.get(0).get(peptide).get(charge);
                    List<Feature> features2 = peptidesInSetsChargeMap.get(1).get(peptide).get(charge);
                    if (features1 == null || features2 == null)
                        continue;
                    List<Float> values1 = new ArrayList<Float>();
                    for (Feature feature : features1)
                        values1.add((float) getRelevantValueFromFeature(feature));
                    List<Float> values2 = new ArrayList<Float>();
                    for (Feature feature : features2)
                        values2.add((float) getRelevantValueFromFeature(feature));    

                    float mean1 = (float)BasicStatistics.mean(values1);
                    float mean2 = (float)BasicStatistics.mean(values2);

                    set1Values.add(mean1);
                    set2Values.add(mean2);
//if (peptide.equals("DCQDIANK")) System.err.println("****" + mean1 + ", " + mean2);
                    if (mean1==mean2)
                        numEqual++;
                }
            }
            else
            {
                //use highest-probability feature
                Double set1Value = 0.0;
                if (comparablePeptidesInSets.get(0).contains(peptide))
                {
                    Feature feature1 = null;
                    float maxProb1 = -1;
                    for (Feature feature : peptidesInSetsMap.get(0).get(peptide))
                        if (MS2ExtraInfoDef.getPeptideProphet(feature) > maxProb1)
                        {
                            feature1 = feature;
                            maxProb1 = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                        }
                    set1Value=getRelevantValueFromFeature(feature1);
                }

                Double set2Value = 0.0;
                if (comparablePeptidesInSets.get(1).contains(peptide))
                {
                    Feature feature2 = null;
                    float maxProb2 = -1;
                    for (Feature feature : peptidesInSetsMap.get(1).get(peptide))
                        if (MS2ExtraInfoDef.getPeptideProphet(feature) > maxProb2)
                        {
                            feature2 = feature;
                            maxProb2 = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                        }
                    set2Value=getRelevantValueFromFeature(feature2);
//if (Math.abs(set1Value-set2Value) > .005) System.err.println(peptide + "set 1: " + set1Value + " set 2: " + set2Value);    

                }

                set1Values.add(set1Value.floatValue());
                set2Values.add(set2Value.floatValue());
                if (set1Value.floatValue() == set2Value.floatValue())
                    numEqual++;
            }
            i++;
        }

        double[] set1LogValues = new double[set1Values.size()];
        double[] set2LogValues = new double[set1Values.size()];
//System.err.println("CALCING LOGS");
        for (int j=0; j<set1Values.size(); j++)
        {
            set1LogValues[j] = Math.log(Math.max(set1Values.get(j), 0.001));
            set2LogValues[j] = Math.log(Math.max(set2Values.get(j), 0.001));
//if (Double.isNaN(set1LogValues[j]) || Double.isInfinite(set1LogValues[j])) System.err.println("BAD1!!!! " + (set1Values.get(j) + ", " + set1LogValues[j]));

//if (Double.isNaN(set2LogValues[j]) || Double.isInfinite(set2LogValues[j])) System.err.println("BAD2!!!! " + (set2Values.get(j) + ", " + set2LogValues[j]));
//System.err.println(set2LogValues[j]);            
        }


        System.err.println("Comparable peptides in set 1: " + comparablePeptidesInSets.get(0).size() +
                           ", set 2: " + comparablePeptidesInSets.get(1).size() + ", in common: " +
                           commonComparablePeptides.size());
        System.err.println("Comparable peptides unique to 1: " +
                           (comparablePeptidesInSets.get(0).size() -commonComparablePeptides.size()) +
                           ", unique to 2: " +
                           (comparablePeptidesInSets.get(1).size() -commonComparablePeptides.size())); 

        System.err.println("Peptides in common: " + commonPeptides.size());
        System.err.println("Comparable features: " + set1LogValues.length);
        System.err.println("# with same value: " + numEqual);

        if (featureSets != null && showCharts)
        {
            PanelWithScatterPlot spd =
                    new PanelWithScatterPlot(set1Values, set2Values, "Common Peptide Attributes");

            List<Float> set1ValuesLog = new ArrayList<Float>();
            List<Float> set2ValuesLog = new ArrayList<Float>();
            for (int j=0; j<set1Values.size(); j++) {
                float set1 = 0;
                float set2 = 0;
                try {
                    set1 =  (float) Math.log(set1Values.get(j));
                    set2 = (float) Math.log(set2Values.get(j));
                    set1ValuesLog.add(set1);
                    set2ValuesLog.add(set2);
                } catch (Exception e) {

                }
            }
            PanelWithScatterPlot spdLog =
                    new PanelWithScatterPlot(set1ValuesLog, set2ValuesLog, "Common Peptide Attributes (log)");
            switch (mode)
        {
                case MODE_PLOT_INTENSITIES:
                    spd.setAxisLabels("Set 1 intensity", "Set 2 intensity");
                    break;
                case MODE_PLOT_TIMES:
                    spd.setAxisLabels("Set 1 time", "Set 2 time");
                    break;
                case MODE_PLOT_PEPTIDEPROPHET:
                    spd.setAxisLabels(featureSetsToHandle[0].getSourceFile().getName() + " Prob",
                            featureSetsToHandle[1].getSourceFile().getName() + " Prob");
                    break;
                case MODE_PLOT_EXPECT:
                    spd.setAxisLabels(featureSetsToHandle[0].getSourceFile().getName() + " Expect",
                            featureSetsToHandle[1].getSourceFile().getName() + " Expect");
                    break;
                case MODE_PLOT_FVAL:
                    spd.setAxisLabels("Set 1 PeptideProphet fval",
                            "Set 2 PeptideProphet fval");
                    break;
                case MODE_PLOT_KSCORE_OR_XCORR:
                    spd.setAxisLabels(set1HasKscore ? "K score" : "XCorr",
                            set2HasKscore ? "K score" : "XCorr");
                    break;
                case MODE_PLOT_RATIOS:
                    spd.setAxisLabels(featureSetsToHandle[0].getSourceFile().getName() + " ratio",
                            featureSetsToHandle[1].getSourceFile().getName() + " ratio");
                    break;
            }
            if (xAxisLabel != null && yAxisLabel != null)
                spd.setAxisLabels(xAxisLabel, yAxisLabel);

            spd.displayInTab();

            spdLog.setAxisLabels(spd.getXAxis().getLabel() + " (log)", spd.getYAxis().getLabel() + " (log)");
            spdLog.displayInTab();
        }
List<Float> valueDifferences = new ArrayList<Float>();
for (int j=0; j<set1Values.size(); j++)
{
    float diff = set2Values.get(j) - set1Values.get(j);
    if (!Float.isInfinite(diff) && !Float.isNaN(diff))
        valueDifferences.add(diff);
}
ApplicationContext.infoMessage("Mean difference of values (set 2 - set 1): " + BasicStatistics.mean(valueDifferences) +
        ".  Median: " + BasicStatistics.median(valueDifferences));        
        ApplicationContext.infoMessage("Correlation Coefficient: " +
                BasicStatistics.correlationCoefficient(set1Values, set2Values));


        //add uncomparable values

        float minLog1 = (float) BasicStatistics.min(set1LogValues);
        float minLog2 = (float) BasicStatistics.min(set2LogValues);
        List<Float> uncompLog1 = new ArrayList<Float>();
        List<Float> minLog1List = new ArrayList<Float>();
        List<Float> minLog2List = new ArrayList<Float>();

        for (float val : uncomparableValuesInSets.get(0))
        {
            float logVal = (float)Math.log(val);
            if (!Float.isNaN(logVal) && !Float.isInfinite(logVal))
            {
                uncompLog1.add(logVal);
                minLog2List.add(minLog2);
            }

        }
        List<Float> uncompLog2 = new ArrayList<Float>();
        for (float val : uncomparableValuesInSets.get(1))
        {
            float logVal = (float)Math.log(val);
            if (!Float.isNaN(logVal) && !Float.isInfinite(logVal))
            {
                uncompLog2.add(logVal);
                minLog1List.add(minLog1);
            }

        }




        if (showCharts && (mode == MODE_PLOT_INTENSITIES || mode == MODE_PLOT_TOTAL_INTENSITIES ||
                           mode == MODE_PLOT_RATIOS) &&
                featureSets != null)
        {
            PanelWithScatterPlot spdLog =
                    new PanelWithScatterPlot(set1LogValues, set2LogValues, "Common Peptide values (log)");
            spdLog.addData(uncompLog1, minLog2List, "Uncomparable in 1");
            spdLog.addData(minLog1List, uncompLog2, "Uncomparable in 2");

            spdLog.setAxisLabels("Set 1 Log", "Set 2 Log");
            spdLog.displayInTab();
        }
        //if we're doing a directory-level thing
        if (featureSets == null)
        {
            for (int j=0; j<set1Values.size(); j++)
            {
                allSet1Values.add(set1Values.get(j));
                allSet2Values.add(set2Values.get(j));
            }
            System.err.println("all plot: " + allSet1Values.size());
            for (int j=0; j<set1LogValues.length; j++)
            {
                if (!Double.isNaN(set1LogValues[j]) && !Double.isInfinite(set1LogValues[j]) &&
                        !Double.isNaN(set2LogValues[j]) && !Double.isInfinite(set2LogValues[j] ))
                {
                    allSet1LogValues.add((float) set1LogValues[j]);
                    allSet2LogValues.add((float) set2LogValues[j]);
                }
            }
        }

        if (featureSets != null)
        {
            for (int j=0; j<featureSetsToHandle.length; j++)
            {
                if (uncomparableValuesInSets.get(j).isEmpty())
                    continue;
                PanelWithHistogram uncomparableHist =
                        new PanelWithHistogram(uncomparableValuesInSets.get(j), "Uncomparable " + (j+1));
                uncomparableHist.displayInTab();
            }
            for (int j=0; j<featureSetsToHandle.length; j++)
            {
                if (uncomparableValuesInSets.get(j).isEmpty())
                    continue;
                List<Float> uncompLog = new ArrayList<Float>();
                for (float val : uncomparableValuesInSets.get(j))
                     uncompLog.add((float)Math.log(val));
                PanelWithHistogram uncomparableHist =
                        new PanelWithHistogram(uncompLog, "Log Uncomparable " + (j+1));
                uncomparableHist.displayInTab();
            }
        }


    }

    protected double getRelevantValueFromFeature(Feature feature)
    {

        double result = 0;

        switch(mode)
        {
            case MODE_PLOT_INTENSITIES:
                if (feature != null) result = (double) feature.getIntensity();
                break;
            case MODE_PLOT_TOTAL_INTENSITIES:
                if (feature != null) result = (double) feature.getTotalIntensity();
                break;
            case MODE_PLOT_TIMES:
                if (feature != null) result = (double) feature.getTime();
                break;
            case MODE_PLOT_PEPTIDEPROPHET:
                if (feature != null) result = MS2ExtraInfoDef.getPeptideProphet(feature);
//if (Math.abs(set1Value-result) > 0.1) System.err.println(peptide + ", " + set1Value + ", " + result);

                break;
            case MODE_PLOT_EXPECT:
                if (feature != null) result = Double.parseDouble(MS2ExtraInfoDef.getSearchScore(feature, "expect"));
//if (Math.abs(set1Value-result) > 0.1) System.err.println(peptide + ", " + set1Value + ", " + result);

                break;
            case MODE_PLOT_FVAL:
                if (feature != null) {
//                    System.err.println(MS2ExtraInfoDef.getFval(feature));
                    result = MS2ExtraInfoDef.getFval(feature);
                }
                break;
            case MODE_PLOT_KSCORE_OR_XCORR:
                if (feature != null && MS2ExtraInfoDef.getSearchScore(feature,"dotproduct") != null)
                    result = MS2ExtraInfoDef.getDoubleSearchScore(feature,"dotproduct");
                else
                    if (feature != null) result = MS2ExtraInfoDef.getDoubleSearchScore(feature,"xcorr");
                break;
            case MODE_PLOT_RATIOS:
                    if (feature != null) result = IsotopicLabelExtraInfoDef.getRatio(feature);
                break;
            case MODE_PLOT_LIGHT_AREAS:
                      if (feature != null)
                          result = IsotopicLabelExtraInfoDef.getLightIntensity(feature);
                break;
        }
        return result;
    }



    /**
     * Note: _last_ featureset specified takes precedence in terms of feature location
     * @throws CommandLineModuleExecutionException
     */
    protected void createFeatureFileForPeptideUnion(FeatureSet[] featureSetsToHandle)
            throws CommandLineModuleExecutionException
    {
        Map<String, Feature> unionPeptideFeatures = new HashMap<String, Feature>();

        int i=1;
        for (FeatureSet featureSet : featureSetsToHandle)
        {
            Feature[] thisSetFeatures = featureSet.getFeatures();
            ApplicationContext.setMessage("Adding features from set with " +
                                          thisSetFeatures.length + " features");
            for (Feature feature : thisSetFeatures)
            {
                List<String> peptideList = MS2ExtraInfoDef.getPeptideList(feature);
                if (peptideList != null)
                {
                    for (String peptide : peptideList)
                    {
                        unionPeptideFeatures.put(peptide, feature);
                    }
                }
            }
        }

        Feature[] unionFeatures =
                unionPeptideFeatures.values().toArray(new Feature[unionPeptideFeatures.size()]);
        ApplicationContext.infoMessage("Creating combined featureset with " + unionFeatures.length + " features");
        FeatureSet outSet = new FeatureSet(unionFeatures);
        try
        {
            outSet.save(outFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

    protected double[][] calcOverlaps(FeatureSet[] featureSetsToHandle)
    {
        int numFeatureSets = featureSetsToHandle.length;

        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet((float)minPeptideProphet);

        Set<String>[] peptideSets = (Set<String>[]) new Set[numFeatureSets];
        Set<String>[] quantPeptideSets = (Set<String>[]) new Set[numFeatureSets];

        double[][] overlaps = new double[numFeatureSets][numFeatureSets];

        Set<String> allPeptides = new HashSet<String>();
        Set<String> allQuantPeptides = new HashSet<String>();

        for (int i=0; i<numFeatureSets; i++)
        {
            featureSetsToHandle[i] = featureSetsToHandle[i].filter(sel);

            Set<String> currentSetPeptides = new HashSet<String>();
            Set<String> currentSetQuantPeptides = new HashSet<String>();

            for (Feature feature : featureSetsToHandle[i].getFeatures())
            {
                List<String> peptideList = MS2ExtraInfoDef.getPeptideList(feature);
                if (peptideList != null)
                {
                    for (String peptide : peptideList)
                    {
                        currentSetPeptides.add(peptide);
                        allPeptides.add(peptide);
                        if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                        {
                            currentSetQuantPeptides.add(peptide);
                            allQuantPeptides.add(peptide);
                        }
                    }
                }
            }
            peptideSets[i] = currentSetPeptides;
            quantPeptideSets[i] = currentSetQuantPeptides;

        }





        for (int i=0; i<numFeatureSets; i++)
        {
            for (int j=0; j<numFeatureSets; j++)
            {
                if (i==j)
                {
                    overlaps[i][j] = peptideSets[i].size();
                    continue;
                }
                if (overlaps[j][i] > 0)
                {
                    overlaps[i][j] = overlaps[j][i];
                    continue;
                }
                int overlap = 0;
                for (String peptide : peptideSets[i])
                {
                    if (peptideSets[j].contains(peptide))
                    {
                        _log.debug(peptide);
                        overlap++;
                    }
                }
                overlaps[i][j] = overlap;
            }
        }
        return overlaps;
    }


    /**
     * Return the set of peptides in second but not first
     * @param featureSetsToHandle
     * @return
     */
    protected Pair<Set<String>,Set<String>> showOverlap(FeatureSet[] featureSetsToHandle)
            throws CommandLineModuleExecutionException
    {
        int numFeatureSets = featureSetsToHandle.length;

        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet((float)minPeptideProphet);

        Set<String>[] peptideSets = (Set<String>[]) new Set[numFeatureSets];
        Set<String>[] quantPeptideSets = (Set<String>[]) new Set[numFeatureSets];

        double[][] overlaps = new double[numFeatureSets][numFeatureSets];

        Set<String> allPeptides = new HashSet<String>();
        Set<String> allQuantPeptides = new HashSet<String>();

        for (int i=0; i<numFeatureSets; i++)
        {
            featureSetsToHandle[i] = featureSetsToHandle[i].filter(sel);

            Set<String> currentSetPeptides = new HashSet<String>();
            Set<String> currentSetQuantPeptides = new HashSet<String>();

            for (Feature feature : featureSetsToHandle[i].getFeatures())
            {
                List<String> peptideList = MS2ExtraInfoDef.getPeptideList(feature);
                if (peptideList != null)
                {
                    for (String peptide : peptideList)
                    {
                        currentSetPeptides.add(peptide);
                        allPeptides.add(peptide);
                        if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                        {
                            currentSetQuantPeptides.add(peptide);
                            allQuantPeptides.add(peptide);
                        }
                    }
                }
            }
            peptideSets[i] = currentSetPeptides;
            quantPeptideSets[i] = currentSetQuantPeptides;
        }

        for (int i=0; i<numFeatureSets; i++)
        {
            for (int j=0; j<numFeatureSets; j++)
            {
                if (i==j)
                {
                    overlaps[i][j] = peptideSets[i].size();
                    continue;
                }
                if (overlaps[j][i] > 0)
                {
                    overlaps[i][j] = overlaps[j][i];
                    continue;
                }
                int overlap = 0;
                for (String peptide : peptideSets[i])
                {
                    if (peptideSets[j].contains(peptide))
                    {
                        _log.debug(peptide);
                        overlap++;
                    }
                }
                overlaps[i][j] = overlap;
            }
        }

        double[] peptideCounts = new double[numFeatureSets];
        double[] quantPeptideCounts = new double[numFeatureSets];

        for (int i=0; i<featureSetsToHandle.length; i++)
        {
            ApplicationContext.infoMessage("Set " + (i+1) + " (" + featureSetsToHandle[i].getSourceFile() +
                    ") unique peptides: " + peptideSets[i].size() + ", unique quant: " + quantPeptideSets[i].size());
            peptideCounts[i] = peptideSets[i].size();
            quantPeptideCounts[i] = quantPeptideSets[i].size();

        }

        Set<String> pepsIn2ButNot1 = new HashSet<String>();
        Set<String> quantPepsIn2ButNot1 = new HashSet<String>();


        if (peptideSets.length == 2)
        {
            for (String peptide2 : peptideSets[1])
                if (!peptideSets[0].contains(peptide2))
                    pepsIn2ButNot1.add(peptide2);
            for (String quantPeptide2 : quantPeptideSets[1])
                if (!quantPeptideSets[0].contains(quantPeptide2))
                    quantPepsIn2ButNot1.add(quantPeptide2);

            int numIn2butnot1 = pepsIn2ButNot1.size();
            Map<String, Float> in2ButNot1PeptideProbabilityMap = new HashMap<String,Float>();
            for (Feature feature2 : featureSetsToHandle[1].getFeatures())
            {
                String peptide2 = MS2ExtraInfoDef.getFirstPeptide(feature2);
                if (peptide2 != null && !peptideSets[0].contains(peptide2))
                {
                    Float existingProb = in2ButNot1PeptideProbabilityMap.get(peptide2);
                    if (existingProb == null)
                        existingProb = 0f;
                    float newProb = (float) MS2ExtraInfoDef.getPeptideProphet(feature2);
                    if (newProb > existingProb)
                        in2ButNot1PeptideProbabilityMap.put(peptide2, newProb);
                }
            }

            List<Float> in2ButNot1Probabilities = new ArrayList<Float>();
            in2ButNot1Probabilities.addAll(in2ButNot1PeptideProbabilityMap.values());
            System.err.println("Peptides in 2 but not 1: " + numIn2butnot1 +
                    " (" + Rounder.round(100 * (double) numIn2butnot1 / (double) peptideSets[0].size(),1) +
                    "% of original set 1 size)");
            float meanprob =  (float) BasicStatistics.mean(in2ButNot1Probabilities);
            System.err.println("\tmean probability: " +
                    BasicStatistics.mean(in2ButNot1Probabilities));
            float expectedtrue = (in2ButNot1Probabilities.size() * meanprob);
            System.err.println("\texpected true: " + expectedtrue + " (" +
                    Rounder.round(100 * expectedtrue / peptideSets[0].size(),1) + "% of set 1 peptides)");
        }



//        PaintScale paintScale = PanelWithHeatMap.createPaintScale(0, 100, Color.BLUE, Color.RED);
        StringBuffer htmlStringBuf = new StringBuffer("<table border=\"1\"><tr><td></td>");
        for (int i=0; i<numFeatureSets; i++)
        {
            htmlStringBuf.append("<th>" + (i+1) + "</th>");
        }
        htmlStringBuf.append("</tr>\n");
        for (int i=0; i<numFeatureSets; i++)
        {
            htmlStringBuf.append("<tr><th>" + (i+1) + "</th>");
            StringBuffer rowText = new StringBuffer();
            for (int j=0; j<numFeatureSets; j++)
            {
                if (j>0)
                    rowText.append("\t");
                int percentOverlap = 0;
                if (peptideSets[i].size() > 0)
//                    percentOverlap = (100 * overlaps[i][j] / (peptideSets[i].size() + peptideSets[j].size() - overlaps[i][j]));
                percentOverlap = 100 * (int) overlaps[i][j] / peptideSets[i].size();
                rowText.append((int) overlaps[i][j] + " (" + percentOverlap + "%)");
int red = percentOverlap * 256 / 100 - 1;
int blue = (100 - percentOverlap) * 256 / 100 - 1;
//System.err.println(percentOverlap + ", " + red + ", " + blue + "; " + Integer.toString(red, 16) + ", " + Integer.toString(blue, 16));                
String colorAsHex = Integer.toString(red, 16) + "00" + Integer.toString(blue, 16);
                String text = (int) overlaps[i][j] + " (" + percentOverlap + "%)";
//                if (i>j)
//                {
//                    colorAsHex = "FFFFFF";
//                    text = "&nbsp;";
//                }

                htmlStringBuf.append("<td bgcolor=\"#" + colorAsHex + "\">" + text + "</td>");
            }
            ApplicationContext.infoMessage(rowText.toString());
            htmlStringBuf.append("</tr>\n");
        }
        htmlStringBuf.append("</table>");
//System.err.println("**" + Integer.toString(256, 16));
        if (showCharts)
        {
            PanelWithRPerspectivePlot plot3D = new PanelWithRPerspectivePlot();
            double[] plot3DAxisPoints = new double[numFeatureSets];
            for (int i=0; i<numFeatureSets; i++)
                plot3DAxisPoints[i] = i;
            plot3D.setRotationAngle(45);
            plot3D.plot(plot3DAxisPoints, plot3DAxisPoints, overlaps);
            plot3D.displayInTab();
//            try
//            {
//                BrowserController.navigateOrPanelTempFileWithContents(htmlStringBuf.toString(), "tablehtml.html", this);
//            }
//            catch (Exception e)
//            {
//                throw new CommandLineModuleExecutionException(e);
//            }
        }

        Set<String> commonAcrossAll = new HashSet<String>();
        List<Float> numbersOfFilesPeptidesOccurIn = new ArrayList<Float>();
        List<Float> numbersOfFilesQuantPeptidesOccurIn = new ArrayList<Float>();

        int numSingleFilePeptides = 0;
        for (String peptide : allPeptides)
        {
            int numSetsThisPeptide = 0;
            int numQuantSetsThisPeptide = 0;

            for (int i = 0; i < numFeatureSets; i++)
            {
                if (peptideSets[i].contains(peptide))
                {
                    numSetsThisPeptide++;
                }
                if (quantPeptideSets[i].contains(peptide))
                {
                    numQuantSetsThisPeptide++;
                }
            }
            if (numSetsThisPeptide == 1)
                numSingleFilePeptides++;
            boolean inAll = (numSetsThisPeptide == numFeatureSets);

            numbersOfFilesPeptidesOccurIn.add((float) numSetsThisPeptide);
            if (numQuantSetsThisPeptide>0)
                numbersOfFilesQuantPeptidesOccurIn.add((float) numQuantSetsThisPeptide);

            if (inAll)
               commonAcrossAll.add(peptide);

        }




        ApplicationContext.infoMessage("Total distinct peptides in all runs: " + allPeptides.size());
        ApplicationContext.infoMessage("Total distinct QUANTIFIED peptides in all runs: " + allQuantPeptides.size());


        System.err.println("Average number of peptides per set: " + BasicStatistics.mean(peptideCounts));
        System.err.println("Average number of quant peptides per set: " + BasicStatistics.mean(quantPeptideCounts));

        System.err.println("Average number of sets per peptide: " + BasicStatistics.mean(numbersOfFilesPeptidesOccurIn) +
                      ", median: " + BasicStatistics.median(numbersOfFilesPeptidesOccurIn) + ", sum: " +
                BasicStatistics.sum(numbersOfFilesPeptidesOccurIn));
        if (allPeptides.size() > 0)
            System.err.println("Number of single-file peptides: " + numSingleFilePeptides +
                    "(" +(100 * numSingleFilePeptides / allPeptides.size()) + "%)");
        System.err.println("Average number of sets per QUANTITATED peptide: " + BasicStatistics.mean(numbersOfFilesQuantPeptidesOccurIn) +
                      ", median: " + BasicStatistics.median(numbersOfFilesQuantPeptidesOccurIn) + ", sum: " +
                      BasicStatistics.sum(numbersOfFilesQuantPeptidesOccurIn));

        
        ApplicationContext.infoMessage("Common peptides across all runs: " + commonAcrossAll.size());
        if (shouldListAllCommon)
            for (String commonPep : commonAcrossAll)
                ApplicationContext.infoMessage("\t" + commonPep);
        if (_log.isDebugEnabled())
        {
            for (String peptide : commonAcrossAll)
                System.err.println("\t" + peptide);
        }
        double averageConcordanceAllPairs = 0;
        for (int i=0; i<overlaps.length; i++)
        {
            for (int j=0; j<overlaps.length; j++)
            {
                if (i != j)
                    averageConcordanceAllPairs += overlaps[i][j];
            }
        }
        averageConcordanceAllPairs /= (numFeatureSets * numFeatureSets - numFeatureSets);
        System.err.println("Average concordance of all non-identity pairs of featuresets: " + averageConcordanceAllPairs);

        if (numFeatureSets == 2 && showCharts)
        {
            PanelWithVenn vennPanel =
                    new PanelWithVenn(peptideSets[0].size(), peptideSets[1].size(), commonAcrossAll.size());
            vennPanel.setSize(new Dimension(400,300));
            JDialog vennDialog = new JDialog();
            vennDialog.add(vennPanel);
            vennDialog.setSize(vennPanel.getSize());
            vennDialog.setVisible(true);
        }
        return new Pair<Set<String>,Set<String>>(pepsIn2ButNot1,quantPepsIn2ButNot1);
//scaff
//double minDistance = Double.MAX_VALUE;
//Pair<String, String> minDistancePeptides = null;
//for (String commonPeptide1 : commonAcrossAll)
//{
//    double t1 = -1;
//    for (Feature feature : featureSets[0].getFeatures())
//        if (commonPeptide1.equals(MS2ExtraInfoDef.getFirstPeptide(feature)))
//        {
//            t1 = feature.getTime();
//            break;
//        }
//    if (t1 <= 0) continue;
//    for (String commonPeptide2 : commonAcrossAll)
//    {
//        if (commonPeptide1.equals(commonPeptide2))
//            continue;
//        double t2 = -1;
//        for (Feature feature : featureSets[0].getFeatures())
//            if (commonPeptide2.equals(MS2ExtraInfoDef.getFirstPeptide(feature)))
//            {
//                t2 = feature.getTime();
//                break;
//            }
//        if (t2 <= 0) continue;
//
//        double distance = Math.abs(t1 - t2);
//
//        double calcHDistance = Math.abs(AmtUtilities.calculateNormalizedHydrophobicity(commonPeptide1) - AmtUtilities.calculateNormalizedHydrophobicity(commonPeptide2));
//        if (distance < minDistance && calcHDistance > 0.3)
//        {
//            minDistance = distance;
//            minDistancePeptides = new Pair<String,String>(commonPeptide1, commonPeptide2);
//System.err.println("minDistance: " + distance + ", calcHDistance: " + calcHDistance + ", hs: " + AmtUtilities.calculateNormalizedHydrophobicity(commonPeptide1) + ", " + AmtUtilities.calculateNormalizedHydrophobicity(commonPeptide2));
//        }
//    }
//}
//System.err.println("Min distance peptides: " + minDistancePeptides + " -- distance: " + minDistance);
//end scaff

        
    }



}
