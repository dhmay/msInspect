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
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.APMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;


/**
 * Commandline module for FDR calculation and filtering
 *
 * This allows the user to take a feature file (which can be a pepXML file) and do an FDR calculation based on
 * reverse database hits (rev_ prefix), using some single peptide quality score (PeptideProphet probability, or
 * any arbitrary search_score).  The results are then filtered according to some maximum FDR threshold and output
 * to a new file, with all filtered results having a PeptideProphet probability of some specified value. 
 */
public class CalculateFDRCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalculateFDRCLM.class);

    public static final int MIN_FEATURES_FOR_FDR_CALC = 50;

    protected File[] featureFiles;
    protected String searchScoreName = null;
    protected File outFile = null;
    protected File outDir = null;

    protected boolean higherIsBetter=false;
    protected boolean showCharts=false;

    protected float maxFDRToKeep = .05f;
    protected float passingFeaturePeptideProphetValue = .95f;

    protected boolean setPeptideProphet1MinusFDR = false;

    protected float targetDecoyDBSizeRatio = 1.0f;

    protected int scoreType = MODE_PPROPHET;

    protected int outFormat = OUT_FORMAT_INPUT;

    protected boolean byCharge = true;
    protected boolean calcAllRunsFDRTogether = true;

    public static final String DEFAULT_REV_PROTEIN_PREFIX = "rev_";
    protected String reverseProteinPrefix = DEFAULT_REV_PROTEIN_PREFIX;

    protected File saveChartsDir = null;


    public static final String DEFAULT_SEARCH_SCORE_NAME = "expect";



    protected final static String[] scoreTypeStrings =
            {
                    "pprophet",
                    "searchscore"
            };

    protected final static String[] scoreTypeExplanations =
            {
                    "Use PeptideProphet probability",
                    "Use a search score (name must be provided)"
            };

    protected static final int MODE_PPROPHET = 0;
    protected static final int MODE_SEARCH_SCORE = 1;


    protected final static String[] outFormatStrings =
            {
                    "pepxml",
                    "msinspect",
                    "apml",
                    "input"
            };

    protected final static String[] outFormatExplanations =
            {
                    "PepXML",
                    "msInspect .tsv",
                    "APML",
                    "Same as input format (PepXML or msInspect)"
            };

    protected static final int OUT_FORMAT_PEPXML = 0;
    protected static final int OUT_FORMAT_MSINSPECT = 1;
    protected static final int OUT_FORMAT_APML = 2;
    protected static final int OUT_FORMAT_INPUT = 3;


    protected boolean shouldWriteOutput = false;


    public CalculateFDRCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "calcfdr";
        mShortDescription = "Use reverse database hits to calculate False Discovery Rate and filter.";
        mHelpMessage = "Use reverse database hits to calculate the 'target' False Discovery Rate (#decoy/#target) " +
                "for each peptide identification, based on some score (either PeptideProphet probability, or another " +
                "search_score specified explicitly).  Calculate q-values and use those instead of raw FDR.  " +
                "Filter the results " +
                "using FDR cutoff 'maxfdr' and assign all results the same arbitrary PeptideProphet " +
                "probability (argument 'pprophetvalue')";

        CommandLineArgumentDefinition[] argDefs =
                {
                       new EnumeratedValuesArgumentDefinition("scoretype", false, scoreTypeStrings,
                                                          scoreTypeExplanations, "searchscore"),
                        createUnnamedSeriesFileArgumentDefinition(true, "Input feature file(s)"),
                        new StringArgumentDefinition("searchscorename", false,
                                "Name of the search score to use (for 'searchscore' mode)",
                                                       DEFAULT_SEARCH_SCORE_NAME),
                        new FileToWriteArgumentDefinition("out", false, "Output file (for single file processing)"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
                        new BooleanArgumentDefinition("higherisbetter", false,
                                "Is a higher value better, for this score (for 'searchscore' mode)?", higherIsBetter),
                        new BooleanArgumentDefinition("showcharts", false,
                                "Plot an ROC curve?", showCharts),
                        new DecimalArgumentDefinition("pprophetvalue", false,
                                "Set the PeptideProphet score of every passing feature to this value",
                                passingFeaturePeptideProphetValue),
                        new DecimalArgumentDefinition("maxfdr", false,
                                "Maximum FDR to keep in output file",
                                maxFDRToKeep),
                        new EnumeratedValuesArgumentDefinition("outformat", false, outFormatStrings, outFormatExplanations,
                                                           "input"),
                        new DecimalArgumentDefinition("targetdecoydbsizeratio", false,
                                "Ratio of the number of peptides in the target search database to the number " +
                                "of peptides in the decoy search database.",
                                targetDecoyDBSizeRatio),
                        new BooleanArgumentDefinition("bycharge", false,
                                "Calculate FDR separately by charge state? Charge states with too few " +
                                        "identifications will be dropped", byCharge),
                        new BooleanArgumentDefinition("together", false,
                                "Calcualte FDR on all fractions together?  Otherwise, calculate FDR separately for each" +
                                "run", calcAllRunsFDRTogether),
                        new DirectoryToReadArgumentDefinition("savechartsdir", false, "Directory to save charts to"),
                        new BooleanArgumentDefinition("setpprophet1minusfdr", false,
                                "Set PeptideProphet score to 1 - FDR?", setPeptideProphet1MinusFDR),
                        new StringArgumentDefinition("revproteinprefix", false,
                                "Prefix in the FASTA file for proteins that are reversed sequences",
                                DEFAULT_REV_PROTEIN_PREFIX)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        scoreType = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("scoretype")).getIndexForArgumentValue(
                getStringArgumentValue("scoretype"));

        featureFiles = this.getUnnamedSeriesFileArgumentValues();

        higherIsBetter = getBooleanArgumentValue("higherisbetter");
        
        reverseProteinPrefix = getStringArgumentValue("revproteinprefix");

        switch(scoreType)
        {
            case MODE_SEARCH_SCORE:
                searchScoreName = getStringArgumentValue("searchscorename");
                ApplicationContext.infoMessage("Building FDR from search score '" + searchScoreName + "'");
                break;
            case MODE_PPROPHET:
                assertArgumentAbsent("searchscorename","mode");
                assertArgumentAbsent("higherisbetter");
                ApplicationContext.infoMessage("Building FDR from PeptideProphet probability");                
                higherIsBetter = true;
                break;
        }

        maxFDRToKeep = getFloatArgumentValue("maxfdr");
        ApplicationContext.infoMessage("Max FDR to keep: " + maxFDRToKeep);
        passingFeaturePeptideProphetValue = getFloatArgumentValue("pprophetvalue");

        setPeptideProphet1MinusFDR = getBooleanArgumentValue("setpprophet1minusfdr");
        if (hasArgumentValue("pprophetvalue") && setPeptideProphet1MinusFDR)
        {
            throw new ArgumentValidationException("'pprophetvalue' argument can't be specified if " +
                    "'setpprophet1minusfdr' argument is 'true'");
        }

        targetDecoyDBSizeRatio = getFloatArgumentValue("targetdecoydbsizeratio");

        showCharts = getBooleanArgumentValue("showcharts");
        saveChartsDir = getFileArgumentValue("savechartsdir");



        byCharge = getBooleanArgumentValue("bycharge");
        calcAllRunsFDRTogether = getBooleanArgumentValue("together");



        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        shouldWriteOutput = (outFile != null || outDir != null);

        if (featureFiles.length > 1 && hasArgumentValue("out"))
        {
            throw new ArgumentValidationException("ERROR: Multiple inputs, one output.");
        }
        if (outFile != null)
            assertArgumentAbsent("outdir");


        outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(
                getStringArgumentValue("outformat"));

    }

    protected boolean isInputPepXML(File featureFile)
            throws IOException
    {
        boolean inputIsPepXML = false;
        if (PepXMLFeatureFileHandler.getSingletonInstance().canHandleFile(featureFile))
            inputIsPepXML = true;
        return inputIsPepXML;
    }

    protected List<FeatureSet> loadFeatureSetsFromFile(File featureFile)
            throws Exception
    {

        List<FeatureSet> featureSetsInFile = new ArrayList<FeatureSet>();

        if (isInputPepXML(featureFile))
            featureSetsInFile = PepXMLFeatureFileHandler.getSingletonInstance().loadAllFeatureSets(featureFile);
        else
            featureSetsInFile.add(new FeatureSet(featureFile));
        return featureSetsInFile;
    }

    protected int calcOutputFormat(File featureFile)
            throws IOException
    {
        int thisFileOutFormat = outFormat;
        if (thisFileOutFormat == OUT_FORMAT_INPUT)
        {
            thisFileOutFormat = isInputPepXML(featureFile) ? OUT_FORMAT_PEPXML :
                    (APMLFeatureFileHandler.getSingletonInstance().canHandleFile(featureFile) ?
                            OUT_FORMAT_APML :  OUT_FORMAT_MSINSPECT);
        }
        return thisFileOutFormat;
    }

    protected Map<Integer, List<Feature>> splitFeaturesByCharge(Feature[] features)
    {
        Map<Integer, List<Feature>> result = new HashMap<Integer, List<Feature>>();
        for (Feature feature : features)
        {
            int charge = feature.getCharge();
            List<Feature> featuresThisCharge = result.get(charge);
            if (featuresThisCharge == null)
            {
                featuresThisCharge = new ArrayList<Feature>();
                result.put(charge, featuresThisCharge);
            }
            featuresThisCharge.add(feature);
        }
        return result;
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (calcAllRunsFDRTogether)
        {
            ApplicationContext.infoMessage("Calculating FDR of all runs together");
            List<Feature> allFeaturesAllRuns = new ArrayList<Feature>();
            Map<File, List<FeatureSet>> fileFeatureSetListMap = new HashMap<File, List<FeatureSet>>();

            for (File featureFile : featureFiles)
            {
                try
                {
                    List<FeatureSet> featureSetsInFile = loadFeatureSetsFromFile(featureFile);

                    for (FeatureSet featureSet : featureSetsInFile)
                    {
                        for (Feature feature : featureSet.getFeatures())
                            allFeaturesAllRuns.add(feature);
                    }
                    fileFeatureSetListMap.put(featureFile, featureSetsInFile);

                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Failed to load feature file " +
                            featureFile.getAbsolutePath());
                }
            }
            if (byCharge)
            {
                Map<Integer, List<Feature>> featuresByCharge =
                        splitFeaturesByCharge(allFeaturesAllRuns.toArray(new Feature[allFeaturesAllRuns.size()]));
                for (int charge=0; charge<10; charge++)
                {
                    if (featuresByCharge.containsKey(charge))
                    {
                        ApplicationContext.infoMessage("Calculating FDR for charge " + charge);
                        List<Feature> featuresThisCharge = featuresByCharge.get(charge);
                        try
                        {
                            calcFDROnFeatures(featuresThisCharge.toArray(new Feature[featuresThisCharge.size()]),
                                    "_charge"+charge);
                        }
                        catch (CommandLineModuleExecutionException e)
                        {
                            ApplicationContext.infoMessage("FDR calc for charge " + charge +
                                    " failed, those features will not get FDR assignments");
                        }
                    }
                }
            }
            else
            {
                calcFDROnFeatures(allFeaturesAllRuns.toArray(new Feature[allFeaturesAllRuns.size()]), "");
            }

            if (shouldWriteOutput)
            {
                for (File featureFile : fileFeatureSetListMap.keySet())
                {
                    File outputFile = outFile;
                    if (outFile == null)
                        outputFile = new File(outDir, featureFile.getName());
                    for (FeatureSet featureSet : fileFeatureSetListMap.get(featureFile))
                    {
                        List<Feature> featuresToKeep = new ArrayList<Feature>();
                        for (Feature feature : featureSet.getFeatures())
                        {
                            //a hit with any forward-database proteins is considered a forward hit
                            boolean foundForwardProtein = false;
                            for (String protein : MS2ExtraInfoDef.getProteinList(feature))
                                if (!protein.contains(reverseProteinPrefix))
                                {
                                    foundForwardProtein = true;
                                    break;
                                }
                            if (foundForwardProtein &&
                                    MS2ExtraInfoDef.hasFalseDiscoveryRate(feature) &&
                                    MS2ExtraInfoDef.getFalseDiscoveryRate(feature) <= maxFDRToKeep)
                                featuresToKeep.add(feature);
                        }

                        featureSet.setFeatures(featuresToKeep.toArray(new Feature[featuresToKeep.size()]));
                    }
                    try
                    {
                        writeFile(fileFeatureSetListMap.get(featureFile), outputFile, calcOutputFormat(featureFile));
                    }
                    catch(IOException e)
                    {
                        throw new CommandLineModuleExecutionException("Failure writing output file " +
                                outputFile.getAbsolutePath(),e);
                    }
                }
            }
        }
        else
        {
            ApplicationContext.infoMessage("Calculating FDR of each run separately");
            for (File featureFile : featureFiles)
            {
                try
                {
                    List<FeatureSet> featureSetsInFile = loadFeatureSetsFromFile(featureFile);

                    for (FeatureSet featureSet : featureSetsInFile)
                    {
                        featureSet.setFeatures(calcFDROnFeatures(featureSet.getFeatures(), "_" + featureFile.getName()));
                    }

                    if (shouldWriteOutput)
                    {
                        File outputFile = outFile;
                        if (outFile == null)
                            outputFile = new File(outDir, featureFile.getName());
                        writeFile(featureSetsInFile, outputFile, calcOutputFormat(featureFile));
                    }
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Failed to handle file " +
                            featureFile.getAbsolutePath(),e);
                }
            }
        }

    }

    protected void writeFile(List<FeatureSet> featureSets, File outputFile, int thisFileOutFormat)
            throws CommandLineModuleExecutionException
    {                                  
        try
        {
            switch (thisFileOutFormat)
            {
                case OUT_FORMAT_PEPXML:
                    PepXMLFeatureFileHandler.getSingletonInstance().saveFeatureSets(featureSets, outputFile);
                    break;
                default:
                    if (featureSets.size() > 1)
                        throw new CommandLineModuleExecutionException("Can't save multiple featuresets to " +
                                "msInspect .tsv file " + outputFile.getAbsolutePath());
                    featureSets.get(0).save(outputFile);
                    break;
            }
            ApplicationContext.infoMessage("Saved file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failure writing feature file " +
                    outputFile.getAbsolutePath(),e);
        }
    }


    /**
     *
     * @param features
     * @return forward features passing FDR cutoff
     * @throws CommandLineModuleExecutionException
     */
    protected Feature[] calcFDROnFeatures(Feature[] features, String chartNameSuffix)
            throws CommandLineModuleExecutionException
    {
        Comparator<Feature> comparatorDescendingGoodness =
                new SearchScoreComparator(higherIsBetter);

        ApplicationContext.infoMessage("Total features: " + features.length);
        if (features.length < MIN_FEATURES_FOR_FDR_CALC)
        {
            throw new CommandLineModuleExecutionException("Too few features to calculate FDR");
        }
        Feature[] sortedFeaturesDescGoodness = features;
        Arrays.sort(sortedFeaturesDescGoodness, comparatorDescendingGoodness);

        int numForwardHits = 0;
        int numReverseHits = 0;

        double cutoffSearchScore = 0;

        List<Feature> forwardFeaturesList = new ArrayList<Feature>();

        double lowestScore = Double.MAX_VALUE;
        double highestScore = Double.MIN_VALUE;

        float[] fdrs = new float[sortedFeaturesDescGoodness.length];

        for (int i=0; i<sortedFeaturesDescGoodness.length; i++)
        {
            Feature feature = sortedFeaturesDescGoodness[i];
            double thisFeatureScore = getSearchScoreValue(feature);
            if (thisFeatureScore < lowestScore)
                lowestScore = thisFeatureScore;
            if (thisFeatureScore > highestScore)
                highestScore = thisFeatureScore;

            //a hit with any forward-database proteins is considered a forward hit
            boolean foundForwardProtein = false;
            for (String protein : MS2ExtraInfoDef.getProteinList(feature))
                if (!protein.contains(reverseProteinPrefix))
                    foundForwardProtein = true;

            if (foundForwardProtein)
            {
                numForwardHits++;
                forwardFeaturesList.add(feature);
            }
            else
                numReverseHits++;

            //Calculate the "target" false discovery rate: the expected proportion of target hits that are bad.
            //Number of reverse hits must be scaled by the ratio of database sizes.  Ideally this ratio is the
            //ratio of the number of /peptides/ in the target vs. decoy databases, but number of /proteins/
            //is a reasonable approximation that's good enough for government work.
            //If no forward hits, FDR = 1, to avoid division by 0.
            double fdr = 1.0;
            if (numForwardHits > 0)
                fdr = (double) numReverseHits * targetDecoyDBSizeRatio / (double) numForwardHits;
//            fdr = Rounder.round(fdr,3);
            fdrs[i] = (float) fdr;

            MS2ExtraInfoDef.setFalseDiscoveryRate(feature, (float) fdr);
        }

        if (numReverseHits == 0)
            ApplicationContext.infoMessage("WARNING: no reverse proteins found in search results!  Please be " +
                    "sure that you ran your search against a FASTA database in which reversed proteins have " +
                    "identifiers beginning with '" + reverseProteinPrefix + "'.  If your database uses " +
                    "a different prefix, specify it with the 'revproteinprefix' argument.");

        float[] qvals = new float[sortedFeaturesDescGoodness.length];

        boolean foundIt = false;
        for (int i=sortedFeaturesDescGoodness.length-1; i>=0; i--)
        {
            double worstPossibleValue = i==qvals.length-1 ? Float.MAX_VALUE : qvals[i+1];
            qvals[i] = (float) Math.min(worstPossibleValue, fdrs[i]);

            if (!foundIt && qvals[i] <= maxFDRToKeep)
            {
                foundIt = true;
                cutoffSearchScore = getSearchScoreValue(sortedFeaturesDescGoodness[i]);
            }
//if (qvals[i] < fdrs[i]) System.err.println(getSearchScoreValue(sortedFeaturesDescGoodness[i]) + ", " + fdrs[i] + ", " + qvals[i]);

//if (i<qvals.length-1 && qval > qvals[i+1]) System.err.println("Bad! " + qval + ", " + qvals[i+1]);
//if (i<qvals.length-1 && getSearchScoreValue(sortedFeaturesDesc[i]) < getSearchScoreValue(sortedFeaturesDesc[i+1])) System.err.println("Bad! " + getSearchScoreValue(sortedFeaturesDesc[i]) + ", " + getSearchScoreValue(sortedFeaturesDesc[i+1]));
        }

        ApplicationContext.infoMessage("Cutoff value: " + cutoffSearchScore);


        List<Feature> passingForwardFeatures = new ArrayList<Feature>();
        List<Feature> passingFeatures = new ArrayList<Feature>();

        for (Feature feature : forwardFeaturesList)
        {
            double score = getSearchScoreValue(feature);
            if ( (higherIsBetter && score >= cutoffSearchScore) ||
                 (!higherIsBetter && score <= cutoffSearchScore))
                passingForwardFeatures.add(feature);
        }

        for (Feature feature : sortedFeaturesDescGoodness)
        {
            double score = getSearchScoreValue(feature);
            if ( (higherIsBetter && score >= cutoffSearchScore) ||
                 (!higherIsBetter && score <= cutoffSearchScore))
                passingFeatures.add(feature);
        }

        for (Feature passingFeature : passingForwardFeatures)
        {
            float featurePepProphetValue = passingFeaturePeptideProphetValue;
            if (setPeptideProphet1MinusFDR)
                featurePepProphetValue = 1 - MS2ExtraInfoDef.getFalseDiscoveryRate(passingFeature);
            MS2ExtraInfoDef.setPeptideProphet(passingFeature, featurePepProphetValue);
            MS2ExtraInfoDef.setAllNttProb(passingFeature,
                    "(" + featurePepProphetValue + "," + featurePepProphetValue + "," + featurePepProphetValue + ")");
        }

        if (showCharts || saveChartsDir != null)
        {
            float[] searchScores = new float[fdrs.length];
            for (int i=0; i<sortedFeaturesDescGoodness.length; i++)
                searchScores[i] = (float) getSearchScoreValue(sortedFeaturesDescGoodness[i]);
            PanelWithLineChart rocPanel = new PanelWithLineChart(searchScores, fdrs, "FDR (y) vs Score (x)");
            rocPanel.setAxisLabels("Score","FDR (q-value)");
            rocPanel.addData(searchScores, qvals, "q-value vs. score");
            if (showCharts)
                rocPanel.displayInTab();
            if (saveChartsDir != null)
            try
            {
                File chartFile = new File(saveChartsDir,"FDR_vs_score" + chartNameSuffix + ".png");
                rocPanel.saveChartToImageFile(chartFile);
                ApplicationContext.infoMessage("Saved chart file " + chartFile);            

            }
            catch (IOException e)
            {
                ApplicationContext.errorMessage("Failed to save chart",e);
            }
ApplicationContext.infoMessage("WARNING! Due to high precision of score values and FDR values, this chart is messed up");            

//            PanelWithLineChart rocPanel2 = new PanelWithLineChart(fdrs, qvals, "FDR vs q-score");
//            rocPanel2.setAxisLabels("Score","FDR (q-value)");
//            rocPanel2.displayInTab();
        }

        ApplicationContext.infoMessage("Forward features passing cutoff: " + passingForwardFeatures.size());
        return passingForwardFeatures.toArray(new Feature[passingForwardFeatures.size()]);
    }

    public void writeFeatureSet(FeatureSet featureSet, File outputFile, int outputFormat)
            throws CommandLineModuleExecutionException
    {

        try
        {
            switch (outputFormat)
            {
                case OUT_FORMAT_PEPXML:
                    featureSet.savePepXml(outputFile);
                    break;
                case OUT_FORMAT_MSINSPECT:
                    featureSet.save(outputFile);
                    break;
            }
            ApplicationContext.setMessage("Saved output file " + outputFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }


        if (featureSet.getFeatures().length == 0)
        {
            ApplicationContext.infoMessage("No number of these features gives the requested FDR! Wrote empty file.");
            return;
        }

        ApplicationContext.infoMessage("Total # of FORWARD hits passing FDR: " + featureSet.getFeatures().length);
    }

    protected double getSearchScoreValue(Feature feature)
    {
        if (scoreType == MODE_PPROPHET)
            return MS2ExtraInfoDef.getPeptideProphet(feature);
        else
            {
                double val = 0;
                String scoreString = MS2ExtraInfoDef.getSearchScore(feature, searchScoreName);
                try {
                    val = Double.parseDouble(scoreString);

                } catch (Exception e) {
                    val = Float.parseFloat(scoreString);
                }
                return val;
            }
    }

    public class SearchScoreComparator implements Comparator<Feature>
    {
        protected boolean descendingOrder = true;

        public SearchScoreComparator(boolean descendingOrder)
        {
            this.descendingOrder = descendingOrder;
        }

        public int compare(Feature o1, Feature o2)
        {
            double o1Value = getSearchScoreValue(o1);
            double o2Value = getSearchScoreValue(o2);

            if (o1Value == o2Value)
                return 0;

            if (descendingOrder)
                return o1Value < o2Value ? 1 : -1;
            else
                return o1Value > o2Value ? 1 : -1;

        }
    }

}
