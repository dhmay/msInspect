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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.filehandler.APMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;


/**
 * Commandline module for FDR calculation and filtering
 *
 * This allows the user to take a feature file (which can be a pepXML file) and do an FDR calculation based on
 * reverse database hits (rev_ prefix), using some single peptide quality score (PeptideProphet probability, or
 * any arbitrary search_score).  The results are then filtered according to some maximum FDR threshold and output
 * to a new file, with all filtered results having a PeptideProphet probability of some specified value. 
 */
public class CalculateFDRCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalculateFDRCLM.class);

    protected File[] featureFiles;
    protected String searchScoreName = null;
    protected File outFile = null;
    protected File outDir = null;

    protected boolean higherIsBetter=false;
    protected boolean showCharts=false;

    protected float maxFDRToKeep = .05f;
    protected float passingFeaturePeptideProphetValue = .95f;

    protected float targetDecoyDBSizeRatio = 1.0f;

    protected int scoreType = MODE_PPROPHET;

    protected int outFormat = OUT_FORMAT_INPUT;

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
                       createEnumeratedArgumentDefinition("scoretype", false, scoreTypeStrings,
                                                          scoreTypeExplanations, "searchscore"),
                        createUnnamedSeriesArgumentDefinition(ArgumentDefinitionFactory.FILE_TO_READ,true,
                                "Input feature file(s)"),
                        createStringArgumentDefinition("searchscorename", false,
                                "Name of the search score to use (for 'searchscore' mode)",
                                                       DEFAULT_SEARCH_SCORE_NAME),
                        createFileToWriteArgumentDefinition("out", false, "Output file (for single file processing)"),
                        createDirectoryToReadArgumentDefinition("outdir", false, "Output directory"),
                        createBooleanArgumentDefinition("higherisbetter", false,
                                "Is a higher value better, for this score (for 'searchscore' mode)?", higherIsBetter),
                        createBooleanArgumentDefinition("showcharts", false,
                                "Plot an ROC curve?", showCharts),
                        createDecimalArgumentDefinition("pprophetvalue", false,
                                "Set the PeptideProphet score of every passing feature to this value",
                                passingFeaturePeptideProphetValue),
                        createDecimalArgumentDefinition("maxfdr", false,
                                "Maximum FDR to keep in output file",
                                maxFDRToKeep),
                        createEnumeratedArgumentDefinition("outformat", false, outFormatStrings, outFormatExplanations,
                                                           "input"),
                        createDecimalArgumentDefinition("targetdecoydbsizeratio", false,
                                "Ratio of the number of peptides in the target search database to the number " +
                                "of peptides in the decoy search database.",
                                targetDecoyDBSizeRatio),
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
        

        switch(scoreType)
        {
            case MODE_SEARCH_SCORE:
                searchScoreName = getStringArgumentValue("searchscorename");
                break;
            case MODE_PPROPHET:
                assertArgumentAbsent("searchscorename","mode");
                assertArgumentAbsent("higherisbetter");
                higherIsBetter = true;
                break;
        }

        maxFDRToKeep = (float) getDoubleArgumentValue("maxfdr");
        passingFeaturePeptideProphetValue = (float) getDoubleArgumentValue("pprophetvalue");

        targetDecoyDBSizeRatio = (float) getDoubleArgumentValue("targetdecoydbsizeratio");

        showCharts = getBooleanArgumentValue("showCharts");


        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (featureFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }
        if (outFile != null)
            assertArgumentAbsent("outdir");

        outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(
                getStringArgumentValue("outformat"));

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Comparator<Feature> comparatorDescendingGoodness = 
                    new SearchScoreComparator(higherIsBetter);

        for (File featureFile : featureFiles)
        {
            File outputFile = outFile;
            if (outFile == null)
                outputFile = new File(outDir, featureFile.getName());
            handleFile(featureFile, outputFile, comparatorDescendingGoodness);
        }
    }

    protected void handleFile(File featureFile, File outputFile, Comparator<Feature> comparatorDescendingGoodness)
            throws CommandLineModuleExecutionException
    {
        try
        {
            int thisFileOutFormat = outFormat;
            boolean inputIsPepXML = false;
            if (PepXMLFeatureFileHandler.getSingletonInstance().canHandleFile(featureFile))
                inputIsPepXML = true;
            if (thisFileOutFormat == OUT_FORMAT_INPUT)
            {
                thisFileOutFormat = inputIsPepXML ? OUT_FORMAT_PEPXML :
                        (APMLFeatureFileHandler.getSingletonInstance().canHandleFile(featureFile) ?
                                OUT_FORMAT_APML :  OUT_FORMAT_MSINSPECT);
            }

            List<FeatureSet> featureSetsInFile = new ArrayList<FeatureSet>();

            if (inputIsPepXML)
                featureSetsInFile = PepXMLFeatureFileHandler.getSingletonInstance().loadAllFeatureSets(featureFile);
            else
                featureSetsInFile.add(new FeatureSet(featureFile));

            File thisSetOutputFile = outputFile;
            for (FeatureSet featureSet : featureSetsInFile)
            {
                if (featureSetsInFile.size() > 1)
                {
                    thisSetOutputFile = new File(outputFile.getParentFile(),
                            MS2ExtraInfoDef.getFeatureSetBaseName(featureSet) + ".pep.xml");
                }
                handleFeatureSet(featureSet,thisSetOutputFile, comparatorDescendingGoodness, thisFileOutFormat);
            }

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failure loading feature file " +
                    featureFile.getAbsolutePath(),e);
        }

    }


    protected void handleFeatureSet(FeatureSet featureSet, File outputFile,
                                     Comparator<Feature> comparatorDescendingGoodness,
                                     int outputFormat)
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Total peptide IDs: " + featureSet.getFeatures().length);
        Feature[] sortedFeaturesDescGoodness = featureSet.getFeatures();
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
                if (!protein.startsWith("rev_"))
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
        }

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
            MS2ExtraInfoDef.setPeptideProphet(passingFeature, passingFeaturePeptideProphetValue);

        if (showCharts)
        {
            float[] searchScores = new float[fdrs.length];
            for (int i=0; i<sortedFeaturesDescGoodness.length; i++)
                searchScores[i] = (float) getSearchScoreValue(sortedFeaturesDescGoodness[i]);
            PanelWithLineChart rocPanel = new PanelWithLineChart(searchScores, fdrs, "FDR (q-value) vs Score");
            rocPanel.setAxisLabels("Score","FDR (q-value)");
            rocPanel.addData(searchScores, qvals, "q-value vs. score");
            rocPanel.displayInTab();
ApplicationContext.infoMessage("WARNING! Due to high precision of score values and FDR values, this chart is messed up");            

//            PanelWithLineChart rocPanel2 = new PanelWithLineChart(fdrs, qvals, "FDR vs q-score");
//            rocPanel2.setAxisLabels("Score","FDR (q-value)");
//            rocPanel2.displayInTab();
        }


        featureSet.setFeatures(passingForwardFeatures.toArray(new Feature[passingForwardFeatures.size()]));
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


        if (!foundIt)
        {
            ApplicationContext.infoMessage("No number of these features gives the requested FDR! Wrote empty file.");
            return;
        }

        ApplicationContext.infoMessage("Cutoff value: " + cutoffSearchScore);
        ApplicationContext.infoMessage("Total # of hits passing FDR: " + passingFeatures.size());
        ApplicationContext.infoMessage("# of FORWARD hits passing FDR: " + passingForwardFeatures.size());
    }

    protected double getSearchScoreValue(Feature feature)
    {
        if (scoreType == MODE_PPROPHET)
            return MS2ExtraInfoDef.getPeptideProphet(feature);
        else
            return Double.parseDouble(MS2ExtraInfoDef.getSearchScore(feature, searchScoreName));
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
