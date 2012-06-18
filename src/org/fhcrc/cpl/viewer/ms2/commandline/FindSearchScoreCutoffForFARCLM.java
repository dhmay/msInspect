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
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.Comparator;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for feature finding
 */
public class FindSearchScoreCutoffForFARCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FindSearchScoreCutoffForFARCLM.class);

    protected FeatureSet featureSet;
    protected double maxFalseAssignmentRate = 0.05;
    protected String searchScoreName = null;
    protected File outFile = null;
    protected boolean higherIsBetter=true;
    protected boolean plotROC=false;
    protected boolean setPeptideProphetTo1MinusFAR = false;

    protected int mode = MODE_PPROPHET;


    protected final static String[] modeStrings =
            {
                    "pprophet",
                    "searchscore"
            };

    protected final static String[] modeExplanations =
            {
                    "Use PeptideProphet probability",
                    "Use a search score (name must be provided)"
            };

    protected static final int MODE_PPROPHET = 0;
    protected static final int MODE_SEARCH_SCORE = 1;


    public FindSearchScoreCutoffForFARCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "searchscorecutoff";
        mShortDescription = "Use reverse database hits to determine the performance of a search score";
        mHelpMessage = "Use reverse database hits to determine the performance of a search score.  " +
                       "Can be used on PeptideProphet probability or on any search_score.  " +
                       "You specify a maximum false-to-true ratio, and this will give you the score cutoff " +
                       "that gives you that ratio, and optionally an output feature file containing just those " +
                       "features (with false-to-true ratio for each feature in the description)";

        CommandLineArgumentDefinition[] argDefs =
                {
                       new EnumeratedValuesArgumentDefinition("mode",true,modeStrings,
                                                          modeExplanations),
                        createUnnamedFeatureFileArgumentDefinition(true, null),
                        new DecimalArgumentDefinition("maxfar", true, "maximum false assignment rate"),
                        new StringArgumentDefinition("searchscorename", false, "Name of the search score to use"),
                        new FileToWriteArgumentDefinition("out", false, "Output file"),
                        new BooleanArgumentDefinition("higherisbetter", false,
                                "Is a higher value better, for this score?", higherIsBetter),
                        new BooleanArgumentDefinition("plotroc", false,
                                "Plot an ROC curve?", plotROC),
                        new BooleanArgumentDefinition("setpprophet", false,
                                "Set the PeptideProphet score of every passing feature to 1 - the FAR at that score",setPeptideProphetTo1MinusFAR),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        featureSet = this.getFeatureSetArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);


        switch(mode)
        {
            case MODE_SEARCH_SCORE:
                assertArgumentPresent("searchscorename","mode");
                searchScoreName = getStringArgumentValue("searchscorename");
                break;
            case MODE_PPROPHET:
                assertArgumentAbsent("searchscorename","mode");
                assertArgumentAbsent("higherisbetter");
                //for some reason there are sometimes peptideprophet scores of -3
                FeatureSet.FeatureSelector pprophetFeatureSelector = new FeatureSet.FeatureSelector();
                pprophetFeatureSelector.setMinPProphet(0);
                featureSet = featureSet.filter(pprophetFeatureSelector);
                higherIsBetter = true;
                break;
        }

        maxFalseAssignmentRate = getDoubleArgumentValue("maxfar");

        Feature firstFeature = featureSet.getFeatures()[0];

        plotROC = getBooleanArgumentValue("plotroc");

        switch(mode)
        {
            case MODE_SEARCH_SCORE:
                String searchScoreValueString = MS2ExtraInfoDef.getSearchScore(firstFeature, searchScoreName);
                if (searchScoreValueString == null)
                    throw new ArgumentValidationException("Search score named " + searchScoreName + " not found");
                try
                {
                    Double.parseDouble(searchScoreValueString);
                }
                catch (Exception e)
                {
                    throw new ArgumentValidationException("Non-decimal value found for search score " + searchScoreName);
                }
                break;
            case MODE_PPROPHET:
                try
                {
                    MS2ExtraInfoDef.getPeptideProphet(firstFeature);
                }
                catch (Exception e)
                {
                    throw new ArgumentValidationException("No peptideprophet scores found");
                }
                break;
        }



        higherIsBetter = getBooleanArgumentValue("higherisbetter");
        setPeptideProphetTo1MinusFAR = getBooleanArgumentValue("setpprophet");

        outFile = getFileArgumentValue("out");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Comparator<Feature> comparator = null;
        {
            switch(mode)
            {
                case MODE_PPROPHET:
                    comparator =
                            new PeptideProphetComparatorDesc();
                    break;
                case MODE_SEARCH_SCORE:
                    comparator =
                            new SearchScoreComparator(searchScoreName, higherIsBetter);
                    break;
            }
        }
        ApplicationContext.infoMessage("Total features: " + featureSet.getFeatures().length);
        Feature[] sortedFeatures = featureSet.getFeatures();
        Arrays.sort(sortedFeatures, comparator);

        int numForwardHits = 0;
        int numReverseHits = 0;
        boolean foundIt = false;
        int indexOfLastGoodFalseAssignmentRate = -1;
        int numForwardHitsAtLastGoodFAR = 0;
        int numAllHitsAtLastGoodFAR = 0;
        double cutoffSearchScore = 0;

        List<Feature> forwardFeaturesList = new ArrayList<Feature>();


        double lowestScore = Double.MAX_VALUE;
        double highestScore = Double.MIN_VALUE;
        for (int i=0; i<sortedFeatures.length; i++)
        {
            double thisFeatureScore = 0;
            switch(mode)
            {
                case MODE_SEARCH_SCORE:
                    thisFeatureScore = Double.parseDouble(MS2ExtraInfoDef.getSearchScore(sortedFeatures[i], searchScoreName));
                    break;
                case MODE_PPROPHET:
                    thisFeatureScore = MS2ExtraInfoDef.getPeptideProphet(sortedFeatures[i]);
                    break;
            }

            if (thisFeatureScore < lowestScore)
                lowestScore = thisFeatureScore;
            if (thisFeatureScore > highestScore)
                highestScore = thisFeatureScore;

            if (MS2ExtraInfoDef.getFirstProtein(sortedFeatures[i]).startsWith("rev_"))
            {
                numReverseHits++;
//if (i<10) System.err.println("rev\t"+thisFeatureScore);
            }
            else
            {
                forwardFeaturesList.add(sortedFeatures[i]);
                numForwardHits++;
//if (i<10) System.err.println("fwd\t" + thisFeatureScore);

            }


            double falseAssignmentRate = 1;
            if (numForwardHits != 0)
                falseAssignmentRate = (double) numReverseHits / (double) numForwardHits;

            if (falseAssignmentRate <= maxFalseAssignmentRate)
            {
//System.err.println(falseAssignmentRate);                
                foundIt = true;
                indexOfLastGoodFalseAssignmentRate = i;
                numForwardHitsAtLastGoodFAR = numForwardHits;
                numAllHitsAtLastGoodFAR = numReverseHits + numForwardHits;
                cutoffSearchScore = thisFeatureScore;
            }
//            else
//            {
//                System.err.println("No good, score=" + thisFeatureScore + ", bad=" + numReverseHits + ", good=" + numForwardHits + ", ratio=" + falseAssignmentRate);
//            }
        }

        int numBins = 100;
        double[] falseInBins = new double[numBins];
        double[] trueInBins = new double[numBins];

        double binSize = (highestScore - lowestScore) / (double) numBins;
        ApplicationContext.infoMessage("ROC: number of bins: " + numBins + ", bin size: " + binSize);
        ApplicationContext.infoMessage("ROC: min value: " + lowestScore + ", max value: " + highestScore);

        List<Feature> passingForwardFeatures = new ArrayList<Feature>();
        for (int i=0; i<sortedFeatures.length; i++)
        {
            double thisFeatureScore = getSearchScoreValue(sortedFeatures[i]);
            int bin = Math.max(0, (int) ((thisFeatureScore - lowestScore) / binSize) - 1);
            if (MS2ExtraInfoDef.getFirstProtein(sortedFeatures[i]).startsWith("rev_"))
                falseInBins[bin]++;
            else
            {
                passingForwardFeatures.add(sortedFeatures[i]);
                trueInBins[bin]++;
            }
        }

        double[] binValues = new double[numBins];
        double[] sensitivities = new double[numBins];

        int cumulativeFalse = 0;
        int cumulativeTrue = 0;

        for (int i=0; i<numBins; i++)
        {
            int index = i;
            if (higherIsBetter)
                index = (numBins-1) - i;
            binValues[index] = lowestScore + (binSize * index);
            cumulativeFalse += falseInBins[index];
            cumulativeTrue += trueInBins[index];
            sensitivities[index] = (double) cumulativeFalse / (double) cumulativeTrue;
//System.err.println(i + ", " + sensitivities[i] + ", " + cumulativeFalse + ", " + cumulativeTrue);                
        }


        if (plotROC)
        {
            PanelWithLineChart rocPanel = new PanelWithLineChart(binValues, sensitivities, "False/True At or Better Than Score");
            rocPanel.setAxisLabels("Score","False / True");
            new ChartDialog(rocPanel).setVisible(true);
        }

        if (!foundIt)
        {
            ApplicationContext.infoMessage("No number of these features gives the requested FAR!");
            return;
        }

        ApplicationContext.infoMessage("Cutoff value: " + cutoffSearchScore);
        ApplicationContext.infoMessage("Total # of hits passing FAR: " + numAllHitsAtLastGoodFAR);
        ApplicationContext.infoMessage("# of FORWARD hits passing FAR: " + numForwardHitsAtLastGoodFAR);

        if (outFile != null)
        {
            for (Feature passingFeature : passingForwardFeatures)
            {
                double searchScore = getSearchScoreValue(passingFeature);
                int binIndex = 0;
                for (binIndex=0; binIndex<binValues.length; binIndex++)
                {
                    if (binValues[binIndex] >= searchScore)
                        break;
                }
                if (binIndex == binValues.length)
                {
                    binIndex--;
                }
                //1 - the false/true ratio at this score is an approximation of probability.
                //However, this value can be negative, if a local series of negative hits intrudes.
                //Hence, max with 0
                double probability = Math.max(1.0-sensitivities[binIndex],0);
                if (setPeptideProphetTo1MinusFAR)
                    MS2ExtraInfoDef.setPeptideProphet(passingFeature,probability);
                else
                    passingFeature.setDescription("" + probability);
            }


            featureSet.setFeatures(passingForwardFeatures.toArray(new Feature[passingForwardFeatures.size()]));
            try
            {
                featureSet.save(outFile);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

    }

    public static class PeptideProphetComparatorDesc implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            double o1Value = MS2ExtraInfoDef.getPeptideProphet(o1);
            double o2Value = MS2ExtraInfoDef.getPeptideProphet(o2);

            if (o1Value == o2Value)
                return 0;
            return o1Value < o2Value ? 1 : -1;

        }

    }

    protected double getSearchScoreValue(Feature feature)
    {
        if (mode == MODE_PPROPHET)
            return MS2ExtraInfoDef.getPeptideProphet(feature);
        else
            return Double.parseDouble(MS2ExtraInfoDef.getSearchScore(feature, searchScoreName));
    }

    public static class SearchScoreComparator implements Comparator<Feature>
    {
        protected String searchScoreName = null;
        protected boolean descendingOrder = true;

        public SearchScoreComparator(String searchScoreType, boolean descendingOrder)
        {
            searchScoreName = searchScoreType;
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

        protected double getSearchScoreValue(Feature feature)
        {
            return Double.parseDouble(MS2ExtraInfoDef.getSearchScore(feature, searchScoreName));
        }
    }

}
