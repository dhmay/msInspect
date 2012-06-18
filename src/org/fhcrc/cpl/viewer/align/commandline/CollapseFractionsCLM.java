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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.viewer.align.QuantileRegressionAligner;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.MatrixUtil;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;


import java.util.*;
import java.util.List;
import java.io.File;
import java.io.IOException;
import java.awt.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class CollapseFractionsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CollapseFractionsCLM.class);

    protected boolean showCharts = true;

    protected PeptideArrayAnalyzer arrayAnalyzer;

    protected boolean normalizeIntensities = true;

    protected File outFile;

    protected int minFeaturesInFracAndUnfrac = 15;
    protected int minFeaturesForRegress = 15;

    public CollapseFractionsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "collapsefractions";
        mShortDescription = "Collapse multiple fractions into a single featureset";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true,
                                "peptide array in which first run is unfractionated, all the rest fractionated"),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new BooleanArgumentDefinition("normalize", false, "Normalize fraction intensities?", normalizeIntensities),
                        new FileToWriteArgumentDefinition("out", true, "output file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        showCharts = getBooleanArgumentValue("showcharts");
        normalizeIntensities = getBooleanArgumentValue("normalize");
        try
        {
            arrayAnalyzer = new PeptideArrayAnalyzer(getUnnamedFileArgumentValue());
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load peptide array",e);
        }
        outFile = getFileArgumentValue("out");
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<Double[]> intensitiesAllRuns = new ArrayList<Double[]>();

        for (int i=0; i<arrayAnalyzer.getRunCount(); i++)
            intensitiesAllRuns.add(arrayAnalyzer.getRunIntensities(arrayAnalyzer.getRunNames().get(i)));
        Double[] firstRunIntensities = intensitiesAllRuns.get(0);

        if (normalizeIntensities)
            normalizeIntensities(intensitiesAllRuns);
        int nrow = arrayAnalyzer.getRowMaps().length;

        Double[] featureFracIntensitySums = new Double[nrow];
        //Do not add a sum value if there's a missing value between two non-missing values
        int numFoundMissingBracketedByPresent = 0;
        List<Integer> fractionCounts = new ArrayList<Integer>();
        for (int i=0; i<nrow; i++)
        {
            double rowFracSum = 0;
            boolean foundNonMissingFrac = false;
            boolean foundMissingAfterNonMissing = false;
            boolean foundNonMissingAfterMissingAfterNonMissing = false;
            int fractionCount = 0;
            for (int j=1; j<arrayAnalyzer.getRunCount(); j++)
            {
                if (intensitiesAllRuns.get(j)[i] != null)
                {
                    fractionCount++;
                    foundNonMissingFrac = true;
                    if (foundMissingAfterNonMissing)
                        foundNonMissingAfterMissingAfterNonMissing = true;
                    rowFracSum += intensitiesAllRuns.get(j)[i];

//                    else ApplicationContext.infoMessage("Present, missing, present: " + intensitiesAllRuns.get(1)[i] + ", " + intensitiesAllRuns.get(2)[i] + ", " + intensitiesAllRuns.get(3)[i] + ", " + intensitiesAllRuns.get(4)[i]);
                }
                else
                    if (foundNonMissingFrac)
                        foundMissingAfterNonMissing = true;
            }
            if (foundNonMissingAfterMissingAfterNonMissing)
                numFoundMissingBracketedByPresent++;
            if (rowFracSum > 0 && !foundNonMissingAfterMissingAfterNonMissing)
            {
                featureFracIntensitySums[i] = rowFracSum;
                fractionCounts.add(fractionCount);
            }
            else fractionCounts.add(0);
        }
        ApplicationContext.infoMessage("Number discarded because found missing in between present: " + numFoundMissingBracketedByPresent);

        Map<Integer, List<Double>> unfracFractionInCommonLogIntensitiesMap = new HashMap<Integer, List<Double>>();
        Map<Integer, List<Double>> fracFractionInCommonLogIntensitiesMap = new HashMap<Integer, List<Double>>();

        for (int i=0; i<nrow; i++)
        {
            if (featureFracIntensitySums[i] != null && firstRunIntensities[i] != null)
            {
                int numFractionsUsed = fractionCounts.get(i);
                List<Double> unfracLogIntensitiesThisFractionCount =
                        unfracFractionInCommonLogIntensitiesMap.get(numFractionsUsed);
                List<Double> fracLogIntensitiesThisFractionCount =
                        fracFractionInCommonLogIntensitiesMap.get(numFractionsUsed);
                if (unfracLogIntensitiesThisFractionCount == null)
                {
                    unfracLogIntensitiesThisFractionCount = new ArrayList<Double>();
                    unfracFractionInCommonLogIntensitiesMap.put(numFractionsUsed,
                            unfracLogIntensitiesThisFractionCount);
                    fracLogIntensitiesThisFractionCount = new ArrayList<Double>();
                    fracFractionInCommonLogIntensitiesMap.put(numFractionsUsed,
                            fracLogIntensitiesThisFractionCount);
                }
                unfracLogIntensitiesThisFractionCount.add(Math.log(firstRunIntensities[i]));
                fracLogIntensitiesThisFractionCount.add(Math.log(featureFracIntensitySums[i]));
            }
        }

        //todo: technically this should only need to be done for runCount > 1

        for (int i=1; i<=arrayAnalyzer.getRunCount(); i++)
        {
            if (unfracFractionInCommonLogIntensitiesMap.containsKey(i))
            {
                System.err.println("Found " + unfracFractionInCommonLogIntensitiesMap.get(i).size() +
                        " points in " + i + " fractions");
                if (unfracFractionInCommonLogIntensitiesMap.get(i).size() < minFeaturesForRegress)
                {
                    ApplicationContext.infoMessage("\tToo few for regression, so removing instead of adjusting");
                    int numRemoved = 0;

                    for (int j=0; j<featureFracIntensitySums.length; j++)
                    {
                        if (fractionCounts.get(j) == i && featureFracIntensitySums[j] != null)
                        {
                            featureFracIntensitySums[j] = null;
                            numRemoved++;
                        }

                    }
                    ApplicationContext.infoMessage("\tRemoved " + numRemoved + " intensities found in " + i + " fractions");
                    continue;
                }

                //normalize intensities per fraction
                double[] coeffs = calcGoodRegressionCoeffs(fracFractionInCommonLogIntensitiesMap.get(i),
                        unfracFractionInCommonLogIntensitiesMap.get(i)) ;

                if (normalizeIntensities)
                {
                    for (int j=0; j<featureFracIntensitySums.length; j++)
                    {

                        if (fractionCounts.get(j) == i && featureFracIntensitySums[j] != null)
                        {
                            featureFracIntensitySums[j] =
                                    RegressionUtilities.mapValueUsingCoefficients(coeffs, featureFracIntensitySums[j]);
                        }

                    }
                }


                if (showCharts)
                {
                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(fracFractionInCommonLogIntensitiesMap.get(i),
                            unfracFractionInCommonLogIntensitiesMap.get(i), "in"+i+"fracs");
                    pwsp.setAxisLabels("Fractionated","Unfractionated");
                    pwsp.addLineOrCurve(coeffs);
                    pwsp.addLineOrCurve(new double[]{0,1});
                    pwsp.displayInTab();
                }

                if (normalizeIntensities)
                {
                    for (int j=0; j<fracFractionInCommonLogIntensitiesMap.get(i).size(); j++)
                        fracFractionInCommonLogIntensitiesMap.get(i).set(j,
                                RegressionUtilities.mapValueUsingCoefficients(coeffs,
                                        fracFractionInCommonLogIntensitiesMap.get(i).get(j)));
                }
                else
                    ApplicationContext.infoMessage("Note: Skipping intensity normalization for features in " + i + " fractions");


                if (showCharts)
                {
                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(fracFractionInCommonLogIntensitiesMap.get(i),
                            unfracFractionInCommonLogIntensitiesMap.get(i), "FracUnfrac"+i);
                    pwsp.addLineOrCurve(new double[] {0,1});
                    pwsp.displayInTab();
                    ApplicationContext.infoMessage("Log-intensity correlation between fractionated and unfrac for peptides in " + i + " fractions: " +
                            BasicStatistics.correlationCoefficient(unfracFractionInCommonLogIntensitiesMap.get(i),
                                    fracFractionInCommonLogIntensitiesMap.get(i)));
                }
            }
        }

        //Now we write out a feature file with features with appropriate intensities
        Map<Integer, Map<String, List<Feature>>> detailsMap = null;
        try
        {
            ApplicationContext.infoMessage("Loading details file...");
            detailsMap = arrayAnalyzer.loadDetailsRowMapsById();
            ApplicationContext.infoMessage("Loaded " + detailsMap.size() + " array rows' details from details file");

        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to process details file",e); 
        }

        List<Feature> outFeatures = new ArrayList<Feature>();
        for (int i=0; i<nrow; i++)
        {
            int id = i+1;
            if (featureFracIntensitySums[i] == null)
                continue;
            double newIntensity = featureFracIntensitySums[i];   
            String firstAvailableRun = detailsMap.get(id).keySet().iterator().next();
            Feature firstDetailFeature = detailsMap.get(id).get(firstAvailableRun).get(0);
            double intensityScaleFactor = newIntensity / firstDetailFeature.getIntensity();
            firstDetailFeature.setIntensity((float) newIntensity);
            firstDetailFeature.setTotalIntensity((float) intensityScaleFactor * firstDetailFeature.getTotalIntensity());

            outFeatures.add(firstDetailFeature);
        }
        ApplicationContext.infoMessage("Created " + outFeatures.size() + " features to save");

        FeatureSet outFeatureSet = new FeatureSet();
        outFeatureSet.setFeatures(outFeatures.toArray(new Feature[outFeatures.size()]));
        try
        {
            outFeatureSet.save(outFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to save output file",e);
        }

        if (showCharts)
        {
            new PanelWithHistogram(fractionCounts, "FractionCount").displayInTab();

            List<Double> unfracInCommonLogIntensities = new ArrayList<Double>();
            List<Double> fracInCommonLogIntensities = new ArrayList<Double>();

            for (int i=0; i<nrow; i++)
            {
                if (featureFracIntensitySums[i] != null && firstRunIntensities[i] != null)
                {
                    unfracInCommonLogIntensities.add(Math.log(firstRunIntensities[i]));
                    fracInCommonLogIntensities.add(Math.log(featureFracIntensitySums[i]));
                }
            }

            PanelWithScatterPlot fracUnfracScatter =
                    new PanelWithScatterPlot(unfracInCommonLogIntensities, fracInCommonLogIntensities, "Frac vs unfrac");
            fracUnfracScatter.setAxisLabels("Unfractionated","Fractionated");
            fracUnfracScatter.addLineOrCurve(new double[] {0,1});
            fracUnfracScatter.displayInTab();

            for (int i=2; i<=arrayAnalyzer.getRunCount(); i++)
                if (unfracFractionInCommonLogIntensitiesMap.containsKey(i))
                {


                }
//            fracUnfracScatterByFracCount.setAxisLabels("Unfractionated","Fractionated");
//            fracUnfracScatterByFracCount.setPointSize(2);
//            fracUnfracScatterByFracCount.setShowLegend(true);
//            fracUnfracScatterByFracCount.setSeriesColors(new Color[]{Color.BLUE, Color.RED, Color.GREEN, Color.ORANGE});
//            fracUnfracScatterByFracCount.setSeriesColor(0, Color.BLUE);
//            fracUnfracScatterByFracCount.setSeriesColor(1, Color.RED);
//            fracUnfracScatterByFracCount.setSeriesColor(2, Color.GREEN);
//            fracUnfracScatterByFracCount.setSeriesColor(3, Color.ORANGE);
//            fracUnfracScatterByFracCount.addLineOrCurve(new double[] {0,1});
//
//            fracUnfracScatterByFracCount.displayInTab();

            ApplicationContext.infoMessage("Log-intensity correlation between fractionated and unfrac: " +
                    BasicStatistics.correlationCoefficient(unfracInCommonLogIntensities, fracInCommonLogIntensities));
            ApplicationContext.infoMessage("In-common median intensity, frac: " +
                    BasicStatistics.median(fracInCommonLogIntensities) + ", no frac: " +
                    BasicStatistics.median(unfracInCommonLogIntensities));

            List<Double> unfracNotNullLogIntensities = new ArrayList<Double>();
            for (Double intensity : firstRunIntensities)
                if (intensity != null) unfracNotNullLogIntensities.add(Math.log(intensity));
            new PanelWithHistogram(unfracNotNullLogIntensities, "Unfractionated").displayInTab();
            ApplicationContext.infoMessage("Unfractionated intensities: n=" + unfracNotNullLogIntensities.size() +
                    ", median=" + BasicStatistics.median(unfracNotNullLogIntensities));

            List<Double> fracNotNullLogIntensities = new ArrayList<Double>();
            for (Double intensity : featureFracIntensitySums)
                if (intensity != null)
                {
                    fracNotNullLogIntensities.add(Math.log(intensity));
                }
            new PanelWithHistogram(fracNotNullLogIntensities, "Fractionated").displayInTab();
            ApplicationContext.infoMessage("Fractionated intensities: n=" + fracNotNullLogIntensities.size() +
                    ", median=" + BasicStatistics.median(fracNotNullLogIntensities));
        }
    }

    protected void normalizeIntensities(List<Double[]> intensitiesAllRuns)
             throws CommandLineModuleExecutionException
    {
        Double[] firstRunIntensities = intensitiesAllRuns.get(0);

        PanelWithScatterPlot pwspAllRunsNorm = null;
        if (showCharts)
        {
            pwspAllRunsNorm = new PanelWithScatterPlot();
            pwspAllRunsNorm.setAxisLabels("Original","Mapped");
            pwspAllRunsNorm.setName("RunNormalization");
        }


        ApplicationContext.infoMessage("Rows in array: " + firstRunIntensities.length);
        for (int i=1; i<arrayAnalyzer.getRunCount(); i++)      
        {
            ApplicationContext.infoMessage("Normalizing fraction " + i);
            Double[] intensitiesThisRun = intensitiesAllRuns.get(i);
            List<Double> firstRunIntensitiesJustThisRun = new ArrayList<Double>();
            List<Double> thisRunIntensitiesJustThisRun = new ArrayList<Double>();

            for (int j=0; j<firstRunIntensities.length; j++)
            {
                if (intensitiesThisRun[j] == null || firstRunIntensities[j] == null)
                    continue;
                boolean foundInAnotherRun = false;
                for (int k=1; k<arrayAnalyzer.getRunCount(); k++)
                {
                    if (k!=i && intensitiesAllRuns.get(k)[j] != null)
                    {
                        foundInAnotherRun = true;
                        break;
                    }
                }
                if (!foundInAnotherRun)
                {
                    firstRunIntensitiesJustThisRun.add(firstRunIntensities[j]);
                    thisRunIntensitiesJustThisRun.add(intensitiesThisRun[j]);
                }

            }
            ApplicationContext.infoMessage("values in common between ONLY unfrac and fraction " + i + ":" +
                    firstRunIntensitiesJustThisRun.size());
            List<Double> logFirst = new ArrayList<Double>();
            List<Double> logThis = new ArrayList<Double>();
            double minx = Double.MAX_VALUE;
            double maxx = Double.MIN_VALUE;
            for (int l=0; l<firstRunIntensitiesJustThisRun.size(); l++)
            {
                logFirst.add(Math.log(firstRunIntensitiesJustThisRun.get(l)));
                logThis.add(Math.log(thisRunIntensitiesJustThisRun.get(l)));
                minx = Math.min(Math.log(thisRunIntensitiesJustThisRun.get(l)), minx);
                maxx = Math.max(Math.log(thisRunIntensitiesJustThisRun.get(l)), maxx);
            }

            if (logFirst.size() < minFeaturesInFracAndUnfrac)
            {
                ApplicationContext.infoMessage("WARNING: unable to normalize run " + i + ", not enough in common with unfractionated run");
                continue;
            }


            double[] coeffs = null;
//                try
//                {
//                    //todo: un-hardcode
////                    coeffs = RegressionUtilities.modalRegression(logFirst, logThis, 1, 6, 2.0);
//                }
//                catch (IOException e)
//                {
//                    throw new CommandLineModuleExecutionException("Failure in modal regression",e);
//                }
            //todo: un-hardcode
            Pair<double[], double[]> lowLeverageLowResPoints =
                    RegressionUtilities.selectValuesWithLowLeverageAndStudentizedResidual(logThis, logFirst,
                            17, 1.8, false, 1, false, true);
//try
//{
//            coeffs = RegressionUtilities.modalRegression(
//                    lowLeverageLowResPoints.first, lowLeverageLowResPoints.second);
//}
//catch (IOException e)
//{throw new CommandLineModuleExecutionException(e);}
            coeffs = calcGoodRegressionCoeffs(logThis, logFirst);


            if (showCharts)
            {
                PanelWithScatterPlot pwsp = new PanelWithScatterPlot(lowLeverageLowResPoints.first,
                        lowLeverageLowResPoints.second, "logintensities" + i);
                pwsp.addData(logThis, logFirst, "used datapoints");
                pwsp.setAxisLabels("current fraction","unfractionated run");
                pwsp.addLineOrCurve(coeffs, minx, maxx);
                pwsp.displayInTab();
            }

            List<Double> origNotNullIntensities = null;
            if (showCharts)
            {
                origNotNullIntensities = new ArrayList<Double>();
                for (int l=0; l<intensitiesThisRun.length; l++)
                {
                    if (intensitiesThisRun[l] != null)
                    {
                        origNotNullIntensities.add(Math.log(intensitiesThisRun[l]));
                    }
                }
            }

            for (int l=0; l<intensitiesThisRun.length; l++)
            {
                if (intensitiesThisRun[l] != null)
                    intensitiesThisRun[l] = Math.exp(RegressionUtilities.mapValueUsingCoefficients(coeffs,
                            Math.log(intensitiesThisRun[l])));
            }

            if (showCharts)
            {
                List<Double> mappedNotNullIntensities = new ArrayList<Double>();

                for (int l=0; l<intensitiesThisRun.length; l++)
                {
                    if (intensitiesThisRun[l] != null)
                    {
                        mappedNotNullIntensities.add(Math.log(intensitiesThisRun[l]));
                    }
                }

                pwspAllRunsNorm.addData(origNotNullIntensities,
                        mappedNotNullIntensities, "mappedintensities" + i);
            }

        }
        if (showCharts)
            pwspAllRunsNorm.displayInTab();


    }

    protected double[] calcGoodRegressionCoeffs(List<Double> xvalsList, List<Double> yvalsList)
    {
        double[] xvals = new double[xvalsList.size()];
        double[] yvals = new double[xvalsList.size()];
        for (int i=0; i<xvals.length; i++)
        {
            xvals[i] = xvalsList.get(i);
            yvals[i] = yvalsList.get(i);
        }
        return calcGoodRegressionCoeffs(xvals, yvals);
    }

    protected double[] calcGoodRegressionCoeffs(double[] xvals, double[] yvals)
    {
            Pair<double[], double[]> lowLeverageLowResPoints =
                    RegressionUtilities.selectValuesWithLowLeverageAndStudentizedResidual(xvals, yvals,
                            17, 1.8, false, 1, false, true);
//try
//{
//            coeffs = RegressionUtilities.modalRegression(
//                    lowLeverageLowResPoints.first, lowLeverageLowResPoints.second);
//}
//catch (IOException e)
//{throw new CommandLineModuleExecutionException(e);}
            return RegressionUtilities.symmetricalRegression(
                    lowLeverageLowResPoints.first, lowLeverageLowResPoints.second);
    }


}
