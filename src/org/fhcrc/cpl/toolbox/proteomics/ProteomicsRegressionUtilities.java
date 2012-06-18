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
package org.fhcrc.cpl.toolbox.proteomics;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.MatrixUtil;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;

import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * Utilities for regression analysis specific to msInspect -- mostly making use of Feature.
 * These utilities were originally in the class that became toolbox.RegressionUtilities
 */
public class ProteomicsRegressionUtilities
{
    protected static Logger _log = Logger.getLogger(ProteomicsRegressionUtilities.class);

    //modes for matching features by scan or by time
    public static final int REGRESSION_MODE_SCAN = 0;
    public static final int REGRESSION_MODE_TIME = 1;
    //time-based matching has been shown to be far superior
    public static final int DEFAULT_REGRESSION_MODE = REGRESSION_MODE_TIME;

    public static double[] modalRegressionIterateDegree(double[] xset, double[] yset,
                                                        int maxDegree, double leverageNumerator,
                                                        double maxStudentizedResidual,
                                                        boolean showCharts)
            throws IOException
    {
        Feature[] dummyFeatures = new Feature[xset.length];
        double[] xsetOrig = xset;
        double[] ysetOrig = yset;
        for (int j=0; j<xset.length; j++)
        {
            dummyFeatures[j] = new Feature();
            dummyFeatures[j].setTime((float) xset[j]);
            AmtExtraInfoDef.setObservedHydrophobicity(dummyFeatures[j], yset[j]);
        }
        for (int degree=1; degree<maxDegree; degree++)
        {
            _log.debug("Applying degree " + degree);
            boolean modal = true;
            if (degree==1) modal = false;
            dummyFeatures = selectFeaturesWithLowLeverageAndStudentizedResidual(
                    dummyFeatures,
                    xset, yset,
                    leverageNumerator,
                    maxStudentizedResidual, modal, degree, showCharts, false);
            xset = new double[dummyFeatures.length];
            yset = new double[dummyFeatures.length];
            for (int j=0; j<dummyFeatures.length; j++)
            {
                xset[j] = dummyFeatures[j].getTime();
                yset[j] = AmtExtraInfoDef.getObservedHydrophobicity(dummyFeatures[j]);
            }
            _log.debug("Remaining pairs: " + xset.length);
        }

        double[] regressionResult = RegressionUtilities.modalRegression(xset, yset, maxDegree);

        if (showCharts)
        {
            int maxX = 0;
            int minX = Integer.MAX_VALUE;

            for (double x : xset)
            {
                if (x > maxX)
                    maxX = (int) x;
                if (x < minX)
                    minX = (int) x;
            }
            ScatterPlotDialog spd = new ScatterPlotDialog();

            int numDotsOnChart = (maxX-minX+1) / 2;
            double[] chartXVals = new double[numDotsOnChart];
            double[] chartYVals = new double[numDotsOnChart];

            for (int j=0; j<numDotsOnChart; j++)
            {
                chartXVals[j] = minX + (2 * j);
                chartYVals[j] =
                        RegressionUtilities.mapValueUsingCoefficients(regressionResult, chartXVals[j]);
            }
            spd.addData(xset, yset,
                        "Matches with low stud. res.");
            spd.addData(xsetOrig, ysetOrig, "all mass matches");
            spd.addData(chartXVals, chartYVals, "regression function");
            spd.setAxisLabels("X","Y");
            spd.setVisible(true);
        }

        return regressionResult;

    }

    /**
     * return an array containing the subset of the allFeatures array in which all
     * features have leverage < leverageNumerator/n
     *
     * Used in AmtDatabaseBuilder, and possibly also in selecting MS1 features with
     * no embedded MS2 for matching
     * @param allFeatures
     * @param leverageNumerator
     * @param scanOrTimeMode
     * @return
     */
    public static Feature[] selectFeaturesWithLowLeverage(Feature[] allFeatures,
                                                          double leverageNumerator,
                                                          int scanOrTimeMode)
    {
        int n = allFeatures.length;
        double[] timeValues = new double[n];
        for (int i=0; i<n; i++)
            timeValues[i] = scanOrTimeMode == REGRESSION_MODE_TIME ?
                allFeatures[i].getTime() : allFeatures[i].getScan();
        double[] leverages = BasicStatistics.leverages(timeValues);

        double maxLeverage = leverageNumerator / (double) n;
//System.err.println("max leverage: " + maxLeverage + ", n: " + n);
int passed=0;
        return selectFeaturesWithLowAbsoluteSomething(allFeatures, leverages, maxLeverage);
    }


    public static Feature[] selectFeaturesWithLowAbsoluteSomething(Feature[] allFeatures,
                                                          double[] somethings,
                                                          double maxSomething)
    {
        int[] indexes = RegressionUtilities.selectIndexesWithLowAbsoluteSomething(somethings, maxSomething);
        Feature[] result = new Feature[indexes.length];
        for (int i=0; i<indexes.length; i++)
            result[i] = allFeatures[indexes[i]];

        return result;
    }


    /**
     *
     * @param features
     * @param xValues
     * @param yValues
     * @param maxStudentizedResidual
     * @param modalRegression
     * @param degree
     * @param showCharts
     * @param assumePositiveCorr Assume positive correlation between x and y?  If so, we can try to adjust for bias
     * @return
     */
    public static Feature[] selectFeaturesWithLowStudentizedResidual(
            Feature[] features, double[] xValues, double[] yValues,
            double maxStudentizedResidual, boolean modalRegression, int degree, boolean showCharts,
            boolean assumePositiveCorr)
    {
        Feature[] result =  selectFeaturesWithLowLeverageAndStudentizedResidual(
            features, xValues, yValues, Integer.MAX_VALUE, maxStudentizedResidual,
                modalRegression, degree, showCharts, assumePositiveCorr);
        return result;
    }

    /**
     * Select low-leverage features, where leverage is defined by mass
     * @param features
     * @param maxLeverage
     * @return
     */
    public static List<Feature> selectFeaturesWithLowAbsMassLeverage(Feature[] features, float maxLeverage)
    {
        double[] masses = new double[features.length];
        for (int i=0; i<features.length; i++)
            masses[i] = features[i].getMass();
        double[] leverages = BasicStatistics.leverages(masses);
        List<Feature> result = new ArrayList<Feature>();
        for (int i=0; i<features.length; i++)
            if (Math.abs(leverages[i]) < maxLeverage)
                result.add(features[i]);
        return result;
    }

    /**
     * Calculate studentized residuals
      * @param xValues
     * @param yValues
     * @return
     */
    public static double[] calculateStudentizedResiduals(double[] xValues, double[] yValues)
    {
        double[] regressionResult = MatrixUtil.linearRegression(xValues, yValues);
        double[] leverages = BasicStatistics.leverages(xValues);


        double[] residuals = new double[xValues.length];
        for (int i=0; i<xValues.length; i++)
        {
            double predictedValue =
                    RegressionUtilities.mapValueUsingCoefficients(regressionResult,
                                                           xValues[i]);
            residuals[i] = yValues[i] - predictedValue;
        }
        return BasicStatistics.studentizedResiduals(xValues,
                        residuals, leverages);
    }


    /**
     *
     * @param features
     * @param xValues
     * @param yValues
     * @param leverageNumerator
     * @param maxStudentizedResidual
     * @param modalRegression
     * @param degree
     * @param showCharts
     * @param assumePositiveCorr Assume positive correlation between x and y?  If so, we can try to adjust for bias
     * @return
     */
    public static Feature[] selectFeaturesWithLowLeverageAndStudentizedResidual(
            Feature[] features, double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts, boolean assumePositiveCorr)
    {
        int[] indexes = RegressionUtilities.selectIndexesWithLowLeverageAndStudentizedResidual(
                xValues, yValues, leverageNumerator, maxStudentizedResidual,
                modalRegression, degree, showCharts, assumePositiveCorr);
        Feature[] result = new Feature[indexes.length];
        for (int i=0; i<result.length; i++)
            result[i] = features[indexes[i]];
        return result;
    }


}
