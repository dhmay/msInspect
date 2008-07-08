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
package org.fhcrc.cpl.viewer.util;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.viewer.gui.util.ScatterPlotDialog;
import org.labkey.common.tools.MatrixUtil;
import org.labkey.common.tools.BasicStatistics;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.io.IOException;
import java.io.InputStream;


/**
 * Utilities for regression analysis
 */
public class RegressionUtilities
{
    protected static Logger _log = Logger.getLogger(RegressionUtilities.class);

    //key definitions for accessing regression results
    public static final String REGRESSION_SLOPE_KEY = "REGRESSION_SLOPE_KEY";
    public static final String REGRESSION_INTERCEPT_KEY = "REGRESSION_INTERCEPT_KEY";
    public static final String REGRESSION_SIGMA_KEY = "REGRESSION_SIGMA_KEY";
    //modes for matching features by scan or by time
    public static final int REGRESSION_MODE_SCAN = 0;
    public static final int REGRESSION_MODE_TIME = 1;
    //time-based matching has been shown to be far superior
    public static final int DEFAULT_REGRESSION_MODE = REGRESSION_MODE_TIME;

    public static final int DEFAULT_MAX_MILLIS_FOR_ROBUST_REGRESSION = 180000;
    


    public static double[] robustRegression(double[] xset, double[] yset)
    {
        return robustRegression(xset, yset, DEFAULT_MAX_MILLIS_FOR_ROBUST_REGRESSION);
    }

    public static double[] robustRegression(double[] xset, double[] yset, int millis)
    {
        Map<String,double[]> variableValueMap = new HashMap<String,double[]>(2);
        variableValueMap.put("x",xset);
        variableValueMap.put("y",yset);
        //robust regression can be a time-consuming business, need to allow
        //plenty of time for R to return
        //TODO: parameterize time to wait?
        return RInterface.processRCoefficientResponse(RInterface.evaluateRExpression("coefficients(rlm((y~x)));",
                variableValueMap, new String[] {"MASS"}, millis));
    }

    /**
     * Given a slope and intercept of a line relating x and y, predict x from y
     * @param slope
     * @param intercept
     * @param x
     * @return
     */
    public static double predictYFromX(double slope, double intercept, double x)
    {
        return slope * x + intercept;
    }

    /**
     * Given a slope and intercept of a line relating x and y, predict y from x
     * @param slope
     * @param intercept
     * @return
     */
    public static double predictXFromY(double slope, double intercept, double y)
    {
        return (y - intercept) / slope;
    }
    


    /**
     * Default: two coefficients
     * @param xset
     * @param yset
     * @return
     * @throws java.io.IOException
     */
    public static double[] modalRegression(double[] xset, double[] yset)
            throws IOException
    {
        return modalRegression(xset, yset, 1);
    }

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
            dummyFeatures = selectFeaturesWithLowStudentizedResidual(
                    dummyFeatures,
                    xset, yset,
                    leverageNumerator,
                    maxStudentizedResidual, modal, degree, showCharts);
            xset = new double[dummyFeatures.length];
            yset = new double[dummyFeatures.length];
            for (int j=0; j<dummyFeatures.length; j++)
            {
                xset[j] = dummyFeatures[j].getTime();
                yset[j] = AmtExtraInfoDef.getObservedHydrophobicity(dummyFeatures[j]);
            }
            _log.debug("Remaining pairs: " + xset.length);
        }

        double[] regressionResult = modalRegression(xset, yset, maxDegree);

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
                        mapValueUsingCoefficients(regressionResult, chartXVals[j]);
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
     * Call Yan's Modal Regression R code.
     * This code is dependent on the "quantreg" library being installed.  If that package
     * isn't installed, it will throw IOException.
     *
     * Steps to install:
     *
     * source("http://bioconductor.org/biocLite.R")
     * biocLite(c("quantreg"))
     *
     * @param xset
     * @param yset
     * @param degree The degree of the polynomial.  Minimum 1
     * @return
     * @throws IOException
     */
    public static double[] modalRegression(double[] xset, double[] yset, int degree)
            throws IOException
    {
        if (degree < 1)
            throw new RuntimeException("Expected degree parameter >=1");

        //first write the source, then run the function
        StringBuffer modalRegressionSourceCodeBuf = new StringBuffer();
        InputStream in = RegressionUtilities.class.getResourceAsStream("modal_regression.R");
        int inByte;
        while ((inByte = in.read()) != -1)
            modalRegressionSourceCodeBuf.append((char) inByte);
        in.close();

        Map<String,double[]> variableValueMap = new HashMap<String,double[]>(2);
        variableValueMap.put("x",xset);
        variableValueMap.put("y",yset);

        //for higher-degree regression
        String expandXsetCommand = "";
        if (degree > 1)
        {
            String expansionString = "";
            for (int i=1; i<=(degree); i++)
            {
                expansionString = expansionString + "x^" + i;
                if (i < degree)
                    expansionString = expansionString + ", ";
            }
            expandXsetCommand =
                    "x<-matrix(c(" + expansionString + "), nrow=length(x), byrow=F);";
        }

        String rCommand = modalRegressionSourceCodeBuf.toString() +
                expandXsetCommand + " modal_regress(x,y)";

        //modal regression can be a time-consuming business, need to allow
        //plenty of time for R to return
        //TODO: parameterize time to wait?

        String rResponse = null;

        try
        {
            rResponse = RInterface.evaluateRExpression(rCommand, variableValueMap,
                    new String[] {"quantreg"},
                    150000);
        }
        catch (RuntimeException e)
        {
            throw new IOException("Failure running R for modal regression.  Type: " +
                    e.getClass().getName() + ", Message: " + e.getMessage() +
                    ".  Make sure quantreg package is installed");
        }
        return RInterface.processRCoefficientResponse(rResponse);
    }

    /**
     * Damon's poor man's version of modal regression.
     * This is a poor man's attempt to remove the effect of noise bias from the regression.
     * First I do a regular regression.  This will have a low-slope, high-intercept
     * bias because of the random false mass matches all over the place.
     * Next, I do a regression to predict Time in terms of H.  This should have exactly
     * the OPPOSITE bias, if the distribution of crap is uniform.
     * Then I find the intersection point of the two lines and return a line through
     * that point with the average slope of the two regression lines.
     *
     * Works pretty well in many cases, but it has a major assumption, namely that the
     * distribution of the crap is uniform around the good stuff.
     * @param xvals
     * @param yvals
     * @return
     */
    public static double[] poorMansModalRegression(double[] xvals, double[] yvals)
    {
        //calculate the regression line predicting H from time
        double[] regressionResult =
                MatrixUtil.linearRegression(xvals, yvals);

        //calculate the regression line predicting time from H
        double[] reverseRegressionResult =
                MatrixUtil.linearRegression(yvals, xvals);

        //solve to determine intersection point
        //forward line: y = ax + b
        //reverse line: x = cy + d
        double a = regressionResult[1];
        double b = regressionResult[0];

        //reverse regression gives x in terms of y
        double cBackward = reverseRegressionResult[1];
        double dBackward = reverseRegressionResult[0];

        //Turn it around, y in terms of x
        double c = 1.0 / cBackward;
        double d = - dBackward / cBackward;

        //find intersection point
        double xIntersect = (d - b) / (a - c);
        double yIntersect = a * xIntersect + b;

        //Correct slope is average of the two slopes
        double averageSlope = (a + c) / 2;
        //find intercept
        double averageIntercept = yIntersect -(averageSlope * xIntersect);

        double[] result = new double[2];
        result[0] = averageIntercept;
        result[1] = averageSlope;

        _log.debug("poorMansModalRegression: Forward slope and intercept: " + a + ", " + b);
        _log.debug("poorMansModalRegression: Reverse slope and intercept: " + c + ", " + d);
        _log.debug("poorMansModalRegression: Intersection: " + xIntersect + ", " + yIntersect);

        return result;
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
        Integer[] indexes = selectIndexesWithLowAbsoluteSomething(somethings, maxSomething);
        Feature[] result = new Feature[indexes.length];
        for (int i=0; i<indexes.length; i++)
            result[i] = allFeatures[indexes[i]];

        return result;
    }

    public static Integer[] selectIndexesWithLowAbsoluteSomething(double[] somethings,
                                                      double maxSomething)
    {
        List<Integer> resultList = new ArrayList<Integer>();

        for (int i=0; i<somethings.length; i++)
        {
            if (Math.abs(somethings[i]) < maxSomething)
            {
                resultList.add(i);
            }
        }

        return resultList.toArray(new Integer[0]);
    }

    public static Feature[] selectFeaturesWithLowLeverageAndStudentizedResidual(
            Feature[] features, double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts
            )
    {
        Feature[] result = selectFeaturesWithLowLeverageAndOrStudentizedResidual(
            features, xValues, yValues, leverageNumerator, maxStudentizedResidual,
                modalRegression, degree, false, showCharts);
        return result;
    }

    public static Feature[] selectFeaturesWithLowStudentizedResidual(
            Feature[] features, double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts)
    {
        Feature[] result =  selectFeaturesWithLowLeverageAndOrStudentizedResidual(
            features, xValues, yValues, leverageNumerator, maxStudentizedResidual,
                modalRegression, degree, true, showCharts);
        return result;
    }

    public static Feature[] selectFeaturesWithLowLeverageAndOrStudentizedResidual(
            Feature[] features, double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean keepHighLeverage, boolean showCharts)
    {
        int nAll = features.length;
        double maxLeverage = leverageNumerator / (double) nAll;
        _log.debug("selectFeaturesWithLowLeverageAndOrStudentizedResidual, starting with " + nAll +
                   ", maxLeverageNumerator=" + leverageNumerator + ", maxLeverage=" + maxLeverage +
                   ", maxStudRes=" + maxStudentizedResidual);

        double[] leverages = BasicStatistics.leverages(xValues);
//for (double leverage : leverages)
//    System.err.println("" + leverage);
//PanelWithHistogram pwh = new PanelWithHistogram(leverages);
//pwh.displayDialog("leverages");
        Integer[] indexesWithLowLeverage = selectIndexesWithLowAbsoluteSomething(
                leverages, maxLeverage);

        int nLowLeverage = indexesWithLowLeverage.length;
        double[] xValuesWithLowLeverage = new double[nLowLeverage];
        double[] yValuesWithLowLeverage = new double[nLowLeverage];
        Feature[] featuresWithLowLeverage = new Feature[nLowLeverage];
        for (int i=0; i<nLowLeverage; i++)
        {
            xValuesWithLowLeverage[i] = xValues[indexesWithLowLeverage[i]];
            yValuesWithLowLeverage[i] = yValues[indexesWithLowLeverage[i]];
            featuresWithLowLeverage[i] = features[indexesWithLowLeverage[i]];
        }
//        double[] regressionResult = MatrixUtil.linearRegression(xValuesWithLowLeverage,
//                yValuesWithLowLeverage);
//        double[] regressionResult = AmtUtilities.robustRegression(xValuesWithLowLeverage,
//                yValuesWithLowLeverage);

        _log.debug("selectFeaturesWithLowLeverageAndOrStudentizedResidual, low leverage: " + nLowLeverage);

        double[] regressionResult = null;
        if (!modalRegression)
            regressionResult = poorMansModalRegression(xValuesWithLowLeverage,
                                                       yValuesWithLowLeverage);
        else
        {
            try
            {
                regressionResult = modalRegression(xValuesWithLowLeverage, yValuesWithLowLeverage, degree);
            }
            catch (IOException e)
            {
                throw new RuntimeException(e);
            }
        }


        double[] residuals = null;
        double[] xValuesForStudRes = null;
        double[] yValuesForStudRes = null;
        Feature[] featuresForStudResCutoff = null;
        if (keepHighLeverage)
        {
            residuals = new double[xValues.length];
            xValuesForStudRes = xValues;
            yValuesForStudRes = yValues;
            featuresForStudResCutoff = features;
        }
        else
        {
            residuals = new double[nLowLeverage];
            xValuesForStudRes = xValuesWithLowLeverage;
            yValuesForStudRes = yValuesWithLowLeverage;
            featuresForStudResCutoff = featuresWithLowLeverage;
        }


        for (int i=0; i<xValuesForStudRes.length; i++)
        {
            double predictedValue =
                    mapValueUsingCoefficients(regressionResult,
                                                           xValuesForStudRes[i]);
            residuals[i] = yValuesForStudRes[i] - predictedValue;
        }

        double[] studentizedResiduals =
                BasicStatistics.studentizedResiduals(xValuesForStudRes,
                        residuals, leverages);
        Feature[] featuresWithLowStudentizedResidual =
                selectFeaturesWithLowAbsoluteSomething(featuresForStudResCutoff,
                        studentizedResiduals, maxStudentizedResidual);
        _log.debug("selectFeaturesWithLowLeverageAndOrStudentizedResidual, result: " + featuresWithLowStudentizedResidual.length);

        if (showCharts)
        {
            int maxX = 0;
            int minX = Integer.MAX_VALUE;
            double[] xValuesWithLowStudRes = new double[featuresWithLowStudentizedResidual.length];
            double[] yValuesWithLowStudRes = new double[featuresWithLowStudentizedResidual.length];

            for (int i=0; i<featuresWithLowStudentizedResidual.length; i++)
            {
                xValuesWithLowStudRes[i] = featuresWithLowStudentizedResidual[i].getTime();
                yValuesWithLowStudRes[i] = AmtExtraInfoDef.getObservedHydrophobicity(featuresWithLowStudentizedResidual[i]);
            }

            for (double x : xValuesWithLowStudRes)
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
                        mapValueUsingCoefficients(regressionResult, chartXVals[j]);
            }
            spd.addData(xValuesWithLowStudRes, yValuesWithLowStudRes,
                        "Matches with low stud. res.");
            spd.addData(xValues, yValues, "all mass matches");
            spd.addData(chartXVals, chartYVals, "regression function");
            spd.setAxisLabels("X","Y");
            spd.setVisible(true);
        }
        return featuresWithLowStudentizedResidual;
    }

    public static double mapValueUsingCoefficients(double[] coefficients, double valueToMap)
    {
        double result = 0;
        for (int j=0; j<coefficients.length; j++)
            result += coefficients[j] * Math.pow(valueToMap, j);
        return result;
    }

    
}
