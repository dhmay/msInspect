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
package org.fhcrc.cpl.toolbox.statistics;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.datastructure.Pair;

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

    public static final int DEFAULT_MAX_MILLIS_FOR_ROBUST_REGRESSION = 180000;


    /**
     * Cover method.  Inefficiently converts to doubles
     * @param categories these should be distinct values -- practically speaking, ints would have been better.
     * double[] for simplicity
     * @param outcomes
     * @return
     */
    public static AnovaResult oneWayAnova(List<? extends Number> categories, List<? extends Number> outcomes)
    {
        double[] categoriesDouble = new double[categories.size()];
        double[] outcomesDouble = new double[categories.size()];

        for (int i=0; i<categories.size(); i++)
        {
            categoriesDouble[i] =  categories.get(i).doubleValue();
            outcomesDouble[i] = outcomes.get(i).doubleValue();
        }

        return oneWayAnova(categoriesDouble, outcomesDouble);
    }

    /**
     *
     * @param categories these should be distinct values -- practically speaking, ints would have been better.
     * double[] for simplicity
     * @param outcomes
     * @return
     */
    public static AnovaResult oneWayAnova(double[] categories, double[] outcomes)
    {
        Map<String,double[]> variableValueMap = new HashMap<String,double[]>(2);
        variableValueMap.put("x",categories);
        variableValueMap.put("y",outcomes);
        String rResponse = RInterface.evaluateRExpression("summary(aov(y~x))",
                variableValueMap, new String[] {"MASS"});
        AnovaResult result = new AnovaResult();

        String textLine2 = rResponse.split("\n")[1];
        String[] stringChunks = textLine2.split("\\s+");
        result.setFValue(Float.parseFloat(stringChunks[4]));
        //If there's an extremely small p-value, "<" will appear in position 5
        if ("<".equals(stringChunks[5]))
            result.setPValue(0f);
        else
            result.setPValue(Float.parseFloat(stringChunks[5]));

        return result;
    }


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

    public static double[] modalRegression(List<? extends Number> xset, List<? extends Number> yset, int degree)
            throws IOException
    {
        double[] xsetArray = new double[xset.size()];
        double[] ysetArray = new double[xset.size()];
        for (int i=0; i<xsetArray.length; i++)
        {
            xsetArray[i] = xset.get(i).doubleValue();
            ysetArray[i] = yset.get(i).doubleValue();

        }
        return modalRegression(xsetArray, ysetArray, degree);
    }

    public static double[] modalRegression(List<? extends Number> xset, List<? extends Number> yset, int degree,
                                           double maxLeverageNumerator, double maxStudentizedResidual)
            throws IOException
    {
        double[] xsetArray = new double[xset.size()];
        double[] ysetArray = new double[xset.size()];
        for (int i=0; i<xsetArray.length; i++)
        {
            xsetArray[i] = xset.get(i).doubleValue();
            ysetArray[i] = yset.get(i).doubleValue();

        }
        Pair<double[], double[]> passingXYValues =
                selectValuesWithLowLeverageAndStudentizedResidual(xsetArray, ysetArray, maxLeverageNumerator,
                maxStudentizedResidual, false, 1, false, true);
        return modalRegression(passingXYValues.first, passingXYValues.second, degree);
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
     * todo: require some minum cardinality (5?), and throw an IllegalArg if not enough datapoints. Because R will fail
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


    public static int[] selectIndexesWithLowAbsoluteSomething(double[] somethings,
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

        int[] result = new int[resultList.size()];
        for (int i=0; i<resultList.size(); i++)
            result[i] = resultList.get(i);
        return result;
    }


    public static double mapValueUsingCoefficients(double[] coefficients, double valueToMap)
    {
        double result = 0;
        for (int j=0; j<coefficients.length; j++)
            result += coefficients[j] * Math.pow(valueToMap, j);
        return result;
    }

    public static List<Double> mapValuesUsingCoefficients(double[] coefficients, List<? extends Number> valuesToMap)
    {
        double[] valuesToMapArray = new double[valuesToMap.size()];

        for (int i=0; i<valuesToMap.size(); i++)
        {
            valuesToMapArray[i] = valuesToMap.get(i).doubleValue();
        }
        double[] resultArray = mapValuesUsingCoefficients(coefficients, valuesToMapArray);
        List<Double> result = new ArrayList<Double>();
        for (double x : resultArray)
            result.add(x);
        return result;
    }

    public static double[] mapValuesUsingCoefficients(double[] coefficients, double[] valuesToMap)
    {
        double[] result = new double[valuesToMap.length];
        for (int i=0; i<result.length; i++)
            result[i] = mapValueUsingCoefficients(coefficients, i);
        return result;
    }


    public static int[] selectIndexesWithLowLeverageAndStudentizedResidual(
            List<? extends Number> xValues, List<? extends Number> yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts, boolean assumePositiveCorr)
    {
        double[] xArray = new double[xValues.size()];
        double[] yArray = new double[xValues.size()];

        for (int i=0; i<xValues.size(); i++)
        {
            xArray[i] = xValues.get(i).doubleValue();
            yArray[i] = yValues.get(i).doubleValue();
        }
        return selectIndexesWithLowLeverageAndStudentizedResidual(xArray, yArray, leverageNumerator,
                maxStudentizedResidual, modalRegression, degree, showCharts, assumePositiveCorr);
    }

    public static Pair<double[], double[]> selectValuesWithLowLeverageAndStudentizedResidual(
            List<? extends Number> xValues, List<? extends Number> yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts, boolean assumePositiveCorr)
    {
        double[] xArray = new double[xValues.size()];
        double[] yArray = new double[xValues.size()];

        for (int i=0; i<xValues.size(); i++)
        {
            xArray[i] = xValues.get(i).doubleValue();
            yArray[i] = yValues.get(i).doubleValue();
        }
        return selectValuesWithLowLeverageAndStudentizedResidual(xArray, yArray,
             leverageNumerator,  maxStudentizedResidual,
             modalRegression,  degree,  showCharts,  assumePositiveCorr);
    }

    public static Pair<double[], double[]> selectValuesWithLowLeverageAndStudentizedResidual(
            double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts, boolean assumePositiveCorr)
    {
        int[] passingIndexes = selectIndexesWithLowLeverageAndStudentizedResidual(xValues, yValues,
             leverageNumerator,  maxStudentizedResidual,
             modalRegression,  degree,  showCharts,  assumePositiveCorr);
        double[] newX = new double[passingIndexes.length];
        double[] newY = new double[passingIndexes.length];
        for (int i=0; i<passingIndexes.length; i++)
        {
            newX[i] = xValues[passingIndexes[i]];
            newY[i] = yValues[passingIndexes[i]];
        }
        return new Pair<double[], double[]>(newX, newY);
    }

    /**
     * Does what it says on the tin
     * @param xValues
     * @param yValues
     * @param leverageNumerator
     * @param maxStudentizedResidual
     * @param modalRegression
     * @param degree
     * @param showCharts
     * @param assumePositiveCorr Assume positive correlation?  If so, we can try to eliminate regression bias
     * @return
     */
    public static int[] selectIndexesWithLowLeverageAndStudentizedResidual(
            double[] xValues, double[] yValues,
            double leverageNumerator, double maxStudentizedResidual,
            boolean modalRegression, int degree, boolean showCharts, boolean assumePositiveCorr)
    {
        int nAll = xValues.length;
        double maxLeverage = leverageNumerator / (double) nAll;
        _log.debug("selectFeaturesWithLowLeverageAndStudentizedResidual, starting with " + nAll +
                   ", maxLeverageNumerator=" + leverageNumerator + ", maxLeverage=" + maxLeverage +
                   ", maxStudRes=" + maxStudentizedResidual);

        double[] leverages = BasicStatistics.leverages(xValues);
//for (double leverage : leverages)
//    System.err.println("" + leverage);
//PanelWithHistogram pwh = new PanelWithHistogram(leverages);
//pwh.displayDialog("leverages");
        int[] indexesWithLowLeverage = selectIndexesWithLowAbsoluteSomething(
                leverages, maxLeverage);

        int nLowLeverage = indexesWithLowLeverage.length;
        double[] xValuesWithLowLeverage = new double[nLowLeverage];
        double[] yValuesWithLowLeverage = new double[nLowLeverage];
        for (int i=0; i<nLowLeverage; i++)
        {
            xValuesWithLowLeverage[i] = xValues[indexesWithLowLeverage[i]];
            yValuesWithLowLeverage[i] = yValues[indexesWithLowLeverage[i]];
        }
        double[] simpleRegressionResult = MatrixUtil.linearRegression(xValuesWithLowLeverage,
                                                           yValuesWithLowLeverage);
        _log.debug("selectFeaturesWithLowLeverageAndStudentizedResidual, low leverage: " + nLowLeverage);
        if (showCharts && _log.isDebugEnabled())
        {
            PanelWithScatterPlot pwsp = new PanelWithScatterPlot(xValuesWithLowLeverage, yValuesWithLowLeverage, "lowleverage");
            pwsp.addLine(simpleRegressionResult[1], simpleRegressionResult[0],
                    BasicStatistics.min(xValuesWithLowLeverage), BasicStatistics.max(xValuesWithLowLeverage));
            pwsp.displayInTab();
        }
        double[] regressionResult = null;
        if (!modalRegression)
        {
            //if not doing modal regression, we may not actually use simple regression, either.  We can correct for
            //bias toward the lower right quadrant by performing the regression twice, once inverted, and averaging
            //the two sets of coefficients
            _log.debug("forward: " + simpleRegressionResult[1] + "x + " + simpleRegressionResult[0] + " = y");
            if (!assumePositiveCorr)
                regressionResult = simpleRegressionResult;
            else
            {
                double[] regressionResultInverse = MatrixUtil.linearRegression(yValuesWithLowLeverage, xValuesWithLowLeverage);
                _log.debug("backward: " + regressionResultInverse[1] + "x + " + regressionResultInverse[0] + " = y");

                regressionResult = new double[2];

                double b1Inverse = 1.0 / regressionResultInverse[1];
                double b0Inverse = -regressionResultInverse[0] / regressionResultInverse[1];
                _log.debug("backward translated: " + b1Inverse + "x + " + b0Inverse + " = y");

                regressionResult[1] = (simpleRegressionResult[1] + b1Inverse) / 2;
                regressionResult[0] = (simpleRegressionResult[0] + b0Inverse) / 2;
                _log.debug("average: " + regressionResult[1] + "x + " + regressionResult[0] + " = y");
            }
        }
        else
        {
            try
            {
                regressionResult =
                        modalRegression(xValuesWithLowLeverage,
                                yValuesWithLowLeverage, degree);
            }
            catch (IOException e)
            {
                throw new RuntimeException(e);
            }
        }


        double[] residuals = new double[nLowLeverage];

        for (int i=0; i<xValuesWithLowLeverage.length; i++)
        {
            double predictedValue =
                    mapValueUsingCoefficients(regressionResult, xValuesWithLowLeverage[i]);
            residuals[i] = yValuesWithLowLeverage[i] - predictedValue;
        }
        double[] studentizedResiduals =
                BasicStatistics.studentizedResiduals(xValuesWithLowLeverage,
                        residuals, leverages);
        int[] indexesWithLowStudentizedResidual =
                selectIndexesWithLowAbsoluteSomething(studentizedResiduals, maxStudentizedResidual);
        _log.debug("selectFeaturesWithLowLeverageAndStudentizedResidual, result: " +
                indexesWithLowStudentizedResidual.length);

        int nResult = indexesWithLowStudentizedResidual.length;
        int[] result = new int[nResult];
        for (int i=0; i<nResult; i++)
            result[i] = indexesWithLowLeverage[indexesWithLowStudentizedResidual[i]];

        if (showCharts)
        {
            int maxX = 0;
            int minX = Integer.MAX_VALUE;
            double[] xValuesWithLowStudRes = new double[indexesWithLowStudentizedResidual.length];
            double[] yValuesWithLowStudRes = new double[indexesWithLowStudentizedResidual.length];

            for (int i =0; i<indexesWithLowStudentizedResidual.length; i++)
            {
                xValuesWithLowStudRes[i] = xValuesWithLowLeverage[indexesWithLowStudentizedResidual[i]];
                yValuesWithLowStudRes[i] = yValuesWithLowLeverage[indexesWithLowStudentizedResidual[i]];
            }

            for (double x : xValuesWithLowStudRes)
            {
                if (x > maxX)
                    maxX = (int) x;
                if (x < minX)
                    minX = (int) x;
            }
            PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
            pwsp.setName("LevAndStudRes");
            int numDotsOnChart = (maxX-minX+1) / 2;
            double[] chartXVals = new double[numDotsOnChart];
            double[] chartYVals = new double[numDotsOnChart];

            for (int j=0; j<numDotsOnChart; j++)
            {
                chartXVals[j] = minX + (2 * j);
                chartYVals[j] =
                        mapValueUsingCoefficients(regressionResult, chartXVals[j]);
            }
            pwsp.addData(xValuesWithLowStudRes, yValuesWithLowStudRes,
                        "Matches with low stud. res.");
            pwsp.addData(xValues, yValues, "all mass matches");
            pwsp.addData(chartXVals, chartYVals, "regression function");
            pwsp.setAxisLabels("X","Y");
            pwsp.displayInTab();
        }
        return result;
    }

    public static double[] symmetricalRegression(List<? extends Number> xset, List<? extends Number> yset)
    {
        double[] xsetArray = new double[xset.size()];
        double[] ysetArray = new double[xset.size()];
        for (int i=0; i<xsetArray.length; i++)
        {
            xsetArray[i] = xset.get(i).doubleValue();
            ysetArray[i] = yset.get(i).doubleValue();

        }
        return symmetricalRegression(xsetArray, ysetArray);
    }

    /**
     * If we assume a positive correlation between variables, we can correct for
     *bias toward the lower right quadrant by performing the regression twice, once inverted, and averaging
     *the two sets of coefficients
     * @param xValues
     * @param yValues
     * @return
     */
    public static double[] symmetricalRegression(double[] xValues, double[] yValues)
    {
        double[] simpleRegressionResult = MatrixUtil.linearRegression(xValues,
                yValues);
        double[] regressionResultInverse = MatrixUtil.linearRegression(yValues, xValues);
        _log.debug("backward: " + regressionResultInverse[1] + "x + " + regressionResultInverse[0] + " = y");

        double[]        regressionResult = new double[2];

        double b1Inverse = 1.0 / regressionResultInverse[1];
        double b0Inverse = -regressionResultInverse[0] / regressionResultInverse[1];
        _log.debug("backward translated: " + b1Inverse + "x + " + b0Inverse + " = y");

        regressionResult[1] = (simpleRegressionResult[1] + b1Inverse) / 2;
        regressionResult[0] = (simpleRegressionResult[0] + b0Inverse) / 2;
        _log.debug("average: " + regressionResult[1] + "x + " + regressionResult[0] + " = y");
        return regressionResult;
    }

    public static class AnovaResult
    {
        protected float fValue;
        protected float pValue;

        public float getFValue()
        {
            return fValue;
        }

        public void setFValue(float fValue)
        {
            this.fValue = fValue;
        }

        public float getPValue()
        {
            return pValue;
        }

        public void setPValue(float pValue)
        {
            this.pValue = pValue;
        }
    }
}
