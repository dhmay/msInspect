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
package org.fhcrc.cpl.toolbox;

import org.apache.log4j.Logger;

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


    public static double mapValueUsingCoefficients(double[] coefficients, double valueToMap)
    {
        double result = 0;
        for (int j=0; j<coefficients.length; j++)
            result += coefficients[j] * Math.pow(valueToMap, j);
        return result;
    }

    
}
