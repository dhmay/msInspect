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
package org.fhcrc.cpl.viewer.align;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseMatcher;

import java.io.*;

/**
 * Aligner using quantile regression
 */
public class QuantileRegressionAligner extends Aligner
{
    private static Logger _log = Logger.getLogger(QuantileRegressionAligner.class);

    protected int nonlinearMappingPolynomialDegree =
            AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;



    protected double[] modalRegressionCoefficients = null;

    public QuantileRegressionAligner()
    {
    }

    /**
     * Cover method. In the absence of filename information, use a filename that's
     * likely to be unique
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public double[] alignPairs(Pair<Double,Double>[] pairs, int maxValueToWarp)
    {
        return alignPairs(pairs, maxValueToWarp, "length" + pairs.length);
    }
    
    /**
     * Generic method for creating a warping based on pairs of values.
     * Create temp files with names based on tempFileNameStart, in case the user
     * wants to peruse them
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public double[] alignPairs(Pair<Double,Double>[] pairs, double maxValueToWarp,
                               String tempFileNameStart)
    {
        PrintWriter pw = null;
        double[] result = null;
        if (pairs.length < 4)
            throw new RuntimeException("QuantileRegressionAligner.alignPairs: at least 4 pairs are necessary for " +
                    "alignment, only " + pairs.length + " were provided");
        _log.debug("alignPairs, #pairs=" + pairs.length);
        try
        {
            double[] baseTimes = new double[pairs.length];
            double[] toAlignTimes = new double[pairs.length];
            
            for (int j=0; j<pairs.length; j++)
            {
                Pair<Double,Double> matchedPair = pairs[j];
                baseTimes[j] = matchedPair.first;
                toAlignTimes[j] = matchedPair.second;
            }
            try
            {
                modalRegressionCoefficients = RegressionUtilities.modalRegression(baseTimes,
                        toAlignTimes, nonlinearMappingPolynomialDegree);
                _log.debug("Regression complete");

                result = new double[(int) (maxValueToWarp+1)];
                for (int j=0; j<=maxValueToWarp; j++)
                {
                    result[j] =
                        RegressionUtilities.mapValueUsingCoefficients(modalRegressionCoefficients, j);

                }

            }
            catch (IOException e)
            {
                e.printStackTrace(System.err);
                //if we failed, it could be because of timeout, or because of a missing
                //quantreg package, or...?
                //todo: move text to bundle
                throw new RuntimeException("Failure calling R for modal regression.  R may have timed out.\n" +
                        "This may also be because the required \"quantreg\" package is not installed.\n"+
                        "  To install this package, run the following lines in R:\n" +
                        "source(\"http://bioconductor.org/biocLite.R\")\n" +
                        "biocLite(c(\"quantreg\"))");
            }

        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Failed alignment: " + x.getMessage(), x);
            return null;
        }
        finally
        {
            if (pw != null)
                pw.close();
        }

        return result;
    }


    public double[] getModalRegressionCoefficients()
    {
        return modalRegressionCoefficients;
    }


    public int getNonlinearMappingPolynomialDegree()
    {
        return nonlinearMappingPolynomialDegree;
    }

    public void setNonlinearMappingPolynomialDegree(int nonlinearMappingPolynomialDegree)
    {
        this.nonlinearMappingPolynomialDegree = nonlinearMappingPolynomialDegree;
    }
}
