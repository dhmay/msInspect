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
package org.fhcrc.cpl.viewer.amt;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.MatrixUtil;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Utilities related to AMT
 */
public class AmtUtilities
{
    private static Logger _log = Logger.getLogger(AmtUtilities.class);

    //name and version of the hydrophobicity algorithm currently used.
    //Values should be stored in normalized form, but it's still important
    //to keep track of what version of the algorithm they came from.
    public static final String HYDROPHOBICITY_ALGORITHM_NAME = "krokhin";
    public static final double HYDROPHOBICITY_ALGORITHM_VERSION = 3.0;

    //the search_engine string to use for pepXML files that contain AMT results
    public static final String AMT_SEARCH_ENGINE_CODE_FOR_PEPXML = "msInspect/AMT";

    public static void recordHydrophobicities(FeatureSet featureSet,
                                              Map<String,Double> regressionLine,
                                              int scanOrTimeMode)
    {
        for (Feature feature : featureSet.getFeatures())
            AmtExtraInfoDef.setObservedHydrophobicity(feature,
                    RegressionUtilities.predictYFromX(getSlopeFromRegressionLine(regressionLine),
                                               getInterceptFromRegressionLine(regressionLine),
                                               (scanOrTimeMode== ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN ?
                                                 feature.getScan() :
                                                 feature.getTime())));
    }

    public static void recordMS1Hydrophobicities(FeatureSet featureSet,
                                                 Map<String,Double> regressionLine,
                                                 int scanOrTimeMode)
    {
        recordHydrophobicities(featureSet, regressionLine, scanOrTimeMode);
    }

    /**
     * Pick numProteins proteins at "random" from a specified list of proteins.
     * It would sure be nice if we didn't have to load all the proteins in order to
     * pick a random sample, but since FastaLoader just provides an iterator, this is
     * about all we've got
     * @param fullProteinArray
     * @param numProteinsToLoad
     * @return
     */
    public static ArrayList<Protein> pickRandomProteins(ArrayList<Protein> fullProteinArray,
                                                        int numProteinsToLoad)
    {
        //cowardly refuse to pull out a random subset larger than the original set
        int fullProteinArraySize = fullProteinArray.size();
        if (fullProteinArraySize < numProteinsToLoad)
            return null;
        ArrayList<Protein> randomProteinArray = new ArrayList<Protein>(numProteinsToLoad);

        Random random = new Random();
        int[] chosenIndexes = new int[numProteinsToLoad];
        Arrays.fill(chosenIndexes,-1);
        for (int i=0; i<numProteinsToLoad; i++)
        {
            boolean foundOne = false;
            int nextIndex = -1;
            //prevent duplicate proteins
            while (!foundOne)
            {
                nextIndex = random.nextInt(fullProteinArraySize);
                if (Arrays.binarySearch(chosenIndexes,nextIndex) < 0)
                    foundOne = true;
            }
            randomProteinArray.add(fullProteinArray.get(nextIndex));
        }
        return randomProteinArray;
    }




    public static double[] getTimesForSortedScanArray(MSRun run,
                                                      int[] sortedScanArray)
    {
        double[] result = new double[sortedScanArray.length];
        if (run == null)
        {
            _log.error("Null run, can't calculate times for scans");
            return null;
        }
        int ms2IndexInRun = 0;
        MSRun.MSScan[] ms2Scans = run.getMS2Scans();
        for (int i = 0; i < sortedScanArray.length; i++)
        {
            int thisFeatureScanNum = sortedScanArray[i];
            //this relies on the fact that featuresForRegression is sorted by scan
            MSRun.MSScan thisFeatureScan = null;
            while (true)
            {
                if (ms2IndexInRun >= ms2Scans.length)
                {
                    //This will throw a NullPointerException down below when we
                    //try to access the retention time.
                    break;
                }
                int evalScanNum = ms2Scans[ms2IndexInRun].getNum();
                if (evalScanNum == thisFeatureScanNum)
                {
                    thisFeatureScan = ms2Scans[ms2IndexInRun];
                    break;
                }
                else if (evalScanNum > thisFeatureScanNum)
                {
                    _log.error("Error!  Read scan number " + evalScanNum + " while attempting to find scan number " + thisFeatureScanNum + ".");
                    _log.error("I.e., scan number " + thisFeatureScanNum + " was not found");
                    _log.error("Probably the wrong mzXml file was provided for pepXml file.");
                    break;
                }
                else
                    ms2IndexInRun++;
            }
            assert(thisFeatureScan != null);
            result[i] = thisFeatureScan.getDoubleRetentionTime();
        }
        return result;
    }

    /**
     * Workhorse method for relating time to hydrophobicity using linear regression,
     * can be used to relate time in terms of hydrophobicity or vice versa.
     *
     * optionally use robust regression, which takes a lot longer
     *
     * @return a hashmap containing three elements, whose keys are defined in constants: the slope,
     * intercept and sigma of the regression
     *
     * Note:  this regression predicts hydrophobicity in terms of time.  This relationship is
     * not reciprocal, due to all sorts of arcane weirdness related to linear regression.
     */
    public static Map<String,Double> calculateTimeHydroRelationshipEitherWay(
            Feature[] featuresForRegression, int scanOrTimeMode,
            boolean predictTimeFromHydro,
            boolean robustRegression)
    {


        Map<String,Double> result = new HashMap<String,Double>();

        //clone because we're sorting
        Feature[] featuresForRegressionClone =
                featuresForRegression.clone();

        if (featuresForRegressionClone.length == 0)
        {
            throw new RuntimeException("AmtUtilities.calculateTimeHydroRelationshipEitherWay: " +
            "regression failed, no eligible features to use for regression.");
        }

        Arrays.sort(featuresForRegressionClone, new Feature.ScanAscComparator());
        double[] hydrophobicities = new double[featuresForRegressionClone.length];
        double[] scansOrTimes = new double[hydrophobicities.length];

        for (int i=0; i<featuresForRegressionClone.length; i++)
        {
            Feature feature = featuresForRegressionClone[i];
            hydrophobicities[i] = calculateNormalizedHydrophobicity(
                    MS2ExtraInfoDef.getFirstPeptide(feature));
            scansOrTimes[i] =
                    (scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME) ?
                           feature.getTime() : feature.getScan();
        }


        try
        {
            double regressionLineSlope, regressionLineIntercept, sigma;


            double[] regressionResult;

            if (robustRegression)
                regressionResult = RegressionUtilities.robustRegression(
                    predictTimeFromHydro ? hydrophobicities : scansOrTimes,
                    predictTimeFromHydro ? scansOrTimes : hydrophobicities);
            else
                regressionResult = MatrixUtil.linearRegression(
                        predictTimeFromHydro? hydrophobicities : scansOrTimes,
                        predictTimeFromHydro? scansOrTimes : hydrophobicities);

            regressionLineSlope = regressionResult[1];
            regressionLineIntercept = regressionResult[0];
                        //not currently using sigma for anything, but I have a hunch it'll be crucial
            sigma = MatrixUtil.sigma(hydrophobicities, scansOrTimes, regressionResult);
            result.put(RegressionUtilities.REGRESSION_SLOPE_KEY, regressionLineSlope);
            result.put(RegressionUtilities.REGRESSION_INTERCEPT_KEY, regressionLineIntercept);
            result.put(RegressionUtilities.REGRESSION_SIGMA_KEY, sigma);
        }
        catch (Exception e)
        {
            _log.error("Failure to calculate hydrophobicity->time prediction",e);
        }

        return result;
    }

    /**
     * @return a HashMap containing three elements, whose keys are defined in constants: the slope,
     * intercept and sigma of the regression
     *
     * Note:  this regression predicts time in terms of hydrophobicity.  This relationship is
     * not reciprocal, due to all sorts of arcane weirdness related to linear regression.
     */
    public static Map<String,Double> calculateScanOrTimeHydrophobicityRelationship(
            Feature[] featuresForRegression, int scanOrTimeMode,
            boolean robustRegression)
    {
        return calculateTimeHydroRelationshipEitherWay(featuresForRegression, scanOrTimeMode,
                                                       false,
                                                       robustRegression);
    }

    /**
     * @return a HashMap containing three elements, whose keys are defined in constants: the slope,
     * intercept and sigma of the regression
     *
     * Note:  this regression predicts hydrophobicity in terms of time.  This relationship is
     * not reciprocal, due to all sorts of arcane weirdness related to linear regression.
     */
    public static Map<String,Double> calculateHydrophobicityScanOrTimeRelationship(
            Feature[] featuresForRegression, int scanOrTimeMode,
            boolean robustRegression)
    {
        return calculateTimeHydroRelationshipEitherWay(featuresForRegression, scanOrTimeMode,
                                                       true,
                                                       robustRegression);
    }

    /**
     * Given a peptide, and a line describing the relationship between hydrophobicity and
     * scan number, calculate the hydrophobicity of the peptide and return the predicted scan
     * for that hydrophobicity
     * @param intercept
     * @param peptide
     * @return
     */
    protected static double predictScanOrTime(double slope, double intercept, Peptide peptide)
    {
        return predictScanOrTime(slope, intercept,
                                 calculateNormalizedHydrophobicity(peptide));
    }

    /**
     * Given a hydrophobicity, and a line describing the relationship between hydrophobicity and
     * scan number, return the predicted scan for a given hydrophobicity.
     * @param slope
     * @param intercept
     * @param hydrophobicity
     * @return
     */
    public static double predictScanOrTime(double slope, double intercept, double hydrophobicity)
    {
        double predictedScanOrTime = RegressionUtilities.predictYFromX(slope, intercept, hydrophobicity);
        predictedScanOrTime = Math.max(predictedScanOrTime,0);
        return predictedScanOrTime;
    }

    public static double getSlopeFromRegressionLine(Map<String,Double> regressionLine)
    {
        return regressionLine.get(RegressionUtilities.REGRESSION_SLOPE_KEY);
    }

    public static double getInterceptFromRegressionLine(Map<String,Double> regressionLine)
    {
        return regressionLine.get(RegressionUtilities.REGRESSION_INTERCEPT_KEY);
    }

    /**
     * Predict the scan or time of a feature (depending on what the regression line represents)
     * @param regressionLine
     * @param feature
     * @return
     */
    public static double predictScanOrTime(Map<String,Double> regressionLine,
                                           Feature feature)
    {
        double slope = getSlopeFromRegressionLine(regressionLine);
        double intercept = getInterceptFromRegressionLine(regressionLine);
        return predictScanOrTime(slope, intercept, feature);
    }

    /**
     * Predict the scan or time of a feature (depending on what the regression parameters represent)
     * @param feature
     * @return
     */
    public static double predictScanOrTime(double slope, double intercept, Feature feature)
    {
        //this is lame.  Really peptide should have a constructor that doesn't require a protein
        Protein fakeProtein =
                new Protein("", MS2ExtraInfoDef.getFirstPeptide(feature).getBytes());
        Peptide currentPeptide = new Peptide(fakeProtein,0,fakeProtein.getBytes().length);
        //predict either the scan or the time, based on hydrophobicity.
        //the result will depend on which mode we're in
        return predictScanOrTime(slope, intercept, currentPeptide);
    }

    public static double convertDeltaScanOrTimeToHydro(Map<String,Double> regressionLineMap,
                                                       double scanOrTimeValueToConvert)
    {
        double slope = AmtUtilities.getSlopeFromRegressionLine(regressionLineMap);
        double intercept = AmtUtilities.getInterceptFromRegressionLine(regressionLineMap);
        return convertDeltaScanOrTimeToHydro(slope, intercept, scanOrTimeValueToConvert);
    }

    public static double convertDeltaScanOrTimeToHydro(double slope, double intercept,
                                                       double deltaScanOrTime)
    {
        return RegressionUtilities.predictYFromX(slope, intercept, deltaScanOrTime) -
               RegressionUtilities.predictYFromX(slope, intercept, 0);
    }





    /**
     * Create a peptide object from a peptide sequence.  Useful for finding mass, etc.
     * @param peptideSequence
     * @return
     */
    public static Peptide createPeptideFromSequence(String peptideSequence)
    {
        //this is lame.  Really peptide should have a constructor that doesn't require a protein
        Protein fakeProtein = new Protein("",peptideSequence.getBytes());
        return new Peptide(fakeProtein,0,fakeProtein.getBytes().length);
    }



    //Methods related to calculation and normalization of peptide hydrophobicity.
    //It's all about convenience.


    public static double calculateRawHydrophobicity(String peptideSequence)
    {
        return calculateRawHydrophobicity(createPeptideFromSequence(peptideSequence));
    }

    /**
     * Cover method for Peptide.getHydrophobicity / getHydrophobicity3.
     * Use this to switch between versions of the hydrophobicity calculator
     * @param peptide
     * @return
     */
    public static double calculateRawHydrophobicity(Peptide peptide)
    {
        //return peptide.getHydrophobicity();
        return peptide.getHydrophobicity3();
    }

    public static double normalizeHydrophobicity(double inputHydrophobicity)
    {
        return HydrophobicityNormalizer.normalize(inputHydrophobicity,
                                                  HYDROPHOBICITY_ALGORITHM_NAME,
                                                  HYDROPHOBICITY_ALGORITHM_VERSION);
    }

    public static double calculateNormalizedHydrophobicity(Peptide peptide)
    {
        return normalizeHydrophobicity(calculateRawHydrophobicity(peptide));
    }

    public static double calculateNormalizedHydrophobicity(String peptideSequence)
    {
        return normalizeHydrophobicity(calculateRawHydrophobicity(peptideSequence));
    }

    /**
     * Warning about this one: if there's more than one peptide associated with
     * the feature, it'll take the first 
     * @param feature
     * @return
     */
    public static double calculateNormalizedHydrophobicity(Feature feature)
    {
        return calculateNormalizedHydrophobicity(MS2ExtraInfoDef.getFirstPeptide(feature));
    }




    /**
     * Assumes features have times populated if necessary
     *
     * Given a line that should predict time from  hydrophobicity, indicate by how much
     * (in units of scans or seconds) each feature appears to deviate from that line
     * @param features
     * @param scanOrTimeMode
     * @return
     */
    public static void addHydrophobicityToFeatures(Feature[] features,
                                                   int scanOrTimeMode,
                                                   double[] timeToHCoefficients)
    {
        try
        {
            for (Feature feature : features)
            {
                AmtExtraInfoDef.setObservedHydrophobicity(feature,
                        RegressionUtilities.mapValueUsingCoefficients(timeToHCoefficients,
                            (scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME) ?
                                feature.getTime() : feature.getScan()));
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("Error adding hydrophobicity observations to features",e);
        }
    }


    /**
     * Clusterable for clustering based on feature mz and hydrophobicity.
     * Allows access to the original Feature
     */
    public static class FeatureMzHydroClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMzHydroClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mz;
        }
        public double getDimension2Value()
        {
            return AmtExtraInfoDef.getObservedHydrophobicity(parentFeature);
        }
    }

    /**
     * Clusterable for clustering based on feature mass and hydrophobicity.
     * Allows access to the original Feature
     */
    public static class FeatureMassHydroClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMassHydroClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mass;
        }
        public double getDimension2Value()
        {
            return AmtExtraInfoDef.getObservedHydrophobicity(parentFeature);
        }
    }
}
