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
package org.fhcrc.cpl.viewer.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.gui.util.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.RegressionUtilities;
import org.fhcrc.cpl.toolbox.*;


import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;


/**
 * Utilities for recalibrating and characterizing the mass calibration of sets of
 * features, based on an approach that uses theoretical peptide mass clusters.
 *
 * Peptides should have masses with certain predictable qualities.  In particular, the
 * distances between masses should cluster around multiples of a certain number, called
 * the "mass wavelength".  This number would be 1 (the mass of a proton in a Carbon atom),
 * except that different elements' nucleon masses are different, due to the effect
 * of the mass defect.  So the clustering isn't perfect, but it is pretty good.  You can
 * determine how the deviation from predicted cluster changes as the distance between
 * masses increases.
 *
 * This deviation can be used to calibrate masses in a given run, and it can be used to
 * toss out wildly deviating masses from runs known to be calibrated correctly.
 *
 * For more information, and for the theoretical basis of pretty much everything in here,
 * check out:
 * "Analytical model of peptide mass cluster centres with applications"
 * Wolski, Farrow, et al.
 * http://www.proteomesci.com/content/4/1/18
 */
public class MassCalibrationUtilities
{
    protected static Logger _log = Logger.getLogger(MassCalibrationUtilities.class);

    //the weighted average of all peptides' mass defect per dalton, at monoisotopic mass.
    //That is, for the distribution of elements within the average peptide,
    //the weighted average of the mass defect, divided by the atomic weight
    //of the light (common) version of the element

//    public static final double THEORETICAL_MASS_WAVELENGTH = 1.000482;
//public static final double THEORETICAL_MASS_WAVELENGTH = 1.000467;
public static final double DEFAULT_THEORETICAL_MASS_WAVELENGTH = 1.000476; //empirically derived



    //The intercept for the line with slope THEORETICAL_MASS_WAVELENGTH, that
    //maps mass to mass cluster deviation.
    //I have no great understanding of why this is not zero, but it is not zero.
    public static final double THEORETICAL_MASS_DEFECT_INTERCEPT = .034;
//public static final double THEORETICAL_MASS_DEFECT_INTERCEPT = 0.013;
//public static final double THEORETICAL_MASS_DEFECT_INTERCEPT = .0;


    //maximum number of pairs for use in leverage calculation, when doing
    //robust regression to determine a run's characteristics.
    //This is purely a performance consideration.  Higher is better, but much higher
    //than this default value is really, really expensive
    public static final int DEFAULT_MAX_PAIRS_FOR_LEVERAGE_CALC = 500000;

    //Parameters related to robust regression
    public static final double DEFAULT_MAX_LEVERAGE_NUMERATOR_FOR_REGRESSION = 4.0;
    public static final double DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_REGRESSION = 2.3;

    //default maximum deviation from theoretical mass cluster, for use in filtering.
    //This default value is derived from the Wolski et al. paper.  I think we could
    //probably go tighter.
    public static int DEFAULT_MAX_DEVIATION_PPM = 200;


    public MassCalibrationUtilities()
    {
    }

    /**
     * Convenience method for a single feature set
     * @param features
     * @return
     */
    public static PanelWithScatterPlot plotMassDefectDeviation(Feature[] features,
                                                            double theoreticalMassWavelength,
                                                            boolean shouldDisplayPPMLine,
                                                            double ppmForLine)
    {
        List<Feature[]> featureList = new ArrayList<Feature[]>(1);
        featureList.add(features);
        return plotMassDefectDeviation(featureList, theoreticalMassWavelength,
                shouldDisplayPPMLine, ppmForLine,true);
    }

    /**
     * Plot the mass defect deviations of multiple featuresets, in different colors.
     * Optionally displays lines indicating deviations of ppmForLine parts per million
     * from theoretical cluster center.
     * @param featureArrays
     * @return
     */
    public static PanelWithScatterPlot plotMassDefectDeviation(List<Feature[]> featureArrays,
                                                            double theoreticalMassWavelength,
                                                            boolean shouldDisplayPPMLine,
                                                            double ppmForLine,
                                                            boolean useMz)
    {
        PanelWithScatterPlot pwsp = new PanelWithScatterPlot();

        int setNum=1;
        for (Feature[] features : featureArrays)
        {
            double[][] scatterPlotData = new double[2][features.length];

            for (int i=0; i<features.length; i++)
            {
                double mass = features[i].getMass();
                if (useMz)
                    mass = features[i].getMz();

                double massOffset = (mass -  THEORETICAL_MASS_DEFECT_INTERCEPT)
                        % theoreticalMassWavelength;

                             
                if (massOffset >= 0.5)
                    massOffset -= 1;
                scatterPlotData[0][i] = (float) (mass);
                scatterPlotData[1][i] = (float) massOffset;
            }
            pwsp.addData(scatterPlotData[0], scatterPlotData[1], "Set " + setNum++);
        }

        if (shouldDisplayPPMLine)
        {
            double minMass = Double.MAX_VALUE;
            double maxMass = Double.MIN_VALUE;
            for (Feature feature : featureArrays.get(0))
            {
                float featureMass = feature.getMass();
                if (useMz)
                    featureMass = feature.getMz();
                if (featureMass < minMass)
                    minMass = featureMass;
                if (featureMass > maxMass)
                    maxMass = featureMass;
            }
            int numPoints = (int) (maxMass - minMass + 1);
            double[] linesXValues = new double[2 * numPoints];
            double[] linesYValues = new double[2 * numPoints];
            for (int i=0; i<numPoints; i++)
            {
                linesXValues[i] = i;
                linesXValues[numPoints + i] = i;
                linesYValues[i] = ((double)i) * ppmForLine / 1000000.0;
                linesYValues[numPoints+i] = -linesYValues[i];
            }
            pwsp.addData(linesXValues, linesYValues, ppmForLine + "PPM Deviation");
        }

        pwsp.setAxisLabels("Mass","Mass Cluster Deviation");
        return pwsp;
    }

    /**
     * master method for mass correction.
     *
     * Masses can be corrected separately for different ranges of scans within the file,
     * since calibration can drift throughout the run.
     *
     * @param features
     * @param maxPairsForLeverageCalc
     * @param numPartitions
     * @param theoreticalMassWavelength
     * @param showCharts
     */
    public static void calibrateMasses(Feature[] features,
                                       int maxPairsForLeverageCalc,
                                       int numPartitions,
                                       double theoreticalMassWavelength,
                                       int initialMassFilterPPM,
                                       boolean showCharts)
    {
        Feature[] featuresByScan = new Feature[features.length];
        System.arraycopy(features, 0, featuresByScan, 0, features.length);
        Arrays.sort(featuresByScan, new Feature.ScanAscComparator());

        //if features for mass calibration calculation must first undergo an initial mass filter:
        //1.  do the initial filter
        //2.  calculate all the partition parameters
        //3.  apply those parameters to ALL the features
        if (initialMassFilterPPM > 0)
        {
            Feature[] initialFilteredFeatures =
                    MassCalibrationUtilities.filterFeaturesByMassDefectDeviation(features,
                            initialMassFilterPPM);
            Feature[] initialFilteredFeaturesByScan = new Feature[initialFilteredFeatures.length];
            System.arraycopy(initialFilteredFeatures, 0, initialFilteredFeaturesByScan, 0,
                    initialFilteredFeatures.length);
            Arrays.sort(initialFilteredFeaturesByScan, new Feature.ScanAscComparator());
            Pair<Integer, Pair<Double, Double>>[] partitionParameters =
                    MassCalibrationUtilities.calculateWavelengthsAndOffsetsMultiplePartitions(
                            initialFilteredFeaturesByScan, maxPairsForLeverageCalc, numPartitions,
                            theoreticalMassWavelength, showCharts);

            calibrateMasses(featuresByScan, partitionParameters, theoreticalMassWavelength, showCharts);
        }
        else
        {
            //Calculate the calibration parameters for all partitions at the outset
            Pair<Integer, Pair<Double, Double>>[] partitionParameters =
                    calculateWavelengthsAndOffsetsMultiplePartitions(
                            featuresByScan, maxPairsForLeverageCalc, numPartitions,
                            theoreticalMassWavelength, showCharts);
            calibrateMasses(featuresByScan, partitionParameters, theoreticalMassWavelength, showCharts);
        }
    }

    public static void calibrateMasses(Feature[] featuresByScan,
                                       Pair<Integer, Pair<Double, Double>>[] partitionParameters,
                                       double theoreticalMassWavelength,
                                       boolean showCharts)
    {
        int numPartitions = partitionParameters.length;

        int paramIndex = 0;
        Pair<Double,Double> currentParameters = null;
        double currentWavelength =0;
        double currentOffset =0;
        int nextPartitionStartScan = -1;
        List<Feature> currentPartitionFeatureList = new ArrayList<Feature>();

        List<Feature[]> partitionFeaturesList = new ArrayList<Feature[]>();

        for (Feature feature : featuresByScan)
        {
            if (feature.getScan() >= nextPartitionStartScan)
            {
                currentParameters = partitionParameters[paramIndex++].getValue();
                currentWavelength = currentParameters.first;
                currentOffset = currentParameters.second;
                ApplicationContext.infoMessage("Partition " + (paramIndex) + ": wavelength=" +
                           currentWavelength + ", offset=" + currentOffset);                
                if (numPartitions > paramIndex)
                    nextPartitionStartScan = partitionParameters[paramIndex].getKey();
                else
                {
                    //last partition
                    nextPartitionStartScan = Integer.MAX_VALUE;
                }
                partitionFeaturesList.add(currentPartitionFeatureList.toArray(new Feature[currentPartitionFeatureList.size()]));
                currentPartitionFeatureList = new ArrayList<Feature>();
            }

            float newMass = feature.getMass();
            newMass += feature.getMass() * (theoreticalMassWavelength - currentWavelength) - currentOffset;

            feature.setMass(newMass);
            if (feature.getCharge() > 0)
                feature.updateMz();
            currentPartitionFeatureList.add(feature);

        }
        //last partition
        partitionFeaturesList.add(currentPartitionFeatureList.toArray(new Feature[currentPartitionFeatureList.size()]));
        
        if (showCharts)
        {
            PanelWithScatterPlot pwsp =
                    plotMassDefectDeviation(partitionFeaturesList, theoreticalMassWavelength, false, 0,true);
            pwsp.setName("after");
            pwsp.displayInTab();
        }

    }


    /**
     * Returns an array of pairs.  The first item in the pair is a _scan number_ indicating
     * where the region begins.  The second item is a Double,Double pair containing the
     * wavelength and intercept to be applied in the region.
     * @param featuresByScan
     * @param maxPairsForLeverageCalc
     * @param numPartitions
     * @param theoreticalMassWavelength
     * @param showCharts
     * @return
     */
    public static Pair<Integer, Pair<Double,Double>>[]
        calculateWavelengthsAndOffsetsMultiplePartitions(Feature[] featuresByScan,
                                       int maxPairsForLeverageCalc,
                                       int numPartitions,
                                       double theoreticalMassWavelength,
                                       boolean showCharts)
    {


        int partitionSize = featuresByScan.length / numPartitions;

        List<Feature[]> partitionFeaturesList = new ArrayList<Feature[]>();

        Pair<Integer, Pair<Double,Double>>[] result =
                (Pair<Integer, Pair<Double,Double>>[]) new Pair[numPartitions];

        for (int i=0; i<numPartitions; i++)
        {
            int partitionStart = i * partitionSize;
            int partitionEnd = ((i+1) * partitionSize) - 1;

            if (i==numPartitions-1)
                partitionEnd = featuresByScan.length-1;

            int thisPartitionSize = (partitionEnd - partitionStart) - 1;

            Feature[] featuresInPartition = new Feature[thisPartitionSize];
            System.arraycopy(featuresByScan, partitionStart,
                             featuresInPartition, 0, thisPartitionSize);
            _log.debug("Partition " + (i+1) + ": scans " + partitionStart + "-" + partitionEnd +
                       ".  " + featuresInPartition.length + " features");
            partitionFeaturesList.add(featuresInPartition);
        }
        if (showCharts)
        {
            PanelWithScatterPlot pwsp =
                    plotMassDefectDeviation(partitionFeaturesList, theoreticalMassWavelength, false, 0,true);
            pwsp.setName("Before");
            pwsp.displayInTab();
        }

        for (int i=0; i<partitionFeaturesList.size(); i++)
        {
            Feature[] featuresInPartition = partitionFeaturesList.get(i);
            Pair<Double,Double> wavelengthAndOffset =
                    calculateWavelengthAndOffset(featuresInPartition, maxPairsForLeverageCalc,
                                                 theoreticalMassWavelength);

            result[i] =
                    new Pair<Integer, Pair<Double,Double>>(
                            featuresByScan[i * partitionSize].getScan(),
                            wavelengthAndOffset);

        }
        return result;
    }


    /**
     * Calibrate masses in a single partition.
     *
     * First we calculate the mass wavelength, and the offset from theoretical mass
     * intercept. Then we adjust each mass.
     * @param features
     * @param maxPairsForLeverageCalc
     * @param theoreticalMassWavelength
     */
    protected static void calibrateMassesSinglePartition(Feature[] features,
                                                       int maxPairsForLeverageCalc,
                                                       double theoreticalMassWavelength)
    {
        Pair<Double,Double> wavelengthAndOffset =
                calculateWavelengthAndOffset(features, maxPairsForLeverageCalc, theoreticalMassWavelength);
        double wavelength = wavelengthAndOffset.first;
        double offset = wavelengthAndOffset.second;
        double[] deltaMasses = new double[features.length];
        double[] deltaMassesPPM = new double[features.length];

        for (int i=0; i<features.length; i++)
        {
            Feature feature = features[i];
            float newMass = feature.getMass();
            newMass += feature.getMass() * (theoreticalMassWavelength - wavelength) - offset;

            deltaMasses[i] = newMass - feature.getMass();
            deltaMassesPPM[i] = deltaMasses[i] * (1000000.0/feature.getMass());
            feature.setMass(newMass);
            if (feature.getCharge() > 0)
                feature.updateMz();
        }
        if (_log.isDebugEnabled())
        {                       
            double meanMassChange = BasicStatistics.mean(deltaMasses);
            double meanMassChangePPM = BasicStatistics.mean(deltaMassesPPM);

            for (int i=0; i<deltaMasses.length; i++)
                deltaMasses[i] = Math.abs(deltaMasses[i]);
            double meanAbsoluteMassChange = BasicStatistics.mean(deltaMasses);
            _log.debug("Average mass change: " + meanMassChange + " (" + meanMassChangePPM + " ppm)");
            _log.debug("Average _absolute_ mass change: " + meanAbsoluteMassChange);             
        }
    }

    /**
     * Use robust regression techniques to map mass to mass offset from theoretical cluster
     * @param features
     * @param maxPairsForLeverageCalc more is better, but more costs more time
     * @return a Pair containing slope, intercept of the regression line.
     */
    public static Pair<Double,Double> calculateWavelengthAndOffset(Feature[] features,
                                                              int maxPairsForLeverageCalc,
                                                              double theoreticalMassWavelength)
    {
        int n = features.length;
        double[] allMasses = new double[n];
        for (int i=0; i<features.length; i++)
            allMasses[i] = features[i].getMass();

        //calculate the mass defect wavelength.  That is, the distance between clusters.
        //This can be thought of as the slope of the regression line
        double massDefectWavelength =
                calculateMassDefectWavelength(features, maxPairsForLeverageCalc, theoreticalMassWavelength);

        //calculate mass defect intercept.  Once we've got the wavelength, this is just
        //a matter of taking the average of (the feature masses (mod wavelength)) and
        //subtracting the theoretical intercept.
        //At least according to the Wolski paper

        double[] massDefectIntercepts = new double[features.length];
        for (int i=0; i<features.length; i++)
        {
            double intercept =
//                    (features[i].getMass()) % massDefectWavelength;

                    (features[i].getMass() %
                          massDefectWavelength) - THEORETICAL_MASS_DEFECT_INTERCEPT;
            if (intercept >= 0.5)
                intercept = intercept - 1.0;
            massDefectIntercepts[i] = intercept;
        }
//ScatterPlotDialog spd = new ScatterPlotDialog(allMasses, massDefectIntercepts, "intercepts");
//spd.setVisible(true);

        double meanMassDefectIntercept = BasicStatistics.mean(massDefectIntercepts);

        _log.debug("wavelength and offset: " + massDefectWavelength + ", " + meanMassDefectIntercept);
        return new Pair<Double,Double>(massDefectWavelength,meanMassDefectIntercept);
    }

    /**
     * Calculate the wavelength; that is, the distance between mass clusters.  I.e., the slope
     * of the line that maps mass to mass defect.  
     * First, find all pairs of features (in reality, this is capped at maxPairsForLeverageCalc,
     * due to performance considerations -- the problem is just too large.  We ensure a representative
     * sample).  Calculate the distances between each pair, and the wavelength for each pair.
     * Then, toss out the ones with high leverage.
     * Then, perform a linear regression to map distance to wavelength.  Calculate all the residuals,
     * and toss out the ones with high studentized residuals.  Do the regression again and return
     * the slope, which is the wavelength.  The intercept has no meaning in this context.
     */
    public static double calculateMassDefectWavelength(Feature[] features,
                                                       int maxPairsForLeverageCalc,
                                                       double theoreticalMassWavelength)
    {
        Pair<double[],double[]> distancesAndWavelengthsForRegression =
                findDistancesAndWavelengthResidualsForRegression(features, maxPairsForLeverageCalc,
                        theoreticalMassWavelength);
        _log.debug("assembled distances, wavelengths: " +
                distancesAndWavelengthsForRegression.first.length);

        double[] distancesForRegression = distancesAndWavelengthsForRegression.first;
        double[] wavelengthResidualsForRegression = distancesAndWavelengthsForRegression.second;



//ScatterPlotDialog spd = new ScatterPlotDialog(distancesForRegression, wavelengthsForRegression, "after leverage cutoff");
//spd.setVisible(true);
        double[] regressionResult =
                MatrixUtil.linearRegression(distancesForRegression,
                        wavelengthResidualsForRegression);

//ScatterPlotDialog spd = new ScatterPlotDialog(distancesForRegression, wavelengthsForRegression, "points after stud res cutoff");
//double[] regressionLineXVals = new double[2000];
//double[] regressionLineYVals = new double[2000];
//for (int i=0; i<regressionLineXVals.length; i++)
//{
//    regressionLineXVals[i] = i;
//    regressionLineYVals[i] = regressionResult[1] * i + regressionResult[0];
//}
//spd.addData(regressionLineXVals, regressionLineYVals, "First Regression line");
//spd.setVisible(true);
        _log.debug("first reg: wavelength: " + regressionResult[1]);

        int n = distancesForRegression.length;

        double[] massDefectResiduals = new double[n];

        _log.debug("calculateRegressionLine: Pairs after leverage cutoff: " + n);

        for (int i=0; i<n; i++)
        {
            double predictedMassDefect =
                    RegressionUtilities.predictYFromX(
                            regressionResult[1],
                            regressionResult[0],
                            distancesForRegression[i]);
            massDefectResiduals[i] =
                    Math.abs(wavelengthResidualsForRegression[i] - predictedMassDefect);
        }

        double[] studentizedResiduals =
                BasicStatistics.studentizedResiduals(distancesForRegression,
                        massDefectResiduals);

        boolean[] pairsInSecondRegression = new boolean[n];

        int numInSecondRegression = 0;

        for (int i=0; i<n; i++)
        {
            //todo: parameterize
            if (studentizedResiduals[i] < DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_REGRESSION)
            {
                pairsInSecondRegression[i] = true;
                numInSecondRegression++;
            }
            else
                pairsInSecondRegression[i] = false;
        }

        double[] distancesForSecondRegression = new double[numInSecondRegression];
        double[] wavelengthResidualsForSecondRegression = new double[numInSecondRegression];

        int index=0;
        for (int i=0; i<n; i++)
        {
            if (pairsInSecondRegression[i])
            {
                distancesForSecondRegression[index] = distancesForRegression[i];
                wavelengthResidualsForSecondRegression[index] =
                        wavelengthResidualsForRegression[i];
                index++;
            }
        }

        _log.debug("calculateRegressionLine: pairs after stud.res. cutoff: " +
                numInSecondRegression);



        regressionResult =
                MatrixUtil.linearRegression(distancesForSecondRegression,
                        wavelengthResidualsForSecondRegression);

double[] secondRegressionLineXVals = new double[2000];
double[] secondRegressionLineYVals = new double[secondRegressionLineXVals.length];
for (int i=0; i<secondRegressionLineXVals.length; i++)
{
    secondRegressionLineXVals[i] = i;
    secondRegressionLineYVals[i] = regressionResult[1] * i + regressionResult[0];
}
//ScatterPlotDialog spd2 =
//        new ScatterPlotDialog(secondRegressionLineXVals, secondRegressionLineYVals, "Second Regression line");
//spd2.addData(distancesForSecondRegression, wavelengthResidualsForSecondRegression, "points after stud res cutoff");
//spd2.setVisible(true);

        _log.debug("After second regression, wavelength: " + (theoreticalMassWavelength + regressionResult[1]));

        return theoreticalMassWavelength + regressionResult[1];
    }

    /**
     * First, find all pairs of features... or, actually, a representative sample of 
     * maxPairsForLeverageCalc pairs, due to performance considerations -- the problem 
     * is just too large.  
     * We ensure a fairly random sample, by figuring out by what factor there are more 
     * pairs than we can handle (x) and then taking every xth pair.  This won't ensure
     * that we actually get maxPairs... pairs, only that we get at least maxPairs.../2.  
     * Anything more complicated gets too intensive, though, due to the sheer number of
     * pairs... even choosing, say, 10 million random numbers takes way too long.  
     *
     * Next, calculate the distances between each pair, and the wavelength for each pair.
     *
     * Then, calculate the leverage of each (based on distance) and toss out the ones 
     * with high leverage.
     *
     * Return what's left
     *
     * @param features
     * @param maxPairsForLeverageCalc
     * @return
     */
    protected static Pair<double[],double[]> findDistancesAndWavelengthResidualsForRegression(
            Feature[] features, int maxPairsForLeverageCalc, double theoreticalMassWavelength)
    {
        int numFeatures = features.length;

        int numPairs = 0;
        
        numPairs = (numFeatures*(numFeatures-1))/2;

        ApplicationContext.infoMessage("Total pairs: " + numPairs);

        int pairInterval = 1;
        if (numPairs > maxPairsForLeverageCalc)
        {
            ApplicationContext.infoMessage("Too many pairs (" + numPairs + "), reducing...");

            pairInterval = (Math.round(numPairs / maxPairsForLeverageCalc));
            if (pairInterval == 1)
                pairInterval = 2;

            //todo: this is stupid.  I should be able to calculate this.  Why am I so lazy??!
            _log.debug("Counting reduced pairs");
            int numKeptPairs=0;
            int allPairsIndex=0;
            for (int i=0; i<numFeatures; i++)
            {
                for (int j=0; j<numFeatures; j++)
                {
                    if (i<j && (allPairsIndex % pairInterval) == 0)
                    {
                        numKeptPairs++;
                    }
                    allPairsIndex++;
                }
            }

            numPairs = numKeptPairs;
            _log.debug("Done Counting reduced pairs: " + numPairs);
        }

        float[] pairWavelengthResiduals = new float[numPairs];
        float[] pairDistances = new float[numPairs];

        int allPairsIndex=0;
        int keptPairsIndex=0;

        for (int i=0; i<numFeatures; i++)
        {
            for (int j=0; j<numFeatures; j++)
            {
                if (i<j && (allPairsIndex % pairInterval) == 0)
                {
                    pairDistances[keptPairsIndex] =
                            Math.abs(features[i].getMass() - features[j].getMass());

                    double thisPairWavelengthResidual =
                            pairDistances[keptPairsIndex] %
                                    theoreticalMassWavelength;
                    if (thisPairWavelengthResidual >= 0.5)
                        thisPairWavelengthResidual -= 1;
                    pairWavelengthResiduals[keptPairsIndex++] = (float) thisPairWavelengthResidual;
                }
                allPairsIndex++;
            }
        }

//ScatterPlotDialog spd = new ScatterPlotDialog(pairDistances, pairWavelengths, "all points");
//spd.setVisible(true);

        double maxLeverage =
                DEFAULT_MAX_LEVERAGE_NUMERATOR_FOR_REGRESSION / (double) numPairs;
        float[] leverages = BasicStatistics.leverages(pairDistances);

        ApplicationContext.setMessage("Calculated leverages");

        List<Integer> pairsInFirstRegression = new ArrayList<Integer>();
        int numPairsInRegression = 0;
        for (int i=0; i<numPairs; i++)
        {
            if (leverages[i] < maxLeverage)
            {
                pairsInFirstRegression.add(i);
                numPairsInRegression++;
            }
        }

        double[] regressionDistances = new double[numPairsInRegression];
        double[] regressionWavelengthResiduals = new double[numPairsInRegression];
        int index=0;
        for (int pairIndex : pairsInFirstRegression)
        {
                regressionDistances[index] = pairDistances[pairIndex];
                regressionWavelengthResiduals[index] = pairWavelengthResiduals[pairIndex];
                index++;
        }

        return new Pair<double[],double[]>(regressionDistances, regressionWavelengthResiduals);
    }

    /**
     * Not the entry point for the --filter command
     * @param featureSet
     */
    public static void filterFeatureSetByMassDefectDeviation(FeatureSet featureSet,
                                                             double theoreticalMassWavelength)
    {
        featureSet.setFeatures(filterFeaturesByMassDefectDeviation(featureSet.getFeatures(),
                theoreticalMassWavelength));
    }

    /**
     * Not the entry point for the --filter command
      * @param featuresList
     * @return
     */
    public static List<Feature> filterFeaturesByMassDefectDeviation(List<Feature> featuresList,
                                                                double maxDeviationPPM)
    {
        Feature[] features = featuresList.toArray(new Feature[featuresList.size()]);
        Feature[] filtered = filterFeaturesByMassDefectDeviation(features,
                maxDeviationPPM,
                DEFAULT_THEORETICAL_MASS_WAVELENGTH);
        List<Feature> result = new ArrayList<Feature>(filtered.length);
        for (Feature feature : filtered)
            result.add(feature);
        return result;
    }

    /**
     * Not the entry point for the --filter command
      * @param features
     * @return
     */
    public static Feature[] filterFeaturesByMassDefectDeviation(Feature[] features,
                                                                double maxDeviationPPM)
    {
        return filterFeaturesByMassDefectDeviation(features,
                maxDeviationPPM,
                DEFAULT_THEORETICAL_MASS_WAVELENGTH);
    }

    /**
     * Not the entry point for the --filter command
     * @param features
     * @param maxDeviationPPM
     * @return
     */
    public static Feature[] filterFeaturesByMassDefectDeviation(Feature[] features,
                                                    double maxDeviationPPM,
                                                    double theoreticalMassWavelength)
    {
        List<Feature> resultList = new ArrayList<Feature>();
        for (Feature feature : features)
        {
            if (Math.abs(calculateMassDefectDeviationPPM(feature, theoreticalMassWavelength)) <= maxDeviationPPM)
                resultList.add(feature);
        }

        return resultList.toArray(new Feature[resultList.size()]);
    }

    /**
     * Convenience method for PPM conversion
     * @param feature
     * @return
     */
    public static double calculateMassDefectDeviationPPM(Feature feature,
                                                         double theoreticalMassWavelength)
    {
        double massDefectDaltons = calculateMassDefectDeviation(feature, theoreticalMassWavelength);
        return ((massDefectDaltons * 1000000.0) / feature.getMass());
    }

    /**
     * Calculate the deviation from the closest theoretical mass cluster.
     * Assumes masses are already calibrated correctly
     * @param feature
     * @return
     */
    public static double calculateMassDefectDeviation(Feature feature,
                                                      double theoreticalMassWavelength)
    {
        double mass = feature.getMass();

        //mass offset is the mass mod the weighted average of all peptides'
        //actual mass per dalton
        double predictedMassDefectDeviation =
                (mass - THEORETICAL_MASS_DEFECT_INTERCEPT) % (theoreticalMassWavelength);

        //massOffset is now divided... half the peptides will be above the average,
        //and those will be just above 0.  The other half will be below the average,
        //and those will be just below 1.
        //This step moves the ones just below 1 so they're just below 0
        if (predictedMassDefectDeviation >= 0.5)
            predictedMassDefectDeviation -= 1;

        return predictedMassDefectDeviation;
    }

}
