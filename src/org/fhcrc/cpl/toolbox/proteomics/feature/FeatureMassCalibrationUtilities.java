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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;


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
public class FeatureMassCalibrationUtilities
{
    protected static Logger _log = Logger.getLogger(FeatureMassCalibrationUtilities.class);

    public FeatureMassCalibrationUtilities()
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
        List<float[]> massArrays = new ArrayList();
        for (Feature[] featureArray : featureArrays)
        {
            float[] massArray = new float[featureArray.length];
            for (int i=0; i<featureArray.length; i++)
            {
                massArray[i] = useMz ? featureArray[i].getMz() : featureArray[i].getMass();
            }
            massArrays.add(massArray);
        }
        return MassCalibrationUtilities.plotMassDefectDeviation(massArrays,
                theoreticalMassWavelength, shouldDisplayPPMLine, ppmForLine);
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
                    FeatureMassCalibrationUtilities.filterFeaturesByMassDefectDeviation(features,
                            initialMassFilterPPM);
            Feature[] initialFilteredFeaturesByScan = new Feature[initialFilteredFeatures.length];
            System.arraycopy(initialFilteredFeatures, 0, initialFilteredFeaturesByScan, 0,
                    initialFilteredFeatures.length);
            Arrays.sort(initialFilteredFeaturesByScan, new Feature.ScanAscComparator());
            Pair<Integer, Pair<Double, Double>>[] partitionParameters =
                    FeatureMassCalibrationUtilities.calculateWavelengthsAndOffsetsMultiplePartitions(
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
        float[] masses = new float[features.length];
        for (int i=0; i<features.length; i++)
            masses[i] = features[i].getMass();
        float[] newMasses = MassCalibrationUtilities.calibrateMassesSinglePartition(masses, maxPairsForLeverageCalc, 
                theoreticalMassWavelength);
        for (int i=0; i<features.length; i++)
        {
            features[i].setMass(masses[i]);
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
        float[] allMasses = new float[features.length];
        for (int i=0; i<features.length; i++)
            allMasses[i] = features[i].getMass();
        return MassCalibrationUtilities.calculateWavelengthAndOffset(allMasses, maxPairsForLeverageCalc,
                theoreticalMassWavelength);
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
                MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH);
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
                MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH);
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
            if (Math.abs(MassCalibrationUtilities.calculateMassDefectDeviationPPM(feature.getMass(),
                    theoreticalMassWavelength)) <= maxDeviationPPM)
                resultList.add(feature);
        }

        return resultList.toArray(new Feature[resultList.size()]);
    }

}
