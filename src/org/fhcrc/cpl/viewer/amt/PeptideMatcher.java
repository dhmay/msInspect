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

import java.util.*;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.apache.log4j.Logger;

public class PeptideMatcher
{
    private static Logger _log = Logger.getLogger(PeptideMatcher.class);

    //number of times the standard deviation by which a peptide's intensity must differ from the
    //protein's average intensity in order for it to be tossed out
//Not implemented
//    protected static final double INTENSITY_OUTLIER_STD_DEV_MULTIPLIER = 1;


    //Two different modes for time tolerance when matching peptides.
    //Fixed:  always use the same time tolerance
    //Adaptive:
    //  For unobserved peptides, UNOBSERVED_PEPTIDE_GRADIENT_PERCENT_TOLERANCE percent of the gradient
    //  For peptides observed once before, deltaTime
    //  For peptides observed multiple times before, a multiple of the standard deviation of the observations
    public static final int PEPTIDE_TIME_TOLERANCE_MODE_FIXED=0;
    public static final int PEPTIDE_TIME_TOLERANCE_MODE_ADAPTIVE=1;
    public static final int DEFAULT_PEPTIDE_TIME_TOLERANCE_MODE=PEPTIDE_TIME_TOLERANCE_MODE_FIXED;

    public static final int UNOBSERVED_PEPTIDE_GRADIENT_PERCENT_TOLERANCE=5;

    public static final Feature.MassAscComparator featureMassAscComp = new Feature.MassAscComparator();


    /**
     * Calculate the standard deviation of the absolute value of hydrophobicity differences in
     * peptides identified between two pepxml files
     * @param ms2Features1
     * @param ms2Features2
     * @param scanOrTimeMode
     * @param run1
     * @param run2
     * @param minPeptideProphet
     * @return
     */
    public static double calculateHydroDiffStandardDeviation(FeatureSet ms2Features1,
                                                             FeatureSet ms2Features2,
                                                             int scanOrTimeMode,
                                                             MSRun run1, MSRun run2,
                                                             double minPeptideProphet,
                                                             boolean robustRegression)
    {
        double[] hydrophobicityDifferences =
                calculateHydrophobicityDifferences(ms2Features1, ms2Features2, scanOrTimeMode,
                                                   run1, run2, minPeptideProphet,
                                                   robustRegression);
//for (double diff : hydrophobicityDifferences) System.err.println(" " + diff);
        return BasicStatistics.standardDeviation(hydrophobicityDifferences);
    }

    /**
     * Calculate the median of the absolute value of hydrophobicity differences in
     * peptides identified between two pepxml files
     * @param ms2Features1
     * @param ms2Features2
     * @param scanOrTimeMode
     * @param run1
     * @param run2
     * @param minPeptideProphet
     * @return
     */
    public static double calculateHydroDiffMedian(FeatureSet ms2Features1,
                                                             FeatureSet ms2Features2,
                                                             int scanOrTimeMode,
                                                             MSRun run1, MSRun run2,
                                                             double minPeptideProphet,
                                                             boolean robustRegression)
    {
        double[] hydrophobicityDifferences =
                calculateHydrophobicityDifferences(ms2Features1, ms2Features2, scanOrTimeMode,
                                                   run1, run2, minPeptideProphet, robustRegression);
        return BasicStatistics.median(hydrophobicityDifferences);
    }

    /**
     * return an array of the absolute values of differences between hydrophobicities of
     * peptides observed in both runs
     * @param ms2Features1
     * @param ms2Features2
     * @param scanOrTimeMode
     * @param run1
     * @param run2
     * @param minPeptideProphet
     * @return
     */
    public static double[] calculateHydrophobicityDifferences(FeatureSet ms2Features1,
                                                             FeatureSet ms2Features2,
                                                             int scanOrTimeMode,
                                                             MSRun run1, MSRun run2,
                                                             double minPeptideProphet,
                                                             boolean robustRegression)
    {
        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet((float)minPeptideProphet);
        ms2Features1 = ms2Features1.filter(sel);
        ms2Features2 = ms2Features2.filter(sel);


        ms2Features1.populateTimesForMS2Features(run1);
        ms2Features2.populateTimesForMS2Features(run2);



        //Calculate the time-hydro relationships for each run
        Map<String, Double> lineMap =
                AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(ms2Features1.getFeatures(),
                        scanOrTimeMode, robustRegression);
        double slope1 = lineMap.get(RegressionUtilities.REGRESSION_SLOPE_KEY);
        double intercept1 = lineMap.get(RegressionUtilities.REGRESSION_INTERCEPT_KEY);

        lineMap = AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(ms2Features2.getFeatures(),
                scanOrTimeMode, robustRegression);
        double slope2 = lineMap.get(RegressionUtilities.REGRESSION_SLOPE_KEY);
        double intercept2 = lineMap.get(RegressionUtilities.REGRESSION_INTERCEPT_KEY);

        //Keep track of the peptides found in set 1 and their features
        HashMap<String, Feature> peptidesFromFeatureSet1 = new HashMap<String, Feature>();
        for (Feature set1Feature : ms2Features1.getFeatures())
            peptidesFromFeatureSet1.put(MS2ExtraInfoDef.getFirstPeptide(set1Feature), set1Feature);

        ArrayList<Double> hydroDifferencesList = new ArrayList<Double>();

        //For features in the two sets with the same peptide, calculate their hydrophobicities
        //and store the difference between them
        for (Feature set2Feature : ms2Features2.getFeatures())
        {
            if (peptidesFromFeatureSet1.containsKey(MS2ExtraInfoDef.getFirstPeptide(set2Feature)))
            {
                //get time/scan for set 1 feature
                double time1;
                Feature set1Feature = peptidesFromFeatureSet1.get(MS2ExtraInfoDef.getFirstPeptide(set2Feature));
                if (scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN)
                    time1 = set1Feature.getScan();
                else
                    time1 = run1.getMS2Scans()[run1.getIndexForMS2ScanNum(set1Feature.getScan())].getDoubleRetentionTime();

                //get time/scan for set 2 feature
                double time2;
                if (scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN)
                    time2 = set2Feature.getScan();
                else
                    time2 = run2.getMS2Scans()[run2.getIndexForMS2ScanNum(set2Feature.getScan())].getDoubleRetentionTime();

                //translate time to hydrophobicity
                double hydro1 =
                        RegressionUtilities.predictYFromX(
                                slope1, intercept1, time1);
                double hydro2 =
                        RegressionUtilities.predictYFromX(
                                slope2, intercept2, time2);
                //store the absolute value of the difference
                hydroDifferencesList.add(Math.abs(hydro1 - hydro2));
            }
        }

        //Convert the ArrayList to an array.  It would sure be nice if there were a better way
        double[] hydrophobicityDifferences = new double[hydroDifferencesList.size()];
        for (int i = 0; i < hydrophobicityDifferences.length; i++)
        {
            hydrophobicityDifferences[i] = Math.abs(hydroDifferencesList.get(i));
        }

        return hydrophobicityDifferences;
    }



    //Utility methods for getting what we need from a map representing a regression line


    /**
     * For a given feature, try all possible mass modifications to match against a featureset.
     * That is, for each modification, try the low value for _all matching residues_, then
     * try the high value.  Does _not_ try combinations of low and high.
     * So the number of times this method will be run is at most 2^(number of var mods).
     * If it finds a match sooner, it'll return.
     * The return value is a cloned version of the matched feature, with the peptide added.
     * Note: the mass table is altered during execution of this method
     * @param feature
     * @param currentPeptide
     * @param ms1Features
     * @param comp
     * @param varModArray
     * @param startIndex
     * @param massTable
     * @return a cloned version of the matched feature, with the peptide added
     */
    public static Feature recursivelyMatchAllMods(Feature feature, Peptide currentPeptide,
                                                  Feature[] ms1Features,
                                                  ClusteringFeatureSetMatcher featureSetMatcher,
                                                  Comparator<Feature> comp,
                                                  MS2Modification[] varModArray,
                                                  int startIndex, double[] massTable)
    {
        if (startIndex >= varModArray.length)
        {
            feature.setMass((float) currentPeptide.getMass(massTable));

            List<Feature> matches = null;// featureSetMatcher.findMatchingFeatures(ms1Features, feature,
                                           //                               comp);
            Feature result = null;
            if (matches != null && !matches.isEmpty())
            {
                //clone the feature and mark it with the peptide
                //TODO: here we toss out all matches but the first.  Is that ok?
                result = (Feature) matches.get(0).clone();
                MS2ExtraInfoDef.setSinglePeptide(result,new String(currentPeptide.getChars()));
            }
            return result;
        }

        Feature matchedFeature = recursivelyMatchAllMods(feature, currentPeptide, ms1Features,
                                                     featureSetMatcher,
                                                     comp, varModArray, startIndex+1, massTable);
        if (matchedFeature != null)
            return matchedFeature;

        MS2Modification varMod = varModArray[0];
        int indexToModify = (int) varMod.getAminoAcid().charAt(0);
        float massDiff = varMod.getMassDiff();
        massTable[indexToModify] = massTable[indexToModify] + massDiff;
        return recursivelyMatchAllMods(feature, currentPeptide, ms1Features,
                                       featureSetMatcher,
                                       comp, varModArray, startIndex+1, massTable);
    }

}
