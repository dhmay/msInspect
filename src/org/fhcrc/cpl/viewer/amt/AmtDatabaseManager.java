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
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;

/**
 * Manages AMT databases -- prunes outliers of various kinds
 *
 */
public class AmtDatabaseManager
{
    static Logger _log = Logger.getLogger(AmtDatabaseManager.class);

    //default parameters for median alignment
    public static final int DEFAULT_MIN_PEPTIDES_FOR_ALIGNMENT_REGRESSION=30;
    public static final int DEFAULT_MIN_OBSERVATIONS_FOR_ALIGNMENT_REGRESSION=3;

    //parameters for database self-alignment
    public static final int DEFAULT_MATCHING_DEGREE_FOR_DB_ALIGNMENT = AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;


    /**
     * Remove entries with only one observation, where that observation differs from
     * predicted H by more than significantHydroDifference
     * @param amtDB
     * @param significantHydroDifference
     */
    public static void removePredictedHOutliers(AmtDatabase amtDB,
                                                float significantHydroDifference)
    {
        int excludedEntryCount = 0;
        for (AmtPeptideEntry peptideEntry : amtDB.getEntries())
        {
            //if there's only one observation for this peptide entry, and we're conditionally
            //excluding single-observation entries...
            if (peptideEntry.getNumObservations() < 2)
            {
                //if the single obseration is too far off from the predicted hydrophobicity...
                if (Math.abs(peptideEntry.getPredictedHydrophobicity() -
                        peptideEntry.getMedianObservedHydrophobicity()) >
                        significantHydroDifference)
                {
                    amtDB.removeEntry(peptideEntry.getPeptideSequence());
                    excludedEntryCount++;
                }
            }
        }
        ApplicationContext.infoMessage("Removed " + excludedEntryCount +
                " entries with suspicious observed H values");
    }

    /**
     * Adjust AMT database entries to account for the effect of acrylamide
     * @param amtDB Altered in place
     * @param fromAcrylamideToNot If true, adjust acrylamide-containing entries to remove acrylamide's effect.
     * If false, adjust non-acrylamide-containing entries to add the effect of acrylamide
     * @return
     */
    public static AmtDatabase adjustEntriesForAcrylamide(AmtDatabase amtDB,
                                                         boolean fromAcrylamideToNot,
                                                         boolean showCharts)
    {
        List<Double> oldHVals = new ArrayList<Double>();
        List<Double> newHVals = new ArrayList<Double>();

        for (AmtPeptideEntry peptideEntry : amtDB.getEntries())
        {
            for (AmtPeptideEntry.AmtPeptideObservation obs : peptideEntry.getObservations())
            {
                double oldHVal = obs.getObservedHydrophobicity();
                double newHydrophobicity;
                if (fromAcrylamideToNot)
                    newHydrophobicity = AcrylamideHydrophobicityAdjuster.removeAcrylamideEffect(
                            peptideEntry.getPeptideSequence(), oldHVal);
                else
                    newHydrophobicity = AcrylamideHydrophobicityAdjuster.adjustHForAcrylamide(
                            peptideEntry.getPeptideSequence(), oldHVal);
                obs.setObservedHydrophobicity(newHydrophobicity);
                oldHVals.add(oldHVal);
                newHVals.add(newHydrophobicity);
            }
            peptideEntry.recalculateStats();
        }

        if (showCharts)
        {
            ScatterPlotDialog spd = new ScatterPlotDialog(oldHVals, newHVals, "New vs. Old hydrohpobicity values");
            spd.setVisible(true);
        }                                      

        return amtDB;
    }


    protected static Set<String> intersection(Set<String> overlapSet1, Set<String> overlapSet2)
    {
        Set<String> intersSet = new HashSet<String>(overlapSet1);
        intersSet.retainAll(overlapSet2);
        return intersSet;
    }

    /**
     *
     **
     * @param amtDB
     * @param minMatchedPeptides
     */
    public static AmtDatabase alignAllRunsUsingCommonPeptides(AmtDatabase amtDB,
                                           int minMatchedPeptides, int matchingdegree,
                                           float maxLeverageNumerator, float maxStudentizedResidual, boolean showCharts)
    {

        Map<AmtRunEntry, Set<String>> peptidesInRuns =
                new HashMap<AmtRunEntry, Set<String>>(amtDB.numRuns());

        int numAllRuns = amtDB.numRuns();

        ApplicationContext.infoMessage("Choosing best runs to start with...");

        for (AmtRunEntry run : amtDB.getRuns())
        {
            Set<String> peptidesThisRun = new HashSet<String>();
            for (AmtPeptideEntry peptideEntry1 : amtDB.getPeptideEntriesForRun(run))
                peptidesThisRun.add(peptideEntry1.getPeptideSequence());
            peptidesInRuns.put(run, peptidesThisRun);
        }

        Pair<AmtRunEntry, AmtRunEntry> closestRuns = null;
        int numMatchedBetweenClosest = 0;
        Set<String> closestPeptideOverlap = null;
        for (AmtRunEntry run1 : amtDB.getRuns())
        {
            Set<String> run1Peptides = peptidesInRuns.get(run1);
            for (AmtRunEntry run2 : amtDB.getRuns())
            {
                Set<String> run2Peptides = peptidesInRuns.get(run2);

                //if it can't possibly be the best one, skip it
                if (run2Peptides.size() < numMatchedBetweenClosest)
                    continue;

                Set<String> matchedPeptides = intersection(run1Peptides, run2Peptides);
                if (matchedPeptides.size() > numMatchedBetweenClosest)
                {
                    closestPeptideOverlap = matchedPeptides;
                    numMatchedBetweenClosest = matchedPeptides.size();
                    closestRuns = new Pair<AmtRunEntry, AmtRunEntry>(run1, run2);
                }
            }
        }

        if (numMatchedBetweenClosest < minMatchedPeptides)
            throw new RuntimeException("No two runs had enough in common to align");
        ApplicationContext.infoMessage("Adding initial runs...");

        AmtDatabase result = new AmtDatabase();
        AmtRunEntry baseRun = closestRuns.first;
        addRun(result, amtDB, baseRun);


        Map<AmtRunEntry, double[]> runTHMappingCoefficientMap =
                new HashMap<AmtRunEntry, double[]>();

        runTHMappingCoefficientMap.put(closestRuns.second,
                alignRun(result, amtDB, closestRuns.second, closestPeptideOverlap, matchingdegree,
                        maxLeverageNumerator, maxStudentizedResidual, showCharts));
        addRun(result, amtDB, closestRuns.second);


        Set<String> peptidesAddedToDatabase = new HashSet<String>();
        peptidesAddedToDatabase.addAll(peptidesInRuns.get(baseRun));
        peptidesAddedToDatabase.addAll(peptidesInRuns.get(closestRuns.second));

        Collection<AmtRunEntry> runsToAlign = new ArrayList<AmtRunEntry>();
        for (AmtRunEntry run : amtDB.getRuns())
        {
            if (run != baseRun && run != closestRuns.second)
                runsToAlign.add(run);
        }

        ApplicationContext.infoMessage("Adding the rest of the runs...");

        //keep track of the peptides already in the database that also exist in all runs
        Map<AmtRunEntry, Set<String>> databasePeptidesInRuns =
                new HashMap<AmtRunEntry, Set<String>>(amtDB.numRuns());
        for (AmtRunEntry runEntry : runsToAlign)
            databasePeptidesInRuns.put(runEntry, intersection(peptidesAddedToDatabase, peptidesInRuns.get(runEntry)));

        PanelWithChart runningChartPanel = null;
        while (runsToAlign.size() > 0)
        {
            ApplicationContext.setMessage("Runs remaining: " + runsToAlign.size() +
                    " (" + Rounder.round(100.0 * (float)(numAllRuns-runsToAlign.size())/(float)numAllRuns, 2) +
                    "% complete)");

            AmtRunEntry nextRun = null;
            int maxOverlapSize = -1;
            for (AmtRunEntry possibleNextRun : runsToAlign)
            {
                Set<String> possibleNextRunOverlap = databasePeptidesInRuns.get(possibleNextRun);
                possibleNextRunOverlap.addAll(intersection(peptidesAddedToDatabase, peptidesInRuns.get(possibleNextRun)));

                if (possibleNextRunOverlap.size() > maxOverlapSize)
                {
                    maxOverlapSize = possibleNextRunOverlap.size();
                    nextRun = possibleNextRun;
                }
            }

            if (maxOverlapSize < minMatchedPeptides)
            {
                ApplicationContext.infoMessage("Unable to align remaining " + runsToAlign.size() + " runs, excluding them");
                break;
            }

            runTHMappingCoefficientMap.put(nextRun,
                    alignRun(result, amtDB, nextRun, databasePeptidesInRuns.get(nextRun), matchingdegree,
                            maxLeverageNumerator, maxStudentizedResidual, showCharts));
            addRun(result, amtDB, nextRun);
            runsToAlign.remove(nextRun);
            
            //cleanup
            databasePeptidesInRuns.remove(nextRun);


        }

        if (showCharts)
        {
            double maxTimeAnyRun = 0;
            for (AmtRunEntry run : amtDB.getRuns())
            {
                for (AmtPeptideEntry.AmtPeptideObservation obs : amtDB.getObservationsForRun(run))
                    if (obs.getTimeInRun() > maxTimeAnyRun)
                        maxTimeAnyRun = obs.getTimeInRun();
            }

            PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
            pwsp.setName("Realignment of DB runs");
            for (double[] thMapCoefficients : runTHMappingCoefficientMap.values())
            {
                int numDotsOnChart = ((int) maxTimeAnyRun+1) /2;
                double[] mappingXVals = new double[numDotsOnChart];
                double[] mappingYVals = new double[numDotsOnChart];
                for (int j=0; j<numDotsOnChart; j++)
                {
                    mappingXVals[j] = 2 * j;
                    mappingYVals[j] = RegressionUtilities.mapValueUsingCoefficients(thMapCoefficients, j);
                }
                pwsp.addData(mappingXVals, mappingYVals, "");
            }

//            spd.addData(reverseRegressionLineXVals, reverseRegressionLineYVals, "reverse regression line");
//            spd.addData(regressionLineXVals, regressionLineYVals, "regression line");
//            spd.addData(lowStudResMs1Times, lowStudResMs2H, "low studentized residual");

            pwsp.setAxisLabels("MS1 time","AMT Hydrophobicity");
            pwsp.displayInTab();
        }

        return result;
    }

    protected static double[] alignRun(AmtDatabase newDatabase, AmtDatabase sourceDatabase,
                          AmtRunEntry runToAdd, Set<String> peptideOverlap, int matchingDegree,
                          float maxLeverageNumerator, float maxStudentizedResidual, boolean showCharts)
    {
        List<Pair<Feature, Feature>> peptideMatches = new ArrayList<Pair<Feature, Feature>>();

        for (String commonPeptide : peptideOverlap)
        {
            AmtPeptideEntry newDatabaseEntry = newDatabase.getEntry(commonPeptide);
            AmtPeptideEntry.AmtPeptideObservation runToAddObservation =
                    sourceDatabase.getEntry(commonPeptide).getObservationForRun(runToAdd);
            Feature newDatabaseFeature = new Feature();
            AmtExtraInfoDef.setObservedHydrophobicity(newDatabaseFeature,
                    newDatabaseEntry.getMedianObservedHydrophobicity());
            Feature runToAddFeature = new Feature();
            runToAddFeature.setTime((float) runToAddObservation.getTimeInRun());
            peptideMatches.add(new Pair<Feature,Feature>(runToAddFeature, newDatabaseFeature));
        }
        AmtDatabaseMatcher matcher = new AmtDatabaseMatcher();
        matcher.setMaxRegressionStudentizedResidual(maxStudentizedResidual);
        matcher.setMaxRegressionLeverageNumerator(maxLeverageNumerator);
        matcher.setBuildCharts(true);

        float maxTimeForMap = 0;
        for (AmtPeptideEntry.AmtPeptideObservation obs : sourceDatabase.getObservationsForRun(runToAdd))
        {
            if (obs.getTimeInRun() > maxTimeForMap)
                maxTimeForMap = (float) obs.getTimeInRun() + 1;
        }
        double[] timeHydMapCoefficients = matcher.calculateTHMapCoefficientsWithMatchedFeatures(
                peptideMatches.toArray((Pair<Feature,Feature>[]) new Pair[peptideMatches.size()]),
                matchingDegree);
        runToAdd.setTimeHydMapCoefficients(timeHydMapCoefficients);
        for (AmtPeptideEntry.AmtPeptideObservation newRunObs : sourceDatabase.getObservationsForRun(runToAdd))
        {
            double newHydrophobicity =
                    RegressionUtilities.mapValueUsingCoefficients(timeHydMapCoefficients,
                            newRunObs.getTimeInRun());
            newRunObs.setObservedHydrophobicity(newHydrophobicity);
        }

        if (showCharts)
        {
            if (MultiChartDisplayPanel.getSingletonInstance().getNumCharts() > 0)
                MultiChartDisplayPanel.getSingletonInstance().removeAllCharts();
            PanelWithChart pwc = new PanelWithChart(matcher.getTimeHydrophobicityMappingChart());
            pwc.setName("Alignment for run " + newDatabase.numRuns());
            pwc.displayInTab();
        }

        return timeHydMapCoefficients;
    }



    protected static void addRun(AmtDatabase newDatabase, AmtDatabase sourceDatabase,
                          AmtRunEntry runToAdd)
    {
            Map<MS2Modification, MS2Modification> modificationMap =
                    newDatabase.addRunEntry(runToAdd);


            for (AmtPeptideEntry peptideEntry : sourceDatabase.getPeptideEntriesForRun(runToAdd))
            {
                AmtPeptideEntry newEntry = null;
//                    newEntry.setPeptideSequence(peptideEntry.getPeptideSequence());
//                    newEntry.setPredictedHydrophobicity(AmtUtilities.calculateNormalizedHydrophobicity(peptideEntry.getPeptideSequence()));
                for (AmtPeptideEntry.AmtPeptideModificationStateEntry modState : peptideEntry.getModificationStateEntries())
                {
                    AmtPeptideEntry.AmtPeptideObservation obsForRun =
                            modState.getObservationForRun(runToAdd);
                    if (obsForRun != null)
                    {
                        if (newEntry == null)
                        {
                            newEntry = AmtPeptideEntry.createEntryFromObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                        }
                        else
                            newEntry.addObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                    }
                }
                if (newEntry != null)
                    newDatabase.addObservationsFromEntry(newEntry, modificationMap);
            }
            ApplicationContext.infoMessage("Adding run " + runToAdd.getPepXmlFilename());
    }

    /**
     * Pick obserations with a low leverage value.
     *
     * Since there's a linear map from time to hydrophobicity, it doesn't matter
     * whether we base the leverage on time or hydrophobicity. Using hydrophobicity
     *
     * @param observations
     * @return
     */
    protected static AmtPeptideEntry.AmtPeptideObservation[]
        pickObservationsWithLowLeverage(AmtPeptideEntry.AmtPeptideObservation[] observations)
    {
        List<AmtPeptideEntry.AmtPeptideObservation> resultList =
                new ArrayList<AmtPeptideEntry.AmtPeptideObservation>();
        int n = observations.length;
        double[] xValues = new double[n];
        for (int i = 0; i < n; i++)
            xValues[i] = observations[i].getObservedHydrophobicity();
        double[] leverages = BasicStatistics.leverages(xValues);
//todo: parameterize?
        double maxLeverage = AmtDatabaseBuilder.DEFAULT_MAX_LEVERAGE_NUMERATOR / (double) n;
        List<Double> timesForRegressionList = new ArrayList<Double>();
        for (int i = 0; i < n; i++)
        {
            if (leverages[i] < maxLeverage)
            {
                resultList.add(observations[i]);
            }
        }
        return resultList.toArray(new AmtPeptideEntry.AmtPeptideObservation[0]);
    }

    /**
     * Given a set of Strings and an array of candidate sets, find the candidate set
     * with the most overlap with the base.  Ignore any sets in setsToIgnore
     *
     * This would be trivial if HashSets actually implemented the stuff that they
     * should
     * @param baseSet
     * @param candidateSets
     * @param setsToIgnore
     * @return
     */
    protected static int findSetWithMostOverlap(Set<String> baseSet,
                                                Set<String>[] candidateSets,
                                                List<Integer> setsToIgnore)
    {
        int maxOverlap = -1;
        int result = 0;
        for (int i=0; i<candidateSets.length; i++)
        {
            if (setsToIgnore.contains(i))
                continue;
            int currentOverlap = 0;
            for (String baseString : baseSet)
                if (candidateSets[i].contains(baseString))
                    currentOverlap++;
            if (currentOverlap > maxOverlap)
            {
                maxOverlap = currentOverlap;
                result = i;
            }
        }
        return result;
    }

    /**
     * Throw out observations that are above a given multiple of the standard deviation
     * of all observations for that peptide.
     *
     * REMEMBER: for low degrees of freedom (few peptide observations), standard deviation
     * gets wonky.  You really want to compare against a T value instead.
     * But according to Yan, a cutoff of 3 standard deviations is still pretty darn
     * safe.  For 3 peptide observations, it's still equivalent to tossing out <5% of
     * observations. It gets more conservative as you go up in df.
     *
     * So, IDEALLY, the parameter would be percentToThrowAway, rather than
     * hydroStdDevMultipleCutoff, but that's complicated, and this is safe enough.
     *
     * @param amtDB
     * @param hydroStdDevMultipleCutoff
     * @return
     */
    public static AmtPeptideEntry.AmtPeptideObservation[] removeHydrophobicityOutliers(
            AmtDatabase amtDB, double hydroStdDevMultipleCutoff)
    {
        List<AmtPeptideEntry.AmtPeptideObservation> resultList =
                new ArrayList<AmtPeptideEntry.AmtPeptideObservation>();

        //for debug logging
        int[] outliersPerRun = new int[amtDB.numRuns()];

        for (AmtPeptideEntry peptideEntry : amtDB.getEntries())
        {
            double deltaHydroCutoff = hydroStdDevMultipleCutoff *
                                      peptideEntry.getHydrophobicityStandardDeviation();
            if (peptideEntry.getNumObservations() < 3)
                continue;
            for (AmtPeptideEntry.AmtPeptideObservation obs :
                    peptideEntry.getObservations())
            {
                double deltaHydro =
                        Math.abs(peptideEntry.getMedianObservedHydrophobicity() -
                                 obs.getObservedHydrophobicity());
                if (deltaHydro > deltaHydroCutoff)
                {
                    resultList.add(obs);
                    //for debug logging
                    outliersPerRun[amtDB.getSequenceForRun(obs.getRunEntry()) - 1]++;
                    peptideEntry.removeObservation(obs);
                }
            }
        }
        
        //for debug logging
        if (_log.isDebugEnabled())
        {
            int totalOutliers = 0;
            for (int i=0; i<amtDB.numRuns(); i++)
            {
                _log.debug("Run " + (i+1) + ": " + outliersPerRun[i] + " outliers");
                totalOutliers += outliersPerRun[i];
            }
            _log.debug("total outliers: " + totalOutliers);
        }

        return resultList.toArray(new AmtPeptideEntry.AmtPeptideObservation[0]);
    }

    /**
     * This is a utility method to determine minimum and maximum masses for matching individual
     * features by mass.  Returns min and then max masses, in an array.
     * @param features
     * @param massMatchDeltaMass
     * @param massMatchDeltaMassType
     * @return
     */
    public static float[][] findMinAndMaxMassesForMatch(Feature[] features,
                                                        float massMatchDeltaMass,
                                                        int massMatchDeltaMassType)
    {
        float[] minMassesForMatch = new float[features.length];
        float[] maxMassesForMatch = new float[features.length];

        for (int j=0; j<features.length; j++)
        {
            float ms1FeatureMass = features[j].getMass();
            minMassesForMatch[j] = ms1FeatureMass - massMatchDeltaMass;
            maxMassesForMatch[j] = ms1FeatureMass + massMatchDeltaMass;
            if (massMatchDeltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_PPM)
            {
                float absDeltaMass = (massMatchDeltaMass * minMassesForMatch[j] / 1000000f);
                minMassesForMatch[j] = ms1FeatureMass - absDeltaMass;
                maxMassesForMatch[j] = ms1FeatureMass + absDeltaMass;
            }
        }
        return new float[][] { minMassesForMatch,maxMassesForMatch };
    }

    /**
     * Return a database containing only entries from the subset of runs
     * from the passed-in AMT database that match, by mass only, at least
     * minMassMatchPercent of their entries to ms1Features.
     *
     * Add the best runs first, and stop if db gets too big
     * @param amtDatabase
     * @param ms1FeaturesOrig
     * @param minMassMatchPercent
     * @param massMatchDeltaMass
     * @param massMatchDeltaMassType
     * @param maxNewDBEntries
     * @param modificationsInMS1
     * @param showCharts
     * @return
     */
    public static AmtDatabase removeRunsWithoutMassMatches(AmtDatabase amtDatabase,
                                                       Feature[] ms1FeaturesOrig,
                                                       int minMassMatchPercent,
                                                       float massMatchDeltaMass,
                                                       int massMatchDeltaMassType,
                                                       int maxNewDBEntries,
                                                       int maxNewDBRuns,
                                                       MS2Modification[] modificationsInMS1,
                                                       boolean showCharts)
    {
        Feature[] ms1FeaturesCopy = new Feature[ms1FeaturesOrig.length];
        System.arraycopy(ms1FeaturesOrig, 0, ms1FeaturesCopy, 0, ms1FeaturesCopy.length);
        Arrays.sort(ms1FeaturesCopy, new Feature.MassAscComparator());

        AmtDatabaseMatcher amtDatabaseMatcher = new AmtDatabaseMatcher();

        AmtDatabase result = new AmtDatabase();
        float[] percentMatchedPerRun = new float[amtDatabase.numRuns()];
        AmtRunEntry[] runEntries = amtDatabase.getRuns();

        float[][] minAndMaxMassesForMatch = 
                findMinAndMaxMassesForMatch(ms1FeaturesCopy, massMatchDeltaMass, massMatchDeltaMassType);
        float[] minMassesForMatch = minAndMaxMassesForMatch[0];
        float[] maxMassesForMatch = minAndMaxMassesForMatch[1];

        final Map<AmtRunEntry, Float> runMassMatchPercentMap =
                new HashMap<AmtRunEntry, Float>();

        int numPassing = 0;
        for (int i=0; i<runEntries.length; i++)
        {
            if (i % (Math.max(runEntries.length/10,1)) == 0)
            {
                ApplicationContext.setMessage(Rounder.round(((double) i * 100.0 / (double) runEntries.length),0) +
                        "% done evaluating runs; passing: " + numPassing + " / " + i +
                        " (" + Rounder.round(((double) numPassing  * 100.0) / (double) i, 0) + "%)...");
            }
            AmtRunEntry runEntry = runEntries[i];
            AmtDatabaseFeatureSetGenerator featureGen =
                    new AmtDatabaseFeatureSetGenerator(amtDatabase);
            FeatureSet runFeatureSet =
                    featureGen.createFeatureSetForRun(runEntry, modificationsInMS1);
            Feature[] runFeatures = runFeatureSet.getFeatures();
            Arrays.sort(runFeatures, new Feature.MassAscComparator());

            int numMassMatchedFeatures = 0;

            int ms1FeaturesIndex = 0;
            int runFeaturesIndex = 0;
            while (runFeaturesIndex < runFeatures.length &&
                   ms1FeaturesIndex < ms1FeaturesCopy.length)
            {
                while (runFeaturesIndex < runFeatures.length &&
                       runFeatures[runFeaturesIndex].getMass() <=
                        minMassesForMatch[ms1FeaturesIndex])
                    runFeaturesIndex++;
                if (runFeaturesIndex >= runFeatures.length)
                    break;
                if (runFeatures[runFeaturesIndex].getMass() <=
                        maxMassesForMatch[ms1FeaturesIndex])
                    numMassMatchedFeatures++;
                ms1FeaturesIndex++;
            }
            float massMatchesPercent =
                    (((float) numMassMatchedFeatures) / ((float) runFeatureSet.getFeatures().length) * 100);
//            System.err.println("Matches this run: " + numMassMatchedFeatures + " (" + massMatchesPercent + "%)");
            percentMatchedPerRun[i] = massMatchesPercent;

            runMassMatchPercentMap.put(runEntry, massMatchesPercent);

            if (massMatchesPercent >= minMassMatchPercent)
                numPassing++;
        }

        if (numPassing == 0)
        {
            ApplicationContext.infoMessage("No runs matched enough MS1 features.  Quitting");
            return result;
        }

        AmtRunEntry[] runEntriesSorted = new AmtRunEntry[runEntries.length];
        for (int i=0; i<runEntries.length; i++)
            runEntriesSorted[i] = runEntries[i];

        Comparator<AmtRunEntry> runComparatorByMassMatchPercentDesc =
                new Comparator<AmtRunEntry> ()
        {
            public int compare(AmtRunEntry o1, AmtRunEntry o2)
            {
                float o1Matches = runMassMatchPercentMap.get(o1);
                float o2Matches = runMassMatchPercentMap.get(o2);

                return o1Matches < o2Matches ? 1 : o1Matches > o2Matches ? -1 : 0;
            }
        };

        Arrays.sort(runEntriesSorted, runComparatorByMassMatchPercentDesc);
        ApplicationContext.infoMessage("Runs passing: " + numPassing +
                ".  Adding runs to database...");
        for (AmtRunEntry runEntry : runEntriesSorted)
        {
            float massMatchesPercent = runMassMatchPercentMap.get(runEntry);

            //since runs are sorted, stop if this run doesn't pass
            if (massMatchesPercent < minMassMatchPercent)
                break;

            int numPeptidesIfRunAdded = result.numEntries();
            for (AmtPeptideEntry peptideEntry : amtDatabase.getPeptideEntriesForRun(runEntry))
                if (!result.contains(peptideEntry.getPeptideSequence()))
                    numPeptidesIfRunAdded++;

            if (numPeptidesIfRunAdded > maxNewDBEntries)
            {
                ApplicationContext.setMessage("Too many entries in passing runs, stopping early to keep DB small");
                ApplicationContext.setMessage("Adding " + result.numRuns() + " runs, out of " + numPassing + " mass-matched");
                break;
            }

            Map<MS2Modification, MS2Modification> modificationMap =
                    result.addRunEntry(runEntry);

            for (AmtPeptideEntry peptideEntry : amtDatabase.getPeptideEntriesForRun(runEntry))
            {
                AmtPeptideEntry newEntry = null;
//                    newEntry.setPeptideSequence(peptideEntry.getPeptideSequence());
//                    newEntry.setPredictedHydrophobicity(AmtUtilities.calculateNormalizedHydrophobicity(peptideEntry.getPeptideSequence()));
                for (AmtPeptideEntry.AmtPeptideModificationStateEntry modState : peptideEntry.getModificationStateEntries())
                {
                    AmtPeptideEntry.AmtPeptideObservation obsForRun =
                            modState.getObservationForRun(runEntry);
                    if (obsForRun != null)
                    {
                        if (newEntry == null)
                        {
                            newEntry = AmtPeptideEntry.createEntryFromObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                        }
                        else
                            newEntry.addObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                    }
                }
                if (newEntry != null)
                    result.addObservationsFromEntry(newEntry, modificationMap);


            }
            _log.debug("Adding run " + runEntry.getPepXmlFilename() +
                    ", mass matches: " + massMatchesPercent + "%");

            if (result.numRuns() >= maxNewDBRuns)
            {
                ApplicationContext.setMessage("Too many passing runs, stopping early to keep DB small");
                ApplicationContext.setMessage("Adding " + result.numRuns() + " runs, out of " + numPassing + " mass-matched");
                break;
            }
        }

        if (result.numRuns() == 0)
        {
            ApplicationContext.infoMessage("Too many entries in a single run!  Quitting");
            return result;
        }

        ApplicationContext.infoMessage(result.numRuns() + " out of " + amtDatabase.numRuns() + " runs kept");
        ApplicationContext.infoMessage("Database summary:  " + result);
        return result;
    }


    public static AmtDatabase removeRunsByStructure(AmtDatabase amtDatabase,
                                                       int ms1RunCol, int ms1RunRow, int ms1RunExp,
                                                       int maxNewDBEntries,
                                                       int maxNewDBRuns,
                                                       Fractionation2DUtilities.FractionatedAMTDatabaseStructure amtDatabaseStructure,
                                                       boolean showCharts)
    {
        final Map<AmtRunEntry, Float> runDistanceMap =
                new HashMap<AmtRunEntry, Float>();
        AmtRunEntry[] runEntries = amtDatabase.getRuns();

        for (int i=0; i<runEntries.length; i++)
        {
            Pair<Integer, int[]> expAndPos = amtDatabaseStructure.calculateExperimentAndPosition(i);
            int runEntryCol = expAndPos.second[0];
            int runEntryRow = expAndPos.second[1];

            //manhattan distance
            float distance = Math.abs(runEntryCol - ms1RunCol) + Math.abs(runEntryRow - ms1RunRow);
            runDistanceMap.put(runEntries[i], distance);
        }

        AmtRunEntry[] runEntriesSorted = new AmtRunEntry[runEntries.length];
        for (int i=0; i<runEntriesSorted.length; i++)
            runEntriesSorted[i] = runEntries[i];

        Comparator<AmtRunEntry> runComparatorByDistanceAsc =
                new Comparator<AmtRunEntry> ()
        {
            public int compare(AmtRunEntry o1, AmtRunEntry o2)
            {
                float o1Matches = runDistanceMap.get(o1);
                float o2Matches = runDistanceMap.get(o2);

                return o1Matches < o2Matches ? -1 : o1Matches > o2Matches ? 1 : 0;
            }
        };

        Arrays.sort(runEntriesSorted, runComparatorByDistanceAsc);


        return removeRunsInOrder(amtDatabase, maxNewDBEntries, maxNewDBRuns, amtDatabaseStructure,
                                 runEntriesSorted, showCharts, null, null, null, null);
    }

    protected static Set<String> createPeptideSetFromFeatures(Feature[] features)
    {
        Set<String> peptides = new HashSet<String>();
        for (Feature feature : features)
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide != null)
                peptides.add(peptide);
        }
        return peptides;
    }

    /**
     * Remove runs with insufficient peptide matches to the MS/MS features passed in
     * @param amtDatabase
     * @param ms2Features
     * @param minPeptideMatchPercent
     * @param maxNewDBEntries
     * @param maxNewDBRuns
     * @param showCharts
     * @return
     */
    public static AmtDatabase removeRunsWithoutPeptideMatches(AmtDatabase amtDatabase,
                                                       Feature[] ms2Features,
                                                       int minPeptideMatchPercent,
                                                       int maxNewDBEntries,
                                                       int maxNewDBRuns,
                                                       Fractionation2DUtilities.FractionatedAMTDatabaseStructure amtDatabaseStructure,
                                                       boolean showCharts,
                                                       Feature[] ms1Features,
                                                       MS2Modification[] ms2ModificationsForMatching,
                                                       AmtDatabaseMatcher dbMatcher)
    {
        float[] percentMatchedPerRun = new float[amtDatabase.numRuns()];
        AmtRunEntry[] runEntries = amtDatabase.getRuns();

        final Map<AmtRunEntry, Integer> runPeptideMatchNumberMap =
                new HashMap<AmtRunEntry, Integer>();
        final Map<AmtRunEntry, Float> runPeptideMatchPercentMap =
                new HashMap<AmtRunEntry, Float>();

        Set<String> ms2Peptides = createPeptideSetFromFeatures(ms2Features);

        int numPassing = 0;
        for (int i=0; i<runEntries.length; i++)
        {
            if (i % (Math.max(runEntries.length/10,1)) == 0)
            {
                ApplicationContext.setMessage(Rounder.round(((double) i * 100.0 / (double) runEntries.length),0) +
                        "% done evaluating runs; passing: " + numPassing + " / " + i +
                        " (" + Rounder.round(((double) numPassing  * 100.0) / (double) i, 0) + "%)...");
            }
            AmtRunEntry runEntry = runEntries[i];

            int numPeptidesThisRun = 0;
            int numPeptidesInCommon = 0;
            for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
            {
                String peptideSequence = peptideEntry.getPeptideSequence();
                AmtPeptideEntry.AmtPeptideObservation obs =
                        peptideEntry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    numPeptidesThisRun++;
                    if (ms2Peptides.contains(peptideSequence))
                        numPeptidesInCommon++;
                }
            }

            float peptideMatchesPercent =
                    (((float) numPeptidesInCommon) /
                            ((float) numPeptidesThisRun) * 100);
//            System.err.println("Matches this run: " + numMassMatchedFeatures + " (" + massMatchesPercent + "%)");
            percentMatchedPerRun[i] = peptideMatchesPercent;

            runPeptideMatchPercentMap.put(runEntry, peptideMatchesPercent);
            runPeptideMatchNumberMap.put(runEntry, numPeptidesInCommon);

            if (peptideMatchesPercent >= minPeptideMatchPercent)
                numPassing++;
        }

        if (numPassing == 0)
        {
            ApplicationContext.infoMessage("No runs matched enough MS1 features.  Quitting");
            return new AmtDatabase();
        }

        if (showCharts && amtDatabaseStructure != null)
        {
            double[] peptideMatchesPerRun = new double[amtDatabase.numRuns()];
            for (int i=0; i<peptideMatchesPerRun.length; i++)
                peptideMatchesPerRun[i] = runPeptideMatchNumberMap.get(amtDatabase.getRunBySequence(i+1));

            Fractionation2DUtilities.showHeatMapChart(
                    amtDatabaseStructure,
                    peptideMatchesPerRun, "AMT Database Run Peptides In Common with This Run",
                    true);

            double[] peptidePercentMatchesPerRun = new double[amtDatabase.numRuns()];
            for (int i=0; i<peptidePercentMatchesPerRun.length; i++)
                peptidePercentMatchesPerRun[i] = runPeptideMatchPercentMap.get(amtDatabase.getRunBySequence(i+1));

            Fractionation2DUtilities.showHeatMapChart(
                    amtDatabaseStructure,
                    peptidePercentMatchesPerRun, "AMT Database Run Peptide % In Common with This Run",
                    true);
        }

        List<AmtRunEntry> passingRunEntriesList = new ArrayList<AmtRunEntry>(runEntries.length);
        for (AmtRunEntry runEntry : runEntries)
            if (runPeptideMatchPercentMap.get(runEntry) >= minPeptideMatchPercent)
                passingRunEntriesList.add(runEntry);

        AmtRunEntry[] runEntriesSorted = new AmtRunEntry[passingRunEntriesList.size()];
        for (int i=0; i<runEntriesSorted.length; i++)
            runEntriesSorted[i] = passingRunEntriesList.get(i);

        Comparator<AmtRunEntry> runComparatorByPeptideMatchPercentDesc =
                new Comparator<AmtRunEntry> ()
        {
            public int compare(AmtRunEntry o1, AmtRunEntry o2)
            {
                float o1Matches = runPeptideMatchPercentMap.get(o1);
                float o2Matches = runPeptideMatchPercentMap.get(o2);

                return o1Matches < o2Matches ? 1 : o1Matches > o2Matches ? -1 : 0;
            }
        };

        Comparator<AmtRunEntry> runComparatorByPeptideMatchNumberDesc =
                new Comparator<AmtRunEntry> ()
        {
            public int compare(AmtRunEntry o1, AmtRunEntry o2)
            {
                float o1Matches = runPeptideMatchNumberMap.get(o1);
                float o2Matches = runPeptideMatchNumberMap.get(o2);

                return o1Matches < o2Matches ? 1 : o1Matches > o2Matches ? -1 : 0;
            }
        };

        Arrays.sort(runEntriesSorted, runComparatorByPeptideMatchPercentDesc);

        ApplicationContext.infoMessage("Runs passing: " + numPassing +
                ".  Adding runs to database...");
        return removeRunsInOrder(amtDatabase, maxNewDBEntries, maxNewDBRuns, amtDatabaseStructure,
                                 runEntriesSorted, showCharts, ms1Features, ms2ModificationsForMatching,
                                 ms2Features, dbMatcher);

    }

    /**
     * Trims down an AMT database by adding runs one by one from runEntriesSorted, stopping when
     * adding another run would violate maxNewDBEntries or maxNewDBRuns. 
     * @param amtDatabase
     * @param maxNewDBEntries
     * @param maxNewDBRuns
     * @param amtDatabaseStructure
     * @param runEntriesSorted  sorted in _descending_ order of goodness
     * @param showCharts
     * @return
     */
    public static AmtDatabase removeRunsInOrder(AmtDatabase amtDatabase,
                                                int maxNewDBEntries,
                                                int maxNewDBRuns,
                                                Fractionation2DUtilities.FractionatedAMTDatabaseStructure amtDatabaseStructure,
                                                AmtRunEntry[] runEntriesSorted,
                                                boolean showCharts,
                                                Feature[] ms1Features,
                                                MS2Modification[] ms2ModificationsForMatching,
                                                Feature[] ms2Features,
                                                AmtDatabaseMatcher databaseMatcher)
    {
        //if true, we will engage in an extremely costly plot of FAR vs. number of runs kept in the database.
        //This was done in order to produce that plot for the HUPO conference.  Not normally useful, as
        //it's very hand-wavy
        boolean plotFarVsNumRuns = false;

        AmtDatabase result = new AmtDatabase();

        int numSortedRuns = runEntriesSorted.length;

        //variables for db density vs. FAR chart
        List<Double> farsWithNumbersOfRuns = new ArrayList<Double>();
        List<Double> numEntriesWithNumbersOfRuns = new ArrayList<Double>();
        

        List<AmtRunEntry> runsAdded = new ArrayList<AmtRunEntry>();
        for (AmtRunEntry runEntry : runEntriesSorted)
        {
            int numPeptidesIfRunAdded = result.numEntries();
            for (AmtPeptideEntry peptideEntry : amtDatabase.getPeptideEntriesForRun(runEntry))
                if (!result.contains(peptideEntry.getPeptideSequence()))
                    numPeptidesIfRunAdded++;

            if (numPeptidesIfRunAdded > maxNewDBEntries)
            {
                ApplicationContext.setMessage("Too many entries in passing runs, stopping early to keep DB small");
                ApplicationContext.setMessage("Adding " + result.numRuns() + " runs, out of " + numSortedRuns);
                break;
            }
            if (result.numRuns() + 1 >= maxNewDBRuns)
            {
                ApplicationContext.setMessage("Too many passing runs, stopping early to keep DB small");
                ApplicationContext.setMessage("Adding " + result.numRuns() + " runs, out of " + numSortedRuns);
                break;
            }

            Map<MS2Modification, MS2Modification> modificationMap =
                    result.addRunEntry(runEntry);

            for (AmtPeptideEntry peptideEntry : amtDatabase.getPeptideEntriesForRun(runEntry))
            {
                AmtPeptideEntry newEntry = null;
                for (AmtPeptideEntry.AmtPeptideModificationStateEntry modState :
                        peptideEntry.getModificationStateEntries())
                {
                    AmtPeptideEntry.AmtPeptideObservation obsForRun =
                            modState.getObservationForRun(runEntry);
                    if (obsForRun != null)
                    {
                        if (newEntry == null)
                        {
                            newEntry = AmtPeptideEntry.createEntryFromObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                        }
                        else
                            newEntry.addObservation(peptideEntry.getPeptideSequence(),
                                    modState.getModifications(), obsForRun);
                    }
                }
                if (newEntry != null)
                    result.addObservationsFromEntry(newEntry, modificationMap);


            }
            runsAdded.add(runEntry);
            _log.debug("Adding run " + runEntry.getPepXmlFilename() + " (#" + result.numRuns() + ")");
        }

        if (result.numRuns() == 0)
        {
            ApplicationContext.infoMessage("Too many entries in a single run!  Quitting");
            return result;
        }

        if (showCharts && amtDatabaseStructure != null)
        {
            double[] matchedRuns = new double[amtDatabase.numRuns()];
            for (AmtRunEntry runEntry : amtDatabase.getRuns())
            {
                if (runsAdded.contains(runEntry))
                {
                    matchedRuns[amtDatabase.getSequenceForRun(runEntry)-1]++;
                }
            }

            Fractionation2DUtilities.showHeatMapChart(
                    amtDatabaseStructure,
                    matchedRuns, "AMT Runs Kept For Matching",
                    false);
        }

        //Only for creating an extremely costly plot        
        if (showCharts && ms1Features != null && ms2ModificationsForMatching != null && plotFarVsNumRuns)
        {
            List<Double> numbersOfRuns = new ArrayList<Double>();
            for (int i=0; i<farsWithNumbersOfRuns.size(); i++)
                numbersOfRuns.add((double)(i+1));
            System.err.println("Line chart.  Size: " + farsWithNumbersOfRuns);
            PanelWithLineChart pwlc = new PanelWithLineChart(numbersOfRuns, farsWithNumbersOfRuns,
                    "FAR vs. Number of Fractions");
            ChartDialog cd = new ChartDialog(pwlc);
            cd.setTitle("Line chart: FAR vs. number of fractions");
            cd.setVisible(true);

            PanelWithLineChart pwlc2 = new PanelWithLineChart(numEntriesWithNumbersOfRuns, farsWithNumbersOfRuns,
                    "FAR vs. Number of Entries");
            ChartDialog cd2 = new ChartDialog(pwlc2);
            cd2.setTitle("Line chart: FAR vs. number of entries");
            cd2.setVisible(true);
        }

        ApplicationContext.infoMessage(result.numRuns() + " out of " + amtDatabase.numRuns() + " runs kept");
        ApplicationContext.infoMessage("Database summary:  " + result);
        return result;
    }



}
