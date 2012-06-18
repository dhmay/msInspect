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
import java.util.List;
import java.io.IOException;
import java.io.File;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.BaseFeatureSetMatcherImpl;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.Window2DFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.jfree.chart.JFreeChart;

import javax.swing.*;


/**
 * This class performs simple matching between an AMT database and
 * a feature files.  It constructs a probability model for the matches
 * based on mass and NRT deviation.
 */
public class AmtDatabaseMatcher
{
    static Logger _log = Logger.getLogger(AmtDatabaseMatcher.class);

    //These defaults are for robust regression
    //default value for the cutoff numerator in the leverage cutoff
    public static final double DEFAULT_LEVERAGE_NUMERATOR = 6;
    public static final double DEFAULT_MAX_LEVERAGE_NUMERATOR = 12;
    public static final double DEFAULT_MAX_STUDENTIZED_RESIDUAL = 3.0;

    public static final float DEFAULT_MASS_MATCH_DELTA_MASS = 5;
    public static final int DEFAULT_MASS_MATCH_DELTA_MASS_TYPE =
            FeatureSetMatcher.DELTA_MASS_TYPE_PPM;
    public static final float DEFAULT_2D_MATCH_DELTA_MASS = 10;
    public static final int DEFAULT_2D_MATCH_DELTA_MASS_TYPE =
            FeatureSetMatcher.DELTA_MASS_TYPE_PPM;
    public static final float DEFAULT_2D_MATCH_DELTA_ELUTION = 0.15f;

    //minimum number of matches required to perform a regression.  This is probably not conservative enough
    public static final int MIN_MATCHED_FEATURES_FOR_REGRESSION = 8;

    //minimum number of feature pairs that must be passed to Yan's
    //quantile regression code to get it to work.
    public static final int DEFAULT_QUANTILE_REG_MIN_FEATURES = 85;
    protected int quantileRegressionMinFeatures = DEFAULT_QUANTILE_REG_MIN_FEATURES;

    //default degree of the polynomial to fit when mapping T to H nonlinearly
    public static final int DEFAULT_NONLINEAR_MAPPING_DEGREE = 5;

    protected int nonlinearMappingPolynomialDegree = DEFAULT_NONLINEAR_MAPPING_DEGREE;

    //what's a 'significant' difference in hydrophobicity, in terms of figuring
    //out whether a single-observation peptide should be used for matching?        feat
    public static final float DEFAULT_SIGNIFICANT_HYDRO_DIFFERENCE = .4f;

    //These variables are for the mass-only match that we use in order to come up
    //with (ms1,amt) pairs for alignment if we have no embedded MS2
    protected float massMatchDeltaMass = DEFAULT_MASS_MATCH_DELTA_MASS;
    protected int massMatchDeltaMassType = DEFAULT_MASS_MATCH_DELTA_MASS_TYPE;

    //These parameters are used in the actual ('loose') match
    protected float realMatchDeltaMass = DEFAULT_2D_MATCH_DELTA_MASS;
    protected int realMatchDeltaMassType = DEFAULT_2D_MATCH_DELTA_MASS_TYPE;
    protected float realMatchDeltaElution = DEFAULT_2D_MATCH_DELTA_ELUTION;

    //Controls the number of iterations for EM model
    protected int minEMIterations = AmtMatchProbabilityAssigner.DEFAULT_MIN_EM_ITERATIONS;
    protected int maxEMIterations = AmtMatchProbabilityAssigner.DEFAULT_MAX_EM_ITERATIONS;

    protected int maxRProbAssignmentMillis = AmtMatchProbabilityAssigner.DEFAULT_MAX_EM_ITERATIONS;


    //Should we use MS1 times for alignment?  This requires mass-and-time matching between MS1 and MS2
    public static final boolean DEFAULT_USE_MS1_TIMES_FOR_ALIGNMENT = true;
    protected boolean useMs1TimesForAlignment = DEFAULT_USE_MS1_TIMES_FOR_ALIGNMENT;

    //tolerances for matching MS1 features to MS2 for alignment
    protected float ms1Ms2MassTolerancePPM = AmtDatabaseBuilder.DEFAULT_MS1_MS2_MASS_TOLERANCE_PPM;
    protected float ms1Ms2TimeToleranceSeconds = AmtDatabaseBuilder.DEFAULT_MS1_MS2_TIME_TOLERANCE_SECONDS;
    

    //if this is true, we do a dummy match instead of a real match
    protected boolean doDecoyMatch = false;

    //Parameters for robust regression
    protected double maxRegressionLeverageNumerator =
            DEFAULT_MAX_LEVERAGE_NUMERATOR;
    protected double maxRegressionStudRes =
            DEFAULT_MAX_STUDENTIZED_RESIDUAL;

    protected float minMatchProbabilityToKeep = AmtMatchProbabilityAssigner.DEFAULT_MIN_MATCH_PROBABILITY;
    protected float maxMatchFDRToKeep = AmtMatchProbabilityAssigner.DEFAULT_MAX_MATCH_FDR;


    //Hard maximum on second-best probability
    protected float maxSecondBestProbability = AmtMatchProbabilityAssigner.DEFAULT_MAX_SECONDBEST_PROBABILITY;
    //Minimum difference between second-best probability and best probability
    protected float minSecondBestProbabilityDifference =
            AmtMatchProbabilityAssigner.DEFAULT_MIN_SECONDBEST_PROBABILITY_DIFFERENCE;
    

    //For creation of decoy database, used in determining probabilities
    public static final int DEFAULT_DECOY_DB_MASS_ADJUSTMENT_DA = 11;

    //Persisting the calculated mapping coefficients so they can be accessed elsewhere
    public double[] timeHydMapCoefficients = null;

    //control whether to display useful charts when matching to databases
    protected boolean buildCharts = false;

    //Stores the result of matching
    public FeatureSetMatcher.FeatureMatchingResult featureMatchingResult = null;

    //Persist charts related to matching, for later retrieval
    protected JFreeChart massMatchScatterplot = null;
    protected JFreeChart timeHydrophobicityMappingChart = null;
    protected JFreeChart massCalibrationChart = null;
    protected PanelWithRPerspectivePlot massTimeErrorPerspectivePlot = null;
    public JFreeChart massDeltaMassScatterPlot = null;

    //amount by which to adjust AMT feature masses.  For dummy matching
    protected double amtFeatureMassAdjustment = 0;

    //Stores information about the AMT database's structure.  For matching
    //fractionated data to fractionated data
    protected boolean amtDBDimensionsDefined = false;
    protected Fractionation2DUtilities.FractionatedAMTDatabaseStructure
            amtDatabaseStructure;

    protected AmtMatchProbabilityAssigner probabilityAssigner = null;


    public void matchWithFDRCalc(AmtDatabase amtDatabaseThisMatch,
                                Feature[] amtDBFeaturesThisMatch,
                                FeatureSet ms1FeatureSetToMatch,
                                FeatureSet embeddedMs2FeatureSet,
                                float minEmbeddedMs2PeptideProphet,
                                MS2Modification[] ms2ModificationsArray,
                                File matchingOutputFile,
                                boolean removeFractions,
                                int minRunsToKeep, int maxRunsToKeep,
                                boolean calibrateMassesUsingMatches,
                                boolean showCharts)
            throws IOException
    {
        float portionDecoy = .5f;

        int numTotalDBEntries = amtDatabaseThisMatch.numEntries();
        int numDecoy = Math.round(portionDecoy * numTotalDBEntries);

        AmtPeptideEntry[] allEntries = amtDatabaseThisMatch.getEntries();

        Set<String> decoyPeptides1 = new HashSet<String>();

        //select decoy peptides
        while (decoyPeptides1.size() < numDecoy)
        {
            int entryIndex = (int) Math.round(Math.random() * (numTotalDBEntries-1));
            AmtPeptideEntry entry =  allEntries[entryIndex];
            String peptide = entry.getPeptideSequence();
            if (!decoyPeptides1.contains(peptide))
            {
                decoyPeptides1.add(peptide);
            }
        }

        Set<String> decoyPeptides2 = new HashSet<String>();
        for (AmtPeptideEntry peptideEntry : allEntries)
        {
            if (!decoyPeptides1.contains(peptideEntry.getPeptideSequence()))
                decoyPeptides2.add(peptideEntry.getPeptideSequence());
        }



        Map<Feature, Float> firstHalfMatchesWithFDRs = matchWithPortionDecoy(
                amtDatabaseThisMatch,
                amtDBFeaturesThisMatch,
                ms1FeatureSetToMatch,
                embeddedMs2FeatureSet,
                minEmbeddedMs2PeptideProphet,
                ms2ModificationsArray,
                matchingOutputFile,
                removeFractions,
                minRunsToKeep, maxRunsToKeep,
                calibrateMassesUsingMatches,
                showCharts,
                decoyPeptides1, 1f);


    }

    /**
     * Change the masses of a random portion (portionDecoy) of the AMT
     * database entries to make them decoy entries, perform the match,
     * and report on the specificity of the match
     * @param amtDatabaseThisMatch
     * @param amtDBFeaturesThisMatch
     * @param ms1FeatureSetToMatch
     * @param embeddedMs2FeatureSet
     * @param minEmbeddedMs2PeptideProphet
     * @param ms2ModificationsArray
     * @param matchingOutputFile
     * @param removeFractions
     * @param minRunsToKeep
     * @param maxRunsToKeep
     * @param calibrateMassesUsingMatches
     * @param showCharts
     * @param decoyPeptides
     */
    public Map<Feature, Float> matchWithPortionDecoy(AmtDatabase amtDatabaseThisMatch,
                                Feature[] amtDBFeaturesThisMatch,
                                FeatureSet ms1FeatureSetToMatch,
                                FeatureSet embeddedMs2FeatureSet,
                                float minEmbeddedMs2PeptideProphet,
                                MS2Modification[] ms2ModificationsArray,
                                File matchingOutputFile,
                                boolean removeFractions,
                                int minRunsToKeep, int maxRunsToKeep,
                                boolean calibrateMassesUsingMatches,
                                boolean showCharts,
                                Set<String> decoyPeptides,
                                float targetDecoyRatio)
            throws IOException
    {
        List<Feature> targetDecoyDBFeatures = new ArrayList<Feature>();
        for (Feature feature : amtDBFeaturesThisMatch)
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

            if (decoyPeptides.contains(peptide))
            {
                Feature decoyFeature = new Feature(feature);
                decoyFeature.setMass(feature.getMass() + DEFAULT_DECOY_DB_MASS_ADJUSTMENT_DA);
                targetDecoyDBFeatures.add(decoyFeature);
            }
            else
                targetDecoyDBFeatures.add(feature);
        }

        List<Feature> matchedMS1Features = matchAgainstMs1(amtDatabaseThisMatch,
                targetDecoyDBFeatures.toArray(new Feature[targetDecoyDBFeatures.size()]),
                ms1FeatureSetToMatch,
                embeddedMs2FeatureSet,
                minEmbeddedMs2PeptideProphet,
                ms2ModificationsArray,
                matchingOutputFile,
                removeFractions,
                minRunsToKeep, maxRunsToKeep,
                calibrateMassesUsingMatches,
                showCharts);

        Collections.sort(matchedMS1Features, new PeptideProphetComparatorDesc());
        int numTargetMatches = 0;
        int numDecoyMatches = 0;
        float[] fdrs = new float[matchedMS1Features.size()];
        float[] probabilities = new float[matchedMS1Features.size()];
        float[] sums1MinusProbs = new float[matchedMS1Features.size()];
        float[] sumsProbs = new float[matchedMS1Features.size()];
        float[] theoreticalFDRs = new float[matchedMS1Features.size()];





        Map<Feature, Float> result = new HashMap<Feature, Float>();
        float sum1MinusProbs = 0f;
        float sumProbs = 0f;

        for (int i=0; i<matchedMS1Features.size(); i++)
        {
            Feature matchedMS1Feature = matchedMS1Features.get(i);
            boolean decoyMatch = decoyPeptides.contains(MS2ExtraInfoDef.getFirstPeptide(matchedMS1Feature));
            if (decoyMatch)
                numDecoyMatches++;
            else
                numTargetMatches++;
            probabilities[i] = (float) MS2ExtraInfoDef.getPeptideProphet(matchedMS1Feature);
            float fdr = 1f;
            if (numTargetMatches > 0)
                fdr = targetDecoyRatio * numDecoyMatches / numTargetMatches;
            fdrs[i] = fdr;

            sum1MinusProbs += 1 - probabilities[i];
            sums1MinusProbs[i] = sum1MinusProbs;
            sumProbs += probabilities[i];
            sumsProbs[i] = sumProbs;

            theoreticalFDRs[i] = (sum1MinusProbs / (i+1));

            if (!decoyMatch)
                result.put(matchedMS1Feature, fdr);    
        }

        if (showCharts)
        {
            PanelWithLineChart pwlc = new PanelWithLineChart(probabilities, fdrs,
                    "Prob vs FDR");
            pwlc.setAxisLabels("Match Probability", "FDR (#decoy/#target)");            
            pwlc.displayInTab();

            PanelWithLineChart pwlc2 = new PanelWithLineChart(theoreticalFDRs, fdrs,
                    "Theoretical FDR vs FDR");
            pwlc2.setAxisLabels("Theoretical FDR: sum(1-p)/sum(1)", "Actual FDR (#decoy/#target)");
            pwlc2.displayInTab();
        }

        return result;
    }



    public static class PeptideProphetComparatorDesc implements Comparator<Feature>
    {
        public int compare(Feature o1, Feature o2)
        {
            double o1Value = MS2ExtraInfoDef.getPeptideProphet(o1);
            double o2Value = MS2ExtraInfoDef.getPeptideProphet(o2);

            if (o1Value == o2Value)
                return 0;
            return o1Value < o2Value ? 1 : -1;

        }

    }

    /**
     * For each peptide, per modification state, remove all but one observation.  That observation
     * gets the median retention time of all of them
     * @param featureSet
     * @return
     */
    public static void representPeptidesWithMedianTimePerPeptidePerMod(FeatureSet featureSet)
    {

        //map peptides to features
        HashMap<String,List<Feature>> peptideFeatureListMap =
                new HashMap<String,List<Feature>>();

        Feature[] features = featureSet.getFeatures();

        int initialNumFeatures = features.length;

        for (Feature feature : features)
        {
            String peptide =
                    MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            List<Feature> thisPeptideFeatureList = peptideFeatureListMap.get(peptide);
            if (thisPeptideFeatureList == null)
            {
                thisPeptideFeatureList = new ArrayList<Feature>();
                peptideFeatureListMap.put(peptide, thisPeptideFeatureList);
            }
            thisPeptideFeatureList.add(feature);
        }

        List<Feature> resultFeatureList = new ArrayList<Feature>();

        //for debugging
        int sumNumModStatesPerPeptide = 0;

        for (String peptide : peptideFeatureListMap.keySet())
        {
            List<Feature> featuresThisPeptide = peptideFeatureListMap.get(peptide);

            //only exact mod state matches fold together
            Map<String, List<Feature>> modStateFeatureListMap =
                    new HashMap<String, List<Feature>>();
            for (Feature featureThisPeptide : featuresThisPeptide)
            {
                MS2ExtraInfoDef.getModifiedAminoAcids(featureThisPeptide);
                Map<Integer, List<ModifiedAminoAcid>> positionModListMap =
                        (Map<Integer, List<ModifiedAminoAcid>>)
                                featureThisPeptide.getProperty("modifiedaminoacids");
                //special handling for no modifications
                String modStateString = "";
                if (positionModListMap != null && positionModListMap.size() > 0)
                    modStateString = MS2ExtraInfoDef.getSingletonInstance().convertToString("modifiedaminoacids", positionModListMap);
                List<Feature> featuresThisModState = modStateFeatureListMap.get(modStateString);
                if (featuresThisModState == null)
                {
                    featuresThisModState = new ArrayList<Feature>();
                    modStateFeatureListMap.put(modStateString, featuresThisModState);
                }
                featuresThisModState.add(featureThisPeptide);
            }

            if (_log.isDebugEnabled())
                sumNumModStatesPerPeptide += modStateFeatureListMap.size();
            //add the first for each mod state
            for (String modStateString : modStateFeatureListMap.keySet())
            {
                List<Feature> featuresThisModState =
                        modStateFeatureListMap.get(modStateString);
                List<Float> retentionTimesThisModState = new ArrayList<Float>();
                for (Feature feature : featuresThisModState)
                     retentionTimesThisModState.add(feature.getTime());

                float medianTime = (float) BasicStatistics.median(retentionTimesThisModState);

                int lastScan = 0;
                Feature firstFeature = null;

                for (Feature feature : featuresThisModState)
                {
                    if (feature.getScan() > lastScan)
                        lastScan = feature.getScan();
                    if (firstFeature == null || feature.getScan() < firstFeature.getScan())
                        firstFeature = feature;
                }
                firstFeature.setScanLast(lastScan);
                firstFeature.setTime(medianTime);
                resultFeatureList.add(firstFeature);
            }
        }

        //sort all the earliest features for each peptide & set them as the feature array
        //for this featureset
        Collections.sort(resultFeatureList, new Feature.MzScanAscComparator());
        _log.debug("representPeptidesWithMedianTimePerPeptidePerMod, initial features: " + initialNumFeatures + 
                ", resulting features: " + resultFeatureList.size());

        featureSet.setFeatures(resultFeatureList.toArray(new Feature[0]));

    }


    /**
     * Match an AMT database against a single MS1 feature file. This is the master matching method
     *
     * @return a list of matched MS1 features
     */
    public List<Feature> matchAgainstMs1(AmtDatabase amtDatabaseThisMatch,
                                Feature[] amtDBFeaturesThisMatch,
                                FeatureSet ms1FeatureSetToMatch,
                                FeatureSet embeddedMs2FeatureSet,
                                float minEmbeddedMs2PeptideProphet,
                                MS2Modification[] ms2ModificationsArray,
                                File matchingOutputFile,
                                boolean removeFractions,
                                int minRunsToKeep, int maxRunsToKeep,
                                boolean calibrateMassesUsingMatches,
                                boolean showCharts)
            throws IOException
    {
        boolean recordDBRunsInFeatureFile = false;

        ApplicationContext.setMessage("Matching against file " + ms1FeatureSetToMatch.getSourceFile().getName());

        Feature[] guideFeaturesForAlignment = null;
        if (embeddedMs2FeatureSet != null)
        {
            FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
            sel.setMinPProphet(minEmbeddedMs2PeptideProphet);
            _log.debug("Embedded MS2: before filter, " + embeddedMs2FeatureSet.getFeatures().length + " features");
            embeddedMs2FeatureSet = embeddedMs2FeatureSet.filter(sel);
            _log.debug("Embedded MS2: after filter, " + embeddedMs2FeatureSet.getFeatures().length + " features");

            if (embeddedMs2FeatureSet.getFeatures().length < MIN_MATCHED_FEATURES_FOR_REGRESSION)
                throw new IllegalArgumentException("ERROR: after filter, too few MS2 features (" +
                    embeddedMs2FeatureSet.getFeatures().length + ") for alignment.  You can try again with a " +
                    "less restrictive cutoff.");

            boolean nonzeroMs2TimesExist = false;
            for (Feature feature : embeddedMs2FeatureSet.getFeatures())
            {
                if (feature.getTime() > 0)
                {
                    nonzeroMs2TimesExist = true;
                    break;
                }
            }
            if (!nonzeroMs2TimesExist)
            {
                throw new IllegalArgumentException("ERROR! MS2 feature retention times are all zero!  " +
                    "Please populate feature " +
                    "times, using the 'populatems2times' command, or provide associated mzXML file(s) to " +
                    "use for populating scan times.");
            }

            if (useMs1TimesForAlignment)
            {                
                representPeptidesWithMedianTimePerPeptidePerMod(embeddedMs2FeatureSet);

                Window2DFeatureSetMatcher featureSetMatcher =
                        new Window2DFeatureSetMatcher();
                featureSetMatcher.setMassDiffType(FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
                featureSetMatcher.setMaxMassDiff(ms1Ms2MassTolerancePPM);
                featureSetMatcher.setMinMassDiff(-ms1Ms2MassTolerancePPM);
                featureSetMatcher.setMaxElutionDiff(ms1Ms2TimeToleranceSeconds);
                featureSetMatcher.setMinElutionDiff(-ms1Ms2TimeToleranceSeconds);
                featureSetMatcher.setElutionMode(BaseFeatureSetMatcherImpl.ELUTION_MODE_TIME);

                FeatureSetMatcher.FeatureMatchingResult ms1MS2MatchingResult =
                        featureSetMatcher.matchFeatures(ms1FeatureSetToMatch, embeddedMs2FeatureSet);
                List<Feature> singlyMatchedMS1FeatureCopies = new ArrayList<Feature>();
                for (Feature feature : ms1MS2MatchingResult.getMasterSetFeatures())
                {
                    List<Feature> ms2MatchedFeatures = ms1MS2MatchingResult.getSlaveSetFeatures(feature);
                    Set<String> peptides = new HashSet<String>();
                    for (Feature ms2MatchedFeature : ms2MatchedFeatures)
                        peptides.add(MS2ExtraInfoDef.getFirstPeptide(ms2MatchedFeature));
                    if (peptides.size() == 1)
                    {
                        Feature featureCopy = new Feature(feature);
                        MS2ExtraInfoDef.setSinglePeptide(featureCopy, peptides.iterator().next());
                        singlyMatchedMS1FeatureCopies.add(featureCopy);
                    }
                }
                guideFeaturesForAlignment = singlyMatchedMS1FeatureCopies.toArray(new Feature[singlyMatchedMS1FeatureCopies.size()]);
                ApplicationContext.infoMessage("Using MS1 features for alignment, matching to MS2.  " +
                        guideFeaturesForAlignment.length + " out of " +
                        ms1FeatureSetToMatch.getFeatures().length + " MS1 features used for alignment");
            }
            else
            {
                _log.debug("Using MS2 features for alignment");
                MS2ExtraInfoDef.removeAllButFirstFeatureForEachPeptide(embeddedMs2FeatureSet);
                guideFeaturesForAlignment = embeddedMs2FeatureSet.getFeatures();
            }

            _log.debug("Guide features for alignment: " + guideFeaturesForAlignment.length);

        }


        FeatureSet amtDatabaseFeatureSet =
                new FeatureSet(amtDBFeaturesThisMatch);

        calculateFeatureHydrophobicities(
                ms1FeatureSetToMatch.getFeatures(),
                guideFeaturesForAlignment,
                amtDatabaseFeatureSet, nonlinearMappingPolynomialDegree);



        if (showCharts)
        {
            JFreeChart nonlinearChart = getTimeHydrophobicityMappingChart();
            PanelWithChart pwc = new PanelWithChart(nonlinearChart);
            pwc.setName("T->H Map");
            pwc.displayInTab();
        }

        //mess with the masses of all AMT database entries
        if (doDecoyMatch)
        {
            FeatureSet dummyDBFeatureSet = new FeatureSet();
            List<Feature> dummyDBFeatures = new ArrayList<Feature>();
            for (Feature oldDBFeature : amtDatabaseFeatureSet.getFeatures())
            {
                Feature newFeature = new Feature(oldDBFeature);
                newFeature.setMass(newFeature.getMass() + DEFAULT_DECOY_DB_MASS_ADJUSTMENT_DA);
                dummyDBFeatures.add(newFeature);
            }
            dummyDBFeatureSet.setFeatures(dummyDBFeatures.toArray(new Feature[0]));
            amtDatabaseFeatureSet = dummyDBFeatureSet;
        }

        ApplicationContext.infoMessage("Loose matching with tolerances " +
                realMatchDeltaMass + ", " + realMatchDeltaElution);

        FeatureSetMatcher.FeatureMatchingResult looseMatchingResult =
                callWindowMatcher(ms1FeatureSetToMatch, amtDatabaseFeatureSet,
                        -realMatchDeltaMass, realMatchDeltaMass,
                        -realMatchDeltaElution, realMatchDeltaElution);
//annotateFractionConcordanceForMatches(amtDatabaseThisMatch, embeddedMs2Features, looseMatchingResult);
        if (removeFractions)
        {
            amtDatabaseThisMatch =
                    reduceDatabaseByRunSimilarity(amtDatabaseThisMatch,
                            guideFeaturesForAlignment, looseMatchingResult,
                            minRunsToKeep, maxRunsToKeep, showCharts);
            AmtDatabaseFeatureSetGenerator featureGenerator =
                new AmtDatabaseFeatureSetGenerator(amtDatabaseThisMatch);
            amtDBFeaturesThisMatch =
                    featureGenerator.createFeaturesForModifications(ms2ModificationsArray);
            amtDatabaseFeatureSet =
                new FeatureSet(amtDBFeaturesThisMatch);

            calculateFeatureHydrophobicities(
                    ms1FeatureSetToMatch.getFeatures(),
                    guideFeaturesForAlignment,
                    amtDatabaseFeatureSet, nonlinearMappingPolynomialDegree);
            if (showCharts)
            {
                JFreeChart nonlinearChart = getTimeHydrophobicityMappingChart();
                PanelWithChart pwc = new PanelWithChart(nonlinearChart);
                pwc.setName("T->H Map after fraction removal");
                pwc.displayInTab();
            }
            looseMatchingResult =
                callWindowMatcher(ms1FeatureSetToMatch, amtDatabaseFeatureSet,
                        -realMatchDeltaMass, realMatchDeltaMass,
                        -realMatchDeltaElution, realMatchDeltaElution);
        }

        if (calibrateMassesUsingMatches)
        {
            ApplicationContext.setMessage("Loose match before calibrate: " + looseMatchingResult.size() +
                    " matches. Calibrating...");

            calibrateMS1FeaturesWithMatches(
                    ms1FeatureSetToMatch.getFeatures(),
                    looseMatchingResult,
                    showCharts);

            looseMatchingResult =
                    callWindowMatcher(ms1FeatureSetToMatch, amtDatabaseFeatureSet,
                            -realMatchDeltaMass, realMatchDeltaMass,
                            -realMatchDeltaElution, realMatchDeltaElution);

        }

        if (showCharts)
        {
            List<Float> massErrors = new ArrayList<Float>();
            List<Float> hErrors = new ArrayList<Float>();

            for (Feature masterSetFeature : looseMatchingResult.getMasterSetFeatures())
            {
                List<Feature> matchesThisFeature = looseMatchingResult.get(masterSetFeature);
                for (Feature matchedFeature : matchesThisFeature)
                {
                    double deltaMass =  (masterSetFeature.getMass() - matchedFeature.getMass()) *
                                         (1000000 / masterSetFeature.getMass());
                    double deltaH = (AmtExtraInfoDef.getObservedHydrophobicity(masterSetFeature) -
                            (AmtExtraInfoDef.getObservedHydrophobicity(matchedFeature)));
                    massErrors.add((float)deltaMass);
                    hErrors.add((float)deltaH);
                }
            }
            PanelWithScatterPlot pspError =
                    new PanelWithScatterPlot(hErrors, massErrors, "Loose match error data");
            pspError.setAxisLabels("deltaNRT (NRT units)", "deltaMass (ppm)");
            pspError.setPointSize(2);
            pspError.displayInTab();
        }

        Feature[] amtDecoyFeatures = new Feature[amtDBFeaturesThisMatch.length];
        for (int i=0; i<amtDecoyFeatures.length; i++)
        {
            amtDecoyFeatures[i] = new Feature(amtDBFeaturesThisMatch[i]);
            amtDecoyFeatures[i].setMass(amtDecoyFeatures[i].getMass() +
                                        AmtDatabaseMatcher.DEFAULT_DECOY_DB_MASS_ADJUSTMENT_DA);
        }
        FeatureSet amtDecoyFeatureSet = new FeatureSet(amtDecoyFeatures);
        FeatureSetMatcher.FeatureMatchingResult decoyMatchingResult =
                callWindowMatcher(ms1FeatureSetToMatch, amtDecoyFeatureSet,
                        -realMatchDeltaMass, realMatchDeltaMass,
                        -realMatchDeltaElution, realMatchDeltaElution);


        ApplicationContext.setMessage("Loose match: " + looseMatchingResult.size() + " matches.  " +
                "Decoy: " + decoyMatchingResult.size() + " matches");




        probabilityAssigner =
                new AmtMatchProbabilityAssigner(
                        -realMatchDeltaMass, realMatchDeltaMass,
                        -realMatchDeltaElution, realMatchDeltaElution,
                        minMatchProbabilityToKeep, maxMatchFDRToKeep);
        probabilityAssigner.setMaxSecondBestProbability(maxSecondBestProbability);
        probabilityAssigner.setMinSecondBestProbabilityDifference(minSecondBestProbabilityDifference);
        probabilityAssigner.setMinEMIterations(minEMIterations);
        probabilityAssigner.setMaxEMIterations(maxEMIterations);
        probabilityAssigner.setMaxRProbAssignmentMillis(maxRProbAssignmentMillis);
        

        List<Feature> matchedMS1Features = probabilityAssigner.assignMatchesAndProbabilities(
                looseMatchingResult, decoyMatchingResult,
                showCharts);
        ms1FeatureSetToMatch.addExtraInformationType(AmtExtraInfoDef.getSingletonInstance());


        if (matchingOutputFile != null)
            _log.debug("\toutput file " + matchingOutputFile.getName());
        _log.debug("temp dir = " + TempFileManager.getTmpDir().getAbsolutePath());


        Set<String> thisRunMatchedPeptides = new HashSet<String>();
        for (Feature feature : ms1FeatureSetToMatch.getFeatures())
        {
            List<String> peptidesThisFeature = MS2ExtraInfoDef.getPeptideList(feature);
            if (peptidesThisFeature == null)
                continue;
            for (String peptide : peptidesThisFeature)
            {
                thisRunMatchedPeptides.add(peptide);
            }
        }

        ApplicationContext.infoMessage("    (matched " + thisRunMatchedPeptides.size() + " distinct peptides)");
        ApplicationContext.infoMessage((100 * (double) thisRunMatchedPeptides.size() / (double) amtDatabaseThisMatch.numEntries()) +
                " % of database");

//        allMatchedPeptides.addAll(thisRunMatchedPeptides);

        File amtDBFile = amtDatabaseThisMatch.getAmtDBSourceFile();
        String amtDBFileName = "";
        if (amtDBFile != null)
            amtDBFileName = amtDBFile.getName();
        AmtExtraInfoDef.setFeatureSetMatchedDatabaseName(ms1FeatureSetToMatch,
                                                         amtDBFileName);

        ms1FeatureSetToMatch.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());
        if (ms2ModificationsArray != null)
        {
            MS2ExtraInfoDef.setFeatureSetModifications(ms1FeatureSetToMatch,
                                                       ms2ModificationsArray);
        }
        if (recordDBRunsInFeatureFile)
        {
            AmtRunEntry[] runs = amtDatabaseThisMatch.getRuns();
            List<String> runsMatchedList = new ArrayList<String>(runs.length);
            for (AmtRunEntry run : runs)
                runsMatchedList.add(run.getPepXmlFilename());
            AmtExtraInfoDef.setFeatureSetRunsMatched(ms1FeatureSetToMatch, runsMatchedList);
        }
        if (matchingOutputFile != null)
        {
            try
            {
                ms1FeatureSetToMatch.save(matchingOutputFile);
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Error saving output featureset",e);
            }
        }

        return matchedMS1Features;
    }





    /**
     * Perform a mass-only match between MS1 features and AMT database.  This is
     * for determining alignment parameters when no embedded MS2 are available
     * @param amtDatabaseFeatureSet
     * @param ms1Features
     * @param massMatchDeltaMass
     * @param massMatchDeltaMassType
     * @return
     */
    public Pair<Feature,Feature>[] matchOnMassOnly(FeatureSet amtDatabaseFeatureSet,
           Feature[] ms1Features,
           float massMatchDeltaMass, int massMatchDeltaMassType
           )
    {
        //deltaElution has to be set to a hugely high number.  And it must not be set to
        //900000 (five zeroes), because that is a sentinel value.  Booo.
        ClusteringFeatureSetMatcher massOnlyFeatureSetMatcher =
                new ClusteringFeatureSetMatcher(massMatchDeltaMass,
                        massMatchDeltaMassType,
                        9000000);
        massOnlyFeatureSetMatcher.setElutionMode(BaseFeatureSetMatcherImpl.ELUTION_MODE_TIME);
        massOnlyFeatureSetMatcher.setMassBucketIncrement(1);
        massOnlyFeatureSetMatcher.setNumMassBuckets(
                Math.round(massMatchDeltaMass) - 1);
        massOnlyFeatureSetMatcher.setNumElutionBuckets(1);

        FeatureSetMatcher.FeatureMatchingResult massMatchingResult =
                massOnlyFeatureSetMatcher.matchFeatures(
                        new FeatureSet(ms1Features),
                        amtDatabaseFeatureSet);

        Set<Feature> ms1MatchedFeatures = massMatchingResult.getMasterSetFeatures();

        int i=0;
        Pair<Feature,Feature>[] result =
            (Pair<Feature,Feature>[]) new Pair[ms1MatchedFeatures.size()];

        for (Feature ms1Feature : ms1MatchedFeatures)
        {
            result[i++] =
                    new Pair<Feature,Feature>(ms1Feature,
                            massMatchingResult.getBestMatch(ms1Feature));
        }


        return result;
    }


    /**
     * Cover method
     * @param amtDatabaseFeatureSet
     * @param ms1Features
     * @param degree
     * @return
     */
    public double[] calculateTHMapCoefficientsWithMassMatching(
            FeatureSet amtDatabaseFeatureSet,
            Feature[] ms1Features, int degree)
    {
         return calculateTHMapCoefficientsWithMassMatching(amtDatabaseFeatureSet,
                 ms1Features, degree, false);
    }


    /**
     * Use mass-only matching to calculate a time-hydrophobicity map.
     *
     * First we perform a simple mass-only matching.  Plotted on a time-hydrophobicity
     * chart, this will generally produce a uniformly distributed, diffuse noise of
     * false positives, with a high-density line of true positives somewhere in the
     * middle.  We then look for that line.
     *
     * Simple regression won't do, because of all the noise -- it will tend to get
     * dragged down unpredictably. We use Yan's modal regression to determine the
     * high-density line.  This requires an R installation, with the quantreg package
     * installed.
     *
     * if onlyUseSamePeptideMatches is true, then we're doing same-peptide matching,
     * but still checking mass (because of potential effects on RT due to modifications)
     *
     * @param amtDatabaseFeatureSet
     * @param features
     * @param degree
     * @return
     */
    public double[] calculateTHMapCoefficientsWithMassMatching(
           FeatureSet amtDatabaseFeatureSet,
           Feature[] features, int degree,
           boolean onlyUseSamePeptideMatches)
    {
        _log.debug("calculateTHMapCoefficientsWithMassMatching, db features: " +
                amtDatabaseFeatureSet.getFeatures().length + ", other features: " +
                features.length);
        Pair<Feature,Feature>[] massMatchedFeatures =
            matchOnMassOnly(amtDatabaseFeatureSet, features,
                            massMatchDeltaMass, massMatchDeltaMassType);
        if (buildCharts)
        {
////This chart is maybe overkill            
//            double[] deltaMassesPPM = new double[massMatchedFeatures.length];
//            for (int i=0; i<deltaMassesPPM.length; i++)
//            {
//                deltaMassesPPM[i] = (massMatchedFeatures[i].second.getMass() -
//                        massMatchedFeatures[i].first.getMass()) *
//                        1000000 / massMatchedFeatures[i].first.getMass();
//            }
//            PanelWithHistogram deltaMassHist = new PanelWithHistogram(deltaMassesPPM);
//            deltaMassHist.setName("massMatchDeltaPPM");
//            deltaMassHist.displayInTab();
        }
        _log.debug("calculateTHMapCoefficientsWithMassMatching, mass matches = " + massMatchedFeatures.length);
        if (onlyUseSamePeptideMatches)
        {
            List<Pair<Feature,Feature>> restrictedMassMatches =
                    new ArrayList<Pair<Feature,Feature>>();

            Comparator<Pair<Feature,Feature>> ms2TimeFeaturePairComparatorAsc =
                    new Comparator<Pair<Feature,Feature>>()
            {
                public int compare(Pair<Feature,Feature> o1,Pair<Feature,Feature> o2)
                {
                    float scan1 = o1.first.getScan();
                    float scan2 = o2.first.getScan();

                    return (scan1 > scan2 ? 1 : scan2 > scan1 ? -1 : 0);
                }
            };

            Arrays.sort(massMatchedFeatures, ms2TimeFeaturePairComparatorAsc);

            Set<String> alreadyMatchedPeptides = new HashSet<String>();
            for (Pair<Feature,Feature> featurePair : massMatchedFeatures)
            {
                String firstPeptide = MS2ExtraInfoDef.getFirstPeptide(featurePair.first);
                if (firstPeptide == null)
                    continue;

                if (firstPeptide.equals(MS2ExtraInfoDef.getFirstPeptide(featurePair.second)))
                {
                    if (alreadyMatchedPeptides.contains(firstPeptide))
                        continue;
                    alreadyMatchedPeptides.add(firstPeptide);

                    restrictedMassMatches.add(featurePair);
                }
            }

            massMatchedFeatures = restrictedMassMatches.toArray(new Pair[restrictedMassMatches.size()]);
            _log.debug("Restricted matches by peptide");
        }

        _log.debug("calculateTHMapWithMassMatching: mass matches: " +
                massMatchedFeatures.length);

        return calculateTHMapCoefficientsWithMatchedFeatures(massMatchedFeatures, degree);
    }

    /**
     * Given pairs of matched features, map RT to H
     * @param matchedFeatures
     * @param degree
     * @return
     */
    public double[] calculateTHMapCoefficientsWithMatchedFeatures(
           Pair<Feature,Feature>[] matchedFeatures,
           int degree)
    {
        if (matchedFeatures.length < MIN_MATCHED_FEATURES_FOR_REGRESSION)
        {
            ApplicationContext.infoMessage("ERROR: Insufficient mass-matched features ( " +
                    matchedFeatures.length + ") to build T->H map. If you're " +
                    "seeing this error, it means that mass-matching features to the database in order to build the " +
                    "map didn't go well.  You can try increasing your mass tolerance with the 'massmatchdeltamass' " +
                    "parameter.  Please also make sure you have specified the appropriate modifications for your " +
                    "MS1 features using the 'modifications' argument.  Best of all would be to provide embedded " +
                    "MS/MS search results in PepXML format, with the 'embeddedms2' argument.");
            return null;
        }

        //create a t->H mapping for this MS1 set based on the matching results
        double[] ms1Times = new double[matchedFeatures.length];
        double[] amtHydrophobicities = new double[matchedFeatures.length];



        Feature[] dummyFeatures = new Feature[matchedFeatures.length];
        int i=0;
        boolean nonzeroMs1TimesExist = false;
        for (Pair<Feature,Feature> matchedPair : matchedFeatures)
        {
            Feature ms1Feature = matchedPair.first;
            Feature matchedMs2Feature = matchedPair.second;
            Feature dummyFeature = (Feature) ms1Feature.clone();
            ms1Times[i] = ms1Feature.getTime();
            if (ms1Times[i] > 0)
                nonzeroMs1TimesExist = true;
            amtHydrophobicities[i] = AmtExtraInfoDef.getObservedHydrophobicity(matchedMs2Feature);
            AmtExtraInfoDef.setObservedHydrophobicity(dummyFeature,
                    amtHydrophobicities[i]);
            dummyFeatures[i] = dummyFeature;
            i++;
        }
        if (!nonzeroMs1TimesExist)
            throw new IllegalArgumentException("ERROR! MS1 feature retention times are all zero!  " +
                    "Please populate feature " +
                    "times, using the 'populatems2times' command, or provide associated mzXML file(s) to " +
                    "use for populating scan times.");
//ScatterPlotDialog spdt = new ScatterPlotDialog(ms1Times, amtHydrophobicities, "mass-only matches");
//spdt.setVisible(true);
        

        double[] resultCoefficients;
        Feature[] featuresForRegression = dummyFeatures;
        double[] ms1TimesForRegression = null;
        double[] hydrophobicitiesForRegression = null;

        //The last argument has a big effect on the initial regression results, and thus on what
        //datapoints get excluded
        featuresForRegression =
                ProteomicsRegressionUtilities.selectFeaturesWithLowLeverageAndStudentizedResidual(
                        dummyFeatures,
                        ms1Times, amtHydrophobicities,
                        maxRegressionLeverageNumerator,
                        maxRegressionStudRes, false, 1, false, true);
        ApplicationContext.infoMessage("Using " + featuresForRegression.length +
                " features (out of " + matchedFeatures.length +
                " mass-matched) for regression");


        ms1TimesForRegression = new double[featuresForRegression.length];
        hydrophobicitiesForRegression = new double[featuresForRegression.length];

        for (int j=0; j<featuresForRegression.length; j++)
        {
            ms1TimesForRegression[j] = featuresForRegression[j].getTime();
            hydrophobicitiesForRegression[j] =
                    AmtExtraInfoDef.getObservedHydrophobicity(featuresForRegression[j]);
        }

        try
        {
            resultCoefficients = RegressionUtilities.modalRegression(ms1TimesForRegression,
                    hydrophobicitiesForRegression,
                    degree);
        }
        catch (IOException e)
        {
            e.printStackTrace(System.err);
            //if we failed, it could be because of timeout, or because of a missing
            //quantreg package, or...?
            //todo: move text to bundle
            throw new RuntimeException("ERROR: Failure calling R for modal regression.  R may have timed out.\n" +
                    "This may also be because the required \"quantreg\" package is not installed.\n"+
                    "  To install this package, run the following lines in R:\n" +
                    "source(\"http://bioconductor.org/biocLite.R\")\n" +
                    "biocLite(c(\"quantreg\"))");
        }


        _log.debug("Mapping coefficients:");
        for (int j=0; j<resultCoefficients.length; j++)
        {
            _log.debug("\tDegree " + j + ": " + resultCoefficients[j]);
        }

        if (buildCharts)
        {
            int maxTime = 0;
            int minTime = Integer.MAX_VALUE;
            for (double ms1Time : ms1Times)
            {
                if (ms1Time > maxTime)
                    maxTime = (int) ms1Time;
                if (ms1Time < minTime)
                    minTime = (int) ms1Time;
            }
            PanelWithScatterPlot psp = new PanelWithScatterPlot();

            for (int j=0; j<featuresForRegression.length; j++)
            {
                Feature feature = featuresForRegression[j];
                ms1TimesForRegression[j] = feature.getTime();
                hydrophobicitiesForRegression[j] = AmtExtraInfoDef.getObservedHydrophobicity(feature);
            }

            int numDotsOnChart = (maxTime-minTime+1) / 2;
            double[] chartXVals = new double[numDotsOnChart];
            double[] chartYVals = new double[numDotsOnChart];

            for (int j=0; j<numDotsOnChart; j++)
            {
                chartXVals[j] = minTime + (2 * j);
                chartYVals[j] =
                        RegressionUtilities.mapValueUsingCoefficients(resultCoefficients, chartXVals[j]);
            }
            psp.addData(ms1TimesForRegression, hydrophobicitiesForRegression,
                        "Matches used in regression");            
            psp.addData(ms1Times, amtHydrophobicities, "all mass matches");
            psp.addData(chartXVals, chartYVals, "regression function");

            psp.setAxisLabels("MS1 time","AMT Hydrophobicity");
            timeHydrophobicityMappingChart = psp.getChart();
        }

        return resultCoefficients;
    }



    /**
     * Outermost method for calculating feature hydrophobicities.  Delegator.
     * @param ms1Features
     * @param alignmentGuideFeatures
     * @param amtDatabaseFeatureSet
     * @param degree
     */
    public void calculateFeatureHydrophobicities(
            Feature[] ms1Features, Feature[] alignmentGuideFeatures,
            FeatureSet amtDatabaseFeatureSet, int degree)
    {
        if (alignmentGuideFeatures == null)
        {
            _log.debug("Using mass-only matching to map time to hydrophobicity.  Deltamass: " +
                    massMatchDeltaMass + ", delta mass type: " + massMatchDeltaMassType);
            timeHydMapCoefficients =
                    calculateTHMapCoefficientsWithMassMatching(
                            amtDatabaseFeatureSet, ms1Features, degree,
                            false);
        }
        else
        {
            timeHydMapCoefficients =
                    calculateTHMapCoefficientsWithMassMatching(
                            amtDatabaseFeatureSet, alignmentGuideFeatures, degree,
                            true);
        }        

        int maxTimePlusOne = 0;
        for (Feature feature : ms1Features)
        {
            maxTimePlusOne = Math.max(maxTimePlusOne, (int) (feature.getTime() + 1));
        }

        for (Feature feature : ms1Features)
        {
            AmtExtraInfoDef.setObservedHydrophobicity(feature,
                    RegressionUtilities.mapValueUsingCoefficients(timeHydMapCoefficients, feature.getTime()));
        }
    }


    /**
     * Calls Window2DFeatureSetMatcher to do a dead-simple 2D match between
     * two featuresets
     * @param ms1FeatureSet
     * @param amtFeatureSet
     * @param lowMassTolerance
     * @param highMassTolerance
     * @param lowHTolerance
     * @param highHTolerance
     * @return
     */
    public FeatureSetMatcher.FeatureMatchingResult
            callWindowMatcher(              FeatureSet ms1FeatureSet,
                                            FeatureSet amtFeatureSet,
                                            float lowMassTolerance,
                                            float highMassTolerance,
                                            float lowHTolerance,
                                            float highHTolerance)
    {
        Window2DFeatureSetMatcher window2DFeatureSetMatcher =
                new Window2DFeatureSetMatcher();
        window2DFeatureSetMatcher.setMatchingParameters(
                lowMassTolerance, highMassTolerance,
                lowHTolerance, highHTolerance,
                FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
        return window2DFeatureSetMatcher.matchFeatures(ms1FeatureSet, amtFeatureSet);
    }

    /**
     * Calibrate MS1 features based on an initial match to the AMT database.
     * Warning!  This will invalidate the hash keys in matchingResult when it updates the
     * features.   Do NOT try to use matchingResult after calling this.
     * @param ms1Features
     * @param matchingResult
     * @param showCharts
     * @return
     */
    public double[] calibrateMS1FeaturesWithMatches(
            Feature[] ms1Features,
            FeatureSetMatcher.FeatureMatchingResult matchingResult,
            boolean showCharts)
    {
        _log.debug("calibrateMS1FeaturesWithMatches 1");
        double[] result = calculateMassCalibrationParameters(matchingResult, showCharts);
        double slope = result[1];
        double intercept = result[0];

        if (showCharts)
        {
            int numMatches = matchingResult.size();
            double[] ms1FeatureMasses = new double[numMatches];
            double[] massErrorData = new double[numMatches];

            int i=0;            
            for (Feature ms1Feature : matchingResult.getMasterSetFeatures())
            {

                ms1FeatureMasses[i] =  ms1Feature.getMass() -
                    (float) (ms1Feature.getMass() * slope + intercept);
                List<Feature> matchedFeatures = matchingResult.get(ms1Feature);
                Feature matchedAmtFeature =  matchedFeatures.get(0);
                massErrorData[i] = ms1FeatureMasses[i] - matchedAmtFeature.getMass();
                i++;
            }

            //scatterplot of mass vs. deltaMass
            PanelWithScatterPlot psp = new PanelWithScatterPlot(ms1FeatureMasses, massErrorData,
                    "MS1 feature mass vs. (signed) match mass error");
            psp.setName("After calibration");
            psp.setAxisLabels("MS1 Feature Mass", "Absolute (Da) error (MS1 - AMT)");
            psp.displayInTab();
        }

        //Warning!  This will invalidate the hash keys in matchingResult
        for (Feature ms1Feature : ms1Features)
        {
            ms1Feature.setMass(ms1Feature.getMass() -
                    (float) (ms1Feature.getMass() * slope + intercept));
            ms1Feature.updateMz();
        }

        _log.debug("calibrateMS1FeaturesWithMatches, recalibrated");

        return result;
    }

    /**
     * Calculate linear mass calibration parameters by doing robust linear regression of
     * mass error in the matches vs. MS1 feature mass.
     * @param matchingResult
     * @param showCharts
     * @return
     */
    public double[] calculateMassCalibrationParameters(
            FeatureSetMatcher.FeatureMatchingResult matchingResult,
            boolean showCharts)
    {
        int numMatches = matchingResult.size();

        double[] ms1FeatureMasses = new double[numMatches];
        double[] massErrorData = new double[numMatches];

        int i=0;
        double minMass = Double.MAX_VALUE;
        double maxMass = Double.MIN_VALUE;
        for (Feature ms1Feature : matchingResult.getMasterSetFeatures())
        {
            ms1FeatureMasses[i] = ms1Feature.getMass();
            if (ms1FeatureMasses[i] < minMass)
                minMass = ms1FeatureMasses[i];
            if (ms1FeatureMasses[i] > maxMass)
                maxMass = ms1FeatureMasses[i];

            Feature matchedAmtFeature =  matchingResult.get(ms1Feature).get(0);

            massErrorData[i] = ms1Feature.getMass() - matchedAmtFeature.getMass();

            i++;
        }

        double[] result = RegressionUtilities.robustRegression(ms1FeatureMasses, massErrorData);
        _log.debug("calculateMassCalibrationParameters, slope=" + result[1] + ", intercept=" + result[0]);


        if (showCharts)
        {
            //scatterplot of mass vs. deltaMass
            PanelWithScatterPlot psp = new PanelWithScatterPlot(ms1FeatureMasses, massErrorData,
                                                          "MS1 feature mass vs. (signed) match mass error");
            psp.addLine(result[1], result[0], minMass, maxMass);
            psp.setName("Before calibration");
            psp.setAxisLabels("MS1 Feature Mass", "Absolute (Da) error (MS1 - AMT)");
            psp.displayInTab();

        }

        return result;
    }


    /**
     * Separating this out so the interesting code flows better
     * @param matchingResult
     */
    public void createMassTimeErrorPlots(
            FeatureSetMatcher.FeatureMatchingResult matchingResult)
    {

        int numUnambiguousMatches = 0;
        for (Feature ms1Feature : matchingResult.getMasterSetFeatures())
        {
            if (matchingResult.get(ms1Feature).size() == 1)
                numUnambiguousMatches++;
        }
        double[] ms1FeatureMasses = new double[numUnambiguousMatches];
        double[] ms1FeatureHydrophobicities = new double[numUnambiguousMatches];
        double[] massErrorData = new double[numUnambiguousMatches];
        double[] elutionErrorData = new double[numUnambiguousMatches];
        int i=0;
        for (Feature ms1Feature : matchingResult.getMasterSetFeatures())
        {
            ms1FeatureMasses[i] = ms1Feature.getMass();
            ms1FeatureHydrophobicities[i] = AmtExtraInfoDef.getObservedHydrophobicity(ms1Feature);

            if (matchingResult.get(ms1Feature).size() > 1)
                continue;

            Feature matchedAmtFeature =  matchingResult.get(ms1Feature).get(0);
            //convert to ppm
            massErrorData[i] =
                    (ms1Feature.getMass() -
                            matchedAmtFeature.getMass()) *
                            (1000000 / ms1Feature.getMass());
//System.err.println("Error: " + histogramData[i-1] + ", " +ms1Feature.getMass() + ", " + result.get(ms1Feature).get(0).getMass());

            elutionErrorData[i] =
                    (AmtExtraInfoDef.getObservedHydrophobicity(ms1Feature) -
                            (AmtExtraInfoDef.getObservedHydrophobicity(matchedAmtFeature)));
            i++;
        }


        //3D mass-elution histogram
        JDialog perspDialog = new JDialog();
        perspDialog.setSize(1000,800);
        massTimeErrorPerspectivePlot = new PanelWithRPerspectivePlot();
        massTimeErrorPerspectivePlot.setChartHeight(800);
        massTimeErrorPerspectivePlot.setChartWidth(1000);

        massTimeErrorPerspectivePlot.setSize(1000,800);
        massTimeErrorPerspectivePlot.setTiltAngle(25);
        massTimeErrorPerspectivePlot.setRotationAngle(-30);
        massTimeErrorPerspectivePlot.setAxisRVariableNames("Hydrophobicity","Mass_ppm","Matches");
        double xBinSize = .002;
        double yBinSize = 1;
        massTimeErrorPerspectivePlot.plotPointsSummary(elutionErrorData, massErrorData,
                xBinSize, yBinSize);
        perspDialog.add(massTimeErrorPerspectivePlot);
        perspDialog.setVisible(true);

        //scatterplot of mass vs. deltaMass
        ScatterPlotDialog spd = new ScatterPlotDialog(ms1FeatureMasses, massErrorData,
                "MS1 feature mass vs. (signed) match mass error");
        spd.setAxisLabels("MS1 Feature Mass", "PPM error (MS1 - AMT)");
        massDeltaMassScatterPlot = spd.getPanelWithScatterPlot().getChart();
        spd.setVisible(true);

        //scatterplot of hydrophobicity error vs. mass error
        ScatterPlotDialog spdHandM = new ScatterPlotDialog(elutionErrorData, massErrorData,
                "MS1 feature mass error vs. H error");
        spdHandM.setAxisLabels("H error", "PPM error");
        spdHandM.setVisible(true);
    }


    /**
     * Utility method
     * @param features
     * @return
     */
    public static Set<String> createPeptideSetFromFeatures(Feature[] features)
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
     * Remove peptide entries from the AMT database that don't occur in runs that are similar,
     * in peptide content, to ms2Features.
     * @param amtDatabase
     * @param ms2Features
     * @param matchingResult
     * @param showCharts
     * @return
     */
    public AmtDatabase reduceDatabaseByRunSimilarity(
            AmtDatabase amtDatabase,
            Feature[] ms2Features,
            FeatureSetMatcher.FeatureMatchingResult matchingResult,
            int minRunsToKeep, int maxRunsToKeep,
            boolean showCharts)
    {
        ApplicationContext.infoMessage("Removing unlikely entries from database");
        float[] percentMatchedPerRun = new float[amtDatabase.numRuns()];
        AmtRunEntry[] runEntries = amtDatabase.getRuns();

        final Map<AmtRunEntry, Float> runPeptideOverlapPercentMap =
                new HashMap<AmtRunEntry, Float>();
        final Map<AmtRunEntry, Integer> runPeptideOverlapCountMap =
                new HashMap<AmtRunEntry, Integer>();
        final Map<AmtRunEntry, Integer> runMatchedPeptideCountMap =
                new HashMap<AmtRunEntry, Integer>();

        final Map<AmtRunEntry, Set<String>> runMatchedPeptideMap =
                new HashMap<AmtRunEntry, Set<String>>();
        final Map<AmtRunEntry, Set<String>> runPeptideOverlapMap =
                new HashMap<AmtRunEntry, Set<String>>();
        final Map<AmtRunEntry, Set<String>> runPeptideSetMap =
                new HashMap<AmtRunEntry, Set<String>>();

        Set<String> ms2Peptides = createPeptideSetFromFeatures(ms2Features);

        Set<Feature> matchedAmtFeatures = matchingResult.getSlaveSetFeatures();
        Set<String> matchedAmtPeptides = new HashSet<String>();
        for (Feature matchedAmtFeature : matchedAmtFeatures)
        matchedAmtPeptides.add(MS2ExtraInfoDef.getFirstPeptide(matchedAmtFeature));

        for (int i=0; i<runEntries.length; i++)
        {
            if (i % (Math.max(runEntries.length/10,1)) == 0)
            {
                ApplicationContext.setMessage(Rounder.round(((double) i * 100.0 / (double) runEntries.length),0) +
                        "% done evaluating runs");
            }
            AmtRunEntry runEntry = runEntries[i];

            Set<String> matchedPeptidesInThisRun = new HashSet<String>();

            Set<String> peptidesThisRun = new HashSet<String>();
            Set<String> ms2PeptidesInCommon = new HashSet<String>();
            for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
            {
                AmtPeptideEntry.AmtPeptideObservation obs =
                        peptideEntry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    String peptideSequence = peptideEntry.getPeptideSequence();
                    peptidesThisRun.add(peptideSequence);
                    if (ms2Peptides.contains(peptideSequence))
                        ms2PeptidesInCommon.add(peptideSequence);
                    if (matchedAmtPeptides.contains(peptideSequence))
                        matchedPeptidesInThisRun.add(peptideSequence);
                }
            }

            float peptideMatchesPercent =
                    (((float) ms2PeptidesInCommon.size()) /
                            ((float) peptidesThisRun.size()) * 100);
//            System.err.println("Matches this run: " + numMassMatchedFeatures + " (" + massMatchesPercent + "%)");
            percentMatchedPerRun[i] = peptideMatchesPercent;

            runPeptideOverlapPercentMap.put(runEntry, peptideMatchesPercent);
            runPeptideOverlapCountMap.put(runEntry, ms2PeptidesInCommon.size());

            runMatchedPeptideCountMap.put(runEntry, matchedPeptidesInThisRun.size());

            runMatchedPeptideMap.put(runEntry, matchedPeptidesInThisRun);
            runPeptideOverlapMap.put(runEntry, ms2PeptidesInCommon);
            runPeptideSetMap.put(runEntry, peptidesThisRun);
        }

        AmtRunEntry[] runEntriesSorted = new AmtRunEntry[amtDatabase.numRuns()];
        for (int i=0; i<amtDatabase.numRuns(); i++)
            runEntriesSorted[i] = runEntries[i];

        Comparator<AmtRunEntry> runComparatorByPeptideMatchPercentDesc =
                new Comparator<AmtRunEntry> ()
        {
            public int compare(AmtRunEntry o1, AmtRunEntry o2)
            {
                float o1Matches = runPeptideOverlapPercentMap.get(o1);
                float o2Matches = runPeptideOverlapPercentMap.get(o2);

                return o1Matches < o2Matches ? 1 : o1Matches > o2Matches ? -1 : 0;
            }
        };

        Arrays.sort(runEntriesSorted, runComparatorByPeptideMatchPercentDesc);

        float[] runIndexes = new float[runEntriesSorted.length];


        float[] cumOverlapCount = new float[runEntriesSorted.length];
        float[] cumMatchCount = new float[runEntriesSorted.length];
        float[] cumOverlapMatchCountRatios = new float[runEntriesSorted.length];
        Set<String> cumOverlap = new HashSet<String>();
        Set<String> cumMatches = new HashSet<String>();

        Set<String> allPeptidesToKeep = new HashSet<String>();
        boolean stillAddingPeptides = true;

        //first derivative of the match:count ratio progression
        float[] deltaCumOverlapMatchCountRatios = new float[runEntriesSorted.length];
        float[] deltaDelta = new float[runEntriesSorted.length];

        float[] matchAdditions = new float[runEntriesSorted.length];
        for (int i=0; i<runEntriesSorted.length; i++)
        {
            cumOverlap.addAll(runPeptideOverlapMap.get(runEntriesSorted[i]));
            cumMatches.addAll(runMatchedPeptideMap.get(runEntriesSorted[i]));

            runIndexes[i] = i;
            cumOverlapCount[i] = cumOverlap.size();
            cumMatchCount[i] = cumMatches.size();
            if (i>20)
                matchAdditions[i] = cumMatchCount[i] - cumMatchCount[i-1];
            cumOverlapMatchCountRatios[i] = cumOverlapCount[i] / cumMatchCount[i];
            deltaCumOverlapMatchCountRatios[i] = cumOverlapMatchCountRatios[i] -
                    (i == 0 ? 0 : cumOverlapMatchCountRatios[i-1]);
            deltaDelta[i] = deltaCumOverlapMatchCountRatios[i] -
                    (i == 0 ? 0 : deltaCumOverlapMatchCountRatios[i-1]);
            //if, that is, we're adding a smaller proportion of peptides from MS/MS overlap to peptides
            //from AMT matches than we were a second ago...
//            if (i>0 && deltaCumOverlapMatchCountRatios[i] >= -0.01)

        }

        float[] matchAdditionsDifferences = new float[runEntriesSorted.length];
        int maxMatchAdditionsDifference = Integer.MIN_VALUE;
        int maxMatchAdditionsDifferenceIndex = 0;
        int windowSize = 10;
        for (int i=Math.max(windowSize, minRunsToKeep-1); i<Math.min(runEntriesSorted.length-windowSize, maxRunsToKeep-1); i++)
        {
            int matchAdditionsDifference = (int)
                    ((cumMatchCount[i+windowSize] - cumMatchCount[i]) - (cumMatchCount[i] - cumMatchCount[i-windowSize]));
//System.err.println(i + ", " + (cumMatchCount[i+windowSize] - cumMatchCount[i]) + " - " + (cumMatchCount[i] - cumMatchCount[i-windowSize]) + " = " + matchAdditionsDifference);
            matchAdditionsDifferences[i] = matchAdditionsDifference;
            if (matchAdditionsDifference > maxMatchAdditionsDifference)
            {
                maxMatchAdditionsDifference = matchAdditionsDifference;
                maxMatchAdditionsDifferenceIndex = i;
            }
        }
//System.err.println("Max additions difference: " + maxMatchAdditionsDifference + ", index=" + maxMatchAdditionsDifferenceIndex);
        int lastRunIndexToKeep = maxMatchAdditionsDifferenceIndex;

        if (lastRunIndexToKeep < minRunsToKeep-1)
            lastRunIndexToKeep = minRunsToKeep-1;
        if (lastRunIndexToKeep > maxRunsToKeep-1)
            lastRunIndexToKeep = maxRunsToKeep-1;
        for (int i=0; i<=lastRunIndexToKeep; i++)
            allPeptidesToKeep.addAll(runPeptideSetMap.get(runEntriesSorted[i]));




        if (showCharts)
        {
            PanelWithLineChart pwlc2 = new PanelWithLineChart(runIndexes, cumOverlapCount, "cumul overlap count");
            pwlc2.addData(runIndexes, cumMatchCount, "cumul match count");
            pwlc2.addData(runIndexes, matchAdditions, "# of new matches added this run");
            pwlc2.addData(runIndexes, matchAdditionsDifferences, "match addition differences");

            float[] cutoffLineX = new float[2];
            cutoffLineX[0] = lastRunIndexToKeep;
            cutoffLineX[1] = lastRunIndexToKeep;
            float[] cutoffLineY = new float[2];
            cutoffLineY[0] = 0;
            cutoffLineY[1] = 100;
            pwlc2.addData(cutoffLineX, cutoffLineY, "Cutoff");

            ChartDialog cd2 = new ChartDialog(pwlc2);
            cd2.setTitle("cumulative");
            cd2.setVisible(true);
        }

        if (lastRunIndexToKeep == -1)
        {
            ApplicationContext.infoMessage("Keeping all runs!");
            return amtDatabase;
        }
        else
            ApplicationContext.infoMessage("Keeping peptide entries from " + (lastRunIndexToKeep+1) + " out of " + amtDatabase.numRuns() + " runs");
        AmtDatabase result = (AmtDatabase) amtDatabase.waistDeepCopy();
        for (AmtPeptideEntry peptideEntry : result.getEntries())
            if (!allPeptidesToKeep.contains(peptideEntry.getPeptideSequence()))
                result.removeEntry(peptideEntry.getPeptideSequence());
        ApplicationContext.infoMessage("Kept " + result.numEntries() + " peptide entries, out of " + amtDatabase.numEntries());
        ApplicationContext.infoMessage("New database: " + result.toString());
        return result;
    }





    //Getters and Setters

    public float getMassMatchDeltaMass()
    {
        return massMatchDeltaMass;
    }

    public void setMassMatchDeltaMass(float massMatchDeltaMass)
    {
        this.massMatchDeltaMass = massMatchDeltaMass;
    }

    public int getMassMatchDeltaMassType()
    {
        return massMatchDeltaMassType;
    }

    public void setMassMatchDeltaMassType(int massMatchDeltaMassType)
    {
        this.massMatchDeltaMassType = massMatchDeltaMassType;
    }

    public float getRealMatchDeltaMass()
    {
        return realMatchDeltaMass;
    }

    public void setRealMatchDeltaMass(float realMatchDeltaMass)
    {
        this.realMatchDeltaMass = realMatchDeltaMass;
    }

    public int getRealMatchDeltaMassType()
    {
        return realMatchDeltaMassType;
    }

    public void setRealMatchDeltaMassType(int realMatchDeltaMassType)
    {
        this.realMatchDeltaMassType = realMatchDeltaMassType;
    }

    public float getRealMatchDeltaElution()
    {
        return realMatchDeltaElution;
    }

    public void setRealMatchDeltaElution(float realMatchDeltaElution)
    {
        this.realMatchDeltaElution = realMatchDeltaElution;
    }

    public void setBuildCharts(boolean buildCharts)
    {
        this.buildCharts = buildCharts;
    }

    public JFreeChart getMassMatchScatterplot()
    {
        return massMatchScatterplot;
    }

    public JFreeChart getTimeHydrophobicityMappingChart()
    {
        return timeHydrophobicityMappingChart;
    }

    public double getMaxRegressionLeverageNumerator()
    {
        return maxRegressionLeverageNumerator;
    }

    public void setMaxRegressionLeverageNumerator(double leverageNumeratorModalRegression)
    {
        this.maxRegressionLeverageNumerator = leverageNumeratorModalRegression;
    }

    public double getMaxRegressionStudentizedResidual()
    {
        return maxRegressionStudRes;
    }

    public void setMaxRegressionStudentizedResidual(double maxRegressionStudRes)
    {
        this.maxRegressionStudRes = maxRegressionStudRes;
    }


    public double getAmtFeatureMassAdjustment()
    {
        return amtFeatureMassAdjustment;
    }

    public void setAmtFeatureMassAdjustment(double amtFeatureMassAdjustment)
    {
        this.amtFeatureMassAdjustment = amtFeatureMassAdjustment;
    }

    public int getQuantileRegressionMinFeatures()
    {
        return quantileRegressionMinFeatures;
    }

    public void setQuantileRegressionMinFeatures(int quantileRegressionMinFeatures)
    {
        this.quantileRegressionMinFeatures = quantileRegressionMinFeatures;
    }

    public JFreeChart getMassCalibrationChart()
    {
        return massCalibrationChart;
    }

    /**
     * Define the "dimensions" of the AMT database.
     *
     * This only has any meaning for a single multi-fraction database like an IPAS, with one-
     * or (really) two-dimensional fractionation.
     * @param amtDatabaseStructure
     */
    public void defineAMTDBStructure(
            Fractionation2DUtilities.FractionatedAMTDatabaseStructure amtDatabaseStructure)
    {
        this.amtDatabaseStructure = amtDatabaseStructure;
        amtDBDimensionsDefined=true;
    }

    public Fractionation2DUtilities.FractionatedAMTDatabaseStructure getAMTDBStructure()
    {
        return amtDatabaseStructure;
    }


    /**
     * This is for paring down an AMT database by removing runs without many peptides
     * in common with the MS2 features provided.  For getting rid of runs in a heavily
     * fractionated database.  Unnecessary?
     * @param amtDB
     * @param ms2FeatureSetToMatch
     * @param maxEntries
     * @param minRunMassMatchPercent
     * @param maxRuns
     * @param showCharts
     * @param ms1Features
     * @param ms2ModificationsForMatching
     * @return The modified database, or null if it can't be done
     */
    public AmtDatabase buildAmtDatabaseForPeptideMatches(
            AmtDatabase amtDB,
            FeatureSet ms2FeatureSetToMatch,
            int maxEntries,
            int minRunMassMatchPercent,
            int maxRuns,
            boolean showCharts,
            Feature[] ms1Features,
            MS2Modification[] ms2ModificationsForMatching)
    {


            ApplicationContext.infoMessage("Discarding runs with insufficient peptide matches to MS2 features...");
            int oldNumRuns = amtDB.numRuns();

            AmtDatabase result = AmtDatabaseManager.removeRunsWithoutPeptideMatches(
                            amtDB, ms2FeatureSetToMatch.getFeatures(),
                            minRunMassMatchPercent,
                            maxEntries,
                            maxRuns,
                            amtDatabaseStructure,
                            showCharts,
                            ms1Features, ms2ModificationsForMatching, this);
            if (result.numEntries() == 0)
            {
                throw new RuntimeException("FAILURE: No runs in the AMT database match " +
                        minRunMassMatchPercent + "% of peptides from MS2 features in file " +
                        ms2FeatureSetToMatch.getSourceFile().getName());
            }
            else
            {
                ApplicationContext.infoMessage("Retained " + result.numRuns() + " out of " + oldNumRuns + " runs in AMT database that matched at least " + minRunMassMatchPercent + "% of features");
                ApplicationContext.infoMessage("Reduced DB: " + result.toString());
            }
        return result;
     }

    /**
     * This method removes runs from an AMT database based on distance in fraction
     * space from the run we're matching to.  This is purely for demonstration of
     * how well it doesn't work, should be discarded after paper is written.
     * @param amtDB
     * @param ms2FeatureSetToMatch
     * @param maxEntries
     * @param minRunMassMatchPercent
     * @param maxRuns
     * @param showCharts
     * @return  The modified database, or null if it can't be done
     */
    public AmtDatabase buildAmtDatabaseForGeographicRestriction(
            AmtDatabase amtDB,
            FeatureSet ms2FeatureSetToMatch,
            int maxEntries,
            int minRunMassMatchPercent,
            int maxRuns,
            boolean showCharts)
    {
        AmtRunEntry[] runs = amtDB.getRuns();
        String ms2FileName = ms2FeatureSetToMatch.getSourceFile().getName();
        int ms2Index = 0;
        for (ms2Index=0; ms2Index<runs.length; ms2Index++)
        {
            String runFileName =  runs[ms2Index].getPepXmlFilename();
            if (runFileName.substring(0, runFileName.indexOf(".")).equalsIgnoreCase(ms2FileName.substring(0, ms2FileName.indexOf("."))))
                break;
        }

        Pair<Integer, int[]> expAndPos = amtDatabaseStructure.calculateExperimentAndPosition(ms2Index);
        int runEntryCol = expAndPos.second[0];
        int runEntryRow = expAndPos.second[1];



            ApplicationContext.infoMessage("Discarding runs far away from MS1 run...");
            int oldNumRuns = amtDB.numRuns();

            AmtDatabase result = AmtDatabaseManager.removeRunsByStructure(
                            amtDB,
                    runEntryRow, runEntryCol, expAndPos.first,
                            maxEntries,
                            maxRuns,
                            amtDatabaseStructure,
                            showCharts);
            if (result.numEntries() == 0)
            {
                throw new RuntimeException("FAILURE: unable to satisfy maxEntries, maxRuns constraints using runs geographically close to " +
                            ms2FeatureSetToMatch.getSourceFile().getName());
            }
            else
            {
                ApplicationContext.infoMessage("Retained " + result.numRuns() + " out of " + oldNumRuns + " runs in AMT database that matched at least " + minRunMassMatchPercent + "% of features");
                ApplicationContext.infoMessage("Reduced DB: " + result.toString());
            }
        return result;
     }


    public float getMinMatchProbabilityToKeep()
    {
        return minMatchProbabilityToKeep;
    }

    public void setMinMatchProbabilityToKeep(float minMatchProbabilityToKeep)
    {
        this.minMatchProbabilityToKeep = minMatchProbabilityToKeep;
    }


    public AmtMatchProbabilityAssigner getProbabilityAssigner()
    {
        return probabilityAssigner;
    }

    public void setProbabilityAssigner(AmtMatchProbabilityAssigner probabilityAssigner)
    {
        this.probabilityAssigner = probabilityAssigner;
    }


    public boolean isDecoyMatch()
    {
        return doDecoyMatch;
    }

    public void setDecoyMatch(boolean decoyMatch)
    {
        this.doDecoyMatch = decoyMatch;
    }

    public int getMinEMIterations()
    {
        return minEMIterations;
    }

    public void setMinEMIterations(int minEMIterations)
    {
        this.minEMIterations = minEMIterations;
    }

    public int getMaxEMIterations()
    {
        return maxEMIterations;
    }

    public void setMaxEMIterations(int maxEMIterations)
    {
        this.maxEMIterations = maxEMIterations;
    }

    public float getMinSecondBestProbabilityDifference()
    {
        return minSecondBestProbabilityDifference;
    }

    public void setMinSecondBestProbabilityDifference(float minSecondBestProbabilityDifference)
    {
        this.minSecondBestProbabilityDifference = minSecondBestProbabilityDifference;
    }

    public float getMaxSecondBestProbability()
    {
        return maxSecondBestProbability;
    }

    public void setMaxSecondBestProbability(float maxSecondBestProbability)
    {
        this.maxSecondBestProbability = maxSecondBestProbability;
    }

    public boolean isUseMs1TimesForAlignment()
    {
        return useMs1TimesForAlignment;
    }

    public void setUseMs1TimesForAlignment(boolean useMs1TimesForAlignment)
    {
        this.useMs1TimesForAlignment = useMs1TimesForAlignment;
    }

    public float getMs1Ms2MassTolerancePPM()
    {
        return ms1Ms2MassTolerancePPM;
    }

    public void setMs1Ms2MassTolerancePPM(float ms1Ms2MassTolerancePPM)
    {
        this.ms1Ms2MassTolerancePPM = ms1Ms2MassTolerancePPM;
    }

    public float getMs1Ms2TimeToleranceSeconds()
    {
        return ms1Ms2TimeToleranceSeconds;
    }

    public void setMs1Ms2TimeToleranceSeconds(float ms1Ms2TimeToleranceSeconds)
    {
        this.ms1Ms2TimeToleranceSeconds = ms1Ms2TimeToleranceSeconds;
    }

    public int getNonlinearMappingPolynomialDegree()
    {
        return nonlinearMappingPolynomialDegree;
    }

    public void setNonlinearMappingPolynomialDegree(int nonlinearMappingPolynomialDegree)
    {
        this.nonlinearMappingPolynomialDegree = nonlinearMappingPolynomialDegree;
    }

    public float getMaxMatchFDRToKeep()
    {
        return maxMatchFDRToKeep;
    }

    public void setMaxMatchFDRToKeep(float maxMatchFDRToKeep)
    {
        this.maxMatchFDRToKeep = maxMatchFDRToKeep;
    }

    public int getMaxRProbAssignmentMillis() {
        return maxRProbAssignmentMillis;
    }

    public void setMaxRProbAssignmentMillis(int maxRProbAssignmentMillis) {
        this.maxRProbAssignmentMillis = maxRProbAssignmentMillis;
    }
}
