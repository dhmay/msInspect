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

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.apache.log4j.Logger;


import java.util.*;


/**
 * Methods for performing labeled quantitation on AMT data.
 *
 * This is very simple.  After AMT matching, we simply locate all peptides for which we've found
 * both the light and heavy pair (given an isotopic label definition) with the same modification
 * state and charge.  We consider that one feature and record the light area, heavy area, and ratio
 */
public class AmtLabeledQuant
{
    protected static Logger _log = Logger.getLogger(AmtLabeledQuant.class);

    //default acrylamide label
    public static final double ACRYLAMIDE_MASS_DIFF = 3.0101;
    public static final Character ACRYLAMIDE_RESIDUE = 'C';

    protected AnalyzeICAT.IsotopicLabel isotopicLabel = DEFAULT_ISOTOPIC_LABEL;
    public static final AnalyzeICAT.IsotopicLabel  DEFAULT_ISOTOPIC_LABEL =
            new AnalyzeICAT.IsotopicLabel((float) PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[ACRYLAMIDE_RESIDUE],
                    (float) (PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[ACRYLAMIDE_RESIDUE] + ACRYLAMIDE_MASS_DIFF),
                    ACRYLAMIDE_RESIDUE, 3);



    public static final MS2Modification[] DEFAULT_OTHER_MODS =
            new MS2Modification[]
            {
                    ModificationListArgumentDefinition.safeCreateModificationFromString("C57.021"),
                    ModificationListArgumentDefinition.safeCreateModificationFromString("M16V"),
                    //Acrylamide light label weight
                    ModificationListArgumentDefinition.safeCreateModificationFromString("C14.01557"),
            };
    protected MS2Modification[] otherModifications = DEFAULT_OTHER_MODS;


    //TODO: this should probably be expressed in ppm.
    //TODO: come to think of it, mass tolerance really shouldn't come into play.
    //TODO: I should be identifying the modification state in the AMT match.
    protected static final double DEFAULT_MASS_SLOP = 0.1;

    public static final float DEFAULT_MIN_RATIO_HIGH_PROB = 0.7f;
    public static final float DEFAULT_MIN_RATIO_LOW_PROB = 0.05f;


    protected float minRatioHigherProbability = DEFAULT_MIN_RATIO_HIGH_PROB;
    protected float minRatioLowerProbability = DEFAULT_MIN_RATIO_LOW_PROB;

    protected double massSlop = DEFAULT_MASS_SLOP;


    protected boolean showCharts = false;

    protected List<Float> logRatiosAllFiles = new ArrayList<Float>();

    protected boolean perCharge = true;

    float maxRatio = 15f;
    float minRatio = 1f/15f;

    public static final float DEFAULT_MIN_PROB_DIFF_WITHIN_STATE = 0.1f;
    protected float minProbDifferenceWithinState = DEFAULT_MIN_PROB_DIFF_WITHIN_STATE;

    public static final float DEFAULT_MAX_SECOND_BEST_PROB_WITHIN_STATE = 0.4f;
    protected float maxSecondBestProbWithinState = DEFAULT_MAX_SECOND_BEST_PROB_WITHIN_STATE;



    public AmtLabeledQuant()
    {

    }

    /**
     * do the actual work
     */
    public void quantitate(FeatureSet featureSet)
    {
        Map<String, Map<Integer, List<Feature>>> peptideChargeFeatureListMap =
                new HashMap<String, Map<Integer, List<Feature>>>();

        List<Float> possibleHeavyMasses = calculateHeavyLabelResidueMasses();

        for (Feature feature : featureSet.getFeatures())
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            Map<Integer, List<Feature>> thisPeptideChargeFeatureListMap =
                peptideChargeFeatureListMap.get(peptide);
            if (thisPeptideChargeFeatureListMap == null)
            {
                thisPeptideChargeFeatureListMap = new HashMap<Integer, List<Feature>>();
                peptideChargeFeatureListMap.put(peptide, thisPeptideChargeFeatureListMap);
            }

            int charge = feature.getCharge();
            List<Feature> featureList = thisPeptideChargeFeatureListMap.get(charge);
            if (featureList == null)
            {
                featureList = new ArrayList<Feature>();
                thisPeptideChargeFeatureListMap.put(charge, featureList);
            }
            featureList.add(feature);
        }

        List<Feature> featuresWithRatios = new ArrayList<Feature>();
        List<Feature> featuresWithoutRatios = new ArrayList<Feature>();

        List<Float> allLogRatios = new ArrayList<Float>();
        for (String peptide : peptideChargeFeatureListMap.keySet())
        {
            int numExpectedLabels = 1;
            numExpectedLabels = 0;
            for (int i=0; i<peptide.length(); i++)
            {
                if (peptide.charAt(i) == isotopicLabel.getResidue())
                    numExpectedLabels++;
            }
            //if no labeled residues this peptide, skip it
            if (numExpectedLabels == 0)
            {
                Map<Integer, List<Feature>> chargeFeatureMap = peptideChargeFeatureListMap.get(peptide);
                for (List<Feature> featureList : chargeFeatureMap.values())
                    featuresWithoutRatios.addAll(featureList);
                continue;
            }

            List<Feature> featuresWithRatiosThisPeptide = new ArrayList<Feature>();
            for (int charge=1; charge<10; charge++)
            {
                Map<Integer, List<Feature>> thisPeptideChargeFeatureListMap =
                        peptideChargeFeatureListMap.get(peptide);

                List<Feature> featureList = thisPeptideChargeFeatureListMap.get(charge);

                //no features this peptide this charge
                if (featureList == null)
                    continue;

                //Take out all features that don't pass the lower probability threshold
                List<Feature> tooLowProbFeatureList = new ArrayList<Feature>();
                for (Feature feature : featureList)
                    if (MS2ExtraInfoDef.getPeptideProphet(feature) < minRatioLowerProbability)
                        tooLowProbFeatureList.add(feature);
                featuresWithoutRatios.addAll(tooLowProbFeatureList);
                featureList.removeAll(tooLowProbFeatureList);

                if (featureList.size() == 1)
                {
                    featuresWithoutRatios.add(featureList.get(0));
                    continue;
                }

                List<Feature> lowFeatureList = new ArrayList<Feature>();
                List<Feature> highFeatureList = new ArrayList<Feature>();

                for (Feature feature : featureList)
                {
                    int numLabels = 0;

                    List<ModifiedAminoAcid>[] modifiedAcids =
                            MS2ExtraInfoDef.getModifiedAminoAcids(feature);
                    if (modifiedAcids != null)
                    {
                        for (int i=0; i<modifiedAcids.length; i++)
                        {
                            if (peptide.charAt(i) == isotopicLabel.getResidue())
                            {
                                if (modifiedAcids[i] != null)
                                {
                                    for (ModifiedAminoAcid mod : modifiedAcids[i])
                                    {
//System.err.println("\tGot mod: " + mod + ", mass: " + mod.getMass() + ", expected: " + possibleHeavyMasses.get(0));
//System.err.println("\tGot mod: " + mod + ", mass: " + mod.getMass() + ", expected: " + labelHeavyMass);

                                        for (float possibleHeavyMass : possibleHeavyMasses)
                                            if (Math.abs(mod.getMass() - possibleHeavyMass) < massSlop)
                                            {
                                                numLabels++;
                                                break;
                                            }
                                    }
                                }
                            }
                        }
                    }

                    if (numLabels == 0)
                        lowFeatureList.add(feature);
                    else if (numLabels == numExpectedLabels)
                        highFeatureList.add(feature);
                    else featuresWithoutRatios.add(feature);
                }
//System.err.println("\t" + lowFeatureList.size() + ", " + highFeatureList.size());
                if (lowFeatureList.isEmpty() || highFeatureList.isEmpty())
                {
                    featuresWithoutRatios.addAll(lowFeatureList);
                    featuresWithoutRatios.addAll(highFeatureList);
                    continue;
                }

                Feature lowWinner = pickWinnerForModState(lowFeatureList);
                Feature highWinner = pickWinnerForModState(highFeatureList);
//System.err.println("Evaluating " + featureList.size() + " features...");
                if (lowWinner != null && highWinner != null)
                {
//System.err.println("\tRatio!");
                    Feature processedPair = processFeaturePair(lowWinner, highWinner);
                    if (processedPair != null)
                    {
                        double maxProb = Math.max(MS2ExtraInfoDef.getPeptideProphet(lowWinner),
                                MS2ExtraInfoDef.getPeptideProphet(highWinner));
                        MS2ExtraInfoDef.setPeptideProphet(processedPair, maxProb);
                        featuresWithRatiosThisPeptide.add(processedPair);
                        lowFeatureList.remove(lowWinner);
                        highFeatureList.remove(highWinner);
                    }
//                    else ApplicationContext.setMessage("Tossing out extreme ratio");

                    featuresWithoutRatios.addAll(lowFeatureList);
                    featuresWithoutRatios.addAll(highFeatureList);      
                }
                else
                {
                    ApplicationContext.setMessage("Can't calculate ratio for peptide " + peptide +
                            " in charge " + charge + ": " + featureList.size() + " features, " +
                            lowFeatureList.size() + " low and " + highFeatureList.size() + " high.");
                    featuresWithoutRatios.addAll(lowFeatureList);
                    featuresWithoutRatios.addAll(highFeatureList);
                }



/*
                //todo: handle more than two better
                if (featureList.size() == 2)
                {
                    Feature lightFeature = featureList.get(0);
                    Feature heavyFeature = featureList.get(1);

                    boolean addedRatio = false;
                    double lightProb = MS2ExtraInfoDef.getPeptideProphet(lightFeature);
                    double heavyProb = MS2ExtraInfoDef.getPeptideProphet(heavyFeature);

                    double minProb = Math.min(lightProb,heavyProb);
                    double maxProb = Math.max(lightProb,heavyProb);

                    int lightLabelCount = IsotopicLabelExtraInfoDef.getLabelCount(lightFeature);
                    int heavyLabelCount = IsotopicLabelExtraInfoDef.getLabelCount(heavyFeature);
System.err.println("Label counts: " + lightLabelCount + ", " + heavyLabelCount);                    


                    if (    minProb >= minRatioLowerProbability &&
                            maxProb >= minRatioHigherProbability)
                    {

                        double massDiff = heavyFeature.getMass() - lightFeature.getMass();

                        float massDiffError = (float) (massDiff - (numExpectedLabels * separation));

                        if (massDiffError < (massSlop))
                        {
                            Feature processedPair = processFeaturePair(lightFeature, heavyFeature);
                            if (processedPair != null)
                            {
                                MS2ExtraInfoDef.setPeptideProphet(processedPair, maxProb);                                
                                featuresWithRatiosThisPeptide.add(processedPair);
                                addedRatio = true;
                            }
                        }

                    }
                    if (!addedRatio)
                    {
                        featuresWithoutRatios.addAll(featureList);
                    }

                }
                else if (featureList.size() > 2)
                {
                    //OK, we've got too many features.  Still, if the two
                    //highest-probability features are separated, one light, one heavy,
                    //that's a good sign, keep the ratio.

                    //Otherwise, punt

                    float highestProb = -1f;
                    Feature featureWithHighestProb = null;
                    for (Feature feature : featureList)
                    {
                        float prob = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                        if (prob > highestProb)
                        {
                            featureWithHighestProb = feature;
                            highestProb = prob;
                        }
                    }

                    if (highestProb < minRatioHigherProbability)
                    {
                        featuresWithoutRatios.addAll(featureList);
                        continue;
                    }


                    Feature featureWithSecondHighestProb = null;
                    float secondHighestProb = -1f;

                    for (Feature feature : featureList)
                    {
                        if (feature.equals(featureWithHighestProb))
                            continue;
                        float prob = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                        if (prob > secondHighestProb)
                        {
                            featureWithSecondHighestProb = feature;
                            secondHighestProb = prob;
                        }
                    }

                    Feature lightFeature = featureWithHighestProb;
                    Feature heavyFeature = featureWithSecondHighestProb;

                    if (featureWithSecondHighestProb.getMass() < featureWithHighestProb.getMass())
                    {
                        heavyFeature = featureWithHighestProb;
                        lightFeature = featureWithSecondHighestProb;
                    }

                    double massDiff = heavyFeature.getMass() - lightFeature.getMass();

                    float massDiffError = (float) (massDiff - (numExpectedLabels * separation));

                    if (massDiffError < (massSlop))
                    {
                        Feature processedPair = processFeaturePair(lightFeature, heavyFeature);
                        if (processedPair != null)
                            featuresWithRatiosThisPeptide.add(processedPair);
                    }
                    else
                    {

                        ApplicationContext.infoMessage(featureList.size() +
                                " features for peptide " + peptide + " in charge " + charge + ", punting");

                        featuresWithoutRatios.addAll(featureList);
                    }
                }
*/
            }


            if (featuresWithRatiosThisPeptide.size() > 0)
            {
                for (Feature feature : featuresWithRatiosThisPeptide)
                    allLogRatios.add((float) Math.log(IsotopicLabelExtraInfoDef.getRatio(feature)));
                if (perCharge)
                {
                    featuresWithRatios.addAll(featuresWithRatiosThisPeptide);
                    for (Feature feature : featuresWithRatiosThisPeptide)
                          allLogRatios.add((float) Math.log(IsotopicLabelExtraInfoDef.getRatio(feature)));
                }
                else
                {
                    float[] ratiosThisPeptide = new float[featuresWithRatiosThisPeptide.size()];
                    for (int i=0; i<ratiosThisPeptide.length; i++)
                        ratiosThisPeptide[i] = (float) IsotopicLabelExtraInfoDef.getRatio(featuresWithRatiosThisPeptide.get(i));
                    float ratioThisPeptide = BasicStatistics.geometricMean(ratiosThisPeptide);

                    Feature representativeFeature = featuresWithRatiosThisPeptide.get(0);
                    IsotopicLabelExtraInfoDef.setRatio(representativeFeature,ratioThisPeptide);
                    double newLightIntensity = ratioThisPeptide * IsotopicLabelExtraInfoDef.getHeavyIntensity(representativeFeature);
                    IsotopicLabelExtraInfoDef.setLightIntensity(representativeFeature, newLightIntensity);

                    featuresWithRatios.add(representativeFeature);
                    allLogRatios.add((float) Math.log(ratioThisPeptide));
                }

            }
        }

        ApplicationContext.infoMessage("Features with ratios: " + featuresWithRatios.size());

        List<Feature> allFeatures = new ArrayList<Feature>();
        allFeatures.addAll(featuresWithRatios);
        allFeatures.addAll(featuresWithoutRatios);

        featureSet.setFeatures(allFeatures.toArray(new Feature[allFeatures.size()]));
    }

    /**
     * Returns the highest-probability feature in the list, or null if no feature stands
     * out above the rest enough
     * @param featuresInModState
     * @return
     */
    protected Feature pickWinnerForModState(List<Feature> featuresInModState)
    {
        if (featuresInModState.size() == 1)
            return featuresInModState.get(0);

        float highestProb = -1f;
        float secondHighestProb = -1f;
        Feature highestProbFeature = null;

        for (Feature feature : featuresInModState)
        {
            float featureProb = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
            if (featureProb > highestProb)
            {
                secondHighestProb = highestProb;

                highestProb = featureProb;
                highestProbFeature = feature;
            }
            else if (featureProb > secondHighestProb)
            {
                secondHighestProb = featureProb;
            }
        }

        if (secondHighestProb < maxSecondBestProbWithinState &&
            (highestProb - secondHighestProb) > minProbDifferenceWithinState)
            return highestProbFeature;
        else return null;
    }

    protected Feature processFeaturePair(Feature lightFeature, Feature heavyFeature)
    {
        float lightIntensity = lightFeature.getIntensity();
        float heavyIntensity = heavyFeature.getIntensity();
        float ratio = (lightIntensity / heavyIntensity);

        if (ratio < minRatio || ratio > maxRatio)
            return null;

        IsotopicLabelExtraInfoDef.setRatio(lightFeature,ratio);
        IsotopicLabelExtraInfoDef.setLightIntensity(lightFeature, lightIntensity);
        IsotopicLabelExtraInfoDef.setHeavyIntensity(lightFeature, heavyIntensity);
        if (isotopicLabel != null)
        {
            IsotopicLabelExtraInfoDef.setLabel(lightFeature,
                    isotopicLabel);
        }

        return lightFeature;
    }


    /**
     * calculate the masses that represent a heavy label on the relevant residue
     *
     * Warning -- we're not handling the case in which there are THREE variable modifications --
     * the label and two others.  But, come on, that's, you know, crazy.
     * @return
     */
    protected List<Float> calculateHeavyLabelResidueMasses()
    {
        List<Float> lightLabelMasses = new ArrayList<Float>();

        float lightMass = (float) PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[isotopicLabel.getResidue()];

        _log.debug("calculateHeavyLabelResidueMasses. Monoisotopic mass of " +
                isotopicLabel.getResidue() + ": " + lightMass);

        for (MS2Modification mod : otherModifications)
        {
            if (!mod.getVariable() && mod.getAminoAcid().charAt(0) == isotopicLabel.getResidue())
                lightMass += mod.getMassDiff();
        }

        _log.debug("calculateHeavyLabelResidueMasses. Mass of " +
                isotopicLabel.getResidue() + " with static mods: " + lightMass);        

        lightLabelMasses.add(lightMass);

        for (MS2Modification mod : otherModifications)
        {
            //here's where we'd need to add code to handle combinations of multiple var mods
            if (mod.getVariable() && mod.getAminoAcid().charAt(0) == isotopicLabel.getResidue())
                lightLabelMasses.add(lightMass + mod.getMassDiff());
        }

        //I'm lazy.  the var "lightLabelMasses" now will actually represent heavy label masses

        for (int i=0; i<lightLabelMasses.size(); i++)
            lightLabelMasses.set(i, lightLabelMasses.get(i) +
                                 (isotopicLabel.getHeavy() - isotopicLabel.getLight()));
        _log.debug("calculateHeavyLabelResidueMasses. Heavy masses: ");
        for (float mass : lightLabelMasses)
            _log.debug("\t" + mass);

        return lightLabelMasses;
    }


    public AnalyzeICAT.IsotopicLabel getIsotopicLabel()
    {
        return isotopicLabel;
    }

    public void setIsotopicLabel(AnalyzeICAT.IsotopicLabel isotopicLabel)
    {
        this.isotopicLabel = isotopicLabel;
    }

    public double getMassSlop()
    {
        return massSlop;
    }

    public void setMassSlop(double massSlop)
    {
        this.massSlop = massSlop;
    }


    public float getMinRatioHigherProbability()
    {
        return minRatioHigherProbability;
    }

    public void setMinRatioHigherProbability(float minRatioHigherProbability)
    {
        this.minRatioHigherProbability = minRatioHigherProbability;
    }

    public float getMinRatioLowerProbability()
    {
        return minRatioLowerProbability;
    }

    public void setMinRatioLowerProbability(float minRatioLowerProbability)
    {
        this.minRatioLowerProbability = minRatioLowerProbability;
    }

    public MS2Modification[] getOtherModifications()
    {
        return otherModifications;
    }

    public void setOtherModifications(MS2Modification[] otherModifications)
    {
        this.otherModifications = otherModifications;
    }
}
