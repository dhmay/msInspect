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


import org.fhcrc.cpl.toolbox.gui.chart.PanelWithRPerspectivePlot;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;

import java.util.*;
import org.apache.log4j.Logger;

import javax.swing.*;

/**
 * This class handles deconvolution of features; that is, identifying features
 * that we've found in multiple charge states and collapsing them down.  Collapsed
 * features contain the details of the original feature with highest intensity
 * (except for the intensity itself, which is the sum of all collapsed features').
 *
 * Optimization attempts to find the best mass and time tolerance values for performing
 * deconvolution.  It does this by trying a large number of possible settings and,
 * for each combination of parameters, assigning some kind of quality score.  See the
 * methods below for details.
 */
public class Deconvoluter implements Cloneable
{
    static Logger _log = Logger.getLogger(Deconvoluter.class);

    public static final double DEFAULT_DELTA_TIME = 30;
    public static final double DEFAULT_DELTA_MASS = 15;

    protected double deltaTime = DEFAULT_DELTA_TIME;
    protected double deltaMass = DEFAULT_DELTA_MASS;

    //the boundaries for the optimization parameters.
    //ideally the MAXes should be multiples of the MINs.
    public static final double DEFAULT_MAX_DELTA_MASS_PPM = 40;
    public static final double DEFAULT_MAX_DELTA_TIME = 50;
    public static final double DEFAULT_MIN_DELTA_MASS_PPM = 5;
    public static final double DEFAULT_MIN_DELTA_TIME = 10;

    //Increments between different values of the parameters to try
    protected static final double MASS_INCREMENT_PPM = 2;
    protected static final double TIME_INCREMENT_SECONDS = 5;    

    protected FeatureGrouper grouper;

    protected FeatureSet featureSet;

    protected boolean showCharts=false;

    public Deconvoluter(FeatureSet featureSet)
    {
        this.featureSet = featureSet;
        grouper = new FeatureGrouper();
        grouper.addSet(featureSet);
        grouper.setGroupByMass(true);
        setMassType(FeatureClusterer.DELTA_MASS_TYPE_PPM);
        setDeltaMass(DEFAULT_DELTA_MASS);
        setDeltaTime(DEFAULT_DELTA_TIME);
        setElutionMode(FeatureClusterer.ELUTION_MODE_TIME);
    }

    /**
     * Optimize both the mass and time matching tolerances.  As described above,
     * this is based on the ratio of multi-to-1 deconvolutions to 2-to-1 deconvolutions.
     */
    public void optimizeParameters()
    {
        int numMassIncrements = (int) ((DEFAULT_MAX_DELTA_MASS_PPM - DEFAULT_MIN_DELTA_MASS_PPM)
                                        / MASS_INCREMENT_PPM) + 1;
        int numTimeIncrements = (int) ((DEFAULT_MAX_DELTA_TIME - DEFAULT_MIN_DELTA_TIME)
                                        / TIME_INCREMENT_SECONDS) + 1;

        double[] massSettings = new double[numMassIncrements];
        double[] timeSettings = new double[numTimeIncrements];
        double[][] ratioMultiTo2EachSetting = new double[numMassIncrements][numTimeIncrements];

        for (int timeIndex = 0; timeIndex < numTimeIncrements; timeIndex++)
        {
            timeSettings[timeIndex] = timeIndex * TIME_INCREMENT_SECONDS +
                    DEFAULT_MIN_DELTA_TIME;
        }
        for (int massIndex = 0; massIndex < numMassIncrements; massIndex++)
        {
            massSettings[massIndex] = massIndex * MASS_INCREMENT_PPM +
                    DEFAULT_MIN_DELTA_MASS_PPM;
        }

        double bestDeltaMass = 0;
        double bestDeltaTime = 0;
        double highestMultiTo2Ratio = -1;
//System.err.println("deltaMass\tdeltaTime\tnumFeatures\tnumBuckets\tnum2\tsameCharge2\tdiffCharge2\tnum3\tbucketsWith3with3SameCharge\tbucketsWith3with2SameCharge\tnumMore\tbucketsWithMoreWith4OrMoreSameCharge\tbucketsWithMoreWith3SameCharge\tbucketsWithMoreWith2SameCharge");
        for (int massIndex = 0; massIndex < numMassIncrements; massIndex++)
        {
            for (int timeIndex = 0; timeIndex < numTimeIncrements; timeIndex++)
            {
                setDeltaMass(massSettings[massIndex]);
                setDeltaTime(timeSettings[timeIndex]);
                ratioMultiTo2EachSetting[massIndex][timeIndex] =
//this is for debugging                        
//                       collectBucketInfo();
//this is the best one
                        calcMultiToSingleChargeStateDifference();
//this doesn't work well
//                        calcMultiTo2ratio();
//this is downright silly
//deconvolute().getFeatures().length;
                if (ratioMultiTo2EachSetting[massIndex][timeIndex] > highestMultiTo2Ratio)
                {
                    highestMultiTo2Ratio = ratioMultiTo2EachSetting[massIndex][timeIndex];
                    bestDeltaMass = massSettings[massIndex];
                    bestDeltaTime = timeSettings[timeIndex];   
                }
            }
        }

        setDeltaMass(bestDeltaMass);
        setDeltaTime(bestDeltaTime);

        _log.debug("Optimal deconvolution parameters: deltaMass=" + bestDeltaMass + ", deltaTime=" + bestDeltaTime);

        //This plot is pretty obscure to read; mainly for development purposes
        if (showCharts)
        {
            PanelWithRPerspectivePlot plotPanel =
                    new PanelWithRPerspectivePlot();
            plotPanel.setRotationAngle(225);
            plotPanel.plot(massSettings, timeSettings,
                           ratioMultiTo2EachSetting);

            JDialog imageDialog = new JDialog();

            imageDialog.setTitle("X = deltaMass, Y = deltaTime, Z = Ratio of 3-feature collapses to 2-feature");
            imageDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            imageDialog.add(plotPanel);
            plotPanel.repaint();
            imageDialog.setSize(PanelWithRPerspectivePlot.DEFAULT_CHART_DIALOG_WIDTH,
                                PanelWithRPerspectivePlot.DEFAULT_CHART_DIALOG_HEIGHT);
            imageDialog.setVisible(true);
        }
    }


    /**
     * Here's where we calculate the ratio of multiple-features-collapsed-to-one to
     * 2-features-collapsed-to-one, for a given set of matching parameters.  We do this
     * by performing the relevant bits of deconvolution and then analyzing the clusters.
     * @deprecated
     * @return
     */
//    protected double calcMultiTo2ratio()
//    {
//        grouper.split2D(deltaMass, deltaTime);
//        Clusterer2D.BucketSummary[] buckets = grouper.summarize();
//        int numMorethan2s=0;
//        int num2s=0;
//
//        for (int i = 0; i < buckets.length; i++)
//        {
//            Clusterer2D.BucketSummary bucket = buckets[i];
//            if (bucket.featureCount > 1)
//            {
//                if (bucket.featureCount > 2)
//                    numMorethan2s++;
//                else if (bucket.featureCount ==2)
//                    num2s++;
//            }
//        }
//
//        return (double) numMorethan2s / (double) num2s;
//    }

    /**
     *
     * @return
     */
    protected double collectBucketInfo()
    {
        grouper.split2D(deltaMass, deltaTime);
        Clusterer2D.BucketSummary[] buckets = grouper.summarize();

        int numFeatures = deconvolute().getFeatures().length;
        int numBuckets = buckets.length;

        int bucketsWith2 = 0;
        int bucketsWith3 = 0;
        int bucketsWithMoreThan3 = 0;

        int bucketsWith2SameCharge = 0;
        int bucketsWith2DiffCharge = 0;
        int bucketsWith3AllSameCharge = 0;
        int bucketsWith3TwoSameCharge = 0;
        int bucketsWithMoreFourOrMoreSameCharge = 0;

        int bucketsWithMoreThreeSameCharge = 0;
        int bucketsWithMoreTwoSameCharge = 0;


        for (Clusterer2D.BucketSummary bucket : buckets)
        {
            Feature[] featuresInBucket = FeatureGrouper.getFeatures(bucket);

            int[] numInEachCharge = new int[10];
            for (int j=0; j<featuresInBucket.length; j++)
            {
                numInEachCharge[featuresInBucket[j].getCharge()]++;   
            }

            switch(bucket.featureCount)
            {
                case 0:
                    break;
                case 1:
                    break;
                case 2:
                    bucketsWith2++;
                    boolean sameCharge = false;

                    for (int i=0; i<10; i++)
                        if (numInEachCharge[i] == 2)
                        {
                            sameCharge = true;
                            break;
                        }

                    if (sameCharge)
                        bucketsWith2SameCharge++;
                    else
                        bucketsWith2DiffCharge++;
                    break;
                case 3:
                    bucketsWith3++;
                    boolean threeSameCharge = false;
                    boolean twoSameCharge = false;

                    for (int i=0; i<10; i++)
                    {
                        if (numInEachCharge[i] == 3)
                        {
                            threeSameCharge = true;
                            break;
                        }
                        else if (numInEachCharge[i] == 2)
                        {
                            twoSameCharge = true;
                            break;
                        }
                    }

                    if (twoSameCharge)
                        bucketsWith3TwoSameCharge++;
                    else if (threeSameCharge)
                        bucketsWith3AllSameCharge++;
                    break;
                default:
                    bucketsWith3++;
                    boolean moreThanThreeSameCharge = false;
                    boolean threeSameChargeAgain = false;
                    boolean twoSameChargeAgain = false;

                    for (int i=0; i<10; i++)
                    {
                        if (numInEachCharge[i] == 2)
                        {
                            twoSameCharge = true;
                        }
                        else if (numInEachCharge[i] == 3)
                        {
                            threeSameCharge = true;
                            twoSameCharge = false;
                        }
                        else if (numInEachCharge[i] > 3)
                        {
                            moreThanThreeSameCharge = true;
                            twoSameCharge = false;
                            threeSameCharge = false;                              
                            break;
                        }
                    }

                    if (moreThanThreeSameCharge)
                        bucketsWithMoreFourOrMoreSameCharge++;
                    else if (threeSameChargeAgain)
                        bucketsWithMoreThreeSameCharge++;
                    else if (twoSameChargeAgain)
                        bucketsWithMoreTwoSameCharge++;

                    break;
            }
        }
System.err.println(deltaMass + "\t" + deltaTime + "\t" + numFeatures + "\t" + numBuckets + "\t" +
        bucketsWith2 + "\t" + bucketsWith2SameCharge + "\t" + bucketsWith2DiffCharge + "\t" +
        bucketsWith3 + "\t" + bucketsWith3AllSameCharge + "\t" + bucketsWith3TwoSameCharge + "\t" +
        bucketsWithMoreThan3 + "\t" + bucketsWithMoreFourOrMoreSameCharge + "\t" + bucketsWithMoreThreeSameCharge + "\t" + bucketsWithMoreTwoSameCharge);

        return 0;
    }


    /**
     * Here's where we calculate the number of multi-feature buckets for which
     * there are no features with the same charge, and subract the number of buckets
     * for which there are two or more features with the same charge (even if there
     * are other features with different charges.  This gives us
     * an estimate of how good deconvolution looks.
     * 
     * We do this
     * by performing the relevant bits of deconvolution and then analyzing the clusters.
     * @return
     */
    protected double calcMultiToSingleChargeStateDifference()
    {
        grouper.split2D(deltaMass, deltaTime);
        Clusterer2D.BucketSummary[] buckets = grouper.summarize();
        int numMultipleCharges=0;
        int numSameCharges=0;

        for (Clusterer2D.BucketSummary bucket : buckets)
        {
            if (bucket.featureCount > 1)
            {
                Feature[] featuresInBucket = FeatureGrouper.getFeatures(bucket);
                int firstChargeState = featuresInBucket[0].getCharge();
                boolean sameChargesInBucket = false;
                for (int j=1; j<featuresInBucket.length; j++)
                {
                    if (featuresInBucket[j].getCharge() == firstChargeState)
                    {
                        sameChargesInBucket = true;
                        break;
                    }
                }
                if (sameChargesInBucket)
                    numSameCharges++;
                else
                    numMultipleCharges++;
            }
        }

        return (double) numMultipleCharges - (double) numSameCharges;
    }


    /**
     * Combine features that represent the same peptide.  Features
     * must be within a deltaMass, deltaTime window to be considered same.
     * <p/>
     * This removes redunancy caused by multiple charge states.  It
     * also combines features of the same mass/charge that are close
     * together (e.g. a small feature that is split by ion competition)
     * <p/>
     * CONSIDER(mbellew) the mass (mz) tolerance could be tighter
     * for features w/ same charge
     *
     * @return
     */
    public FeatureSet deconvolute()
    {
        _log.debug("Deconvoluting.  Delta Mass: " + deltaMass + "ppm, delta time: " + deltaTime);
        _log.debug("Before deconvolution, #features: " + featureSet.getFeatures().length);
        grouper.split2D(deltaMass, deltaTime);
        Clusterer2D.BucketSummary[] buckets = grouper.summarize();
        Feature[] deconvoluted = new Feature[buckets.length];

        int numPeptideConflicts = 0;
        int numPreservedPeptides = 0;

        for (int i = 0; i < buckets.length; i++)
        {
            Clusterer2D.BucketSummary bucket = buckets[i];
            Feature deconvolutedFeature = null;
            if (bucket.featureCount == 1)
            {
                deconvolutedFeature = (Feature) FeatureGrouper.getFeatures(bucket)[0].clone();
                deconvolutedFeature.setChargeStates(1);
            }
            else
            {
                Feature[] bucketFeatures = FeatureGrouper.getFeatures(bucket);
                Feature best = bucketFeatures[0];
                float totalIntensity = 0.0f;
                String description = "";

                for (Feature f : bucketFeatures)
                {
                    if (description.length() > 0)
                        description += ", ";
                    description += f.charge;
                    if (null != f.getDescription())
                        description += " (" + f.getDescription() + ")";
                    totalIntensity += f.intensity;
                    if (f.totalIntensity > best.totalIntensity)
                        best = f;
                }

                deconvolutedFeature = (Feature) best.clone();

                deconvolutedFeature.setIntensity(totalIntensity);
                deconvolutedFeature.setChargeStates(bucket.featureCount);
                deconvolutedFeature.setDescription(description);


                //if there's MS2 data in this FeatureSet, then we need to make sure
                //that the collapsed feature contains the peptide and protein ID carried
                //by its components.
                //If there are conflicts, we leave the existing ID on the collapsed feature
                //alone, or if it had none, don't assign
                //TODO: somehow move this to MS2ExtraInfoDef?
                if (featureSet.hasExtraInformationType(MS2ExtraInfoDef.getSingletonInstance()))
                {
                    Set<String> featurePeptides = new HashSet<String>();
                    Set<String> featureProteins = new HashSet<String>();

                    for (Feature f : bucketFeatures)
                    {
                        String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(f);
                        if (featurePeptide != null)
                        {
                            featurePeptides.add(featurePeptide);

                            String featureProtein = MS2ExtraInfoDef.getFirstProtein(f);
                            if (featureProtein != null)
                                featureProteins.add(featureProtein);
                        }
                    }

                    if (featurePeptides.size() == 1 &&
                            MS2ExtraInfoDef.getFirstPeptide(deconvolutedFeature) == null)
                    {
                        MS2ExtraInfoDef.setSinglePeptide(deconvolutedFeature,
                                featurePeptides.iterator().next());
                        numPreservedPeptides++;

                        if (featureProteins.size() == 1 &&
                                MS2ExtraInfoDef.getFirstProtein(deconvolutedFeature) == null)
                            MS2ExtraInfoDef.addProtein(deconvolutedFeature,
                                    featureProteins.iterator().next());
                    }
                    else
                    {
                        if (featurePeptides.size() > 1)
                            numPeptideConflicts++;
                    }
                }
            }
            deconvolutedFeature.comprised = FeatureGrouper.getFeatures(bucket);
            deconvoluted[i] = deconvolutedFeature;
        }

        //reporting on peptides preserved and conflicts
        if (featureSet.hasExtraInformationType(MS2ExtraInfoDef.getSingletonInstance()))
        {
            _log.debug("deconvolute: peptides actively preserved: " + numPreservedPeptides);
            _log.debug("deconvolute: peptide conflicts: " + numPeptideConflicts);
        }

        FeatureSet fs = (FeatureSet) featureSet.clone();
//Make map modifiable & set properties.
        Map props = new HashMap();
        props.putAll(featureSet.getProperties());
        if (null != featureSet.getSourceFile())
            props.put("origSourceFile", featureSet.getSourceFile());
        props.put("deconvoluteScanDiff", String.valueOf(deltaTime));
        props.put("deconvoluteMassDiff", String.valueOf(deltaMass));

        fs.setFeatures(deconvoluted);
        fs.setProperties(props);

        _log.debug("After deconvolution, #features: " + deconvoluted.length);

        return fs;
    }

    public int getMassType()
    {
        return grouper.getMassType();
    }

    public void setMassType(int massType)
    {
        grouper.setMassType(massType);
    }


    public double getDeltaTime()
    {
        return deltaTime;
    }

    public void setDeltaTime(double deltaTime)
    {
        this.deltaTime = deltaTime;
    }

    public double getDeltaMass()
    {
        return deltaMass;
    }

    public void setDeltaMass(double deltaMass)
    {
        this.deltaMass = deltaMass;
    }

    public int getElutionMode()
    {
        return grouper.getElutionMode();
    }

    public void setElutionMode(int _elutionMode)
    {
        grouper.setElutionMode(_elutionMode);
    }

    public boolean shouldShowCharts()
    {
        return showCharts;
    }

    public void setShowCharts(boolean showCharts)
    {
        this.showCharts = showCharts;
    }
}
