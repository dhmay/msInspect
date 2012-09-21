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

package org.fhcrc.cpl.viewer.metabologna;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.IOException;
import java.io.File;

/**
 * If multiple Adducts have the same chemical formula, and some of them have modifications, any modified Adducts
 * with that chemical formula will be removed from matching results.
 *
 * Change feature masses to be mz * charge, rather than (mz - H) * charge
 */
public class MetaboliteDatabaseMatcher
{
    protected static Logger _log = Logger.getLogger(MetaboliteDatabaseMatcher.class);

    public static final float DEFAULT_INITIAL_MATCH_PROPORTION = 0.1f;

    public static final float DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM = 5;
    public static final float DEFAULT_CALIBRATION_MAX_STUD_RES = 1;
    public static final float DEFAULT_CALIBRATION_MAX_LEVERAGE = 6;
    public static final int DEFAULT_CALIBRATION_REGRESSION_DEGREES_FREEDOM = 1;

    protected float calMassTolerancePPM = DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM;
    protected float calMaxLeverage = DEFAULT_CALIBRATION_MAX_LEVERAGE;
    protected double calMaxStudRes = DEFAULT_CALIBRATION_MAX_STUD_RES;
    protected int calRegressionDoF = DEFAULT_CALIBRATION_REGRESSION_DEGREES_FREEDOM;

    protected Map<ChemicalFormula, List<Adduct>> allFormulaAdductsMap;

    protected List<ChemicalCompound> databaseCompounds = null;
    protected boolean showCharts = false;

    protected boolean useUnmodifiedAdduct = true;

    //When we're calculating the ppm error that contains 'most' of the points in the 2-peak error distribution,
    //how many standard deviations is 'most'?
    protected float sigmaMultipleForErrorDist = 2;

    //Chemical modifications to apply.  Adducts will be created with /one/ of these and with /none/ of these
    //todo: implment rules about applying different modifications to types of adducts
    protected List<ChemicalModification> chemicalModifications = new ArrayList<ChemicalModification>();

    /**
     * Return mass calibration coefficients.
     * Also, do the calibration.
     * @param features
     * @return calibration coefficients
     * @throws IOException
     */
    public double[] calibrateFeatureMasses(Feature[] features) throws IOException
    {
        _log.debug("calibrating feature masses. " + features.length + " features, mass tolerance: " +
                calMassTolerancePPM + "ppm");
        //First by time, then by mass
        Map<Feature, Map<ChemicalFormula, List<Adduct>>> featureMassMatchMap =
                massMatchFull(features, calMassTolerancePPM, 1);
        List<Float> origMasses = new ArrayList<Float>();
        _log.debug(featureMassMatchMap.size() + " features matched");

        /*
        //this code was used for calibrating feature masses by retention time.  I think this was misguided, although
        //it seemed to do something useful.
        double[] rtCalibrationCoefficients = null;
        if (calWithRTFirst) {
            _log.debug("Calibrating using RT");
            if (showCharts)
            {
                showMassDeltaMassChart(featureMassMatchMap, "DeltaMassVsMass_beforeRTCal");
            }
            rtCalibrationCoefficients = calcMassCalibrationUsingMatches(featureMassMatchMap, true, true);

            if (_log.isDebugEnabled())
            {
                System.err.print("RT calibration coefficients: ");
                for (double coeff : rtCalibrationCoefficients)
                    System.err.print(" " + coeff);
                System.err.println();
            }
            List<Float> featureAbsMassChangesPPM = new ArrayList<Float>();
            for (Feature feature : features)
            {
                float featureMass = feature.getMz() * feature.getCharge();
                origMasses.add(featureMass);
                float massChangePPM = (float) RegressionUtilities.mapValueUsingCoefficients(rtCalibrationCoefficients,
                        feature.getTime());
                float massChangeDa = MassUtilities.convertPPMToDa(massChangePPM, featureMass);
                float newMass = featureMass - massChangeDa;
                featureAbsMassChangesPPM.add(Math.abs(massChangePPM));

                feature.setMass(newMass);
                feature.setMz(newMass / feature.getCharge());
                if (feature.comprised != null)
                {
                    for (Spectrum.Peak peak : feature.comprised)
                    {
                        if (peak == null)
                            continue;
                        peak.setMz(peak.mz - (massChangeDa / Math.abs(feature.charge)));
//System.err.println("old=" + oldPeakMass + ", new=" + newPeakMass + ", ppmdiff: " + MassUtilities.convertDaToPPM(newPeakMass-oldPeakMass, oldPeakMass));
                    }
                }
            }
            _log.debug("Mean absolute mass change by RT: " + Rounder.round(BasicStatistics.mean(featureAbsMassChangesPPM),1) + "ppm");
        }
        */

        _log.debug("Calibrating by mass");
        featureMassMatchMap =
                massMatchFull(features, calMassTolerancePPM, 1);
        if (showCharts)
        {
            showMassDeltaMassChart(featureMassMatchMap, "DeltaMassVsMass_beforeMassCal"   );
        }

        double[] calibrationCoefficientsByMass = calcMassCalibrationUsingMatches(featureMassMatchMap, false, false);
        if (_log.isDebugEnabled())
        {
            System.err.print("Mass calibration coefficients: ");
            for (double coeff : calibrationCoefficientsByMass)
                System.err.print(" " + coeff);
            System.err.println();
        }
        List<Float> featureAbsMassChangesByMassPPM = new ArrayList<Float>();
        for (Feature feature : features)
        {
            //paranoia -- just in case feature.getMass() hasn't been corrected for non-peptide use
            float featureMass = feature.getMz() * feature.getCharge();
            float massDiffDa = (float) RegressionUtilities.mapValueUsingCoefficients(calibrationCoefficientsByMass,
                            featureMass);
            float newMass = featureMass - massDiffDa;
            featureAbsMassChangesByMassPPM.add(MassUtilities.convertDaToPPM(Math.abs(massDiffDa), featureMass));
            feature.setMass(newMass);
            feature.setMz(newMass / feature.getCharge());
            if (feature.comprised != null)
            {
                for (Spectrum.Peak peak : feature.comprised)
                {
                    if (peak == null)
                        continue;
                    float peakMass = peak.mz * feature.getCharge();
                    float massDiffDaThisPeak = (float) RegressionUtilities.mapValueUsingCoefficients(
                            calibrationCoefficientsByMass,
                            peakMass);
                    peak.setMz(peak.mz - (massDiffDaThisPeak/Math.abs(feature.charge)));
               
//System.err.println("old=" + oldPeakMass + ", new=" + newPeakMass + ", ppmdiff: " + MassUtilities.convertDaToPPM(newPeakMass-oldPeakMass, oldPeakMass));
                }
            }
        }
        _log.debug("Mean absolute mass change by mass: " + Rounder.round(BasicStatistics.mean(featureAbsMassChangesByMassPPM),1) + "ppm");

        if (_log.isDebugEnabled())
        {
            List<Float> totalMassChangesPPM = new ArrayList<Float>();
            for (int i=0; i<features.length; i++)
                totalMassChangesPPM.add(MassUtilities.convertDaToPPM(features[i].mass - origMasses.get(i), features[i].mass));
            _log.debug("Mean absolute mass change with both calibrations: " + Rounder.round(BasicStatistics.mean(totalMassChangesPPM),1) + "ppm");
            if (showCharts)
            {
                new PanelWithHistogram(totalMassChangesPPM, "Cal PPM Changes").displayInTab();
            }
        }

        featureMassMatchMap =
                massMatchFull(features, getOrCreateFormulaAdductsMap(), calMassTolerancePPM, 1);
        if (showCharts)
        {
            showMassDeltaMassChart(featureMassMatchMap, "DeltaMassVsMass_afterMassCal");

        }

        return calibrationCoefficientsByMass;
    }

    protected void showMassDeltaMassChart(Map<Feature, Map<ChemicalFormula, List<Adduct>>> featureMassMatchMap,
                                          String name)
    {
            List<Double> masses = new ArrayList<Double>();
            List<Double> deltaMassesPPM = new ArrayList<Double>();
            for (Feature feature : featureMassMatchMap.keySet())
            {
                float featureMass = feature.getMz() * feature.getCharge();
                for (ChemicalFormula formula : featureMassMatchMap.get(feature).keySet())
                {
                    masses.add(formula.getCommonestIsotopeMass());
                    deltaMassesPPM.add((double)MassUtilities.convertDaToPPM(featureMass -
                            (float) formula.getCommonestIsotopeMass(), (float)formula.getCommonestIsotopeMass()));
                }
            }
            new PanelWithScatterPlot(masses, deltaMassesPPM,name ).displayInTab();
    }

    protected double[] calcMassCalibrationUsingMatches(Map<Feature, Map<ChemicalFormula, List<Adduct>>> featureMassMatchMap,
                                                       boolean byTime, boolean usePPM) throws IOException
    {
        List<Double> massDiffsListPreCal = new ArrayList<Double>();
        List<Double> featureXvalListPreCal = new ArrayList<Double>();

        String chartSuffix = byTime ? "_RT" : "_MASS";
        chartSuffix = chartSuffix + (usePPM ? "_PPM" : "DA");
        String xAxisLabel = byTime ? "time" : "mass";

        for (Feature feature : featureMassMatchMap.keySet())
        {
            float featureMass = feature.getMz() * feature.getCharge();

            for (ChemicalFormula formula : featureMassMatchMap.get(feature).keySet())
            {
                double deltaMassDa = featureMass - formula.getCommonestIsotopeMass();
                massDiffsListPreCal.add(usePPM ?
                        MassUtilities.convertDaToPPM((float)deltaMassDa, (float) formula.getCommonestIsotopeMass()) : deltaMassDa);
                featureXvalListPreCal.add(byTime ? (double) feature.getTime() : featureMass);
            }
        }
        double[] featureXvals = new double[featureXvalListPreCal.size()];
        double[] massDiffs = new double[massDiffsListPreCal.size()];
        for (int i=0; i<featureXvals.length; i++)
        {
            featureXvals[i] = featureXvalListPreCal.get(i);
            massDiffs[i] = massDiffsListPreCal.get(i);
        }

        int[] indexesLowStuff = RegressionUtilities.selectIndexesWithLowLeverageAndStudentizedResidual(
                featureXvals, massDiffs, calMaxLeverage,
                calMaxStudRes, false, 1, false, false);

        double[] featureXvalsForRegression = new double[indexesLowStuff.length];
        double[] massDiffsForRegression = new double[indexesLowStuff.length];
        for (int i=0; i<indexesLowStuff.length; i++)
        {
            featureXvalsForRegression[i] = featureXvals[indexesLowStuff[i]];
            massDiffsForRegression[i] = massDiffs[indexesLowStuff[i]];
        }
        double[] resultCoefficients = null;
        resultCoefficients = RegressionUtilities.modalRegression(featureXvalsForRegression,
                massDiffsForRegression, calRegressionDoF);

        if (showCharts)
        {
            double minXval = BasicStatistics.min(featureXvals);
            double maxXval = BasicStatistics.max(featureXvals);

            PanelWithScatterPlot psp = new PanelWithScatterPlot();
            String chartName = "MassCal" + chartSuffix;
            psp.setName(chartName);
            psp.addLineOrCurve(resultCoefficients, minXval, maxXval);

            psp.addData(featureXvalsForRegression, massDiffsForRegression,
                    "Matches used in regression");
            psp.addData(featureXvals, massDiffs, "all mass matches");

            String massDiffLabel = "mass difference" + (usePPM ? " (ppm)" : " (da)");
            psp.setAxisLabels(xAxisLabel,massDiffLabel);
            psp.displayInTab();
        }
        return resultCoefficients;
    }

    public Map<Feature, Map<ChemicalFormula, List<Adduct>>> massMatchBestOnly(Feature[] features,
                                                           Map<ChemicalFormula, List<Adduct>> formulaAdductsMap,
                                                           float massTolPPM)
    {
        Map<Feature, Map<ChemicalFormula, List<Adduct>>> result =
                 massMatchFull(features, formulaAdductsMap, massTolPPM, 1);
        for (Feature feature : result.keySet())
        {
            double smallestDeltaMassPPM = Double.MAX_VALUE;
            ChemicalFormula bestMatch = null;
            for (ChemicalFormula formula : result.get(feature).keySet())
            {
                double deltaMassPPM = MassUtilities.convertDaToPPM(
                        (float)(feature.mass - formula.getCommonestIsotopeMass()), feature.mass);
                if (Math.abs(deltaMassPPM) < smallestDeltaMassPPM)
                {
                    smallestDeltaMassPPM = deltaMassPPM;
                    bestMatch = formula;
                }
            }
            Map<ChemicalFormula, List<Adduct>> bestMap = new HashMap<ChemicalFormula, List<Adduct>>();
            bestMap.put(bestMatch, result.get(feature).get(bestMatch));
            result.put(feature, bestMap);
        }
        return result;
    }


    protected Map<ChemicalFormula, List<Adduct>> getOrCreateFormulaAdductsMap()
    {
        if (allFormulaAdductsMap == null)
        {
            System.err.println("Generating formula->adducts map");
            allFormulaAdductsMap = createFormulaAdductsMap(databaseCompounds);
//for (List<Adduct> adductList : allFormulaAdductsMap.values()) { for (Adduct adduct : adductList) System.err.println(
//        adduct.getCompoundNameAndIonTypeString() + ", mass=" + adduct.getCommonestIsotopeMass()); }            
        }
        return allFormulaAdductsMap;
    }


    /**
     * Apply all chemical modifications to all compounds to which they can be performed, and assemble
     * the results into a map from formulas to lists of adducts with that formula.
     *
     * After generating full map, clean it up.  For each chemical formula, if there are multiple
     * adducts with the same SMILES string, remove any that have modifications
     * @param compounds
     * @return
     */
    protected Map<ChemicalFormula, List<Adduct>> createFormulaAdductsMap(List<ChemicalCompound> compounds)
    {
        Map<String, ChemicalFormula> formulaStringFormulaMap = new HashMap<String, ChemicalFormula>();
        Map<ChemicalFormula, List<Adduct>> allFormulaAdductsMap = new HashMap<ChemicalFormula, List<Adduct>>();
        Map<ChemicalModification, List<ChemicalCompound>> modCompoundsMap =
                new HashMap<ChemicalModification, List<ChemicalCompound>>();
        _log.debug("Creating formula-adducts map for " + compounds.size() + " compounds.");
        if (chemicalModifications.isEmpty())
            _log.debug("\tNo modifications to compounds");
        else
        {
            _log.debug("\tModifications to compounds: ");
            for (ChemicalModification mod : chemicalModifications)
            {
                _log.debug("\t\t" + mod.getName());
            }
        }
        for (int i=0; i<compounds.size(); i++)
        {
            ChemicalCompound compound = compounds.get(i);
            if (i % (compounds.size() / 20) == 0)
                   _log.debug("Created adducts for " + i + " compounds...");

            List<Adduct> adductsThisCompound = new ArrayList<Adduct>();

            Adduct baseAdduct = new Adduct(compound);

            if (useUnmodifiedAdduct) {
                //add the unmodified mass
                adductsThisCompound.add(new Adduct(baseAdduct));
            }

            if (!chemicalModifications.isEmpty())
            {
                for (ChemicalModification mod : chemicalModifications)
                {
                    if (mod.canPerform(baseAdduct))
                    {
                        Adduct modAdduct = new Adduct(compound);
                        mod.perform(modAdduct);
                        adductsThisCompound.add(modAdduct);
//System.err.println(modAdduct.getIonTypeString());                        
//System.err.println(new SmilesGenerator().createSMILES(compound.getCdkMolecule()));                        

                        //accounting
                        List<ChemicalCompound> compoundsThisMod = modCompoundsMap.get(mod);
                        if (compoundsThisMod == null)
                        {
                            compoundsThisMod = new ArrayList<ChemicalCompound>();
                            modCompoundsMap.put(mod, compoundsThisMod);
                        }
                        compoundsThisMod.add(compound);
                    }
                }
            }

            for (Adduct adduct : adductsThisCompound)
            {
                ChemicalFormula formula = adduct.getFormula();
                if (formulaStringFormulaMap.containsKey(formula.toString()))
                        formula = formulaStringFormulaMap.get(formula.toString());
                else formulaStringFormulaMap.put(formula.toString(), formula);

//if (uniqueFormulaStrings.contains(formula.toString()))
//    System.err.println("DUPE!!!! " + formula.toString());
//uniqueFormulaStrings.add(formula.toString());

                List<Adduct> adductsThisFormula = allFormulaAdductsMap.get(formula);
                if (adductsThisFormula == null)
                {
                    adductsThisFormula = new ArrayList<Adduct>();
                    allFormulaAdductsMap.put(formula, adductsThisFormula);
                }
                adductsThisFormula.add(adduct);
            }
        }

        //todo: this is expensive but useful.  Need to put this back in if possible
        if (!chemicalModifications.isEmpty())
        {
            _log.debug("******NOT CHECKING FOR MODS WITH DUPE STRUCTURE!!!! TOO EXPENSIVE!! PUT BACK IN IF POSSIBLE****");

//            _log.debug("Checking for duplicate SMILES strings...");
//            //Remove modified duplicate adducts
//            for (List<Adduct> adductsForAFormula : allFormulaAdductsMap.values())
//                removeModifiedWhenDuplicateSMILES(adductsForAFormula);
        }

        if (!chemicalModifications.isEmpty())
        {
            ApplicationContext.infoMessage("Modifications and numbers of compounds applied to (out of " +
                    compounds.size() + "):");
            for (ChemicalModification mod : chemicalModifications)
            {
                int numComps = 0;
                if (modCompoundsMap.containsKey(mod))
                    numComps = modCompoundsMap.get(mod).size();
                ApplicationContext.infoMessage(mod.getName() + ": " + numComps);
            }
        }
        _log.debug("Created adducts for all compounds");
        return allFormulaAdductsMap;
    }

    /**
     * Remove any modified Adducts for which an unmodified Adduct exists with the same SMILES string
     * @param adducts
     */
    public void removeModifiedWhenDuplicateSMILES(List<Adduct> adducts)
    {
        List<String> unmodifiedSMILESStrings = new ArrayList<String>();
        List<Adduct> modifiedAdducts = new ArrayList<Adduct>();
        for (Adduct adduct : adducts)
        {
            if (adduct.getModifications() == null || adduct.getModifications().isEmpty())
                 unmodifiedSMILESStrings.add(ChemCalcs.createSMILESString(adduct.getMolecule()));
            else
                modifiedAdducts.add(adduct);
        }
        for (Adduct adduct : modifiedAdducts)
        {
            if (unmodifiedSMILESStrings.contains(ChemCalcs.createSMILESString(adduct.getMolecule())))
                adducts.remove(adduct);
        }
    }

    /**
     * This one uses the default formula-adducts map
     * @param featuresSortedMassAsc
     * @param massTolPPM
     * @param numPeaksMustMatch
     * @return
     */
    public Map<Feature, Map<ChemicalFormula, List<Adduct>>> massMatchFull(Feature[] featuresSortedMassAsc,
                                                              float massTolPPM, int numPeaksMustMatch)
    {
        return massMatchFull(featuresSortedMassAsc, getOrCreateFormulaAdductsMap(), massTolPPM, numPeaksMustMatch);
    }



    /**
     * This version takes a map from formulas to adducts
     * @param featuresSortedMassAsc
     * @param massTolPPM
     * @param numPeaksMustMatch at least numPeaksMustMatch peaks ALL must match within tolerance
     * @return
     */
    public Map<Feature, Map<ChemicalFormula, List<Adduct>>> massMatchFull(Feature[] featuresSortedMassAsc,
                                                              Map<ChemicalFormula, List<Adduct>> formulaAdductsMap,
                                                              float massTolPPM, int numPeaksMustMatch)
    {
        //make absolutely sure feature masses are fixed.
        for (Feature feature : featuresSortedMassAsc)
            feature.setMass(feature.getMz() * feature.getCharge());
        List<ChemicalFormula> formulasByMassAsc = new ArrayList<ChemicalFormula>(allFormulaAdductsMap.keySet());
        Collections.sort(formulasByMassAsc, new ChemicalFormula.ComparatorMassAsc());
        _log.debug("Unique chemical formulas in database, with adducts: " + formulasByMassAsc.size());
        if (_log.isDebugEnabled())
        {
            List<Integer> adductCountsByFormula = new ArrayList<Integer>();
            for (List<Adduct> formulas : formulaAdductsMap.values())
                adductCountsByFormula.add(formulas.size());
            _log.debug("Mean adducts per chemical formula: " + BasicStatistics.mean(adductCountsByFormula));
//            if (showCharts)
//                new PanelWithHistogram(adductCountsByFormula, "# of Adducts Per Formula").displayInTab();
        }

        int minPossibleMassIndex = 0;
        Map<Feature, Map<ChemicalFormula, List<Adduct>>> result =
                new HashMap<Feature, Map<ChemicalFormula, List<Adduct>>>();
        boolean newFeature = true;
        for (Feature feature :featuresSortedMassAsc)
        {
            //not using calculated feature mass, because that presumes M+H and subtracts the H
            float featureMass = feature.getMz() * feature.getCharge();
            float massToleranceDa = MassUtilities.calculateAbsoluteDeltaMass(featureMass, massTolPPM,
                            FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
            float minMatchMass = featureMass - massToleranceDa;
            float maxMatchMass = featureMass + massToleranceDa;

//            while ((formulasByMassAsc.get(minPossibleMassIndex).getCommonestIsotopeMass() + 1.0078125) < minMatchMass &&
            while (formulasByMassAsc.get(minPossibleMassIndex).getCommonestIsotopeMass() < minMatchMass &&
                   minPossibleMassIndex < formulasByMassAsc.size()-1)
            {
                minPossibleMassIndex++;
            }

            List<ChemicalFormula> formulasMatchedThisFeature = new ArrayList<ChemicalFormula>();
            for (int i=minPossibleMassIndex ; i< formulasByMassAsc.size(); i++)
            {
                ChemicalFormula formula = formulasByMassAsc.get(i);

                float compoundMass = (float) formula.getCommonestIsotopeMass();

                if (compoundMass < minMatchMass)
                    continue;
                else if (maxMatchMass < compoundMass)
                    break;
                else
                {
                    if (newFeature)
                    {
                        newFeature = false;
                    }
                    if (numPeaksMustMatch > 1)
                    {
                        float minPeaksFeatureOrCompound = Math.min(feature.comprised.length,
                                formula.getPeakMasses().length);
                        if (minPeaksFeatureOrCompound < numPeaksMustMatch)
                            continue;
                        int numMatched = 1;
                        for (int peakInd = 1; peakInd < minPeaksFeatureOrCompound; peakInd++)
                        {
                            float peakMass = Feature.convertMzToMass(feature.comprised[peakInd].mz, feature.charge);
                            float deltaMassPPM =
                                    MassUtilities.convertDaToPPM((float) formula.getPeakMasses()[peakInd] - peakMass,
                                                                 peakMass);
                            if (Math.abs(deltaMassPPM) > massTolPPM)
                            {
                                numMatched++;
                            }
                        }
                        if (numMatched >= numPeaksMustMatch)
                            formulasMatchedThisFeature.add(formula);
                    }
                    else
                        formulasMatchedThisFeature.add(formula);
                }
            }


            if (!formulasMatchedThisFeature.isEmpty())
            {
                Map<ChemicalFormula, List<Adduct>> resultThisFeature = new HashMap<ChemicalFormula, List<Adduct>>();
                for (ChemicalFormula formula : formulasMatchedThisFeature)
                    resultThisFeature.put(formula, formulaAdductsMap.get(formula));
                result.put(feature, resultThisFeature);
            }
        }
        return result;
    }

    /**
     *
     * @param featureSet
     * @param looseMassTolerancePPM
     * @return
     * @throws IOException
     */
    public double calcRunMassError(FeatureSet featureSet, float looseMassTolerancePPM) throws IOException
    {
        List<ChemicalCompound> compoundsList = new ArrayList<ChemicalCompound>();
        List<ChemicalCompound> compoundsWith2PeaksList = new ArrayList<ChemicalCompound>();

//        List<Float> dbAllMassPeakDistances = new ArrayList<Float>();
//        List<Float> dbAllMassPeakDistancesPPMDist = new ArrayList<Float>();
//        List<Float> dbAllMasses = new ArrayList<Float>();
        for (ChemicalCompound compound : databaseCompounds)
        {
            try
            {
                compoundsList.add(compound);
                if (compound.getPeakMasses().length == 1)
                {
                    _log.debug("Single-peak compound " + compound);
                }
                else
                {
                    float peakDist = (float)(compound.getPeakMasses()[1] - compound.getPeakMasses()[0]);
                    if (Math.abs(peakDist) < 1.5f)
                    {
//                        dbAllMasses.add((float) compound.getMass());
//                        dbAllMassPeakDistances.add(peakDist);
//                        if (Math.abs(MassUtilities.convertDaToPPM(peakDist - 1.0035f, (float) compound.getMass())) < 50)
//                            dbAllMassPeakDistancesPPMDist.add(MassUtilities.convertDaToPPM(peakDist - 1.0035f, (float) compound.getMass()));
                        compoundsWith2PeaksList.add(compound);
                    }
                    else
                    {
                        _log.debug("Weird second-peak mass: " + compound);
                    }
                }
            }
            catch (IllegalArgumentException e)
            {
                ApplicationContext.setMessage("Skipping bad compound: " + e.getMessage());
            }

        }

//        if (showCharts)
//        {
//            new PanelWithHistogram(dbAllMassPeakDistances, "db peak distances").displayInTab();
//            new PanelWithHistogram(dbAllMassPeakDistancesPPMDist, "db peak spread ppm").displayInTab();
//            new PanelWithScatterPlot(dbAllMasses, dbAllMassPeakDistances, "db mass vs peakdist").displayInTab();
//        }


        Collections.sort(compoundsWith2PeaksList, new ChemicalCompound.ComparatorMassAsc());
        Collections.sort(compoundsList, new ChemicalCompound.ComparatorMassAsc());


        Feature[] features = featureSet.getFeatures();
        List<Feature> featuresWithMultiplePeaksList = new ArrayList<Feature>();
        for (Feature feature : features)
            if (feature.peaks > 1)
                 featuresWithMultiplePeaksList.add(feature);
        Feature[] featuresWithMultiplePeaks = featuresWithMultiplePeaksList.toArray(new Feature[featuresWithMultiplePeaksList.size()]);

        ApplicationContext.infoMessage("Processing file " + featureSet.getSourceFile().getName() + " with " + features.length +
                " features");
        Arrays.sort(featuresWithMultiplePeaks, new Feature.MassAscComparator());

        //disabling, because this is useless
        if (showCharts && false)
        {
            List<Float> featureMasses = new ArrayList<Float>();
            List<Float> featurePeakDistances = new ArrayList<Float>();
            for (Feature feature : featuresWithMultiplePeaks)
            {
                float featureMass = feature.getMz() * feature.getCharge();
                featureMasses.add(featureMass);
                featurePeakDistances.add(feature.comprised[1].mz - feature.comprised[0].mz);
            }
            new PanelWithScatterPlot(featureMasses, featurePeakDistances, "feature_mass_peakdist").displayInTab();
        }


        Map<Feature, Map<ChemicalFormula, List<Adduct>>> featureMassMatchMap =
                massMatchBestOnly(featuresWithMultiplePeaks, createFormulaAdductsMap(compoundsWith2PeaksList),
                        looseMassTolerancePPM);
        List<Float> deltaMassPPMList = new ArrayList<Float>();
        List<Float> deltaMassPPMSecondPeakList = new ArrayList<Float>();
        List<Float> compoundMassesFor2PeakDiffList = new ArrayList<Float>();
        List<Float> peakDeltaMassDaDifferencesList = new ArrayList<Float>();
        List<Float> peakDeltaMassPPMDifferencesList = new ArrayList<Float>();
        List<Float> secondPeakDeltaMassPPMList = new ArrayList<Float>();


        for (Feature feature : featureMassMatchMap.keySet())
        {
            float featureMass = feature.getMz() * feature.getCharge();
            Map<ChemicalFormula, List<Adduct>> thisFeatureMatches = featureMassMatchMap.get(feature);
            for (ChemicalFormula formula : thisFeatureMatches.keySet())
            {
                double compoundMass = formula.getCommonestIsotopeMass();
                float massCompoundSecondPeak = (float) formula.getPeakMasses()[1];

//if (MassUtilities.convertDaToPPM(Math.abs((float) mass - trackedMass), trackedMass) < 1)
// System.err.println("MATCHED feature mass " + feature.getMass() + " to compound " + compound);
                float deltaMassPPM = MassUtilities.convertDaToPPM(featureMass - (float) compoundMass, (float) compoundMass);
                deltaMassPPMList.add(deltaMassPPM);

                float featureSecondPeakMass = feature.comprised[1].mz / feature.charge;
//System.err.println(mass + ", " + massCompoundSecondPeak + ", " + featureSecondPeakMass);
                float deltaMassSecondPeakPPM =
                        MassUtilities.convertDaToPPM(featureSecondPeakMass - massCompoundSecondPeak,
                                massCompoundSecondPeak);
                deltaMassPPMSecondPeakList.add(deltaMassSecondPeakPPM);
//System.err.println("****" + deltaMassPPM + ", " + deltaMassSecondPeakPPM);
                //todo: 30 is arbitrary
                if (Math.abs(deltaMassSecondPeakPPM - deltaMassPPM) < 30)  {
                    compoundMassesFor2PeakDiffList.add((float) compoundMass);
                    peakDeltaMassPPMDifferencesList.add(deltaMassSecondPeakPPM - deltaMassPPM);
                    peakDeltaMassDaDifferencesList.add(featureSecondPeakMass - featureMass);
                    secondPeakDeltaMassPPMList.add(deltaMassSecondPeakPPM);
                }
            }
        }
        System.err.println("Median mass difference between first and second peak: " + BasicStatistics.median(peakDeltaMassDaDifferencesList));


        if (peakDeltaMassPPMDifferencesList.size() < 10)
            throw new IOException("Not enough 2-peak matches (" + 
                    peakDeltaMassPPMDifferencesList.size() + ") to analyze distribution");

        if (showCharts)
        {
            new PanelWithHistogram(secondPeakDeltaMassPPMList, "peak2_deltappm").displayInTab();            
            new PanelWithHistogram(peakDeltaMassPPMDifferencesList, "2peak_diff").displayInTab();
            new PanelWithScatterPlot(compoundMassesFor2PeakDiffList, peakDeltaMassPPMDifferencesList, "mass_vs_2peak_diff",
                    "Compound Mass", "2-peak delta").displayInTab();            
        }

        ApplicationContext.setMessage("Calculating mixed distribution...");
        Collections.sort(deltaMassPPMList);
        Pair<float[], float[]> paramsAndProbs =
                calculateDistParamsAndProbsEM(peakDeltaMassPPMDifferencesList, DEFAULT_INITIAL_MATCH_PROPORTION);
        float[] distParams = paramsAndProbs.first;
        float[] probabilities = paramsAndProbs.second;

        //todo: implicitly assuming mu==0.  Should be true after calibration.  But I don't think this is actually true
        float mu = distParams[0];
        float sigma = distParams[1];
        float proportion = distParams[2];
        float area = distParams[3];

        float ppmContainingMostPoints = sigmaMultipleForErrorDist * sigma;
        ApplicationContext.setMessage("Done calculating mixed distribution.  sigma=" + sigma +
                ".  Most points within " + ppmContainingMostPoints + "ppm");
        ApplicationContext.setMessage("Mu: " + mu);
        float massErrorPPM = (float) (ppmContainingMostPoints / Math.sqrt(2));
        ApplicationContext.setMessage("Run mass error: " + massErrorPPM);
//        float parity = 0.5f;
//
//        int densParityIndex = Math.abs(Arrays.binarySearch(probabilities, parity));
//        float ppmAtParity = deltaMassPPMList.get(densParityIndex);
//        _log.debug("Optimizing ppm tolerance, staring with " + ppmAtParity + "ppm");
//        double closestProbability = Float.MAX_VALUE;
//        float increment = 0.2f;
//        while (ppmAtParity > 0)
//        {
//            double proportionTimesNormdist = proportion * BasicStatistics.calcNormalCumDensity(mu, sigma, ppmAtParity);
//            double probability = proportionTimesNormdist / (proportionTimesNormdist + ((1-proportion) / area) );
//System.err.println("PPM: " + ppmAtParity + ", prob: " + probability);
//            if (densParityIndex == 0 || probability >= closestProbability)
//                break;
//            ppmAtParity -= increment;
//        }


//        float[] fdrs = AmtMatchProbabilityAssigner.calculateFDRsForProbabilities(probabilities);
//        new PanelWithHistogram(fdrs, "FDRs").displayInTab();

        if (showCharts)
        {
            List<Float> deltaMassesKeptList = new ArrayList<Float>();
            List<Float> deltaMassPPMSecondPeakKeptList = new ArrayList<Float>();
            List<Float> deltaMassesSecondUnderList = new ArrayList<Float>();
            List<Float> deltaMassPPMSecondPeakSecondUnderList = new ArrayList<Float>();
            for (int i=0; i<deltaMassPPMList.size(); i++)
            {
                if (Math.abs(deltaMassPPMSecondPeakList.get(i)) < massErrorPPM)
                {
                    deltaMassesSecondUnderList.add(deltaMassPPMList.get(i));
                    deltaMassPPMSecondPeakSecondUnderList.add(deltaMassPPMSecondPeakList.get(i));
                    if (Math.abs(deltaMassPPMList.get(i)) < massErrorPPM)
                    {
                        deltaMassesKeptList.add(deltaMassPPMList.get(i));
                        deltaMassPPMSecondPeakKeptList.add(deltaMassPPMSecondPeakList.get(i));                                  
                    }
                }
            }
            PanelWithScatterPlot pwspTwoPeakError =
                    new PanelWithScatterPlot(deltaMassesKeptList, deltaMassPPMSecondPeakKeptList, "2peakerror");
            pwspTwoPeakError.addData(deltaMassesSecondUnderList, deltaMassPPMSecondPeakSecondUnderList, "peak 2 under");
            pwspTwoPeakError.addData(deltaMassPPMList, deltaMassPPMSecondPeakList, "all points");
            pwspTwoPeakError.setAxisLabels("Peak 1 Delta Mass", "Peak 2 Delta Mass");
            pwspTwoPeakError.displayInTab();
            PanelWithHistogram pwhDiffs = new PanelWithHistogram(peakDeltaMassPPMDifferencesList, "peak_dist_diff");
            pwhDiffs.displayInTab();
            float[] peakDeltaMassPPMDifferencesArray = new float[peakDeltaMassPPMDifferencesList.size()];
            for (int i=0; i<peakDeltaMassPPMDifferencesArray.length; i++)
                peakDeltaMassPPMDifferencesArray[i] = peakDeltaMassPPMDifferencesList.get(i);
            PanelWithScatterPlot deltaMassProbPlot = new PanelWithScatterPlot(peakDeltaMassPPMDifferencesArray,
                    probabilities, "PeakDeltaMassVsProb");
            deltaMassProbPlot.setAxisLabels("Peak Pair DeltaMass (PPM)","Probability Matching Same Ion");
            deltaMassProbPlot.displayInTab();
        }

        return massErrorPPM;
    }




    /**
      * Use the Expectation Maximization algorithm to determine the parameters of the target normal distribution
     * and the proportion of datapoints in that distribution vs. the uniform background.
      * @param targetMassErrorDataList
      * @param proportionTrue
      * @return a float[] containing mu, sigma, and proportion in the normal dist
      * @throws IOException
      */
     public Pair<float[],float[]> calculateDistParamsAndProbsEM(List<Float> targetMassErrorDataList,
                                          float proportionTrue)
             throws IOException
     {
         //uniquify masses, so duplicate masses don't mess up distribution
         Map<Integer, List<Integer>> uniqueOrigIndexesMap = new HashMap<Integer, List<Integer>>();
         List<Float> uniqueMassesList = new ArrayList<Float>();

         for (int i=0; i<targetMassErrorDataList.size(); i++)
         {
             Float mass = targetMassErrorDataList.get(i);
             int uniqueIndex = uniqueMassesList.indexOf(mass);
             if (uniqueIndex == -1)
             {
                 uniqueMassesList.add(mass);
                 uniqueIndex = uniqueMassesList.size()-1;
             }

             List<Integer> indexes = uniqueOrigIndexesMap.get(uniqueIndex);
             if (indexes == null)
             {
                 indexes = new ArrayList<Integer>();
                 uniqueOrigIndexesMap.put(uniqueIndex, indexes);
             }
             indexes.add(i);
         }

         if (uniqueMassesList.size() < 10)
            throw new IllegalArgumentException("calculateDistParamsAndProbsEM: Too few mass matches (" +
                    uniqueMassesList.size() + ") to calculate parameters");

         double[] uniqueMasses = new double[uniqueOrigIndexesMap.size()];
         for (int i=0; i<uniqueMasses.length; i++)
            uniqueMasses[i] = uniqueMassesList.get(i);


         int numPoints = uniqueMasses.length;



         Map<String,Object> rScalarVarMap = new HashMap<String,Object>();
         Map<String,double[]> rVectorVarMap = new HashMap<String,double[]>();

         rVectorVarMap.put("targetx",uniqueMasses);

         _log.debug("calculateDistParamsAndProbsEM, targetx: " + uniqueMasses.length);


         rScalarVarMap.put("proportion",(double) proportionTrue);
         double area = BasicStatistics.max(uniqueMasses) - BasicStatistics.min(uniqueMasses);
         rScalarVarMap.put("area",area);
         rScalarVarMap.put("miniterations",(double) 30);
         rScalarVarMap.put("maxiterations",(double) 80);
         rScalarVarMap.put("max_deltap_proportion_for_stable",(double) .005f);
         rScalarVarMap.put("max_deltap_proportion",(double) .005f + 0.001);

         rScalarVarMap.put("iters_stable_for_converg",(double) 1);

         rScalarVarMap.put("showcharts",showCharts ? "TRUE" : "FALSE");
         File outChartFile = null;
         File normTestOutChartFile = null;

         if (showCharts)
         {
             outChartFile = TempFileManager.createTempFile("em_plots.jpg", this);
             normTestOutChartFile = TempFileManager.createTempFile("em_normtest_plots.jpg", this);

             rScalarVarMap.put("chart_out_file", "'" + RInterface.generateRFriendlyPath(outChartFile) + "'");
             rScalarVarMap.put("normtest_chart_out_file", "'" +
                     RInterface.generateRFriendlyPath(normTestOutChartFile) + "'");

         }

         String calculateProbabilitiesCommand =
           RInterface.readResourceFile("/org/fhcrc/cpl/viewer/metabologna/assign_massmatch_probabilities_em.R");

         //this timeout is arbitrary and dangerous.  Woo hoo!
         int timeoutMilliseconds=900000;
         String rResult = RInterface.evaluateRExpression(calculateProbabilitiesCommand,
                 rScalarVarMap, rVectorVarMap, null, null, timeoutMilliseconds);

         Map<String, String> varResultStringMap =
                 RInterface.extractVariableStringsFromListOutput(rResult);
         float[] errorProbabilities = new float[numPoints];

         int probsIndex = 0;
         String[] probChunks = varResultStringMap.get("probs").split("\\s");
         for (String probChunk : probChunks)
         {
             if (!probChunk.contains("[") && probChunk.length() > 0)
             {
                 errorProbabilities[probsIndex++] = Float.parseFloat(probChunk);
             }
         }
         if (probsIndex != numPoints)
             throw new IOException("FAILED to read probabilities correctly back from R!");

         boolean converged = varResultStringMap.get("converged").contains("TRUE");
         String numIterString = varResultStringMap.get("num_iterations");
         numIterString = numIterString.substring(numIterString.indexOf("]") + 1).trim();
         int num_iterations = Integer.parseInt(numIterString);
         if (converged)
         {

             ApplicationContext.setMessage("EM estimation converged after " +
                     num_iterations + " iterations");
         }
         else
         {
             ApplicationContext.infoMessage("WARNING!!! EM estimation failed to converge after " +
                     num_iterations + " iterations");
         }

         String[] ksChunks = varResultStringMap.get("ksresults").split("\\s");
         float ks_score_x = Float.parseFloat(ksChunks[2]);

         //todo: move parsing of vector results into RInterface
         String[] distParamValues = varResultStringMap.get("dist_params").split("\\s");
         List<Float> paramVals = new ArrayList<Float>();
         for (String paramVal : distParamValues)
             if (paramVal != null && paramVal.length() > 1 && !paramVal.contains("["))
                 paramVals.add(Float.parseFloat(paramVal));

         float mu_x = paramVals.get(0);
         float sigma_x = paramVals.get(1);
         float proportion = paramVals.get(2);


         _log.debug("Distribution params: mu_=" + mu_x+
                    ", sigma=" + sigma_x + ", proportion=" + proportion);


         if (ks_score_x < 0.005f )
         {
             _log.debug("WARNING!!!!  Kolmogorov-Smirnov test indicates that these matching data " +
                 "may not conform to a bivariate normal distribution intermingled with a uniform " +
                 "distribution.  If this key assumption fails, match probabilities will not be " +
                 "accurate.  Please consider re-running this analysis in non-parametric mode. " +
                 "KS p-value = " + ks_score_x);
         }
         else
         {
             ApplicationContext.setMessage("KS normality test passed.  KS value = " + ks_score_x );
         }

         if (showCharts)
         {
             try
             {
                 PanelWithBlindImageChart pwbic =
                         new PanelWithBlindImageChart(outChartFile,"EM Parameters");
                 pwbic.displayInTab();

             }
             catch (Exception e)
             {
                 ApplicationContext.errorMessage("Error displaying error cutoff chart images, file: "
                               + outChartFile.getAbsolutePath(),e);
             }
             try
             {
                 PanelWithBlindImageChart pwbic =
                         new PanelWithBlindImageChart(normTestOutChartFile,"EM dist analysis");
                 pwbic.displayInTab();
             }
             catch (Exception e)
             {
                 ApplicationContext.infoMessage("Error displaying error cutoff chart images, file: "
                               + normTestOutChartFile.getAbsolutePath() + ", error: " + e.getMessage());
             }
         }
         TempFileManager.deleteTempFiles(this);

         float[] distParams = new float[] { mu_x, sigma_x, proportion, (float) area };

         float[] matchProbabilities = new float[targetMassErrorDataList.size()];
         for (int i=0; i<errorProbabilities.length; i++)
         {
             List<Integer> indexes = uniqueOrigIndexesMap.get(i);
             for (int matchIndex : indexes)
             {
                matchProbabilities[matchIndex] = errorProbabilities[i];
             }
         }

         return new Pair<float[], float[]>(distParams, matchProbabilities);
     }



    public float getCalMassTolerancePPM()
    {
        return calMassTolerancePPM;
    }

    public void setCalMassTolerancePPM(float calMassTolerancePPM)
    {
        this.calMassTolerancePPM = calMassTolerancePPM;
    }

    public float getCalMaxLeverage()
    {
        return calMaxLeverage;
    }

    public void setCalMaxLeverage(float calMaxLeverage)
    {
        this.calMaxLeverage = calMaxLeverage;
    }

    public double getCalMaxStudRes()
    {
        return calMaxStudRes;
    }

    public void setCalMaxStudRes(double calMaxStudRes)
    {
        this.calMaxStudRes = calMaxStudRes;
    }

    public List<ChemicalCompound> getDatabaseCompounds()
    {
        return databaseCompounds;
    }

    public void setDatabaseCompounds(List<ChemicalCompound> databaseCompoundsByMass)
    {
        this.databaseCompounds = databaseCompoundsByMass;
    }

    public boolean isShowCharts()
    {
        return showCharts;
    }

    public void setShowCharts(boolean showCharts)
    {
        this.showCharts = showCharts;
    }

    public int getCalRegressionDoF() {
        return calRegressionDoF;
    }

    public void setCalRegressionDoF(int calRegressionDoF) {
        this.calRegressionDoF = calRegressionDoF;
    }

    public List<ChemicalModification> getChemicalModifications() {
        return chemicalModifications;
    }

    public void setChemicalModifications(List<ChemicalModification> chemicalModifications) {
        this.chemicalModifications = chemicalModifications;
    }

    public boolean isUseUnmodifiedAdduct() {
        return useUnmodifiedAdduct;
    }

    public void setUseUnmodifiedAdduct(boolean useUnmodifiedAdduct) {
        this.useUnmodifiedAdduct = useUnmodifiedAdduct;
    }
}
