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
package org.fhcrc.cpl.viewer.metabologna.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.Window2DFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.viewer.metabologna.*;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.apache.log4j.Logger;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.exception.CDKException;

import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;


/**
 * Match metabolites to a metabolite database.  Note that this requires changing the feature masses of both MS1
 * and MS2 features to be mz * charge, rather than (mz - H) * charge
 */
public class MetaboliteDBMatcherCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MetaboliteDBMatcherCLM.class);

    protected File[] featureFiles;
    protected File massesFile;
    protected float matchMassTolerancePPM = 10f;

    protected int regressionDegreesOfFreedom = MetaboliteDatabaseMatcher.DEFAULT_CALIBRATION_REGRESSION_DEGREES_FREEDOM;

    protected boolean showCharts = false;

    protected FeatureSet ms2FeatureSet = null;
    protected File ms2FeaturesDir = null;

    protected File outFile = null;         
    protected File outFeatureFile = null;

    protected File outDir = null;
    protected File outFeatureDir = null;

    protected int currentFileIndex = 0;

    //for testing.  A single mass to track and report details on
    protected float trackedMass = 165.078979f;

    protected boolean shouldCalibrateMasses = false;
    protected boolean shouldWriteAllFeatures = false;

    protected double maxStudRes = 1;


    protected List<ChemicalCompound> databaseCompounds = null;

    protected MetaboliteDatabaseMatcher metaboliteDBMatcher = new MetaboliteDatabaseMatcher();

    protected int minPeaksMustMatch = 1;

    protected boolean shouldCalcMassError = false;
    protected double massErrorPPM = 2;

    protected boolean shouldGuessFormulas = false;

    protected boolean shouldAdd2H = false;
    protected boolean shouldAddWater = false;
    protected boolean shouldAddH = true;

    protected float looseMatchMassTolerancePPM = 7;




    public MetaboliteDBMatcherCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matchmetabolitedb";

        mHelpMessage ="matchmetabolitedb";
        mShortDescription = "matchmetabolitedb";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedSeriesFileArgumentDefinition(true, "feature file(s)"),
                       new FileToReadArgumentDefinition("db", true, "File full of masses to check.  Can have other columns"),
                       new DecimalArgumentDefinition("masstoleranceppm", false, "mass tolerance for matching, ppm", matchMassTolerancePPM),
                       new BooleanArgumentDefinition("showcharts", false,
                               "show charts?", showCharts),
                       new FileToReadArgumentDefinition("ms2scanfile", false, "MS2 scan feature file"),
                       new DirectoryToReadArgumentDefinition("ms2dir", false, "MS2 scan feature file directory"),

                       new FileToWriteArgumentDefinition("out", false, "Output spreadsheet File"),
                       new FileToWriteArgumentDefinition("outfeatures", false, "Output Feature File"),

                       new DirectoryToWriteArgumentDefinition("outdir", false, "Output Directory for mass spreadsheets"),

                       new DirectoryToWriteArgumentDefinition("outfeaturedir", false, "Output Directory for feature " +
                               "files containing all (or only matched) features"),
                       new BooleanArgumentDefinition("writeallfeatures", false, "Should all features be written, or only matches?",
                               shouldWriteAllFeatures),
                       new BooleanArgumentDefinition("calibratemasses", false,
                               "Calibrate feature masses based on loose (20ppm) mass match?", shouldCalibrateMasses),
                       new IntegerArgumentDefinition("regressiondof", false,
                               "Degrees of freedom to use for regression (mass calibration)", regressionDegreesOfFreedom),
                       new DecimalArgumentDefinition("maxstudres", false,
                               "maximum studentized residual for calibration regression", maxStudRes),
                       new DecimalArgumentDefinition("calmasstoleranceppm", false, "mass tolerance for calibration, ppm",
                               MetaboliteDatabaseMatcher.DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM),
                       new DecimalArgumentDefinition("calmaxleverage", false, "maximum leverage numerator for calibration",
                               MetaboliteDatabaseMatcher.DEFAULT_CALIBRATION_MAX_LEVERAGE),
                       new IntegerArgumentDefinition("minpeaksmatch", false,
                               "Minimum number of peaks that must match within mass tolerance", minPeaksMustMatch),
                       new DecimalArgumentDefinition("masserrorppm", false,
                               "Mass error for run (in PPM. Specifying this turns off mass error calculation)", massErrorPPM),
                       new StringListArgumentDefinition("addformulas", false,
                               "Formula snippets that can appear as additions to formulas in the database"),
                       new StringListArgumentDefinition("subtractformulas", false,
                               "Formula snippets that can be subtractions from formulas in the database"),
                       new BooleanArgumentDefinition("guessformulas", false,
                               "Populate all possible chemical formulas for each feature?", shouldGuessFormulas),
                       new BooleanArgumentDefinition("add2h", false,
                               "Should we try to convert C=C to C-C and add H2?", shouldAdd2H),
                       new BooleanArgumentDefinition("addh", false,
                               "Should we add an H in an undefined location?", shouldAddH),
                       new BooleanArgumentDefinition("addwater", false,
                               "Should we try to convert C=C to C-C and add H2O?", shouldAddWater),
                       new BooleanArgumentDefinition("usebaseadduct", false,
                               "Should we include the [M] adduct (no additions)?", false),
                       new DecimalArgumentDefinition("loosemasstolppm", false,
                               "Loose matching mass tolerance (ppm)", looseMatchMassTolerancePPM),
                       new BooleanArgumentDefinition("calcmasserror", false,
                               "Should we calculate mass error specifically for this run?", shouldCalcMassError),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        regressionDegreesOfFreedom = getIntegerArgumentValue("regressiondof");
        metaboliteDBMatcher.setCalRegressionDoF(regressionDegreesOfFreedom);
        metaboliteDBMatcher.setUseUnmodifiedAdduct(getBooleanArgumentValue("usebaseadduct"));
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        matchMassTolerancePPM =  getFloatArgumentValue("masstoleranceppm");
        metaboliteDBMatcher.setCalMassTolerancePPM(getFloatArgumentValue("calmasstoleranceppm"));
        metaboliteDBMatcher.setCalMaxLeverage(getFloatArgumentValue("calmaxleverage"));

        shouldCalcMassError = getBooleanArgumentValue("calcmasserror");
        massErrorPPM = getFloatArgumentValue("masserrorppm");

        if (shouldCalcMassError && hasArgumentValue("masserrorppm"))
            throw new ArgumentValidationException("Mass error supplied, but also indicated mass error should be calculated.");

        shouldGuessFormulas = getBooleanArgumentValue("guessformulas");
        shouldAdd2H = getBooleanArgumentValue("add2h");
        shouldAddH = getBooleanArgumentValue("addh");

        shouldAddWater = getBooleanArgumentValue("addwater");

        looseMatchMassTolerancePPM = getFloatArgumentValue("loosemasstolppm");


        List<String> formulaStringsToAdd = getStringListArgumentValue("addformulas");
        List<ChemicalModification> chemicalModsToAdd = new ArrayList<ChemicalModification>();
        if (formulaStringsToAdd != null && !formulaStringsToAdd.isEmpty())
            for (String formulaString : formulaStringsToAdd)
            {
                ChemicalFormula formula = new ChemicalFormula(formulaString);

                ApplicationContext.infoMessage("Formula to add: " + formula + ", mass=" +
                        formula.getCommonestIsotopeMass());
                chemicalModsToAdd.add(new SimpleChemicalAddition(formula));
            }

        if (shouldAdd2H)
                chemicalModsToAdd.add(new ReduceDoubleBondAdd2HMod());
        if (shouldAddH) {
            System.err.println("Using M+H mod");
                chemicalModsToAdd.add(new UnknownHAdditionMod());
        }
        if (shouldAddWater)
                chemicalModsToAdd.add(new ReduceDoubleBondAddWaterMod());

//for (ChemicalModification mod : chemicalModsToAdd)
//         System.err.println("***" + mod.getName());

//        List<String> formulaStringsToSubtract = getStringListArgumentValue("subtractformulas");
//        if (formulaStringsToSubtract != null && !formulaStringsToSubtract.isEmpty())
//            for (String formulaString : formulaStringsToSubtract)
//            {
//                ChemicalFormula formula = new ChemicalFormula(formulaString);
//
//                ApplicationContext.infoMessage("Formula to subtract: " + formula + ", mass=" +
//                        formula.getCommonestIsotopeMass());
//                chemicalFormulasToAdd.add(new SimpleChemicalSubtraction(formula));
//            }

        metaboliteDBMatcher.setChemicalModifications(chemicalModsToAdd);

        minPeaksMustMatch = getIntegerArgumentValue("minpeaksmatch");
        try
        {
            //populate 2 theoretical peaks for each compound
            databaseCompounds =
                    ChemicalCompound.loadCompoundsFromFile(getFileArgumentValue("db"), 2);
            ApplicationContext.infoMessage("Loaded " + databaseCompounds.size() + " compounds from database");
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            throw new ArgumentValidationException("Failed to load compound database " +
                    getFileArgumentValue("db").getAbsolutePath() + ". Type: " +
                    e.getClass().getName() + ", message: " + e.getMessage(),e);
        }
//for (ChemicalCompound compound : databaseCompoundsByMass) System.err.println(compound.getCommonestIsotopeMass());
        metaboliteDBMatcher.setDatabaseCompounds(databaseCompounds);

        showCharts = getBooleanArgumentValue("showcharts");
        metaboliteDBMatcher.setShowCharts(showCharts);
        shouldWriteAllFeatures = getBooleanArgumentValue("writeallfeatures");

        outFile = getFileArgumentValue("out");
        outFeatureFile = getFileArgumentValue("outfeatures");

        outDir = getFileArgumentValue("outdir");
        outFeatureDir = getFileArgumentValue("outfeaturedir");

        shouldCalibrateMasses = getBooleanArgumentValue("calibratemasses");

        maxStudRes = getDoubleArgumentValue("maxstudres");

        if (featureFiles.length == 0 && (outFile != null || outDir == null))
        {
            throw new ArgumentValidationException("If more than one feature file specified, must specify output dir");
        }

        if (hasArgumentValue("ms2scanfile"))
        {
            try
            {
                ms2FeatureSet = new FeatureSet(getFileArgumentValue("ms2scanfile"));
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Failed to load MS2 scan features",e);
            }
        }

        ms2FeaturesDir = getFileArgumentValue("ms2dir");

        if (featureFiles.length > 1)
            assertArgumentAbsent("ms2scanfile");
        if (ms2FeaturesDir != null)
            assertArgumentAbsent("ms2scanfile","ms2dir");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        for (int i=0; i<featureFiles.length; i++)
        {
            File featureFile = featureFiles[i];
            ApplicationContext.infoMessage("Processing file " + featureFile.getAbsolutePath());
            FeatureSet ms2FeatureSetThisFile = ms2FeatureSet;
            if (ms2FeaturesDir != null)
            {
                try
                {
                    File ms2FeatureFile =
                            CommandLineModuleUtilities.findFileLikeFile(featureFiles[i], ms2FeaturesDir, "ms2.tsv");
                    ApplicationContext.infoMessage("Found MS2 file " + ms2FeatureFile.getAbsolutePath());

                    ms2FeatureSetThisFile = new FeatureSet(ms2FeatureFile);
                }
                catch (FileNotFoundException e)
                {
                    throw new CommandLineModuleExecutionException("No MS2 file found for file " +
                            featureFiles[i].getAbsolutePath());
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Failed to open MS2 file associated with file " +
                            featureFiles[i].getAbsolutePath());
                }
            }
            FeatureSet featureSet = null;
            try
            {
                featureSet = new FeatureSet(featureFile);
//for playing with FDR                
//for (Feature feature : featureSet.getFeatures()) feature.setMass(feature.getMass() + (float)(5*1.00035));                
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to open feature file " +
                        featureFile.getAbsolutePath());
            }

            if (outFile != null)
                processFile(featureSet, outFile, outFeatureFile, ms2FeatureSetThisFile);
            else
            {

                File outputFile = null;
                File outputFeatureFile = null;
                if (outDir != null)
                    outputFile = CommandLineModuleUtilities.createOutputFile(featureFile, "matched.tsv", outDir);
                if (outFeatureDir != null)
                    outputFeatureFile = CommandLineModuleUtilities.createOutputFile(featureFile, "matched.features.tsv",
                            outFeatureDir);
                processFile(featureSet, outputFile, outputFeatureFile ,ms2FeatureSetThisFile);
            }
            currentFileIndex++;
            if (i < featureFiles.length-1)
                MultiChartDisplayPanel.destroySingleton();
        }
    }

    /**
     * do the actual work
     */
    public void processFile(FeatureSet featureSet, File outputFile, File outputFeatureFile, FeatureSet ms2FeatureSetThisFile)
            throws CommandLineModuleExecutionException
    {
        Arrays.sort(featureSet.getFeatures(), new Feature.MassAscComparator());

        double[] massCalCoefficients = null;
        if (shouldCalibrateMasses)
        {
            try
            {
                massCalCoefficients = metaboliteDBMatcher.calibrateFeatureMasses(featureSet.getFeatures());
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failure during mass calibration", e);
            }
        }

        if (showCharts)
        {
            Map<Feature, Map<ChemicalFormula, List<Adduct>>> looseMatch =
                    metaboliteDBMatcher.massMatchFull(featureSet.getFeatures(), looseMatchMassTolerancePPM, 1);
            List<Float> finalMatchedMassDiffsPPM = new ArrayList<Float>();
            for (Feature feature : looseMatch.keySet())
            {
                for (ChemicalFormula formula : looseMatch.get(feature).keySet())
                {
                     finalMatchedMassDiffsPPM.add(MassUtilities.convertDaToPPM(
                             feature.getMass() - (float) formula.getCommonestIsotopeMass(),
                             (float) formula.getCommonestIsotopeMass()));
                }
            }
            new PanelWithHistogram(finalMatchedMassDiffsPPM, "All Mass Errors (ppm)").displayInTab();
        }

        if (shouldCalcMassError)
        {
            try
            {
                massErrorPPM = metaboliteDBMatcher.calcRunMassError(featureSet, looseMatchMassTolerancePPM);
            }
            catch (IOException e)
            {
//                throw new CommandLineModuleExecutionException("Failed to calculate mass error in run " + featureSet.getSourceFile() + ": " + e.getMessage());
                massErrorPPM = looseMatchMassTolerancePPM;
                System.err.println("Failed to calculate run mass error!!!!  Using " +
                        looseMatchMassTolerancePPM + ".  Error message: " + e.getMessage());
            }
        }

        String headerLine =
                "id\tFeatureMass\tIntensity\tCharge\tScan\tPeaks\tDeltaPPM\tCompound\tAdductMass\tCompoundFormula\t" +
                "AdductFormula\tIonType\tSMILES\tCompoundClass";
        FeatureSetMatcher.FeatureMatchingResult ms2MatchingResult = null;
        if (ms2FeatureSetThisFile != null)
        {

            if (shouldCalibrateMasses)
            {
                for (Feature feature : ms2FeatureSetThisFile.getFeatures())
                {
                    float massDiffDa = (float) RegressionUtilities.mapValueUsingCoefficients(massCalCoefficients,
                            feature.getMass());
                    feature.setMass(feature.getMass() - massDiffDa);
                    feature.updateMz();
                }
            }

            //important: to mass-match MS1 features, ms2 feature masses must be changed to mz * charge
            for (Feature ms2Feature : ms2FeatureSetThisFile.getFeatures())
                ms2Feature.setMass(ms2Feature.getMz() * ms2Feature.getCharge());

            headerLine = headerLine + "\tMS2Scans";
            Window2DFeatureSetMatcher fsm = new Window2DFeatureSetMatcher();

            //todo: be more scientific about this.  Just tripling the MS1 mass tolerance. elution diff is baseless
            float ms2MassTolerancePPM = 3 * (float) massErrorPPM;
            fsm.setMaxMassDiff(ms2MassTolerancePPM);
            fsm.setMinMassDiff(-ms2MassTolerancePPM);
            fsm.setMassDiffType( FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
            fsm.setMaxElutionDiff(30);
            fsm.setMinElutionDiff(-30);
            fsm.setElutionMode(Window2DFeatureSetMatcher.ELUTION_MODE_TIME);
            fsm.setMatchWithinChargeOnly(true);

            ms2MatchingResult = fsm.matchFeatures(featureSet, ms2FeatureSetThisFile);
            ApplicationContext.infoMessage("MS2 matches: " + ms2MatchingResult.size());
        }

        if (shouldGuessFormulas)
        {
            ApplicationContext.infoMessage("Guessing chemical formulas for " + featureSet.getFeatures().length + " features...");
            List<Integer> formulasPerFeature = new ArrayList<Integer>();
            List<Float> featureMasses = new ArrayList<Float>();

            int i=0;
            for (Feature feature : featureSet.getFeatures())
            {
                i++;
                if (i % 50 != 0)
                    continue;
                if (feature.mass > 300 && i % 200 != 0)
                    continue;
//                if (feature.mass > 400 && i % 500 != 0)
//                    continue;
//                if (feature.mass > 500 && i % 1000 != 0)
//                    continue;

                if (feature.mass > 450)
                    continue;
                if (i % (featureSet.getFeatures().length / 40) == 0)
                    ApplicationContext.infoMessage(100 * (float) i / featureSet.getFeatures().length + "% complete.... (" + i + " features)");
                try
                {
                    List<ChemicalFormula> formulas = ChemCalcs.CDKFormulaSet2ChemFormList(ChemCalcs.calcMass2Formulas(feature.mass, massErrorPPM));
                    MetaboliteExtraInfoDef.setFormulaList(feature,formulas);
//Set<String> formulaStrings = new HashSet<String>();
//                    for (ChemicalFormula formula : formulas)
//                    {
//                        if (formulaStrings.contains(formula.toString()))
//                            System.err.println("DUPE!!!! " + formula.toString());
//                        formulaStrings.add(formula.toString());
//                        double deltaMass = MassUtilities.convertDaToPPM((float)formula.getCommonestIsotopeMass() - feature.mass, feature.mass);
//                        if (Math.abs(deltaMass) > 2)
//                            System.err.println("\tMASS! " + deltaMass);
//                    }
                    int numFormulas = MetaboliteExtraInfoDef.getFormulaList(feature).size();
                    formulasPerFeature.add(numFormulas);

                    featureMasses.add(feature.mass);
                }
                catch (Exception e)
                {
                    _log.debug("Formula restriction failure");
                }
            }
            ApplicationContext.infoMessage("Done guessing chemical formulas");
            if (showCharts)
            {
                new PanelWithHistogram(formulasPerFeature, "Formulas per feature").displayInTab();
                new PanelWithScatterPlot(featureMasses, formulasPerFeature, "Mass vs. #formulas").displayInTab();

            }
        }

        PrintWriter outPW = new PrintWriter(System.out);
        if (outputFile != null)
        {
            try
            {
                outPW = new PrintWriter(outputFile);
            }
            catch (FileNotFoundException e)
            {
                throw new CommandLineModuleExecutionException("Failed to open output file",e);
            }
        }

        outPW.println(headerLine);

        Map<Feature, Map<ChemicalFormula, List<Adduct>>> featureMassMatchMap =
                metaboliteDBMatcher.massMatchFull(featureSet.getFeatures(),
                        (float) massErrorPPM, minPeaksMustMatch);
        List<Float> finalMatchedMassDiffsPPM = new ArrayList<Float>();
        int featureId = 0;
        List<Feature> outMatchedFeatures = new ArrayList<Feature>();
        for (Feature feature : featureMassMatchMap.keySet())
        {
            boolean newFeature = true;
            Map<ChemicalFormula, List<Adduct>> thisFeatureMatches = featureMassMatchMap.get(feature);
            for (ChemicalFormula formula : thisFeatureMatches.keySet())
            {
                float mass = (float) formula.getCommonestIsotopeMass();

                float deltaMassPPM = MassUtilities.convertDaToPPM(feature.getMass() - mass, mass);
                finalMatchedMassDiffsPPM.add(deltaMassPPM);
                if (newFeature)
                {
                    featureId++;
                    newFeature = false;
                }

                for (Adduct adduct : thisFeatureMatches.get(formula))
                {

                    String line = featureId + "\t" + feature.getMass() + "\t" + feature.getIntensity() + "\t" +
                            feature.getCharge()+ "\t" + feature.getScan() + "\t" +
                            feature.peaks + "\t" + deltaMassPPM + "\t" + adduct.getCompound().getName() + "\t" +
                            adduct.getCommonestIsotopeMass() + "\t" + adduct.getCompound().getFormula() + "\t" +
                            adduct.getFormula() + "\t" + adduct.getIonTypeString() + "\t" +
                            ChemCalcs.createSMILESString(adduct.getMolecule()) + "\t" +
                            adduct.getCompound().getCompoundClass();
                    if (ms2FeatureSetThisFile != null)
                    {
                        line = line + "\t";
                        List<Feature> matchedMS2Features = ms2MatchingResult.get(feature);
                        if (matchedMS2Features != null)
                        {
                            for (int j=0; j<matchedMS2Features.size(); j++)
                            {
                                if (j>0)
                                    line = line + ";";
                                line = line + matchedMS2Features.get(j).getScan();
                            }
                        }
                    }
                    outPW.println(line);
                    outPW.flush();
                }
            }

            List<Adduct> equivalentMassesWithNames = null;
            for (List<Adduct> adducts : thisFeatureMatches.values())
            {
                for (Adduct adduct : adducts)
                {
                    if (equivalentMassesWithNames == null)
                    {
                        equivalentMassesWithNames = new ArrayList<Adduct>();
                        equivalentMassesWithNames.add(adduct);
                        continue;
                    }
                    float massDiff = Math.abs((float) adduct.getCommonestIsotopeMass() - feature.getMass());
                    float oldMassDiff =  Math.abs((float) equivalentMassesWithNames.get(0).getCommonestIsotopeMass() - feature.getMass());
                    if (massDiff < oldMassDiff)
                    {
                        equivalentMassesWithNames = new ArrayList<Adduct>();
                        equivalentMassesWithNames.add(adduct);
                    }
                    else if (massDiff == oldMassDiff)
                    {
                        equivalentMassesWithNames.add(adduct);
                    }
                }
            }

            String massName = "";
            Collections.sort(equivalentMassesWithNames, new Adduct.ComparatorNameAndIonTypeAsc());
            for (int i=0; i<equivalentMassesWithNames.size(); i++)
            {
                if (i>0)
                {
//System.err.println("***" + massName);                        
                    massName = massName + "|";
                }
                massName = massName + equivalentMassesWithNames.get(i).getCompound().getName();
            }
            //feature.setDescription(massName);
            MS2ExtraInfoDef.addPeptide(feature, massName);
            outMatchedFeatures.add(feature);

        }

        ApplicationContext.setMessage("Matched features within " + massErrorPPM + "ppm: " + outMatchedFeatures.size());


        if (showCharts)
        {
//            new PanelWithScatterPlot(featureTimesList, deltaMassPPMList, "RTvsError").displayInTab();

            new PanelWithHistogram(finalMatchedMassDiffsPPM, "Retained Mass Errors (ppm)").displayInTab();
//            new PanelWithScatterPlot(matchMassList, deltaMassPPMList, dmvmName).displayInTab();

            //build matched compound class piechart if possible
            Map<String, Float> pieChartData = new HashMap<String, Float>();
            int totalClassesFound = 0;
            for (Map<ChemicalFormula, List<Adduct>> formulaMatchesMap: featureMassMatchMap.values())
            {
                Set<String> classesThisMatch = new HashSet<String>();
                for (List<Adduct> adducts: formulaMatchesMap.values())
                {
                    for (Adduct adduct : adducts)
                    {
                        if (adduct.getCompound().getCompoundClass() != null)
                            classesThisMatch.add(adduct.getCompound().getCompoundClass());
                    }
                }
                for (String classThisMatch : classesThisMatch)
                {
                    pieChartData.put(classThisMatch,
                            pieChartData.containsKey(classThisMatch) ?
                                 pieChartData.get(classThisMatch) + 1 :
                                    1);
                    totalClassesFound++;
                }
            }

            if (!pieChartData.isEmpty())
            {
                    PanelWithPieChart pwpc = new PanelWithPieChart(pieChartData, "Classes Found");
                    pwpc.displayInTab();
                    ApplicationContext.infoMessage("Classes found:");
                    for (String classFound : pieChartData.keySet())
                    {
                        ApplicationContext.infoMessage("\t" + classFound + ": " + pieChartData.get(classFound).intValue());
                    }
                    ApplicationContext.infoMessage("Total: " + totalClassesFound);
            }
            else
            {
                ApplicationContext.infoMessage("No compound class data found, can't build piechart");
            }
        }

        if (outputFile != null)
            outPW.close();
        if (outputFeatureFile != null)
        {
            featureSet.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());
            if (!shouldWriteAllFeatures)
                featureSet.setFeatures(outMatchedFeatures.toArray(new Feature[outMatchedFeatures.size()]));

            try
            {
            featureSet.save(outputFeatureFile, false, "APML");
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to save matched features to file "  +
                        outputFeatureFile.getAbsolutePath());
            }
            ApplicationContext.infoMessage("Saved matched feature file " +
                    outputFeatureFile.getAbsolutePath() + " with " +
                    featureSet.getFeatures().length + " features");
        }
    }

    public IMolecule createCOOH()
    {
        IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
        IMolecule molecule = builder.newMolecule();
        molecule.addAtom(builder.newAtom("C"));
        molecule.addAtom(builder.newAtom("O"));
        molecule.addBond(0, 1, IBond.Order.DOUBLE);
        molecule.addAtom(builder.newAtom("O"));
        molecule.addBond(0, 2, IBond.Order.SINGLE);

        try {

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
            LonePairElectronChecker lpcheck = new LonePairElectronChecker();
            lpcheck.saturate(molecule);

        }
        catch (CDKException e)
        {
            e.printStackTrace();
        }
//        Image image = MoleculeRenderer2D.renderMolecule(molecule, width, height);
//        new PanelWithBlindImageChart((BufferedImage) image, "0=C-OH").displayInTab();
        return molecule;
    }
}






