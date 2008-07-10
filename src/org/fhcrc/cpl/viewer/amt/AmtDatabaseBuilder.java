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
package org.fhcrc.cpl.viewer.amt;

import java.util.*;
import java.io.File;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.toolbox.RegressionUtilities;
import org.fhcrc.cpl.viewer.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;

/**
 * Builds AMT databases, in various ways
 *
 * If you want to build an AMT database in some other way, put a method in here to do it
 * your way.
 */
public class AmtDatabaseBuilder
{
    static Logger _log = Logger.getLogger(AmtDatabaseBuilder.class);

    //default cutoff values for tossing out observations with high studentized residuals
    //Marty suggested 2.3 for this
    public static final double DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_REGRESSION = 2.0;
    //Marty suggested 3.0 for this
    public static final double DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_INCLUSION = 2.0;

    //default parameters for median alignment
    public static final int DEFAULT_MIN_PEPTIDES_FOR_ALIGNMENT_REGRESSION=8;
    public static final int DEFAULT_MIN_OBSERVATIONS_FOR_ALIGNMENT_REGRESSION=3;

    //Default maximum leverage numerator (X/n) to keep for regression
    public static final double DEFAULT_MAX_LEVERAGE_NUMERATOR = 4;


    /**
     * Create the default N-Terminal modifications that are used by X!Tandem.
     * In later X!Tandem versions, these will be declared explicitly.
     * If they're not, need to add them.
     * @return
     */
    public static MS2Modification[] createDefaultNTerminalModifications()
    {
        MS2Modification eMod = new MS2Modification();
        eMod.setAminoAcid("E");
        eMod.setMassDiff(-18.0106f);
        eMod.setVariable(true);

        MS2Modification qMod = new MS2Modification();
        qMod.setAminoAcid("Q");
        qMod.setMassDiff(-17.0265f);
        qMod.setVariable(true);

        return new MS2Modification[]
                {
                        eMod,
                        qMod
                };
    }

    /**
     * calculate observed hydrophobicity for all features and
     * add them all as observations to a new database.
     *
     * This one's the workhorse.  Everything else calls it eventually
     *
     * @param scanOrTimeMode
     * @return
     */
    public static AmtDatabase createAmtDatabaseForRun(FeatureSet ms2FeatureSet,
                                                      int scanOrTimeMode,
                                                      double[] timeToHCoefficients,
                                                      boolean calculateHydrophobicities,
                                                      Map<String, Integer> spectralCountsMap)
    {
        //clone the featureset and strip out peptides that fall below minimum pprophet
        FeatureSet featureSetClone = (FeatureSet) ms2FeatureSet.clone();
        FeatureSet.FeatureSelector fs = new FeatureSet.FeatureSelector();
        featureSetClone = featureSetClone.filter(fs);

        //build a map of known offsets from the features
        Feature[] features = featureSetClone.getFeatures();

        if (calculateHydrophobicities)
            AmtUtilities.addHydrophobicityToFeatures(
                        featureSetClone.getFeatures(), scanOrTimeMode,
                        timeToHCoefficients);

        //Now that we've got hydrophobicity on the features, create a db

        AmtDatabase result = new AmtDatabase();
        AmtRunEntry runEntry = new AmtRunEntry(timeToHCoefficients,
                                               null);

        MS2Modification[] runModifications = MS2ExtraInfoDef.getFeatureSetModifications(ms2FeatureSet);


        //check for explicitly-declared E and Q modifications.  If they're not declared, add them
        boolean hasEorQMod = false;
        if (runModifications != null)
        {
            for (MS2Modification mod : runModifications)
            {
                String aa = mod.getAminoAcid();
                if ("E".equals(aa) || "Q".equals(aa))
                {
                    hasEorQMod=true;
                    break;
                }
            }
        }
        if (!hasEorQMod)
        {
            _log.debug("Adding default E and Q modifications, not specified explicitly in file");

            int existingModsLength=0;
            if (runModifications != null)
                existingModsLength = runModifications.length;
            MS2Modification[] defaultEQMods = createDefaultNTerminalModifications();
            MS2Modification[] newModifications =
                    new MS2Modification[existingModsLength + defaultEQMods.length];
            if (runModifications != null)
                System.arraycopy(runModifications, 0, newModifications, 0, runModifications.length);
            for (int i=0; i<defaultEQMods.length; i++)
                newModifications[existingModsLength+i] = defaultEQMods[i];
            runModifications = newModifications;
        }

        runEntry.setModifications(runModifications);
        if (_log.isDebugEnabled())
        {
            if (MS2ExtraInfoDef.getFeatureSetModifications(ms2FeatureSet) != null)
                for (MS2Modification mod : runModifications)
                    _log.debug("Adding modification from run: " + mod);
        }
        File sourceFile = ms2FeatureSet.getSourceFile();
        if (sourceFile != null)
            runEntry.setPepXmlFilename(ms2FeatureSet.getSourceFile().getName());
        result.addRunEntry(runEntry);

        for (Feature feature : features)
            result.resolveModsAndAddObservation(
                    MS2ExtraInfoDef.getFirstPeptide(feature),
                    MS2ExtraInfoDef.getModifiedAminoAcids(feature),
                    MS2ExtraInfoDef.getPeptideProphet(feature),
                    AmtExtraInfoDef.getObservedHydrophobicity(feature), runEntry,
                    spectralCountsMap,
                    feature.getTime());
        return result;
    }

    protected static Map<String, Double> performInitialRegression(Feature[] allFeatures,
                                                                  int scanOrTimeMode)
    {
        Feature[] featuresForFirstRegression =
                RegressionUtilities.selectFeaturesWithLowLeverage(allFeatures,
                        DEFAULT_MAX_LEVERAGE_NUMERATOR, scanOrTimeMode);

        //Do an initial regression
        return AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                        featuresForFirstRegression,
                        scanOrTimeMode,
                        false);
    }

    /**
     * Calculate "studentized residuals" for each feature.
     *
     * This is based on the residual of the feature, and on its leverage. It's a measure
     * of the certainty that the feature is messed up w.r.t. the rest of the features.
     *
     * @param allFeatures
     * the residual squared must be in order to be dropped from the set
     * @return
     */
    public static double[] calculateStudentizedResiduals(Feature[] allFeatures,
                                                         int scanOrTimeMode)
    {
        //number of features
        int n = allFeatures.length;

        if (scanOrTimeMode == RegressionUtilities.REGRESSION_MODE_TIME)
        {
            boolean nonzeroTimesExist = false;
            for (Feature feature : allFeatures)
            {
                if (feature.getTime() > 0)
                {
                    nonzeroTimesExist = true;
                    break;
                }
            }
            if (!nonzeroTimesExist)
            {
                throw new IllegalArgumentException("ERROR! Feature retention times are all zero!  " +
                        "Please populate feature " +
                        "times (e.g., using the 'populatems2times' command), or provide associated mzXML file(s) to " +
                        "use for populating scan times.");
            }
        }

        //Do an initial regression using all features
        Map<String, Double> scanOrTimeToHydroLine =
                performInitialRegression(allFeatures, scanOrTimeMode);
        //residuals based on the initial regression (y - yhat)
        double[] hydrophobicityResiduals = new double[n];
        //scan/time values for all features (x)
        double[] scanOrTimeValues = new double[n];

        //gather the information we need about each feature
        for (int i=0; i<n; i++)
        {
            //record all x values
            scanOrTimeValues[i] = ((scanOrTimeMode == RegressionUtilities.REGRESSION_MODE_TIME) ?
                            allFeatures[i].getTime() :
                            allFeatures[i].getScan());
            //this value was already implicitly calculated above, when regression was
            //performed, but we don't have easy access to it.  Anyway, fast.
            //This is yhat
            double predictedH =
                    RegressionUtilities.predictYFromX(
                            AmtUtilities.getSlopeFromRegressionLine(scanOrTimeToHydroLine),
                            AmtUtilities.getInterceptFromRegressionLine(scanOrTimeToHydroLine),
                            scanOrTimeValues[i]);
            //calculated H: y
            double calculatedH =
                             AmtUtilities.calculateNormalizedHydrophobicity(
                                     MS2ExtraInfoDef.getFirstPeptide(allFeatures[i]));

            //residual is y - yhat
            hydrophobicityResiduals[i] = calculatedH - predictedH;
        }

        //sigmaHydroError requires a denomenator of n-2, so can't use
        //BasicStatistics.standardDeviation
        double sumSquaredResiduals=0;
        for (double residual : hydrophobicityResiduals)
            sumSquaredResiduals += Math.pow(residual, 2);
        double sigmaHydroError = Math.sqrt(sumSquaredResiduals / (n - 2));

        //calculate the leverages of all values
        double[] leverages = BasicStatistics.leverages(scanOrTimeValues);

        double[] studentizedResiduals = new double[n];

        //1 + 1/n.  Same for all inputs, compute once
        double onePlusOneOverN = 1 + 1/n;

        //Find the "studentized residual" of every feature:
        for (int i=0; i < n; i++)
        {
            double scaledResidual = hydrophobicityResiduals[i] / sigmaHydroError;
            studentizedResiduals[i] = scaledResidual /
                    Math.sqrt(onePlusOneOverN + leverages[i]);
        }
        _log.debug("studentized residual mean: " + BasicStatistics.mean(studentizedResiduals) + ", stddev: " + BasicStatistics.standardDeviation(studentizedResiduals));

        return studentizedResiduals;
    }

    /**
     * Given an array of features, choose just the ones with residual <=
     * of maxStudentizedResidual.  This is done twice, so it's a method.
     * @param allFeatures
     * @param studentizedResiduals
     * @param maxStudentizedResidual
     * @return
     */
    public static Feature[]
            chooseFeaturesWithMaxStudentizedResidual(Feature[] allFeatures,
                                                     double[] studentizedResiduals,
                                                     double maxStudentizedResidual)
    {
        //square the max studentized residual, compare with square of studentized
        //residual (since studentized residual is signed)
        double maxStudentizedResidualSquared = Math.pow(maxStudentizedResidual,2);
        List<Feature> passingFeatures = new ArrayList<Feature>();
        for (int i=0; i<allFeatures.length; i++)
        {
            if (Math.pow(studentizedResiduals[i],2) <= maxStudentizedResidualSquared)
                passingFeatures.add(allFeatures[i]);
        }
        return passingFeatures.toArray(new Feature[0]);
    }

    /**
     * Count the spectra with each peptide identification.  Assumes filtering
     * has already been done such that these features are the ones we consider good
     * @param features
     * @return
     */
    protected static Map<String,Integer> countSpectraForPeptides(Feature[] features)
    {
        Map<String, Integer> result = new HashMap<String, Integer>();

        for (Feature feature : features)
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            if (result.containsKey(peptide))
                result.put(peptide, result.get(peptide) + 1);
            else
                result.put(peptide, 1);
        }

        return result;
    }

    /**
     * Create an AMT database for a single run, doing those arcane things we do
     * to cast out the unworthy features.  Here's the sequence:
     * 1.  Calculate the leverage of all features
     * 2.  Toss out all features with leverage > 4/n, perform a regression to
     * map time to normalized H (first cut)
     * 3.  Calculate the studentized residual of all features from the original set
     * 4.  Take out features with studentized residual < 2.0, perform a regression
     * to map time to normalized H (final version)
     * 5.  Take features with studentized residual < 2.0, insert those peptides
     * into the database 
     *
     *
     *
     * TODO: parameterize residual cutoffs
     * @param ms2FeatureSet this method assumes features have times populated
     * @param scanOrTimeMode
     * @param robustRegression
     * @return
     */
    public static AmtDatabase createAmtDatabaseForRun(FeatureSet ms2FeatureSet,
                                                      int scanOrTimeMode,
                                                      boolean robustRegression,
                                                      Pair<Integer, Integer> outNumFeaturesChosen,
                                                      double maxStudentizedResidualForRegression,
                                                      double maxStudentizedResidualForInclusion)
    {
        Map<String, Integer> peptideSpectralCountsMap =
                countSpectraForPeptides(ms2FeatureSet.getFeatures());

        //one feature per modification state
        MS2ExtraInfoDef.removeAllButFirstFeatureForEachPeptide(ms2FeatureSet,
                                                               true);

        Feature[] allMs2Features = ms2FeatureSet.getFeatures();
        double[] studentizedResiduals =
                calculateStudentizedResiduals(allMs2Features,
                                              scanOrTimeMode);

        Feature[] featuresForRegression =
                chooseFeaturesWithMaxStudentizedResidual(allMs2Features,
                        studentizedResiduals,
                        maxStudentizedResidualForRegression);

        Map<String, Double> timeHydroLine =
                AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                        featuresForRegression,
                        scanOrTimeMode,
                        robustRegression);
        double[] timeToHCoefficients = new double[2];
        timeToHCoefficients[1] = AmtUtilities.getSlopeFromRegressionLine(timeHydroLine);
        timeToHCoefficients[0] = AmtUtilities.getInterceptFromRegressionLine(timeHydroLine);


        if (AmtUtilities.getSlopeFromRegressionLine(timeHydroLine) < 0)
                ApplicationContext.infoMessage("WARNING: negative regression line slope!");

        //the features for inclusion in the database may be chosen based on a different
        //cutoff from those chosen for regression.
        //TODO: Should we be recalculating the studentized residuals based on the prediction line
        //we calculated?
        Feature[] featuresForInclusion =
                chooseFeaturesWithMaxStudentizedResidual(allMs2Features,
                        studentizedResiduals,
                        maxStudentizedResidualForInclusion);

        outNumFeaturesChosen.first = featuresForRegression.length;
        outNumFeaturesChosen.second = featuresForInclusion.length;
        _log.debug("Chose " + featuresForRegression.length + "/" + allMs2Features.length +
                " features for regression");
        _log.debug("Chose " + featuresForInclusion.length + "/" + allMs2Features.length +
                " features for inclusion");

        FeatureSet featureSetForInclusion = new FeatureSet(featuresForInclusion);
        featureSetForInclusion.setSourceFile(ms2FeatureSet.getSourceFile());
        MS2Modification[] modifications = MS2ExtraInfoDef.getFeatureSetModifications(ms2FeatureSet);
        _log.debug("Setting modifications: " + modifications);
        MS2ExtraInfoDef.setFeatureSetModifications(featureSetForInclusion,
                modifications);

        AmtDatabase amtDatabase =
                createAmtDatabaseForRun(featureSetForInclusion, scanOrTimeMode,
                        timeToHCoefficients,
                        true, peptideSpectralCountsMap);
        if (_log.isDebugEnabled())
        {
            int numTotalObservations = 0;
            for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
                numTotalObservations += peptideEntry.getNumObservations();
            _log.debug("Number of total observations: " + numTotalObservations);
        }
        return amtDatabase;
    }


    /**
     * Given a directory full of pepXml files and a directory full of the
     * corresponding mzXML files, create an AMT database from all the pepXml files,
     * using the mzXML files to grab retention times.
     *
     * Highly dependent on naming conventions.  If you're starting with an all.pep.xml
     * file, use the method below instead.
     *
     * The nice thing about this is that it can check for existence of all the
     * mzXml files before starting work.
     * @param pepXmlDir
     * @param mzXmlDir
     * @param scanOrTimeMode
     * @return
     * @throws Exception
     */
    public static AmtDatabase createAmtDBFromDirectories(File pepXmlDir, File mzXmlDir,
                                                         int scanOrTimeMode,
                                                         FeatureSet.FeatureSelector featureSelector,
                                                         boolean robustRegression,
                                                         double maxStudentizedResidualForRegression,
                                                         double maxStudentizedResidualForInclusion,
                                                         boolean align)
            throws Exception
    {
        boolean getTimesFromMzxml = (mzXmlDir != null);
        
        AmtDatabase resultDB = new AmtDatabase();
        List<File> ms2FeatureFiles = new ArrayList<File>();
        Map<File,File> pepXmlMzXmlFileMap = new HashMap<File,File>();
        String[] potentialFiles = pepXmlDir.list();
        Arrays.sort(potentialFiles);
        for (String pepXmlFileName : potentialFiles)
        {
            File ms2FeatureFile =
                    new File(pepXmlDir.getAbsolutePath() + File.separatorChar +
                             pepXmlFileName);
            if (ms2FeatureFile.exists() && !ms2FeatureFile.isDirectory())
                ms2FeatureFiles.add(ms2FeatureFile);
        }
        if (getTimesFromMzxml)
        {
            for (File pepXmlFile : ms2FeatureFiles)
            {
                File mzXmlFile = CommandLineModuleUtilities.findFileLikeFile(pepXmlFile, mzXmlDir,"mzXML");

//                String pepXmlFileName = pepXmlFile.getName();
//                String mzXmlFileName =
//                        (pepXmlFileName.substring(0, pepXmlFileName.indexOf(".")) + ".mzXML");
//                File mzXmlFile = null;
//
//                for (String potentialMzXmlFilename : mzXmlDir.list())
//                {
//
//                    if (potentialMzXmlFilename.equalsIgnoreCase(mzXmlFileName))
//                    {
//                        mzXmlFile = new File(mzXmlDir.getAbsolutePath() + File.separatorChar +
//                                potentialMzXmlFilename);
//                        break;
//                    }
//                }
//                if (mzXmlFile == null)
//                    throw new FileNotFoundException("Could not find an mzXML file that's a match for pepXml file " + pepXmlFileName);

                pepXmlMzXmlFileMap.put(pepXmlFile, mzXmlFile);
            }
        }

        //for each pair of pepXml and mzXml files, create an amt database
        //and add its observations, runs and mods to the cumulative db
        int numFeaturesChosenForRegression = 0;
        int numFeaturesChosenForInclusion = 0;

        Pair<Integer,Integer> numFeaturesDummy = new Pair<Integer,Integer>(0,0);

        for (File ms2FeatureFile : ms2FeatureFiles)
        {
            ApplicationContext.infoMessage("Adding features from file " + ms2FeatureFile.getAbsolutePath());

            FeatureSet thisFeatureSet = new FeatureSet(ms2FeatureFile);
            thisFeatureSet = thisFeatureSet.filter(featureSelector);
            if (getTimesFromMzxml)
            {
                File mzXmlFile = pepXmlMzXmlFileMap.get(ms2FeatureFile);
                _log.debug("Populating times for MS2 file " + ms2FeatureFile + " from mzXML file " + mzXmlFile);
                MSRun thisRun = MSRun.load(mzXmlFile.getAbsolutePath());
                thisFeatureSet.populateTimesForMS2Features(thisRun);
            }

            if (align)
            {
                addRunToAmtDatabase(resultDB, thisFeatureSet,
                        scanOrTimeMode, robustRegression, maxStudentizedResidualForRegression,
                        maxStudentizedResidualForInclusion, 2, .05, 10);
            }
            else
            {
                AmtDatabase thisRunDB =
                        createAmtDatabaseForRun(thisFeatureSet,
                                scanOrTimeMode, robustRegression,
                                numFeaturesDummy, maxStudentizedResidualForRegression,
                                maxStudentizedResidualForInclusion);
                numFeaturesChosenForRegression += numFeaturesDummy.first;
                numFeaturesChosenForInclusion += numFeaturesDummy.second;

                resultDB.addObservationsFromAnotherDatabase(thisRunDB);
                _log.debug("Features chosen for regression: " + numFeaturesChosenForRegression);
                _log.debug("Features chosen for inclusion: " + numFeaturesChosenForInclusion);
            }
        }

        return resultDB;
    }

    /**
     * workinprogress
     *
     * Add a run to an AMT database.
     *
     * After linear H calculation, if we have significant agreeing peptide overlap (at least
     * minMatchedPeptides that match database entries with at least minDBObservationsForAlignmentMatch
     * observations, within maxHDiffForAlignmentMatch H units) then use those agreeing
     * peptides as points for alignment, and set all the rest of the H values according
     * to that alignment. 
     *
     * @param amtDatabase
     * @param newRunFeatureSet
     * @param scanOrTimeMode
     * @param robustRegression
     * @param maxStudentizedResidualForRegression
     * @param maxStudentizedResidualForInclusion
     * @param maxHDiffForAlignmentMatch
     * @param minMatchedPeptidesForAlignment
     */
    public static void addRunToAmtDatabase(AmtDatabase amtDatabase, FeatureSet newRunFeatureSet,
                                           int scanOrTimeMode,
                                           boolean robustRegression,
                                           double maxStudentizedResidualForRegression,
                                           double maxStudentizedResidualForInclusion,
                                           double minDBObservationsForAlignmentMatch,
                                           double maxHDiffForAlignmentMatch,
                                           int minMatchedPeptidesForAlignment)
    {
        Pair<Integer, Integer> outNumFeaturesChosen = new Pair<Integer, Integer>(0,0);
        AmtDatabase thisRunDB = createAmtDatabaseForRun(
                newRunFeatureSet,
                scanOrTimeMode,
                robustRegression,
                outNumFeaturesChosen,
                maxStudentizedResidualForRegression,
                maxStudentizedResidualForInclusion);
        amtDatabase.addObservationsFromAnotherDatabase(thisRunDB);
    }

    /**
     * Create an AMT database from an all.pep.xml file -- pass null
     * to the workhorse method so that it knows we need to look in the file
     * from a reference to the mzXML file
     * @param allPepXmlFile
     * @param scanOrTimeMode
     * @param featureSelector
     * @param robustRegression
     * @return
     * @throws Exception
     */
    public static AmtDatabase createAmtDBFromAllPepXml(File allPepXmlFile,
                                                       int scanOrTimeMode,
                                                       FeatureSet.FeatureSelector featureSelector,
                                                       boolean robustRegression,
                                                       double maxStudentizedResidualForRegression,
                                                       double maxStudentizedResidualForInclusion)
            throws Exception
    {
        return createAmtDBFromPepXml(allPepXmlFile, null, scanOrTimeMode,
                                     featureSelector, robustRegression,
                                     maxStudentizedResidualForRegression,
                                     maxStudentizedResidualForInclusion);
    }

    /**
     * Create an AMT database from a pep.xml file.
     *
     * If no mzXml file is specified, look in the pepXml fractions themselves
     * to find mzXML files
     *
     * Warning:
     * Since PepXmlLoader runs through the pepXml file sequentially, we have to
     * load the databases one at a time... so if an mzXml file is missing
     * at the very end, we won't find it until we've done most of the work.
     *
     * TODO: currently this will fail if you're on an OS that uses a different
     * file separator character from the one the all.pep.xml file was created on,
     * because fraction.getSpectrumPath() will have the wrong characters.  Should be converted
     * @param allPepXmlFile
     * @param scanOrTimeMode
     * @return
     * @throws Exception
     */
    public static AmtDatabase createAmtDBFromPepXml(File allPepXmlFile,
                                                    MSRun run,
                                                    int scanOrTimeMode,
                                                    FeatureSet.FeatureSelector featureSelector,
                                                    boolean robustRegression,
                                                    double maxStudentizedResidualForRegression,
                                                    double maxStudentizedResidualForInclusion)
            throws Exception
    {
        AmtDatabase resultDB = new AmtDatabase();
        PepXmlLoader pepXmlLoader = new PepXmlLoader(allPepXmlFile, _log);

        PepXmlLoader.FractionIterator fractionIterator =
                pepXmlLoader.getFractionIterator();
        File pepXmlFileDir = allPepXmlFile.getParentFile();

        boolean shouldLoadRuns =
                (run == null && scanOrTimeMode == RegressionUtilities.REGRESSION_MODE_TIME);

        int numFeaturesChosenForRegression=0;
        int numFeaturesChosenForInclusion=0;
        Pair<Integer,Integer> numFeaturesDummy = new Pair<Integer,Integer>(0,0);
        int numFeaturesEvaluated = 0;
        //for each fraction, create a featureset and locate the mzXml file,
        //create an amt database, and add all the observations, runs and mods
        PepXMLFeatureFileHandler pepXmlFileHandler = new PepXMLFeatureFileHandler();
        while (fractionIterator.hasNext())
        {
            PepXmlLoader.PepXmlFraction fraction =
                    (PepXmlLoader.PepXmlFraction) fractionIterator.next();
            _log.debug("fraction " + fraction.getDataBasename());
            FeatureSet fractionFeatureSet =
                    pepXmlFileHandler.createFeatureSetFromPepXMLFraction(fraction, pepXmlLoader);

            fractionFeatureSet = fractionFeatureSet.filter(featureSelector);
            //remove all but the first occurrence of each peptide

            MS2ExtraInfoDef.setFeatureSetModifications(fractionFeatureSet,
                    fraction.getModifications().toArray(new MS2Modification[0]));
            if (shouldLoadRuns)
            {
                File mzXmlFile = new File(pepXmlFileDir.getAbsolutePath() +
                        File.separatorChar +
                        fraction.getSpectrumPath());
                _log.debug("mzxml file: " + mzXmlFile.getAbsolutePath());

                run = MSRun.load(mzXmlFile.getAbsolutePath());
            }

            if (scanOrTimeMode == RegressionUtilities.REGRESSION_MODE_TIME)
                fractionFeatureSet.populateTimesForMS2Features(run);

            numFeaturesEvaluated += fractionFeatureSet.getFeatures().length;

            AmtDatabase thisRunDB =
                    createAmtDatabaseForRun(fractionFeatureSet,
                            scanOrTimeMode, robustRegression, numFeaturesDummy,
                            maxStudentizedResidualForRegression,
                            maxStudentizedResidualForInclusion);
            _log.debug("Created fraction database: " + thisRunDB);

            numFeaturesChosenForRegression += numFeaturesDummy.first;
            numFeaturesChosenForInclusion += numFeaturesDummy.second;

            resultDB.addObservationsFromAnotherDatabase(thisRunDB);
            thisRunDB.getRunBySequence(1).setMzXmlFilename(run.getFileName());
            thisRunDB.getRunBySequence(1).setPepXmlFilename(fraction.getDataBasename() + ".pep.xml");
        }
        _log.debug("Features chosen for final regression: " + numFeaturesChosenForRegression);
        _log.debug("Features chosen for inclusion: " + numFeaturesChosenForInclusion);
        _log.debug("Total features evaluated: " + numFeaturesEvaluated);
        _log.debug("Resulting database: " + resultDB);

        return resultDB;
    }






    /**
     * Create a feature set based on this AMT database, mapped to a particular run.
     * Make a BUNCH of assumptions.
     * @param run
     * @param regressionMS2FeatureSet
     * @return
     */
    public static FeatureSet buildFeatureSet(AmtDatabase amtDB,
                                             MSRun run, FeatureSet regressionMS2FeatureSet,
                                             int chargeForAllFeatures, double minPeptideProphet,
                                             boolean robustRegression)
    {
        ArrayList<Feature> resultFeaturesList = new ArrayList<Feature>();

        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet((float)minPeptideProphet);

        regressionMS2FeatureSet = regressionMS2FeatureSet.filter(sel);

        regressionMS2FeatureSet.populateTimesForMS2Features(run);

        Map<String,Double> regressionMap =
            AmtUtilities.calculateHydrophobicityScanOrTimeRelationship(
                    regressionMS2FeatureSet.getFeatures(),
                    RegressionUtilities.REGRESSION_MODE_TIME,
                    robustRegression);

        double slope = regressionMap.get(RegressionUtilities.REGRESSION_SLOPE_KEY);
        double intercept = regressionMap.get(RegressionUtilities.REGRESSION_INTERCEPT_KEY);

        for (AmtPeptideEntry peptideEntry : amtDB.getEntries())
        {
            int scan = (int) AmtUtilities.predictScanOrTime(slope, intercept,
                                             peptideEntry.getMedianObservedHydrophobicity());
            Feature newFeature =
                    MS2ExtraInfoDef.createMS2Feature(scan, (float) peptideEntry.getMass(), chargeForAllFeatures,
                                                     peptideEntry.getPeptideSequence(), null, null);
            AmtExtraInfoDef.setObservedHydrophobicity(newFeature,peptideEntry.getMedianObservedHydrophobicity());

            //fake support for our fake charge
            if (chargeForAllFeatures>1)
                newFeature.setPeaks(2);

            //fake intensity
            newFeature.setIntensity(200);

            resultFeaturesList.add(newFeature);
        }

        return new FeatureSet(resultFeaturesList.toArray(new Feature[0]));
    }


}
