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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.MzXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.Window2DFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureMassCalibrationUtilities;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFindingBroker;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.*;

/**
 * Command line module for correcting mzXML precursor masses
 */
public class CorrectPrecursorMzCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CorrectPrecursorMzCLM.class);

    protected File outFile = null;
    protected File outDir;

    protected File[] inFiles;


    FeatureSet featureSet = null;

    File featuresDir;

    protected File outFeatureFile = null;

    //calibration
    protected boolean calibrateSpectra=false;
    protected boolean calibrate = false;
    protected int numPartitions=1;
    protected int initialMassFilterPPM = 0;

    protected float fractionalMzTolerancePPM = 100f;

    protected boolean showCharts = false;

    protected boolean useInaccurateMs1Mz = false;


    //how many peaks above and below the precursor should we check for MS1 features?
    int numPeaksToCheckBelowPrecursor = 3;
    int numPeaksToCheckAbovePrecursor = 1;


    public CorrectPrecursorMzCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "correctprecursormz";

        mHelpMessage =
                "Corrects the precursor m/z values in mzXML spectra, using the features found by msInspect " +
                "(or any compatible feature-finder) to override the m/z values provided by the instrument. " +
                "Optionally, also performs calibration of the precursor masses, and of the spectra themselves. " +
                "Can be called using an existing feature file, or, if none is specified, will perform feature-finding. " +
                "Can be called on one file at a time, or on multiple at once.";

        mShortDescription = "Perform mass correction on mzXML precursor masses";

        CommandLineArgumentDefinition[] argDefs =
            {
                    this.createUnnamedSeriesFileArgumentDefinition(true, "Input mzXML file(s)"),
                    new FileToWriteArgumentDefinition("out", false, "Output File"),
                    new DirectoryToWriteArgumentDefinition("outdir", false, "Output Directory (for multiple inputs)"),
                    new FeatureFileArgumentDefinition("features", true,
                            "Feature file to use in correction"),
                    new DirectoryToReadArgumentDefinition("featuresdir", false,
                            "Directory of feature files to use in correction"),
                    new FileToWriteArgumentDefinition("outfeatures", false,
                            "Output recalibrated feature file"),
                    new BooleanArgumentDefinition("calibratespectra", false,
                            "Calibrate spectra, as well as precursor masses (only valid with calibrate option)?",
                            calibrateSpectra),
                    new IntegerArgumentDefinition("partitions", false,
                            "Number of partitions, by scan, to divide the run into, for calibration", numPartitions),
                    new IntegerArgumentDefinition("initialfilterppm", false,
                            "Initial ppm value used as a pre-calibration cutoff.  Features deviating from theoretical " +
                            "clusters (BEFORE calibration) will be filtered out during calibration.  " +
                            "However, those features WILL appear in the recalibrated featureset, with corrected masses." +
                            "Default = no filter",
                            initialMassFilterPPM),
                    new BooleanArgumentDefinition("calibrate", false, "Calibrate precursor masses?", calibrate),
                    new DecimalArgumentDefinition("fracmztolerance", false,
                            "fractional M/Z tolerance, in PPM, for associating MS1 features with precursor scans",
                            fractionalMzTolerancePPM),
                    new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                    new BooleanArgumentDefinition("useinaccuratems1mz", false,
                            "Use an MS1 feature if m/z not \"accurate\"?", useInaccurateMs1Mz),
                    new IntegerArgumentDefinition("peaksaboveprecursor", false,
                            "Number of peaks above the precursor m/z to check for MS1 features",
                            numPeaksToCheckAbovePrecursor),
                    new IntegerArgumentDefinition("peaksbelowprecursor", false,
                            "Number of peaks below the precursor m/z to check for MS1 features",
                            numPeaksToCheckBelowPrecursor),

            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFiles = this.getUnnamedSeriesFileArgumentValues();

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");
        featuresDir = getFileArgumentValue("featuresdir");

        if (inFiles.length == 1)
        {
            assertArgumentPresent("out");
            assertArgumentAbsent("outdir");
            assertArgumentAbsent("featuresdir");
            outFile = getFileArgumentValue("out");
        }
        else
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
            outDir = getFileArgumentValue("outdir");
        }

        outFeatureFile = getFileArgumentValue("outfeatures");

        featureSet = getFeatureSetArgumentValue("features");

        numPeaksToCheckBelowPrecursor = getIntegerArgumentValue("peaksbelowprecursor");
        numPeaksToCheckAbovePrecursor = getIntegerArgumentValue("peaksaboveprecursor");            


        numPartitions = getIntegerArgumentValue("partitions");

        //This is because I'm lame.  I could spend a little time and get this working
        //for multiple partitions
        if (numPartitions > 1 && hasArgumentValue("outfeatures"))
            throw new ArgumentValidationException("Output feature file is only allowed if there is only one partition");

        calibrate = getBooleanArgumentValue("calibrate");          
        calibrateSpectra = getBooleanArgumentValue("calibratespectra");

        if (calibrateSpectra && !calibrate)
            throw new ArgumentValidationException("Can't calibrate spectra if we're not recalibrating precursors");

        initialMassFilterPPM = getIntegerArgumentValue("initialfilterppm");
        fractionalMzTolerancePPM = getFloatArgumentValue("fracmztolerance");

        useInaccurateMs1Mz = getBooleanArgumentValue("useinaccuratems1mz");

        showCharts = getBooleanArgumentValue("showcharts");

    }

    public void execute() throws CommandLineModuleExecutionException
    {
        for (File inFile : inFiles)
        {
            try
            {
                MSRun run = MSRun.load(inFile.getAbsolutePath());
                if (run == null)
                    throw new CommandLineModuleExecutionException(TextProvider.getText("ERROR_LOADING_FILE"));
                File outputFile = outFile;
                if (outFile == null)
                    outputFile = CommandLineModuleUtilities.createOutputFile(inFile, "corrected.mzXML", outDir);
                FeatureSet featureSetForCorrection = featureSet;

                if (inFiles.length > 1)
                {
                    ApplicationContext.infoMessage("Correcting input file " + inFile.getName() + " to file " + outputFile.getName());

                    if (featuresDir != null)
                    {
                        File featureFileForCorrection =
                                CommandLineModuleUtilities.findFileLikeFile(inFile, featuresDir, "tsv");
                        ApplicationContext.infoMessage("Found feature file " + featureFileForCorrection.getName());                        
                        featureSetForCorrection = new FeatureSet(featureFileForCorrection);
                    }
                }

                handleFile(run, outputFile, featureSetForCorrection);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Error processing file", e);
            }
        }
    }


    /**
     * do the actual work
     */
    public void handleFile(MSRun run, File outputFile, FeatureSet featureSetForCorrection)
            throws CommandLineModuleExecutionException
    {
        if (featureSetForCorrection == null)
        {
            try
            {
                ApplicationContext.setMessage("Finding features...");
                FloatRange range = FeatureExtractor.getMzExtractionRange(run);
                float thisRunMinMz = range.min;
                float thisRunMaxMz = range.max;
                int thisRunScanCount = run.getScanCount();
//                FeatureExtractor find =
//                        FeatureExtractor.getDefault   (run, 1, thisRunScanCount, 6,
//                                new FloatRange(thisRunMinMz, thisRunMaxMz), 2.0);
//
//                if (null != outputFile)
//                {
//                    find.setStatusListener(new FeatureExtractor.StatusListener()
//                    {
//                        public void progress(float percent)
//                        {
//                            ApplicationContext.infoMessage(String.valueOf(percent) +
//                                    "% complete.");
//                        }
//                    });
//                }
//                featureSetForCorrection = find.analyze();

                featureSetForCorrection = FeatureFindingBroker.findPeptides(
                        run, 0, thisRunScanCount,
                        PeakCombiner.DEFAULT_MAX_CHARGE,
                        new FloatRange(thisRunMinMz, thisRunMaxMz),
                        0, FeatureFinder.DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS,
                        FeatureFinder.DEFAULT_FEATURE_FINDING_CLASS,
                        true, false, false);
            }
            catch (Exception e)
            {
                ApplicationContext.infoMessage("Error while finding features");
                throw new CommandLineModuleExecutionException(e);
            }
            ApplicationContext.setMessage("Done finding features.");
        }
        ApplicationContext.setMessage("Finding wavelength and offset...");

        List<Feature> ms1FeaturesByScan = new ArrayList<Feature>();
        for (Feature feature : featureSetForCorrection.getFeatures())
            ms1FeaturesByScan.add(feature);

        Collections.sort(ms1FeaturesByScan, new Feature.ScanAscComparator());
        if (initialMassFilterPPM > 0)
        {
            ms1FeaturesByScan =
                    FeatureMassCalibrationUtilities.filterFeaturesByMassDefectDeviation(
                            ms1FeaturesByScan, initialMassFilterPPM);
        }
        //now we've got a featureset, one way or another.  Determine the adjustment

        Pair<Integer, Pair<Double,Double>>[] calibrationParameters = null;
        if (calibrate)
        {
            calibrationParameters = calibrate(ms1FeaturesByScan.toArray(new Feature[ms1FeaturesByScan.size()]),
                                              featureSetForCorrection, run);
        }

        FeatureSet ms2ScanFeatureSet = run.getTandemFeatureSet(2);
        for (Feature feature : ms2ScanFeatureSet.getFeatures())
        {
            feature.setScanFirst(feature.getScan());
            feature.setScanLast(feature.getScan());
        }



        List<Integer> precursorAdjustmentsToTry = new ArrayList<Integer>();
        for (int i=0; i<=numPeaksToCheckBelowPrecursor; i++)
            precursorAdjustmentsToTry.add(-i);
        for (int i=1; i<=numPeaksToCheckAbovePrecursor; i++)
            precursorAdjustmentsToTry.add(i);

        List<Feature> remainingMs2ScanFeatures = new ArrayList<Feature>();
        for (Feature feature : ms2ScanFeatureSet.getFeatures())
            remainingMs2ScanFeatures.add(feature);
        List<Feature> remainingMs1Features = new ArrayList<Feature>();
        for (Feature feature : ms1FeaturesByScan)
            remainingMs1Features.add(feature);


        List<Float> oldMs2PrecursorMasses = new ArrayList<Float>();
        List<Float> newMs2PrecursorMasses = new ArrayList<Float>();
        List<Feature> matchedMs2ScanFeatures = new ArrayList<Feature>();
        Map<Integer,Feature> scanNumAdjustedFeatureMap = new HashMap<Integer,Feature>();

        for (Integer precursorAdjustment : precursorAdjustmentsToTry)
        {
            ApplicationContext.setMessage("Trying precursor adjustment: " + precursorAdjustment + "Da");
            List<Feature> adjustedMs2ScanFeatures = new ArrayList<Feature>();
            Map<Feature, Feature> adjustedOrigPrecursorFeatureMap = new HashMap<Feature,Feature>();
            for (Feature origFeature : remainingMs2ScanFeatures)
            {
                Feature adjustedFeature = new Feature(origFeature);
                adjustedFeature.setMass(origFeature.getMass() + precursorAdjustment);
                adjustedFeature.updateMz();

                adjustedMs2ScanFeatures.add(adjustedFeature);
                adjustedOrigPrecursorFeatureMap.put(adjustedFeature, origFeature);
            }

            Map<Feature, Feature> adjustedPrecursorMs1FeatureMatch =
                    assignMatches(run, outputFile, remainingMs1Features, adjustedMs2ScanFeatures);
            ApplicationContext.setMessage("\tMatches made: " + adjustedPrecursorMs1FeatureMatch.size());

            for (Feature adjustedPrecursorToFix : adjustedPrecursorMs1FeatureMatch.keySet())
            {
                Feature precursorToAdjust = adjustedOrigPrecursorFeatureMap.get(adjustedPrecursorToFix);

                Feature matchedMs1Feature = adjustedPrecursorMs1FeatureMatch.get(adjustedPrecursorToFix);
                oldMs2PrecursorMasses.add(precursorToAdjust.getMass());
                newMs2PrecursorMasses.add(matchedMs1Feature.getMass());
                matchedMs2ScanFeatures.add(precursorToAdjust);

                precursorToAdjust.setMz(matchedMs1Feature.getMz());
                precursorToAdjust.updateMass();
                scanNumAdjustedFeatureMap.put(precursorToAdjust.getScan(), precursorToAdjust);

                remainingMs2ScanFeatures.remove(precursorToAdjust);
            }

            remainingMs1Features.removeAll(adjustedPrecursorMs1FeatureMatch.values());

            ApplicationContext.setMessage("\tremaining precursors: " + remainingMs2ScanFeatures.size() + ", MS1: " +
                    remainingMs1Features.size());
        }

        //mark MS2 scans as altered
        for (MSRun.MSScan ms2Scan : run.getMS2Scans())
        {
            if (scanNumAdjustedFeatureMap.containsKey(ms2Scan.getNum()))
            {
                ms2Scan.setPrecursorMzCorrected(true);
                ms2Scan.setPrecursorMz(scanNumAdjustedFeatureMap.get(ms2Scan.getNum()).getMz());
            }
        }

        ApplicationContext.setMessage("Done matching.  " + scanNumAdjustedFeatureMap.size() + " out of " + 
                                      ms2ScanFeatureSet.getFeatures().length + " precursor m/z values updated");


        if (showCharts)
        {
            Feature[] ms2ScanFeatures = ms2ScanFeatureSet.getFeatures();
            float[] ms2Scans = new float[ms2ScanFeatures.length];
            float[] ms2Mzs = new float[ms2Scans.length];
            for (int i=0; i<ms2Scans.length; i++)
            {
                ms2Scans[i] = ms2ScanFeatures[i].getScan();
                ms2Mzs[i] = ms2ScanFeatures[i].getMz();
            }

            float[] ms1Scans = new float[ms1FeaturesByScan.size()];
            float[] ms1Mzs = new float[ms1FeaturesByScan.size()];
            for (int i=0; i<ms1FeaturesByScan.size(); i++)
            {
                ms1Scans[i] = ms1FeaturesByScan.get(i).getScan();
                ms1Mzs[i] = ms1FeaturesByScan.get(i).getMz();
            }

            PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
//            pwsp.setPointSize(6);
            pwsp.setShowLegend(true);
            if (!matchedMs2ScanFeatures.isEmpty())
            {
                float[] matchedScans = new float[matchedMs2ScanFeatures.size()];
                float[] matchedMzs = new float[matchedMs2ScanFeatures.size()];
                for (int i=0; i<matchedMs2ScanFeatures.size(); i++)
                {
                    matchedScans[i] = matchedMs2ScanFeatures.get(i).getScan();
                    matchedMzs[i] = matchedMs2ScanFeatures.get(i).getMz();
                }
                pwsp.addData(matchedScans, matchedMzs, "Matches");
            }
            pwsp.addData(ms2Scans, ms2Mzs, "MS2 scan precursors");
            pwsp.addData(ms1Scans, ms1Mzs, "MS1 features");

            pwsp.setName("Features");
            pwsp.displayInTab();


            if (!oldMs2PrecursorMasses.isEmpty())
            {
                float[] fracMzChangeDataPPM = new float[oldMs2PrecursorMasses.size()];
                float[] mzChangeDataDa = new float[oldMs2PrecursorMasses.size()];

                for (int i=0; i<oldMs2PrecursorMasses.size(); i++)
                {
                    float mzDiff = newMs2PrecursorMasses.get(i) - oldMs2PrecursorMasses.get(i);

                    mzChangeDataDa[i] = (mzDiff);
                    fracMzChangeDataPPM[i] = (float)((mzDiff + .5) % 1) - .5f;
                    if (fracMzChangeDataPPM[i] < -.5)
                        fracMzChangeDataPPM[i] += 1;
                    fracMzChangeDataPPM[i] *= 1000000f / oldMs2PrecursorMasses.get(i);
                }
                PanelWithHistogram mzChangeHist = new PanelWithHistogram(fracMzChangeDataPPM,
                        "Fractional mass adjustments (ppm)");
                mzChangeHist.displayInTab();

                PanelWithHistogram mzChangeHistDa = new PanelWithHistogram(mzChangeDataDa,
                        "Mass adjustments (Da)", 200);
                mzChangeHistDa.setName("Mass adjustments (Da)");
                mzChangeHistDa.displayInTab();
            }
        }


        try
        {
            ApplicationContext.setMessage("Writing corrected file...");
            MzXmlWriter mzXmlWriter = new MzXmlWriter(run);
            mzXmlWriter.setMassCalibrationParameters(calibrationParameters);
            mzXmlWriter.setShouldCalibrateSpectra(calibrateSpectra);
            mzXmlWriter.setShouldCalibratePrecursorMasses(false);
            mzXmlWriter.write(outputFile);
            ApplicationContext.infoMessage("Done.  Wrote file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

    }

    public Map<Feature, Feature> assignMatches(MSRun run, File outputFile, List<Feature> featuresForCorrection,
                                               List<Feature>  ms2ScanFeatures)
            throws CommandLineModuleExecutionException
    {
        Window2DFeatureSetMatcher featureSetMatcher = new Window2DFeatureSetMatcher();
        featureSetMatcher.setMassDiffType(FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
        featureSetMatcher.setMaxMassDiff(fractionalMzTolerancePPM);
        featureSetMatcher.setMinMassDiff(-fractionalMzTolerancePPM);
        featureSetMatcher.setElutionMode(Window2DFeatureSetMatcher.ELUTION_MODE_SCAN);
        featureSetMatcher.setElutionRangeMode(Window2DFeatureSetMatcher.ELUTION_RANGE_MODE_RANGE);
        featureSetMatcher.setMinElutionDiff(0);
        featureSetMatcher.setMaxElutionDiff(0);
        featureSetMatcher.setMatchWithinChargeOnly(true);

        FeatureSetMatcher.FeatureMatchingResult matchingResult =
                featureSetMatcher.matchFeatures(new FeatureSet(ms2ScanFeatures.toArray(new Feature[0])),
                                                new FeatureSet(featuresForCorrection.toArray(new Feature[0])));
        Set<Feature> matchedMs2ScanFeatures = matchingResult.getMasterSetFeatures();
//List<Float> massErrors = new ArrayList<Float>();
//List<Float> scanErrors = new ArrayList<Float>();
//for (Feature ms2Feature : matchedMs2ScanFeatures)
//{
//    for (Feature ms1Feature : matchingResult.getSlaveSetFeatures(ms2Feature))
//    {
//
//
//        massErrors.add(ms2Feature.getMz() - ms1Feature.getMz());
//        scanErrors.add((float) (ms2Feature.getScan() - ms1Feature.getScan()));
//    }
//}
//new PanelWithScatterPlot(scanErrors, massErrors, "Error").displayInTab();


        Map<Feature, Feature> result = new HashMap<Feature, Feature>();

        int numInaccurateMz = 0;
        int numTooManyMatches = 0;

        for (Feature matchedMs2ScanFeature : matchedMs2ScanFeatures)
        {
            _log.debug("Matched MS2 feature...");
            List<Feature> matchedMs1Features = matchingResult.getSlaveSetFeatures(matchedMs2ScanFeature);
            if (matchedMs1Features.size() == 1)
            {
                Feature matchedMs1Feature = matchedMs1Features.iterator().next();
                if (useInaccurateMs1Mz || matchedMs1Feature.isAccurateMZ())
                {
                    _log.debug("\tYEP! Matching: \n\t\t" + matchedMs2ScanFeature + "\n\t\t" + matchedMs1Feature);
                    result.put(matchedMs2ScanFeature, matchedMs1Feature);
                }
                else
                {
                    _log.debug("\tNOPE.  Inaccurate MS1 mass");
                    numInaccurateMz++;
                }

            }
            else
            {
                _log.debug("\tNOPE.  Too many matches");
                numTooManyMatches++;
            }
        }

        ApplicationContext.setMessage("Adjusted M/Z of " + result.size() + " out of " +
                                       ms2ScanFeatures.size() + " precursors");
        ApplicationContext.setMessage("\tNumber with too many matches: " + numTooManyMatches);
        if (!useInaccurateMs1Mz)
            ApplicationContext.setMessage("\tNumber with inaccurate MS1 mass: " + numInaccurateMz);
        return result;
    }

    /**
     * Calibrates the features, MS2 precursor masses
     * @param featuresByScan
     * @param ms1FeatureSetForCorrection
     * @param run
     * @throws CommandLineModuleExecutionException
     */
    protected Pair<Integer, Pair<Double,Double>>[] calibrate(Feature[] featuresByScan, FeatureSet ms1FeatureSetForCorrection,
                                               MSRun run)
            throws CommandLineModuleExecutionException
    {
        Pair<Integer, Pair<Double,Double>>[] calibrationParameters =
                FeatureMassCalibrationUtilities.calculateWavelengthsAndOffsetsMultiplePartitions(
                        featuresByScan,
                        MassCalibrationUtilities.DEFAULT_MAX_PAIRS_FOR_LEVERAGE_CALC,
                        numPartitions,
                        MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH,
                        showCharts);
            
        if (_log.isDebugEnabled())
        {
            for (Pair<Integer, Pair<Double,Double>> partitionParameters : calibrationParameters)
                _log.debug("**partition parameters: "+ partitionParameters.first + ",   " + partitionParameters.second.first + ", " + partitionParameters.second.second);
        }
        if (calibrationParameters.length==1)
        {
            double wavelength = calibrationParameters[0].second.first;
            double offset = calibrationParameters[0].second.second;

            ApplicationContext.setMessage("Wavelength: " + wavelength +
                    ", offset: " + offset);
        }

        ApplicationContext.setMessage("Adjusting features...");
        Feature[] ms1Features = ms1FeatureSetForCorrection.getFeatures();
        for (Feature feature : ms1Features)
        {
            float newMass = feature.getMass();
            double wavelength = calibrationParameters[0].second.first;
            double offset = calibrationParameters[0].second.second;
            newMass += feature.getMass() * (MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH - wavelength) -
                    offset;

            feature.setMass(newMass);
            if (feature.getCharge() > 0)
                feature.updateMz();
        }

        if (outFeatureFile != null)
        {

            try
            {
                ms1FeatureSetForCorrection.save(outFeatureFile);
                ApplicationContext.infoMessage("Done writing adjusted features to file " + outFeatureFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
        //TODO: make this work for partitions
        for (MSRun.MSScan ms2Scan : run.getMS2Scans())
        {
            int precursorCharge = ms2Scan.getPrecursorCharge();
            if (precursorCharge < 1)
                precursorCharge = 1;
            float oldPrecursorMz = ms2Scan.getPrecursorMz();
            float oldPrecursorMass = oldPrecursorMz * precursorCharge;
            float newPrecursorMass = oldPrecursorMass + (float)
                    (oldPrecursorMass *
                            (MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH -
                                    calibrationParameters[0].second.first) -
                            calibrationParameters[0].second.second);
            ms2Scan.setPrecursorMz(newPrecursorMass / precursorCharge);
        }

        if (showCharts)
        {
            PanelWithScatterPlot pwsp =
                    FeatureMassCalibrationUtilities.plotMassDefectDeviation(
                            featuresByScan,
                            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH, false, 0);
            pwsp.setName("After");
            pwsp.displayInTab();
        }

        return calibrationParameters;
    }

}
