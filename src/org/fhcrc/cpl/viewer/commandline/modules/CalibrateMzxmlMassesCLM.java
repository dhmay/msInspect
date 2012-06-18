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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.MzXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureMassCalibrationUtilities;
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
import java.io.FileNotFoundException;
import java.util.Arrays;

/**
 * Command linemodule for Saving pieces of mzXML files
 */
public class CalibrateMzxmlMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalibrateMzxmlMassesCLM.class);

    protected File outFile = null;
    protected File outDir;

    protected File[] inFiles;
    protected Pair<Integer, Pair<Double,Double>>[] calibrationParameters;

    boolean shouldCalculateAdjustment = false;

    FeatureSet featureSet = null;
    FeatureSet scanChargeFeatureSet = null;

    File featuresDir;
    File scanChargeFeaturesDir;

    protected File outFeatureFile = null;

    protected boolean ms2PrecursorMassesOnly=false;

    protected int numPartitions=1;

    protected int initialMassFilterPPM = 0;


    public CalibrateMzxmlMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "calibratemzxmlmasses";

        mHelpMessage =
                "Calibrates the masses in mzXML spectra, based on supplied parameters, or on the calibration properties of a specified feature file\n" +
                "This should work well when mass miscalibration is a function of M/Z (e.g., in some TOF data).  If it's a function of mass that behaves differently for different charge states (e.g., some LTQFT data), this will not give good results.";

        mShortDescription = "Perform mass calibration on mzXML spectra";

        CommandLineArgumentDefinition[] argDefs =
            {
                    this.createUnnamedSeriesFileArgumentDefinition(true, "Input mzXML file(s)"),
                    new FileToWriteArgumentDefinition("out", false, "Output File"),
                    new DirectoryToWriteArgumentDefinition("outdir", false, "Output Directory (for multiple inputs)"),
                    new DecimalArgumentDefinition("wavelength", false, "Wavelength"),
                    new DecimalArgumentDefinition("offset", false, "Offset"),
                    new FeatureFileArgumentDefinition("features", false,
                            "Feature file to use in recalibration"),
                    new FeatureFileArgumentDefinition("scanchargefeatures", false,
                            "Feature file to use to assign charges to MS/MS scans"),
                    new DirectoryToReadArgumentDefinition("featuresdir", false,
                            "Directory of feature files to use in recalibration"),
                    new DirectoryToReadArgumentDefinition("scanchargefeaturesdir", false,
                            "Directory of feature files to use for determining charge states"),
                    new FileToWriteArgumentDefinition("outfeatures", false,
                            "Output recalibrated feature file"),
                    new BooleanArgumentDefinition("onlyms2precursormasses", false,
                            "Only recalibrate MS2 precursor masses", ms2PrecursorMassesOnly),
                    new IntegerArgumentDefinition("partitions", false,
                            "Number of partitions by scan", numPartitions),
                    new IntegerArgumentDefinition("initialfilterppm", false,
                            "Initial ppm value used as a pre-calibration cutoff.  Features deviating from theoretical clusters (BEFORE calibration) will be filtered out during calibration.  However, those features WILL appear in the recalibrated featureset, with corrected masses.  Default = no filter",
                            initialMassFilterPPM)
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
        scanChargeFeaturesDir = getFileArgumentValue("scanchargefeaturesdir");



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


        if (hasArgumentValue("wavelength"))
        {
            assertArgumentPresent("offset");
            assertArgumentAbsent("features");
            double wavelength = getFloatArgumentValue("wavelength");
            double offset = getFloatArgumentValue("offset");

            //each pair contains an integer that's the first scan number of the partition,
            //and a pair of doubles representing the parameters for that partition
            calibrationParameters = (Pair<Integer, Pair<Double,Double>>[]) new Pair[1];
            calibrationParameters[0] = new Pair<Integer, Pair<Double,Double>>(0,
                    new Pair<Double,Double>(wavelength,offset));
        }
        else
        {
            assertArgumentAbsent("wavelength");
            if (inFiles.length == 1)
                assertArgumentPresent("features");
            else
                assertArgumentPresent("featuresdir");

            shouldCalculateAdjustment = true;
        }

        featureSet = getFeatureSetArgumentValue("features");
        scanChargeFeatureSet = getFeatureSetArgumentValue("scanchargefeatures");

        numPartitions = getIntegerArgumentValue("partitions");

        //This is because I'm lame.  I could spend a little time and get this working
        //for multiple partitions
        if (numPartitions > 1 && hasArgumentValue("outfeatures"))
            throw new ArgumentValidationException("Output feature file is only allowed if there is only one partition");

        ms2PrecursorMassesOnly = getBooleanArgumentValue("onlyms2precursormasses");

        initialMassFilterPPM = getIntegerArgumentValue("initialfilterppm");

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
                    outputFile = CommandLineModuleUtilities.createOutputFile(inFile, "calibrated.mzXML", outDir);
                FeatureSet featuresForCalibration = featureSet;
                FeatureSet scanChargeFeatures = scanChargeFeatureSet;

                if (inFiles.length > 1)
                {
                    ApplicationContext.infoMessage("Calibrating input file " + inFile.getName() + " to file " + outputFile.getName());

                    if (shouldCalculateAdjustment && featuresDir != null)
                    {
                        File featureFileForCalibration =
                                CommandLineModuleUtilities.findFileLikeFile(inFile, featuresDir, "tsv");
                        ApplicationContext.infoMessage("Found feature file " + featureFileForCalibration.getName());                        
                        featuresForCalibration = new FeatureSet(featureFileForCalibration);
                    }
                    if (shouldCalculateAdjustment && scanChargeFeaturesDir != null)
                    {

                        File featureFileForCalibration = null;
                        try
                        {
                             featureFileForCalibration = CommandLineModuleUtilities.findFileLikeFile(inFile, scanChargeFeaturesDir, "tsv");
                        }
                        catch (FileNotFoundException e)
                        {
                            featureFileForCalibration = CommandLineModuleUtilities.findFileLikeFile(inFile, scanChargeFeaturesDir, "xml");

                        }
                        ApplicationContext.infoMessage("Found scan-charge feature file " + featureFileForCalibration.getName());
                        scanChargeFeatures = new FeatureSet(featureFileForCalibration);
                    }
                }

                handleFile(run, outputFile, featuresForCalibration, scanChargeFeatures);
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
    public void handleFile(MSRun run, File outputFile, FeatureSet featureSetForCalibration,
                           FeatureSet scanChargeFeatures)
            throws CommandLineModuleExecutionException
    {
        if (shouldCalculateAdjustment)
        {
            if (featureSetForCalibration == null)
            {
                try
                {
                    ApplicationContext.setMessage("Finding features...");
                    FloatRange range = FeatureExtractor.getMzExtractionRange(run);
                    float thisRunMinMz = range.min;
                    float thisRunMaxMz = range.max;
                    int thisRunScanCount = run.getScanCount();
                    featureSetForCalibration = FeatureFindingBroker.findPeptides(
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

            Feature[] inputFeatures = featureSetForCalibration.getFeatures();
            if (initialMassFilterPPM > 0)
            {
                 inputFeatures =
                         FeatureMassCalibrationUtilities.filterFeaturesByMassDefectDeviation(
                                 inputFeatures, initialMassFilterPPM);
            }

            Feature[] featuresByScan = new Feature[inputFeatures.length];
            System.arraycopy(featureSetForCalibration.getFeatures(), 0, featuresByScan, 0,
                    inputFeatures.length);
            Arrays.sort(featuresByScan, new Feature.ScanAscComparator());
            //now we've got a featureset, one way or another.  Determine the adjustment

            calibrationParameters =
                    FeatureMassCalibrationUtilities.calculateWavelengthsAndOffsetsMultiplePartitions(
                            featuresByScan,
                            MassCalibrationUtilities.DEFAULT_MAX_PAIRS_FOR_LEVERAGE_CALC,
                            numPartitions,
                            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH,
                            false);




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
        }

        //TODO: get this working for multiple partitions
        if (outFeatureFile != null)
        {
            ApplicationContext.setMessage("Adjusting features...");
            Feature[] features = featureSetForCalibration.getFeatures();
            for (int i=0; i<features.length; i++)
            {
                Feature feature = features[i];
                float newMass = feature.getMass();
                double wavelength = calibrationParameters[0].second.first;
                double offset = calibrationParameters[0].second.second;
                newMass += feature.getMass() * (MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH - wavelength) -
                                offset;

                feature.setMass(newMass);
                if (feature.getCharge() > 0)
                    feature.updateMz();
            }
            try
            {
                featureSetForCalibration.save(outFeatureFile);
                ApplicationContext.infoMessage("Done writing adjusted features to file " + outFeatureFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }


        if (scanChargeFeatures != null)
        {
            int numCharge0 = 0;
            int numCharge1 = 0;
            int numCharge2 = 0;
            int numCharge3 = 0;

            for (Feature ms2Feature : scanChargeFeatures.getFeatures())
            {
                MSRun.MSScan ms2Scan = null;
                try
                {
                    ms2Scan = run.getMS2Scan(ms2Feature.getScan());
                }
                catch (ArrayIndexOutOfBoundsException e)
                {}
                
                if (ms2Scan != null)
                {
                    int charge = ms2Feature.getCharge();

                    if (charge > 0)
                        ms2Scan.setPrecursorCharge(charge);
                    switch (charge)
                    {
                        case 0:
                            numCharge0++;
                            break;
                        case 1:
                            numCharge1++;
                            break;
                        case 2:
                            numCharge2++;
                            break;
                        case 3:
                            numCharge3++;
                            break;
                    }
                }
            }
            ApplicationContext.infoMessage("Corrected MS2 scan charges.  Charge breakdown: 0=" +
                          numCharge0 + ", 1=" + numCharge1 + ", 2=" + numCharge2 + ", 3=" + numCharge3);
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


        try
        {
            ApplicationContext.setMessage("Writing calibrated file...");
            MzXmlWriter mzXmlWriter = new MzXmlWriter(run);  
//            mzXmlWriter.setMassCalibrationParameters(calibrationParameters);
//            mzXmlWriter.setShouldCalibrateSpectra(!ms2PrecursorMassesOnly);
            mzXmlWriter.write(outputFile);
            ApplicationContext.infoMessage("Done.  Wrote file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }


    }

}
