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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureMassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteDatabaseMatcher;
import org.apache.log4j.Logger;


import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for mass calibration
 */
public class CalibrateFeatureMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalibrateFeatureMassesCLM.class);

    protected File[] featureSetFiles;

    protected File outFile;
    protected File outDir;

    protected boolean showCharts = false;

    protected int maxPairsForLeverageCalc =
            MassCalibrationUtilities.DEFAULT_MAX_PAIRS_FOR_LEVERAGE_CALC;

    protected int numPartitions=1;

    protected double theoreticalMassWavelength =
            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH;


    protected int initiaMassFilterPPM = 0;


    public CalibrateFeatureMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "calibratefeaturemasses";
        mShortDescription = "Perform mass calibration on a set of high-mass-accuracy peptide features";
        mHelpMessage = "This tool is used to correct high-accuracy features for mass miscalibration, based on theoretical properties of the distribution of peptide masses.  Based on work by Wolski et al., 2006";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "Input feature file(s)"),
                        new FileToWriteArgumentDefinition("out",false,
                                "Output file (for single input files)"),
                        new DirectoryToWriteArgumentDefinition("outdir",false,
                                "Output directory (for multiple input files)"),
                        new BooleanArgumentDefinition("showcharts", false,
                                "Show charts", showCharts),
                        new IntegerArgumentDefinition("maxpairs", false,
                                "Maximum number of mass pairs to be considered (higher numbers may increase accuracy, but will take much longer)",
                                maxPairsForLeverageCalc),
                        new IntegerArgumentDefinition("partitions", false,
                                "Number of partitions to break the file into (higher numbers are more expensive)",
                                numPartitions),
                        new DecimalArgumentDefinition("theoreticalwavelength",false,
                                "Theoretical mass cluster wavelength",
                                theoreticalMassWavelength),
                        new IntegerArgumentDefinition("initialfilterppm", false,
                                "Initial ppm value used as a pre-calibration cutoff.  Features deviating from theoretical clusters (BEFORE calibration) will be filtered out during calibration.  However, those features WILL appear in the recalibrated featureset, with corrected masses.  Default = no filter",
                                initiaMassFilterPPM),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureSetFiles = this.getUnnamedSeriesFileArgumentValues();

        if (featureSetFiles.length == 1)
        {
            if (hasArgumentValue("out"))
                outFile = getFileArgumentValue("out");
            else
            {
                assertArgumentPresent("outdir");
                String outFileName = featureSetFiles[0].getName();
                outFile = new File(getFileArgumentValue("outdir"), outFileName); 
            }
        }
        else
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
            outDir = getFileArgumentValue("outdir");
        }

        showCharts = getBooleanArgumentValue("showcharts");
        maxPairsForLeverageCalc = getIntegerArgumentValue("maxpairs");

        numPartitions = getIntegerArgumentValue("partitions");
        theoreticalMassWavelength = getDoubleArgumentValue("theoreticalwavelength");

        initiaMassFilterPPM = getIntegerArgumentValue("initialfilterppm");
        
    }


    /**
     * do the actual work
     */
    public void execute()
            throws CommandLineModuleExecutionException
    {
        if (featureSetFiles.length == 1)
        {
            try
            {
                FeatureSet featureSet = new FeatureSet(featureSetFiles[0]);
                adjustFeatureFile(featureSet, outFile);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
        else
        {
            for (File featureSetFile : featureSetFiles)
            {
                try
                {
                    FeatureSet featureSet = new FeatureSet(featureSetFile);
                    String sourceFileName = featureSet.getSourceFile().getName();
                    ApplicationContext.setMessage("Processing file " + sourceFileName + " (" + featureSet.getFeatures().length + " features) ...");
                    adjustFeatureFile(featureSet, new File(outDir, sourceFileName));
                    ApplicationContext.setMessage("Done.");
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
            }
        }
    }

    public void adjustFeatureFile(FeatureSet featureSet, File outputFile)
            throws CommandLineModuleExecutionException
    {
        Feature[] featureArrayCopy = null;
        if (showCharts)
        {
            featureArrayCopy = new Feature[featureSet.getFeatures().length];
            for (int i=0; i<featureSet.getFeatures().length; i++)
                featureArrayCopy[i] = new Feature(featureSet.getFeatures()[i]);
        }


        FeatureMassCalibrationUtilities.calibrateMasses(featureSet.getFeatures(),
                maxPairsForLeverageCalc, numPartitions, theoreticalMassWavelength,
                initiaMassFilterPPM, showCharts);

        try
        {
            featureSet.save(outputFile);
            System.err.println("Saved output featureset");
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        if (showCharts)
        {
            List<Feature[]> featureArrayList = new ArrayList<Feature[]>();
            featureArrayList.add(featureArrayCopy);
            featureArrayList.add(featureSet.getFeatures());
            PanelWithScatterPlot pwsp =
                    FeatureMassCalibrationUtilities.plotMassDefectDeviation(
                            featureArrayList,
                            theoreticalMassWavelength, true, 200, true);
            pwsp.setName("Before/after");
            pwsp.displayInTab();
        }

    }

}
