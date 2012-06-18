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
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureMassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.APMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteDatabaseMatcher;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteUtilities;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 */
public class CalibrateMetMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalibrateMetMassesCLM.class);

    protected File[] featureSetFiles;
    protected MetaboliteDatabaseMatcher metDBMatcher;
    protected File outFile;
    protected File outDir;

    protected List<ChemicalCompound> databaseCompoundsByMass = null;

    protected TabLoader loader;

    protected boolean shouldRemoveImprobable = false;
    protected boolean showCharts = false;

    protected List<ChemicalCompound> databaseCompounds = null;

    protected float calMassTolerancePPM = MetaboliteDatabaseMatcher.DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM;

    protected float numStandardDeviationsUnitMassOffsetToKeep = 4;

    public CalibrateMetMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "calibratemetmasses";
        mShortDescription = "Calibrate metabolite feature masses based on a database match";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedSeriesFileArgumentDefinition(true, "Feature file(s)"),
                        new FileToReadArgumentDefinition("metdb", true, "Metabolite database file"),
                        new FileToWriteArgumentDefinition("out", false, "output file"),
                        new BooleanArgumentDefinition("showcharts", false,
                                "Show charts", showCharts),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "output directory"),
                        new BooleanArgumentDefinition("removeimprobable", false, "Remove improbable masses?", shouldRemoveImprobable),
                        new DecimalArgumentDefinition("masstoleranceppm", false, "Loose mass tolerance for calibration", calMassTolerancePPM),
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
        metDBMatcher = new MetaboliteDatabaseMatcher();
        metDBMatcher.setShowCharts(showCharts);

        calMassTolerancePPM = getFloatArgumentValue("masstoleranceppm");
        metDBMatcher.setCalMassTolerancePPM(calMassTolerancePPM);

        try
        {
            //populate 2 theoretical peaks for each compound
            databaseCompounds =
                    ChemicalCompound.loadCompoundsFromFile(getFileArgumentValue("metdb"), 2);
            ApplicationContext.infoMessage("Loaded " + databaseCompounds.size() + " compounds from database");
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            throw new ArgumentValidationException("Failed to load compound database " +
                    getFileArgumentValue("metdb").getAbsolutePath(),e);
        }
//for (ChemicalCompound compound : databaseCompoundsByMass) System.err.println(compound.getCommonestIsotopeMass());
        metDBMatcher.setDatabaseCompounds(databaseCompounds);


        showCharts = getBooleanArgumentValue("showcharts");
        shouldRemoveImprobable = getBooleanArgumentValue("removeimprobable");
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
                calibrateFeatureFile(featureSet, outFile);
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
                    calibrateFeatureFile(featureSet, new File(outDir, sourceFileName));
                    ApplicationContext.setMessage("Done.");
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException(e);
                }
            }
        }
    }

    protected void calibrateFeatureFile(FeatureSet featureSet, File outputFile)
            throws CommandLineModuleExecutionException
    {
        Arrays.sort(featureSet.getFeatures(), new Feature.MassAscComparator());
        
        try
        {
            metDBMatcher.calibrateFeatureMasses(featureSet.getFeatures());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed calibration for file " +
                    featureSet.getSourceFile().getName(),e);
        }

        if (shouldRemoveImprobable)
        {
            List<Double> databaseUnitMassOffsets = new ArrayList<Double>();
            List<Double> databaseMasses = new ArrayList<Double>();

            for (ChemicalCompound compound : databaseCompounds)
            {
                databaseUnitMassOffsets.add(MetaboliteUtilities.calcUnitMassDefectOffset(compound.getCommonestIsotopeMass()));
                databaseMasses.add(compound.getCommonestIsotopeMass());
            }
            double minInDatabase = BasicStatistics.min(databaseUnitMassOffsets);
            double maxInDatabase = BasicStatistics.max(databaseUnitMassOffsets);

            double meanInDatabase = BasicStatistics.mean(databaseUnitMassOffsets);
            double sdInDatabase = BasicStatistics.standardDeviation(databaseUnitMassOffsets);

            ApplicationContext.infoMessage("Database unit mass offsets: min=" + minInDatabase +
                    ", max=" + maxInDatabase + ", mean=" + meanInDatabase + ", SD=" + sdInDatabase);



            if (showCharts)
                new PanelWithScatterPlot(databaseMasses, databaseUnitMassOffsets, "DB mass vs unit mass offset").displayInTab();

            double minToKeep = meanInDatabase - (numStandardDeviationsUnitMassOffsetToKeep * sdInDatabase);
            double maxToKeep = meanInDatabase + (numStandardDeviationsUnitMassOffsetToKeep * sdInDatabase);
            ApplicationContext.infoMessage("Keeping offsets within " +
                    numStandardDeviationsUnitMassOffsetToKeep + "SD of mean: min=" + minToKeep +
                    ", max=" + maxToKeep);

            List<Feature> featuresToKeep = new ArrayList<Feature>();
            List<Double> allFeatureMasses = new ArrayList<Double>();
            List<Double> allFeatureUnitMassOffsets = new ArrayList<Double>();
            List<Double> keptFeatureMasses = new ArrayList<Double>();
            List<Double> keptFeatureUnitMassOffsets = new ArrayList<Double>();
            for (Feature feature : featureSet.getFeatures())
            {
                float mass = feature.mass;
                double unitMassOffset = MetaboliteUtilities.calcUnitMassDefectOffset(mass);
                allFeatureMasses.add((double) mass);
                allFeatureUnitMassOffsets.add(unitMassOffset);

                if (unitMassOffset >= minToKeep && unitMassOffset < maxToKeep)
                {
                    featuresToKeep.add(feature);
                    keptFeatureMasses.add((double) mass);
                    keptFeatureUnitMassOffsets.add(unitMassOffset);
                }
            }

            ApplicationContext.infoMessage("Kept " + featuresToKeep.size() + " out of " +
                    featureSet.getFeatures().length + " features");

            featureSet.setFeatures(featuresToKeep.toArray(new Feature[featuresToKeep.size()]));

            if (showCharts)
            {
                PanelWithScatterPlot pwspUnitMassOffsets = new PanelWithScatterPlot(
                        keptFeatureMasses, keptFeatureUnitMassOffsets, "mass vs mass offset");
                pwspUnitMassOffsets.addData(allFeatureMasses,allFeatureUnitMassOffsets, "all features");
                pwspUnitMassOffsets.displayInTab();
            }

        }


        try
        {
            APMLFeatureFileHandler.getSingletonInstance().saveFeatureSet(featureSet, outputFile);
            System.err.println("Saved output featureset " + outputFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }



    }
}


