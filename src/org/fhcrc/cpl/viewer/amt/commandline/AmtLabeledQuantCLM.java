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
package org.fhcrc.cpl.viewer.amt.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.viewer.amt.AmtLabeledQuant;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;


import java.io.File;
import java.util.*;


/**
 * Command linemodule for
 */
public class AmtLabeledQuantCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AmtLabeledQuantCLM.class);

    protected File[] featureFiles = null;

    protected File outFile;
    protected File outDir;


    protected boolean showCharts = false;

    protected List<Float> logRatiosAllFiles = new ArrayList<Float>();

    protected boolean perCharge = true;

    float maxRatio = 15f;
    float minRatio = 1f/15f;

    AmtLabeledQuant amtLabeledQuant = null;


    public AmtLabeledQuantCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "amtlabeledquant";
        mShortDescription = "Performs labeled quantitation by mining AMT results for isotopic pairs.";
        mHelpMessage = "Performs labeled quantitation by mining AMT results for isotopic pairs.  " +
                       "Calculates ratios for each peptide within each charge state and then takes " +
                       "the geometric mean.";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "AMT matching result file(s)"),
                        new FileToWriteArgumentDefinition("out",false, null),
                        new DirectoryToWriteArgumentDefinition("outdir",false, null),
                        new DecimalArgumentDefinition("separation", false,
                                "light-heavy separation", AmtLabeledQuant.ACRYLAMIDE_MASS_DIFF),
                        new StringArgumentDefinition("residue", false,
                                "Labeled residue (leave blank for n-terminal)", "" + AmtLabeledQuant.ACRYLAMIDE_RESIDUE),
                        new DecimalArgumentDefinition("minratiohighpprophet", false,
                                "Minimum AMT probability for consideration in a ratio: the higher of the two " +
                                        "probabilities must be at least this high",
                                AmtLabeledQuant.DEFAULT_MIN_RATIO_HIGH_PROB),
                        new DecimalArgumentDefinition("minratiolowpprophet", false,
                                "Minimum AMT probability for consideration in a ratio: the lower of the two " +
                                        "probabilities must be at least this high",
                                AmtLabeledQuant.DEFAULT_MIN_RATIO_LOW_PROB),
                        new BooleanArgumentDefinition("percharge", false,
                                "Calculate ratios per charge (vs. all charges together)?", perCharge),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new DecimalArgumentDefinition("minratio", false, "Minimum ratio to keep", minRatio),
                        new DecimalArgumentDefinition("maxratio", false, "Maximum ratio to keep", maxRatio),
                        new ModificationListArgumentDefinition("othermods", false,
                                "a list of other modifications applied to sample",
                                AmtLabeledQuant.DEFAULT_OTHER_MODS),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();

        if (featureFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }
        else
        {
            assertArgumentPresent("out");
            assertArgumentAbsent("outdir");
        }

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");


        double separation = getDoubleArgumentValue("separation");


        String residue = getStringArgumentValue("residue");
        if (residue != null && residue.length() > 1)
            throw new ArgumentValidationException("Invalid residue " + residue);
        AnalyzeICAT.IsotopicLabel isotopicLabel = AmtLabeledQuant.DEFAULT_ISOTOPIC_LABEL;
        if (residue != null)
        {
             isotopicLabel = new AnalyzeICAT.IsotopicLabel((float) PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[residue.charAt(0)],
                                (float) (PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[residue.charAt(0)] + separation),
                                residue.charAt(0), 5);
        }

        ApplicationContext.infoMessage("Using separation " + separation + "Da, residue " + residue);

        perCharge = getBooleanArgumentValue("percharge");

        showCharts = getBooleanArgumentValue("showcharts");

        amtLabeledQuant = new AmtLabeledQuant();
        amtLabeledQuant.setIsotopicLabel(isotopicLabel);
        amtLabeledQuant.setMinRatioHigherProbability(getFloatArgumentValue("minratiohighpprophet"));
        amtLabeledQuant.setMinRatioLowerProbability(getFloatArgumentValue("minratiolowpprophet"));

        if (hasArgumentValue("othermods"))
        {
            MS2Modification[] specifiedMods =
                    getModificationListArgumentValue("othermods");
            if (specifiedMods != null)
            {
                for (MS2Modification specifiedMod : specifiedMods)
                {
                    _log.debug("Including user-specified modification: " + specifiedMod);
                }
                amtLabeledQuant.setOtherModifications(specifiedMods);
            }
        }

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (featureFiles.length == 1)
        {
            handleFile(featureFiles[0], outFile);
        }
        else
        {
            for (File featureFile : featureFiles)
            {
                File outputFile = new File(outDir, featureFile.getName() + ".amtquant.pep.xml");
                ApplicationContext.infoMessage("Processing file " + featureFile.getName() + "...");
                handleFile(featureFile, outputFile);
            }

            if (showCharts)
            {
                PanelWithHistogram pwh = new PanelWithHistogram(logRatiosAllFiles,"Log ratios");
                ChartDialog cd = new ChartDialog(pwh);
                cd.setVisible(true);

            }
        }
    }


    protected void handleFile(File featureFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        FeatureSet featureSet = null;
        try
        {
            featureSet = new FeatureSet(featureFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to load feature file " + featureFile.getAbsolutePath(),e);
        }

        ApplicationContext.setMessage("Read " + featureSet.getFeatures().length +
                " features from file " + featureFile.getAbsolutePath());

        amtLabeledQuant.quantitate(featureSet);

        try
        {
            featureSet.savePepXml(outputFile);
            ApplicationContext.setMessage("Wrote " + featureSet.getFeatures().length + " features to file " +
                    outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        if (showCharts)
        {
//            if (allLogRatios.size() > 0)
//            {
//                if (featureFiles.length == 1)
//                {
//                    PanelWithHistogram pwh = new PanelWithHistogram(allLogRatios,"Log ratios");
//                    ChartDialog cd = new ChartDialog(pwh);
//                    cd.setVisible(true);
//                }
//                else
//                    logRatiosAllFiles.addAll(allLogRatios);
//            }
//            else ApplicationContext.infoMessage("Failed to build histogram, no ratios to plot");
        }
    }
}
