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
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;


/**
 * Flip peptide ratios and light/heavy areas
 */
public class FlipPepXMLRatiosCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FlipPepXMLRatiosCLM.class);

    protected File[] pepXmlFiles;
    protected File outDir;
    protected File outFile;

    public FlipPepXMLRatiosCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "flippepxmlratios";

        mHelpMessage ="Flip peptide ratios and light/heavy areas";
        mShortDescription = "Flip peptide ratios and light/heavy areas";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedSeriesFileArgumentDefinition(true, "pepXML file(s) to flip"),
                       new FileToWriteArgumentDefinition("out", false, "Output file"),
                       new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        pepXmlFiles = getUnnamedSeriesFileArgumentValues();

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (outFile == null && outDir == null)
            throw new ArgumentValidationException("Either out or outdir is required");
        if (outFile != null && outDir != null)
            throw new ArgumentValidationException("out and outdir cannot both be specified");
        if (pepXmlFiles.length > 1 && outFile != null)
            throw new ArgumentValidationException(
                    "If multiple pepXML files are specified, an output directory must be given");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File pepXmlFile : pepXmlFiles)
        {
            handleFile(pepXmlFile);
        }
    }

    protected void handleFile(File pepXmlFile) throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Processing file " + pepXmlFile.getAbsolutePath());

        File outputFile = outFile;
        if (outputFile == null)
            outputFile = new File(outDir, pepXmlFile.getName());
        FeatureSet featureSet = null;

        try
        {
            featureSet = new FeatureSet(pepXmlFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failure opening file " +
                    pepXmlFile.getAbsolutePath(), e);
        }

        int numRatiosFlipped = 0;
        for (Feature feature : featureSet.getFeatures())
        {
            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
                numRatiosFlipped++;
                double oldHeavy = IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
                IsotopicLabelExtraInfoDef.setHeavyIntensity(feature,
                        IsotopicLabelExtraInfoDef.getLightIntensity(feature));
                IsotopicLabelExtraInfoDef.setLightIntensity(feature, oldHeavy);
                IsotopicLabelExtraInfoDef.setRatio(feature,
                        1.0 / IsotopicLabelExtraInfoDef.getRatio(feature));
            }
        }
        ApplicationContext.infoMessage("\tFlipped " + numRatiosFlipped + " ratios");

        try
        {
            featureSet.savePepXml(outputFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failure writing file " +
                    outputFile.getAbsolutePath(), e);
        }

        ApplicationContext.infoMessage("\tWrote output file " + outputFile.getAbsolutePath());
    }
}
