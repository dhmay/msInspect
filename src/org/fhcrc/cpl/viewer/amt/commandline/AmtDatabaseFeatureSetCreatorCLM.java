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
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;


/**
 * Command linemodule for creating featuresets that represent AMT databases
 */
public class AmtDatabaseFeatureSetCreatorCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AmtDatabaseFeatureSetCreatorCLM.class);

    protected File outFile;
    protected float featureMassAdjustment = 0.0f;
    protected boolean adjustFeatureMasses = false;

    protected AmtDatabase amtDB;
    protected List<MS2Modification> ms2Modifications;


    public AmtDatabaseFeatureSetCreatorCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "createamtfeatureset";
        mShortDescription = "A tool for creating feature sets from amt databases";
        mHelpMessage = "This tool takes an AMT database as input and produces a feature set with one feature per peptide entry.\n"
                +"Intensity is set to an arbitrary value, and scan number is set to a value that relates linearly to median observed hydrophobicity";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, "AMT database for matching"),
                        new FileToWriteArgumentDefinition("out",false, "output filepath"),
                        new ModificationListArgumentDefinition("modifications", false,
                                "a list of modifications to match on"),
                        new DecimalArgumentDefinition("featuremassadjustment", false,
                                "Adjust the masses of all AMT database features  by this amount (in Daltons; " +
                                        "for false positive testing)", featureMassAdjustment),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {

        ms2Modifications = new ArrayList<MS2Modification>();
        MS2Modification[] specifiedMods =
                getModificationListArgumentValue("modifications");
        if (specifiedMods != null)
        for (MS2Modification specifiedMod : specifiedMods)
        {
            ms2Modifications.add(specifiedMod);
            _log.debug("Including user-specified modification: " + specifiedMod);            
        }

        _log.debug("Total modifications: " + ms2Modifications.size());

        File dbFile = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        outFile = getFileArgumentValue("out");
                                                   
        try
        {
            amtDB = new AmtXmlReader(dbFile).getDatabase();
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }
        adjustFeatureMasses = hasArgumentValue("featuremassadjustment");
        if (adjustFeatureMasses)
            featureMassAdjustment = getFloatArgumentValue("featuremassadjustment");
        
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        AmtDatabaseMatcher amtDatabaseMatcher = new AmtDatabaseMatcher();
        MS2Modification[] modArray =
                ms2Modifications.toArray(new MS2Modification[ms2Modifications.size()]);
        AmtDatabaseFeatureSetGenerator featureGenerator = new AmtDatabaseFeatureSetGenerator(amtDB);
        Feature[] amtDBFeatures =
                featureGenerator.createFeaturesForModifications(modArray);
        if (adjustFeatureMasses)
            for (Feature feature : amtDBFeatures)
            {
                feature.setMass(feature.getMass() + featureMassAdjustment);
            }

        FeatureSet outFeatures = new FeatureSet(amtDBFeatures);
        MS2ExtraInfoDef.setFeatureSetModifications(outFeatures, modArray);
        try
        {
            outFeatures.save(outFile);
            ApplicationContext.infoMessage("Wrote " + amtDBFeatures.length +
                                          " features to file " + outFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
