/* 
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for feature finding
 */
public class FilterPeptideRegexpCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FilterPeptideRegexpCLM.class);

    protected File[] files;
    protected File outFile;
    protected File outDir;

    protected String regexp;

    public FilterPeptideRegexpCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "filterpeptideregexp";
        mHelpMessage = "Example: [^C]* will keep all peptides without a Cysteine residue";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(false, null),
                        new StringArgumentDefinition("regexp",true, "Regular expression that every peptide must match"),

                        new FileToWriteArgumentDefinition("out",false, null),
                        new DirectoryToWriteArgumentDefinition("outdir",false, null)

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        files = getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");
        regexp = getStringArgumentValue("regexp");

        if (files.length == 1)
        {
            assertArgumentPresent("out");
            assertArgumentAbsent("outdir");
        }
        else
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (files.length == 1)
        {
            processFile(files[0], outFile);
            ApplicationContext.infoMessage("Wrote output file " + outFile.getAbsolutePath());
        }
        else
        {
            for (File file : files)
            {
                File thisOutFile = new File(outDir, file.getName());
                processFile(file, thisOutFile);
                ApplicationContext.infoMessage("Wrote output file " + thisOutFile.getAbsolutePath());                
            }
        }
    }

    protected void processFile(File file, File currentOutFile)
            throws CommandLineModuleExecutionException
    {
        FeatureSet featureSet = null;
        try
        {
            featureSet = new FeatureSet(file);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        List<Feature> outFeatureList = new ArrayList<Feature>();

        ApplicationContext.infoMessage("Before filter: " + featureSet.getFeatures().length);

        for (Feature feature : featureSet.getFeatures())
        {
            boolean exclude = false;

            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                exclude=true;
            else
            {
                if (!peptide.matches(regexp))
                    exclude=true;
            }
            if (!exclude)
                outFeatureList.add(feature);
        }
        featureSet.setFeatures(outFeatureList.toArray(new Feature[outFeatureList.size()]));
        ApplicationContext.infoMessage("After filter: " + featureSet.getFeatures().length);

        featureSet.setSourceFile(currentOutFile);
        try
        {
            featureSet.save(currentOutFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
