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
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToReadArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.DecimalArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class FindIndistinguishableProteinsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FindIndistinguishableProteinsCLM.class);

    File protXmlFile;
    protected File protListFile;

    protected float minProteinProphet = 0.0f;

    public FindIndistinguishableProteinsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "findindistinguishableproteins";
        mShortDescription = "findindistinguishableproteins";
        mHelpMessage = "findindistinguishableproteins";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, "protxml file"),
                        new FileToReadArgumentDefinition("protfile", true, "File with list of proteins to look for"),
                        new DecimalArgumentDefinition("minpprophet", false,
                                "Minimum ProteinProphet score.", minProteinProphet),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFile = getUnnamedFileArgumentValue();

        protListFile = getFileArgumentValue("protfile");

        minProteinProphet = getFloatArgumentValue("minpprophet");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {


        try
        {
            Set<String> proteinNamesToLookFor = new HashSet<String>();
            BufferedReader br = new BufferedReader(new FileReader(protListFile));

            String line = null;

            while ((line = br.readLine()) != null)
            {
                proteinNamesToLookFor.add(line.trim());
            }
            ApplicationContext.infoMessage("Proteins to look for: " + proteinNamesToLookFor.size());

            List<ProteinGroup> proteinGroups = ProteinUtilities.loadProteinGroupsFromProtXML(protXmlFile);

            Set<String> foundProteins = new HashSet<String>();
            Set<String> notFoundProteins = new HashSet<String>(proteinNamesToLookFor);
            for (ProteinGroup group : proteinGroups)
            {
                if (group.getProbability() < minProteinProphet)
                    continue;
                Set<String> protNamesThisGroup = new HashSet<String>();
                for (ProtXmlReader.Protein protein : group.getProteins())
                {
                    protNamesThisGroup.add(protein.getProteinName());
                }
                List<String> justFoundProteins = new ArrayList<String>();
                for (String protName : notFoundProteins)
                {
                    if (protNamesThisGroup.contains(protName))
                    {
                        foundProteins.add(protName);
                        justFoundProteins.add(protName);
                    }
                }
                notFoundProteins.removeAll(justFoundProteins);


            }

            ApplicationContext.infoMessage("Found: " + foundProteins.size());
            ApplicationContext.infoMessage("Not Found: " + notFoundProteins.size());

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
