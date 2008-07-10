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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.gui.util.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class SummarizeProtXmlCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SummarizeProtXmlCLM.class);

    File protXmlFile;

    public SummarizeProtXmlCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "summarizeprotxml";
        mShortDescription = "summarizeprotxml";
        mHelpMessage = "summarizeprotxml";
        CommandLineArgumentDefinition[] argDefs =
                {
                    createUnnamedFileArgumentDefinition(true, "protxml file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFile = getUnnamedFileArgumentValue();
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

            List<Float> probabilityList = new ArrayList<Float>();

            Map<String, Float> proteinRatioMap = new HashMap<String, Float>();

            Iterator<ProteinGroup> iterator = protXmlReader.iterator();
            while (iterator.hasNext())
            {
                ProteinGroup pg = iterator.next();
                probabilityList.add(pg.getProbability());
            }

            PanelWithHistogram pwh = new PanelWithHistogram(probabilityList, "Probabilities of Recovered Proteins");
            pwh.displayInTab();

            int numAbovePoint5 = 0;
            int numAbovePoint75 = 0;
            int numAbovePoint9 = 0;
            int numAbovePoint95 = 0;

            for (float probability : probabilityList)
            {
                if (probability >= .5f)
                {
                    numAbovePoint5++;
                    if (probability >= .75f)
                    {
                        numAbovePoint75++;
                        if (probability >= .9f)
                        {
                            numAbovePoint9++;
                            if (probability >= .95f)
                            {
                                numAbovePoint95++;
                            }
                        }
                    }
                }
            }

            ApplicationContext.infoMessage("Total protein groups: " + probabilityList.size()) ;
            ApplicationContext.infoMessage("\tProb >= .5: " + numAbovePoint5);
            ApplicationContext.infoMessage("\tProb >= .75: " + numAbovePoint75);
            ApplicationContext.infoMessage("\tProb >= .9: " + numAbovePoint9);
            ApplicationContext.infoMessage("\tProb >= .95: " + numAbovePoint95);


        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
