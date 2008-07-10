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
package org.fhcrc.cpl.viewer.qa.commandline;

import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.*;


/**
 * test
 */
public class QAProtXMLCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(QAProtXMLCLM.class);

    protected File allProtXmlFile;
    protected File protGeneFile;
    protected File outFile;
    protected float minProteinProphet=0f;


    public QAProtXMLCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "qaprotxml";

        mHelpMessage ="Perform QA analysis on a single protXML file";
        mShortDescription = "Perform QA analysis on a single protXML file";

        CommandLineArgumentDefinition[] argDefs =
               {
                       this.createUnnamedFileArgumentDefinition(true, "all.prot.xml filepath"),
                       createFileToReadArgumentDefinition("protgenefile", true,
                               "File associating gene symbols with protein accession numbers"),
                       createFileToWriteArgumentDefinition("out", true,
                               "Output File"),
                       createDecimalArgumentDefinition("minpprophet", false,
                               "Minimum ProteinProphet score to keep in output", minProteinProphet)
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        allProtXmlFile = this.getUnnamedFileArgumentValue();
        protGeneFile = getFileArgumentValue("protgenefile");
        outFile = getFileArgumentValue("out");
        minProteinProphet=(float) getDoubleArgumentValue("minpprophet");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            QAUtilities.createAllProtText(allProtXmlFile, protGeneFile, outFile, minProteinProphet);
            ApplicationContext.infoMessage("Wrote file " + outFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
