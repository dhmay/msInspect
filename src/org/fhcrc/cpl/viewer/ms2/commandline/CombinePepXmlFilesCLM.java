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
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToWriteArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.List;
import java.util.ArrayList;


/**
 * test
 */
public class CombinePepXmlFilesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CombinePepXmlFilesCLM.class);

    protected File[] inputFiles;
    protected File outFile;


    public CombinePepXmlFilesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "combinepepxmlfiles";

        mHelpMessage ="combinepepxmlfiles";
        mShortDescription = "combinepepxmlfiles";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedSeriesFileArgumentDefinition(true, "input files"),
                       new FileToWriteArgumentDefinition("out", true, "Output file"),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inputFiles = getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");

        for (File file : inputFiles)
            if (file.getAbsolutePath().equals(outFile.getAbsolutePath()))
                throw new ArgumentValidationException("ERROR: output file is also specified as an input file.  " +
                        "Quitting");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<File> inputFilesList = new ArrayList<File>();
        for (File inputFile : inputFiles)
            inputFilesList.add(inputFile);
        try
        {
            new PepXMLFeatureFileHandler().combinePepXmlFiles(inputFilesList, outFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
/*
        FileOutputStream outStream = null;
        try
        {
            outStream = new FileOutputStream(outFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to open output file",e);
        }

        for (int i=0; i<inputFiles.length; i++)
        {
            boolean includeHeader = false;
            boolean includeFooter = false;

            if (i==0)
                includeHeader = true;
            if (i==inputFiles.length-1)
                includeFooter = true;
            try
            {
                FileInputStream fis = new FileInputStream(inputFiles[i]);
                DataInputStream dis = new DataInputStream(fis);
                String line = null;

                ApplicationContext.setMessage("Processing file " + inputFiles[i].getAbsolutePath());

                boolean inHeader = true;
                boolean inFooter = false;
                while ((line = dis.readLine()) != null)
                {
                    if (line.contains("/msms_run_summary"))
                        inFooter = true;
                    else if (line.contains("msms_run_summary"))
                        inHeader = false;
                    if ( (!inHeader || includeHeader) && (!inFooter || includeFooter) )
                    {
                        outStream.write(line.getBytes());
                        outStream.write('\n');
                        outStream.flush();
                    }
                }

                dis.close();
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to open input file "
                        + inputFiles[i].getAbsolutePath(),e);
            }
        }

        try
        {
          outStream.close();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to close output file",e);
        }
*/
    }

}
