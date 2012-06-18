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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.BooleanArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.Date;


/**
 */
public class CreateIndexCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CreateIndexCommandLineModule.class);

    protected MSRun run1 = null;
    protected File[] files = null;

    protected boolean shouldForce = false;


    public CreateIndexCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "createindex";

        addArgumentDefinition(createUnnamedSeriesFileArgumentDefinition(true,
                                "A series of feature files to index"));
        addArgumentDefinition(new BooleanArgumentDefinition("force", false,
                "Should force index re-creation if index already exists?", shouldForce));
    }



    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        files = getUnnamedSeriesFileArgumentValues();
        shouldForce = getBooleanArgumentValue("force");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            for (File file : files)
            {
                Date beforeDate = new Date();
                if (shouldForce)
                {
                    String indexFileName = MSRun._indexName(file.getAbsolutePath());
                    File indexFile = new File(indexFileName);
                    if (indexFile.exists())
                    {
                        ApplicationContext.infoMessage("Deleting existing index file " + indexFileName);
                        indexFile.delete();
                    }
                }
                run1 = MSRun.load(file.getAbsolutePath());
                int secondsToCreate = (int) ((new Date().getTime() - beforeDate.getTime()) / 1000f);

                ApplicationContext.infoMessage("Created Index for file " + file.getAbsolutePath() + " in " +
                        secondsToCreate + " seconds");
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

    }
}
