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

package org.fhcrc.cpl.toolbox.test.commandline;

import org.fhcrc.cpl.toolbox.test.BaseCommandTest;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;

import java.util.Map;
import java.util.HashMap;

import junit.framework.TestResult;

/**
 * Abstract base class for all msInspect commandline tests.  Provides convenient ways to invoke
 * command line modules
 */
public abstract class BaseCommandLineTest extends BaseCommandTest
{
    //map of argument values to invoke the module with
    Map<String, String> argumentValues = null;

    //store the commandlinemodule that we need to invoke
    protected CommandLineModule commandLineModule = null;

    /**
     * Each subclass should implement a method to return an instance of the commandlinemodule it needs to invoke
     * @return
     */
    protected abstract CommandLineModule createCommandLineModule();


    public BaseCommandLineTest()
    {

    }

    /**
     * Set an argument value in the argumentValues map, overriding what's there
     * @param argumentName
     * @param value
     */
    protected void setArgumentValue(String argumentName, String value)
    {
        if (argumentValues == null)
            clearArgumentValues();
        argumentValues.put(argumentName.toLowerCase(), value);
    }

    /**
     * Unet an argument value in the argumentValues map
     * @param argumentName
     */
    protected void unsetArgumentValue(String argumentName)
    {
        if (argumentValues == null)
            return;
        argumentValues.remove(argumentName.toLowerCase());
    }

    /**
     * Set the value of the unnamed argument
     * @param value
     */
    protected void setUnnamedArgumentValue(String value)
    {
        setArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,
                         value);
    }

    /**
     * Set the value of the unnamed argument
     * @param value
     */
    protected void setUnnamedSeriesArgumentValue(String value)
    {
        setArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                         value);
    }


    /**
     * Unet the unnamed argument value in the argumentValues map
     */
    protected void unsetUnnamedArgumentValue()
    {
        unsetArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }       

    /**
     * Blow away all argument values
     */
    protected void clearArgumentValues()
    {
        argumentValues = new HashMap<String,String>();
    }

    /**
     * Get a reference to the commandline module, creating if necessary
     * @return
     */
    protected CommandLineModule getCommandLineModule()
    {
        if (commandLineModule == null)
            renewCommandLineModule();
        return commandLineModule;
    }

    /**
     * Create a new instance of the commandline module
     */
    protected void renewCommandLineModule()
    {
        commandLineModule = createCommandLineModule();
    }

    /**
     * Execute the commandline module and add any exceptions to the test result (convenience method)
     * @param result
     * @return true iff there were no exceptions
     */
    protected boolean executeCommandLineModule(TestResult result)
    {
        try
        {
            getCommandLineModule().invoke(argumentValues);
        }
        catch (Exception e)
        {
            result.addError(this, e);
            return false;
        }
        return true;
    }
}
