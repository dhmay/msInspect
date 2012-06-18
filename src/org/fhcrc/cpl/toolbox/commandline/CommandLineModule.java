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

package org.fhcrc.cpl.toolbox.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;

import java.util.Map;

/**
 * All command line modules must implement this interface in order to extend msInspect commandline functionality.
 * The flow is as follows:
 *  -Application finds all CommandLineModule-implementing classes within a specified list of packages
 *  -For each such class, it calls getCommandName() to determine the command that should be used to invoke
 *  the module.  It creates an instance of the module using the no-arg constructor
 *  -If the user enters this command, getArgumentDefinitions() is called to find out what arguments are allowed
 *  -The user's arguments are parsed.  If they all pass basic validation, assignArgumentValues() is called
 *  -If assignArgumentValues() succeeds, execute() is called
 *
 * Although it is not required, we recommend that CommandLineModule classes extend BaseCommandLineModuleImpl,
 * which provides some convenience method implementations.
 */
public interface CommandLineModule
{
    //sentinel value to tell Application that this commandline module has no explicit usage message
    //defined, so Application should use the ArgumentDefinitions defined for the modules to build
    //an appropriate message
    public static final String MODULE_USAGE_AUTOMATIC = "MODULE_USAGE_AUTOMATIC";
    //sentinel value to tell Application that this commandline module has no explicit help message
    //defined, so Application should use the ArgumentDefinitions defined for the modules to build
    //an appropriate message
    public static final String MODULE_HELP_AUTOMATIC = "MODULE_HELP_AUTOMATIC";

    //String to separate unnamed arguments in a series
    public static final String UNNAMED_ARG_SERIES_SEPARATOR = "AAAAARGSEPAAAA";

    /**
     *
     * @return a String telling the user how to invoke this command
     */
    public String getUsage();

    /**
     * a String giving the user detailed help on invoking this command
     * @return
     */
    public String getHelpMessage();

    /**
     * Returns the full help message to be displayed by --help <command>
     * @return
     */
    public String getFullHelp();

    /**
     * Returns an HTML fragment containing full help information for this module
     * @return
     */
    public String getHtmlHelpFragment();

    /**
     *
     * @return a String that is the command that a user should type in order to invoke this
     * CommandLineModule.  This String should not collide with the Command Name of another module
     */
    public String getCommandName();

    /**
     *
     * @return a short (one-line, preferably) description of this functionality
     */
    public String getShortDescription();

    /**
     *
     * @return an array of CommandLineArgumentDefinitions.  This array tells msInspect
     * which required and
     * optional parameters are used by your module.  No particular order
     */
    public CommandLineArgumentDefinition[] getArgumentDefinitions();

    /**
     *
     * @return an array of CommandLineArgumentDefinitions.  This array tells msInspect
     * which BASIC required and
     * optional parameters are used by your module.  No particular order
     */
    public CommandLineArgumentDefinition[] getBasicArgumentDefinitions();

    /**
     *
     * @return an array of CommandLineArgumentDefinitions.  This array tells msInspect
     * which ADVANCED required and
     * optional parameters are used by your module.  No particular order
     */
    public CommandLineArgumentDefinition[] getAdvancedArgumentDefinitions();

    /**
     *
     * @return @return an array of CommandLineArgumentDefinitions, sorted
     * however this module wants them to be sorted for display purposes.  May
     * be more computation-intensive than getArgumentDefinitions().
     */
    public CommandLineArgumentDefinition[] getArgumentDefinitionsSortedForDisplay();

    /**
     *
     * @return a CommandLineArgumentDefinition matching the specified name (case-insensitive).
     * Null if not found
     */
    public CommandLineArgumentDefinition getArgumentDefinition(String argumentName);


    
    /**
     *
     * @param argumentValueMap
     */
    public void digestArguments(Map<String, String> argumentValueMap)
            throws ArgumentValidationException;

    /**
     * Return the original argument name-value pairs that were passed into this module via digestArguments()
     * @return
     */
    public Map<String, String> getArgumentValueStrings();    

    /**
     * the first step in invoking your module.  The values assigned to the various arguments by the user
     * are passed to your module for storage and additional validation.  Any communication with the user about
     * their argument values should be done by this method.
     */
    public void assignArgumentValues()
           throws ArgumentValidationException;

    /**
     * Called by msInspect after assignArgumentValues.  Executes your functionality.  Any
     * communication with the user about execution status should be done by this method.
     */
    public void execute() throws CommandLineModuleExecutionException;

    /**
     * Invoke the CommandLineModule, using the argument values specified in
     * the argumentValues array.  Essentially this should just call digestArguments()
     * and then execute()
     * @param argumentValues
     */
    public void invoke(Map<String, String> argumentValues)
            throws ArgumentValidationException, CommandLineModuleExecutionException;

}
