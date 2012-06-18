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
package org.fhcrc.cpl.viewer;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandLineModuleDiscoverer;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Map;
import java.util.HashMap;

/**
 * Provides macro-running functionality in MSInspect.  Knows how to read commands and
 * arguments form a macro file in the following format:
 * -Comment lines begin with # and are ignored
 * -Commands begin immediately on a given line.  Anything other than a comment that begins
 * immediately at the start of a line is considered a command
 * -Arguments start with a tab character.  Anything other than a comment that begins with
 * a tab is considered an argument.  Arguments are of the format NAME=VALUE.  For this reason,
 * '=' is not currently allowed in arguments.  For unnamed arguments, simply us a tab
 * followed by the argument value
 * -Arguments for a given command should follow directly after that command
 *
 * All handling of individual macro commands is spelled out here.  In general it's desirable
 * to have these macros follow the same code path as the UI, to the greatest extent possible
 */
public class CommandFileRunner
{
    private static Logger _log = Logger.getLogger(CommandFileRunner.class);

    /**
     * Parse a macro file.  For each command, with arguments, invoke one of the known
     * command handlers if there is one for this command
     * @param macroFile
     */
    public static void runMacroFile(File macroFile)
    {
        if (macroFile == null || !macroFile.exists())
        {
            ApplicationContext.infoMessage(TextProvider.getText(
                    "FILE_FILE_DOES_NOT_EXIST",
                    macroFile.getAbsolutePath()));
            return;
        }

        BufferedReader reader = null;
        try
        {
            reader = new BufferedReader(new FileReader(macroFile));
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }

        String command = null;

        Map<String, CommandLineModule> moduleMap =
                ViewerCommandLineModuleDiscoverer.getSingletonInstance().findAllCommandLineModules();


        //all macro command methods return true if success, false if failure.
        //If failure, stop executing further commands
        boolean keepGoing = true;
        while (keepGoing && (command = getNextCommand(reader)) != null)
        {
            if (command.toLowerCase().startsWith("echo") &&
                    (command.length() == "echo".length() ||
                Character.isWhitespace(command.charAt("echo".length()))))
            {
                System.err.println(command.substring(Math.min("echo".length() + 1,
                                                              command.length())));
                continue;
            }

            command = command.toLowerCase();

            CommandLineModule module = moduleMap.get(command);
            if (module == null)
                throw new RuntimeException("Unknown command " + command);
            Map<String,String> argNameValueMap = getArguments(reader, module);
            try
            {
                module.digestArguments(argNameValueMap);
            }
            catch (ArgumentValidationException e)
            {
                ApplicationContext.infoMessage(TextProvider.getText("FAILED_ARGUMENT_VALIDATION") + "\n" +
                        TextProvider.getText("ERROR") + ": " + e.getMessage());
                return;
            }


            try
            {
                module.execute();
                ApplicationContext.infoMessage(TextProvider.getText("COMMAND_COMPLETE", module.getCommandName()));
            }
            catch (CommandLineModuleExecutionException e)
            {
                ApplicationContext.errorMessage(TextProvider.getText("ERROR_RUNNING_COMMAND_COMMAND", module.getCommandName()),e);
            }

            //short sleep just to space out commands, to allow for cleanup
            try
            {
                Thread.sleep(500);
            }
            catch (Exception e) {}
        }
    }

    /**
     * skips past any argument lines to the next command, in my funky file format
     * TODO: this should probably be reimplemented as an Iterator
     * @param reader
     * @return
     */
    public static String getNextCommand(BufferedReader reader)
    {
        String result = null;
        try{
            String line = null;
            while ((line = reader.readLine()) != null)
            {
                if (line.length() > 0 && !line.startsWith("\t") && !line.startsWith("#"))
                {
                    result = line.trim();
                    break;
                }
            }
        }
        catch (Exception e) { e.printStackTrace(System.err); }
        return result;
    }

    /**
     * Gobbles a list of arguments, stopping when it hits eof or
     * a protein definition.  If it hits a protein definition, it backs up
     * @param reader
     * @return
     */
    public static Map<String,String> getArguments(BufferedReader reader,
                                                  CommandLineModule module)
    {
        Map<String,String> result = new HashMap<String,String>();
        try
        {
            String line = null;
            reader.mark(2000);
            while ((line = reader.readLine()) != null)
            {
                if (line.startsWith("\t"))
                {
                    Pair<String,String> argument = processArgument(line, module);
                    if (argument != null)
                        result.put(argument.first, argument.second);
                    reader.mark(2000);
                }
                else if (line.startsWith("#"))
                {
                    continue;
                }
                else
                {
                    //follow the mark back to after the last command block
                    reader.reset();
                    break;
                }
            }
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        return result;
    }

    /**
     * split an argument around an = sign
     * TODO: add error handling
     * @param argument
     * @return the argument name and value, or null if splitting around = gives something
     * other than two strings
     */
    protected static Pair<String,String> processArgument(String argument,
                                                         CommandLineModule module)
    {
        argument = argument.trim();
        String[] param = argument.split("=");
        Pair<String,String> result = new Pair<String,String>("","");
        if (param.length == 1)
        {
            if (module.getArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT) != null)
                result.first = CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT;
            else
                result.first = CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT;
            result.second = param[0];
        }
        else if (param.length == 2)
        {
            result.first = param[0];
            result.second = param[1];
        }
        else
            return null;

        if (result.first.equals(
                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
        {
            _log.debug("Unnamed series.  Before:\n" + result.second);
            result.second = result.second.trim();
            result.second = result.second.replaceAll(" ", CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
            _log.debug("After:\n" + result.second);                    

        }

        return result;
    }

    /**
     * print out a generic error message
     */
    protected static void macroError()
    {
        ApplicationContext.infoMessage("FAILED_TO_RUN_MACRO");
    }

    /**
     * print out a missing-arguments error message
     */
    protected static void macroMissingArguments()
    {
        ApplicationContext.infoMessage(TextProvider.getText("MACRO_COMMAND_MISSING_ARGUMENTS"));
    }

    /**
     * Utility method to create a command, for use in a command file,
     * based on a command name and argument map
     * @param commandName
     * @param argumentMap
     */
    public static String createCommandFileEntry(String commandName,
                                                Map<String,String> argumentMap)
    {
        StringBuffer resultBuf = new StringBuffer();
        resultBuf.append(commandName + "\n");
        for (String argumentName : argumentMap.keySet())
        {
            resultBuf.append("\t");
            if (argumentName.equalsIgnoreCase(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT)
                || argumentName.equalsIgnoreCase(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT))
            {
                //do not append anything for the argument name or equals sign
            }
            else
            {
                resultBuf.append(argumentName + "=");
            }
            resultBuf.append(argumentMap.get(argumentName) + "\n");
        }
        return resultBuf.toString();
    }



    /**
     * action for the menu item that kicks macro-running off
     */
    public static class CommandFileRunnerAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            WorkbenchFileChooser chooser = new WorkbenchFileChooser();
            int chooserStatus = chooser.showOpenDialog(ApplicationContext.getFrame());
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;            
            File file = chooser.getSelectedFile();

            runMacroFile(file);
        }

    }


}
