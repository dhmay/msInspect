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

import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.BooleanArgumentDefinition;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

/**
 * Utility class for commandline modules
 */
public class CommandLineModuleUtilities
{
    public static Logger _log = Logger.getLogger(CommandLineModuleUtilities.class);

    public static String createOutputFilename(String inputFileName, String newSuffix)
    {
        String outputPrefix = null;
        if (inputFileName.contains("."))
            outputPrefix = inputFileName.substring(0, inputFileName.indexOf("."));
        else
            outputPrefix = inputFileName;

        return outputPrefix + "." + newSuffix;
    }

    public static File createOutputFile(File inputFile, String newSuffix, File outputDir)
    {
        String inputFileName = inputFile.getName();
        return new File(outputDir, createOutputFilename(inputFileName, newSuffix));
    }

    public static File findFileLikeFile(File inputFile, File directory, String suffix)
            throws FileNotFoundException
    {
        String inputFileName = inputFile.getName();
        String inputFilePrefix = null;
        try
        {
            inputFilePrefix = inputFileName.substring(0, inputFileName.indexOf("."));
            _log.debug("Searching directory '" + directory.getName() + "' for file with prefix '" + inputFilePrefix +
                    "' and suffix '" + suffix + "'");
        }
        catch (StringIndexOutOfBoundsException e)
        {
            throw new FileNotFoundException("File with no extension cannot be matched with another file: " + inputFileName);
        }

        File result = findFileWithPrefix(inputFilePrefix, directory, suffix);

        if (result == null)
            throw new FileNotFoundException("Failed to find matching file for input file " +
                    inputFile.getAbsolutePath() + " in directory " + directory.getAbsolutePath());
        return result;
    }

    public static File findFileWithPrefix(String prefix, File directory, String suffix)
            throws FileNotFoundException
    {
        File result = null;
        for (File potentialFile : directory.listFiles())
        {
            String potentialFileName = potentialFile.getName();
            if (suffix != null && !potentialFileName.endsWith(suffix))
                continue;
            if (potentialFileName.startsWith(prefix))
            {
                result = potentialFile;
                break;
            }
        }

        if (result == null)
            throw new FileNotFoundException("Failed to find matching file for prefix " +
                    prefix + " in directory " + directory.getAbsolutePath());
        return result;
    }

    /**
     *
     * @param module
     * @param exception
     * @param isLogEnabled
     * @param logFile
     * @return
     * @throws IOException
     */
    public static File createFailureReport(CommandLineModule module, Exception exception,
                                           boolean isLogEnabled, File logFile, String headerText)
            throws IOException
    {
        File reportFile = TempFileManager.createTempFile("failure_" + module.getCommandName() + ".txt", module);
        PrintWriter pw = new PrintWriter(reportFile);
        pw.println(headerText);
        pw.println("revision=" + ApplicationContext.getProperty("REVISION"));
        DateFormat dateFormat =
            new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
        pw.println("date=" + dateFormat.format(new Date()));
        pw.println("java_version=" + String.valueOf(System.getProperty("java.version")));
        pw.println("free_heap=" + (Runtime.getRuntime().freeMemory() / 1024 / 1024) + "MB");
        pw.println("total_current_heap=" + (Runtime.getRuntime().totalMemory() / 1024 / 1024) + "MB");
        pw.println("max_heap=" + (Runtime.getRuntime().maxMemory() / 1024 / 1024) + "MB");

        pw.println("OS=" + String.valueOf(System.getProperty("os.name")));
        pw.println("command=" + module.getCommandName());
        pw.println("arguments:");
        pw.println("**********");
        pw.println(createFailureReportEntry(module.getCommandName(), module.getArgumentValueStrings()));
        pw.println("**********");
        pw.flush();

        pw.println("exception_message=" + exception.getMessage());
        if (exception instanceof CommandLineModuleExecutionException)
        {
            Exception nestedException = ((CommandLineModuleExecutionException) exception).getNestedException();
            pw.println("nested_exception_type=" + (nestedException == null? "null" : nestedException.getClass().getName()));
        }
        pw.println("stack_trace:");
        exception.printStackTrace(pw);

        pw.flush();

        if (isLogEnabled)
        {
            pw.println("**********\nsession_log\n**********");            

            FileReader logReader = new FileReader(logFile);
            BufferedReader br = new BufferedReader(logReader);

            String line = null;
            while ((line = br.readLine()) != null)
            {
                pw.println(line);
                pw.flush();
            }
        }


        pw.close();

        return reportFile;
    }

    /**
     * Create an entry in the failure report for an invocation of a command
     * @param commandName
     * @param argumentMap
     * @return
     */
    public static String createFailureReportEntry(String commandName,
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

    public static String createFailureReportAndPrompt(CommandLineModule module,
                                                      Exception exception, boolean logEnabled, File logFile,
                                                      String additionalErrorText, String failureReportHeaderText)
    {
        String result = null;
        try
        {
            File failureReportFile = CommandLineModuleUtilities.createFailureReport(module, exception,
                    logEnabled, logFile, failureReportHeaderText);
            result = "\nA failure report has been generated in the file:\n  " +
                    failureReportFile.getAbsolutePath() + "\n" + additionalErrorText + "\n";

            if (!logEnabled)
                result += "For more detailed logging of this failure, re-run your command from the command line \n" +
                          "with the '--log' option.";
        }
        catch (IOException e)
        {
            result = "Failed to create failure report!";
        }
        return result;
    }

    public static boolean isUnnamedSeriesArgument(CommandLineArgumentDefinition argDef)
    {
        return argDef.getArgumentName().equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
    }

    public static boolean isUnnamedArgument(CommandLineArgumentDefinition argDef)
    {
        return argDef.getArgumentName().equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    /**
     * Does this set of argument definitions allow the unnamed (single) argument?
     * @return
     */
    public static CommandLineArgumentDefinition getUnnamedArgumentDefinition(CommandLineModule module)
    {
        return module.getArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    public static CommandLineArgumentDefinition getUnnamedArgumentSeriesDefinition(CommandLineModule module)
    {
        return module.getArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
    }


    /**
     * Parse the raw arguments passed in by the user, in the context of the supplied module.  Return a
     * map from argument name to argument value.  If undefined arguments are supplied, throw an
     * IllegalArgumentException
     * @param module
     * @param args
     * @return
     * @throws IllegalArgumentException if undefined arguments are supplied
     */
    public static Map<String,String> parseRawArguments(CommandLineModule module, String[] args)
            throws IllegalArgumentException
    {
        HashMap<String, String> argNameValueMap = new HashMap<String, String>();
        boolean alreadyProcessedNakedArgument = false;
        List<String> unnamedArgumentSeriesList = new ArrayList<String>();

        for (int i = 0; i < args.length; i++)
        {
            _log.debug("Processing argument: " + args[i]);
            boolean ignoreThisArg = false;
            String[] param = args[i].split("=");

            //if there's something really funky with the argument, punt
            if (param.length < 1 || param.length > 2)
            {
                throw new IllegalArgumentException("Unable to parse argument " + args[i]);
            }
            if (param.length == 1)
            {
                //an argument with no = sign.  First check for any BooleanArgumentDefinitions
                //that match the name
                if (param[0].startsWith("--"))
                {
                    CommandLineArgumentDefinition argDef =
                            module.getArgumentDefinition(param[0].substring(2).toLowerCase());
                    if (argDef != null)
                    {
                        //if it's a boolean argument definition, set value to true.  Otherwise, bad
                        if (argDef instanceof BooleanArgumentDefinition)
                        {
                            param = new String[2];
                            param[0] = argDef.getArgumentName();
                            param[1] = "true";
                        }
                        else
                        {
                            //make translatable?
                            ApplicationContext.infoMessage("No argument value specified for argument " + argDef.getArgumentName());
                        }
                    }
                }
                // Check for short args (e.g. "-m0.75")
                else if (param[0].startsWith("-"))
                {
                    String paramName = param[0].substring(1,2);
                    String paramVal = param[0].substring(2);

                    // Command-line implementation explicitly ignores case. This can be deadly for
                    // short args, so currently disallow all upper-case args.
                    if (!paramName.equalsIgnoreCase(paramName))
                    {
                        ApplicationContext.infoMessage("Upper case short-form arguments are disallowed; found '" +
                                paramName + "'");
                    }

                    if ("".equals(paramVal))
                        paramVal = checkBooleanDefault(module, paramName);

                    if (null != paramVal)
                    {
                        param = new String[2];
                        param[0] = paramName;
                        param[1] = paramVal;
                    }
                }
                //Not an unnamed BooleanArgumentDefinition... try the naked argument
                else
                {
                    //if we've got an argument with no = sign, then check to see if this
                    //command allows that.  If so, handle it that way.
                    if (CommandLineModuleUtilities.getUnnamedArgumentDefinition(module) != null)
                    {

                        //two naked arguments, die
                        if (alreadyProcessedNakedArgument)
                        {
                            throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                                    args[i]));
                        }
                        param = new String[2];
                        param[0] = CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT;
                        param[1] = args[i];
                        alreadyProcessedNakedArgument = true;
                    }
                    //ok, maybe this module allows a series of unnamed arguments....
                    else if (CommandLineModuleUtilities.getUnnamedArgumentSeriesDefinition(module) != null)
                    {
                        unnamedArgumentSeriesList.add(args[i]);
                        ignoreThisArg = true;
                    }
                    else
                    {
                        throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                                args[i]));
                    }
                }
            }
            if (ignoreThisArg)
                continue;

            String paramName = param[0];
            if (paramName.startsWith("--"))
                paramName = paramName.substring(2);
            String paramVal = null;
            if (param.length > 1)
                paramVal = param[1];

            _log.debug("Got arg " + paramName + " = '" + paramVal + "'");

            if (module.getArgumentDefinition(paramName) != null)
            {
                if (argNameValueMap.containsKey(paramName.toLowerCase()))
                {
                    throw new IllegalArgumentException("Argument " + paramName +
                            " is specified more than once.  Quitting.\n");
                }
                argNameValueMap.put(paramName.toLowerCase(), paramVal);
            }
            else
            {
                throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                        paramName));
            }
        } // End of arg loop

        if (CommandLineModuleUtilities.getUnnamedArgumentSeriesDefinition(module) != null &&
                !unnamedArgumentSeriesList.isEmpty())
        {
            StringBuffer combinedUnnamedArguments = new StringBuffer();
            boolean firstArg=true;
            for (String unnamedArgumentValue : unnamedArgumentSeriesList)
            {
                //use '*ARGSEP*' as arg separator
                if (!firstArg)
                    combinedUnnamedArguments.append(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                combinedUnnamedArguments.append(unnamedArgumentValue);
                firstArg = false;
            }
            argNameValueMap.put(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                    combinedUnnamedArguments.toString());
        }
        return argNameValueMap;
    }

    /**
     * Return the default value of a Boolean argument as a String
     * @param module
     * @param paramName
     * @return
     */
    private static String checkBooleanDefault(CommandLineModule module, String paramName)
    {
        CommandLineArgumentDefinition argDef =
            module.getArgumentDefinition(paramName);

        // Argument not defined
        if (argDef == null)
            return null;

        // Not a boolean argument
        if (!(argDef instanceof BooleanArgumentDefinition))
        {
            ApplicationContext.infoMessage("No argument value specified for argument " + argDef.getArgumentName());
            return null;
        }

        return "true";
    }
}
