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

package org.fhcrc.cpl.viewer.commandline;

import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.viewer.CommandFileRunner;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import java.io.*;
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

/**
 * Utility class for commandline modules
 */
public class CommandLineModuleUtilities
{
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
        }
        catch (StringIndexOutOfBoundsException e)
        {
            throw new FileNotFoundException("File with no extension cannot be matched with another file: " + inputFileName);
        }
        File result = null;
        for (File potentialFile : directory.listFiles())
        {
            String potentialFileName = potentialFile.getName();
            if (suffix != null && !potentialFileName.endsWith(suffix))
                continue;
            if (potentialFileName.startsWith(inputFilePrefix))
            {
                result = potentialFile;
                break;
            }
        }

        if (result == null)
            throw new FileNotFoundException("Failed to find matching file for input file " +
                    inputFile.getAbsolutePath() + " in directory " + directory.getAbsolutePath());
        return result;
    }

    public static File createFailureReport(CommandLineModule module, Exception exception)
            throws IOException
    {
        File reportFile = TempFileManager.createTempFile("failure_" + module.getCommandName() + ".txt", module);
        PrintWriter pw = new PrintWriter(reportFile);
        pw.println("#msInspect failure report");
        pw.println("#please report this failure and attach this report on the support forum at " +
                   WorkbenchFrame.SUPPORT_URL);
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
        pw.println(CommandFileRunner.createCommandFileEntry(module.getCommandName(), module.getArgumentValueStrings()));
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

        if (Application.isLogEnabled())
        {
            pw.println("**********\nsession_log\n**********");            

            FileReader logReader = new FileReader(Application.getLogFile());
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

    public static String createFailureReportAndPrompt(CommandLineModule module,
                                                      Exception exception)
    {
        String result = null;
        try
        {
            File failureReportFile = CommandLineModuleUtilities.createFailureReport(module, exception);
            result = "\nA failure report has been generated in the file:\n  " +
                    failureReportFile.getAbsolutePath() + "\n" +
                    "If you would like help with this issue, please report the problem in detail\n" +
                    "to the msInspect support forum at:\n  " +
                    WorkbenchFrame.SUPPORT_URL +
                    "\nand attach the failure report file.\n";
            if (!Application.isLogEnabled())
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
}
