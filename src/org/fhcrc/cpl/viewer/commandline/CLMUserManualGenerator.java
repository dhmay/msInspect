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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;


/**
 * test
 */
public class CLMUserManualGenerator
{
    protected static Logger _log = Logger.getLogger(CLMUserManualGenerator.class);

    /**
     * generate a full user manual
     */
    public static void generateFullManual(Writer outW)
            throws IOException
    {
        outW.write("<html>\n<body>\n");
        outW.flush();

        ApplicationContext.setMessage("Creating full HTML manual...");

        writeIntro(outW);
        outW.write("</hr>");
        writeTOCAndManualBody(outW);

        outW.write("</body>\n</html>\n");
        outW.flush();

        ApplicationContext.setMessage("Done creating manual.");
    }

    /**
     * Write a table of contents
     * @param outW
     * @throws IOException
     */
    public static void writeTOCAndManualBody(Writer outW)
            throws IOException
    {
        Map<String, Map<String, CommandLineModule>> packageNameModuleMap =
                CommandLineModuleDiscoverer.findAllCommandLineModulesByPackage();

        //Write out the TOC entries immediately, and buffer the manual body
        //TODO: if this gets too huge, would be better to write it to a file and read it back
        StringBuffer manualBodyBuf = new StringBuffer();



        //arrange known package names in order
        List<String> packageNameList = new ArrayList<String>();
        for (int i=0; i<CommandLineModuleDiscoverer.modulePackageNames.length; i++)
        {
            String packageName = CommandLineModuleDiscoverer.modulePackageNames[i];
            if (packageNameModuleMap.containsKey(packageName))
                packageNameList.add(packageName);
        }
        //add any unknown packages
        for (String packageName : packageNameModuleMap.keySet())
        {
            if (!packageNameList.contains(packageName))
                packageNameList.add(packageName);
        }

        outW.write("<h2>Table of Contents</h2>The available commands are divided into sections: ");
        for (int i=0; i<packageNameList.size(); i++)
        {
            if (i>0) outW.write(", ");
            outW.write("<b>" + CommandLineModuleDiscoverer.getPackageShortDescription(packageNameList.get(i)) + "</b>");
        }
        outW.write("<table border=\"1\">");

        for (String packageName : packageNameList)
        {
            String packageidentifier = CommandLineModuleDiscoverer.getPackageIdentifier(packageName);
            
            String packageShortDescription =  CommandLineModuleDiscoverer.getPackageShortDescription(packageName);
            String packageLongDescription =  CommandLineModuleDiscoverer.getPackageLongDescription(packageName);

            outW.write("<tr><td><a href=\"#package_" + packageidentifier + "\"><h3>" + packageShortDescription + "</h3></a>");
            outW.write(packageLongDescription + "<p>");
            outW.write("</td></tr><tr><td><table border=\"1\">\n<tr><th>Command</th><th>Summary</th></tr>\n");

            Map<String, CommandLineModule> nameModuleMap = packageNameModuleMap.get(packageName);
            List<String> commandList = new ArrayList<String>(nameModuleMap.keySet());
            Collections.sort(commandList);

            manualBodyBuf.append("<hr><a name=\"package_" + packageidentifier + "\"></a><h2>" +
                    packageShortDescription + "</h2><h4>" + packageLongDescription + "</h4>");

            //print table of contents entries for section
            for (String moduleName : commandList)
            {
                CommandLineModule module = nameModuleMap.get(moduleName);
                String shortDescription = nameModuleMap.get(moduleName).getShortDescription();
                if (shortDescription == null || shortDescription.length() < 1)
                    shortDescription = "(no description)";
                outW.write("<tr><td><strong><a href=\"#" + moduleName + "\">" + moduleName + "</a></strong></td><td>" +
                        shortDescription + "</td></tr>\n");

                manualBodyBuf.append("<p/><hr>\n" + module.getHtmlHelpFragment() + "\n");
            }

            //close TOC section
            outW.write("</table></td></tr>\n\n");
            outW.flush();
        }

        //close table of contents
        outW.write("</table>\n<p><hr>\n");
        outW.flush();

        //write manual body
        outW.write(manualBodyBuf.toString());
        outW.flush();
    }

    public static void writeIntro(Writer outW)
            throws IOException
    {
        outW.write("<h1>msInspect Commandline Functions Manual</h1>\n<p>\n");
        outW.write("<h2>Introduction</h2>\n<p>\n");                   
        outW.write("This manual describes all of the commandline functionality available in msInspect.  " +
                "To access any of these commands, type the java command that you use to start msInspect " +
                "(e.g., 'java -Xmx1G -jar viewerApp.jar'), followed by ' --&lt;command&gt;' (the name of the " +
                "command), followed by any arguments.  Argument names are not case-sensitive.\n<p>\n");
        outW.write("To access a graphical interface for entering arguments for any command (with graphical " +
                "help for choosing files, etc.), use the command " +
                "'--interactive', followed by the name of the command you wish to execute.\n<p>\n");
        outW.write("All commandline functions may also be accessed through the msInspect graphical user " +
                "interface using the 'Run Command' menu item under the 'File' menu.\n<p>");
        outW.write("For more information about msInspect, as well as help with the graphical user interface, " +
                "go to the official msInspect webpage at " +
                "<a href=\"http://proteomics.fhcrc.org/CPL/msinspect.html\">" +
                "http://proteomics.fhcrc.org/CPL/msinspect.html</a>\n<p>\n");
        outW.write("Required arguments for each command are indicated in <b>bold text</b> in the arguments " +
                "table. In some cases, additional arguments may be required if certain arguments are " +
                "specified.\n<p>\n");
        outW.write("In addition to the arguments discussed below for each particular command, msInspect accepts " +
                "two special arguments for all commands:\n");
        outW.write("<ol><li><strong>--log</strong> : Turn on logging of all output messages to a log file " +
                "(log file location can be specified with --log=&lt;filepath&gt;)</li>" +
                "<li><strong>--debug</strong> : Turn on full debug logging of all Java classes (individual class " +
                "names can also be specified with --debug=&lt;full_class_name&gt;,&lt;full_class_name&gt;..." +
                "</li></ol>\n<p>\n");
        outW.write("<i>This document automatically generated on " +
                new SimpleDateFormat("MMMM d, yyyy").format(new Date()) +
                " by msInspect revision " + ApplicationContext.getProperty("REVISION") +
                ".  If you are running a newer version of msInspect, you may generate the current version of " +
                "this document with the 'usermanual' command.</i>\n");
    }

    /**
     * Generate manual entry for a single command, in its own doc
     * @param moduleToDocument
     * @param outW
     * @throws IOException
     */
    public static void generateCommandManualEntry(CommandLineModule moduleToDocument, Writer outW)
            throws IOException
    {
        outW.write("<html>\n<body>\n");
        outW.flush();

        outW.write(moduleToDocument.getHtmlHelpFragment());
        outW.flush();

        outW.write("</body>\n</html>\n");
        outW.flush();
    }

}
