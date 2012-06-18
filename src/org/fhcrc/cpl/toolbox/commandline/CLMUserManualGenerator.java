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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleDiscoverer;
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
    public void generateFullManual(Writer outW)
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

    public CommandLineModuleDiscoverer getCommandLineModuleDiscoverer()
    {
        return new CommandLineModuleDiscoverer();
    }

    /**
     * Write a table of contents
     * @param outW
     * @throws IOException
     */
    public void writeTOCAndManualBody(Writer outW)
            throws IOException
    {
        Map<String, Map<String, CommandLineModule>> packageNameModuleMap =
                getCommandLineModuleDiscoverer().findAllCommandLineModulesByPackage();

        //Write out the TOC entries immediately, and buffer the manual body
        //TODO: if this gets too huge, would be better to write it to a file and read it back
        StringBuffer manualBodyBuf = new StringBuffer();



        //arrange known package names in order
        List<String> packageNameList = new ArrayList<String>();
        for (int i=0; i< getCommandLineModuleDiscoverer().modulePackageNames.length; i++)
        {
            String packageName = getCommandLineModuleDiscoverer().modulePackageNames[i];
            if (packageNameModuleMap.containsKey(packageName))
                packageNameList.add(packageName);
        }
        //add any unknown packages
        for (String packageName : packageNameModuleMap.keySet())
        {
            if (!packageNameList.contains(packageName))
                packageNameList.add(packageName);
        }

        if (packageNameList.size() > 1)
        {
            outW.write("<h2>Table of Contents</h2>The available commands are divided into sections: ");
            for (int i=0; i<packageNameList.size(); i++)
            {
                if (i>0) outW.write(", ");
                outW.write("<b>" + getCommandLineModuleDiscoverer().getPackageShortDescription(packageNameList.get(i)) + "</b>");
            }
            outW.write("<table border=\"1\">");
        }

        for (String packageName : packageNameList)
        {
            String packageidentifier = getCommandLineModuleDiscoverer().getPackageIdentifier(packageName);
            
            String packageShortDescription =  getCommandLineModuleDiscoverer().getPackageShortDescription(packageName);
            String packageLongDescription =  getCommandLineModuleDiscoverer().getPackageLongDescription(packageName);

            if (packageNameList.size() > 1)
            {
                outW.write("<tr><td><a href=\"#package_" + packageidentifier + "\"><h3>" + packageShortDescription + "</h3></a>");
                outW.write(packageLongDescription + "<p>");
                outW.write("</td></tr>");
            }
            outW.write("<tr><td><table border=\"1\">\n<tr><th>Command</th><th>Summary</th></tr>\n");

            Map<String, CommandLineModule> nameModuleMap = packageNameModuleMap.get(packageName);
            List<String> commandList = new ArrayList<String>(nameModuleMap.keySet());
            Collections.sort(commandList);

            if (packageNameList.size() > 1)
            {
                manualBodyBuf.append("<hr><a name=\"package_" + packageidentifier + "\"></a><h2>" +
                        packageShortDescription + "</h2><h4>" + packageLongDescription + "</h4>");
            }

            //print table of contents entries for section
            for (String moduleName : commandList)
            {
                CommandLineModule module = nameModuleMap.get(moduleName);
                String shortDescription = nameModuleMap.get(moduleName).getShortDescription();
                if (shortDescription == null || shortDescription.length() < 1)
                    shortDescription = "(no description)";
                outW.write("<tr><td><strong><a href=\"#" + moduleName + "\">" + moduleName + "</a></strong></td><td>" +
                        shortDescription + "</td></tr>\n");

                manualBodyBuf.append("<p><hr>\n" + module.getHtmlHelpFragment() + "\n");
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

    /**
     * Write a manual introduction.  This should be overridden by a class specific to the application
     * @param outW
     * @throws IOException
     */
    public void writeIntro(Writer outW)
            throws IOException
    {
        outW.write("<h1>Commandline Functions Manual</h1>\n<p>\n");
        outW.write("<h2>Introduction</h2>\n<p>\n");                   
        outW.write("This manual describes all of the commandline functionality available.  " +
                "To access any of these commands, type the java command that you use to start the application, " +
                "followed by ' --&lt;command&gt;' (the name of the " +
                "command), followed by any arguments.  Argument names are not case-sensitive.\n<p>\n");
        outW.write("Required arguments for each command are indicated in <b>bold text</b> in the arguments " +
                "table. In some cases, additional arguments may be required if certain arguments are " +
                "specified.  Some commands may have 'advanced' arguments, which are never required; it is not " +
                "recommended to specify values for those arguments unless you know exactly what you're doing.\n<p>\n");
        outW.write("In addition to the arguments discussed below for each particular command, there are " +
                "two special arguments for all commands:\n");
        outW.write("<ol><li><strong>--log</strong> : Turn on logging of all output messages to a log file " +
                "(log file location can be specified with --log=&lt;filepath&gt;)</li>" +
                "<li><strong>--debug</strong> : Turn on full debug logging of all Java classes (individual class " +
                "names can also be specified with --debug=&lt;full_class_name&gt;,&lt;full_class_name&gt;..." +
                "</li></ol>\n<p>\n");
        outW.write("<i>This document automatically generated on " +
                new SimpleDateFormat("MMMM d, yyyy").format(new Date()) + ".");
    }

    /**
     * Generate manual entry for a single command, in its own doc
     * @param moduleToDocument
     * @param outW
     * @throws IOException
     */
    public void generateCommandManualEntry(CommandLineModule moduleToDocument, Writer outW)
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
