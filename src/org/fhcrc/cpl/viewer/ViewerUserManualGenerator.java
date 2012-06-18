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
import org.fhcrc.cpl.toolbox.commandline.CLMUserManualGenerator;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleDiscoverer;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandLineModuleDiscoverer;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;


/**
 * User Manual generation specific to msInspect
 */
public class ViewerUserManualGenerator extends CLMUserManualGenerator
{
    protected static Logger _log = Logger.getLogger(ViewerUserManualGenerator.class);

    public void writeIntro(Writer outW)
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
                "specified.  Some commands have 'advanced' arguments, which are never required; it is not " +
                "recommended to specify values for those arguments unless you know exactly what you're doing.\n<p>\n");
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

    public CommandLineModuleDiscoverer getCommandLineModuleDiscoverer()
    {
        return new ViewerCommandLineModuleDiscoverer();
    }
}
