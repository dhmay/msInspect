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

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.ViewerUserManualGenerator;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.text.html.HTMLFrameHyperlinkEvent;
import javax.swing.text.html.HTMLDocument;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.HyperlinkEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.IOException;


/**
 * Create and display a user manual
 */
public class UserManualCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(UserManualCLM.class);

    protected File outFile = null;

    protected boolean useGUI = true;

    protected CommandLineModule moduleToDocument = null;

    public UserManualCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "usermanual";

        mHelpMessage ="Automatically create an HTML user manual for msInspect commandline functions, using the " +
                "self-documentation features of the commandline modules themselves.";
        mShortDescription = "Create a user manual for commandline functions";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToWriteArgumentDefinition("out", false, "Output file"),
                        new BooleanArgumentDefinition("gui", false,
                                "Open in a GUI window? (otherwise HTML output is simply written to a file or to " +
                                        "standard output", useGUI),
                        new StringArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,
                                false,
                                "Single command to document (leave blank for full manual)"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        outFile = getFileArgumentValue("out");
        useGUI = getBooleanArgumentValue("gui");

        if (this.hasUnnamedArgumentValue())
        {
            String command = getUnnamedArgumentValue().toString();
            try
            {
                moduleToDocument = ViewerCommandLineModuleDiscoverer.getSingletonInstance().getCommandLineModule(command);
            }
            catch (FileNotFoundException e)
            {
                throw new ArgumentValidationException("Unknown command " + command);
            }
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PrintWriter outPW = null;


        if (useGUI && outFile == null)
            outFile = TempFileManager.createTempFile("help_" + getCommandName() + ".html", this);


        if (outFile == null)
        {
            outPW = new PrintWriter(System.out);
        }
        else
        {
            try
            {
                outPW = new PrintWriter(outFile);
            }
            catch (FileNotFoundException e)
            {
                throw new CommandLineModuleExecutionException("Can't write file " + outFile.getAbsolutePath());
            }
        }

        try
        {

            if (moduleToDocument != null)
            {
                new ViewerUserManualGenerator().generateCommandManualEntry(moduleToDocument, outPW);
                outPW.flush();
            }
            else
            {
                new ViewerUserManualGenerator().generateFullManual(outPW);
                outPW.flush();

            }
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to write output document", e);
        }

        if (useGUI)
        {
            try
            {
                HtmlViewerPanel.showFileInDialog(outFile, "Manual");
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to open browser, HTML is in " +
                        outFile.getAbsolutePath(), e);
            }
        }
    }

     class MyHyperLinkListener implements HyperlinkListener
     {
         public void hyperlinkUpdate(HyperlinkEvent e) {
             if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                 JEditorPane pane = (JEditorPane) e.getSource();
                 if (e instanceof HTMLFrameHyperlinkEvent) {
                     HTMLFrameHyperlinkEvent evt = (HTMLFrameHyperlinkEvent)e;
                     HTMLDocument doc = (HTMLDocument)pane.getDocument();
                     doc.processHTMLFrameHyperlinkEvent(evt);
                 } else {
                     try {
                         pane.setPage(e.getURL());
                     } catch (Throwable t) {
                         t.printStackTrace();
                     }
                 }
             }
         }
     }

}
