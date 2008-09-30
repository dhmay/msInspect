/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.commandline.arguments.ViewerArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.CommandFileRunner;
import org.fhcrc.cpl.viewer.ViewerUserManualGenerator;
import org.fhcrc.cpl.toolbox.commandline.CLMUserManualGenerator;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.toolbox.gui.InteractiveModuleFrame;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.event.*;
import java.util.Map;
import java.util.HashMap;
import java.util.prefs.Preferences;
import java.io.*;

/**
 * Dialog window class for specifying command-line arguments for a given module
 */
public class ViewerInteractiveModuleFrame extends InteractiveModuleFrame
{
    static Logger _log = Logger.getLogger(ViewerInteractiveModuleFrame.class);

    //if true, disable the buttons related to commandfile stuff, and the user's notification
    //that the command is complete
    public boolean guiOnly = false;

    public ViewerInteractiveModuleFrame(CommandLineModule module, Map<String, String> moduleArgMap)
    {
        this(module, false, moduleArgMap);
    }

    /**
     * @param module
     */
    public ViewerInteractiveModuleFrame(CommandLineModule module, boolean guiOnly, Map<String, String> moduleArgMap)
    {
        super(module, moduleArgMap);

        this.guiOnly = guiOnly;


        JButton buttonShowCommand = new JButton(TextProvider.getText("SHOW_COMMAND"));
        JButton buttonSaveCommandFile = new JButton(TextProvider.getText("SAVE_TO_FILE"));
        helper.addListener(buttonShowCommand, "buttonShowCommand_actionPerformed");
        helper.addListener(buttonSaveCommandFile, "buttonSaveCommandFile_actionPerformed");

        if (!guiOnly)
        {
            buttonPanel.add(buttonSaveCommandFile, buttonGBC);
            buttonPanel.add(buttonShowCommand, buttonGBC);
        }        

        prefs = Preferences.userNodeForPackage(ViewerInteractiveModuleFrame.class);
        userManualGenerator = new ViewerUserManualGenerator();
        try
        {
            Image iconImage = ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif"));
            setIconImage(iconImage);
        }
        catch (Exception e)
        {
        }
    }


    protected JPanel createArgumentFieldPanel(CommandLineArgumentDefinition argDef,
                                                  boolean firstArg,
                                                  ListenerHelper helper,
                                                  String fieldValue)
    {
        switch (argDef.getDataType())
        {
            case ViewerArgumentDefinitionFactory.DELTA_MASS:
            case ViewerArgumentDefinitionFactory.FASTA_FILE:
            case ViewerArgumentDefinitionFactory.FEATURE_FILE:
                break;
            default:
                return super.createArgumentFieldPanel(argDef, firstArg, helper, fieldValue);
        }
        
        JPanel fieldPanel = new JPanel();
        GridBagConstraints fieldGBC = new GridBagConstraints();
        fieldGBC.gridwidth=GridBagConstraints.REMAINDER;
        fieldGBC.anchor = GridBagConstraints.LINE_START;

        boolean shouldAddFileChooser = false;
        JComponent argComponent = null;

        boolean fieldHasValue = (fieldValue != null && fieldValue.length() > 0);
        //special handling for each data type
        switch (argDef.getDataType())
        {
            case ViewerArgumentDefinitionFactory.DELTA_MASS:
                JTextField deltaMassTextField = new JTextField();
                deltaMassTextField.setPreferredSize(new Dimension(80, 20));
                deltaMassTextField.setMinimumSize(new Dimension(80, 20));

                if (fieldHasValue)
                    deltaMassTextField.setText(fieldValue);
                argComponent = deltaMassTextField;
                break;
            case ViewerArgumentDefinitionFactory.FASTA_FILE:
            case ViewerArgumentDefinitionFactory.FEATURE_FILE:
                JTextField fileTextField = new JTextField();
                fileTextField.setPreferredSize(new Dimension(225, 20));
                fileTextField.setMinimumSize(new Dimension(225, 20));

                if (fieldHasValue)
                    fileTextField.setText(fieldValue);
                argComponent = fileTextField;
                shouldAddFileChooser = true;
                break;
            default:
                throw new IllegalArgumentException("Unknown argument type " + argDef.getArgumentName() +
                        " of type " + argDef.getDataType());
        }

        GridBagConstraints argComponentGBC = new GridBagConstraints();
        argComponentGBC.anchor = GridBagConstraints.LINE_START;
        if (shouldAddFileChooser)
            argComponentGBC.gridwidth=GridBagConstraints.RELATIVE;
        else
            argComponentGBC.gridwidth=GridBagConstraints.REMAINDER;
        if (firstArg)
        {
            argComponentGBC.insets = new Insets(10, 0, 0, 0);
        }
        fieldPanel.add(argComponent, argComponentGBC);
        //add a file chooser that drives off of and populates the text field,
        //if this is a file data type
        if (shouldAddFileChooser)
            addFileChooser(argDef.getArgumentName(), false, helper, fieldPanel);

        argComponentMap.put(argDef, argComponent);
        

        return fieldPanel;
    }

    public void buttonSaveCommandFile_actionPerformed(ActionEvent event)
    {
        WorkbenchFileChooser wfc = new WorkbenchFileChooser();
        int chooserStatus = wfc.showOpenDialog(this);
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return;

        final File file = wfc.getSelectedFile();
        if (null != file)
        {
            Map<String,String> argNameValueMap = parseFieldArguments();
            String commandFileText =
                    CommandFileRunner.createCommandFileEntry(module.getCommandName(),
                            argNameValueMap);
            PrintWriter outPW = null;
            try
            {
                outPW = new PrintWriter(file);
                outPW.print(commandFileText);
                infoMessage(TextProvider.getText("FILE_FILENAME_SAVED",
                            file.getAbsolutePath()));
            }
            catch (FileNotFoundException e)
            {
                infoMessage(TextProvider.getText("ERROR") + ":" + e.getMessage());
            }
            finally
            {
                if (outPW != null)
                    outPW.close();
            }
        }
        notifyDone(event);
    }


    public void buttonShowCommand_actionPerformed(ActionEvent event)
    {
        Map<String,String> argNameValueMap = parseFieldArguments();

        StringBuffer commandLineCommand =
                new StringBuffer("msinspect --" + module.getCommandName());
        for (String argName : argNameValueMap.keySet())
        {
            commandLineCommand.append(" ");

            String argValue = argNameValueMap.get(argName);
            //change separator strings to spaces
            if (argName.equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
                argValue = argValue.replaceAll(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR, " ");

            //show the arg name and an equals sign iff it's a named argument
            if (!argName.equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT) &&
                !argName.equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
                commandLineCommand.append(argName + "=");

            commandLineCommand.append(argValue);
        }
        System.err.println("Command:");
        System.err.println(commandLineCommand.toString());
    }

}

