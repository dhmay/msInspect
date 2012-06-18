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

package org.fhcrc.cpl.toolbox.commandline.arguments;

import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.io.File;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public abstract class FileArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    static Logger _log = Logger.getLogger(FileArgumentDefinition.class);


    //known types of files
    public static final int FILE_TYPE_UNKNOWN = 0;
    public static final int FILE_TYPE_FEATURE = 1;
    public static final int FILE_TYPE_IMAGE = 2;

    protected int fileType = FILE_TYPE_UNKNOWN;

    public FileArgumentDefinition(String argumentName)
    {
        super(argumentName);

    }
    public FileArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public FileArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public FileArgumentDefinition(String argumentName, boolean required, String help, int fileType)
    {
        this(argumentName, required, help);
        setFileType(fileType);
    }

    /**
     * Do some processing on a raw file path to make it suitable to try to create a File from it.
     * I.e., handle escaped spaces and beginning and trailing quotes.  This should be called by all
     * subclasses, rather than just newing up a file using filePath
     * @param filePath
     * @return
     */
    public static File createFileFromRawPath(String filePath)
    {
        _log.debug("Raw filePath: " + filePath);
        filePath = filePath.replaceAll("\\\\ "," ");
        if (filePath.startsWith("\"") && filePath.endsWith("\""))
        {
            filePath = filePath.substring(1, filePath.length() -1);
        }
        _log.debug("Adjusted filePath: " + filePath);

        return new File(filePath);
    }

    public String getValueDescriptor()
    {
        return "<filepath>";
    }

    public static class GUIFileChooserButtonListener implements ActionListener
    {
        protected JTextField argTextField;
        protected JDialog parentDialog;
        protected boolean isMulti;
        protected boolean isDir;

        public GUIFileChooserButtonListener(JTextField argTextField, JDialog parentDialog,
                                            boolean isMulti, boolean isDir)
        {
            this.argTextField = argTextField;
            this.parentDialog = parentDialog;
            this.isMulti = isMulti;
            this.isDir = isDir;
        }

        public void actionPerformed(ActionEvent event)
        {
            JFileChooser fc = new JFileChooser();
            fc.setMultiSelectionEnabled(false);
            String currentFieldValue = argTextField.getText();
            File directory = null;
            if (currentFieldValue != null &&
                    currentFieldValue.length() > 0)
            {
                //if multiple files selected, get the first one
                if (currentFieldValue.contains(" "))
                    currentFieldValue = currentFieldValue.substring(0, currentFieldValue.indexOf(" "));
                File currentFile = new File(currentFieldValue);
                fc.setSelectedFile(currentFile);
                directory = currentFile.getParentFile();
            }
            else
            {
                directory = new File (".");
            }
            if (directory != null && directory.exists())
                fc.setCurrentDirectory(directory);
            if (isDir)
                fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

            int chooserStatus = fc.showOpenDialog(parentDialog);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;

            File[] files = fc.getSelectedFiles();
            if (!isMulti)
                files = new File[] { fc.getSelectedFile() };
            if (null != files && files.length > 0)
            {
                StringBuffer newFileTextBuf = new StringBuffer();
                for (int i=0; i<files.length; i++)
                {
                    if (i>0)
                        newFileTextBuf.append(" ");
                    newFileTextBuf.append(files[i].getAbsolutePath());
                }
                if (isMulti && argTextField.getText() != null && argTextField.getText().length() > 0)
                    argTextField.setText(argTextField.getText() + " " + newFileTextBuf.toString());
                else
                    argTextField.setText(newFileTextBuf.toString());
            }
        }
    }

    /**
     * Same as base method, but resize the text field
     * @param parent
     * @param parentDialog
     * @param defaultValue
     * @return
     */
    public JComponent addComponentsForGUI(Container parent, JDialog parentDialog, String defaultValue,
                                          boolean isMulti, boolean isDir)
    {
        JPanel fieldPanel = new JPanel();

        JTextField argTextField = new JTextField();
        argTextField.setPreferredSize(new Dimension(225, 20));
        argTextField.setMinimumSize(new Dimension(225, 20));

        if (defaultValue != null && defaultValue.length() > 0)
            argTextField.setText(defaultValue);

        GridBagConstraints argComponentGBC = new GridBagConstraints();
        argComponentGBC.anchor = GridBagConstraints.LINE_START;
        argComponentGBC.gridwidth = GridBagConstraints.RELATIVE;
        argComponentGBC.insets = new Insets(5,0,0,0);

        fieldPanel.add(argTextField, argComponentGBC);

        argComponentGBC.gridwidth = GridBagConstraints.REMAINDER;
        JButton chooserButton = new JButton(TextProvider.getText("BROWSE_DOTDOTDOT"));
        chooserButton.setActionCommand(getArgumentName());
        chooserButton.addActionListener(new GUIFileChooserButtonListener(
                argTextField, parentDialog, isMulti, isDir));

        fieldPanel.add(chooserButton, argComponentGBC);

        parent.add(fieldPanel, argComponentGBC);
        return argTextField;
    }

    public JComponent addComponentsForGUISeries(Container parent, JDialog parentDialog, String defaultValue,
                                                boolean isDir)
    {
        //if a series of files to read (most likely case), handle super-specially.  Else, a big text field
        JTextField argTextField = new JTextField();
        argTextField.setPreferredSize(new Dimension(220, 20));
        argTextField.setMinimumSize(new Dimension(220, 20));

        boolean fieldHasValue = (defaultValue != null && defaultValue.length() > 0);

        if (fieldHasValue)
        {
            defaultValue = defaultValue.replaceAll(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR, " ");
            argTextField.setText(defaultValue);
        }

        GridBagConstraints textFieldGBC = new GridBagConstraints();
        textFieldGBC.anchor = GridBagConstraints.LINE_START;
        textFieldGBC.gridwidth=GridBagConstraints.RELATIVE;
        textFieldGBC.insets = new Insets(5, 0, 0, 0);
        parent.add(argTextField, textFieldGBC);

        JButton chooserButton = new JButton(TextProvider.getText("BROWSE_DOTDOTDOT"));
        chooserButton.setActionCommand(getArgumentName());
        chooserButton.addActionListener(new GUIFileChooserButtonListener(
                argTextField, parentDialog, true, isDir));

        return argTextField;
    }

    public int getFileType()
    {
        return fileType;
    }

    public void setFileType(int fileType)
    {
        this.fileType = fileType;
    }
}
