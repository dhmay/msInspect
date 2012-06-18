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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.MzXmlWriter;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import java.io.File;

/**
 * Dialog window class for saving a partial mzXML file
 */
public class SavePartialMzxmlDialog extends JDialog
{
    static Logger _log = Logger.getLogger(SavePartialMzxmlDialog.class);


    public JTextField textStartScan;
    public JTextField textEndScan;
    public JTextField textLowMz;
    public JTextField textHighMz;
    public JTextField textFilename;

    public JButton buttonBrowse;
    public JButton buttonSave;
    public JButton buttonCancel;

    public SavePartialMzxmlDialog()
    {
        super(ApplicationContext.getFrame(), TextProvider.getText("SAVE_PARTIAL_MZXML_DOTDOTDOT"));

        //graphical stuff
        Container contentPanel = null;
        try
        {
            contentPanel = Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/SavePartialMzxmlDialog.xml", this);
            setContentPane(contentPanel);
            pack();
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
            throw new RuntimeException(x);
        }

        ListenerHelper helper = new ListenerHelper(this);
        helper.addListener(buttonBrowse, "buttonBrowse_actionPerformed");
        helper.addListener(buttonCancel, "buttonCancel_actionPerformed");
        helper.addListener(buttonSave, "buttonSave_actionPerformed");

        getRootPane().setDefaultButton(buttonSave);

        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        textStartScan.setText(Integer.toString(run.getScan(0).getNum()));
        textEndScan.setText(Integer.toString(run.getScan(run.getScanCount() - 1).getNum()));
        textLowMz.setText(Float.toString(run.getMzRange().min));
        textHighMz.setText(Float.toString(run.getMzRange().max));
    }

    /**
     * Programmatically populate the textfields
     * @param startScan
     * @param endScan
     * @param lowMz
     * @param highMz
     */
    public void setRegionToSave(int startScan, int endScan, float lowMz, float highMz)
    {
        textStartScan.setText(Integer.toString(startScan));
        textEndScan.setText(Integer.toString(endScan));
        textLowMz.setText(Float.toString(lowMz));
        textHighMz.setText(Float.toString(highMz));
    }

    public void buttonBrowse_actionPerformed(ActionEvent event)
    {
        WorkbenchFileChooser chooser = new WorkbenchFileChooser();
        chooser.setDialogTitle(TextProvider.getText("FILE"));
        int chooserStatus = chooser.showOpenDialog(ApplicationContext.getFrame());
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
             return;        
        File file = chooser.getSelectedFile();
        if (file != null)
            textFilename.setText(file.getAbsolutePath());
    }

    public void buttonSave_actionPerformed(ActionEvent event)
    {
        boolean problems = false;
        int startScanNum = 0, endScanNum = 0;
        float lowMz = 0, highMz = 0;

        //Check the entered parameters
        try
        {
            startScanNum = Integer.parseInt(textStartScan.getText());
            endScanNum = Integer.parseInt(textEndScan.getText());
            lowMz = (float) Double.parseDouble(textLowMz.getText());
            highMz = (float) Double.parseDouble(textHighMz.getText());
        }
        catch (Exception ex)
        {
            problems = true;
        }
        if (textFilename.getText() == null || textFilename.getText().length() < 1)
            problems = true;
        if (startScanNum > endScanNum || lowMz > highMz)
            problems = true;

        if (problems)
        {
            ApplicationContext.infoMessage(TextProvider.getText("PLEASE_CHECK_VALUES_AND_TRY_AGAIN"));
        }
        else
        {
            _log.debug("Saving mzxml file, boundaries: " + startScanNum + "-" + endScanNum +
                       ", " + lowMz + "-" + highMz);
            MzXmlWriter mzXmlWriter = new MzXmlWriter((MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN));

            try
            {
                mzXmlWriter.writeSubregion(new File(textFilename.getText()), startScanNum,
                        endScanNum, lowMz, highMz);
            }
            catch (IOException ioe)
            {
                ApplicationContext.errorMessage(TextProvider.getText("ERROR_WRITING_PARTIAL_MZXML_FILE"), ioe);
            }
            this.setVisible(false);
            this.dispose();
        }
    }

    public void buttonCancel_actionPerformed(ActionEvent e)
    {
        this.setVisible(false);
        this.dispose();
    }

    public static class SavePartialMzxmlAction extends AbstractAction
    {
        public SavePartialMzxmlAction()
        {
            super("Save Partial mzXML file...");
        }

        public void actionPerformed(ActionEvent e)
        {
            SavePartialMzxmlDialog saveDialog = new SavePartialMzxmlDialog();
            saveDialog.setVisible(true);
        }
    }

}

