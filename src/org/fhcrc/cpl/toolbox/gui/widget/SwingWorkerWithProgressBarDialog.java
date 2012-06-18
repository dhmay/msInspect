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
package org.fhcrc.cpl.toolbox.gui.widget;


import javax.swing.*;
import java.awt.*;


/**
 * This is a SwingWorker that pops up a ProgressBar in a dialog box and updates the status.
 * Implementing classes will have to wire up something to do the actual status update, using
 * updateLabelText() and progressBar.setValue().
 *
 * In addition to the progress bar, the dialog contains a status label.  Text is specified in
 * labelTextExpression, with two optional tokens, for the current and maximum value.  
 *
 * See ProteinSummarySelectorFrame for an implementing class
 */
public abstract class SwingWorkerWithProgressBarDialog<T,V> extends
        SwingWorker<T,V>
{
    protected static final int PROGRESS_BAR_WIDTH = 250;
    protected static final int PROGRESS_BAR_HEIGHT = 50;

    protected JProgressBar progressBar;
    protected JDialog progressDialog;
    protected JDialog parent;
    protected JLabel statusLabel;
    protected String labelTextExpression;

    public static final String CURRENT_VALUE_TOKEN = "%%CURRENT_VALUE%%";
    public static final String MAX_VALUE_TOKEN = "%%MAX_VALUE%%";


    public SwingWorkerWithProgressBarDialog(JDialog parent,
                                     int progressBarMin, int progressBarMax, int progressBarStartingValue,
                                     String labelTextExpression, String title)
    {
        this.parent = parent;
        this.labelTextExpression = labelTextExpression;

        progressBar = new JProgressBar(progressBarMin, progressBarMax);
        progressBar.setSize(PROGRESS_BAR_WIDTH, PROGRESS_BAR_HEIGHT);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);
        progressDialog = new JDialog(parent);
        progressDialog.setSize(PROGRESS_BAR_WIDTH + 10, PROGRESS_BAR_HEIGHT + 70);
        progressDialog.setPreferredSize(new Dimension(PROGRESS_BAR_WIDTH + 10, PROGRESS_BAR_HEIGHT + 70));
        progressDialog.pack();

        progressDialog.setLocation(
                Math.max((int) (parent.getLocation().getX() + (parent.getWidth() / 2) - progressDialog.getWidth()), 0),
                Math.max((int) (parent.getLocation().getY() + (parent.getHeight() / 2) - progressDialog.getHeight()),0));

        progressDialog.setTitle(title);
        JPanel progressContainer = new JPanel();
        progressContainer.setSize(PROGRESS_BAR_WIDTH + 10, PROGRESS_BAR_HEIGHT + 10);
        progressDialog.setContentPane(progressContainer);

        progressContainer.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(8,0,8,0);

        statusLabel = new JLabel();

        updateLabelText(progressBarStartingValue);

        if (labelTextExpression != null)
            progressContainer.add(statusLabel, gbc);
        progressContainer.add(progressBar, gbc);
        progressDialog.setVisible(true);
    }

    /**
     * Set the text of the status label
     * @param currentValue
     */
    protected void updateLabelText(int currentValue)
    {
        String labelText = labelTextExpression;
        labelText = labelText.replace(CURRENT_VALUE_TOKEN, "" + currentValue);
        labelText = labelText.replace(MAX_VALUE_TOKEN, "" + progressBar.getMaximum());
        statusLabel.setText(labelText);
    }


}
