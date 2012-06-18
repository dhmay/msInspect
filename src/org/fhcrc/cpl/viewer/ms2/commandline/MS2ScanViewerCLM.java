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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.ms2.gui.MS2ScanViewer;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.awt.*;


/**
 * test
 */
public class MS2ScanViewerCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MS2ScanViewerCLM.class);

    protected File runFile;
    protected float mass;
    protected float massTolerancePPM = 20;
    protected List<Integer> scans=  new ArrayList<Integer>();

    JLabel scanInfoLabel = null;
    protected MS2ScanViewer.MultiMS2ScanViewer multiMS2ScanViewer;

    int numPeaksToLabel = 0;


    public MS2ScanViewerCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "ms2scanviewer";

        mHelpMessage ="ms2scanviewer";
        mShortDescription = "ms2scanviewer";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedFileArgumentDefinition(true, "mzXML file"),
                       new StringListArgumentDefinition("scans", false, "Scan(s) to view"),
                       new DecimalArgumentDefinition("mass", false, "Mass to search for scans around"),
                       new DecimalArgumentDefinition("masstoleranceppm", false, "PPM mass tolerance", massTolerancePPM),
                       new IntegerArgumentDefinition("numpeakstolabel", false,
                               "Number of the highest peaks whose m/z values should be labeled", numPeaksToLabel),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        runFile = getUnnamedFileArgumentValue();
        if (hasArgumentValue("mass"))
            mass = getFloatArgumentValue("mass");
        if (hasArgumentValue("scans"))
        {
            List<String> scanStrings = getStringListArgumentValue("scans");
            for (String scanString : scanStrings)
                scans.add(Integer.parseInt(scanString));
        }
        massTolerancePPM = getFloatArgumentValue("masstoleranceppm");
        numPeaksToLabel = getIntegerArgumentValue("numpeakstolabel");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
	scanInfoLabel = new JLabel();
        MSRun run;
        try
        {
            run = MSRun.load(runFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to load run from file " + runFile.getAbsolutePath());
        }

        List<Integer> scanNumbersToView = new ArrayList<Integer>();
        if (!scans.isEmpty())
            scanNumbersToView.addAll(scans);
        else
        {
            for (MSRun.MSScan ms2Scan : run.getMS2Scans())
            {
                float precursorMass = (ms2Scan.getPrecursorMz() - 1.0072766f) * ms2Scan.getPrecursorCharge();
                float deltaMassPPM = (precursorMass - mass) * 1000000 / mass;

                if (Math.abs(deltaMassPPM) <= massTolerancePPM)
                {
                    ApplicationContext.infoMessage("Adding scan " + ms2Scan.getNum() + " with precursor m/z " +
                            ms2Scan.getPrecursorMz() + ", inferred mass " + precursorMass + ", deltaMassPPM=" + deltaMassPPM);
                    scanNumbersToView.add(ms2Scan.getNum());
                }
            }

            if (scanNumbersToView.isEmpty())
            {
                ApplicationContext.infoMessage("No scans with precursors within " + massTolerancePPM +
                        "PPM of " + mass);
                return;
            }

            ApplicationContext.infoMessage(scanNumbersToView.size() + " scans match");
        }

        multiMS2ScanViewer = new MS2ScanViewer.MultiMS2ScanViewer(run, scanNumbersToView, numPeaksToLabel);
        scanInfoLabel = new JLabel("Scan , Precursor m/z: ");
        scanInfoLabel.setVisible(true);
        ChangeListener scanChangeListener = new MS2ScanChangedListener();
        multiMS2ScanViewer.addChangeListener(scanChangeListener);
        JDialog dialog = new JDialog();
        GridBagConstraints fullRowGBC = new GridBagConstraints();
        fullRowGBC.gridwidth = GridBagConstraints.REMAINDER;
        dialog.setLayout(new GridBagLayout());

        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setTitle("MS2 Scans");
        dialog.setSize(new Dimension(800, 600));
        dialog.add(scanInfoLabel, fullRowGBC);
        dialog.add(multiMS2ScanViewer, fullRowGBC);
        dialog.setVisible(true);      


    }

    protected class MS2ScanChangedListener implements ChangeListener
    {
        public void stateChanged(ChangeEvent e)
        {
            MSRun.MSScan scan = multiMS2ScanViewer.getMs2ScanViewer().getScanInViewer();
            scanInfoLabel.setText("Scan " + scan.getNum() + ", Precursor m/z: " + scan.getPrecursorMz());
            scanInfoLabel.updateUI();
        }
    }
}

