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
package org.fhcrc.cpl.viewer.ms2.gui;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.FloatRange;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import java.util.List;
import java.util.ArrayList;
import java.awt.event.ActionEvent;
import java.awt.*;


/**
 * Enables viewing of a single MS2 scan, with the ability to change the scan being viewed
 */
public class MS2ScanViewer extends JPanel
{
    protected static Logger _log = Logger.getLogger(MS2ScanViewer.class);

    protected MSRun.MSScan scanInViewer = null;
    protected PanelWithChart chartInViewer = null;

    public static final int DEFAULT_RESAMPLING_RESOLUTION = 50;
    protected int resamplingResolution = DEFAULT_RESAMPLING_RESOLUTION;

    public MS2ScanViewer()
    {

    }

    protected void buildChart()
    {
        if (scanInViewer == null)
            return;
        if (chartInViewer != null)
        {
            remove(chartInViewer);
        }
        float[][] spectrum = scanInViewer.getSpectrum();
        float[] mzValues = new float[spectrum[0].length];
        float[] intensityValues = new float[spectrum[0].length];

        for (int i=0; i<mzValues.length; i++)
        {
            mzValues[i] = spectrum[0][i];
            intensityValues[i] = spectrum[1][i];
        }

        chartInViewer = new PanelWithPeakChart(mzValues, intensityValues, "ms2scan");
        add(chartInViewer);
        updateUI();
    }


    public MSRun.MSScan getScanInViewer()
    {
        return scanInViewer;
    }

    public void setScanInViewer(MSRun.MSScan scanInViewer)
    {
        this.scanInViewer = scanInViewer;
        buildChart();
    }

    public PanelWithChart getChartInViewer()
    {
        return chartInViewer;
    }

    public void setChartInViewer(PanelWithPeakChart chartInViewer)
    {
        this.chartInViewer = chartInViewer;
    }

    public int getResamplingResolution()
    {
        return resamplingResolution;
    }

    public void setResamplingResolution(int resamplingResolution)
    {
        this.resamplingResolution = resamplingResolution;
    }

    /**
     * View multiple MS/MS scans.  Buttons to navigate through the scans.
     */
    public static class MultiMS2ScanViewer extends JPanel
    {
        protected MS2ScanViewer ms2ScanViewer = null;
        protected int currentScanIndex;

        protected MSRun run;
        protected List<Integer> ms2ScanNumbers;

        JButton forwardButton = null;
        JButton backButton = null;

        List<ChangeListener> changeListeners = new ArrayList<ChangeListener>();


        public MultiMS2ScanViewer()
        {
            ms2ScanViewer = new MS2ScanViewer();

            ListenerHelper helper = new ListenerHelper(this);

            backButton = new JButton("<-");
            helper.addListener(backButton, "backButton_actionPerformed");

            forwardButton = new JButton("->");
            helper.addListener(forwardButton, "forwardButton_actionPerformed");

            JPanel buttonPanel = new JPanel();

            GridBagConstraints fullRowGBC = new GridBagConstraints();
            fullRowGBC.gridwidth = GridBagConstraints.REMAINDER;
            setLayout(new GridBagLayout());

            add(ms2ScanViewer, fullRowGBC);
            buttonPanel.add(backButton);
            buttonPanel.add(forwardButton, fullRowGBC);
            add(buttonPanel, fullRowGBC);
        }

        public MultiMS2ScanViewer(MSRun run, List<Integer> ms2ScanNumbers)
        {
            this();

            this.run = run;
            this.ms2ScanNumbers = ms2ScanNumbers;
            this.currentScanIndex = 0;

            displayScanNumber(ms2ScanNumbers.get(0));
        }

        public void forwardButton_actionPerformed(ActionEvent event)
        {
            currentScanIndex++;
            displayScanIndex(currentScanIndex);
        }

        public void backButton_actionPerformed(ActionEvent event)
        {
            currentScanIndex--;
            displayScanIndex(currentScanIndex);
        }

        protected void updateButtonUI()
        {
            if (currentScanIndex < ms2ScanNumbers.size()-1)
                forwardButton.setEnabled(true);
            else forwardButton.setEnabled(false);
            forwardButton.updateUI();

            if (currentScanIndex > 0)
                backButton.setEnabled(true);
            else backButton.setEnabled(false);
            backButton.updateUI();
        }

        public void displayScanIndex(int scanIndex)
        {
            displayScanNumber(ms2ScanNumbers.get(scanIndex));
        }

        public void displayScanNumber(int scanNumber)
        {
            _log.debug("MultiMS2ScanViewer displaying scan " + scanNumber);

            ms2ScanViewer.setScanInViewer(run.getMS2Scan(run.getIndexForMS2ScanNum(scanNumber)));
            updateButtonUI();
            for (ChangeListener listener : changeListeners)
            {
                listener.stateChanged(new ChangeEvent("changed"));
            }
        }

        public void addChangeListener(ChangeListener listener)
        {
            changeListeners.add(listener);
        }

        public MS2ScanViewer getMs2ScanViewer()
        {
            return ms2ScanViewer;
        }
    }
}
