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
package org.fhcrc.cpl.viewer.ms2.gui;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.XYPlot;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import java.util.Collections;
import java.util.Comparator;
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
    protected PanelWithPeakChart chartInViewer = null;

    public static final int DEFAULT_RESAMPLING_RESOLUTION = 50;
    protected int resamplingResolution = DEFAULT_RESAMPLING_RESOLUTION;

    protected int numHighestPeaksToLabel = 0;

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

        //add labels for highest peaks
        if (numHighestPeaksToLabel > 0) {
            XYPlot xyPlot = chartInViewer.getChart().getXYPlot();
            List<Pair<Float, Float>> peaksAsPairs = new ArrayList<Pair<Float, Float>>();
            for (int i=0; i<mzValues.length; i++)
                peaksAsPairs.add(new Pair<Float, Float>(mzValues[i], intensityValues[i]));
            Collections.sort(peaksAsPairs, new Comparator<Pair<Float, Float>>() {
                @Override
                public int compare(Pair<Float, Float> pair1, Pair<Float, Float> pair2) {
                    float diff = pair1.second - pair2.second;
                    if (diff > 0) return -1;
                    if (diff < 0) return 1;
                    return 0;
                }
            });
            for (int i=0; i<numHighestPeaksToLabel; i++) {
                Pair<Float, Float> pair = peaksAsPairs.get(i);
                xyPlot.addAnnotation(new XYTextAnnotation("" + Rounder.round(pair.first, 3), pair.first, pair.second));
            }
        }

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

        protected MSRun.MSScan[] scans = null;

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

        public MultiMS2ScanViewer(MSRun run, List<Integer> ms2ScanNumbers, int numPeaksToLabel)
        {
            this();

            MSRun.MSScan[] scans = new MSRun.MSScan[ms2ScanNumbers.size()];
            for (int i=0; i< ms2ScanNumbers.size(); i++) {
                scans[i] = run.getMS2Scan(run.getIndexForMS2ScanNum(ms2ScanNumbers.get(i)));
            }

            init(scans, numPeaksToLabel);
        }

        public MultiMS2ScanViewer(MSRun.MSScan[] scans, int numPeaksToLabel)
        {
            this();
            init(scans, numPeaksToLabel);
        }

        protected void init(MSRun.MSScan[] scans, int numPeaksToLabel) {
            ms2ScanViewer.setNumHighestPeaksToLabel(numPeaksToLabel);

            this.scans = scans;
            this.currentScanIndex = 0;

            ms2ScanViewer.setScanInViewer(scans[0]);
        }

        public void forwardButton_actionPerformed(ActionEvent event)
        {
            currentScanIndex++;
            displayScan(scans[currentScanIndex]);
        }

        public void backButton_actionPerformed(ActionEvent event)
        {
            currentScanIndex--;
            displayScan(scans[currentScanIndex]);
        }

        protected void updateButtonUI()
        {
            if (currentScanIndex < scans.length-1)
                forwardButton.setEnabled(true);
            else forwardButton.setEnabled(false);
            forwardButton.updateUI();

            if (currentScanIndex > 0)
                backButton.setEnabled(true);
            else backButton.setEnabled(false);
            backButton.updateUI();
        }


        public void displayScan(MSRun.MSScan scan)
        {
            _log.debug("MultiMS2ScanViewer displaying scan " + scan.getNum());

            ms2ScanViewer.setScanInViewer(scan);
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

    public int getNumHighestPeaksToLabel() {
        return numHighestPeaksToLabel;
    }

    public void setNumHighestPeaksToLabel(int numHighestPeaksToLabel) {
        this.numHighestPeaksToLabel = numHighestPeaksToLabel;
    }
}
