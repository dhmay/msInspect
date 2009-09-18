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
package org.fhcrc.cpl.viewer.quant.gui;

import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.Rounder;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.XYPlot;

import javax.swing.*;
import java.util.List;
import java.util.ArrayList;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class PanelWithLogRatioHistAndFields extends JPanel
{
    protected static Logger _log = Logger.getLogger(PanelWithLogRatioHistAndFields.class);

    protected float maxRatioBound = 20f;
    protected float minRatioBound = 1f / maxRatioBound;


    protected JLabel maxLowRatioLabel;
    protected JLabel minHighRatioLabel;
    protected JLabel numPassingRatiosLabel;

    protected PanelWithHistogram logRatioHistogram;

    //min and max ratio for retention.  Ratio must be > min OR < max, or both
    protected float minHighRatio = 0f;
    protected float maxLowRatio = 999f;

    List<ActionListener> rangeUpdateListeners = new ArrayList<ActionListener>();

    protected List<Float> logRatios;

    protected LogRatioHistMouseListener histMouseListener;

    public PanelWithLogRatioHistAndFields()
    {
        initGUI();
    }

    public PanelWithLogRatioHistAndFields(float maxLowRatio, float minHighRatio)
    {
        this();
        this.maxLowRatio = maxLowRatio;
        this.minHighRatio = minHighRatio;
        updateFieldText();        
    }

    public PanelWithLogRatioHistAndFields(float maxLowRatio, float minHighRatio, List<Float> logRatios)
    {
        this(maxLowRatio, minHighRatio);
        setLogRatios(logRatios);
    }

    protected void initGUI()
    {
        setLayout(new GridBagLayout());
        maxLowRatioLabel = new JLabel("Max Low Ratio: ");
        minHighRatioLabel = new JLabel("Min High Ratio: ");
        numPassingRatiosLabel = new JLabel("Retained: ");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.LINE_START;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 1;
        gbc.weightx = 1;
        gbc.gridwidth = 1;
        add(maxLowRatioLabel, gbc);
        gbc.anchor = GridBagConstraints.LINE_END;
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        add(minHighRatioLabel, gbc);
        gbc.anchor = GridBagConstraints.LINE_END;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        add(numPassingRatiosLabel, gbc);
    }

    public void setLogRatios(List<Float> logRatios)
    {
        if (logRatioHistogram != null)
        {
            remove(logRatioHistogram);
        }
        this.logRatios = logRatios;


        float minLogRatioBound = (float) Math.log(minRatioBound);
        float maxLogRatioBound = (float) Math.log(maxRatioBound);

        List<Float> boundedLogRatios = new ArrayList<Float>(logRatios.size());
        for (float logRatio : logRatios)
        {
                boundedLogRatios.add((float) Math.min(maxLogRatioBound,
                    Math.max(minLogRatioBound, logRatio)));
        }

        logRatioHistogram = new PanelWithHistogram(boundedLogRatios, "Log Ratios", 200);
        Dimension histDimension = new Dimension(300, 80);

        logRatioHistogram.setPreferredSize(histDimension);
        logRatioHistogram.getChart().removeLegend();


        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.LINE_START;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 1;
        gbc.weightx = 1;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        add(logRatioHistogram, gbc);


        histMouseListener = new LogRatioHistMouseListener(logRatioHistogram);
        ChartPanel histChartPanel = logRatioHistogram.getChartPanel();
        histChartPanel.setMouseZoomable(false);
        logRatioHistogram.getChartPanel().addMouseListener(histMouseListener);
        logRatioHistogram.getChartPanel().addMouseMotionListener(histMouseListener);
        histMouseListener.addRangeUpdateListener(new LogRatioHistogramListener(this));

        logRatioHistogram.updateUI();
        //if there are specified minHigh and maxHigh values, and they're valid, draw the initially selected region
        if (minHighRatio > maxLowRatio && minHighRatio > 0 && maxLowRatio > 0)
            histMouseListener.setSelectedRegionWithChartValues((float)Math.log(maxLowRatio),
                    (float)Math.log(minHighRatio));


        //remove axes from chart
        ((XYPlot)logRatioHistogram.getPlot()).getRangeAxis().setVisible(false);
        ((XYPlot)logRatioHistogram.getPlot()).getDomainAxis().setVisible(false);

        logRatioHistogram.updateUI();


    }

    public void setSize(int width, int height)
    {
        super.setSize(width, height);
        super.setPreferredSize(new Dimension(width, height));
        int chartWidth = width;
        int chartHeight = height - 50;
        Dimension histDimension =  new Dimension(chartWidth, chartHeight);
        if (logRatioHistogram != null)
        {
            logRatioHistogram.setPreferredSize(histDimension);
            logRatioHistogram.updateUI();
        }
    }

    public void updateFieldText()
    {
        maxLowRatioLabel.setText("Max Low Ratio: " + Rounder.round(maxLowRatio, 2));
        minHighRatioLabel.setText("Min High Ratio: " + Rounder.round(minHighRatio, 2));
        if (logRatios != null)
        {
            int numRetained = 0;
            float maxLowLogRatio = (float) Math.log(maxLowRatio);
            float minHighLogRatio = (float) Math.log(minHighRatio);

            for (float logRatio : logRatios)
            {
                if (logRatio <= maxLowLogRatio || logRatio >= minHighLogRatio)
                    numRetained++;
//                else if (Math.exp(logRatio) < 0.2) System.err.println("FAIL: logRatio: " + logRatio + ", ratio: " + Math.exp(logRatio) + ", maxLowLog: " + maxLowLogRatio + ", maxLo: " + maxLowRatio);
            }
            numPassingRatiosLabel.setText("Retained: " + numRetained + " / " +
                    logRatios.size() +
                    " (" + Rounder.round(100f * numRetained / logRatios.size(),1) + "%)");
        }
        else numPassingRatiosLabel.setText("Retained: ");


    }

    public float getMinHighRatio()
    {
        return minHighRatio;
    }

    public void setMinHighRatio(float minHighRatio)
    {
        this.minHighRatio = minHighRatio;
        updateFieldText();
    }

    public float getMaxLowRatio()
    {
        return maxLowRatio;
    }

    public void setMaxLowRatio(float maxLowRatio)
    {
        this.maxLowRatio = maxLowRatio;
        updateFieldText();        
    }

    public void addRangeUpdateListener(ActionListener listener)
    {
        rangeUpdateListeners.add(listener);
    }

    /**
     * A chart listener that picks up events indicating changes to the selected area
     */
    protected class LogRatioHistogramListener implements ActionListener
    {
        protected PanelWithLogRatioHistAndFields parent;
        public LogRatioHistogramListener(PanelWithLogRatioHistAndFields parent)
        {
            this.parent = parent;
        }
        public void actionPerformed(ActionEvent e)
        {
            LogRatioHistMouseListener logRatioHistMouseListener = (LogRatioHistMouseListener) e.getSource();
            minHighRatio = (float) Math.exp(logRatioHistMouseListener.getSelectedRealXMaxValue());
            maxLowRatio = (float) Math.exp(logRatioHistMouseListener.getSelectedRealXMinValue());
            updateFieldText();

            for (ActionListener listener : rangeUpdateListeners)
                listener.actionPerformed(new ActionEvent(parent, 0, maxLowRatio + "," + minHighRatio));
        }
    }

    public float getMaxRatioBound()
    {
        return maxRatioBound;
    }

    public void setMaxRatioBound(float maxRatioBound)
    {
        this.maxRatioBound = maxRatioBound;
    }

    public float getMinRatioBound()
    {
        return minRatioBound;
    }

    public void setMinRatioBound(float minRatioBound)
    {
        this.minRatioBound = minRatioBound;
    }
}
