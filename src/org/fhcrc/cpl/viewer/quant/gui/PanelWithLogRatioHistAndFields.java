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

/**
 * A goofy little GUI thing that shoves together a histogram of log-ratios with an overlay for selecting extreme
 * values, and 3 text fields summarizing the selection.  Kind of specific-seeming, but used in 3 different places
 * in Qurate.
 *
 * The layout is very sensitive to tweaks in dimensions.  Careful when sizing.
 *
 * todo: for some reason a line started showing up at x=0 (log of 1:1).  Would like to control that.
 */
public class PanelWithLogRatioHistAndFields extends JPanel
{
    protected static Logger _log = Logger.getLogger(PanelWithLogRatioHistAndFields.class);

    //boundaries for chart values, in non-log space
    protected float maxRatioBound = 40f;
    protected float minRatioBound = 1f / maxRatioBound;

    //Labels for selection summary info
    protected JLabel maxLowRatioLabel;
    protected JLabel minHighRatioLabel;
    protected JLabel numPassingRatiosLabel;

    protected PanelWithHistogram logRatioHistogram;

    //min and max ratio for retention.  Ratio must be > min OR < max, or both
    protected float minHighRatio = 0f;
    protected float maxLowRatio = 999f;

    //for marking a special value
    protected float domainCrosshairValue = Float.NaN;

    //objects to be notified when the user changes the selected range.  Notifications occur while dragging
    //and also when released
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

    /**
     * Lay out GUI components, but NOT the histogram
     */
    protected void initGUI()
    {
        setLayout(new GridBagLayout());
        maxLowRatioLabel = new JLabel("Max Low Ratio: ");
        minHighRatioLabel = new JLabel("Min High Ratio: ");
        numPassingRatiosLabel = new JLabel("Retained: ");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 0;
        gbc.weightx = 1;        
        gbc.gridwidth = 1;
        add(maxLowRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        add(minHighRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        add(numPassingRatiosLabel, gbc);
    }

    /**
     * Set the log ratios, build the histogram and display, removing the old one if there was one.  Nothing gets
     * cleaned up related to the old chart; it'll just hang around intil GC
     * todo: do I need to dispose of the old chart in a better way?
     * @param logRatios
     */
    public void setLogRatios(List<Float> logRatios)
    {
        _log.debug("setLogRatios 1");
        if (logRatioHistogram != null)
        {
            remove(logRatioHistogram);
        }
        this.logRatios = logRatios;

        float minLogRatioBound = (float) Math.log(minRatioBound);
        float maxLogRatioBound = (float) Math.log(maxRatioBound);

        _log.debug("setLogRatios 2");


        List<Float> boundedLogRatios = new ArrayList<Float>(logRatios.size());
        for (float logRatio : logRatios)
            boundedLogRatios.add(Math.min(maxLogRatioBound, Math.max(minLogRatioBound, logRatio)));

        logRatioHistogram = new PanelWithHistogram(boundedLogRatios, "Log Ratios", 200);
        Dimension histDimension = new Dimension(300, 80);
        if (!Float.isNaN(domainCrosshairValue))
        {
            logRatioHistogram.getChart().getXYPlot().setDomainCrosshairValue(domainCrosshairValue, true);
            logRatioHistogram.getChart().getXYPlot().setDomainCrosshairVisible(true);
        }

        _log.debug("setLogRatios 1");


        logRatioHistogram.setPreferredSize(histDimension);
        logRatioHistogram.getChart().removeLegend();
        logRatioHistogram.getChart().getXYPlot().getDomainAxis().setLowerBound(minLogRatioBound);
        logRatioHistogram.getChart().getXYPlot().getDomainAxis().setUpperBound(maxLogRatioBound);


        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.LINE_START;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 10;
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
            updateSelectedRegion();
        
        //remove axes from chart
        ((XYPlot)logRatioHistogram.getPlot()).getRangeAxis().setVisible(false);
        ((XYPlot)logRatioHistogram.getPlot()).getDomainAxis().setVisible(false);

        logRatioHistogram.updateUI();
    }


    /**
     * Resizing the histogram is tricky
     * @param width
     * @param height
     */
    public void setSize(int width, int height)
    {
        super.setSize(width, height);
        super.setPreferredSize(new Dimension(width, height));
        int chartWidth = width;
        int chartHeight = height - 45;
        Dimension histDimension =  new Dimension(chartWidth, chartHeight);
        if (logRatioHistogram != null)
        {
            logRatioHistogram.setPreferredSize(histDimension);
            logRatioHistogram.updateUI();
        }
    }

    /**
     * Update all 3 text fields after selection change
     */
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

    /**
     * Also updates region selected and text of fields
     * @param minHighRatio
     */
    public void setMinHighRatio(float minHighRatio)
    {
        this.minHighRatio = minHighRatio;
        updateFieldText();
        updateSelectedRegion();
    }

    public float getMaxLowRatio()
    {
        return maxLowRatio;
    }

    /**
     * Also updates region selected and text of fields
     * @param maxLowRatio
     */
    public void setMaxLowRatio(float maxLowRatio)
    {
        this.maxLowRatio = maxLowRatio;
        updateFieldText();
        updateSelectedRegion();
    }

    protected void updateSelectedRegion()
    {
        if (histMouseListener != null)
            histMouseListener.setSelectedRegionWithChartValues((float)Math.log(maxLowRatio),
                    (float)Math.log(minHighRatio));
    }

    public void addRangeUpdateListener(ActionListener listener)
    {
        rangeUpdateListeners.add(listener);
    }

    /**
     * A chart listener that picks up events indicating changes to the selected area, recalculates the
     * low and high ratios defined (in not-log space), and passes on the notification to all listeners
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

    public float getDomainCrosshairValue()
    {
        return domainCrosshairValue;
    }

    public void setDomainCrosshairValue(float domainCrosshairValue)
    {
        this.domainCrosshairValue = domainCrosshairValue;
    }
}
