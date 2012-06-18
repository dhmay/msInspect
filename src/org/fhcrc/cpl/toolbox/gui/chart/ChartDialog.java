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

package org.fhcrc.cpl.toolbox.gui.chart;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.Plot;
import org.fhcrc.cpl.toolbox.TextProvider;

import javax.swing.*;
import java.awt.event.ComponentListener;
import java.awt.event.ComponentEvent;
import java.awt.*;
import java.io.File;
import java.io.IOException;

/**
 * Generic render-a-chart-in-a-dialog class.  Doesn't do much to start with,
 * but this could be a good repository for things that we like to do with charts.
 *
 * For instance, maybe it would be fun to add a way to export chart contents
 * to, say, a .tsv file.... It might be hard to overload the chart context menu,
 * but we could add buttons or something
 */
public class ChartDialog extends JDialog
{
    PanelWithChart _panelWithChart = null;

    public ChartDialog()
    {
        super();
        setTitle(TextProvider.getText("CHART"));
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    }

    public ChartDialog(Plot plot)
    {
        this();
        init(plot);
    }

    public ChartDialog(Plot plot, String title)
    {
        this();
        init(plot);
        setTitle(title);
    }

    public ChartDialog(JFreeChart chart)
    {
        this(chart.getPlot());
    }

    public ChartDialog(PanelWithChart panelWithChart)
    {
        this();
        _panelWithChart = panelWithChart;
        init();
//        setContentPane(panelWithChart);
    }

    public ChartDialog(JFreeChart chart, String title)
    {
        this(chart.getPlot(), title);        
    }

    protected void init()
    {
        setSize(800,600);
        initPanelWithChart();
    }

    protected void initPanelWithChart()
    {
        if (_panelWithChart != null)
        {         
            add(_panelWithChart);
            _panelWithChart.setMinimumSize(getSize());
            _panelWithChart.setPreferredSize(getSize());
        }
        addComponentListener(new ResizeListener());
    }

    protected class ResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
            if (_panelWithChart != null)
            {
                int newPanelWidth = getSize().width - 20;
                int newPanelHeight = getSize().height - 50;

                _panelWithChart.setPreferredSize(new Dimension(newPanelWidth, newPanelHeight));
                _panelWithChart.updateUI();
            }
        }

        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}  
    }



    protected void init(Plot plot)
    {
        _panelWithChart = new PanelWithChart(plot);
        init();
    }

    public PanelWithChart getPanelWithChart()
    {
        return _panelWithChart;
    }

    public void setPanelWithChart(PanelWithChart panelWithChart)
    {
        if (_panelWithChart != null)
            remove(_panelWithChart);
        _panelWithChart = panelWithChart;
        initPanelWithChart();
    }

    public void saveChartToImageFile(File outFile) throws IOException
    {
        _panelWithChart.saveChartToImageFile(outFile);
    }

    protected void saveChartDataToCSV(File outFile)
    {
        _panelWithChart.saveChartDataToCSV(outFile);
    }

    protected void saveChartDataToTSV(File outFile)
    {
        _panelWithChart.saveChartDataToTSV(outFile);
    }

    public void setShowLegend(boolean showLegend)
    {
        _panelWithChart.setShowLegend(showLegend);
    }

    public boolean isShowLegend()
    {
        return _panelWithChart.isShowLegend();
    }
}
