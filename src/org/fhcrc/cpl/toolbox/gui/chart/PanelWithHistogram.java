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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYZDataset;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.awt.*;

/**
 * PanelWithChart implementation to make it easy to put out scatterplots.
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithHistogram extends PanelWithChart
{
    public static final int DEFAULT_BREAKS = 100;

    public static final double BAR_OFFSET = .13;

    protected double[] dataValues;
    protected HistogramDataset dataset;

    int breaks = DEFAULT_BREAKS;

    String title = "";

    protected boolean offsetSeries = false;

    public PanelWithHistogram()
    {
        super();
    }

    public PanelWithHistogram(String title)
    {
        this();
        this.title = title;
    }

    public PanelWithHistogram(double[] data)
    {
        this(data,"");
    }

    public PanelWithHistogram(double[] data,
                              String name)
    {
        this(data, name, DEFAULT_BREAKS);
    }

    public PanelWithHistogram(float[] data)
    {
        this(data,"");
    }

    public PanelWithHistogram(List<? extends Number> data, String name)
    {
        this(data, name, DEFAULT_BREAKS);
    }

    public PanelWithHistogram(List<? extends Number> data, String name, int breaks)
    {
        double[] dataArray = new double[data.size()];
        for (int i=0; i<dataArray.length; i++)
            dataArray[i] = data.get(i).doubleValue();
        init(dataArray,name, breaks);
    }

    public PanelWithHistogram(float[] data,
                              String name)
    {
        this(data, name, DEFAULT_BREAKS);
    }

    public PanelWithHistogram(float[] data,
                              String name, int breaks)
    {
        this();
        double[] dataDouble = new double[data.length];
        for (int i=0; i<data.length; i++)
            dataDouble[i] = data[i];

        init(dataDouble,name,breaks);

    }


    public PanelWithHistogram(double[] data,
                              String name,
                              int breaks)
    {
        this();
        init(data,name,breaks);
    }

    protected void init(double[] data, String name, int breaks)
    {
        setName(name);
        this.breaks = breaks;
        addData(data, name);
//        buildChart();
    }

    public void buildChart()
    {
        _chart = ChartFactory.createHistogram(
                title,title,title, dataset,
                PlotOrientation.VERTICAL,
                true,true,false);
        
        //dhmay adding 2009/09/14.  As of jfree 1.0.13, shadows on by default, and gray background
        ((XYBarRenderer)_chart.getXYPlot().getRenderer()).setShadowVisible(false);


        init(_chart);
    }

    public void addData(double[] data, String name)
    {
        if (dataset == null)
            dataset = new HistogramDataset();
        if (offsetSeries)
        {
            int dataSetIndex = dataset.getSeriesCount();
            if (dataSetIndex > 0)
            {
                for (int i=0; i<data.length; i++)
                {
                    data[i] += (dataSetIndex *BAR_OFFSET );
                }
            }
        }
        dataValues = data;
        dataset.addSeries(name,data,breaks);
        buildChart();
    }


    public void addData(List<? extends Number> data, String name)
    {
        double[] dataArray = new double[data.size()];
        for (int i=0; i<dataArray.length; i++)
            dataArray[i] = data.get(i).doubleValue();
        addData(dataArray, name);

    }

    public HistogramDataset getDataset()
    {
        return dataset;
    }


    public boolean isOffsetSeries()
    {
        return offsetSeries;
    }

    public void setOffsetSeries(boolean offsetSeries)
    {
        this.offsetSeries = offsetSeries;
    }

    public int getBreaks()
    {
        return breaks;
    }

    public void setBreaks(int breaks)
    {
        this.breaks = breaks;
    }

    protected void saveChartDataToFile(File outFile, String delimiter)
    {
        _log.debug("PanelWithHistogram saveChartDataToFile 1, *delimiter*=*" + delimiter + "*");

        PrintWriter pw = null;

        try
        {
            pw = new PrintWriter(outFile);

            pw.println("value");

            for (double dataValue : dataValues)
            {
                pw.println("" + dataValue);
                pw.flush();
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_SAVING_CHART_DATA"), e);
        }
        finally
        {
            if (pw != null)
                pw.close();
        }
    }
}
