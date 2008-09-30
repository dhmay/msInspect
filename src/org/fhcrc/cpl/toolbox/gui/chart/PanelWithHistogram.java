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

package org.fhcrc.cpl.toolbox.gui.chart;

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.ChartFactory;
import org.jfree.data.statistics.HistogramDataset;

import java.util.List;

/**
 * PanelWithChart implementation to make it easy to put out scatterplots.
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithHistogram extends PanelWithChart
{
    public static final int DEFAULT_BREAKS = 100;

    public static final double BAR_OFFSET = .13;

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

    public PanelWithHistogram(List<Float> data, String name)
    {
        this(data, name, DEFAULT_BREAKS);
    }

    public PanelWithHistogram(List<Float> data, String name, int breaks)
    {
        double[] dataArray = new double[data.size()];
        for (int i=0; i<dataArray.length; i++)
            dataArray[i] = data.get(i);
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
                true,false,false);

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
        dataset.addSeries(name,data,breaks);

        buildChart();
    }

    public void addData(List<Double> data, String name)
    {
        double[] dataArray = new double[data.size()];
        for (int i=0; i<dataArray.length; i++)
            dataArray[i] = data.get(i);
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
}
