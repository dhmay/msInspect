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

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PiePlot;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.category.CategoryToPieDataset;
import org.jfree.data.general.PieDataset;
import org.jfree.data.general.DefaultPieDataset;

import java.util.List;
import java.util.Map;
import java.awt.*;

/**
 * A simple pie-chart implementation.  Could get more complicated.
 */
public class PanelWithPieChart extends PanelWithChart
{
    public static final int DEFAULT_LABEL_FONT_SIZE = 12;

    protected DefaultPieDataset dataset;


    public PanelWithPieChart()
    {
        init();
    }

    public PanelWithPieChart(Map<String, Float> nameValueMap, String chartName)
    {
        init(nameValueMap);
        setName(chartName);
    }

    public PanelWithPieChart(List<String> xValues, List<Float> yValues, String chartName)
    {
        String[] xValuesArray = new String[xValues.size()];
        double[] yValuesArray = new double[yValues.size()];
        for (int i=0; i<xValues.size(); i++)
        {
            xValuesArray[i] = xValues.get(i);
            yValuesArray[i] = yValues.get(i);
        }
        init(xValuesArray, yValuesArray);
        setName(chartName);

    }

    public PanelWithPieChart(String[] xValues, int[] yValues,
                             String chartName)
    {
        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<yValues.length; i++)
            yValuesDouble[i] = yValues[i];
        init(xValues, yValuesDouble);
        setName(chartName);

    }

    public PanelWithPieChart(String[] xValues, float[] yValues,
                             String chartName)
    {
        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<yValues.length; i++)
            yValuesDouble[i] = yValues[i];
        init(xValues, yValuesDouble);
        setName(chartName);
    }

    public PanelWithPieChart(String[] xValues, double[] yValues,
                             String chartName)
    {
        init(xValues, yValues);
        setName(chartName);
    }

    protected void init()
    {
        dataset = new DefaultPieDataset();
        JFreeChart chart = ChartFactory.createPieChart(null, dataset, true, true, false);
        
        init(chart.getPlot());
        //set font size to default
        setLabelFontSize(DEFAULT_LABEL_FONT_SIZE);
    }

    protected void init(String[] xValues, double[] yValues)
    {
        init();
        addData(xValues, yValues);
    }

    public void init(Map<String, Float> nameValueMap)
    {
        init();
        addData(nameValueMap);
    }

    public void addData(List<String> names, List<Float> values)
    {
        for (int i=0; i<names.size(); i++)
        {
            addData(names.get(i), values.get(i));
        }
    }

    public void addData(Map<String, Float> nameValueMap)
    {
        for (String name : nameValueMap.keySet())
        {
            addData(name, nameValueMap.get(name));
        }
    }

    public void addData(String[] names, double[] values)
    {
        for (int i=0; i<names.length; i++)
        {
            addData(names[i], values[i]);
        }
    }

    public void addData(String dataName, double dataValue)
    {
        dataset.setValue(dataName, dataValue);
    }

    public void setLabelFontSize(int labelFontSize)
    {
        Font labelFont = ((PiePlot) getPlot()).getLabelFont();
        Font newFont = new Font(labelFont.getName(), labelFont.getStyle(), labelFontSize);
        ((PiePlot) getPlot()).setLabelFont(newFont);
    }


}
