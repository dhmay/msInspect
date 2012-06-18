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
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.data.category.DefaultCategoryDataset;

import java.util.List;
import java.util.Map;

/**
 * A simple bar-chart implementation.  Could get more complicated.
 */
public class PanelWithBarChart extends PanelWithChart
{
    protected DefaultCategoryDataset dataset;

    public PanelWithBarChart()
    {
        init();
    }

    public PanelWithBarChart(Map<String, Float> nameValMap, String categoryName)
    {
        String[] xValuesArray = new String[nameValMap.size()];
        double[] yValuesArray = new double[nameValMap.size()];
        int i=0;
        for (String name : nameValMap.keySet())
        {
            xValuesArray[i] = name;
            yValuesArray[i++] = nameValMap.get(name);
        }
        init(xValuesArray, yValuesArray, categoryName);
    }

    public PanelWithBarChart(List<String> xValues, List<Float> yValues, String categoryName)
    {
        String[] xValuesArray = new String[xValues.size()];
        double[] yValuesArray = new double[yValues.size()];
        for (int i=0; i<xValues.size(); i++)
        {
            xValuesArray[i] = xValues.get(i);
            yValuesArray[i] = yValues.get(i);
        }
        init(xValuesArray, yValuesArray, categoryName);
    }

    public PanelWithBarChart(String[] xValues, int[] yValues,
                             String categoryName)
    {
        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<yValues.length; i++)
            yValuesDouble[i] = yValues[i];
        init(xValues, yValuesDouble, categoryName);
    }

    public PanelWithBarChart(String[] xValues, float[] yValues,
                             String categoryName)
    {
        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<yValues.length; i++)
            yValuesDouble[i] = yValues[i];
        init(xValues, yValuesDouble, categoryName);
    }

    public PanelWithBarChart(String[] xValues, double[] yValues,
                             String categoryName)
    {
        init(xValues, yValues, categoryName);
    }

    protected void init()
    {
        dataset = new DefaultCategoryDataset();
        JFreeChart chart = ChartFactory.createBarChart(null, null, null, dataset,
                    PlotOrientation.VERTICAL, true, false, false);

        init(chart.getPlot());
    }

    protected void init(String[] xValues, double[] yValues, String categoryName)
    {
        init();
        setName(categoryName);
        addData(xValues, yValues, categoryName);
    }

    public void addData(String[] xValues, double[] yValues, String categoryName)
    {

        if (xValues.length != yValues.length)
            throw new RuntimeException("PanelWithBarChart: x and y values have different cardinality");
        for (int i=0; i<xValues.length; i++)
        {
            dataset.addValue(yValues[i], xValues[i], categoryName);
        }
    }

    /**
     * Angle the labels on the category (x) axis
     * @param angle in radians
     */
    public void setCategoryLabelAngle(double angle)
    {
        ((CategoryPlot) getPlot()).getDomainAxis().setCategoryLabelPositions(
                    CategoryLabelPositions.createUpRotationLabelPositions(angle));
    }

    /**
     * Angle the labels on the category (x) axis to a nice, jaunty 30 degrees
     */
    public void setCategoryLabelAngle30Degrees()
    {
        setCategoryLabelAngle(Math.PI / 6.0);
    }

}
