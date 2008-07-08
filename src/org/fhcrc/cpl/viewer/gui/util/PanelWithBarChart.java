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

package org.fhcrc.cpl.viewer.gui.util;

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.data.category.DefaultCategoryDataset;

/**
 * A simple bar-chart implementation.  Could get more complicated.
 */
public class PanelWithBarChart extends PanelWithChart
{

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

    protected void init(String[] xValues, double[] yValues, String categoryName)
    {

        if (xValues.length != yValues.length)
            throw new RuntimeException("PanelWithBarChart: x and y values have different cardinality");

        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        JFreeChart chart = ChartFactory.createBarChart(null, null, null, dataset,
                    PlotOrientation.VERTICAL, true, false, false);

        init(chart.getPlot());


        for (int i=0; i<xValues.length; i++)
        {
            dataset.addValue(yValues[i], xValues[i], categoryName);
        }

    }

}
