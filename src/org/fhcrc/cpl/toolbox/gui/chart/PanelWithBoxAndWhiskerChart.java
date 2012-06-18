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

import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.BoxAndWhiskerToolTipGenerator;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import java.util.List;
import java.util.ArrayList;
import java.awt.*;

/**
 * A simple bar-chart implementation.  Could get more complicated.
 */
public class PanelWithBoxAndWhiskerChart extends PanelWithChart
{
    protected DefaultBoxAndWhiskerCategoryDataset dataset = null;


    public PanelWithBoxAndWhiskerChart(double[] values,
                                       String categoryName,
                                       String chartName)
    {
        setName(chartName);
        init(values, categoryName);
    }

    public PanelWithBoxAndWhiskerChart(List<Float> values,
                                       String categoryName,
                                       String chartName)
    {
        setName(chartName);
        init(values, categoryName);
    }

    public PanelWithBoxAndWhiskerChart(String chartName)
    {
        setName(chartName);
        init();
    }


    protected void init()
    {
        dataset = new DefaultBoxAndWhiskerCategoryDataset();

        final BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
        renderer.setFillBox(false);
        renderer.setToolTipGenerator(new BoxAndWhiskerToolTipGenerator());
        final CategoryAxis xAxis = new CategoryAxis("Type");
        final NumberAxis yAxis = new NumberAxis("Value");
        final CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);

        final JFreeChart chart = new JFreeChart(
                getName(),
                new Font("SansSerif", Font.BOLD, 14),
                plot,
                true
        );

        init(chart.getPlot());
    }

    protected void init(List<Float> values, String categoryName)
    {
        init();
        addData(values, categoryName);
    }

    protected void init(double[] values, String categoryName)
    {
        init();
        addData(values, categoryName);
    }

    public void addData(double[] values, String categoryName)
    {
        List<Double> valuesList = new ArrayList<Double>();
        for (double value : values)
            valuesList.add(value);

        dataset.add(valuesList, categoryName, getName());
    }

    public void addData(List<Float> values, String categoryName)
    {
        List<Double> valuesList = new ArrayList<Double>();
        for (double value : values)
            valuesList.add(value);

        dataset.add(valuesList, categoryName, getName());
    }

}
