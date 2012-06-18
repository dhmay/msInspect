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

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYDataset;

import java.awt.geom.Ellipse2D;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * PanelWithChart implementation to make it easy to put out Line Charts
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithLineChart extends PanelWithChart
{
    protected XYSeriesCollection dataset;
    protected NumberAxis xAxis = new NumberAxis("X");
    protected NumberAxis yAxis = new NumberAxis("Y");

    protected StandardXYItemRenderer renderer = null;

    public static Shape defaultShape = new Ellipse2D.Double(-1,-1,3,3);

    protected static Color[] seriesColors = new Color[]{
            Color.BLUE,
            Color.RED,
            Color.CYAN,
            Color.ORANGE,            
            Color.MAGENTA,
            Color.PINK,
            Color.GREEN,
            Color.YELLOW
    };

    public PanelWithLineChart()
    {
        super();
        init();
    }

    public PanelWithLineChart(java.util.List<? extends Number> xValuesList, java.util.List<? extends Number> yValuesList,
                             String dataSetName)
    {
        this();

        addData(xValuesList, yValuesList, dataSetName);
        setName(dataSetName);

    }
    
    public PanelWithLineChart(double[] xValues, double[] yValues,
                             String dataSetName)
    {
        this();
        addData(xValues, yValues, dataSetName);
        setName(dataSetName);

    }

    public PanelWithLineChart(int[] xValues, int[] yValues,
                             String dataSetName)
    {
        this();
        addData(xValues, yValues, dataSetName);
        setName(dataSetName);

    }

    public PanelWithLineChart(float[] xValues, float[] yValues,
                             String dataSetName)
    {
        this();
        addData(xValues, yValues, dataSetName);
        setName(dataSetName);
    }

    public PanelWithLineChart(double[] xValues, double[] yValues,
                             String dataSetName, Shape shape, Color color)
    {
        this();
        addData(xValues, yValues, dataSetName, shape, color);
        setName(dataSetName);

    }

    protected void init()
    {
        dataset = new XYSeriesCollection();
        renderer = new StandardXYItemRenderer();
        renderer.setPlotLines(true);
        renderer.setBaseShapesVisible(false);

        //set all possible series to the default shape
        for (int i=0; i<10; i++)
            renderer.setSeriesShape(i, defaultShape);



        XYToolTipGenerator toolTipGenerator = new XYToolTipGenerator()
        {
            public String generateToolTip(XYDataset xyDataset, int s, int i)
            {
                double X = Math.round(xyDataset.getXValue(s, i) * 1000.0) / 1000.0;
                double Y = Math.round(xyDataset.getYValue(s, i) * 1000.0) / 1000.0;

                return "(" + X + ", " + Y + ")";
            }
        };


        //dhmay adding for jfreechart 1.0.6 upgrade.  If this isn't here, we get a
        //nullPointerException in XYBarRenderer.drawItemLabel
        renderer.setBaseItemLabelGenerator(new NullLabelGenerator());

        renderer.setSeriesItemLabelsVisible(0, true);
        renderer.setBaseToolTipGenerator(toolTipGenerator);

        _chart = ChartFactory.createXYLineChart(null, null, null, dataset,
                    PlotOrientation.VERTICAL, true, true, false);

        _chart.getXYPlot().setRenderer(renderer);
        init(_chart);
    }

    /**
     * This should not be necessary.  dhmay creating this class for the jfreechart 1.0.6 upgrade.
     * We shouldn't have to specify a labelgenerator at all for the barchart
     *
     */
    protected static class NullLabelGenerator implements XYItemLabelGenerator
    {
        public String generateLabel(XYDataset dataset,
                                    int series,
                                    int item)
        {
            return null;
        }
    }




    public void setAxisLabels(String xLabel, String yLabel)
    {
        xAxis.setLabel(xLabel);
        yAxis.setLabel(yLabel);
    }

    protected Color getDefaultColorForSeries(int index)
    {
        if (index >= seriesColors.length)
            index = index % seriesColors.length;
        return seriesColors[index];
    }

    public void addData(java.util.List<? extends Number> xValuesList, java.util.List<? extends Number> yValuesList,
                        String dataSetName)
    {
        double[] xValues = new double[xValuesList.size()];
        double[] yValues = new double[yValuesList.size()];
        for (int i=0; i<xValues.length; i++)
        {
            xValues[i] = xValuesList.get(i).doubleValue();
            yValues[i] = yValuesList.get(i).doubleValue();
        }

        addData(xValues, yValues, dataSetName);

    }

    public void addData(int[] xValues, int[] yValues, String dataSetName)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[yValues.length];

        for (int i=0; i<xValues.length; i++)
        {
            xValuesDouble[i] = xValues[i];
            yValuesDouble[i] = yValues[i];
        }
        addData(xValuesDouble, yValuesDouble, dataSetName, defaultShape,
                getDefaultColorForSeries(dataset.getSeriesCount()));
    }

    public void addData(float[] xValues, float[] yValues, String dataSetName)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[yValues.length];

        for (int i=0; i<xValues.length; i++)
        {
            xValuesDouble[i] = xValues[i];
            yValuesDouble[i] = yValues[i];
        }
        addData(xValuesDouble, yValuesDouble, dataSetName, defaultShape,
                getDefaultColorForSeries(dataset.getSeriesCount()));
    }

    public void addData(double[] xValues, double[] yValues, String dataSetName)
    {
        addData(xValues, yValues, dataSetName, defaultShape,
                getDefaultColorForSeries(dataset.getSeriesCount()));
    }

    public void addData(double[] xValues, double[] yValues,
                        String dataSetName, Shape shape, Color color)
    {
        if (xValues.length != yValues.length)
            throw new RuntimeException("PanelWithLineChart: x and y values have different cardinality (" +
                    -xValues.length + " and " + yValues.length + ")");


        XYSeries series = new XYSeries(dataSetName);

        for (int i=0; i<xValues.length; i++)
        {
            series.add(xValues[i], yValues[i]);
        }
        dataset.addSeries(series);
        setSeriesColor(dataset.getSeriesCount()-1, color);
    }

    public void addDataFloat(float[] xValues, float[] yValues,
                        String dataSetName, Shape shape, Color color)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[yValues.length];

        for (int i=0; i<xValues.length; i++)
        {
            xValuesDouble[i] = xValues[i];
            yValuesDouble[i] = yValues[i];
        }
        addData(xValuesDouble, yValuesDouble, dataSetName, shape, color);
    }

    public void addSeries(XYSeries series)
    {
        dataset.addSeries(series);
    }

    public void setSeriesColor(int i, Color color)
    {
        renderer.setSeriesPaint(i, color);
    }

    XYSeries getDataSeries(int i)
    {
        return dataset.getSeries(i);
    }

    public StandardXYItemRenderer getRenderer()
    {
        return renderer;
    }

    public NumberAxis getXAxis()
    {
        return xAxis;
    }

    public NumberAxis getYAxis()
    {
        return yAxis;
    }
}
