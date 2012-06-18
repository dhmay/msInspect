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

import java.util.List;
import java.util.ArrayList;
import java.awt.*;

/**
 * ChartDialog implementation to make it dead simple to put out scatterplots.
 * Super-convenience class -- if you want to do anything serious with the
 * plots, do it by accessing the PanelWithScatterPlot using
 * getPanelWithScatterPlot()
 */
public class ScatterPlotDialog extends ChartDialog
{
    protected PanelWithScatterPlot _panelWithScatterPlot;


    public ScatterPlotDialog()
    {
        super();
        init();
    }

    public ScatterPlotDialog(List<Double> xValues, List<Double> yValues,
                             String dataSetName)
    {
        this();

        addData(xValues, yValues, dataSetName);
    }

    public ScatterPlotDialog(double[] xValues, double[] yValues,
                             String dataSetName)
    {
        this();
        _panelWithScatterPlot.addData(xValues, yValues, dataSetName);
    }

    public ScatterPlotDialog(float[] xValues, float[] yValues,
                             String dataSetName)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<xValues.length; i++)
        {
            xValuesDouble[i] = xValues[i];
            yValuesDouble[i] = yValues[i];
        }
        init();
        _panelWithScatterPlot.addData(xValuesDouble, yValuesDouble, dataSetName);      
    }

    public ScatterPlotDialog(double[] xValues, double[] yValues,
                             String dataSetName, Shape shape, Color color)
    {
        this();
        _panelWithScatterPlot.addData(xValues, yValues, dataSetName, shape, color);
    }

    protected void init()
    {
        init(true);
    }

    protected void init(boolean showLegend)
    {
        _panelWithScatterPlot = new PanelWithScatterPlot(showLegend);
        _panelWithChart = _panelWithScatterPlot;
        super.init();
    }

    public void setAxisLabels(String xLabel, String yLabel)
    {
        _panelWithScatterPlot.setAxisLabels(xLabel, yLabel);
    }

    public PanelWithScatterPlot getPanelWithScatterPlot()
    {
        return _panelWithScatterPlot;
    }

    public void addData(List<Double> xValues, List<Double> yValues,
                             String dataSetName)
    {
        double[] xValuesArray = new double[xValues.size()];
        double[] yValuesArray = new double[yValues.size()];
        for (int i=0; i<xValues.size(); i++)
        {
            xValuesArray[i] = xValues.get(i);
            yValuesArray[i] = yValues.get(i);
        }
        addData(xValuesArray, yValuesArray, dataSetName);
    }

    public void addDataRedBlueHeatmap(float[] xValues, float[] yValues, float[] zValues,
                                      int numShades)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[xValues.length];
        double[] zValuesDouble = new double[xValues.length];

        for (int i=0; i<xValues.length; i++)
        {
            xValuesDouble[i] = (double) xValues[i];
            yValuesDouble[i] = (double) yValues[i];
            zValuesDouble[i] = (double) zValues[i];

        }

        addDataRedBlueHeatmap(xValuesDouble, yValuesDouble, zValuesDouble, numShades);
    }

    public void addDataRedBlueHeatmap(double[] xValues, double[] yValues, double[] zValues,
                                      int numShades)
    {
        int numPoints = xValues.length;
        double minZ = Double.MAX_VALUE;
        double maxZ = Double.MIN_VALUE;
        for (double zValue : zValues)
        {
            if (zValue < minZ)
                minZ = zValue;
            if (zValue > maxZ)
                maxZ = zValue;
        }

        double zRange = maxZ-minZ;

        for (int j=0; j<numShades; j++)
        {
            double minZValThisGroup = minZ + j * zRange /numShades;
            double maxZValThisGroup = minZ + (j+1) * zRange /numShades;
            int red = (255/numShades) * j;
            int blue = 255 - (255/numShades) * j;
            Color color = new Color(blue, 10,red);

            java.util.List<Double> thisGroupX = new ArrayList<Double>();
            java.util.List<Double> thisGroupY = new ArrayList<Double>();

            for (int k=0; k<numPoints; k++)
            {
                if (zValues[k] <= maxZValThisGroup &&
                        zValues[k] >= minZValThisGroup)
                {
                    thisGroupX.add(xValues[k]);
                    thisGroupY.add(yValues[k]);
//if (Double.isNaN(xValues[k]) || Double.isInfinite(xValues[k]) ||
//        Double.isNaN(yValues[k]) || Double.isInfinite(yValues[k]))System.err.println(xValues[k] + " , " + yValues[k]);
                }
            }
            addData(thisGroupX, thisGroupY, ""+minZValThisGroup);
            _panelWithScatterPlot.setSeriesColor(j, color);
            _panelWithScatterPlot.setPointSize(3);
        }

    }


    public void addData(double[] xValues, double[] yValues, String dataSetName)
    {
        _panelWithScatterPlot.addData(xValues, yValues, dataSetName);
    }

    public void addData(double[] xValues, double[] yValues,
                        String dataSetName, Shape shape, Color color)
    {
        _panelWithScatterPlot.addData(xValues, yValues, dataSetName, shape, color);
    }

    public double[] addRegressionLine()
    {
        return _panelWithScatterPlot.addRegressionLine();
    }

    public void addLine(double slope, double intercept, double minX, double maxX)
    {
        _panelWithScatterPlot.addLine(slope, intercept, minX, maxX);
    }

    public double[] addRegressionLine(boolean robustRegression)
    {
        return _panelWithScatterPlot.addRegressionLine(robustRegression);
    }

    public double[] addRegressionLine(int seriesIndex)
    {
        return _panelWithScatterPlot.addRegressionLine(seriesIndex, false);
    }

    public double[] addRegressionLine(int seriesIndex, boolean robustRegression)
    {
        return _panelWithScatterPlot.addRegressionLine(seriesIndex, robustRegression);
    }


}
