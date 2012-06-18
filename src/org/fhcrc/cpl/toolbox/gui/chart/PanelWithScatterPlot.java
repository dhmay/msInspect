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
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYDataItem;
import org.fhcrc.cpl.toolbox.statistics.MatrixUtil;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.apache.log4j.Logger;

import java.awt.geom.Ellipse2D;
import java.awt.*;
import java.util.*;

/**
 * PanelWithChart implementation to make it easy to put out scatterplots.
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithScatterPlot extends PanelWithChart
{
    protected static Logger _log = Logger.getLogger(PanelWithScatterPlot.class);

    protected XYSeriesCollection dataset;
    protected NumberAxis xAxis = new NumberAxis("X");
    protected NumberAxis yAxis = new NumberAxis("Y");


    protected StandardXYItemRenderer renderer = null;

    protected static Shape defaultShape = new Ellipse2D.Double(-1,-1,3,3);

    protected static Color[] SERIES_COLORS = new Color[]{
            Color.BLUE,
            Color.RED,
            Color.GREEN,            
            Color.CYAN,
            Color.MAGENTA,
            Color.ORANGE,
            Color.PINK,
            Color.YELLOW
    };
    protected Color[] seriesColors;


    public PanelWithScatterPlot()
    {
        super();
        init();
    }

    public PanelWithScatterPlot(boolean showLegend)
    {
        super();
        this.showLegend = showLegend;        
        init();
    }

    public PanelWithScatterPlot(java.util.List<? extends Number> xValues, java.util.List<? extends Number> yValues,
                                String dataSetName, String xAxisLabel, String yAxisLabel)
    {
        this(xValues, yValues, dataSetName);
        setAxisLabels(xAxisLabel, yAxisLabel);
    }
    
    public PanelWithScatterPlot(java.util.List<? extends Number> xValues, java.util.List<? extends Number> yValues,
                                String dataSetName)
    {
        this();
        double[] xValuesArray = new double[xValues.size()];
        double[] yValuesArray = new double[yValues.size()];
        for (int i=0; i<xValues.size(); i++)
        {
            xValuesArray[i] = xValues.get(i).doubleValue();
            yValuesArray[i] = yValues.get(i).doubleValue();
        }
        addData(xValuesArray, yValuesArray, dataSetName);
        setName(dataSetName);
    }

    public PanelWithScatterPlot(float[] xValues, float[] yValues,
                             String dataSetName)
    {
        this();
        addData(xValues, yValues, dataSetName);
        setName(dataSetName);
    }

    public PanelWithScatterPlot(double[] xValues, double[] yValues,
                             String dataSetName)
    {
        this();
        setName(dataSetName);
        addData(xValues, yValues, dataSetName);
    }

    public PanelWithScatterPlot(double[] xValues, double[] yValues,
                             String dataSetName, Shape shape, Color color)
    {
        this();
        setName(dataSetName);
        addData(xValues, yValues, dataSetName, shape, color);
    }

    /**
     * Create a scatterplot with the logs of the values passed in
     * @param xValues
     * @param yValues
     * @param dataSetName
     * @return
     */
    public static PanelWithScatterPlot createPlotForLogValues(java.util.List<? extends Number> xValues,
                                                              java.util.List<? extends Number> yValues,
                                                              String dataSetName)
    {
        java.util.List<Double> xValuesLog = new ArrayList<Double>();
        java.util.List<Double> yValuesLog = new ArrayList<Double>();

        for (int i=0; i<xValues.size(); i++)
        {
            xValuesLog.add(Math.log(xValues.get(i).doubleValue()));
            yValuesLog.add(Math.log(yValues.get(i).doubleValue()));
        }
        return new PanelWithScatterPlot(xValuesLog, yValuesLog, dataSetName);
    }

    protected void init()
    {
        seriesColors = new Color[SERIES_COLORS.length];
        System.arraycopy(SERIES_COLORS, 0, seriesColors, 0, SERIES_COLORS.length);
        dataset = new XYSeriesCollection();
        renderer = new StandardXYItemRenderer();
        renderer.setPlotLines(false);
        renderer.setBaseShapesVisible(true);
        renderer.setShapesFilled(true);
        //set all possible series to the default shape
        for (int i=0; i<10; i++)
            renderer.setSeriesShape(i, defaultShape);

        XYPlot scatterPlot = new XYPlot(dataset, xAxis, yAxis, renderer);
        init(scatterPlot);
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

            java.util.List<Float> thisGroupX = new ArrayList<Float>();
            java.util.List<Float> thisGroupY = new ArrayList<Float>();

            for (int k=0; k<numPoints; k++)
            {
                if (zValues[k] <= maxZValThisGroup &&
                        zValues[k] >= minZValThisGroup)
                {
                    thisGroupX.add((float) xValues[k]);
                    thisGroupY.add((float) yValues[k]);
//if (Double.isNaN(xValues[k]) || Double.isInfinite(xValues[k]) ||
//        Double.isNaN(yValues[k]) || Double.isInfinite(yValues[k]))System.err.println(xValues[k] + " , " + yValues[k]);
                }
            }
            addData(thisGroupX, thisGroupY, ""+minZValThisGroup);
//            setSeriesColor(j, color);
//                setPointSize(3);
        }

    }


    public void addData(float[] xValues, float[] yValues, String dataSetName)
    {
        double[] xValuesDouble = new double[xValues.length];
        double[] yValuesDouble = new double[xValues.length];

        for (int i=0; i<xValues.length; i++)
        {
             xValuesDouble[i] = xValues[i];
            yValuesDouble[i] = yValues[i];

        }
        addData(xValuesDouble, yValuesDouble, dataSetName);
    }

    public void addData(java.util.List<? extends Number> xValues, java.util.List<? extends Number> yValues, String dataSetName)
    {
        double[] xArray=new double[xValues.size()];
        double[] yArray=new double[yValues.size()];

        for (int i=0; i<xValues.size(); i++)
        {
            xArray[i] = xValues.get(i).doubleValue();
            yArray[i] = yValues.get(i).doubleValue();
        }
        addData(xArray, yArray, dataSetName, defaultShape,
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
            throw new RuntimeException("PanelWithScatterPlot: x and y values have different cardinality");


        XYSeries series = new XYSeries(dataSetName);

        for (int i=0; i<xValues.length; i++)
        {
            series.add(xValues[i], yValues[i]);
        }
        dataset.addSeries(series);
        setSeriesColor(dataset.getSeriesCount(), color);
    }

    public void addSeries(XYSeries series)
    {
        dataset.addSeries(series);
    }

    public double[] addRegressionLine()
    {
        return addRegressionLine(0, false);
    }

    public double[] addRegressionLine(boolean robustRegression)
    {
        return addRegressionLine(0, false);
    }    

    /**
     * This method only does anything if there's exactly one series in the dataset.
     *
     * Performs linear regression, and then plots a regression line, from the
     * minimum to the maximum X value of the series.
     *
     * Removes the old series, adds this new one, and then adds the old one again,
     * so the regression line will appear on top. 
     */
    public double[] addRegressionLine(int seriesIndex, boolean robustRegression)
    {
        if (dataset == null || dataset.getSeriesCount() < seriesIndex+1)
            return null;
        XYSeries series = dataset.getSeries(seriesIndex);
        int n = series.getItemCount();
        double[] xValues = new double[n];
        double[] yValues = new double[n];

        double maxX = Double.MIN_VALUE;
        double minX = Double.MAX_VALUE;

        for (int i=0; i<n; i++)
        {
            XYDataItem dataItem = series.getDataItem(i);
            xValues[i] = dataItem.getX().doubleValue();
            yValues[i] = dataItem.getY().doubleValue();

            if (xValues[i] > maxX) maxX = xValues[i];
            if (xValues[i] < minX) minX = xValues[i];
        }

        _log.debug("addRegressionLine, minX = " + minX + ", maxX = " + maxX);

        RegressionUtilities.robustRegression(xValues, yValues);
        double[] regressionCoefficients = null;
        if (robustRegression)
            regressionCoefficients = RegressionUtilities.robustRegression(xValues, yValues);
        else
            regressionCoefficients = MatrixUtil.linearRegression(xValues, yValues);
        _log.debug("addRegressionLine, coeffs = " + regressionCoefficients[0] + ", " +
                regressionCoefficients[1]);


        addLine(regressionCoefficients[1], regressionCoefficients[0], minX, maxX);

        return regressionCoefficients;
    }

    /**
     * adds a curve using a polynomial of degree coefficients.length
     * @param coefficients
     * @param minX
     * @param maxX
     */
    public void addLineOrCurve(double[] coefficients, double minX, double maxX)
    {
        int numLineDots = 1000;
        double[] lineXvals = new double[numLineDots];
        double[] lineYvals = new double[numLineDots];

        for (int i= 0; i< numLineDots; i++)
        {
            lineXvals[i] = minX+ (i * (maxX - minX) / numLineDots);
            lineYvals[i] = 0;
            for (int j=0; j<coefficients.length; j++)
                lineYvals[i] += coefficients[j] * Math.pow(lineXvals[i],j);
        }

        _log.debug("addLine, Y vals:  first = " + lineYvals[0] + ", last = " + lineYvals[lineYvals.length-1]);

        addData(lineXvals, lineYvals, "Line or Curve");
    }

    public void addLineOrCurve(double[] coefficients)
    {

        addLineOrCurve(coefficients, dataset.getDomainLowerBound(true), dataset.getDomainUpperBound(true));
    }


    public void addLine(double slope, double intercept, double minX, double maxX)
    {
        addLineOrCurve(new double[] {intercept, slope}, minX, maxX);
    }

    public void addHorizontalLine(double yValue, double minX, double maxX)
    {
        int numLineDots = 1000;
        double[] lineXvals = new double[numLineDots];
        double[] lineYvals = new double[numLineDots];
        for (int i= 0; i< numLineDots; i++)
        {
            lineXvals[i] = minX+ (i * (maxX - minX) / numLineDots);
            lineYvals[i] = yValue;
        }

        _log.debug("addLine, Y vals:  first = " + lineYvals[0] + ", last = " + lineYvals[lineYvals.length-1]);

        addData(lineXvals, lineYvals, "Line");
    }

    public void addVerticalLine(double xValue, double minY, double maxY)
    {
        int numLineDots = 1000;
        double[] lineXvals = new double[numLineDots];
        double[] lineYvals = new double[numLineDots];
        for (int i= 0; i< numLineDots; i++)
        {
            lineXvals[i] = xValue;
            lineYvals[i] =  minY+ (i * (maxY - minY) / numLineDots);
        }

        _log.debug("addLine, Y vals:  first = " + lineYvals[0] + ", last = " + lineYvals[lineYvals.length-1]);

        addData(lineXvals, lineYvals, "Line");
    }

    public void setPointSize(int pointSize)
    {
        int upLeft = - (pointSize / 2);

        renderer.setShape(new Ellipse2D.Double(-upLeft,-upLeft,pointSize,pointSize));
    }

    public void addCrosshairs(double initialDomainCrosshairValue, double initialRangeCrosshairValue)
    {
        ((XYPlot)getPlot()).setDomainCrosshairVisible(true);
        ((XYPlot)getPlot()).setRangeCrosshairVisible(true);
        ((XYPlot)getPlot()).setDomainCrosshairValue(initialDomainCrosshairValue);
        ((XYPlot)getPlot()).setRangeCrosshairValue(initialRangeCrosshairValue);
    }

    public void addCrosshairsAndListener(CrosshairChangeListener crosshairListener,
                                         double initialDomainCrosshairValue, double initialRangeCrosshairValue)
    {
        addCrosshairs(initialDomainCrosshairValue, initialRangeCrosshairValue);
        getChart().addProgressListener(crosshairListener);
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

    public Color[] getSeriesColors() {
        return seriesColors;
    }

    public void setSeriesColors(Color[] seriesColors) {
        this.seriesColors = seriesColors;
    }
}
