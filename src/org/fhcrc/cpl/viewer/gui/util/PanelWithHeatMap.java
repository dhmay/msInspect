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

package org.fhcrc.cpl.viewer.gui.util;

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.XYZToolTipGenerator;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.*;
import org.jfree.ui.RectangleAnchor;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;

import java.awt.*;

/**
 * PanelWithChart implementation to make it easy to put out Line Charts
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithHeatMap extends PanelWithChart
{
    static Logger _log = Logger.getLogger(PanelWithHeatMap.class);

//this doesn't work -- bug in XYZBlockRenderer
//    protected XYZToolTipGenerator toolTipGenerator;

    protected XYDataset dataset;
    protected NumberAxis xAxis = new NumberAxis("X");
    protected NumberAxis yAxis = new NumberAxis("Y");

    protected String dataSetName = "";

    LookupPaintScale paintScale = null;

    protected XYBlockRenderer renderer = null;

    protected double[] xValues = null;
    protected double[] yValues = null;
    protected double[][] zValues = null;

    protected double lowerZBound = 0;
    protected double upperZBound = 0;

    public static final int PALETTE_BLUE_RED = 0;
    public static final int PALETTE_RED = 1;
    public static final int PALETTE_BLUE = 2;
    public static final int PALETTE_GRAYSCALE = 3;
    public static final int PALETTE_WHITE_BLUE = 4;
    //This palette starts with blue at 0 and scales through the values.  Any negative values get the same red color
    public static final int PALETTE_POSITIVE_WHITE_BLUE_NEGATIVE_RED = 5;



    public static final int DEFAULT_PALETTE = PALETTE_BLUE_RED;

    protected int palette = DEFAULT_PALETTE;

    public static final int DISTINCT_PALETTE_VALUES = 100;


    public PanelWithHeatMap()
    {
        super();
        init();
    }

    public PanelWithHeatMap(String dataSetName)
    {
        super(dataSetName);
        this.dataSetName = dataSetName;
    }




    public PanelWithHeatMap(float[] xValues, float[] yValues, float[][] zValues,
                             String dataSetName)
    {
        this(dataSetName);
        setData(xValues, yValues, zValues);
    }

    public PanelWithHeatMap(double[] xValues, double[] yValues, double[][] zValues,
                             String dataSetName)
    {
        this(dataSetName);
        setData(xValues, yValues, zValues);
    }

    public PanelWithHeatMap(float[][] zValues,
                             String dataSetName)
    {
        this(dataSetName);
        double[][] zValuesDouble = new double[zValues.length][zValues[0].length];
        for (int i=0; i<zValues.length; i++)
            for (int j=0; j<zValues[0].length; j++)
            {
                zValuesDouble[i][j] = zValues[i][j];
            }
        setData(zValuesDouble);
    }

    public PanelWithHeatMap(double[][] zValues,
                            String dataSetName)
    {
        this(dataSetName);
        setData(zValues);
    }

    public PanelWithHeatMap(double[] zValues,
                             int width, int height,
                             int organization,
                             String dataSetName)
    {
        this(dataSetName);
        setData(zValues, width, height, organization);
    }

    public PanelWithHeatMap(double[] zValues,
                             Fractionation2DUtilities.FractionatedExperimentStructure expStructure,
                             String dataSetName)
    {
        this(dataSetName);
        setData(zValues, expStructure.columns, expStructure.rows, expStructure.organization);
    }

    public PanelWithHeatMap(java.util.List<Double> zValues,
                             Fractionation2DUtilities.FractionatedExperimentStructure expStructure,
                             String dataSetName)
    {
        this(dataSetName);
        setData(zValues, expStructure.rows, expStructure.columns, expStructure.organization);
    }



    protected void init()
    {
        dataset = new XYSeriesCollection();
        renderer = new XYBlockRenderer();
        //set all possible series to the default shape


//        init(chart.getXYPlot());
    }

    public void setAxisLabels(String xLabel, String yLabel)
    {
        xAxis.setLabel(xLabel);
        yAxis.setLabel(yLabel);
    }

    public void setData(float[] xValues, float[] yValues, float[][] zValues)
    {
        double[] xValuesDouble = new double[xValues.length];
        for (int i=0; i<xValues.length; i++)
            xValuesDouble[i] = xValues[i];

        double[] yValuesDouble = new double[yValues.length];
        for (int i=0; i<yValues.length; i++)
            yValuesDouble[i] = yValues[i];

        double[][] zValuesDouble = new double[zValues.length][zValues[0].length];
        for (int i=0; i<zValues.length; i++)
            for (int j=0; j<zValues[0].length; j++)
            {
                zValuesDouble[i][j] = zValues[i][j];
            }

        setData(xValuesDouble, yValuesDouble, zValuesDouble);
    }

    public void setData(double[] zValues, int width, int height, int organization)
    {
        Fractionation2DUtilities.FractionatedExperimentStructure expStructure =
                new Fractionation2DUtilities.FractionatedExperimentStructure(width, height, organization);
        double[][] zValuesRedone = expStructure.convertIndicesToPositions(zValues);
        setData(zValuesRedone);
    }

    public void setData(java.util.List<Double> zValues, int width, int height, int organization)
    {
        setData(convertDoubleListToDoubleArray(zValues), width, height, organization);
    }

    public static double[] convertDoubleListToDoubleArray(java.util.List<Double> doubleList)
    {
        double[] result = new double[doubleList.size()];
        for (int i=0; i<doubleList.size(); i++)
            result[i] = doubleList.get(i);
        return result;
    }


    public void setData(double[][] zValues)
    {
        int width = zValues.length;
        int height = zValues[0].length;

        double[] xValues = new double[width];
        for (int i=0; i<width; i++)
            xValues[i] = i+1;

        double[] yValues = new double[height];
        for (int i=0; i<height; i++)
            yValues[i] = i+1;
        _log.debug("Inferring width=" + width + ", height=" + height);
        if (_log.isDebugEnabled())
        {
            for (int i=0; i<width; i++)
            {
                System.err.println("\n");
                for (int j=0; j<height; j++)
                {
                    System.err.print(zValues[i][j] + "\t");
                }
            }
            System.err.println();
        }
        setData(xValues, yValues, zValues);
    }

    public void setData(double[] xValues, double[] yValues, double[][] zValues)
    {
        this.xValues = xValues;
        this.yValues = yValues;
        this.zValues = zValues;

        double minZValue = Double.MAX_VALUE;
        double maxZValue = Double.MIN_VALUE;
        int width = xValues.length;
        int height = yValues.length;
        int numCells = width * height;
        _log.debug("Number of cells in heat map: " + numCells);
        if (zValues.length != width || zValues[0].length != height)
            throw new RuntimeException("PanelWithHeatMap: wrong number of z values for x and y values (" +
                    zValues.length + " vs. " + width + ", " + zValues[0].length + " vs. " + height +
                    ", x/y first, z second)");
        DefaultXYZDataset theDataset = new DefaultXYZDataset();
        double[][] data = new double[3][numCells];
        for(int j=0; j<height; j++){
            for(int i=0; i<width; i++)
            {
                int cellIndex = (j * width) + i;
                data[0][cellIndex]= xValues[i];
                data[1][cellIndex]= yValues[j];
                data[2][cellIndex]= zValues[i][j];
                //keep track of lowest/highest z values
                minZValue = Math.min(zValues[i][j], minZValue);
                maxZValue = Math.max(zValues[i][j], maxZValue);
            }
        }
        lowerZBound = Rounder.round(minZValue,3);
        upperZBound = Rounder.round(maxZValue,3);
        if (lowerZBound == upperZBound)
            upperZBound += .0001;
        _log.debug("low,high values: " + lowerZBound + ", " + upperZBound);
        theDataset.addSeries("Range: " + lowerZBound + "-" + upperZBound,data);

        dataset = theDataset;
        XYBlockRenderer thisRenderer = new XYBlockRenderer();
        renderer = thisRenderer;

        PaintScale paintScale = createPaintScale(palette);
        setPaintScale(paintScale);
        //This is necessary to get everything to line up
        renderer.setBlockAnchor(RectangleAnchor.BOTTOM);

//        if (toolTipGenerator != null)
//        {
//System.err.println("****Setting TTG");
//            renderer.setBaseToolTipGenerator(new StandardXYZToolTipGenerator());
//                    //toolTipGenerator);
//        }

        if (getPlot() != null)
        {
            ((XYPlot) getPlot()).setDataset(dataset);
            ((XYPlot) getPlot()).setRenderer(renderer);

            invalidate();
            return;
        }
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

//        plot.setRenderer(renderer);

//        setPaintScale(createGrayPaintScale(minZValue, maxZValue));


//        setBlockWidthHeight(1f,1f);


        JFreeChart chart = new JFreeChart(dataSetName,JFreeChart.DEFAULT_TITLE_FONT,plot,true);
            
//        chart.addLegend(new LegendTitle(renderer));
//        PaintScaleLegend legend = new PaintScaleLegend(paintScale, xAxis);
//        chart.addLegend(legend);

        
//        LegendItemCollection c1 = new LegendItemCollection();
//
//        LegendItem item1 = new LegendItem("Label", "Description",
//                "ToolTip", "URL", true,
//                new Rectangle2D.Double(1.0, 2.0, 3.0, 4.0), true, Color.red,
//                true, Color.blue, new BasicStroke(1.2f), true,
//                new Line2D.Double(1.0, 2.0, 3.0, 4.0),
//                new BasicStroke(2.1f), Color.green);
//        LegendItem item2 = new LegendItem("Label", "Description",
//                "ToolTip", "URL", true,
//                new Rectangle2D.Double(1.0, 2.0, 3.0, 4.0),
//                true, Color.red, true, Color.blue, new BasicStroke(1.2f), true,
//                new Line2D.Double(1.0, 2.0, 3.0, 4.0), new BasicStroke(2.1f),
//                Color.green);
//        c1.add(item1);
//
//        chart.getLegend().setSources(new LegendItemSource[]{renderer});

        init(chart);
    }

    public PaintScale createPaintScale(int palette)
    {
        PaintScale result;
        switch (palette)
        {
            case PALETTE_BLUE_RED:
                result = createPaintScale(lowerZBound, upperZBound, Color.BLUE, Color.RED);
                break;
            case PALETTE_RED:
                result = createPaintScale(lowerZBound, upperZBound, new Color(70,5,5), new Color(255,5,5));
                break;
            case PALETTE_BLUE:
                result = createPaintScale(lowerZBound, upperZBound, new Color(5,5,70), new Color(5,5,255));
                break;
            case PALETTE_GRAYSCALE:
                result = createPaintScale(lowerZBound, upperZBound, Color.WHITE, new Color(5,5,5));
                break;
            case PALETTE_WHITE_BLUE:
                result = createPaintScale(lowerZBound, upperZBound, Color.WHITE, Color.BLUE);
                break;
            case PALETTE_POSITIVE_WHITE_BLUE_NEGATIVE_RED:
                result = new LookupPaintScale(lowerZBound, upperZBound+0.1, Color.RED);
                addValuesToPaintScale((LookupPaintScale) result, 0, upperZBound, Color.WHITE, Color.BLUE);
                break;
            default:
                result = createPaintScale(lowerZBound, upperZBound, Color.WHITE, new Color(5,5,5));
                break;
        }
        return result;
    }

    protected static void addValuesToPaintScale(LookupPaintScale paintScale, double lowerBound, double upperBound,
                                                Color lowColor, Color highColor)
    {
        int distinctValues = DISTINCT_PALETTE_VALUES;

        if (upperBound <= lowerBound)
            upperBound = lowerBound + .0001;
        double increment = (upperBound - lowerBound) / distinctValues;

        int redDiff = highColor.getRed() - lowColor.getRed();
        int greenDiff = highColor.getGreen() - lowColor.getGreen();
        int blueDiff = highColor.getBlue() - lowColor.getBlue();
        double redIncrement = (redDiff / distinctValues);
        double greenIncrement = (greenDiff / distinctValues);
        double blueIncrement = (blueDiff / distinctValues);

        _log.debug("Palette: ");

        for (int i=0; i<distinctValues; i++)
        {
            int r = (int) (lowColor.getRed() + (i * redIncrement));
            int g = (int) (lowColor.getGreen() + (i * greenIncrement));
            int b = (int) (lowColor.getBlue() + (i * blueIncrement));
            Color incrementColor = new Color(r,g,b);
            double incrementStart = lowerBound + (i * increment);
            paintScale.add(incrementStart, incrementColor);
            _log.debug("\t" + incrementStart + "-" + (incrementStart + increment) + ": " + incrementColor);
        }
    }

    /**
     * Creates a paintscale that has _distinctValues_ even increments between
     * lowColor and highColor, starting with value lowerBound and ending with upperBound
     * @param lowerBound
     * @param upperBound
     * @param lowColor
     * @param highColor
     * @return
     */
    public static LookupPaintScale createPaintScale(double lowerBound,
                                       double upperBound,
                                       Color lowColor, Color highColor)
    {
        //prevent rounding errors that make highest value undefine
        LookupPaintScale result = new LookupPaintScale(lowerBound, upperBound+0.01, lowColor);
        addValuesToPaintScale(result, lowerBound, upperBound, lowColor, highColor);
        return result;
    }







    public void setBlockWidthHeight(float blockWidth, float blockHeight)
    {
        renderer.setBlockWidth(blockWidth);
        renderer.setBlockHeight(blockHeight);
        repaint();
    }

    public void setPaintScale(PaintScale newPaintScale)
    {
        renderer.setPaintScale(newPaintScale);
        repaint();
    }

    public void setPalette(int palette)
    {
        this.palette = palette;
        setPaintScale(createPaintScale(palette));
        repaint();
    }

    public XYItemRenderer getRenderer()
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


    public double[] getXValues()
    {
        return xValues;
    }

    public void setXValues(double[] xValues)
    {
        this.xValues = xValues;
        setData(xValues, yValues, zValues);
    }

    public double[] getYValues()
    {
        return yValues;
    }

    public void setYValues(double[] yValues)
    {
        this.yValues = yValues;
        setData(xValues, yValues, zValues);
    }

    public double[][] getZValues()
    {
        return zValues;
    }

    public void setZValues(double[][] zValues)
    {
        this.zValues = zValues;
        setData(xValues, yValues, zValues);
    }
}
