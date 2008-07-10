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

import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.RInterface;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;

/**
 * This class knows how to ask R to plot perspective 3D charts to a file
 * and then display them in a panel
 */
public class PanelWithRPerspectivePlot extends PanelWithBlindImageChart
{
    protected static Logger _log =
            Logger.getLogger(PanelWithRPerspectivePlot.class);

    public static final int DEFAULT_CHART_WIDTH = 700;
    public static final int DEFAULT_CHART_HEIGHT = 700;

    //because the dialog boxes need to be bigger than the image
    public static final int DEFAULT_CHART_DIALOG_WIDTH = DEFAULT_CHART_WIDTH + 10;
    public static final int DEFAULT_CHART_DIALOG_HEIGHT = DEFAULT_CHART_HEIGHT + 25;


    public static final int DEFAULT_ROTATION_ANGLE = 0;
    public static final int DEFAULT_TILT_ANGLE = 30;
    public static final double DEFAULT_SHADE = .5;


    protected int chartWidth = DEFAULT_CHART_WIDTH;
    protected int chartHeight = DEFAULT_CHART_HEIGHT;
    protected int rotationAngle = DEFAULT_ROTATION_ANGLE;
    protected int tiltAngle = DEFAULT_TILT_ANGLE;
    protected double shade = DEFAULT_SHADE;


    public static final String DEFAULT_BACKGROUND_COLOR_STRING = "white";
    public static final String DEFAULT_FOREGROUND_COLOR_STRING = "lightgreen";

    protected String xAxisName = "x";
    protected String yAxisName = "y";
    protected String zAxisName = "z";



    protected boolean showAxisDetails = true;

    protected String backgroundColorString = DEFAULT_BACKGROUND_COLOR_STRING;
    protected String foregroundColorString = DEFAULT_FOREGROUND_COLOR_STRING;

    protected java.util.List<LineVariables> lineVariableList = new ArrayList<LineVariables>();

    protected int millisToWait = 300000;

    public PanelWithRPerspectivePlot()
    {
    }



    public void plotPointsSummary(java.util.List<Float> xValuesList, java.util.List<Float> yValuesList,
                                  double xBinSize, double yBinSize)
    {
        double[] xValues = new double[xValuesList.size()];
        double[] yValues = new double[yValuesList.size()];

        for (int i=0; i<xValues.length; i++)
        {
            xValues[i] = xValuesList.get(i);
            yValues[i] = yValuesList.get(i);
        }                                       
    }

    public void plotPointsSummary(double[] xValues, double[] yValues,
                                  double xBinSize, double yBinSize)
    {
        _log.debug("plotPointsSummary: xBinSize = " + xBinSize + ", yBinSize = " +  yBinSize +
                   ", values: " + xValues.length);

        double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE, maxY = Double.MIN_VALUE;

        for (int i=0; i<xValues.length; i++)
        {
            if (xValues[i] < minX)
                minX = xValues[i];
            if (xValues[i] > maxX)
                maxX = xValues[i];
            if (yValues[i] < minY)
                minY = yValues[i];
            if (yValues[i] > maxY)
                maxY = yValues[i];            
        }

        _log.debug("  minX=" + minX + ", maxX=" + maxX + ", minY=" + minY + ", maxY=" + maxY);

        int numXBins = countBins(minX, maxX, xBinSize);
        int numYBins = countBins(minY, maxY, yBinSize);

        double[] xBins = enumerateBins(minX, xBinSize, numXBins);
        double[] yBins = enumerateBins(minY, yBinSize, numYBins);

        double[][] zMatrix = new double[numXBins][numYBins];
        for (int i=0; i<xValues.length; i++)
        {
            int xBin = calcBin(xValues[i], minX, xBinSize);
            int yBin = calcBin(yValues[i], minY, yBinSize);
            zMatrix[xBin][yBin]++;
        }

        plot(xBins, yBins, zMatrix);
    }

    protected int calcBin(double val, double globalMin, double binSize)
    {
        return (int)((val-globalMin) / binSize);
    }

    protected int countBins(double minValue, double maxValue, double binSize)
    {
        return (int)((maxValue - minValue) / binSize) + 1;
    }

    protected double[] enumerateBins(double minValue, double binSize, int numBins)
    {
        double[] result = new double[numBins];
        _log.debug("enumerateBins, min=" + minValue + ", binSize=" + binSize + ", numBins=" + numBins);        
        for (int i=0; i<numBins; i++)
        {
           result[i] = minValue + (i * binSize);
        }
        return result;
    }

    public void plot(float[] xArray, float[] yArray, float[][] zMatrix)
    {
        double[][] zMatrixDouble = new double[zMatrix.length][zMatrix[0].length];
        for (int i=0; i<zMatrix.length; i++)
        {
            for (int j=0; j<zMatrix[0].length; j++)
                zMatrixDouble[i][j] = zMatrix[i][j]; 
        }
        plot(xArray, yArray, zMatrixDouble);
    }

    public void plot(float[] xArray, float[]yArray, double[][] zMatrix)
    {
        double[]xArrayDouble = new double[xArray.length];
        double[]yArrayDouble = new double[yArray.length];

        for (int i=0; i<xArray.length; i++)
            xArrayDouble[i] = xArray[i];
        for (int i=0; i<yArray.length; i++)
            yArrayDouble[i] = yArray[i];

        plot(xArrayDouble, yArrayDouble, zMatrix);
    }

    public void addLine(double[] xValues, double[] yValues, double[] zValues, String color)
    {
        if (lineVariableList == null)
            lineVariableList = new ArrayList<LineVariables>();
        lineVariableList.add(new LineVariables(xValues, yValues, zValues, color));
    }

    public void plot(double[] xArray, double[]yArray, double[][] zMatrix)
    {
        Map<String, double[]> vectorVarMap = new HashMap<String, double[]>();
        vectorVarMap.put(xAxisName, xArray);
        vectorVarMap.put(yAxisName, yArray);

        Map<String, double[][]> matrixVarMap = new HashMap<String, double[][]>();
        matrixVarMap.put(zAxisName, zMatrix);

        //need to try to create a unique filename
        String outputFileName = "persp" + xArray.length + "" + yArray.length + "" +
                                BasicStatistics.mean(zMatrix[0]) + ".png";


        String tickType = showAxisDetails? "detailed" : "simple";

        String rExpression = "png(\"" + outputFileName +"\"," + chartWidth +
                "," + chartHeight + "); par(bg = \"" +
                backgroundColorString +
                "\",mar=rep(.5,4)); persp(" + xAxisName + "," + yAxisName + "," + zAxisName +
                ",theta=" + rotationAngle +
                ", phi=" + tiltAngle +
                ", shade=" + shade +
                ", col=\"" + foregroundColorString + "\", ticktype=\"" + tickType +"\") -> res;" +
                "round(res,3); ";
//rExpression = rExpression + " x<-c(-.05, .05, .1);";
//rExpression = rExpression + "lines (trans3d(x, y=10, z= 6 + sin(x), pmat = res), col = 3);";
        if (lineVariableList != null)
        {
            int lineNum = 0;
            for (LineVariables lineVariables : lineVariableList)
            {
                _log.debug("Adding line with " + lineVariables.xValues.length + " points");                
                vectorVarMap.put("linex" + lineNum,lineVariables.xValues);
                vectorVarMap.put("liney" + lineNum,lineVariables.yValues);
                vectorVarMap.put("linez" + lineNum,lineVariables.zValues);
                rExpression = rExpression + "lines(trans3d(linex" + lineNum + ", liney" + lineNum +
                                            ", linez" + lineNum + ", pmat = res), col = \"" +
                                            lineVariables.color + "\");";

                lineNum++;
            }
        }


        rExpression = rExpression + " dev.off();";

        RInterface.evaluateRExpression(rExpression,
                vectorVarMap, matrixVarMap, null, millisToWait);
        File pngFile = TempFileManager.createTempFile(outputFileName, this);
        //can't use TempFileManager, since this object may be long gone when the
        //user closes the window
        pngFile.deleteOnExit();

        try
        {
             setImage(pngFile);
        }
        catch (IOException e)
        {
            throw new RuntimeException("Failed to load image from file " + pngFile.getAbsolutePath(),e);
        }
    }


    //be careful here.  These actually have to be valid R variable names
    public void setAxisRVariableNames(String xName, String yName, String zName)
    {
        xAxisName = xName;
        yAxisName = yName;
        zAxisName = zName;
    }


    public int getRotationAngle()
    {
        return rotationAngle;
    }

    public void setRotationAngle(int rotationAngle)
    {
        this.rotationAngle = rotationAngle;
    }


    public int getTiltAngle()
    {
        return tiltAngle;
    }

    public void setTiltAngle(int tiltAngle)
    {
        this.tiltAngle = tiltAngle;
    }


    public double getShade()
    {
        return shade;
    }

    public void setShade(double shade)
    {
        this.shade = shade;
    }

    public int getChartHeight()
    {
        return chartHeight;
    }

    public void setChartHeight(int chartHeight)
    {
        this.chartHeight = chartHeight;
    }

    public int getChartWidth()
    {
        return chartWidth;
    }

    public void setChartWidth(int chartWidth)
    {
        this.chartWidth = chartWidth;
    }

    public String getBackgroundColorString()
    {
        return backgroundColorString;
    }

    public void setBackgroundColorString(String backgroundColorString)
    {
        this.backgroundColorString = backgroundColorString;
    }

    public String getForegroundColorString()
    {
        return foregroundColorString;
    }

    public void setForegroundColorString(String foregroundColorString)
    {
        this.foregroundColorString = foregroundColorString;
    }


    public boolean getShowAxisDetails()
    {
        return showAxisDetails;
    }

    public void setShowAxisDetails(boolean showAxisDetails)
    {
        this.showAxisDetails = showAxisDetails;
    }


    public int getMillisToWait()
    {
        return millisToWait;
    }

    public void setMillisToWait(int millisToWait)
    {
        this.millisToWait = millisToWait;
    }

    public static class LineVariables
    {
        public double[] xValues;
        public double[] yValues;
        public double[] zValues;
        public String color = "black";

        public LineVariables(double[] x, double[] y, double[] z, String color)
        {
            this.xValues = x;
            this.yValues = y;
            this.zValues = z;
            this.color = color;
        }

    }


    
}
