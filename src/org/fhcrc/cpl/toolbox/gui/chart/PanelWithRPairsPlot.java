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

import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.HashMap;
import java.util.List;

/**
 * This class knows how to ask R to plot a "pairs" plot
 */
public class PanelWithRPairsPlot extends PanelWithBlindImageChart
{
    protected static Logger _log =
            Logger.getLogger(PanelWithRPairsPlot.class);

    public static final int DEFAULT_CHART_WIDTH = 700;
    public static final int DEFAULT_CHART_HEIGHT = 700;

    //because the dialog boxes need to be bigger than the image
    public static final int DEFAULT_CHART_DIALOG_WIDTH = DEFAULT_CHART_WIDTH + 10;
    public static final int DEFAULT_CHART_DIALOG_HEIGHT = DEFAULT_CHART_HEIGHT + 25;

    protected double[][] data = null;

    protected int millisToWait = 300000;

    protected int chartWidth = DEFAULT_CHART_WIDTH;
    protected int chartHeight = DEFAULT_CHART_HEIGHT;

    public PanelWithRPairsPlot()
    {
    }

    public void plot(List<List<Double>> dataRowLists)
    {
        plot(dataRowLists, false);
    }

    public void plot(List<List<Double>> dataRowLists, boolean transpose)
    {
        double[][] data = new double[dataRowLists.size()][dataRowLists.get(0).size()];
        for (int row=0; row<dataRowLists.size(); row++)
        {
            List<Double> currentRow = dataRowLists.get(row);
            for (int column=0; column<dataRowLists.get(0).size(); column++)
            {
                data[row][column] = currentRow.get(column);
            }
        }
        plot(data, transpose);
    }

    public void plot(double[][] data)
    {
        plot(data, false);
    }

    public void plot(double[][] data, boolean transpose)
    {
        this.data = data;

        Map<String, double[][]> matrixVarMap = new HashMap<String, double[][]>();
        matrixVarMap.put("data", data);

        StringBuffer rExpressionBuf = new StringBuffer();
        if (transpose)
               rExpressionBuf.append("data <- t(data); ");
        //need to try to create a unique filename
        String outputFileName = "pairs.png";
        File pngFile = TempFileManager.createTempFile(outputFileName, this);
        //can't use TempFileManager, since this object may be long gone when the
        //user closes the window
        pngFile.deleteOnExit();

        rExpressionBuf.append("png(\"" + outputFileName +"\"," + chartWidth +
                "," + chartHeight + "); pairs(data); ");
        rExpressionBuf.append(" dev.off();");

        _log.debug(rExpressionBuf.toString());
        RInterface.evaluateRExpression(rExpressionBuf.toString(),
                null, matrixVarMap, null, millisToWait);

        try
        {
             setImage(pngFile);
        }
        catch (IOException e)
        {
            throw new RuntimeException("Failed to load image from file " + pngFile.getAbsolutePath(),e);
        }
        TempFileManager.deleteTempFiles(this);
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
}
