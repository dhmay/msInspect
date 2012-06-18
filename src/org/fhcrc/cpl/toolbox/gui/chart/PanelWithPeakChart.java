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

import java.awt.*;

/**
 * Triple up the input data.  Around each input data point, put two zero values.
 * The lines will then look like peaks
 */
public class PanelWithPeakChart extends PanelWithLineChart
{
    public PanelWithPeakChart()
    {
        super();
    }

    public PanelWithPeakChart(double[] xValues, double[] yValues,
                             String dataSetName)
    {
        this();
        setName(dataSetName);
        addData(xValues, yValues, dataSetName);
    }

    public PanelWithPeakChart(float[] xValues, float[] yValues,
                             String dataSetName)
    {
        this();
        setName(dataSetName);

        addData(xValues, yValues, dataSetName);
    }

    public PanelWithPeakChart(double[] xValues, double[] yValues,
                             String dataSetName, Shape shape, Color color)
    {
        this();
        setName(dataSetName);
        
        addData(xValues, yValues, dataSetName, shape, color);
    }

    public void addData(double[] xValues, double[] yValues,
                        String dataSetName, Shape shape, Color color)
    {
        double[] paddedXValues = new double[xValues.length*3];
        double[] paddedYValues = new double[xValues.length*3];

        for (int i=0; i<xValues.length; i++)
        {
            paddedXValues[3*i] = (xValues[i]-.000001f);
            paddedYValues[3*i] = (0f);
            paddedXValues[3*i + 1] = (xValues[i]);
            paddedYValues[3*i + 1] =(yValues[i]);
            paddedXValues[3*i + 2] = (xValues[i]+.000001f);
            paddedYValues[3*i + 2] =(0f);
        }
        super.addData(paddedXValues, paddedYValues, dataSetName, shape, color);
    }

}
