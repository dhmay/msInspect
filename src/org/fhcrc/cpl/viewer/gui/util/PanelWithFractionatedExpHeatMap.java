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

package org.fhcrc.cpl.viewer.gui.util;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHeatMap;


/**
 * PanelWithHeatMap subclass tailored to plotting results from a 2-dimensional fractionation experiment, one
 * value per fraction.
 *
 * These methods were formerly part of PanelWithHeatMap
 */
public class PanelWithFractionatedExpHeatMap extends PanelWithHeatMap
{
    static Logger _log = Logger.getLogger(PanelWithFractionatedExpHeatMap.class);

    public PanelWithFractionatedExpHeatMap(double[] zValues,
                             int width, int height,
                             int organization,
                             String dataSetName)
    {
        super(dataSetName);
        setData(zValues, width, height, organization);
    }

    public PanelWithFractionatedExpHeatMap(double[] zValues,
                             Fractionation2DUtilities.FractionatedExperimentStructure expStructure,
                             String dataSetName)
    {
        super(dataSetName);
        setData(zValues, expStructure.columns, expStructure.rows, expStructure.organization);
    }

    public PanelWithFractionatedExpHeatMap(java.util.List<Double> zValues,
                             Fractionation2DUtilities.FractionatedExperimentStructure expStructure,
                             String dataSetName)
    {
        super(dataSetName);
        setData(zValues, expStructure.rows, expStructure.columns, expStructure.organization);
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
}
