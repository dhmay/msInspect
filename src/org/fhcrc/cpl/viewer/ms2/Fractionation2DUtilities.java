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
package org.fhcrc.cpl.viewer.ms2;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHeatMap;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.LookupPaintScale;

import javax.swing.*;
import java.util.List;
import java.util.ArrayList;
import java.awt.*;


/**
 * Utilities for working with 2D-fractionated MS/MS experiments
 */
public class Fractionation2DUtilities
{
    protected static Logger _log = Logger.getLogger(Fractionation2DUtilities.class);

    public static final int BY_ROW = 0;
    public static final int BY_COLUMN = 1;

    public Fractionation2DUtilities()
    {

    }

    /**
     * Describes the structure of an AMT database containing one or more fractionated
     * experiments
     */
    public static class FractionatedAMTDatabaseStructure
    {
        protected List<FractionatedExperimentStructure> experimentStructures;
        protected List<Integer> beginningIndicesForExperiments;

        public FractionatedAMTDatabaseStructure()
        {
            experimentStructures = new ArrayList<FractionatedExperimentStructure>();
            beginningIndicesForExperiments = new ArrayList<Integer>();
        }

        public FractionatedAMTDatabaseStructure(
                List<FractionatedExperimentStructure> experimentStructures)
        {
            this();
            for (FractionatedExperimentStructure experimentStructure : experimentStructures)
                addExperimentStructure(experimentStructure);
        }

        public FractionatedAMTDatabaseStructure(FractionatedExperimentStructure experimentStructure,
                                                int numExperiments)
        {
            this();
            for (int i=0; i<numExperiments; i++)
                addExperimentStructure(experimentStructure);
        }

        public String toString()
        {
            StringBuffer resultBuf = new StringBuffer("Experiments: " + getNumExperiments()+"\n");
            for (int i=0; i<experimentStructures.size(); i++)
            {
                FractionatedExperimentStructure experimentStruct = experimentStructures.get(i);
                resultBuf.append("Experiment " + i + ": " + experimentStruct + "\n");
            }
            return resultBuf.toString();
        }

        public void addExperimentStructure(FractionatedExperimentStructure experimentStructure)
        {
            beginningIndicesForExperiments.add(getNumFractions());
            experimentStructures.add(experimentStructure);
        }

        public FractionatedExperimentStructure getExperimentStructure(int experimentNumber)
        {
            return experimentStructures.get(experimentNumber);
        }

        public int getNumFractions()
        {
            int numFractions = 0;
            for (FractionatedExperimentStructure experimentStructure : experimentStructures)
            {
                numFractions += experimentStructure.getNumFractions();
            }
            return numFractions;
        }

        /**
         * experiment number is zero-based
         * @param index
         * @return
         */
        public Pair<Integer,int[]> calculateExperimentAndPosition(int index)
        {
            int experimentNumber;
            for (experimentNumber = 0; experimentNumber < getNumExperiments(); experimentNumber++)
            {
                if ((experimentNumber >= getNumExperiments() - 1) ||
                        beginningIndicesForExperiments.get(experimentNumber+1) > index)
                    break;
            }

            int indexInExperiment = index - beginningIndicesForExperiments.get(experimentNumber);
            int[] positionInExperiment =
                    experimentStructures.get(experimentNumber).convertIndexToPosition(indexInExperiment);
            return new Pair<Integer,int[]>(experimentNumber, positionInExperiment);
        }


        public int getNumExperiments()
        {
            return experimentStructures.size();
        }
    }


    public static class FractionatedExperimentStructure
    {
        public int rows;
        public int columns;
        public int organization;

        public FractionatedExperimentStructure(int columns, int rows, int organization)
        {
            this.rows=rows;
            this.columns=columns;
            this.organization=organization;
        }

        public String toString()
        {
            return "cols=" + columns + ", rows=" + rows +
                    ", organization=" + (organization == BY_ROW? "row" : "col");
        }

        public int getNumFractions()
        {
            return rows * columns;
        }

        /**
         * Takes a single index in a 1D array representing the experiment and returns the
         * equivalent 2D location
         * @param index
         * @return a position in two coordinates.  First column, then row
         */
        public int[] convertIndexToPosition(int index)
        {
            int[] result = new int[2];
            switch(organization)
            {
                case BY_ROW:
                    result[0] = index % columns;
                    result[1] = index / columns;
                    break;
                case BY_COLUMN:
                    result[0] = index / rows;
                    result[1] = index % rows;
                    break;
                default:
                    throw new IllegalArgumentException("Illegal organization argument to convertIndexToPosition: " +
                            organization);
            }
            return result;
        }

        /**
         * Takes data representing a value for each index in a 1D array representing the experiment
         * and converts it into a 2D array.
         * @param values
         * @return
         */
        public double[][] convertIndicesToPositions(double[] values)
        {
            double[][] result = new double[columns][rows];
            switch(organization)
            {
                case BY_ROW:
                    for (int i=0; i<rows; i++)
                        for (int j=0; j<columns; j++)
                            result[j][i] = values[(columns * i) + j];
                    break;
                case BY_COLUMN:

                    for (int i=0; i<rows; i++)
                        for (int j=0; j<columns; j++)
                            result[i][j] = values[(columns * i) + j];
                    break;
                default:
                    throw new IllegalArgumentException(
                            "Illegal organization argument to convertIndexToPosition: " +
                                    organization);
            }
            return result;
        }
    }

    /**
     * num rows, num columns, fraction order, number of experiments
     * @param structureArgument
     * @return
     */
    public static FractionatedAMTDatabaseStructure
        digestAmtStructureArgument(String structureArgument)
                        throws ArgumentValidationException
    {
        String[] chunks = structureArgument.split(",");
        if (chunks.length != 4)
            throw new ArgumentValidationException("argument amtdbdimensions must have four comma-separated parts (rows, columns, organization, # experiments)." +
                                                  "  Specified argument has " + chunks.length);
        int lastCommaIndex = structureArgument.lastIndexOf(",");
        FractionatedExperimentStructure experimentStructure =
                digestDimensionArgument(structureArgument.substring(0, lastCommaIndex));

        try
        {
            int numExperiments = Integer.parseInt(chunks[3]);
            return new FractionatedAMTDatabaseStructure(experimentStructure, numExperiments);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }
    }

    /**
     * Takes a String argument defining the structure of an MS/MS experiment
     * (number of rows, number of columns, and whether the runs appear in the order
     * of rows or columns) and parses it
     * @param dimensionArgument
     * @return
     * @throws ArgumentValidationException
     */
    public static FractionatedExperimentStructure
            digestDimensionArgument(String dimensionArgument)
            throws ArgumentValidationException
    {
        int dbRows, dbCols, organization;
      
        String[] chunks = dimensionArgument.split(",");
        if (chunks.length != 3)
            throw new ArgumentValidationException("argument amtdbdimensions must have three comma-separated parts (rows, columns, organization)." +
                                                  "  Specified argument has " + chunks.length);

        try
        {
            dbRows = Integer.parseInt(chunks[0]);
            dbCols = Integer.parseInt(chunks[1]);
        }
        catch (NumberFormatException e)
        {
            throw new ArgumentValidationException("Could not understand the number of rows and columns passed in argument amtdbdimensions");
        }

        if (chunks[2].equalsIgnoreCase("row"))
            organization = BY_ROW;
        else if (chunks[2].equalsIgnoreCase("col"))
            organization = BY_COLUMN;
        else
            throw new ArgumentValidationException("Please specify 'row' or 'col' (without quotes) for the organization of database rows");

        return new FractionatedExperimentStructure(dbCols, dbRows, organization);
    }

    public static void showHeatMapChart(FractionatedAMTDatabaseStructure amtDatabaseStructure,
                                        double[] dataForChart, String chartName, boolean showLegend)
    {
        int chartWidth=1000;
        int chartHeight=1000;

        double globalMinValue = Double.MAX_VALUE;
        double globalMaxValue = Double.MIN_VALUE;

        for (double value : dataForChart)
        {
            if (value < globalMinValue)
                globalMinValue = value;
            if (value > globalMaxValue)
                globalMaxValue = value;
        }

        _log.debug("showHeatMapChart: experiment structures:");
        List<double[][]> amtPeptidesInExperiments = new ArrayList<double[][]>();
        for (int i=0; i<amtDatabaseStructure.getNumExperiments(); i++)
        {
            int experimentWidth = amtDatabaseStructure.getExperimentStructure(i).columns;
            int experimentHeight = amtDatabaseStructure.getExperimentStructure(i).rows;
            _log.debug("\t" + amtDatabaseStructure.getExperimentStructure(i));
            amtPeptidesInExperiments.add(new double[experimentWidth][experimentHeight]);
        }
        for (int i=0; i<dataForChart.length; i++)
        {
            Pair<Integer, int[]> positionInExperiments =
                    amtDatabaseStructure.calculateExperimentAndPosition(i);
            int xpos = positionInExperiments.second[0];
            int ypos = positionInExperiments.second[1];
            int experimentIndex = positionInExperiments.first;
//System.err.println("i, xpos, ypos: " + i + ", " + xpos + ", " + ypos);            
            amtPeptidesInExperiments.get(experimentIndex)[xpos][ypos] =
                              dataForChart[i];

        }

        JDialog cd = new JDialog(ApplicationContext.getFrame(),"Heat Map(s)");
        cd.setSize(chartWidth,chartHeight);
        cd.setPreferredSize(new Dimension(chartWidth,chartHeight));
        cd.setLayout(new FlowLayout());

        int nextPerfectSquareRoot = 0;
        while (true)
        {
            nextPerfectSquareRoot++;
            if ((nextPerfectSquareRoot * nextPerfectSquareRoot) >= amtPeptidesInExperiments.size())
            {
                break;
            }
        }

        int plotWidth = (int) ((double) (chartWidth - 20)/
                (double) nextPerfectSquareRoot);
        int plotHeight = (int) ((double) (chartHeight - 20)/
                (double) nextPerfectSquareRoot);

        _log.debug("Rescaled chart dimensions: " + plotWidth + "x" + plotHeight);


        LookupPaintScale paintScale =
                PanelWithHeatMap.createPaintScale(globalMinValue, globalMaxValue, Color.BLUE, Color.RED);

        for (int i=0; i<amtPeptidesInExperiments.size(); i++)
        {
            PanelWithHeatMap pwhm =
                    new PanelWithHeatMap(amtPeptidesInExperiments.get(i),"Experiment " + (i+1));
//            pwhm.setPalette(PanelWithHeatMap.PALETTE_BLUE_RED);
            pwhm.setPaintScale(paintScale);
            pwhm.setPreferredSize(new Dimension(plotWidth, plotHeight));
            pwhm.setAxisLabels("AX Fraction","RP Fraction");
            if (!showLegend)
                pwhm.getChart().removeLegend();
            cd.add(pwhm);
        }
        cd.setTitle(chartName);
        cd.setVisible(true);
    }

}
