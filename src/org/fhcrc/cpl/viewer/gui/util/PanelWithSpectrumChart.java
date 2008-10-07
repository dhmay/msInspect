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

import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.labels.XYZToolTipGenerator;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.xy.XYDataset;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.toolbox.FloatRange;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;

import java.util.Map;
import java.util.HashMap;

/**
 * PanelWithChart implementation to make it easy to put out Line Charts
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithSpectrumChart extends PanelWithHeatMap
{
    static Logger _log = Logger.getLogger(PanelWithSpectrumChart.class);

    protected MSRun run;

    protected int minScan = 0;
    protected int maxScan = 0;

    protected float minMz = 0;
    protected float maxMz = 0;

    protected int resolution = DEFAULT_RESOLUTION;

    public static final int DEFAULT_RESOLUTION = 100;

    protected int scanLine1 = 0;
    protected int scanLine2 = 0;

    protected boolean generateLineCharts = false;

    protected Map<Integer, PanelWithLineChart> scanLineChartMap = null;


    public PanelWithSpectrumChart()
    {
        super();
    }

    public PanelWithSpectrumChart(MSRun run, int minScan, int maxScan, float minMz, float maxMz, int resolution,
                                  int scanLine1, int scanLine2, boolean generateLineCharts)
    {
        this();

        this.run = run;
        this.minScan = minScan;
        this.maxScan = maxScan;
        this.minMz = minMz;
        this.maxMz = maxMz;
        this.resolution = resolution;
        this.scanLine1 = scanLine1;
        this.scanLine2 = scanLine2;
        this.generateLineCharts = generateLineCharts;

        this.setPalette(PanelWithHeatMap.PALETTE_POSITIVE_WHITE_BLUE_NEGATIVE_RED);
        setAxisLabels("scan", "m/z");

        scanLineChartMap = new HashMap<Integer, PanelWithLineChart>();

        generateChart();
    }

    protected void generateChart()
    {
        int minScanIndex = Math.abs(run.getIndexForScanNum(minScan));
        int maxScanIndex = Math.abs(run.getIndexForScanNum(maxScan));
        int numScans = maxScanIndex - minScanIndex + 1;

        int numMzBins = (int) ((float) resolution * (maxMz - minMz)) + 1;        

        double[] scanValues = new double[numScans];
        double[] scanIndexValues = new double[numScans];

        double[] mzValues = new double[numMzBins];

        for (int i=0; i<numScans; i++)
        {
            scanIndexValues[i] = minScanIndex + i;
            scanValues[i] = run.getScanNumForIndex(minScanIndex + i);
        }
        for (int i=0; i<numMzBins; i++)
            mzValues[i] = minMz + (i / (float) resolution);
        //Resampling adds some slop on each end
        FloatRange mzWindowForResample = new FloatRange(minMz, maxMz);
        double[][] intensityValues = new double[numScans][numMzBins];
        _log.debug("Loading spectrum in range....");



        double maxIntensityOnChart = 0;

        for (int scanArrayIndex = 0; scanArrayIndex < numScans; scanArrayIndex++)
        {
            int scanIndex = minScanIndex + scanArrayIndex;
            MSRun.MSScan scan = run.getScan(scanIndex);

            float[] signal = Spectrum.Resample(scan.getSpectrum(), mzWindowForResample, resolution);
            //todo: this is horrible
            double[] signalAsDouble = new double[signal.length];


            for (int i=0; i<signal.length; i++)
                signalAsDouble[i] = signal[i];

            if (generateLineCharts)
            {
                PanelWithLineChart lineChart =
                        new PanelWithPeakChart(mzValues, signalAsDouble, "Scan " + (int) scanValues[scanArrayIndex]);
                lineChart.setAxisLabels("m/z", "intensity");
                scanLineChartMap.put((int) scanValues[scanArrayIndex], lineChart);
            }
            
            intensityValues[scanArrayIndex] = signalAsDouble;

            
            for (int i=0; i<signal.length; i++)
            {
                if (signalAsDouble[i] > maxIntensityOnChart)
                    maxIntensityOnChart = signalAsDouble[i];
            }


        }
        ApplicationContext.infoMessage("Maximum intensity on chart: " + maxIntensityOnChart);
        if (scanLine1 > 0 || scanLine2 > 0)
        {
            boolean showedScanLine1 = false;
            boolean showedScanLine2 = false;
            for (int i=0; i<scanValues.length; i++)
            {
                boolean shouldShowScanLine = false;
                
                int scanNumber = (int) scanValues[i];
                int nextScanNumber = -1;
                if (i < scanValues.length)
                    nextScanNumber = (int) scanValues[i+1];

                if (scanLine1 > 0 && !showedScanLine1 && nextScanNumber >= scanLine1)
                {
                    showedScanLine1=true;
                    shouldShowScanLine = true;
                }
                if (scanLine2 > 0 && !showedScanLine2 && scanNumber > scanLine2)
                {
                    showedScanLine2=true;
                    shouldShowScanLine = true;
                }

                if (shouldShowScanLine)
                {
                    for (int j=0; j<intensityValues[i].length; j++)
                        intensityValues[i][j] = -1;
                }

                if ((scanLine1 == 0 || showedScanLine1) && (scanLine2 == 0 || showedScanLine2))
                    break;
            }
        }


        int numScansPadded = maxScan - minScan + 1;
        double[] scanValuesPadded;
        double[][] intensityValuesPadded;
        if (numScansPadded == numScans)
        {
            scanValuesPadded = scanValues;
            intensityValuesPadded = intensityValues;
        }
        else
        {
            _log.debug("Padding!");
            scanValuesPadded = new double[numScansPadded];
            intensityValuesPadded = new double[numScansPadded][numMzBins];

            int unPaddedIndex = 0;
            for (int i=0; i<scanValuesPadded.length; i++)
            {
                int scanValue = minScan + i;
                scanValuesPadded[i] = scanValue;
                
                if (scanValue >= scanValues[unPaddedIndex])
                    unPaddedIndex++;
                for (int j=0; j<numMzBins; j++)
                    intensityValuesPadded[i][j] = intensityValues[unPaddedIndex][j];
            }
        }


        _log.debug("Done loading spectrum in range.");

//        xAxis.setVisible(false);
        setData(scanValuesPadded, mzValues, intensityValuesPadded);

        ((XYPlot) _plot).getDomainAxis().setRange(minScan, maxScan);
        ((XYPlot) _plot).getRangeAxis().setRange(minMz, maxMz);
    }

    /*        Tooltips don't work with XYBlockRenderer
    protected class SpectrumToolTipGenerator implements XYZToolTipGenerator
    {
        protected double[] scanValues;
        protected double[] mzValues;

        public SpectrumToolTipGenerator(double[] scanValues, double[] mzValues)
        {
            this.scanValues = scanValues;
            this.mzValues = mzValues;
        }

        public String generateToolTip(XYDataset data, int series, int item)
        {
System.err.println("asdf1");
            return "scan " + scanValues[series] + ", m/z " + mzValues[item];
        }

        public String generateToolTip(XYZDataset data, int series, int item)
        {
System.err.println("asdf2");            

            return "scan " + scanValues[series] + ", m/z " + mzValues[item];
        }
    }
    */

    public int getMinScan()
    {
        return minScan;
    }

    public void setMinScan(int minScan)
    {
        this.minScan = minScan;
    }

    public int getMaxScan()
    {
        return maxScan;
    }

    public void setMaxScan(int maxScan)
    {
        this.maxScan = maxScan;
    }

    public float getMinMz()
    {
        return minMz;
    }

    public void setMinMz(float minMz)
    {
        this.minMz = minMz;
    }

    public float getMaxMz()
    {
        return maxMz;
    }

    public void setMaxMz(float maxMz)
    {
        this.maxMz = maxMz;
    }

    public int getResolution()
    {
        return resolution;
    }

    public void setResolution(int resolution)
    {
        this.resolution = resolution;
    }

    public MSRun getRun()
    {
        return run;
    }

    public void setRun(MSRun run)
    {
        this.run = run;
    }

    public int getScanLine1()
    {
        return scanLine1;
    }

    public void setScanLine1(int scanLine1)
    {
        this.scanLine1 = scanLine1;
    }

    public int getScanLine2()
    {
        return scanLine2;
    }

    public void setScanLine2(int scanLine2)
    {
        this.scanLine2 = scanLine2;
    }

    public Map<Integer, PanelWithLineChart> getScanLineChartMap()
    {
        return scanLineChartMap;
    }

    public void setScanLineChartMap(Map<Integer, PanelWithLineChart> scanLineChartMap)
    {
        this.scanLineChartMap = scanLineChartMap;
    }
}
