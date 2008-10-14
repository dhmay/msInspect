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
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.toolbox.FloatRange;
import org.fhcrc.cpl.toolbox.gui.chart.*;

import javax.imageio.ImageIO;
import java.util.*;
import java.util.List;
import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import java.awt.*;

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

    protected float lightMz = 0;
    protected float heavyMz = 0;

    protected boolean generateLineCharts = false;

    protected Map<Integer, PanelWithLineChart> scanLineChartMap = null;

    protected boolean generate3DChart = false;
    protected PanelWithRPerspectivePlot contourPlot = null;
    protected int contourPlotWidth = DEFAULT_CONTOUR_PLOT_WIDTH;
    protected int contourPlotHeight = DEFAULT_CONTOUR_PLOT_HEIGHT;
    protected int contourPlotRotationAngle = DEFAULT_CONTOUR_PLOT_ROTATION;
    protected int contourPlotTiltAngle = DEFAULT_CONTOUR_PLOT_TILT;

    public static final int DEFAULT_CONTOUR_PLOT_WIDTH = 1000;
    public static final int DEFAULT_CONTOUR_PLOT_HEIGHT = 1000;
    public static final int DEFAULT_CONTOUR_PLOT_ROTATION = 80;
    public static final int DEFAULT_CONTOUR_PLOT_TILT = 20;

    //this is a bit of a hack.  While we're in here, supply the scan level of a particular scan
    //scan to return the level for
    protected int scanToCheckLevel = -1;
    protected boolean specifiedScanFoundMS1 = false;


    public PanelWithSpectrumChart()
    {
        super();
    }

    protected void init()
    {
        super.init();
        setPalette(PanelWithHeatMap.PALETTE_POSITIVE_WHITE_BLUE_NEGATIVE_BLACK_RED);
        setAxisLabels("scan", "m/z");
        scanLineChartMap = new HashMap<Integer, PanelWithLineChart>();
    }

    public PanelWithSpectrumChart(MSRun run, int minScan, int maxScan, float minMz, float maxMz,
                                  int scanLine1, int scanLine2, float lightMz, float heavyMz)
    {
        init();

        this.run = run;
        this.minScan = minScan;
        this.maxScan = maxScan;
        this.minMz = minMz;
        this.maxMz = maxMz;
        this.scanLine1 = scanLine1;
        this.scanLine2 = scanLine2;
        this.lightMz = lightMz;
        this.heavyMz = heavyMz;
    }

    /**
     *
     */
    public void generateChart()
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
        float minMzForRange = (int) minMz;
        if (minMzForRange > minMz) minMzForRange--;
        float maxMzForRange = (int) maxMz;
        if (maxMzForRange < maxMz) maxMzForRange++;


        FloatRange mzWindowForResample = new FloatRange(minMzForRange, maxMzForRange);
        double[][] intensityValues = new double[numScans][numMzBins];
        _log.debug("Loading spectrum in range....");

        double maxIntensityOnChart = 0;

        for (int scanArrayIndex = 0; scanArrayIndex < numScans; scanArrayIndex++)
        {
            int scanIndex = minScanIndex + scanArrayIndex;
            MSRun.MSScan scan = run.getScan(scanIndex);
            if (scan.getNum() == scanToCheckLevel)
                specifiedScanFoundMS1 = true;
            float[] signalFromResample = Spectrum.Resample(scan.getSpectrum(), mzWindowForResample, resolution);
            int firstIndexToKeep = -1;

            for (int i=0; i<signalFromResample.length; i++)
            {
                float mzValueThisIndex = minMzForRange + (i / (float) resolution);
                if (mzValueThisIndex >= minMz && firstIndexToKeep == -1)
                {
                    firstIndexToKeep = i;
                    break;
                }
            }


            //todo: this is horrible
            double[] signal = new double[numMzBins];
            for (int i=0; i<numMzBins; i++)
                signal[i] = signalFromResample[firstIndexToKeep + i];

            if (generateLineCharts)
            {
                PanelWithLineChart lineChart =
                        new PanelWithPeakChart(mzValues, signal, "Scan " + (int) scanValues[scanArrayIndex]);
                lineChart.setAxisLabels("m/z", "intensity");
                scanLineChartMap.put((int) scanValues[scanArrayIndex], lineChart);
            }

            intensityValues[scanArrayIndex] = signal;

            for (int i=0; i<signal.length; i++)
            {
                if (signal[i] > maxIntensityOnChart)
                    maxIntensityOnChart = signal[i];
            }
        }

        if (generateLineCharts)
        {
            for (PanelWithLineChart lineChart : scanLineChartMap.values())
            {
                ((XYPlot) lineChart.getPlot()).getRangeAxis().setRange(0, maxIntensityOnChart);
            }
        }


        int scanLine1Index = -1;
        int scanLine2Index = -1;


        int numScansPadded = maxScan - minScan + 1;
        double[] scanValuesPadded;
        double[][] intensityValuesPadded;

        if (numScansPadded == numScans)
        {
            scanValuesPadded = scanValues;
            intensityValuesPadded = intensityValues;

            if (scanLine1 > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] < scanLine1)
                        scanLine1Index = i;
                    else break;
                }
            }
            if (scanLine2 > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] <= scanLine2+1)
                        scanLine2Index = i;
                    else break;
                }
            }
        }
        else
        {
            _log.debug("Padding! unpadded: " + numScans + ", padded: " + numScansPadded);

            scanValuesPadded = new double[numScansPadded];
            intensityValuesPadded = new double[numScansPadded][numMzBins];

            int unPaddedIndex = 0;

            for (int i=0; i<scanValuesPadded.length; i++)
            {
                int scanValue = minScan + i;
                scanValuesPadded[i] = scanValue;

                if (unPaddedIndex < scanValues.length-1 && scanValue >= scanValues[unPaddedIndex+1])
                    unPaddedIndex++;

                System.arraycopy(intensityValues[unPaddedIndex], 0, intensityValuesPadded[i], 0,
                        intensityValues[unPaddedIndex].length);
            }

            //add the lines for the scanlines, just outside the specified boundaries
            if (scanLine1 > 0)
            {
                scanLine1Index = 0;
                for (int i=0; i<scanValuesPadded.length; i++)
                {
                    if (scanValuesPadded[i] < scanLine1)
                    {
                        scanLine1Index = i;
                    }
                    else break;
                }
            }
            if (scanLine2 > 0)
            {
                scanLine2Index = 0;
                int nextScanAfterLine2 = 0;
                for (int i=0; i<scanValues.length; i++)
                    if (scanValues[i] > scanLine2)
                    {
                        nextScanAfterLine2 = (int) scanValues[i];
                        break;
                    }
                if (nextScanAfterLine2 == 0)
                    scanLine2Index = 0;
                else
                {
                    for (int i=0; i<scanValuesPadded.length; i++)
                    {
                        if (scanValuesPadded[i] >= nextScanAfterLine2)
                        {
                            scanLine2Index = i;
                            break;
                        }
                    }
                }
            }
        }
        if (scanLine1Index > 0)
            for (int j=0; j<intensityValuesPadded[scanLine1Index].length; j++)
            {
                double origValue = intensityValuesPadded[scanLine1Index][j];
                double newValue = -origValue;
                if (newValue == 0)
                    newValue = -0.000001;
                intensityValuesPadded[scanLine1Index][j] = newValue;
            }
        if (scanLine2Index > 0)
            for (int j=0; j<intensityValuesPadded[scanLine2Index].length; j++)
            {
                double origValue = intensityValuesPadded[scanLine2Index][j];
                double newValue = -origValue;
                if (newValue == 0)
                    newValue = -0.000001;
                intensityValuesPadded[scanLine2Index][j] = newValue;
            }

        //tick marks for specified m/z's
        float intensityForTickMark = -0.001f;
        if (lightMz > 0)
        {
            int closestLightMzIndex = Math.abs(Arrays.binarySearch(mzValues, lightMz));
            intensityValuesPadded[0][closestLightMzIndex] = intensityForTickMark;
            intensityValuesPadded[1][closestLightMzIndex] = intensityForTickMark;

            intensityValuesPadded[intensityValuesPadded.length-1][closestLightMzIndex] = intensityForTickMark;
            intensityValuesPadded[intensityValuesPadded.length-2][closestLightMzIndex] = intensityForTickMark;
        }
        if (heavyMz > 0)
        {
            int closestHeavyMzIndex = Math.abs(Arrays.binarySearch(mzValues, heavyMz));
            intensityValuesPadded[0][closestHeavyMzIndex] = intensityForTickMark;
            intensityValuesPadded[1][closestHeavyMzIndex] = intensityForTickMark;

            intensityValuesPadded[intensityValuesPadded.length-1][closestHeavyMzIndex] = intensityForTickMark;
            intensityValuesPadded[intensityValuesPadded.length-2][closestHeavyMzIndex] = intensityForTickMark;

        }

        if (generate3DChart)
        {
            _log.debug("Generating R contour plot...");

//float[][] intensityValuesPaddedAsFloat = new float[intensityValuesPadded.length][intensityValuesPadded[0].length];
//for (int i=0; i<intensityValuesPadded.length; i++)
//    for (int j=0; j<intensityValuesPadded[0].length; j++)
//       intensityValuesPaddedAsFloat[i][j] = (float) intensityValuesPadded[i][j];
//
//JFrame frame = SurfaceFrame.ShowSurfaceFrame(intensityValuesPaddedAsFloat);
//frame.setVisible(true);

            contourPlot = new  PanelWithRPerspectivePlot();
            contourPlot.setRotationAngle(contourPlotRotationAngle);
            contourPlot.setTiltAngle(contourPlotTiltAngle);
            contourPlot.setChartWidth(contourPlotWidth);
            contourPlot.setChartHeight(contourPlotHeight);
            contourPlot.setUseGradientForColor(true);
            contourPlot.setShowBorders(false);

            if (scanLine1 > 0)
            {
                double[] line1XValues = new double[] {scanLine1,scanLine1};
                double[] line1YValues = new double[] {minMz, maxMz};
                double[] line1ZValues = new double[] {0,0};
                contourPlot.addLine(line1XValues, line1YValues, line1ZValues, "red");
            }
            if (scanLine2 > 0)
            {
                double[] line2XValues = new double[] {scanLine2,scanLine2};
                double[] line2YValues = new double[] {minMz, maxMz};
                double[] line2ZValues = new double[] {0,0};
                contourPlot.addLine(line2XValues, line2YValues, line2ZValues, "red");
            }

            double closestLightMz = mzValues[Math.abs(Arrays.binarySearch(mzValues, lightMz))];
            double closestHeavyMz = mzValues[Math.abs(Arrays.binarySearch(mzValues, heavyMz))];

            double[] tickXValues = new double[] { minScan-1, maxScan+1 };

            double[] lightTickYValues = new double[] { closestLightMz, closestLightMz };
            double[] heavyTickYValues = new double[] { closestHeavyMz, closestHeavyMz };

            double[] tickZValues = new double[] {0,0};

            contourPlot.addLine(tickXValues, lightTickYValues, tickZValues, "red");
            contourPlot.addLine(tickXValues, heavyTickYValues, tickZValues, "red");            


            contourPlot.plot(scanValuesPadded, mzValues, intensityValuesPadded);
            _log.debug("Generated R contour plot.");
        }

        _log.debug("Done loading spectrum in range.");

        setData(scanValuesPadded, mzValues, intensityValuesPadded);

        ((XYPlot) _plot).getDomainAxis().setRange(minScan, maxScan);
        ((XYPlot) _plot).getRangeAxis().setRange(minMz, maxMz);
    }




    /**
     *
     * @param imageWidthEachScan
     * @param imageHeightEachScan
     * @param maxTotalImageHeight a hard boundary on the total image height.  If imageHeightEachScan is too big,
     * given the total number of charts and this arg, it gets knocked down
     * @param outputFile
     * @throws java.io.IOException
     */
    public void savePerScanSpectraImage(int imageWidthEachScan, int imageHeightEachScan, int maxTotalImageHeight,
                                        File outputFile)
            throws IOException
    {
        int numCharts = scanLineChartMap.size();

        int widthPaddingForLabels = 50;

        imageHeightEachScan = Math.min(imageHeightEachScan, maxTotalImageHeight / numCharts);

        List<Integer> allScanNumbers = new ArrayList<Integer>(scanLineChartMap.keySet());
        Collections.sort(allScanNumbers);
        List<PanelWithChart> allCharts = new ArrayList<PanelWithChart>();

        for (int scanNumber : allScanNumbers)
        {
            PanelWithLineChart scanChart = scanLineChartMap.get(scanNumber);
            allCharts.add(scanChart);
            scanChart.setSize(imageWidthEachScan - widthPaddingForLabels, imageHeightEachScan);
        }                 

        BufferedImage perScanChartImage = MultiChartDisplayPanel.createImageForAllCharts(allCharts);

        BufferedImage perScanChartImageWithLabels = new BufferedImage(imageWidthEachScan,
                perScanChartImage.getHeight(), BufferedImage.TYPE_INT_RGB);

        Graphics2D g = perScanChartImageWithLabels.createGraphics();
        g.drawImage(perScanChartImage, widthPaddingForLabels, 0, null);
        g.setPaint(Color.WHITE);
        g.drawRect(0, 0, widthPaddingForLabels, perScanChartImage.getHeight());

        for (int i=0; i<allCharts.size(); i++)
        {
            int scanNumber = allScanNumbers.get(i);

            int chartTop = i * imageHeightEachScan;
            int chartMiddle = chartTop + (imageHeightEachScan / 2);

            if (scanLine1 > 0 && scanLine2 > 0)
            {
                if (scanNumber >= scanLine1 && scanNumber <= scanLine2)
                    g.setPaint(Color.GREEN);
                else
                    g.setPaint(Color.RED);
            }
            else g.setPaint(Color.BLACK);

            g.drawString("" + scanNumber, 5, chartMiddle);
        }
        g.dispose();

        ImageIO.write(perScanChartImageWithLabels,"png",outputFile);
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

    public boolean isSpecifiedScanFoundMS1()
    {
        return specifiedScanFoundMS1;
    }

    public void setSpecifiedScanFoundMS1(boolean specifiedScanFoundMS1)
    {
        this.specifiedScanFoundMS1 = specifiedScanFoundMS1;
    }

    public boolean isGenerate3DChart()
    {
        return generate3DChart;
    }

    public void setGenerate3DChart(boolean generate3DChart)
    {
        this.generate3DChart = generate3DChart;
    }

    public PanelWithRPerspectivePlot getContourPlot()
    {
        return contourPlot;
    }

    public void setContourPlot(PanelWithRPerspectivePlot contourPlot)
    {
        this.contourPlot = contourPlot;
    }

    public int getContourPlotWidth()
    {
        return contourPlotWidth;
    }

    public void setContourPlotWidth(int contourPlotWidth)
    {
        this.contourPlotWidth = contourPlotWidth;
    }

    public int getContourPlotHeight()
    {
        return contourPlotHeight;
    }

    public void setContourPlotHeight(int contourPlotHeight)
    {
        this.contourPlotHeight = contourPlotHeight;
    }

    public int getContourPlotRotationAngle()
    {
        return contourPlotRotationAngle;
    }

    public void setContourPlotRotationAngle(int contourPlotRotationAngle)
    {
        this.contourPlotRotationAngle = contourPlotRotationAngle;
    }

    public int getContourPlotTiltAngle()
    {
        return contourPlotTiltAngle;
    }

    public void setContourPlotTiltAngle(int contourPlotTiltAngle)
    {
        this.contourPlotTiltAngle = contourPlotTiltAngle;
    }

    public boolean isGenerateLineCharts()
    {
        return generateLineCharts;
    }

    public void setGenerateLineCharts(boolean generateLineCharts)
    {
        this.generateLineCharts = generateLineCharts;
    }

    public int getScanToCheckLevel()
    {
        return scanToCheckLevel;
    }

    public void setScanToCheckLevel(int scanToCheckLevel)
    {
        this.scanToCheckLevel = scanToCheckLevel;
    }
}
