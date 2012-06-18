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

package org.fhcrc.cpl.viewer.quant.gui;

import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.viewer.quant.QuantEvent;

import javax.imageio.ImageIO;
import java.util.*;
import java.util.List;
import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import java.awt.*;

/**
 * Constructs several charts with information about a quantitative event
 *
 * TODO: no reason to extend PanelWithHeatMap here -- that's legacy.
 */
public class PanelWithSpectrumChart extends PanelWithHeatMap
{
    static Logger _log = Logger.getLogger(PanelWithSpectrumChart.class);

    //Default mass to assume separates each peak.  Properly this should be the mass of a proton,
    //but for historical reasons it's 1.0Da in Q3, do it's 1.0Da here
    public static final float DEFAULT_PEAK_SEPARATION_MASS = 1.f;

    //Default peak mass tolerance, same as in Q3
    public static final float DEFAULT_PEAK_TOLERANCE_PPM = 25f;
    //number of peaks used by Q3
    public static final int MAX_Q3_PEAKS = 5;

    protected MSRun run;

    protected float ratioOnePeak = QuantEvent.MISSING_RATIO;

    //boundaries
    protected int minScan = 0;
    protected int maxScan = 0;
    protected float minMz = 0;
    protected float maxMz = 0;

    //resolution.  Bumping this up makes things take longer, and also makes things look a bit funny
    protected int resolution = DEFAULT_RESOLUTION;
    public static final int DEFAULT_RESOLUTION = 100;

    //info about the quant event
    protected int lightFirstScanLine = 0;
    protected int lightLastScanLine = 0;
    protected int heavyFirstScanLine = 0;
    protected int heavyLastScanLine = 0;
    protected float lightMz = 0;
    protected float heavyMz = 0;
    protected int charge = 0;

    protected double upperZBoundNearPeak = 0;

    protected int numSafePeaks = -1;

    protected float[] peakMzs = null;


    protected float peakSeparationMass = DEFAULT_PEAK_SEPARATION_MASS;
    protected float peakTolerancePPM = DEFAULT_PEAK_TOLERANCE_PPM;


    protected boolean shouldGenerateLineCharts = false;
    protected Map<Integer, PanelWithLineChart> scanLineChartMap = null;

    protected boolean shouldGenerate3DChart = false;
    protected PanelWithRPerspectivePlot contourPlot = null;

    protected PanelWithLineChart intensitySumChart;
    protected PanelWithLineChart intensitySumPeaksChart;

    //3D parameter defaults
    public static final int DEFAULT_CONTOUR_PLOT_WIDTH = 1000;
    public static final int DEFAULT_CONTOUR_PLOT_HEIGHT = 1000;
    public static final int DEFAULT_CONTOUR_PLOT_ROTATION = 80;
    public static final int DEFAULT_CONTOUR_PLOT_TILT = 20;
    public static final boolean DEFAULT_CONTOUR_PLOT_SHOW_AXES = true;

    //3D plot parameters
    protected int contourPlotWidth = DEFAULT_CONTOUR_PLOT_WIDTH;
    protected int contourPlotHeight = DEFAULT_CONTOUR_PLOT_HEIGHT;
    protected int contourPlotRotationAngle = DEFAULT_CONTOUR_PLOT_ROTATION;
    protected int contourPlotTiltAngle = DEFAULT_CONTOUR_PLOT_TILT;
    protected boolean contourPlotShowAxes = DEFAULT_CONTOUR_PLOT_SHOW_AXES;

    //this is a bit of a hack -- we store and return the scan level of the event scan.
    //As long as we're looking at scans in this run
    protected int idEventScan = -1;
    protected float idEventMz = -1;

    protected boolean specifiedScanFoundMS1 = false;

    //information about other events included with this quant event because they're similar
    protected List<Integer> otherEventScans;
    protected List<Float> otherEventMZs;

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
                                  int lightFirstScanLine, int lightLastScanLine,
                                  int heavyFirstScanLine, int heavyLastScanLine,
                                  float lightMz, float heavyMz, int charge)
    {
        init();

        this.run = run;
        this.minScan = minScan;
        this.maxScan = maxScan;
        this.minMz = minMz;
        this.maxMz = maxMz;
        this.lightFirstScanLine = lightFirstScanLine;
        this.lightLastScanLine = lightLastScanLine;
        this.heavyFirstScanLine = heavyFirstScanLine;
        this.heavyLastScanLine = heavyLastScanLine;
        this.lightMz = lightMz;
        this.heavyMz = heavyMz;
        this.charge = charge;
    }

    /**
     *  Generate all the charts
     */
    public void generateCharts()
    {
        int minScanIndex = Math.abs(run.getIndexForScanNum(minScan));
        int maxScanIndex = Math.abs(run.getIndexForScanNum(maxScan));
        int numScans = maxScanIndex - minScanIndex + 1;
        int numMzBins = (int) ((float) resolution * (maxMz - minMz)) + 1;
        double[] scanValues = new double[numScans];

        double[] mzValues = new double[numMzBins];

        for (int i=0; i<numScans; i++)
            scanValues[i] = run.getScanNumForIndex(minScanIndex + i);
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

        double[] sumIntensitiesInQuantRange = new double[numMzBins];
        //carry just the intensities used in quantitation
        double[] sumPeakIntensitiesInQuantRange = new double[numMzBins];
        if (numSafePeaks < 0) 
            numSafePeaks = Math.min(Math.round((heavyMz - lightMz) * charge), MAX_Q3_PEAKS);
        if (heavyMz < lightMz)
            numSafePeaks = 4;
        if (peakMzs != null)
            numSafePeaks = Math.min(numSafePeaks, peakMzs.length);

        List<Float> peakMzsToQuantitate = new ArrayList<Float>();
        List<Float> nonMonoisotopicPeakMzList = new ArrayList<Float>();

        for (int i=0; i<numSafePeaks; i++)
        {
            float lightPeakMz = lightMz + (peakSeparationMass / charge) * i;
            float heavyPeakMz = heavyMz + (peakSeparationMass / charge) * i;
            if (peakMzs != null)
            {
                lightPeakMz = peakMzs[i];
                heavyPeakMz = heavyMz + (lightPeakMz - lightMz);
            }
            peakMzsToQuantitate.add(lightPeakMz);
            peakMzsToQuantitate.add(heavyPeakMz);
            nonMonoisotopicPeakMzList.add(lightPeakMz);
            if (heavyPeakMz > lightMz)
                nonMonoisotopicPeakMzList.add(heavyPeakMz);
        }
        int heavyMonoisotopicPeakIndex = numSafePeaks;

        float[] nonMonoisotopicPeakMzs = new float[nonMonoisotopicPeakMzList.size()];
        for (int i=0; i<nonMonoisotopicPeakMzs.length; i++)
        {
            nonMonoisotopicPeakMzs[i] = nonMonoisotopicPeakMzList.get(i);
        }
        Collections.sort(peakMzsToQuantitate);

        //dhmay adding 20090723
        //For calculating single-peak ratio based on theoretical highest peak
        int maxIdealPeakIndex = Spectrum.calcMaxIdealPeakIndex(lightMz  - Spectrum.HYDROGEN_ION_MASS) * charge;
        float lightMaxIdealPeakIntensitySum = 0f;
        float heavyMaxIdealPeakIntensitySum = 0f;




        for (int scanArrayIndex = 0; scanArrayIndex < numScans; scanArrayIndex++)
        {
            int scanIndex = minScanIndex + scanArrayIndex;
            MSRun.MSScan scan = run.getScan(scanIndex);
            if (scan.getNum() == idEventScan)
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

            //this is horrible.  arraycopy would be better, but need to convert float to double
            double[] signal = new double[numMzBins];
            for (int i=0; i<numMzBins; i++)
            {
                signal[i] = signalFromResample[firstIndexToKeep + i];

            }

            if (shouldGenerateLineCharts)
            {
                double[] signalPeaks = new double[numMzBins];
                int currentPeakIndex = -1;
                float minMzCurrentPeak = -1;
                float maxMzCurrentPeak = -1;
                for (int i=0; i<numMzBins; i++)
                {
                    double mz = mzValues[i];
                    if (mz > maxMzCurrentPeak && currentPeakIndex < peakMzsToQuantitate.size()-1)
                    {
                        currentPeakIndex++;
                        float currentPeakMz = peakMzsToQuantitate.get(currentPeakIndex);
                        float currentPeakMass = (currentPeakMz - Spectrum.HYDROGEN_ION_MASS) * charge;
                        float deltaMz = (peakTolerancePPM * currentPeakMass / 1000000f) / charge;
                        minMzCurrentPeak = currentPeakMz - deltaMz;
                        maxMzCurrentPeak = currentPeakMz + deltaMz;
                    }
                    if (mz >= minMzCurrentPeak && mz <= maxMzCurrentPeak)
                    {
                        signalPeaks[i] += signal[i];
                    }
                }
                upperZBoundNearPeak = Math.max(upperZBoundNearPeak, BasicStatistics.max(signalPeaks));

                PanelWithLineChart lineChart =
                        new PanelWithPeakChart(mzValues, signalPeaks, "Scan " + (int) scanValues[scanArrayIndex],
                                PanelWithLineChart.defaultShape, Color.BLUE);
                lineChart.addData(mzValues, signal, "all",
                        PanelWithLineChart.defaultShape, Color.GRAY);

                lineChart.setAxisLabels("m/z", "intensity");
                scanLineChartMap.put((int) scanValues[scanArrayIndex], lineChart);
            }

            intensityValues[scanArrayIndex] = signal;

            upperZBound = Math.max(upperZBound, BasicStatistics.max(signal));



            if (scan.getNum() >= lightFirstScanLine && scan.getNum() >= heavyFirstScanLine &&
                    scan.getNum() <= lightLastScanLine && scan.getNum() <= heavyLastScanLine)
            {
                int currentPeakIndex = -1;
                float minMzCurrentPeak = -1;
                float maxMzCurrentPeak = -1;
                for (int i=0; i<sumIntensitiesInQuantRange.length; i++)
                {
                    sumIntensitiesInQuantRange[i] += signal[i];
                    float mz = (float) mzValues[i];

                    if (mz > maxMzCurrentPeak && currentPeakIndex < peakMzsToQuantitate.size()-1)
                    {
                        currentPeakIndex++;
                        float currentPeakMz = peakMzsToQuantitate.get(currentPeakIndex);
                        float currentPeakMass = (currentPeakMz - Spectrum.HYDROGEN_ION_MASS) * charge;
                        float deltaMz = (peakTolerancePPM * currentPeakMass / 1000000f) / charge;
                        minMzCurrentPeak = currentPeakMz - deltaMz;
                        maxMzCurrentPeak = currentPeakMz + deltaMz;
                    }

                    if (mz >= minMzCurrentPeak && mz <= maxMzCurrentPeak)
                    {
                        sumPeakIntensitiesInQuantRange[i] += signal[i];
                        //dhmay adding 20090723 for single-peak quantitation
                        if (currentPeakIndex == maxIdealPeakIndex)
                            lightMaxIdealPeakIntensitySum += signal[i];
                        else if (currentPeakIndex == heavyMonoisotopicPeakIndex + maxIdealPeakIndex)
                            heavyMaxIdealPeakIntensitySum += signal[i];
                    }
                }
            }
        }


        ratioOnePeak = lightMaxIdealPeakIntensitySum / Math.max(heavyMaxIdealPeakIntensitySum, 0.001f);


        double maxIntensityOnSumChart = BasicStatistics.max(sumIntensitiesInQuantRange);
        double maxPeakIntensity = BasicStatistics.max(sumPeakIntensitiesInQuantRange);
        if (maxIntensityOnSumChart < 1)
            maxIntensityOnSumChart = maxPeakIntensity;

        float[] monoisotopicMzs = new float[] { lightMz, heavyMz };

        double intensityForSumChartHeight = maxIntensityOnSumChart * 1.5;
        float[] monoisotopicIntensitiesSumChart =
                new float[] { (float) intensityForSumChartHeight, (float) intensityForSumChartHeight};
        intensitySumChart = new PanelWithPeakChart(mzValues, sumPeakIntensitiesInQuantRange,"Peak Intensities");



        intensitySumChart.addData(mzValues, sumIntensitiesInQuantRange,"Intensity Sum",
                PanelWithLineChart.defaultShape, Color.LIGHT_GRAY);
        if (heavyMz > 0)
            intensitySumChart.addDataFloat(monoisotopicMzs, monoisotopicIntensitiesSumChart,
                "Monoisotopic light and heavy", PanelWithLineChart.defaultShape, Color.GREEN);
        else
            intensitySumChart.addDataFloat(new float[] {lightMz}, new float[] { (float) intensityForSumChartHeight },
                "Monoisotopic light and heavy", PanelWithLineChart.defaultShape, Color.GREEN);
        if (nonMonoisotopicPeakMzs.length > 0)
        {
//for (float mz : nonMonoisotopicPeakMzs) System.err.println("***" + mz);
            float[] nonmonoisotopicIntensitiesSumChart =
                    new float[nonMonoisotopicPeakMzs.length];
            Arrays.fill(nonmonoisotopicIntensitiesSumChart, (float) intensityForSumChartHeight);
//System.err.println("ADDING LINES, peakMzs null ? " + (peakMzs == null));            
            if (peakMzs != null)
            {
                intensitySumChart.addDataFloat(nonMonoisotopicPeakMzs, nonmonoisotopicIntensitiesSumChart,
                        "Nonmonoisotopic light and heavy (exact)", PanelWithLineChart.defaultShape, Color.RED);
            }
            else
            {
                intensitySumChart.addDataFloat(nonMonoisotopicPeakMzs, nonmonoisotopicIntensitiesSumChart,
                        "Nonmonoisotopic light and heavy (approx)", PanelWithLineChart.defaultShape, Color.ORANGE);

            }
        }
        intensitySumChart.getChart().removeLegend();
        intensitySumChart.setAxisLabels("m/z","Intensity (Sum)");
        if (maxPeakIntensity > 1)
            ((XYPlot) intensitySumChart.getPlot()).getRangeAxis().setRange(0, maxPeakIntensity * 1.25);
        intensitySumChart.setSize(getWidth(), getHeight());

        if (shouldGenerateLineCharts)
        {
            float rangeMax = (float) upperZBoundNearPeak * 1.25f;
            for (PanelWithLineChart lineChart : scanLineChartMap.values())
            {
                ((XYPlot) lineChart.getPlot()).getRangeAxis().setRange(0, rangeMax);
                    float[] monoisotopicIntensities =
                            new float[] { (float) rangeMax, rangeMax};
                    lineChart.addDataFloat(monoisotopicMzs, monoisotopicIntensities, "Monoisotopic light and heavy",
                            PanelWithLineChart.defaultShape, Color.GREEN);
                if (nonMonoisotopicPeakMzs.length > 0)
                {
                    float[] nonMonoisotopicIntensities =
                            new float[nonMonoisotopicPeakMzs.length];
                    Arrays.fill(nonMonoisotopicIntensities, rangeMax);
                    lineChart.addDataFloat(nonMonoisotopicPeakMzs, nonMonoisotopicIntensities,
                            "Nonmonoisotopic light and heavy", PanelWithLineChart.defaultShape, Color.ORANGE);
                }
                lineChart.getChart().removeLegend();
            }
        }

        int lightFirstScanLineIndex = -1;
        int lightLastScanLineIndex = -1;
        int heavyFirstScanLineIndex = -1;
        int heavyLastScanLineIndex = -1;

        int numScansPadded = maxScan - minScan + 1;
        double[] scanValuesPadded;
        double[][] intensityValuesPadded;
        setPaintScale(createHeatMapPaintScale());

        if (numScansPadded == numScans)
        {
            scanValuesPadded = scanValues;
            intensityValuesPadded = intensityValues;

            if (lightFirstScanLine > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] < lightFirstScanLine)
                        lightFirstScanLineIndex = i;
                    else break;
                }
            }
            if (lightLastScanLine > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] <= lightLastScanLine +1)
                        lightLastScanLineIndex = i;
                    else break;
                }
            }
            if (heavyFirstScanLine > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] < heavyFirstScanLine)
                        heavyFirstScanLineIndex = i;
                    else break;
                }
            }
            if (heavyLastScanLine > 0)
            {
                for (int i=0; i<scanValues.length; i++)
                {
                    if (scanValues[i] <= heavyLastScanLine +1)
                        heavyLastScanLineIndex = i;
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
            if (lightFirstScanLine > 0)
            {
                lightFirstScanLineIndex = 0;
                for (int i=0; i<scanValuesPadded.length; i++)
                {
                    if (scanValuesPadded[i] < lightFirstScanLine)
                    {
                        lightFirstScanLineIndex = i;
                    }
                    else break;
                }
            }
            if (lightLastScanLine > 0)
            {
                lightLastScanLineIndex = 0;
                int nextScanAfterLine2 = 0;
                for (int i=0; i<scanValues.length; i++)
                    if (scanValues[i] > lightLastScanLine)
                    {
                        nextScanAfterLine2 = (int) scanValues[i];
                        break;
                    }
                if (nextScanAfterLine2 == 0)
                    lightLastScanLineIndex = 0;
                else
                {
                    for (int i=0; i<scanValuesPadded.length; i++)
                    {
                        if (scanValuesPadded[i] >= nextScanAfterLine2)
                        {
                            lightLastScanLineIndex = i;
                            break;
                        }
                    }
                }
            }
            if (heavyFirstScanLine > 0)
            {
                heavyFirstScanLineIndex = 0;
                for (int i=0; i<scanValuesPadded.length; i++)
                {
                    if (scanValuesPadded[i] < heavyFirstScanLine)
                    {
                        heavyFirstScanLineIndex = i;
                    }
                    else break;
                }
            }
            if (heavyLastScanLine > 0)
            {
                heavyLastScanLineIndex = 0;
                int nextScanAfterLine2 = 0;
                for (int i=0; i<scanValues.length; i++)
                    if (scanValues[i] > heavyLastScanLine)
                    {
                        nextScanAfterLine2 = (int) scanValues[i];
                        break;
                    }
                if (nextScanAfterLine2 == 0)
                    heavyLastScanLineIndex = 0;
                else
                {
                    for (int i=0; i<scanValuesPadded.length; i++)
                    {
                        if (scanValuesPadded[i] >= nextScanAfterLine2)
                        {
                            heavyLastScanLineIndex = i;
                            break;
                        }
                    }
                }
            }
        }


        int lightMzIndex = 0;
        for (lightMzIndex=0; lightMzIndex<intensityValuesPadded[0].length; lightMzIndex++)
            if (mzValues[lightMzIndex] > lightMz)
                    break;
        int heavyMzIndex = 0;
        for (heavyMzIndex=0; heavyMzIndex<intensityValuesPadded[0].length; heavyMzIndex++)
            if (mzValues[heavyMzIndex] > heavyMz)
                    break;

        if (lightFirstScanLineIndex > 0)
            shadeArea(intensityValuesPadded, 0, 0, lightFirstScanLineIndex, lightMzIndex);

        if (lightLastScanLineIndex > 0)
            shadeArea(intensityValuesPadded, lightLastScanLineIndex, 0, intensityValuesPadded.length-1, lightMzIndex);

        if (heavyFirstScanLineIndex > 0)
            shadeArea(intensityValuesPadded, 0, heavyMzIndex, heavyFirstScanLineIndex,
                    intensityValuesPadded[0].length-1);

        if (lightLastScanLineIndex > 0)
            shadeArea(intensityValuesPadded, heavyLastScanLineIndex, heavyMzIndex, intensityValuesPadded.length-1,
                    intensityValuesPadded[0].length-1);

        //if light and heavy first scan are same, shade the remaining area on the left.  Same with last scan
        if (lightFirstScanLineIndex > 0 && heavyFirstScanLineIndex == lightFirstScanLineIndex)
            shadeArea(intensityValuesPadded, 0, lightMzIndex+1, lightFirstScanLineIndex, heavyMzIndex-1);
        if (lightLastScanLineIndex > 0 && heavyLastScanLineIndex == lightLastScanLineIndex)
            shadeArea(intensityValuesPadded, lightLastScanLineIndex, lightMzIndex+1,
                    intensityValuesPadded.length-1, heavyMzIndex-1);

        //special value for tick intensity
        float intensityForTickMark = -2f;
        float intensityForIdCross = -2f;

        //cross for ID event
        if (idEventScan > 0 && idEventMz > 0)
        {
            drawHeatmapCrossForEvent(scanValuesPadded, mzValues, intensityValuesPadded,
                    idEventScan, idEventMz, intensityForIdCross);
        }

        if (otherEventScans != null && !otherEventScans.isEmpty())
        {
            for (int i=0; i< otherEventScans.size(); i++)
                drawHeatmapCrossForEvent(scanValuesPadded, mzValues, intensityValuesPadded,
                        otherEventScans.get(i), otherEventMZs.get(i), intensityForIdCross);
        }

        //tick marks for specified m/z's
        if (lightMz > 0)
        {
            placeHeatmapTick(mzValues, intensityValuesPadded, lightMz, intensityForTickMark);
        }
        if (heavyMz > 0)
        {
            placeHeatmapTick(mzValues, intensityValuesPadded, heavyMz, intensityForTickMark);
        }

        if (peakMzs != null)
        {
            for (float mz : nonMonoisotopicPeakMzs)
                placeHeatmapTick(mzValues, intensityValuesPadded, mz, intensityForTickMark);
        }

        if (shouldGenerate3DChart)
            generate3DPlot(scanValuesPadded, mzValues, intensityValuesPadded);

        _log.debug("Done loading spectrum in range.");

        setData(scanValuesPadded, mzValues, intensityValuesPadded);
        getChart().removeLegend();
//try {contourPlot.saveAllImagesToFiles(new File("/home/dhmay/temp/charts"));} catch(IOException e) {}
        ((XYPlot) _plot).getDomainAxis().setRange(minScan, maxScan);
        ((XYPlot) _plot).getRangeAxis().setRange(minMz, maxMz);
    }

    protected void placeHeatmapTick(double[] mzValues, double[][]intensityValuesPadded, float mzVal, float intensityForTickMark)
    {
            int closestLightMzIndex = Math.abs(Arrays.binarySearch(mzValues, mzVal));
            int[] tickXVals = new int[] { 0, 1, 2, 3, 0, 1, 2, 3,
                    intensityValuesPadded.length-1, intensityValuesPadded.length-2,
                    intensityValuesPadded.length-3, intensityValuesPadded.length-4,
                    intensityValuesPadded.length-1, intensityValuesPadded.length-2,
                    intensityValuesPadded.length-3, intensityValuesPadded.length-4};

            int indMaybePlus1 = Math.min(closestLightMzIndex+1, intensityValuesPadded[0].length-1);
            int[] tickYVals = new int[] { closestLightMzIndex, closestLightMzIndex, closestLightMzIndex, closestLightMzIndex,
                    indMaybePlus1, indMaybePlus1, indMaybePlus1, indMaybePlus1,
                    closestLightMzIndex, closestLightMzIndex, closestLightMzIndex, closestLightMzIndex,
                    indMaybePlus1, indMaybePlus1, indMaybePlus1, indMaybePlus1};

            for (int i=0; i<tickXVals.length; i++)
            {
                intensityValuesPadded[tickXVals[i]][tickYVals[i]] = intensityForTickMark;
            }
    }

    /**
     * Create a 3D plot for a quantitative event
     * @param scanValuesPadded
     * @param mzValues
     * @param intensityValuesPadded
     */
    protected void generate3DPlot(double[] scanValuesPadded, double[] mzValues, double[][] intensityValuesPadded)
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
        contourPlot.setName("3D");
        contourPlot.setRotationAngle(contourPlotRotationAngle);
        contourPlot.setTiltAngle(contourPlotTiltAngle);
        contourPlot.setChartWidth(contourPlotWidth);
        contourPlot.setChartHeight(contourPlotHeight);
        contourPlot.setUseGradientForColor(true);
        contourPlot.setShowBorders(false);
        contourPlot.setShowBox(contourPlotShowAxes);

//contourPlot.setTiltAngles(Arrays.asList(new Integer[] { 82, 76, 70, 66, 60, 54, 49, 45, 41, 38, 35, 32, 30, 28, 27, 26, 25, 24, 23 }));
// contourPlot.setRotationAngles(Arrays.asList(new Integer[] {0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
//                    50, 55, 60, 65, 70, 75, 80, 85, 85}));


        double[] twoZeroes = new double[] {0,0};
        if (lightFirstScanLine > 0)
        {
            double[] line1XValues = new double[] {lightFirstScanLine, lightFirstScanLine};
            double[] line1YValues = new double[] {minMz, lightMz};
            contourPlot.addLine(line1XValues, line1YValues, twoZeroes, "red", PanelWithRPerspectivePlot.LINE_STYLE_DASHED);
        }
        if (lightLastScanLine > 0)
        {
            double[] line2XValues = new double[] {lightLastScanLine, lightLastScanLine};
            double[] line2YValues = new double[] {minMz, lightMz};
            contourPlot.addLine(line2XValues, line2YValues, twoZeroes, "red", PanelWithRPerspectivePlot.LINE_STYLE_DASHED);
        }
        if (heavyFirstScanLine > 0)
        {
            double[] line1XValues = new double[] {heavyFirstScanLine, heavyFirstScanLine};
            double[] line1YValues = new double[] {heavyMz, maxMz};
            contourPlot.addLine(line1XValues, line1YValues, twoZeroes, "red", PanelWithRPerspectivePlot.LINE_STYLE_DASHED);
        }
        if (heavyLastScanLine > 0)
        {
            double[] line2XValues = new double[] {heavyLastScanLine, heavyLastScanLine};
            double[] line2YValues = new double[] {heavyMz, maxMz};
            contourPlot.addLine(line2XValues, line2YValues, twoZeroes, "red", PanelWithRPerspectivePlot.LINE_STYLE_DASHED);
        }
        //draw a little X on the MS/MS event
        if (idEventScan > 0 && idEventMz > 0)
        {
            drawContourCrossForEvent(idEventScan, idEventMz, "blue");
        }
        if (otherEventScans != null && !otherEventScans.isEmpty())
        {
            for (int i=0; i< otherEventScans.size(); i++)
                drawContourCrossForEvent(otherEventScans.get(i), otherEventMZs.get(i), "yellow");
        }

        double closestLightMz = mzValues[Math.abs(Arrays.binarySearch(mzValues, lightMz))];
        double closestHeavyMz = mzValues[Math.abs(Arrays.binarySearch(mzValues, heavyMz))];

        //light and heavy monoisotopes
        double[] tickXValues = new double[] { minScan-1, maxScan+1 };

        double[] lightTickYValues = new double[] { closestLightMz, closestLightMz };
        double[] heavyTickYValues = new double[] { closestHeavyMz, closestHeavyMz };

        contourPlot.addLine(tickXValues, lightTickYValues, twoZeroes, "red");
        contourPlot.addLine(tickXValues, heavyTickYValues, twoZeroes, "red");

        contourPlot.plot(scanValuesPadded, mzValues, intensityValuesPadded);
        _log.debug("Generated R contour plot.");
    }

    /**
     * Invert the color of a region of intensity values
     * @param intensityValuesPadded
     * @param xMin
     * @param yMin
     * @param xMax
     * @param yMax
     */
    protected void shadeArea(double[][] intensityValuesPadded, int xMin, int yMin, int xMax, int yMax)
    {
        for (int x=xMin; x<=xMax; x++)
        {
            for (int y=yMin; y<=yMax; y++)
            {
                //Set intensity of background (<1% of max) in shaded region to special value
                if (intensityValuesPadded[x][y] < upperZBound / 100f) intensityValuesPadded[x][y] = -1;
            }
        }
    }

    protected void drawHeatmapCrossForEvent(double[] scanValuesPadded, double[] mzValues, double[][] intensityValuesPadded,
                                     int scan, float mz, float intensityForIdCross)
    {
        //cross for ID event
        if (scan > 0 && mz > 0)
        {
            int closestEventMzIndex = Math.abs(Arrays.binarySearch(mzValues, mz));
            int closestEventScanIndex = Math.abs(Arrays.binarySearch(scanValuesPadded, scan));

            conditionallySetValueAt(intensityValuesPadded, closestEventScanIndex, closestEventMzIndex, intensityForIdCross);
            conditionallySetValueAt(intensityValuesPadded, closestEventScanIndex-1, closestEventMzIndex-1, intensityForIdCross);
            conditionallySetValueAt(intensityValuesPadded, closestEventScanIndex+1, closestEventMzIndex+1, intensityForIdCross);
            conditionallySetValueAt(intensityValuesPadded, closestEventScanIndex-1, closestEventMzIndex+1, intensityForIdCross);
            conditionallySetValueAt(intensityValuesPadded, closestEventScanIndex+1, closestEventMzIndex-1, intensityForIdCross);
            
        }
    }

    protected void conditionallySetValueAt(double[][] intensityValuesPadded, int scanIndex, int mzIndex, float intensity)
    {
        int maxScanIndex = intensityValuesPadded.length-1;
        int maxMzIndex = intensityValuesPadded[0].length-1;

        if (scanIndex >= 0 && scanIndex <= maxScanIndex && mzIndex >= 0 && mzIndex <= maxMzIndex)
            intensityValuesPadded[scanIndex][mzIndex] = intensity;
    }

    protected void drawContourCrossForEvent(int scan, float mz, String color)
    {
        double crossTop = mz+0.1;
        double crossBottom = mz-0.1;
        double crossLeft = scan-2;
        double crossRight = scan+2;

        double[] twoZeroes = new double[] {0,0};
        contourPlot.addLine(new double[] {crossLeft, crossRight},
                new double[] {crossTop, crossBottom}, twoZeroes, color);
        contourPlot.addLine(new double[] {crossLeft, crossRight},
                new double[] {crossBottom, crossTop}, twoZeroes, color);
    }


    /**
     * Create a paint scale that cycles from white to blue in the positive range, defines a shade value at 1.0, and
     * defaults to RED for everything else
     * @return
     */
    protected LookupPaintScale createHeatMapPaintScale()
    {
        LookupPaintScale result = new LookupPaintScale(-1, upperZBound*=1.04, Color.RED);
        addValuesToPaintScale(result, 0, upperZBound, Color.WHITE, Color.BLUE);
        //-1 is special value for background level
        Color shadeColor = new Color(245, 245, 245);
        addValuesToPaintScale(result, -1.1, -0.9, shadeColor, shadeColor);
        
        return result;
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

            if (lightFirstScanLine > 0 && lightLastScanLine > 0)
            {
                if (scanNumber >= lightFirstScanLine && scanNumber <= lightLastScanLine)
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

    public int getLightFirstScanLine()
    {
        return lightFirstScanLine;
    }

    public void setLightFirstScanLine(int lightFirstScanLine)
    {
        this.lightFirstScanLine = lightFirstScanLine;
    }

    public int getLightLastScanLine()
    {
        return lightLastScanLine;
    }

    public void setLightLastScanLine(int lightLastScanLine)
    {
        this.lightLastScanLine = lightLastScanLine;
    }

    public int getHeavyFirstScanLine()
    {
        return heavyFirstScanLine;
    }

    public void setHeavyFirstScanLine(int heavyFirstScanLine)
    {
        this.heavyFirstScanLine = heavyFirstScanLine;
    }

    public int getHeavyLastScanLine()
    {
        return heavyLastScanLine;
    }

    public void setHeavyLastScanLine(int heavyLastScanLine)
    {
        this.heavyLastScanLine = heavyLastScanLine;
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
        return shouldGenerate3DChart;
    }

    public void setGenerate3DChart(boolean generate3DChart)
    {
        this.shouldGenerate3DChart = generate3DChart;
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

    public boolean getContourPlotShowAxes()
    {
        return contourPlotShowAxes;
    }

    public void setContourPlotShowAxes(boolean contourPlotShowAxes)
    {
        this.contourPlotShowAxes = contourPlotShowAxes;
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
        return shouldGenerateLineCharts;
    }

    public void setGenerateLineCharts(boolean generateLineCharts)
    {
        this.shouldGenerateLineCharts = generateLineCharts;
    }

    public int getIdEventScan()
    {
        return idEventScan;
    }

    public void setIdEventScan(int idEventScan)
    {
        this.idEventScan = idEventScan;
    }

    public float getIdEventMz()
    {
        return idEventMz;
    }

    public void setIdEventMz(float idEventMz)
    {
        this.idEventMz = idEventMz;
    }

    public List<Integer> getOtherEventScans()
    {
        return otherEventScans;
    }

    public void setOtherEventScans(List<Integer> otherEventScans)
    {
        this.otherEventScans = otherEventScans;
    }

    public List<Float> getOtherEventMZs()
    {
        return otherEventMZs;
    }

    public void setOtherEventMZs(List<Float> otherEventMZs)
    {
        this.otherEventMZs = otherEventMZs;
    }

    public PanelWithLineChart getIntensitySumChart()
    {
        return intensitySumChart;
    }

    public void setIntensitySumChart(PanelWithLineChart intensitySumChart)
    {
        this.intensitySumChart = intensitySumChart;
    }

    public float getPeakSeparationMass()
    {
        return peakSeparationMass;
    }

    public void setPeakSeparationMass(float peakSeparationMass)
    {
        this.peakSeparationMass = peakSeparationMass;
    }

    public float getPeakTolerancePPM()
    {
        return peakTolerancePPM;
    }

    public void setPeakTolerancePPM(float peakTolerancePPM)
    {
        this.peakTolerancePPM = peakTolerancePPM;
    }

    public int getCharge()
    {
        return charge;
    }

    public void setCharge(int charge)
    {
        this.charge = charge;
    }

    public float getRatioOnePeak()
    {
        return ratioOnePeak;
    }

    public int getNumSafePeaks() {
        return numSafePeaks;
    }

    public void setNumSafePeaks(int numSafePeaks) {
        this.numSafePeaks = numSafePeaks;
    }

    public float[] getPeakMzs() {
        return peakMzs;
    }

    public void setPeakMzs(float[] peakMzs) {
        this.peakMzs = peakMzs;
    }
}
