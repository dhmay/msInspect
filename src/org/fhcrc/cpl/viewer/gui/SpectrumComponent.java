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
package org.fhcrc.cpl.viewer.gui;

import modwt.Filter;
import modwt.Transform;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.util.Haar;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.Smooth2D;
import org.fhcrc.cpl.viewer.feature.FeatureStrategyCentroided;
import org.fhcrc.cpl.viewer.feature.extraction.SmootherCreator;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.datastructure.FloatArray;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentEvent;
import java.beans.PropertyChangeEvent;
import java.util.Arrays;
import java.util.LinkedList;


public class SpectrumComponent
	{
	static Logger _log = Logger.getLogger(SpectrumComponent.class);

	public JPanel contentPanel;
	public JTextField timeTextField;
	public JPanel chartArea;
	public JComboBox displayMode;
        public JComboBox spectrumElutionComboBox;
	public JButton copyButton;
	public JSpinner scanSpinner;
	public JSpinner mzSpinner;

	int _selectedScanNum;
	double _selectedMZ;

	JFreeChart _chart = null;
	ChartPanel _chartPanel = null;
	ListenerHelper helper = new ListenerHelper(this);

	float[][] _spectrum = null;

	String _type = "scan";

        //define the regular background color, so we can switch back after the initial white
        //in the chart panel
        protected static final Color _backgroundColor = new java.awt.Color(238,238,238);

        //dhmay adding to record vertical axis range of spectrum chart, for locking the Y axis
        protected Range _lockedYAxisRange = null;
        protected Range _previousYAxisRange = null;

	public SpectrumComponent()
		{
		try
			{
			Localizer.getSwingEngine(this).render("org/fhcrc/cpl/viewer/gui/SpectrumComponentForm.xml");
			Localizer.localizeComponentTree(contentPanel);  
            contentPanel.doLayout();
			}
		catch (Exception x)
			{
			ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
			throw new RuntimeException(x);
			}

		if (null != displayMode)
			{
			displayMode.addItem("raw");
			displayMode.addItem("background centroided");
			displayMode.addItem("clean centroided");
			displayMode.addItem("resampled");
			displayMode.addItem("subtract background");
			displayMode.addItem("smoothed");
			displayMode.addItem("peaks");
			displayMode.addItem("wavelet peaks");
			displayMode.addItem("wavelet decomposition");
			displayMode.addItem("wavelet decomposition 2");
			displayMode.addItem("wavelet multiresolution");
			displayMode.addItem("haar8");
			displayMode.addItem("haar16");
			helper.addListener(displayMode, "display_actionPerformed");
			}

		if (null != copyButton)
			{
			helper.addListener(copyButton, "copy_actionPerformed");
			}

		if (null != scanSpinner)
			{
			scanSpinner.setModel(new SpinnerRunModel());
			helper.addListener(scanSpinner, "scanSpinner_stateChanged");
			}

		if (null != mzSpinner)
			{
			mzSpinner.setModel(new SpinnerNumberModel(0.0, 0.0, 1000000.0, 1.0/36));
			helper.addListener(mzSpinner, "mzSpinner_stateChanged");
			}

                //dhmay changing from radio group to combo box, 12/19/2005
		if (null != spectrumElutionComboBox)
			{
			spectrumElutionComboBox.addItem("spectrum");
			spectrumElutionComboBox.addItem("elution");
			helper.addListener(spectrumElutionComboBox, "spectrumElutionComboBox_actionPerformed");
                        }

		helper.addListener(contentPanel, "_componentResized");
		helper.addListener(Application.getInstance(), "updateChart_propertyChange", "thresholding");
		//helper.addListener(ApplicationContext, "scan_propertyChange", "MSScan");
		helper.addListener(Application.getInstance(), "selectedPoint_propertyChange", SharedProperties.SELECTED_POINT);
		helper.addListener(Application.getInstance(), "run_propertyChange", SharedProperties.MS_RUN);
		}


	public void spectrumElutionComboBox_actionPerformed(ActionEvent e)
		{
                String previousType = _type;
		String spectrumElutionMode = (String)spectrumElutionComboBox.getSelectedItem();
		if ("elution".equals(spectrumElutionMode))
			setType("elution");
		else
			setType("scan");

                //if changing from elution to spectrum or vice versa, unlock Y axis
                if (!(_type.equals(previousType)))
                  {
                    _lockedYAxisRange = null;
                    _previousYAxisRange = null;
                    updateChart(true);
                  }
                }

	public void scanSpinner_stateChanged(ChangeEvent e)
		{
		Integer Num = (Integer)scanSpinner.getValue();
		MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
		if (null == Num || null == run)
			{
			clearChart();
			return;
			}
		int scanNum = Num.intValue();
		if (_selectedScanNum == scanNum)
			return;

		_selectedScanNum = scanNum;
		_log.debug("_selectedScanNum = " + _selectedScanNum);
        updateChart(true);
		}


	public void mzSpinner_stateChanged(ChangeEvent e)
		{
		Double Num = (Double)mzSpinner.getValue();
		if (null == Num)
			{
			clearChart();
			return;
			}
		if (_selectedMZ == Num.doubleValue())
			return;

		_selectedMZ = Num.doubleValue();
		_log.debug("_selectedMZ = " + _selectedMZ);
        updateChart(true);
		}


	public void setType(String type)
		{
		if (!"scan".equals(type) && !"elution".equals(type))
			throw new IllegalArgumentException(type);
		if (type.equals(_type))
			return;
		_type = type;
		redrawChart();
		}


/*	public void scan_propertyChange(PropertyChangeEvent e)
		{
		if (!"scan".equals(_type))
			 return;
		_selectedScan = (MSRun.MSScan)ApplicationContext.getProperty(SharedProperties.MS_SCAN);
		updateChart();
		} */


	public void run_propertyChange(PropertyChangeEvent e)
		{
		clearChart();
		MSRun run = (MSRun)e.getNewValue();
		((SpinnerRunModel)scanSpinner.getModel()).setRun(run);
		}


	public void selectedPoint_propertyChange(PropertyChangeEvent e)
		{
		Spectrum.Peak p = (Spectrum.Peak)e.getNewValue();

		if (null == p)
			{
			_selectedMZ = -1;
			_selectedScanNum = -1;
			clearChart();
			return;
			}

		if (_selectedMZ == p.mz && _selectedScanNum == p.scan)
			return;

		_selectedMZ = p.mz;
		_selectedScanNum = p.scan;
		_log.debug("_selectedMZ = " + _selectedMZ);
		_log.debug("_selectedScanNum = " + _selectedScanNum);
		updateChart(false);
		}


	public void copy_actionPerformed(ActionEvent event)
		{
		copyChartData();
		}


	public void display_actionPerformed(ActionEvent event)
		{
		updateChart(true);
		}


	public void updateChart_propertyChange(PropertyChangeEvent e)
		{
		updateChart(true);
		}


	public void _componentResized(ComponentEvent event)
		{
		redrawChart();
		}


	public Component getComponent()
		{
		return contentPanel;
		}


	protected void copyChartData()
		{
		if (null == _chart || null == _spectrum)
			return;
		XYPlot xyplot;
		if (_chart.getPlot() instanceof XYPlot)
			{
            xyplot = _chart.getXYPlot();
			}
		else
			{
			// UNDONE:
			return;
			}
		ValueAxis domainAxis = xyplot.getDomainAxis();
		Range range = domainAxis.getRange();
		double min = range.getLowerBound();
		double max = range.getUpperBound();

		float[] x = _spectrum[0];
		float[] y = _spectrum[1];
		int start = Arrays.binarySearch(x, (float)min);
		if (start < 0) start = -(start + 1);
		int end = Arrays.binarySearch(x, (float)max);
		if (end < 0) end = -(end + 1);

		StringBuffer sb = new StringBuffer(20 * (end - start));
		for (int i = start; i < end; i++)
			sb.append(x[i]).append('\t').append(y[i]).append('\n');
		StringSelection sel = new StringSelection(sb.toString());
		Toolkit.getDefaultToolkit().getSystemClipboard().setContents(sel, sel);
		}


	/**
	 * redo chart from scratch (e.g. for a resize)
	 */
	protected void redrawChart()
		{
		// UNDONE: just set preferredSize?
		_chart = null;
		_chartPanel = null;
		updateChart(true);
		}


	/**
	 * update the chart with new data
	 */
	public void updateChart(boolean preserveDomain)
		{
		_log.debug("updateChart: start");
		float[][] spectrumIN = null;
		MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
		MSRun.MSScan selectedScan = getSelectedScan();

		FloatRange rangeSpectrum = null;    // restricted range of displayed data (for autoZoom)
		Pair<Spectrum.Peak, Spectrum.Peak> chartRange = null;             // tell other components what we are displaying

		Range domainRange = null;
		if (preserveDomain && _chart != null)
			domainRange = _chart.getXYPlot().getDomainAxis().getRange();

		if (null == run || null == selectedScan)
		{
			clearChart();
			//if no chart to display, show white background
            chartArea.setBackground(java.awt.Color.WHITE);
			return;
		}

        chartArea.setBackground(_backgroundColor);

		if ("scan".equals(_type))
			{
			scanSpinner.setValue(new Integer(selectedScan.getNum()));
			double t = selectedScan.getDoubleRetentionTime();
			timeTextField.setText("" + Math.round(t*10)/10.0);
			scanSpinner.setEnabled(true);
			//mzTextField.setText("");
			//mzTextField.setEnabled(false);
			mzSpinner.setValue(new Double(_selectedMZ));
			mzSpinner.setEnabled(false);

/*			ArrayList l = new ArrayList();
			l.add(selectedScan.getSpectrum());
            spectrumIN = Spectrum.CombineRawSpectra(
                    run.getTemplateSpectrum(),
                    l,
                    new FloatRange(selectedScan.getLowMz(),selectedScan.getHighMz())
                    ); */

			spectrumIN = selectedScan.getSpectrum();

			if (autoZoom())
				{
				// we don't want to filter the analyzed data, in case
				// that affects any of the algorithms
				// we only filter the result Datasets
				// UNDONE: choose range based on size of MSDetailPanel

				// range of 27+ is enough to see 4 ICAT labelled features (9*(4-1))
				rangeSpectrum = new FloatRange((float)_selectedMZ - 10F, (float)_selectedMZ +30F);
				}

			if (null == rangeSpectrum)
				{
				chartRange = new Pair<Spectrum.Peak, Spectrum.Peak>(
				        new Spectrum.Peak(_selectedScanNum, spectrumIN[0][0], 0.0F),
				        new Spectrum.Peak(_selectedScanNum, spectrumIN[0][spectrumIN[0].length-1]));
				}
			else
				{
				chartRange = new Pair<Spectrum.Peak, Spectrum.Peak>(
				        new Spectrum.Peak(_selectedScanNum, rangeSpectrum.min, 0.0F),
				        new Spectrum.Peak(_selectedScanNum, rangeSpectrum.max, 0.0F));
				}
			}
		else // elution
			{
			if (_selectedMZ == -1)
				{
				clearChart();
				return;
				}

			int RESAMPLE = 36;
			float mz = (float)Math.round(_selectedMZ * RESAMPLE) / RESAMPLE;
			_log.debug("updateChart mz=" + mz);

			scanSpinner.setValue(new Integer(selectedScan.getNum()));
			scanSpinner.setEnabled(false);
			double t = selectedScan.getDoubleRetentionTime();
			timeTextField.setText("" + Math.round(t*10)/10.0);
			//mzTextField.setText("" + mz);
			//mzTextField.setEnabled(true);
			mzSpinner.setValue(new Double(_selectedMZ));
			mzSpinner.setEnabled(true);

			FloatArray arrayIntensity = new FloatArray();
			FloatArray arrayScan = new FloatArray();
			int scanIndex = run.getIndexForScanNum(selectedScan.getNum());
			int scanStart = Math.max(0,scanIndex-64);
			int scanEnd = Math.min(run.getScanCount() - 1, scanStart + 128);
			FloatRange r = new FloatRange(mz, mz);
			for (int s = scanStart ; s < scanEnd ; s++)
				{
				MSRun.MSScan scan = run.getScan(s);
				float[][] spectrumRaw = scan.getSpectrum();
				float[] resample = Spectrum.Resample(spectrumRaw, r, RESAMPLE);
				float f = resample[0];
				arrayIntensity.add(f);
				//arrayScan.add((float)scan.getDoubleRetentionTime());
//				arrayScan.add((float)s);
				arrayScan.add((float)scan.getNum());
				}
			spectrumIN = new float[][] {arrayScan.toArray(null), arrayIntensity.toArray(null)};

			chartRange = new Pair<Spectrum.Peak, Spectrum.Peak>(
			        new Spectrum.Peak(run.getScan(scanStart).getNum(), mz, 0.0F),
			        new Spectrum.Peak(run.getScan(scanEnd).getNum(), mz, 0.0F));
			}

		String mode = (String)displayMode.getSelectedItem();

		java.util.List<XYSeriesCollection> listMultipleCharts = null;
		XYSeriesCollection series = new XYSeriesCollection();
		series.setIntervalWidth(0.0);


		//
		// Process source spectrum
		//

		float[][] spectrum = new float[][] {(float[])spectrumIN[0].clone(), (float[])spectrumIN[1].clone()};
		ProcessSpectrum:
			{
			if ("raw".equals(mode))
				{
				break ProcessSpectrum;
				}

			if ("background centroided".equals(mode))
				{
				if (run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
					spectrum = FeatureStrategyCentroided.backgroundSpectrum(spectrum);
				break ProcessSpectrum;
				}

			if ("clean centroided".equals(mode))
				{
                if (run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
					spectrum = FeatureStrategyCentroided.cleanSpectrum(spectrum);
				break ProcessSpectrum;
				}


			// all subsequent processing expects resampling
			// don't resample elution profile
			if ("scan".equals(_type))
				{
				int len = spectrum[0].length;
				spectrum = Spectrum.ResampleSpectrum(spectrumIN,
						new FloatRange((float)Math.floor(spectrum[0][0]), (float)Math.ceil(spectrum[0][len-1])),
						36, false);
				}
			else
				{
				spectrum[1] = Spectrum.MedianSmooth(spectrum[1]);
				Spectrum.SmoothALittle(spectrum[1]);
				}

			if ("resampled".equals(mode))
				break ProcessSpectrum;

/*			if ("compressed".equals(mode)) (data encoding test)
				{
				float[] s = spectrum[1];
				short c[] = new short[s.length];
				double m = 1.0;
				for (int i = 0 ; i<c.length ; i++)
					m = Math.max(m, s[i]);
				double sqrt = Math.sqrt(m+1);
				m = Math.round(Math.log(sqrt))+1;
				double f = Math.floor(0x7fff / Math.log(sqrt));
				m = 0;
				for (int i = 0 ; i<c.length ; i++)
					{
					c[i] = (short)(Math.round(Math.log((s[i]+1)/sqrt) * f));
					m = Math.max(m, Math.abs(c[i]));
					}
				System.err.println("MAX " + m);
				for (int i = 0 ; i<c.length ; i++)
					s[i] = (float)((Math.exp((double)c[i]/f) * sqrt) - 1);
				break ProcessSpectrum;
				} */

			// remove background for further processing
			if (run.getHeaderInfo().getDataProcessing().getCentroided() != 1)
				{
				int window = "scan".equals(_type) ? 72 : 15;
				float x[] = spectrum[1];
				float bg[] = Spectrum.MinimaWindow(x, spectrum[0].length, window, null);
				for (int i=0 ; i<bg.length ; i++)
					bg[i] = Math.max(0, x[i] - bg[i]);
				spectrum = new float[][] {spectrum[0], bg};
				}


			if ("subtract background".equals(mode))
				break ProcessSpectrum;


			if (mode.startsWith("threshold"))
				{
				spectrum[1] = Spectrum.WaveletD3(spectrum[1], null);
				break ProcessSpectrum;
				}

			if ("peaks".equals(mode) || "smoothed".equals(mode))
				{
				if ("scan".equals(_type))
					{
					double s = Smooth2D.smoothYfactor;
					spectrum[1] = Spectrum.FFTsmooth(spectrum[1], s, false);
					}
				else
					{
					//FeatureStrategyPeakClusters.smoothTHRESHOLD sm = new FeatureStrategyPeakClusters.smoothTHRESHOLD();
					spectrum[1] = SmootherCreator._thresholdElution(spectrum[1]);
					}
				break ProcessSpectrum;
				}
			} // ProcessSpectrum:


		//
		// This is the source spectrum
		//

		String name = mode.indexOf("peaks") != -1 ? "-spectrum" : "|spectrum";
		series.addSeries(new SpectrumXYSeries(name, spectrum, rangeSpectrum));
		_spectrum = spectrum;


		//
		// add additional series/charts
		//

		// show min/med for reference
        if (mode.equals("resampled") || mode.equals("subtract background"))
	        {
	        float[] T;
	        if (mode.equals("resampled"))
		        {
		        T = Spectrum.MinimaWindow(spectrum[1], spectrum[0].length, 72, null);
	            series.addSeries(new SpectrumXYSeries("-min", new float[][]{spectrum[0], T}, rangeSpectrum));
		        }
	        T = Spectrum.MedianWindow(spectrum[1], spectrum[0].length, 72, false);
	        series.addSeries(new SpectrumXYSeries("-med", new float[][]{spectrum[0], T}, rangeSpectrum));
	        }

		if (mode.startsWith("wavelet decomposition") || mode.startsWith("wavelet multiresolution"))
			{
			Filter f = new Filter("haar");
			float[] signal = spectrum[1]; // Spectrum.PadToDouble(spectrum[1], 32);
			int N = spectrum[0].length; // signal.length;
			int levels = 5;

			//if ("wavelet decomposition 2".equals(mode)) Spectrum.SmoothALittle(signal);

			float[][] modwt1 = Transform.decompose(signal, N, levels, f, "modwt", "periodic", null);
			float[][] t = modwt1;

			if ("wavelet multiresolution".equals(mode))
				{
				t = Transform.multiresolution(t, N, levels, f, "modwt", "periodic", null);
				}
			else if ("wavelet decomposition 2".equals(mode))
				{
				float[][] mra = Transform.multiresolution(t, N, levels, f, "modwt", "periodic", null);
				float[][] modwt2 = Transform.decompose(t[2], N, levels, f, "modwt", "periodic", null);
				// show original, d1, and d2
				float[] a = modwt1[2];
				float[] b = modwt2[2];
				float[] m = mra[2];
				Spectrum.Rotate(a, -3);
				Spectrum.Rotate(b, -7);
				t = new float[][] {a, /*b,*/ m}; // b is a mirror image of m
				}

			// copy array into data series
			listMultipleCharts = new java.util.LinkedList<XYSeriesCollection>();
			XYSeriesCollection ds = new XYSeriesCollection(new SpectrumXYSeries("spectrum", spectrum, rangeSpectrum));
			ds.setIntervalWidth(0.0);
			listMultipleCharts.add(ds);
			for (int i=0 ; i<t.length ; i++)
				{
				//float[] s = Spectrum.UnpadToFloat(t[i], 32, null);
				float[][] l = new float[][] {spectrum[0], t[i]};
				ds = new XYSeriesCollection(new SpectrumXYSeries("level " + (i+1), l, rangeSpectrum));
				ds.setIntervalWidth(0.0);
				listMultipleCharts.add(ds);
				}
			}
		else if (mode.startsWith("haar"))
			{
			int l = Integer.parseInt(mode.substring(4));
			listMultipleCharts = new LinkedList<XYSeriesCollection>();

			float[] t1 = Haar.transform(spectrum[1], l);
			float[] t2 = Haar.transform(t1, l);
			XYSeriesCollection ds = new XYSeriesCollection(new SpectrumXYSeries("spectrum", spectrum, rangeSpectrum));
			ds.setIntervalWidth(0.0);
			listMultipleCharts.add(ds);
			XYSeriesCollection ds1 = new XYSeriesCollection(new SpectrumXYSeries(mode, new float[][] {spectrum[0], t1}, rangeSpectrum));
			ds1.setIntervalWidth(0.0);
			listMultipleCharts.add(ds1);
			XYSeriesCollection ds2 = new XYSeriesCollection(new SpectrumXYSeries(mode, new float[][] {spectrum[0], t2}, rangeSpectrum));
			ds2.setIntervalWidth(0.0);
			listMultipleCharts.add(ds2);
			}
		else if ("peaks".equals(mode) || "threshold peaks".equals(mode))
			{
			double noise = 0.1; //"peaks".equals(mode) ? 2.0 : 1.0; // Spectrum.Noise(spectrum[1], 0, spectrum[1].length);
			int[] peakList = Spectrum.PickPeakIndexes(spectrum[1], noise);

			float[][] peaks = new float[2][peakList.length];
			for (int i = 0; i < peakList.length; i++)
				{
				int p = peakList[i];
				if (p >= spectrum[0].length) continue;
				peaks[0][i] = spectrum[0][p];
				peaks[1][i] = spectrum[1][p];
				}
			series.addSeries(new SpectrumXYSeries("|peaks", peaks, rangeSpectrum));
			_spectrum = peaks;
			}
		else if ("wavelet peaks".equals(mode))
			{
			//Spectrum.Peak[] peaks = Spectrum.WaveletPeaks(spectrum);
			Spectrum.Peak[] peaks = Spectrum.WaveletPeaksD3(spectrum);
			float[][] peakSpectrum = new float[2][peaks.length];
			for (int p=0 ; p<peaks.length ; p++)
				{
				peakSpectrum[0][p] = peaks[p].mz;
				peakSpectrum[1][p] = peaks[p].intensity;
				}

			// show d3 and thresholded spectrum
			int levels = 3, N = spectrum[0].length;
			Filter f = new Filter("haar");
			float[][] modwt = Transform.decompose(spectrum[1], N, levels, f, "modwt", "periodic", null);
			float[][] mra = Transform.multiresolution(modwt, N, levels, f, "modwt", "periodic", null);
			float[] thresholded = new float[N];
			for (int i=0 ; i<N ; i++)
				thresholded[i] = mra[2][i] + mra[3][i];

			series.removeAllSeries();
			//series.addSeries(new SpectrumXYSeries("-spectrum", new float[][] {spectrum[0],thresholded}, rangeSpectrum));
			series.addSeries(new SpectrumXYSeries("-mra", new float[][] {spectrum[0], mra[2]}, rangeSpectrum));
			series.addSeries(new SpectrumXYSeries("|peaks", peakSpectrum, rangeSpectrum));
			_spectrum = peakSpectrum;
		}

		//
		// now update or create chart
		//

		Color[] colors = series.getSeriesCount() == 3 ?
			        new Color[] {Color.RED, Color.BLUE, Color.BLUE} :
		        series.getSeriesCount() == 2 ?
		            new Color[] {Color.BLUE, Color.RED} :
		            new Color[] {Color.RED};

		// NOTE: the more often we call setDatasets instead of creating new chart the better
		// CONSIDER: if we don't save chart, at least preserve zoom settings
		if (false && _chart != null && _chartPanel != null && !(_chart.getPlot() instanceof CombinedDomainXYPlot) && listMultipleCharts == null)
		{
			SpectrumChartFactory.setColors(_chartPanel, colors);
			_chart.getXYPlot().setDataset(series);
			_chartPanel.updateUI();
		}
		else
		{
			if (listMultipleCharts == null)
				{
				_log.debug("updateChart: series=" + series.getSeriesCount() + " length(0)=" + series.getItemCount(0));
				_chartPanel = SpectrumChartFactory.CreateChartPanel(series, colors);
				}
			else
				{
				_log.debug("updateChart: charts=" + listMultipleCharts.size());
				_chartPanel = SpectrumChartFactory.CreateChartPanel(listMultipleCharts, colors);
				}
			_chart = _chartPanel.getChart();

			// there seem to be mystery margins so give a little extra space
			Dimension size = chartArea.getSize();
			_chartPanel.setPreferredSize(new Dimension(size.width - 10, size.height - 10));
			chartArea.removeAll();
			chartArea.add(_chartPanel);
			chartArea.doLayout();
			}
		if (null != domainRange)
			_chart.getXYPlot().getDomainAxis().setRange(domainRange);

                //dhmay: if the user has locked the Y axis, then try to use the stored Y axis range
                if (isYAxisLocked())
                {
                  //a locked axis value won't be available if the axis was _just_ locked
                  if (null == _lockedYAxisRange)
                  {
                    //if we've already displayed a chart, use that chart's axis range.  Otherwise,
                    //don't force, and record this chart's axis range
                    if (_previousYAxisRange != null)
                      _lockedYAxisRange = _previousYAxisRange;
                    else
                      _lockedYAxisRange = _chart.getXYPlot().getRangeAxis().getRange();
                  }
                  _chart.getXYPlot().getRangeAxis().setRange(_lockedYAxisRange);
                }
                else
                {
                  //if the Y axis isn't locked, then dump the stored range so it isn't used later
                  _lockedYAxisRange = null;
                }

                //always store the previous Y axis range
                if (_chart.getXYPlot() != null
                    && _chart.getXYPlot().getRangeAxis() != null
                    && _chart.getXYPlot().getRangeAxis().getRange() != null) {
                _previousYAxisRange = _chart.getXYPlot().getRangeAxis().getRange();
                }

		ApplicationContext.setProperty(SharedProperties.CHART_RANGE, chartRange);
		}


	/* UNDONE move to Spectrum.java
	private static int removeZeros(float[][] spectrum)
		{
		int dst = 0;
		float[] s = spectrum[1];
		float[] mz = spectrum[0];
		for (int src = 0; src < mz.length; src++)
			{
			if (s[src] != 0.0F)
				{
				mz[dst] = mz[src];
				s[dst] = s[src];
				dst++;
				}
			}
		if (dst < mz.length) Arrays.fill(mz, dst, mz.length-1, 0.0F);
		if (dst < s.length) Arrays.fill(s, dst, s.length-1, 0.0F);
		return dst;
		}*/


	protected void clearChart()
		{
		_log.debug("clearChart");
		scanSpinner.setValue(null);
		scanSpinner.setEnabled(false);
		timeTextField.setText("");
		chartArea.removeAll();
		_chart = null;
		}


	public static class SpinnerRunModel extends SpinnerNumberModel
		{
        MSRun _run = null;
		int _scanIndex = -1;

		SpinnerRunModel()
			{
			}

		void setRun(MSRun run)
			{
			if (_run == run)
				return;
			_run = run;
			fireStateChanged();
			}

		public Object getNextValue()
			{
			if (null == _run)
				return null;
			if (_scanIndex >= _run.getScanCount())
				return null;
			_scanIndex++;
			fireStateChanged();
			MSRun.MSScan scan = _run.getScan(_scanIndex);
			return new Integer(scan.getNum());
			}

		public Object getPreviousValue()
			{
			if (null == _run)
				return new Integer(0);
			if (_scanIndex <= 0)
				return null;
			_scanIndex--;
			fireStateChanged();
			MSRun.MSScan scan = _run.getScan(_scanIndex);
			return new Integer(scan.getNum());
			}

		public Object getValue()
			{
			if (null == _run || _scanIndex < 0 ||_scanIndex >= _run.getScanCount())
				return new Integer(0);
			return new Integer(_run.getScan(_scanIndex).getNum());
			}

		public int getIndex()
			{
			return _scanIndex;
			}

		public void setValue(Object value)
			{
			if (null == _run)
				return;
			if (!(value instanceof Integer))
				return;
			int v = ((Integer)value).intValue();

			Object current = getValue();
			if (null != current && v == ((Integer)current).intValue())
				return;

			int i = _run.getIndexForScanNum(v);
			i = i < 0 ? -(i+1) : i;
			_scanIndex = i;
			fireStateChanged();

			return;
			/*for (int i=0, count=_run.getScanCount() ; i<count ; i++)
				{
				int num = _run.getScan(i).getNum();
				if (num < v) continue;
				if (num > v) break;
                _scanIndex = i;
				fireStateChanged();
				return;
				}*/
			}
		}


	/*
	public void setPoint(Point point)
		{
		if (null == point)
			return;

		if (_chart != null && _chartPanel != null)
			{
			float min = (float)(point.getY() - 10);
			float max = min + 20;
			XYPlot xyplot = _chart.getXYPlot();
			ValueAxis domainAxis = xyplot.getDomainAxis();
			domainAxis.setRange(min, max);

			float[][] spectrum = _selectedScan.getSpectrum();
			double maxIntensity = 0;
			int i = 0;
			for (i = 0; i < spectrum[0].length && spectrum[0][i] < min; i++)
				;

			for (; i < spectrum[0].length && spectrum[0][i] < max; i++)
				if (spectrum[1][i] > maxIntensity)
					maxIntensity = spectrum[1][i];

			ValueAxis rangeAxis = xyplot.getRangeAxis();
			if (null != rangeAxis)
				rangeAxis.setRange(0, maxIntensity * 1.1);
			return;
			}
		}*/


	protected MSRun.MSScan getSelectedScan()
		{
		MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
		if (null == run)
			return null;
		int i = run.getIndexForScanNum(_selectedScanNum);

		if (i < 0)
			{
			// Try getting the previous scan (helps if we're looking at ms2 features)
			if (i == -1)
				return null;
			i = -i - 2;
			_selectedScanNum = run.getScan(i).getNum();
			}
		return run.getScan(i);
		}


	protected boolean autoZoom()
		{
		Boolean B = (Boolean)ApplicationContext.getProperty(SharedProperties.AUTO_ZOOM);
		return null == B ? false : B.booleanValue();
		}

        //dhmay adding 12/19/2005
	protected boolean isYAxisLocked()
		{
		Boolean B = (Boolean)ApplicationContext.getProperty(SharedProperties.LOCK_Y_AXIS);
		return null == B ? false : B.booleanValue();
		}
       
	}
