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
package org.fhcrc.cpl.viewer.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.ExtractEdgeFeatures;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.io.InputStream;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Properties;

/**
 * User: mbellew
 * Date: Nov 2, 2004
 *
 * This is meant to represent the peaks used in FeatureStrategyCombined
 */
public class FeatureStrategyPeaks extends FeatureStrategyUsingWindow // extends FeatureExtractor
	{
	static Logger _log = Logger.getLogger(FeatureStrategyPeaks.class);

	static int _WindowMargin = 64;
	static int _FeatureScanWindowStart;
	static int _FeatureScanWindowWidth;
	static float _FeatureMzWindowStart;
	static float _FeatureMzWindowHeight;
	static float _AverageWindowWidth;

	static CPUTimer timerAnalyze = new CPUTimer("FeatureStrategyCombined.analyze");
	static CPUTimer timerResample = new CPUTimer("FeatureStrategyCombined.get and resample");
	static CPUTimer timerBackground = new CPUTimer("FeatureStrategyCombined.background");
	static CPUTimer timerEdgeFeatures = new CPUTimer("FeatureStrategyCombined.edges");
	static CPUTimer timerExtractPeaks = new CPUTimer("FeatureStrategyCombined.peaks");
	static CPUTimer timerExtractPeptides = new CPUTimer("FeatureStrategyCombined.peptides");

	static
	{
	initProps();
	}

	private static boolean bPropsInitialized = false;


	int _startNum = 0;
	int _endNum = 0;
	Scan[] _scans = null;


	private static void initProps()
		{
		assert false == (bPropsInitialized = false);
		if (bPropsInitialized)
			return;

		try
			{
			InputStream is = ExtractEdgeFeatures.class.getResourceAsStream("features.properties");
			Properties props = new Properties();
			props.load(is);
			_FeatureScanWindowStart = Integer.parseInt(((String)props.get("feature.ScanWindowStart")).trim());
			_FeatureScanWindowWidth = Integer.parseInt(((String)props.get("feature.ScanWindowWidth")).trim());
			_FeatureMzWindowStart = Float.parseFloat(((String)props.get("feature.MzWindowStart")).trim());
			_FeatureMzWindowHeight = Float.parseFloat(((String)props.get("feature.MzWindowHeight")).trim());
			_AverageWindowWidth = Integer.parseInt(((String)props.get("feature.AverageWindowWidth")).trim());
			}
		catch (java.io.IOException x)
			{
			x.printStackTrace();
			}
		}


	public FeatureStrategyPeaks(MSRun run, int scanIndex, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scanIndex, count, maxCharge, range, sn);

		int scanMax = Math.min(scanIndex + count, run.getScanCount());
		count = scanMax - scanIndex;
		_scans = getScans(run, scanIndex, count);

		_startNum = run.getScan(scanIndex).getNum();
		_endNum = run.getScan(scanMax-1).getNum();
		}


	public int getType()
		{
		return TYPE_2D; // lie for detail pane
		}


	public Feature[] _analyze() throws InterruptedException
		{
		//return analyzeScanAtATime(_scans);
		return analyzeWindow(_scans, 256, 32);
		}


	/**
	 * THIS IS THE MAIN FEATURE FINDING ROUTINE
	 */
	protected Collection analyze1D(Scan scan)
		{
		float[][] spectrum = Spectrum.ResampleSpectrum(scan.getSpectrum(), _mzRange,
                SpectrumResampler.getResampleFrequency(), false);

		// HACK try averaging
		int scanIndex=0;
		for (; scanIndex < _scans.length; scanIndex++)
			if (scan == _scans[scanIndex])
				break;
		int c=1;
		if (scanIndex-1 >=0)
			{
			c++;
			float[] t = Spectrum.Resample(_scans[scanIndex-1].getSpectrum(), _mzRange,
                    SpectrumResampler.getResampleFrequency());
			for (int i = 0; i < t.length; i++)
				spectrum[1][i] += t[i];
			}
		if (scanIndex+1 < _scans.length)
			{
			c++;
			float[] t = Spectrum.Resample(_scans[scanIndex+1].getSpectrum(), _mzRange,
                    SpectrumResampler.getResampleFrequency());
			for (int i = 0; i < t.length; i++)
				spectrum[1][i] += t[i];
			}
		for (int i = 0; i < spectrum[1].length; i++)
			spectrum[1][i] /= c;
		// END HACK

		float[] signal = spectrum[1];
		int height = signal.length;

		// separate background and signal components
		// we're pretty conservative about what we call background,
		// we can use background and/or median later to filter peaks
		// further.

		float[] background = Spectrum.RemoveBackground(new float[][] {signal})[0];
		float[] median = Spectrum.MedianWindow(signal, height, 72, false);

		assert timerExtractPeaks.start();
		// .analyze() is destructive, need to copy
		float[] spectraT = (float[])signal.clone();
        float[] d3 = Spectrum.WaveletD3(spectraT, null);
		int[] indexes = Spectrum.PickPeakIndexes(d3, 0F);

		ArrayList peaks = new ArrayList();
		for (int i = 0; i < indexes.length; i++)
			{
			int index = indexes[i];
			float mz = spectrum[0][index];
			float in = spectrum[1][index];
			if (in < median[index] * 2)
				continue;
			Spectrum.Peak p = new Spectrum.Peak(scan.getNum(), mz, in);
			p.setMedian(median[index]);
			p.setBackground(background[index]);
			peaks.add(p);
			}
		Spectrum.Peak[] arr = (Spectrum.Peak[])peaks.toArray(new Spectrum.Peak[0]);
		return PeaksAsFeatures(arr);
		}


	protected Collection analyze2D(Scan[] scans)
		{
		_log.debug("analyze2D " + scans[0].getNum() + "-" + scans[scans.length - 1].getNum());
		assert timerAnalyze.start();

		//
		// Convert data into 2D matrix
		// we will do all processing on this data until the end and
		// then process back to "scan" space
		//
		assert timerResample.start();
		float[][] spectra = new float[scans.length][];
		for (int i = 0; i < scans.length; i++)
			spectra[i] = Spectrum.Resample(scans[i].getSpectrum(), _mzRange, SpectrumResampler.getResampleFrequency());
		int width = spectra.length;
		int height = spectra[0].length;
		// a little smoothing across elution, median eliminate lock spray!
		float[] row = null, t = null;
		for (int i=0 ; i<height ; i++)
			{
			row = Spectrum.getRow(spectra, i, row);
			t = Spectrum.MedianSmooth(row, row.length, t);
			Spectrum.SmoothALittle(t);
			Spectrum.setRow(spectra, i, t);
			}
		assert timerResample.stop();

		_log.debug("analyze2D datasize = " + (width * height * 4));

		// separate background and signal components
		// we're pretty conservative about what we call background,
		// we can use background and/or median later to filter peaks
		// further.

		timerBackground.start();

		float[][] original = new float[width][];
		for (int i=0 ; i<width ; i++)
			original[i] = (float[])spectra[i].clone();

		float[][] background = Spectrum.RemoveBackground(spectra);
		float[][] median = new float[spectra.length][];
		for (int i=0 ; i<width ; i++)
			median[i] = Spectrum.MedianWindow(spectra[i], height, 72, false);
		timerBackground.stop();

		ArrayList peaks = new ArrayList();

		float[][] wavelets = new float[width][];
		Pair tmp = new Pair(null, null);
		for (int s = 0; s < spectra.length; s++)
			{
			float[] spectraT = (float[])spectra[s].clone();
	        float[] d3 = Spectrum.WaveletD3(spectraT, tmp);
			wavelets[s] = d3;

			int[] indexes = Spectrum.PickPeakIndexes(d3, 0F);

			for (int p = 0; p < indexes.length; p++)
				{
				int index = indexes[p];
				float mz = this._mzRange.min + index/ SpectrumResampler.getResampleFrequency();
				float in = spectra[s][index];
				if (in < median[s][index] * 2)
					continue;
				Spectrum.Peak peak = new Spectrum.Peak(_scans[s].getNum(), mz, in);
				peak.setMedian(median[s][index]);
				peak.setBackground(background[s][index]);

				if (peak.intensity < (0.5*peak.background + 2*peak.getMedian() + 3))
					continue;

				peaks.add(peak);
				}
			}

		float debugMapMZ = 0;
		if (debugMapMZ != 0)
			{
			try{
			Writer out = new java.io.StringWriter();
			int middle = Math.round((debugMapMZ - _mzRange.min) * SpectrumResampler.getResampleFrequency());
			int start = Math.max(0,middle-100);
			int end = Math.min(height,middle+100);
			out.write("# scans [" + _scans[0].getNum() + "," + _scans[_scans.length-1].getNum() + "] mz [" +
                    (_mzRange.min + start/ SpectrumResampler.getResampleFrequency()) + "," +
                    (_mzRange.min + (end-1)/ SpectrumResampler.getResampleFrequency()) + "]" + "\n");
			float[][] source = wavelets;// original; //spectra; // wavelets : original;
			for (int z = end-1 ; z>=start ; z--)
				{
				for (int s = 0; s < source.length; s++)
					{
					out.write(Math.max(0,source[s][z]) + "\t");
					}
				out.write("\n");
				}
			out.write("\n");
			StringSelection sel = new StringSelection(out.toString());
			Clipboard clip = Toolkit.getDefaultToolkit().getSystemClipboard();
			clip.setContents(sel,sel);
			} catch (Exception x){
			}
			}

		Spectrum.Peak[] arr = (Spectrum.Peak[])peaks.toArray(new Spectrum.Peak[0]);
		return PeaksAsFeatures(arr);
		}


	Collection PeaksAsFeatures(Spectrum.Peak[] peaks)
		{
		ArrayList l = new ArrayList();
		for (int p = 0; p < peaks.length; p++)
			l.add(new Feature(peaks[p]));
		return l;
		}
	}
