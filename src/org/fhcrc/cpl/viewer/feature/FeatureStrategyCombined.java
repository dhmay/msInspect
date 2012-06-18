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
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.ExtractEdgeFeatures;
import org.fhcrc.cpl.viewer.feature.ExtractMaxima2D;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

import java.io.InputStream;
import java.util.*;

/**
 * User: mbellew
 * Date: Nov 2, 2004
 */
public class FeatureStrategyCombined extends FeatureStrategyUsingWindow // extends FeatureExtractor
	{
	static Logger _log = Logger.getLogger(FeatureStrategyCombined.class);

	static int _WindowMargin = 16;
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


	public FeatureStrategyCombined(MSRun run, int scanIndex, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scanIndex, count, maxCharge, range, sn);

		int c2 = Math.max(128, count + 2 * _WindowMargin);
		scanIndex = Math.max(0, scanIndex - (c2 - count) / 2);
		int scanMax = Math.min(scanIndex + c2, run.getScanCount());
		count = scanMax - scanIndex;
		_scans = getScans(run, scanIndex, count);

		_startNum = run.getScan(scanIndex).getNum();
		_endNum = run.getScan(scanMax-1).getNum();
		}


	public int getType()
		{
		return TYPE_2D;
		}


	public Feature[] _analyze() throws InterruptedException
		{
		Feature[] features = analyzeWindow(_scans, 256, _WindowMargin);
		ArrayList filtered = new ArrayList();
		for (int i = 0; i < features.length; i++)
			{
			Feature feature = features[i];
			if (feature.scan >= _startNum && feature.scan <= _endNum)
				filtered.add(feature);
			}
		return (Feature[])filtered.toArray(new Feature[filtered.size()]);
		}


	/**
	 * THIS IS THE MAIN FEATURE FINDING ROUTINE
	 */
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
		assert timerResample.stop();

		int width = spectra.length;
		int height = spectra[0].length;
		_log.debug("analyze2D datasize = " + (width * height * 4));

		// separate background and signal components
		// we're pretty conservative about what we call background,
		// we can use background and/or median later to filter peaks
		// further.

		timerBackground.start();
		float[][] background = Spectrum.RemoveBackground(spectra);
		float[][] median = new float[spectra.length][];
		for (int s=0 ; s<scans.length ; s++)
			median[s] = Spectrum.MedianWindow(spectra[s], height, 72, false);
		timerBackground.stop();

		//
		// Perform a high-level 2D analysis to detect peptide elution
		//
		// This is done using edge detection technique, this seems to
		// be effective at detecting even low intensity features, while
		// ignoring noise
		//

		timerEdgeFeatures.start();
		Spectrum.Peak[][] listPeaks = ExtractEdgeFeatures.analyze(spectra, _mzRange, 0);
		Spectrum.Peak[] grossFeatures = listPeaks[2];
		// UNDONE: can be make use of the edge matrices computed by ExtractEdgeFeatures
		for (int i = 0; i < grossFeatures.length; i++)
			{
			Spectrum.Peak peak = grossFeatures[i];
			int imz = (int)((peak.mz-_mzRange.min) * SpectrumResampler.getResampleFrequency());
			// TODO intensity of charge 0+ is not comparable to other peaks  
			//peak.intensity = spectra[peak.scan][imz];
			peak.background = background[peak.scan][imz];
			peak.setMedian(median[peak.scan][imz]);
			}
		timerEdgeFeatures.stop();

		//
		// STEP 2
		//
		// Analyse these areas for isotopic distributions, and report
		// peptide features.
		//
		// Ideally, we would do some sort of 2D deconvolution of the
		// elution/isotopic distribution.  However, I don't know how
		// do that.
		//
		// see ExtractPeaks()
		//

		assert timerExtractPeaks.start();

		// .analyze() is destructive, need to copy
		float[][] spectraT = new float[spectra.length][];
		for (int i = 0; i < spectraT.length; i++)
			spectraT[i] = (float[])spectra[i].clone();

		Spectrum.Peak[] peaksAll = ExtractMaxima2D.analyze(
				spectraT, _mzRange.min,  SpectrumResampler.getResampleInterval(),
				new FeatureStrategyWavelet2D.SmoothWavelet(),
		        //new Smooth2D(),
				-Float.MAX_VALUE);

		//
		// crude filtering, eliminate peaks not near grossFeatures
		//

		// OK go back through peaks and set intensity based on spectra
		ArrayList list = new ArrayList();
		for (int i = 0; i < peaksAll.length; i++)
			{
			Spectrum.Peak peak = peaksAll[i];
			int imz = (int)((peak.mz-_mzRange.min) * SpectrumResampler.getResampleFrequency());
			if (imz >= height || peak.scan >= width)
				continue;
			peak.intensity = spectra[peak.scan][imz];
			peak.background = background[peak.scan][imz];
			peak.setMedian(median[peak.scan][imz]);

			// TODO: this is a made up cut-off
			// TODO: still need to consider what a good value would be
			if (peak.intensity < (0.5*peak.background + 2*peak.getMedian() + 3))
				continue;
			list.add(peaksAll[i]);
			}
		peaksAll = (Spectrum.Peak[])list.toArray(new Spectrum.Peak[list.size()]);

		Arrays.sort(peaksAll, Spectrum.comparePeakIntensityDesc);// debug only
		Spectrum.Peak[] peaks = FeatureStrategyUsingWindow2D.FilterPeakSet(grossFeatures, peaksAll, null, 4);
		assert timerExtractPeaks.stop();

		assert timerExtractPeptides.start();
		for (int i = 0; i < peaks.length; i++)
			{
			Spectrum.Peak peak = peaks[i];
			_logDebug(peak.toString());
			}
        ArrayList peptidesAll = ExtractPeptideFeatures(_run, peaks);

		// we know that all peptides returned are near a grossFeature, because
		// we filtered them.  However, there may be gross features that we did not
		// find a good match for.  Add these back as 0+ features.

		Tree2D treeAll = new Tree2D();
		for (Iterator it = peptidesAll.iterator(); it.hasNext();)
			{
			Feature feature = (Feature)it.next();
			treeAll.add(feature.scan, feature.mz, feature);
			}
		for (int i = 0; i < grossFeatures.length; i++)
			{
			Spectrum.Peak f = grossFeatures[i];
			if (!treeAll.containsPoints(f.scan-4, f.mz-1, f.scan+4, f.mz+1))
				peptidesAll.add(f);
			}

		ArrayList fitleredPeptides = new ArrayList();
		// convert back to scanNum, and do a little more filtering
		for (Iterator iterator = peptidesAll.iterator(); iterator.hasNext();)
			{
			Feature feature = (Feature)iterator.next();
            if (feature.intensity < 2*feature.getMedian() + 0.5*feature.background + 1)
				continue;
			int index = feature.scan;
			feature.scan = scans[index].getNum();
			// fix up the retention time to match updated scan number
			feature.setTime((float)scans[index].getDoubleRetentionTime());

            fitleredPeptides.add(feature);
			}
		assert timerExtractPeptides.stop();

		assert timerAnalyze.stop();
		CPUTimer.DumpAllTimers();
		return fitleredPeptides;
		}


	Collection PeaksAsFeatures(Spectrum.Peak[] peaks)
		{
		ArrayList l = new ArrayList();
		for (int p = 0; p < peaks.length; p++)
			l.add(new Feature(peaks[p]));
		return l;
		}
	}
