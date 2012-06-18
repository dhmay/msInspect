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

import org.fhcrc.cpl.viewer.util.Haar;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.datastructure.IntegerArray;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

import java.awt.*;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Properties;

/**
 * User: mbellew
 * Date: Sep 15, 2004
 * Time: 6:42:33 PM
 */
public class ExtractEdgeFeatures
	{

	static int _FilterSizeV;
	static int _FilterSizeH;
	static float _HaarThresholdV;
	static float _HaarMinimumV;
	static float _HaarThresholdH;
	static float _HaarMinimumH;
	static int _FeatureWindowV_width;
	static int _FeatureWindowV_height;
	static int _FeatureWindowV_count;
	static int _FeatureWindowH_width;
	static int _FeatureWindowH_height;
	static int _FeatureWindowH_count;


	static
	{
	initProps();
	}


	private static boolean bPropsInitialized = false;

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
			_FilterSizeV = Integer.parseInt(((String)props.get("edge.FilterSizeV")).trim());
			_FilterSizeH = Integer.parseInt(((String)props.get("edge.FilterSizeH")).trim());
			_HaarThresholdV = Float.parseFloat(((String)props.get("edge.ThresholdV")).trim());
			_HaarMinimumV = Float.parseFloat(((String)props.get("edge.MinimumV")).trim());
			_HaarThresholdH = Float.parseFloat(((String)props.get("edge.ThresholdH")).trim());
			_HaarMinimumH = Float.parseFloat(((String)props.get("edge.MinimumH")).trim());

			String t = (String)props.get("edge.FeatureWindowH");
			String[] a = t.split(",");
			_FeatureWindowH_width = Integer.parseInt(a[0].trim());
			_FeatureWindowH_height = Integer.parseInt(a[1].trim());
			_FeatureWindowH_count = Integer.parseInt(a[2].trim());

			t = (String)props.get("edge.FeatureWindowV");
			a = t.split(",");
			_FeatureWindowV_width = Integer.parseInt(a[0].trim());
			_FeatureWindowV_height = Integer.parseInt(a[1].trim());
			_FeatureWindowV_count = Integer.parseInt(a[2].trim());
			}
		catch (java.io.IOException x)
			{
			x.printStackTrace();
			}
		}


	private ExtractEdgeFeatures()
		{
		}


	/**
	 * scans should have a margin of at least 4 scans on either side, in which we won't detect features
	 *
	 * @param scans
	 * @param rangePeaks
	 * @param threshold
	 * @return
	 */
	public static Spectrum.Peak[][] analyze(Scan[] scans, FloatRange rangePeaks, float threshold)
		{
		threshold = 0; // UNDONE: ignoring passed in value

		// use larger range for resample to avoid edge checking later
		FloatRange range = new FloatRange(rangePeaks.min - .5F, rangePeaks.max + 1.0F);

		float[][] spectra = new float[scans.length][];
		for (int i = 0; i < scans.length; i++)
			{
			spectra[i] = Spectrum.Resample(scans[i].getSpectrum(), range, SpectrumResampler.getResampleFrequency());
			}
		Spectrum.Peak[][] listPeaks =  analyze(spectra, range, threshold);

		// translate back to scanNum
		for (int i = 0; i < listPeaks.length; i++)
			{
			Spectrum.Peak[] peaks = listPeaks[i];
			for (int j = 0; j < peaks.length; j++)
				{
				Spectrum.Peak peak = peaks[j];
				int index = peak.scan;
				peak.scan = scans[index].getNum();
				}
			}
		return listPeaks;
		}


	public static Spectrum.Peak[][] analyze(float[][] spectra, FloatRange range, float threshold)
		{
		initProps();
		int scanMax = spectra.length;
		int imzMax = spectra[0].length;

		//
		// find maxima
		//

		char[][] edgeV = new char[scanMax][imzMax];

		ArrayList peaksV = new ArrayList();
		float[] haar = new float[scanMax];
		float[] elution = new float[scanMax];
		for (int imz = 0; imz < imzMax; imz++)
			{
			// create elution profile array
			Spectrum.getRow(spectra, imz, elution);
			Haar.transform(elution, _FilterSizeV, haar);
			float med = Spectrum.MedianSampled(haar, true);
			float t = med * _HaarThresholdV;
			for (int i = 0; i < haar.length; i++)
				haar[i] -= t;

			t = Math.max(threshold, med * _HaarMinimumV);
			int[] elutionPeaks = ZeroIndexes(haar, t);

			for (int p = 0; p < elutionPeaks.length; p++)
				{
				int s = elutionPeaks[p];
				float mz = range.min + imz * (float) SpectrumResampler.getResampleInterval();
				edgeV[s][imz] = 1;
				peaksV.add(new Spectrum.Peak(s, mz, haar[s]));
				}
			}


		char[][] edgeH = new char[spectra.length][spectra[0].length];
		ArrayList peaksH = new ArrayList();
		haar = new float[spectra[0].length];
		for (int s = 0; s < scanMax; s++)
			{
			haar = Haar.transform(spectra[s], _FilterSizeV);
			float med = Spectrum.MedianSampled(haar, true);
			float t = med * _HaarThresholdH;
			for (int i = 0; i < haar.length; i++)
				haar[i] -= t;

			t = Math.max(threshold, med * _HaarMinimumH);
			int[] scanPeaks = ZeroIndexes(haar, t);
			for (int p = 0; p < scanPeaks.length; p++)
				{
				int imz = scanPeaks[p];
				float mz = range.min + imz * SpectrumResampler.getResampleInterval();
				edgeH[s][imz] = 1;
				peaksH.add(new Spectrum.Peak(s, mz, haar[imz]));
				}
			}

		//
		// CROSSINGS
		//

		// there are a lot of ways to do this,
		// I'm going to assume that edges are very sparce
		// making this strategy reasonable
		char[][] edgeVx = new char[edgeV.length][edgeV[0].length];
		int sFromOffset = - _FeatureWindowV_width/2;
		int sToOffset = _FeatureWindowV_width + sFromOffset;
		int imzFromOffset = sFromOffset;
		int imzToOffset = _FeatureWindowV_height - imzFromOffset;
		for (int s=0; s<scanMax ; s++)
			{
			char[] v = edgeV[s];
			for (int imz = 0; imz < imzMax; imz++)
				{
				if (v[imz] == 0)
					continue;
				// found a point add to edgeVx window
				int sfrom = Math.max(0, s+sFromOffset);
				int sto = Math.min(scanMax, s+sToOffset);
				int imzfrom = Math.max(0, imz+imzFromOffset);
				int imzto = Math.min(imzMax, imz+imzToOffset);
				for (int sadd=sfrom ; sadd<sto ; sadd++)
					for (int imzadd=imzfrom ; imzadd <imzto ; imzadd++)
						edgeVx[sadd][imzadd] += 1;
				}
			}

		char[][] edgeHx = new char[edgeH.length][edgeH[0].length];
		imzFromOffset = -_FeatureWindowH_height/2;
		imzToOffset = _FeatureWindowH_height + imzFromOffset;
		sFromOffset = imzFromOffset;
		sToOffset = _FeatureWindowH_width + sFromOffset;
		for (int s=0; s<scanMax ; s++)
			{
			char[] h = edgeH[s];
			for (int imz = 0; imz < imzMax; imz++)
				{
				if (h[imz] == 0)
					continue;
				// found a point add to edgeVx window
				int sfrom = Math.max(0, s+sFromOffset);
				int sto = Math.min(scanMax, s+sToOffset);
				int imzfrom = Math.max(0, imz+imzFromOffset);
				int imzto = Math.min(imzMax, imz+imzToOffset);
				for (int sadd=sfrom ; sadd<sto ; sadd++)
					for (int imzadd=imzfrom ; imzadd <imzto ; imzadd++)
						edgeHx[sadd][imzadd] += 1;
				}
			}

		edgeV = edgeH = null;

		/*
		for (int s = 0; s < scanMax; s++)
			{
			for (int imz = 0; imz < imzMax; imz++)
				{
				int x = Math.min(s + _FeatureWindowV_width, edgeV.length);
				for (int k = s + 1; k < x; k++)
					edgeV[s][imz] += edgeV[k][imz];
				/*
				int n = s - _FeatureWindowV_width/2;
				int x = n + _FeatureWindowV_width;
				n = Math.max(n, 0);
				x = Math.min(x, edgeV.length);
				for (int k = s + 1; k < x; k++)
					edgeV[s][imz] += edgeV[k][imz];
				* /
				}
			}
		for (int s = 0; s < scanMax; s++)
			{
			for (int imz = 0; imz < imzMax; imz++)
				{
				int x = Math.min(imz + _FeatureWindowV_height, edgeV[0].length);
				for (int k = imz + 1; k < x; k++)
					edgeV[s][imz] += edgeV[s][k];
				}
			}
		for (int s = 0; s < scanMax; s++)
			{
			for (int imz = 0; imz < imzMax; imz++)
				{
				int x = Math.min(s + _FeatureWindowH_width, edgeH.length);
				for (int k = s + 1; k < x; k++)
					edgeH[s][imz] += edgeH[k][imz];
				}
			}
		for (int s = 0; s < scanMax; s++)
			{
			for (int imz = 0; imz < imzMax; imz++)
				{
				int x = Math.min(imz + _FeatureWindowH_height, edgeH[0].length);
				for (int k = imz + 1; k < x; k++)
					edgeH[s][imz] += edgeH[s][k];
				}
			}
		*/

		// combine into one array of 'good' points
		char [][] edgeFeatures = new char[scanMax][imzMax];
		for (int s = 0; s < scanMax; s++)
			for (int imz = 0; imz < imzMax; imz++)
				{
				if (edgeVx[s][imz] < _FeatureWindowV_count || edgeHx[s][imz] < _FeatureWindowH_count)
					continue;
				else
					{
					// this function is not that important,
					// it's just used to break ties in loop below
					edgeFeatures[s][imz] = (char)(edgeVx[s][imz] + edgeHx[s][imz]);
					}
				}

		// combine adjacent points into one (more or less?) point
		for (int s = 1; s < scanMax-1; s++)
			for (int imz = 1; imz < imzMax-1; imz++)
				{
				char v = edgeFeatures[s][imz];
				if (v == 0)
					continue;
				// find max value in this island
				collapseIsland(edgeFeatures, s, imz);
				}

		ArrayList peaksX = new ArrayList();
		for (int s = 1; s < scanMax - 1; s++)
			{
			for (int imz = 1; imz < imzMax - 1; imz++)
				{
				if (0 == edgeFeatures[s][imz])
					continue;
				int mzIndex = imz; //  + _FeatureWindowH_height / 2;
				float mz = range.min + mzIndex * SpectrumResampler.getResampleInterval();
				//if (mz < rangePeaks.min || mz > rangePeaks.max)
				//	continue;
				int scanIndex = s; //  + _FeatureWindowV_width / 2;
				// this is only a rough estimate of intensity,
				// doesn't account for background
				// doesn't find highest peaks in area
				float intensity = (spectra[scanIndex][imz] +
				        spectra[scanIndex][imz + 1] +
				        spectra[scanIndex + 1][imz] +
				        spectra[scanIndex + 1][imz + 1]) / 4;
				//int scanNum = scans[s].getNum();
				Feature f = new Feature(s, mz, intensity);
                f.kl = 0;
				f.setDescription("" + (int)edgeVx[s][imz] + "x" + (int)edgeHx[s][imz]);
				peaksX.add(f);
				}
			}

		Spectrum.Peak[] v = (Spectrum.Peak[])peaksV.toArray(new Spectrum.Peak[peaksV.size()]);
		Spectrum.Peak[] h = (Spectrum.Peak[])peaksH.toArray(new Spectrum.Peak[peaksH.size()]);
		Feature[] x = (Feature[])peaksX.toArray(new Feature[peaksX.size()]);
		return new Spectrum.Peak[][]{v, h, x};
		}


	/** this simple function with blow up if the islands aren't relatively small */
	private static void collapseIsland(char[][] edgeFeatures, int s, int imz)
		{
        Point bestPoint = new Point(s, imz); // out parameter
		char best = edgeFeatures[s][imz];
		//float best = spectra[s][imz];
		best = _collapse(edgeFeatures, best, s, imz, bestPoint);
		edgeFeatures[bestPoint.x][bestPoint.y] = best;
		}


	private static char _collapse(char[][] edgeFeatures, char best, int s, int imz, Point bestPoint)
		{
		if (s < 0 || s >= edgeFeatures.length)
			return best;
		if (imz < 0 || imz >= edgeFeatures[0].length)
			return best;
		char f = edgeFeatures[s][imz];
		if (f == 0)
			return best;
		if (f > best)
			{
			best = f;
			bestPoint.x = s;
			bestPoint.y = imz;
			}
		edgeFeatures[s][imz] = 0;
		best = _collapse(edgeFeatures, best, s-1, imz-1, bestPoint);
		best = _collapse(edgeFeatures, best, s-1, imz, bestPoint);
		best = _collapse(edgeFeatures, best, s-1, imz+1, bestPoint);
		best = _collapse(edgeFeatures, best, s, imz-1, bestPoint);
		best = _collapse(edgeFeatures, best, s, imz+1, bestPoint);
		best = _collapse(edgeFeatures, best, s+1, imz-1, bestPoint);
		best = _collapse(edgeFeatures, best, s+1, imz, bestPoint);
		best = _collapse(edgeFeatures, best, s+1, imz+1, bestPoint);
		return best;
		}


	public static FeatureSet[] EdgeFeatures(Scan[] scans, FloatRange rangePeaks, float threshold)
		{
		Spectrum.Peak[][] p = analyze(scans, rangePeaks, threshold);
		FeatureSet v = new FeatureSet(p[0], new Color(0x00, 0x00, 0xcc));
		FeatureSet h = new FeatureSet(p[1], new Color(0x00, 0xcc, 0x00));
		FeatureSet x = new FeatureSet(p[2], new Color(0xff, 0x00, 0x00));
 		return new FeatureSet[]{v, h, x};
		}


	/* find zero crossings following large peaks */
	public static int[] ZeroIndexes(float[] signal, double minFilter)
		{
		IntegerArray arr = new IntegerArray();
		int len = signal.length;

		float prev = -Float.MAX_VALUE;
		float curr;
		boolean found = false;
		for (int i = 0; i < len;)
			{
			// rising
			for (; i < len && prev <= (curr = signal[i]); i++)
				prev = curr;
			if (prev > minFilter)
				found = true; // true until we find the next zero crossing
			// falling
			for (; i < len && prev >= (curr = signal[i]); i++)
				{
				if (found && curr <= 0)
					{
					arr.add(i - 1);
					found = false;
					}
				prev = curr;
				}
			found = false;
			}
		return arr.toArray(null);
		}
	}
