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

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;

import java.util.ArrayList;
import java.util.Collection;

/**
 * User: mbellew
 * Date: Sep 7, 2004
 * Time: 12:10:06 PM
 * <p/>
 * Obviously this is WAY to slow to run for the whole file,
 * this is to just test an idea.
 */
public class FeatureStrategyUsingWindow2D extends FeatureStrategyUsingWindow //extends FeatureExtractor
	{
	int _startNum = 0;
	int _endNum = 0;

	public FeatureStrategyUsingWindow2D(MSRun run, int scanIndex, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scanIndex, count, maxCharge, range, sn);

		int c2 = Math.max(256, count);
		scanIndex = scanIndex - (c2 - count) / 2;
		scanIndex = Math.max(scanIndex, 0);
		int scanMax = Math.min(scanIndex + c2, run.getScanCount());
		count = scanMax - scanIndex;
		_scans = getScans(run, scanIndex, count);

    	_startNum = run.getScan(scanIndex).getNum();
		_endNum = run.getScan(scanIndex+count-1).getNum();
		}


	public Feature[] _analyze() throws InterruptedException
		{
		Feature[] features = analyzeWindow(_scans, 256, 16);
        ArrayList filtered = new ArrayList();
        for (int i = 0; i < features.length; i++)
            {
            Feature feature = features[i];
            if (feature.scan >= _startNum && feature.scan <= _endNum)
                filtered.add(feature);
            }
        return (Feature[]) filtered.toArray(new Feature[filtered.size()]);
		}


	protected Collection analyze2D(Scan[] scans)
		{
		//
		// Calculate maxima in the region of the scan to analyze
		//
		Spectrum.Peak[] interestingPeaks = null;
		if (true)
			interestingPeaks = InterestingFeatures(scans);

		Spectrum.Peak[] peaks = ExtractPeaks(scans);

		peaks = FilterPeakSet(interestingPeaks, peaks, _run, 4);

		return ExtractPeptideFeatures(_run, peaks);
		}


	public static Spectrum.Peak[] FilterPeakSet(Spectrum.Peak[] interestingPeaks, Spectrum.Peak[] peaks, MSRun run, int scanRange)
		{
		// eliminate peaks not near interesting areas
		if (null != interestingPeaks)
			{
			Tree2D tree = new Tree2D();
			for (int i=0 ; i<interestingPeaks.length ; i++)
				{
				Spectrum.Peak p = interestingPeaks[i];
				// CONSIDER: could use getRetentionTime() instead of getIndexForScanNum()
				int s = (run == null) ? p.scan : run.getIndexForScanNum(p.scan);
				tree.add(s, p.mz, p);
				}
			ArrayList usablePeaks = new ArrayList();
			for (int i=0 ; i<peaks.length ; i++)
				{
                Spectrum.Peak p = peaks[i];
				if (null == p)
					continue;
				int s = (run == null) ? p.scan : run.getIndexForScanNum(p.scan);
				float mz = p.getMz();
				if (tree.containsPoints(s-scanRange, mz-4.1F, s+scanRange, mz+1.1F))
					usablePeaks.add(p);
				}
			peaks = (Spectrum.Peak[])usablePeaks.toArray(new Spectrum.Peak[usablePeaks.size()]);
			}
		return peaks;
		}


	protected Spectrum.Peak[] InterestingFeatures(Scan[] scans)
		{
		Spectrum.Peak[][] t = ExtractEdgeFeatures.analyze(scans, _mzRange, 0);
		return t[2];
		}


	protected Spectrum.Peak[] ExtractPeaks(Scan[] scans)
		{
		return ExtractMaxima2D.analyze(scans, _mzRange);
		}


	static int _closest(Spectrum.Peak[] peaks, float x, int start)
		{
		float dist = Math.abs(x - peaks[start].mz);
		int i = start + 1;
		for (; i < peaks.length; i++)
			{
			float d = Math.abs(x - peaks[i].mz);
			if (d > dist)
				break;
			dist = d;
			}
		return i - 1;
		}

	/*
	 * score a candidate feature f
	 * <p/>
	 * Expecting charge=0 features and peaks
	 *
	 * @param f
	 * @param peaks
	public boolean ScoreFeature(Spectrum.Feature f, Spectrum.Peak[] peaks, double noise)
		{
		double sn = 2.0;
		int charge = f.charge;
		float invCharge = 1.0F / charge;
		// UNDONE: make peak count configurable
		// UNDONE: referenence distribution has 6 peaks
		float[] signal = new float[6];

		// find start position in peak list
		int p = Arrays.binarySearch(peaks, f, Spectrum.comparePeakMassAsc);
		if (p < 0) p = -(p + 1);
		p = Math.max(0, p - 1);
		int pLastFound = p;

		float mzP0 = f.mz;
		float sum = 0.0F;
		float intensityLast = 0.0F;
		boolean missingPeak = false;
		boolean skippedPeak = false;
		for (int Pi = 0; Pi < signal.length; Pi++)
			{
			float mzPi = mzP0 + Pi * invCharge;
			p = _closest(peaks, mzPi, p);
			Spectrum.Peak peakFound = peaks[p];
			float dist = Math.abs(peakFound.mz - mzPi);
			boolean found = !peakFound.excluded && dist < 2 * DELTA &&
			        (missingPeak ? peakFound.mz > sn * noise : peakFound.mz > noise);
			if (found)
				{
				// did we skip any significant peaks???
				if (missingPeak) skippedPeak = true;
				if (!skippedPeak)
					for (int s = pLastFound + 1; s < p; s++)
						{
						Spectrum.Peak skipped = peaks[s];
						if (skipped.intensity > intensityLast / 4)
							skippedPeak = true;
						}
				signal[Pi] = peaks[p].intensity;
				pLastFound = p;
				}
			else
				{
				signal[Pi] = 0.1F;
				missingPeak = true;
				}
			intensityLast = signal[Pi];
			sum += signal[Pi];
			}
		signal[0] /= sum;
		signal[1] /= sum;
		signal[2] /= sum;
		signal[3] /= sum;
		signal[4] /= sum;
		signal[5] /= sum;
		f.intensity = sum;
		f.kl = Spectrum.KLPoissonDistance(f.mz * charge, signal);
		return skippedPeak;
		}


	static double distanceNearestFraction(double x, double denom)
		{
		// scale so we can do distance to nearest integer
		double t = x * denom;
		double delta = Math.abs(t - Math.round(t));
		return delta / denom;
		}
	 */


	public int getType()
		{
		return TYPE_2D;
		}
	}
