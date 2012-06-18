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

import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

/**
 * User: mbellew
 * Date: Sep 7, 2004
 * Time: 12:10:06 PM
 */


public class FeatureStrategySpectrumFit extends FeatureExtractor
	{
	static float[] peakShape = new float[2* SpectrumResampler.getResampleFrequency()-1];

	static
		{
		// CONSIDER: peaks are really more gamma shaped than normal
		double sigma = SpectrumResampler.getResampleFrequency()/15.0;
		double s2 = -1.0/(2.0 * sigma * sigma);
		for (int i=0, x=-(SpectrumResampler.getResampleFrequency()-1); i<peakShape.length ; i++, x++)
			peakShape[i] = (float)Math.exp(x*x*s2);
		}

	int _freq;
	Scan[] _scans;

	public FeatureStrategySpectrumFit(MSRun run, int scan, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scan, count, maxCharge, range, sn);
		this._freq = SpectrumResampler.getResampleFrequency();
		_scans = getScans(run, scan, count);
		}


	public Feature[] _analyze() throws InterruptedException
		{
		return analyzeScanAtATime(_scans);
		}


	public Collection analyze1D(Scan scan)
		{
		ArrayList featureList = new ArrayList();
		float[][] _zeroChargeSpectrum = Spectrum.ResampleSpectrum(scan.getSpectrum(), _mzRange, SpectrumResampler.getResampleFrequency(), true);
		float noise = Spectrum.Noise(_zeroChargeSpectrum[1], 0, _zeroChargeSpectrum[1].length);

		float[][] smoothSpectrum = new float[][]
			{
			_zeroChargeSpectrum[0],
			Spectrum.FFTsmooth(_zeroChargeSpectrum[1], 8.0, true)
			};

		//
		// Generate Peaks
		//
		Spectrum.Peak[] peaks = Spectrum.PickPeaks(smoothSpectrum, _sn * noise);
		Spectrum.Peak[] peaksByIntensity = new Spectrum.Peak[peaks.length];
		System.arraycopy(peaks, 0, peaksByIntensity, 0, peaks.length);
		Arrays.sort(peaksByIntensity, Spectrum.comparePeakIntensityDesc);

		//ArrayList candidates = new ArrayList();
		//float[] primaryPeaks = new float[12];
		Spectrum.Peak findPeak = new Spectrum.Peak(0, 0);

		for (int i = 0; i < peaksByIntensity.length; i++)
			{
			Spectrum.Peak peakHighest = peaksByIntensity[i];
			if (peakHighest.excluded)
				continue;
			if (peakHighest.intensity < _sn * noise)
				break;

			_logDebug("--------------");
			_logDebug("peak " + peakHighest.mz + " " + peakHighest.intensity);


			float mzWindowStart = peakHighest.mz - 2.1F;
			Feature bestFeature = null;

			//
			// generate candidate features
			//

			findPeak.mz = mzWindowStart;
			int windowStartIndex = Arrays.binarySearch(peaks, findPeak, Spectrum.comparePeakMzAsc);
			if (windowStartIndex < 0) windowStartIndex = -(windowStartIndex + 1);
			for (int j = windowStartIndex; j < peaks.length; j++)
				{
				Spectrum.Peak p = peaks[j];
				if (p.mz > peakHighest.mz)
					break;
				if (p.excluded)
					continue;
				// even at high end of mass range (6400amu) leading peak is not expected
				// to be smaller than about 1/6 highest peak, and typically much larger
				if (p.intensity < peakHighest.intensity / 10.0F)
					continue;
				double distance = peakHighest.mz - p.mz;
				for (int charge = _maxCharge; charge >= 1; charge--)
					{
					double r = distanceNearestFraction(distance, charge);
					if (r < 2 * (SpectrumResampler.getResampleInterval() + .001))
						{
						// get index of peak p in _zeroChargeSpectrum
						int pIndex = (int)Math.round((p.mz - _zeroChargeSpectrum[0][0]) * SpectrumResampler.getResampleFrequency());
						assert _zeroChargeSpectrum[0][pIndex] == p.mz;

						Score scoreBest = null;
						float mzBest = p.mz;
						for (int start = pIndex-1; start <= pIndex+1 ; start++)
							{
							if (start < 0 || start >= _zeroChargeSpectrum[0].length)
								continue;
                            Score score = Score(_zeroChargeSpectrum, start, charge);
							if (null == scoreBest || score.score > scoreBest.score)
								{
								scoreBest = score;
								mzBest = _zeroChargeSpectrum[0][start];
								}
							}

						if (null == bestFeature || scoreBest.score > bestFeature.intensity)
							{
							Feature f = new Feature(scan.getNum(), mzBest, (float) scoreBest.score);
							f.kl = (float)scoreBest.diff;
                            f.charge = charge;
							bestFeature = f;
							}
						}
					}
				}

			_logDebug("feature " + bestFeature.mz + " " + bestFeature.charge + "+ " + (float)bestFeature.kl);
			featureList.add(bestFeature);

			//
			// exclude all peaks explained by fBest
			//

			// UNDONE: I DON'T HANDLE OVERLAPPING FEATURES YET
			// UNDONE: just kill all peaks in the area!

			findPeak.mz = bestFeature.mz - 4/bestFeature.charge;
			int startFeature = Arrays.binarySearch(peaks, findPeak, Spectrum.comparePeakMzAsc);
			if (startFeature < 0) startFeature = -(startFeature + 1);
			for (int j = startFeature; j < peaks.length; j++)
				{
				Spectrum.Peak p = peaks[j];
				double diffMz = p.mz - bestFeature.mz;
				int peakFound = (int)Math.round(diffMz * bestFeature.charge);
				if (peakFound > 5)
					break;
				p.excluded = true;
				}
			}

		// adjust back to +1 state
		for (int i = 0; i < featureList.size(); i++)
			((Feature)featureList.get(i)).mz += Spectrum.HYDROGEN_ION_MASS;
		return featureList;
		}


	static final int peaksExpected = 6;
	float[] expectedSignal = new float[peaksExpected * SpectrumResampler.getResampleFrequency() + 1];
	int lenExpectedSignal = 0;
	double massExpected = 0.0;
	int chargeExpected = 0;

	static class Score
		{
		double diff;
		double intensity;
		double score;
		}

	Score Score(float[][] spectrum, int start, int charge)
		{
		assert 0 == SpectrumResampler.getResampleFrequency() % charge;
		assert 1 == expectedSignal.length % 2;
		assert 1 == peakShape.length % 2;

		float mz = spectrum[0][start];
		float mass = mz * charge;

		// recompute expected signal if necessary
		if (10 > Math.abs(mass - massExpected) || charge != chargeExpected)
			{
			float[] poisson = Spectrum.Poisson(mass);
			Arrays.fill(expectedSignal, 0.0F);
            for (int p=0 ; p<peaksExpected ; p++)
				{
				int expectedPeakCenter = p* SpectrumResampler.getResampleFrequency()/charge + SpectrumResampler.getResampleFrequency()/2;
				int peakShapeCenter = peakShape.length / 2;
				int offStart = -Math.min(expectedPeakCenter,peakShapeCenter);
				int offEnd = Math.min(peakShapeCenter, expectedSignal.length - expectedPeakCenter -1);
				float intensity = poisson[p];
				for (int off = offStart ; off <= offEnd ; off++)
					{
					expectedSignal[expectedPeakCenter+off] +=
							peakShape[peakShapeCenter+off] * intensity;
					}
				}
			massExpected = mass;
			chargeExpected = charge;
			lenExpectedSignal = SpectrumResampler.getResampleFrequency() + (peaksExpected-1)* SpectrumResampler.getResampleFrequency()/charge + 1;

			if (false)
				{
				_logDebug("\n----------\n\nmass = " + mass + "\ncharge = " + charge + "\n");
				for (int i = 0 ; i <lenExpectedSignal ; i++)
					_logDebug("" + expectedSignal[i]);
				}
			}

		Score s = SignalDiff(spectrum[1], start - SpectrumResampler.getResampleFrequency()/2, expectedSignal, 0, lenExpectedSignal);
		_logDebug("" +  mz  + "\t" + charge + "+\t" + (float)s.diff + "\t" + (float)s.score);
		return s;
		}


	Score SignalDiff(float[] signal1, int off1, float[] signal2, int off2, int len)
		{
		// handle edge cases -- find overlapping region of signals
		int t = Math.min(off1,off2);
		if (t < 0)
			{
			off1 -= t;
			off2 -= t;
			len += t;
			}
		t = Math.min(signal1.length - (off1+len), signal2.length - (off2+len));
		if (t < 0)
			{
			len += t;
			}

		double sum1 = 0.0, sum2 = 0.0;
		for (int i = 0 ; i < len ; i++)
			{
			sum1 += signal1[off1+i];
			sum2 += signal2[off2+i];
			}

		if (false)	// simple diff
			{
			double r = sum2/sum1;
			double diff = 0.0;
			for (int i = 0 ; i<len ; i++)
				diff += Math.abs(signal1[off1+i]*r - signal2[off2+i]);
			diff /= sum2;
			diff /= 2;

			Score s = new Score();
			s.diff = diff;
			s.intensity = sum1;
			s.score = sum1 * (1-diff);
			return s;
			}
		else
			{
			double diff = Spectrum.KLDistanceSymmetric(signal1, off1, signal2, off2, len);

			Score s = new Score();
			s.diff = diff;
			s.intensity = sum1;
			s.score = sum1 * Math.exp(-diff);
			return s;
			}
		}


	static double distanceNearestFraction(double x, double denom)
		{
		// scale so we can do distance to nearest integer
        double t = x*denom;
		double delta = Math.abs(t - Math.round(t));
		return delta/denom;
		}
	}
