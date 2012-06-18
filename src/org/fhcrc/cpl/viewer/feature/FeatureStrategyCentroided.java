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
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;

import java.util.Arrays;

/**
 * User: mbellew
 * Date: Oct 6, 2004
 * Time: 9:30:21 PM
 */
public class FeatureStrategyCentroided extends FeatureStrategyUsingWindow
	{
	public FeatureStrategyCentroided(MSRun run, int scan, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scan, count, maxCharge, range, sn);
		}


	protected Spectrum.Peak[] _pickPeaks(Scan scan)
		{
		float[][] centroided = scan.getSpectrum();
		Spectrum.Peak[] cleanPeaks = clean(centroided);
		for (int i=0 ; i<cleanPeaks.length ; i++)
			cleanPeaks[i].scan = scan.getNum();
		return cleanPeaks;
		}


	// UNDONE: parameter
	static int MAX_CHARGE = 6;


	/** for display/debugging */
	public static float[][] cleanSpectrum(float[][] centroided)
		{
		int length = centroided[0].length;
		Spectrum.Peak[] rawPeaks = new Spectrum.Peak[length];
		for (int i=0 ; i<length ; i++)
            rawPeaks[i] = new Spectrum.Peak(1, centroided[0][i], centroided[1][i]);

		Spectrum.Peak[] cleanPeaks = clean(centroided);

		float[][] clean = new float[2][cleanPeaks.length];
		for (int i=0 ; i<cleanPeaks.length ; i++)
			{
			clean[0][i] = cleanPeaks[i].mz;
			clean[1][i] = cleanPeaks[i].intensity;
			}
		return clean;
		}


	/** for display/debugging */
	public static float[][] backgroundSpectrum(float[][] centroided)
		{
		FloatRange range = new FloatRange();
		float[][] clean = background(centroided, range);
		return clean;
		}


	static int BUCKET_SIZE = MAX_CHARGE;

	/**
	 * compute background based on 'excess' peaks in the spectrum
	 *
	 * I don't really understand this kind of background, and this is
	 * based on observing only one machine type, so mileage may vary....
	 */
	static float[][] background(float[][] raw, FloatRange range)
		{
		int length = raw[0].length;
		range.min = (float)Math.floor(raw[0][0]);
		range.max = (float)Math.ceil(raw[0][length-1]);

		float[][] bg = new float[2][(int)((range.max-range.min)*BUCKET_SIZE + 1)];
		for (int i=0; i<bg[0].length ; i++)
			bg[0][i] = range.min + i * 1.0F/BUCKET_SIZE;

//		Spectrum.Peak[] byIntensity = (Spectrum.Peak[])raw.clone();
//        Arrays.sort(byIntensity, Spectrum.comparePeakIntensityDesc);

		int start = 0;
		int end = 0;
		float[] tmp = new float[100];
		for ( ; start < length ; start++)
			{
			float mzStart = raw[0][start];
			for ( ; end < length && raw[0][end] < mzStart + 1 ; end++ )
				;
			int count = end - start;

            if (count <= MAX_CHARGE)
	            continue;

			tmp = Spectrum.realloc(tmp, count);
			System.arraycopy(raw[1], start, tmp, 0, count);

			Arrays.sort(tmp, 0, count);
			Spectrum.Reverse(tmp, 0, count);

			// don't want to be too abrupt so use some sort of sliding scale
			float denseCutoff = 0;
			if (count > 2 * MAX_CHARGE)
				denseCutoff = tmp[2*MAX_CHARGE];
			else
				{
				float ratio = (2F*MAX_CHARGE-count) / MAX_CHARGE;
				denseCutoff = Math.max(0,tmp[count-1] * ratio);
				}

			for (int j=start ; j<end ; j++)
				{
				float in = raw[1][j];
				float mz = raw[0][j];
				int b = Math.round(((mz - range.min) * BUCKET_SIZE));
				bg[1][b] = Math.max(bg[1][b], Math.min(in, denseCutoff));
				}
			// in areas of high peak density this can get expensive, skip ahead a bit
			start++;
			}

		Spectrum.SmoothALittle(bg[1]);
		Spectrum.SmoothALittle(bg[1]);
		return bg;
		}


	/** centroided spectra are usually not very clean
	 *
	 * We are going to clean it up by doing our own peak detection
	 * and then matching up our peaks with the originals
	 *
	 * @param raw
	 * @return
	 */
	static Spectrum.Peak[] recentroid(float[][] raw)
		{
		int length = raw[0].length;
		FloatRange range = new FloatRange();
		range.min = (float)Math.floor(raw[0][0]);
		range.max = (float)Math.ceil(raw[0][length-1]);

        float[][] spectrum = Spectrum.ResampleSpectrum(raw, range, SpectrumResampler.getResampleFrequency(), false);
		Spectrum.Peak[] peaks = Spectrum.WaveletPeaks(spectrum);

		float delta = 0.6F/BUCKET_SIZE; // just larger than 1/2*BUCKET_SIZE
		for (int i=0 ; i<peaks.length ; i++)
			{
			Spectrum.Peak p = peaks[i];
			p.intensity = 0;
			int j = Arrays.binarySearch(raw[0], p.mz - delta);
			if (j < 0) j = -(j+1);
            for ( ; j < length && raw[0][j] < p.mz + delta ; j++)
	            {
                if (raw[1][j] > p.intensity)
	                {
	                p.mz = raw[0][j];
	                p.intensity = raw[1][j];
	                }
	            }
//			assert p.intensity > 0;
			}
		return peaks;
		}


	static Spectrum.Peak[] clean(float[][] raw)
		{
		FloatRange range = new FloatRange();
		float[][] bg = background(raw, range);

		Spectrum.Peak[] peaks = recentroid(raw);

		for (int i=0 ; i<peaks.length ; i++)
			{
			Spectrum.Peak p = peaks[i];
			int b = Math.round(((p.mz-range.min) * BUCKET_SIZE));
			p.background = bg[1][b];
			if (true)
				{
				// remove background, filter
				p.intensity -= p.background;
				if (p.intensity <= p.background)
					p.excluded = true;
				}
			if (i > 0)
				{
				// hmm, conceivably we could double up two peaks
				// UNDONE: do something more clever???
				if (p.mz == peaks[i-1].mz)
					p.excluded = true;
				}
			}

		// return non-excluded peaks
		int count = 0;
		for (int i = 0 ; i<peaks.length ; i++)
			if (!peaks[i].excluded)
				count++;
        Spectrum.Peak[] clean = new Spectrum.Peak[count];
		int end = 0;
		for (int i = 0 ; i<peaks.length ; i++)
			if (!peaks[i].excluded)
				clean[end++] = peaks[i];
		return clean;
		}


	static Spectrum.Peak tmp = new Spectrum.Peak();
	static int find(Spectrum.Peak[] a, float mz)
		{
		synchronized(tmp)
			{
			tmp.mz = mz;
			int i = Arrays.binarySearch(a, tmp, Spectrum.comparePeakMzAsc);
			return i < 0 ? -(i+1) : i;
			}
		}
	}
