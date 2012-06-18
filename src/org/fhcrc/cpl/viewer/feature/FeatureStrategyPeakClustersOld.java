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

import modwt.Filter;
import modwt.Transform;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;

import java.io.InputStream;
import java.util.*;

/**
 * User: mbellew
 * Date: Nov 2, 2004
 */
public class FeatureStrategyPeakClustersOld extends FeatureStrategyUsingWindow // extends FeatureExtractor
	{
	static Logger _log = Logger.getLogger(FeatureStrategyPeakClustersOld.class);

	static int _WindowMargin = 64;

	static int _FeatureScanWindowStart;
	static int _FeatureScanWindowWidth;
	static float _FeatureMzWindowStart;
	static float _FeatureMzWindowHeight;
	static float _AverageWindowWidth;

	static CPUTimer timerAnalyze = new CPUTimer("FeatureStrategyPeakClustersOld.analyze");
	static CPUTimer timerResample = new CPUTimer("FeatureStrategyPeakClustersOld.get and resample");
	static CPUTimer timerBackground = new CPUTimer("FeatureStrategyPeakClustersOld.background");
	static CPUTimer timerExtractPeaks = new CPUTimer("FeatureStrategyPeakClustersOld.peaks");
	static CPUTimer timerExtractPeptides = new CPUTimer("FeatureStrategyPeakClustersOld.peptides");

	static
	{
	initProps();
	}

	private static boolean bPropsInitialized = false;

    //start and end scan numbers
	int _startNum = 0;
	int _endNum = 0;


    Scan[] _scans = null;

	public boolean debugReturnPeaks = false;

        /**
         * Load properties from feature.properties
         */
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
			_FeatureScanWindowStart =
                    Integer.parseInt(((String)props.get("feature.ScanWindowStart")).trim());
			_FeatureScanWindowWidth =
                    Integer.parseInt(((String)props.get("feature.ScanWindowWidth")).trim());
			_FeatureMzWindowStart =
                    Float.parseFloat(((String)props.get("feature.MzWindowStart")).trim());
			_FeatureMzWindowHeight =
                    Float.parseFloat(((String)props.get("feature.MzWindowHeight")).trim());
			_AverageWindowWidth =
                    Integer.parseInt(((String)props.get("feature.AverageWindowWidth")).trim());
			}
		catch (java.io.IOException x)
			{
			x.printStackTrace();
			}
		}


	public FeatureStrategyPeakClustersOld(MSRun run, int scanIndex, int count,
                                       int maxCharge, FloatRange range, double sn)
		{
		super(run, scanIndex, count, maxCharge, range, sn);

		int c2 = Math.max(256, count + 2 * _WindowMargin);
		scanIndex = Math.max(0, scanIndex - (c2 - count) / 2);
		int scanMax = Math.min(scanIndex + c2, run.getScanCount());
		count = scanMax - scanIndex;
		_scans = getScans(run, scanIndex, count);

		_startNum = run.getScan(scanIndex).getNum();
		_endNum = run.getScan(scanMax - 1).getNum();
		}


	public int getType()
		{
		return TYPE_2D;
		}


	public static float[] _thresholdElution(float[] elution)
		{
		// UNDONE: cache intermediate arrays
		int K = 3;
		int N = elution.length;
		Pair<float[][], float[][]> tmp = new Pair<float[][], float[][]>(null, null);

        tmp.first = Transform.decompose(elution, N, K, new Filter("haar"), "modwt",
                                        "periodic", tmp.first);
		tmp.second = Transform.multiresolution(tmp.first, N, K, new Filter("haar"), "modwt",
                                               "periodic", tmp.second);
		float[][] mra = tmp.second;
		float[] threshold;

//			{
//			threshold = new float[N];
//			for (int i = 0; i < N; i++)
//				{
//				threshold[i] = mra[2][i] + mra[3][i]; // + mra[4][i] + mra[5][i] + mra[6][i];
//				}
//			}

        // even smoother
        threshold = mra[3];
		return threshold;
		}


	public Feature[] _analyze() throws InterruptedException
		{
		Feature[] features = analyzeWindow(_scans, 256, _WindowMargin);
		List<Feature> filtered = new ArrayList<Feature>();
		for (Feature feature : features)
			{
			if (feature.scan >= _startNum && feature.scan <= _endNum)
				filtered.add(feature);
			}
		return filtered.toArray(new Feature[filtered.size()]);
		}


	//
	// So what's up with all these smoothers???
	//
	// Really this is a problem, the sample rate varies drastically
	// across experimental setups.  Most notably due to the number
	// of MS1 vs. MS2 scans, but also scan rate of the machine, and
	// the elution gradient.  I've settled on smoothTHRESHOLD for now.
	//
	// If the MS1 scan rate is too low, using this type of 2D
	// algorithm will not work well.
	//
	public static Smooth2D smoothTHRESHOLD = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			return spectrum;
			}

		public float[] SmoothElution(float[] elution)
			{
			return _thresholdElution(elution);
			}
		};
	Smooth2D smoothALOT = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			//Spectrum.SmoothALittle(spectrum);
			return spectrum;
			}

		protected float[] SmoothElution(float[] elution)
			{
			return Spectrum.FFTsmooth(elution, 12, false);
			}
		};
	Smooth2D smoothMEDIUM = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			//Spectrum.SmoothALittle(spectrum);
			return spectrum;
			}

		protected float[] SmoothElution(float[] elution)
			{
			return Spectrum.FFTsmooth(elution, 6, false);
			}
		};
	Smooth2D smoothD4 = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			//Spectrum.SmoothALittle(spectrum);
			return spectrum;
			}

		protected float[] SmoothElution(float[] elution)
			{
			return Spectrum.WaveletD4(elution);
			}
		};
	Smooth2D smoothALITTLE = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			//Spectrum.SmoothALittle(spectrum);
			return spectrum;
			}

		protected float[] SmoothElution(float[] elution)
			{
			Spectrum.SmoothALittle(elution);
			return elution;
			}
		};
	Smooth2D smoothNONE = new Smooth2D()
		{
		protected float[] SmoothSpectra(float[] spectrum)
			{
			return spectrum;
			}

		protected float[] SmoothElution(float[] elution)
			{
			return elution;
			}
		};


        /**
         * Resample the spectra, within the M/Z range _mzRange,
         * onto a regular grid, with frequency RESAMPLE_FREQ
         * @param scans
         * @param currentThread
         * @param useMedianSmooth
         * @return
         * @throws InterruptedException
         */
        protected float[][] resampleSpectra(Scan[] scans,
                                            Thread currentThread,
                                            boolean useMedianSmooth)
                throws InterruptedException
        {
            float[][] resampledSpectra = new float[scans.length][];
            for (int i = 0; i < scans.length; i++)
            {
                float[][] raw = scans[i].getSpectrum();
                if (currentThread.isInterrupted())
                    throw new InterruptedException();
                resampledSpectra[i] =
                        Spectrum.Resample(raw, _mzRange, SpectrumResampler.getResampleFrequency());
            }
            int height = resampledSpectra[0].length;
            {
                float[] row = null, s = null;
                for (int i = 0; i < height; i++)
                {
                    row = Spectrum.getRow(resampledSpectra, i, row);
                    if (useMedianSmooth) // this will remove lockspray
                    { // use median
                        s = Spectrum.MedianSmooth(row, row.length, s);
                        Spectrum.SmoothALittle(s);
                        Spectrum.setRow(resampledSpectra, i, s);
                    }
                    else
                    {
                        Spectrum.SmoothALittle(row);
                        Spectrum.setRow(resampledSpectra, i, row);
                    }
                }
            }
            return resampledSpectra;
        }

        /**
         * Calculate median intensity at each point on the grid
         *
         * @param scans
         * @param spectra
         * @param width
         * @param height
         * @return
         */
        protected float[][] calculateMedian(Scan[] scans,
                                             float[][] spectra,
                                             int width, int height)
        {
            float[][] median = new float[spectra.length][];
            for (int i = 0; i < scans.length; i++)
                median[i] = Spectrum.MedianWindow(spectra[i], height, 2 * SpectrumResampler.getResampleFrequency(), false);
            float[] row = null;
            for (int r = 0; r < height; r++)
            {
                row = Spectrum.getRow(spectra, r, row);
                float[] m = Spectrum.MedianWindow(row, width, SpectrumResampler.getResampleFrequency(), false);
                for (int s = 0; s < m.length; s++)
                    median[s][r] = Math.max(median[s][r], m[s]);
            }
            return median;
        }

    /**
	 * THIS IS THE MAIN FEATURE FINDING ROUTINE
     *
     * Structure:
     *   Spectrum.Resample()
     *   Smoothing
     *   Separate background from signal (Spectrum.RemoveBackground)
     *   Extract peaks -- wavelet decomposition, smooth, ExtractMaxima2D.analyze()
     *   Sort by decreasing intensity
     *   Throw out peaks with < 5 scans
     *   Break up what look like double-humped elution profiles by breaking peaks
     *   at valleys in the profile
     *   For each peak:
     *           -Compute the extents
     *                  -Walking down scans, make sure that the peak's intensity
     *                   is greater than the intensities of the two m/z's around it
     *                      -Make sure that intensity is above threshold
     *                      -If you start going back up in intensity, stop
     *                  -Same for end scan
     *                  -If too short (<5 scans), toss it out
     *            -"Integrate" to determine total intensity: for each scan, add a
     *             block of intensity corresponding to the intensity on that scan
     *             multiplied by half the amount of time between the adjoining scans
     *            -Create a feature to represent the individual peak
     *   Exclude all features initially
     *   Filter features: for each feature:
     *           find all neighbors within 5 scans and 1.1 m/z units
     *           if the m/z's are too close (within 1.5 / ExtractMaxima2D.getResampleFrequency(), for some reason), no dice
     *           make sure they're lined up, scan-wise
     *           If not enough scans combined in the two features, no dice
     *           If intensities not high enough, no dice
     *           If we got here, unexclude both
     *   FeatureStrategyUsingWindow.ExtractPeptideFeatures(), to tie features together
     *   Dump windows around feature, if asked for
     *   Change scan numbers, which are currently indexes, to the actual scan numbers
     *   If centroided, call AccurateMassCentroid() to fix mass
     *
	 */
	protected Collection<Feature> analyze2D(Scan[] scans) throws InterruptedException
		{
		boolean useMedianSmooth = false;

		Thread currentThread = Thread.currentThread();

		// debug counters
		int d_rawPeaks = 0;
		int d_tooShort = 0;
		int d_filtered = 0;

		_logDebug("analyze2D " + scans[0].getNum() + "-" + scans[scans.length - 1].getNum());
		assert timerAnalyze.start();

		//
		// Convert data into 2D matrix
		// we will do all processing on this data until the end and
		// then process back to "scan" space
		//
		assert timerResample.start();
		float[][] spectra =
                resampleSpectra(scans, currentThread, useMedianSmooth);
		assert timerResample.stop();

        int width = spectra.length;
        int height = spectra[0].length;
        _logDebug("analyze2D datasize = " + (width * height * 4));

		//
		// separate background and signal components
		// we're pretty conservative about what we call background,
		// we can use background and/or median later to filter peaks
		// further.
		//

		timerBackground.start();
        float[][] background = Spectrum.RemoveBackground(spectra);
        float[][] median = calculateMedian(scans, spectra, width, height);
		timerBackground.stop();

		if (currentThread.isInterrupted())
			throw new InterruptedException();

		//
		// the D3 matrix is used to locate peaks,
		// the D3 processing effectively sharpens the raw data, this does a great job of cleaning
		// up messy data, e.g. wide/overlapping peaks, and stray peaks in centroided data
		//
		// why D3?  This is based on the *default* resampling frequency of 36/Da and our max expected charge
		// which I'm presuming to be 6.  Level three represents a feature width of 2^3=8, which
		// is near to 36/6.  A big change in resampling or desired max charge, might require a
        // change.
        //
        //    todo: now that resampling frequency is parameterized, see if this needs to be optimized
        //

		timerExtractPeaks.start();
		float[][] waveletsD3 = new float[width][];
		Pair<float[][], float[][]> tmp = new Pair<float[][], float[][]>(null, null);
		for (int s = 0; s < spectra.length; s++)
			waveletsD3[s] = Spectrum.WaveletD3(spectra[s], tmp);

		//
		// extract all maxima from a time smoothed version of the wavelet transformed data
		//
		// we have a lot of versions of the data now, but it's helpful to keep it all around
		//

        float[][] smooth = new float[waveletsD3.length][];
		for (int s = 0; s < spectra.length; s++)
			smooth[s] = waveletsD3[s].clone();
		smoothTHRESHOLD.smooth(smooth);
		Spectrum.Peak[] maxD3 =
                ExtractMaxima2D.analyze(smooth, _mzRange.min, SpectrumResampler.getResampleInterval(), null, 0.0F);
		d_rawPeaks = maxD3.length;

		// ???? if no raw peaks to process, bail with empty collection to avoid out
        // of bounds errors
		if ( 0 == d_rawPeaks )
			return new ArrayList<Feature>();

		Arrays.sort(maxD3, Spectrum.comparePeakIntensityDesc);
		if (true)
			{
			_logDebug("raw peaks: " + maxD3.length);
			if (maxD3.length > 0)
				for (int i = 0; i < 10; i++)
					_logDebug("  " + (i) + ": " +
                              maxD3[maxD3.length * i / 100].getIntensity());
			if (maxD3.length > 0)
				for (int i = 1; i < 10; i++)
					_logDebug("  " + (i * 10) + ": " +
                              maxD3[maxD3.length * i / 10].getIntensity());
			_logDebug(" 100: " + maxD3[maxD3.length - 1].getIntensity());
			}


		//
		// We probably have a LOT of peaks at this point.  95% are extremely small.
		// The size<shortPeak check in this loop will throw out most, and I
        // don't need to be too aggressive here
		//
		// NOTE: this loop also tries to handle "double humped elution profiles"
		// by breaking peaks at "valleys" in the elution profile.  This could perhaps
		// be improved upon with a second pass that specifically handles this case and
		// makes sure that in a region the splitting is consistent.
		//


        //TODO: parameterize.  Careful, because this is messed with below
        int shortPeak = 5;
		ArrayList<Feature> waveletPeaks = new ArrayList<Feature>();
		for (int i = 0, max = maxD3.length; i < max; i++)
			{
			Spectrum.Peak peak = maxD3[i];
			int imz = Math.round((peak.mz - _mzRange.min) * SpectrumResampler.getResampleFrequency());
			if (imz < 2 || imz >= height - 2)
				continue;
			int scan = peak.scan;
			float peakIn = spectra[peak.scan][imz];
			float peakBg = background[peak.scan][imz];
			float peakMd = median[peak.scan][imz];

			// compute extent of this feature
			// this is complicated by the fact that we may have
            // double-humped peaks we need to separate
			int scanStart = scan;
			int scanEnd = scan;
			float thresholdIn = 1 + 0.5F * background[scan][imz] +
                                2 * median[scan][imz];
			thresholdIn = Math.max(peakIn / 20, thresholdIn);
			float minIn = Float.MAX_VALUE;
			int scanMinIn = scan;
			for (; scanStart > 0; scanStart--)
				{
				float spectraIn = spectra[scanStart][imz];
				// are we still on the ridge and > 5% of max
				if (spectraIn < thresholdIn ||
                        waveletsD3[scanStart][imz] < waveletsD3[scanStart][imz - 2] ||
                        waveletsD3[scanStart][imz] < waveletsD3[scanStart][imz + 2])
					{
					scanStart++;
					break;
					}
				// have we crossed a 'valley'
				float smoothIn = smooth[scanStart][imz];
				if (smoothIn < minIn)
					{
					minIn = spectraIn;
					scanMinIn = scanStart;
					}
				else if (smoothIn > minIn * 1.1 + thresholdIn) // made up cutoff
					{
					scanStart = scanMinIn;
					break;
					}
				}
			minIn = Float.MAX_VALUE;
			scanMinIn = scan;
			for (; scanEnd < width - 1; scanEnd++)
				{
				float spectraIn = spectra[scanEnd][imz];
				// are we still on the ridge and > 5% of max
				if (spectraIn < thresholdIn ||
                        waveletsD3[scanEnd][imz] < waveletsD3[scanEnd][imz - 2] ||
                        waveletsD3[scanEnd][imz] < waveletsD3[scanEnd][imz + 2])
					{
					scanEnd--;
					break;
					}
				// have we crossed a 'valley'
				float smoothIn = smooth[scanEnd][imz];
				if (smoothIn < minIn)
					{
					minIn = spectraIn;
					scanMinIn = scanEnd;
					}
				else if (smoothIn > minIn * 1.1 + thresholdIn) // made up cutoff
					{
					scanEnd = scanMinIn;
					break;
					}
				}

			if (scanEnd - scanStart + 1 < shortPeak)
				{
				d_tooShort++;
				continue;
				}

			// total intensity calculation
			// "integrate" intensity over time
			float totalIn = 0;
			for (int s = scanStart; s <= scanEnd; s++)
				{
				double in = spectra[s][imz];
				double time = 1.0;
				if (s > 0 && s + 1 < scans.length)
					time = (scans[s + 1].getDoubleRetentionTime() -
                            scans[s - 1].getDoubleRetentionTime()) /
                            2;
				totalIn += time * in;
				}

			Feature f = new Feature(peak.scan, scanStart, scanEnd,
                                    peak.mz, peakIn, 0, -1, 0);
			f.background = peakBg;
			f.median = peakMd;
			f.totalIntensity = totalIn;
			waveletPeaks.add(f);
			}

		if (currentThread.isInterrupted())
			throw new InterruptedException();

		Collections.sort(waveletPeaks, Spectrum.compareFeatureLengthDesc);

		// some stats
		if (_log.isDebugEnabled())
			{
			_logDebug("FEATURE LENGTHS " + waveletPeaks.size() + " peaks");
			if (waveletPeaks.size() > 0)
				{
				for (int i = 0; i < 10; i++)
					_logDebug("  " + (i) + ": " +
                            maxD3[waveletPeaks.size() * i / 100].getIntensity());
				}
			if (waveletPeaks.size() > 0)
				{
				for (int i = 1; i < 10; i++)
					_logDebug("  " + (i * 10) + ": " +
                            maxD3[waveletPeaks.size() * i / 10].getIntensity());
				}
			}

		// OK, up the bar for a feature (min length of at largest peak in the feature)
		shortPeak++;
		_logDebug("SHORT PEAK " + shortPeak);

		//
		// Filtering
		//
		// Sticking with the idea of persistance, eliminate peaks
		// that don't seem to have correllated peaks
		//

		Tree2D tree = new Tree2D();
		for (Feature f : waveletPeaks)
			{
			f.excluded = true;
			tree.add(f.scan, f.mz, f);
			}
		Collections.sort(waveletPeaks, Spectrum.compareFeatureLengthDesc);
//QUESTION: so only features with other related features get included... so how do we
//end up with any single-peak features?
        for (Feature f : waveletPeaks)
			{
			int lenF = f.getScanLast() - f.getScanFirst() + 1;
			ArrayList<Feature> l = (ArrayList<Feature>)
                    tree.getPoints(f.scan - 5, f.mz - 1.1F, f.scan + 5, f.mz + 1.1F);
			for (Feature neighbor : l)
				{
				int lenN = neighbor.getScanLast() - neighbor.getScanFirst() + 1;
				// we want features of different masses
//QUESTION: what is the significance of 1.5/ExtractMaxima2D.getResampleFrequency()?
                if (Math.abs(neighbor.mz - f.mz) < 1.5 / SpectrumResampler.getResampleFrequency())
					continue;
				// are centers aligned
				int s = Math.max(f.getScanFirst(), neighbor.getScanFirst());
				int e = Math.min(f.getScanLast(), neighbor.getScanLast());
				if (f.scan < s || f.scan > e || neighbor.scan < s || neighbor.scan > e)
					continue;
				// if neither is chosen yet, add a small hurdle
				if (f.excluded && neighbor.excluded)
					{
//QUESTION: why checking the total number of scans in both?
                    if (lenF + lenN < shortPeak)
						continue;
					if (f.intensity < 1 + 0.5 * f.background + 3 * Math.max(1, f.median) &&
					        neighbor.intensity < 1 + 0.5 * f.background +
                                                 3 * Math.max(1, neighbor.median))
						continue;
					}
//QUESTION: why aren't we tying these together somehow at this point?
                f.excluded = false;
				neighbor.excluded = false;
				}
			}
		ArrayList<Feature> copy = new ArrayList<Feature>();
		for (Feature f : waveletPeaks)
			{
			if (!f.excluded)
				copy.add(f);
			else
				d_filtered++;
			}
		waveletPeaks = copy;
		assert timerExtractPeaks.stop();

		if (currentThread.isInterrupted())
			throw new InterruptedException();

		if (_log.isDebugEnabled())
			{
			_logDebug("kept " + waveletPeaks.size() + " peaks after filtering");
			Collections.sort(waveletPeaks, Spectrum.comparePeakIntensityDesc);
			for (int i = 0; i < waveletPeaks.size() && i < 100; i++)
				_logDebug(waveletPeaks.get(i).toString());
			}


		if (debugReturnPeaks) // see intermediate results
			{
            for (Feature peak : waveletPeaks)
				{
				// fix up scan
				Scan scan = scans[peak.scan];
				peak.setTime((float)scan.getDoubleRetentionTime());
				peak.scan = scan.getNum();
				peak.scanFirst = scans[peak.scanFirst].getNum();
				peak.scanLast = scans[peak.scanLast].getNum();
				peak.setPeaks(1); // otherwise these will get filtered by default
				}
			return waveletPeaks;
			}


		Feature[] peaks = waveletPeaks.toArray(new Feature[0]);
		Arrays.sort(peaks, Spectrum.comparePeakScanAsc);// debug only
		Arrays.sort(peaks, Spectrum.comparePeakMzAsc);// debug only
		_logDebug(Feature.getFeatureHeader(null));

//too verbose
//		for (Feature peak : peaks)
//			{
//			_logDebug(peak.toString(null));
//			}


		//
		// ExtractPeptideFeatures
		//
		// combine peaks into features representing peptides
		//

		assert timerExtractPeptides.start();
		ArrayList<Feature> peptidesAll = ExtractPeptideFeatures(_run, peaks);

		//
		// Extract a window of intensities around each Feature
		//
		// this is may be used by downstream analsysis programs
		//

		int dumpWindowSize = getDumpWindowSize();
		if ( dumpWindowSize > 0 )
			{
			int nSamples = dumpWindowSize * SpectrumResampler.getResampleFrequency();
			for (Feature f : peptidesAll)
				{
				int mzIndex = (int)( ( f.mz - _mzRange.min) * SpectrumResampler.getResampleFrequency() );
				int scanIndex = f.scan;
				f.intensityLeadingPeaks = dumpWindowSize;
				f.intensityTrailingPeaks = dumpWindowSize;
				f.intensityWindow = new float[2 * nSamples];
				int k = 0;
				for (int j = mzIndex - nSamples; j < mzIndex + nSamples; j++, k++)
					{
					if ( j < 0 || j >= spectra[scanIndex].length )
						f.intensityWindow[k] = 0.f; // pad with zeros if we hit an edge
					else
						f.intensityWindow[k] = spectra[scanIndex][j];
					}
				}
			}

		//
		// fix up scanNum
		//

		for (Feature f : peptidesAll)
			{
			Scan scan = scans[f.scan];
			f.setTime((float)scan.getDoubleRetentionTime());
			f.scan = scan.getNum();
			f.setScanFirst(scans[f.getScanFirst()].getNum());
			f.setScanLast(scans[f.getScanLast()].getNum());
			}

	

		// if data are centroided, or if requested explicitly for profile
		// mode data, attempt to get accurate masses
		if (getAccurateMassAdjustmentScans() > 0 || _run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
			{
			Feature[] a = peptidesAll.toArray(new Feature[0]);
			Arrays.sort(a, Spectrum.comparePeakScanAsc);
			for (int i = 0; i < a.length; i++)
				{
//				float mz = AccurateMassCentroid(_run, a[i]);
				float mz = AccurateMass(_run, a[i], getAccurateMassAdjustmentScans());
				if (mz > 0)
					{
					a[i].setMz(mz);
					a[i].setAccurateMZ(true);
					a[i].updateMass();
					}
				if (currentThread.isInterrupted())
					throw new InterruptedException();
				}
			}

		assert timerExtractPeptides.stop();
		assert timerAnalyze.stop();
		CPUTimer.dumpAllTimers();
		return peptidesAll;
		}


	Collection PeaksAsFeatures(Spectrum.Peak[] peaks)
		{
		ArrayList l = new ArrayList();
		for (int p = 0; p < peaks.length; p++)
			l.add(new Feature(peaks[p]));
		return l;
		}


	public static float AccurateMass(MSRun run, Feature f, int nScans)
		{
		if (run.getHeaderInfo().getDataProcessing().getCentroided() == 1)
			return AccurateMassCentroid(run, f);
		else
			return AccurateMassProfile(run, f, nScans);
		}


	public static float AccurateMassCentroid(MSRun run, Feature f)
		{
		// NOTE about isotopes
		// some heavy isotopes are greater than +1Da (C=1.0033), and some are less (N=0.9971)
		// I use 1.0013 as working number (C13 is more common than N15)
		final double ISOTOPE_FACTOR = 1.0013;

		// CONSIDER: does it make sense to do any sort of averaging across a few scans?
		Scan scan = run.getScan(run.getIndexForScanNum(f.scan));
		float[][] s = scan.getSpectrum();
		double delta = .66 / SpectrumResampler.getResampleFrequency();

		double mzP0 = 0;
		double sumMZ = 0;
		double sumIN = 0;
		for (int i = 0; i < 2 && i < f.comprised.length; i++)
			{
			double mzPi;
			if (f.comprised[i] == null)
				mzPi = f.mz + i * ISOTOPE_FACTOR / f.charge;
			else
				mzPi = f.comprised[i].mz;

			// find biggest close peak
			// could use more complicated "distance" measure
			// but in centroided data there is usually only one 'real' peak
			// in a space this small with maybe some tiny side peaks
			int p = Arrays.binarySearch(s[0], (float)(mzPi - delta));
			if (p < 0) p = -1 * (p + 1);
			double mzBiggest = 0;
			double inBiggest = 0;
			for (; p < s[0].length; p++)
				{
				float mz = s[0][p];
				if (mz > mzPi + delta)
					break;
				float in = s[1][p];
				if (in > inBiggest)
					{
					mzBiggest = mz;
					inBiggest = in;
					}
				}
			if (mzBiggest == 0) // UNDONE: use wider window???
				{
				//System.out.println("missed " + i + "\t" + f.toString());
				return 0;
				}
			if (f.charge == 0)
				return (float)mzBiggest;
			if (i == 0)
				mzP0 = mzBiggest;
			// weighted average of peaks
			sumMZ += inBiggest * (mzBiggest - i * ISOTOPE_FACTOR / f.charge);
			sumIN += inBiggest;
			}
		double avgMZ = sumMZ / sumIN;
		// if we got this all right we'd expect ABS(newMZ-mzP0) to be less
        // than the resolution of the machine
		// I'm going to assume that the centroided means FT (very accurate)
		if (Math.abs(avgMZ - mzP0) > mzP0 * 5.0 / 1000000.0)
			return 0;

		// NOTE: could return avgMZ here, however
		//  a) it's very close to mzP0 (see test)
		//  b) I suspect people would rather see the value of peak from the source data
		return (float)mzP0;
		}

            /**
             * adjustment of profile-mode mass; average results from some number of adjacent scans
             */
            public static float AccurateMassProfile(MSRun run, Feature f, int nScans)
            {
                if (nScans <= 0)
                    return 0.f;

                int scanIndex = run.getIndexForScanNum(f.scan);

                int lowScanIndex  = (int) Math.max(scanIndex - (nScans - 1)/2.0 + .5, 0);
                int highScanIndex = (int) Math.min(scanIndex + (nScans - 1)/2.0 + .5, run.getScanCount() - 1);

                float sumMz = 0.f;
                int n = 0;

                for (int s = lowScanIndex; s <= highScanIndex; s++)
                {
                    Scan scan = run.getScan(s);
                    float maxMz = AccurateMassProfileCenter(scan, f);

                    if (maxMz > 0.f)
                    {
                        sumMz += maxMz;
                        n++;
                    }
                }
                return n > 0 ? sumMz / n : 0.f;
            }

            /**
             * adjustment of profile-mode mass using max peak
             */
            private static float AccurateMassProfileMax(Scan scan, Feature f)
            {
                float[][] s = scan.getSpectrum();
                double delta = .5 / SpectrumResampler.getResampleFrequency();
                
                double lowMz = f.mz - delta;
                double highMz = f.mz + delta;

                int p = Arrays.binarySearch(s[0], (float) lowMz);
                if (p < 0)
                    p = -1 * (p + 1);

                double maxMz = 0;
                double maxInt = 0;

                for (; p < s[0].length; p++)
                {
                    if (s[0][p] > highMz)
                        break;

                    if (s[1][p] > maxInt)
                    {
                        maxMz = s[0][p];
                        maxInt = s[1][p];
                    }
                }

                return (float) maxMz;
            }

            /**
             * adjustment of profile-mode mass using center of mass
             */
            private static float AccurateMassProfileCenter(Scan scan, Feature f)
            {
                float[][] s = scan.getSpectrum();
                double delta = .66666667 / SpectrumResampler.getResampleFrequency();
                
                double lowMz = f.mz - delta;
                double highMz = f.mz + delta;

                int p = Arrays.binarySearch(s[0], (float) lowMz);
                if (p < 0)
                    p = -1 * (p + 1);

                double sumMz = 0;
                double sumInt = 0;

                for (; p < s[0].length; p++)
                {
                    if (s[0][p] > highMz)
                        break;

                    sumMz += s[0][p] * s[1][p]; // mz weighted by intensity
                    sumInt += s[1][p];
                }

                // Couldn't figure out a decent match
                if (sumInt <= 0.0)
                    return 0.f;

                return (float) (sumMz/sumInt);
            }

	protected static void _logDebug(String s)
		{
		_log.debug(s);
		}
	}
