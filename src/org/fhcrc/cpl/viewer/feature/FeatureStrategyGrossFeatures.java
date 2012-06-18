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
import org.fhcrc.cpl.viewer.feature.ExtractEdgeFeatures;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.Scan;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Arrays;

/**
 * User: mbellew
 * Date: Nov 2, 2004
 */
public class FeatureStrategyGrossFeatures extends FeatureStrategyUsingWindow //extends FeatureExtractor
	{
	int _startNum = 0;
	int _endNum = 0;

	public FeatureStrategyGrossFeatures(MSRun run, int scanIndex, int count, int maxCharge, FloatRange range, double sn)
		{
		super(run, scanIndex, count, maxCharge, range, sn);

		int c2 = Math.max(256, count);
		scanIndex = Math.max(0, scanIndex - (c2 - count) / 2);
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
	            {
	            feature.setPeaks(1); // will get filtered by default if ==0
                filtered.add(feature);
	            }
            }
        return (Feature[]) filtered.toArray(new Feature[filtered.size()]);
		}


	protected Collection analyze2D(Scan[] scans)
		{
		//
		// Calculate maxima in the region of the scan to analyze
		//
		Spectrum.Peak[][] t = ExtractEdgeFeatures.analyze(scans, _mzRange, 0);
		Spectrum.Peak[] features = t[2];
		return Arrays.asList(features);

		/*
		ArrayList l = new ArrayList();
		for (int f=0 ; f<features.length ; f++)
			{
			Spectrum.Peak peak = features[f];
			l.add(new Spectrum.Feature(peak.scan, peak.mz, peak.intensity, 0, 0));
			}
        return l; */
		}


	public int getType()
		{
		return TYPE_2D;
		}
	}
