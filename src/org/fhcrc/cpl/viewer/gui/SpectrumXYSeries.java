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

import org.apache.log4j.Logger;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.general.SeriesException;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: Oct 11, 2004
 * Time: 11:27:56 AM
 */
public class SpectrumXYSeries extends XYSeries
	{
	private static Logger _log = Logger.getLogger(XYSeries.class);

	public SpectrumXYSeries(String name, float[][] spectrum, FloatRange r)
		{
		this(name, spectrum, Math.min(spectrum[0].length, spectrum[1].length), r);
		}


	public SpectrumXYSeries(String name, float[][] spectrum, int len, FloatRange r)
		{
		super(name, false,
              true); //allow duplicates for speed.  We're guaranteed no duplicates by the nature of the data

		try
			{
			assert len <= spectrum[0].length;
			if (len == -1)
				len = spectrum[0].length;
			len = Math.min(len, spectrum[0].length);

			float[] X = spectrum[0];
			float[] Y = spectrum[1];
			for (int i = 0 ; i<len ; i++)
				{
				float x = X[i], y = Y[i];
				if (y == 0)
					continue;
				if (r != null && !r.contains(x))
					continue;
				this.add(x, y, false);
				}
			}
		catch (SeriesException x)
			{
			throw x;
			}
		}
	}
