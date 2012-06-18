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

import org.jfree.data.xy.AbstractIntervalXYDataset;

/**
 * User: mbellew
 * Date: Jun 29, 2004
 * Time: 11:44:22 AM
 */
public class SpectrumDataSet extends AbstractIntervalXYDataset
	{
	float[][] peaks;
	int length;
	Double startXValue;
	Double endXValue;
	Double startYValue;
	Double endYValue;


	public SpectrumDataSet(float[][] p, int len)
		{
		peaks = p;
		length = len;
		}


	public SpectrumDataSet(float[][] p)
		{
		this(p, p[0].length);
		}


	public double getStartXValue(int series, int item)
		{
		return (double)peaks[0][item];
		}


	public double getEndXValue(int series, int item)
		{
		return (double)peaks[0][item];
		}


	public double getStartYValue(int series, int item)
		{
		return (double)peaks[1][item];
		}


	public double getEndYValue(int series, int item)
		{
		return (double)peaks[1][item];
		}


	public int getItemCount(int series)
		{
		return length;
		}


	public double getXValue(int series, int item)
		{
		return getEndXValue(series, item);
		}


	public double getYValue(int series, int item)
		{
		return getEndYValue(series, item);
		}


	public int getSeriesCount()
		{
		return 1;
		}

	/**
	 * Deprecated as of JFreeChart 1.0.0?
	public String getSeriesName(int series)
		{
		return "Peaks";
		}
	 */

	public Comparable getSeriesKey(int series)
		{
		return "Peaks";
		}

	public Number getX(int s, int i)
		{
		return new Double(getXValue(s,i));
		}

	public Number getY(int s, int i)
		{
		return new Double(getYValue(s,i));
		}

	public Number getStartX(int s, int i)
		{
		return new Double(getStartXValue(s,i));
		}

	public Number getEndX(int s, int i)
		{
		return new Double(getEndXValue(s,i));
		}

	public Number getStartY(int s, int i)
		{
		return new Double(getStartYValue(s,i));
		}

	public Number getEndY(int s, int i)
		{
		return new Double(getEndYValue(s,i));
		}
	}
