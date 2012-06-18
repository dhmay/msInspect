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
package org.fhcrc.cpl.viewer.util;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: Oct 5, 2004
 * Time: 5:19:58 PM
 */
public class Haar
	{
	public static float[] transform(float[] signal, int l)
		{
		return transform(signal, l, null);
		}

	public static float[] transform(float[] signal, int l, float[] x)
		{
		x = (x==null) ? new float[signal.length] : x;

		// init sum
		double sum = -1 * l * signal[0];
		for (int i=0 ; i<l && i<signal.length ; i++)
			sum += signal[i];

		for (int i = 0; i<signal.length ; i++)
			{
			sum += (i-l >= 0) ? signal[i-l] : signal[0];
			sum -= 2 * signal[i];
			sum += (i+l < signal.length) ? signal[i+l] : signal[signal.length-1];
            x[i] = (float)(sum / l);
			}
		return x;
		}
	}
