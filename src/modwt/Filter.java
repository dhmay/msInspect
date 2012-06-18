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

package modwt;


/**
 * This structure is a convenient way to keep the pieces of a wavelet
 * Filter in one place.
 * <p/>
 * L = length of Filter
 * h = wavelet Filter (dvector of length L)
 * g = scaling Filter (dvector of length L)
 */
public class Filter
	{
	final double[] hhaar = new double[]
		{0.7071067811865475, -0.7071067811865475};
	final double[] ghaar = new double[]
		{0.7071067811865475, 0.7071067811865475};
	final double[] hd4 = new double[]
		{-0.1294095225512603, -0.2241438680420134,
		 0.8365163037378077, -0.4829629131445341};
	final double[] gd4 = new double[]
		{0.4829629131445341, 0.8365163037378077,
		 0.2241438680420134, -0.1294095225512603};
	final double[] hd6 = new double[]
		{0.0352262918857096, 0.0854412738820267,
		 -0.1350110200102546, -0.4598775021184915,
		 0.8068915093110928, -0.3326705529500827};
	final double[] gd6 = new double[]
		{0.3326705529500827, 0.8068915093110928,
		 0.4598775021184915, -0.1350110200102546,
		 -0.0854412738820267, 0.0352262918857096};
	final double[] hd8 = new double[]
		{-0.0105974017850021, -0.0328830116666778,
		 0.0308413818353661, 0.1870348117179132,
		 -0.0279837694166834, -0.6308807679358788,
		 0.7148465705484058, -0.2303778133074431};
	final double[] gd8 = new double[]
		{0.2303778133074431, 0.7148465705484058,
		 0.6308807679358788, -0.0279837694166834,
		 -0.1870348117179132, 0.0308413818353661,
		 0.0328830116666778, -0.0105974017850021};
	final double[] hla8 = new double[]
		{0.03222310060407815, 0.01260396726226383,
		 -0.09921954357695636, -0.29785779560560505,
		 0.80373875180538600, -0.49761866763256290,
		 -0.02963552764596039, 0.07576571478935668};
	final double[] gla8 = new double[]
		{-0.07576571478935668, -0.02963552764596039,
		 0.49761866763256290, 0.80373875180538600,
		 0.29785779560560505, -0.09921954357695636,
		 -0.01260396726226383, 0.03222310060407815, };
	final double[] hla16 = new double[]
		{0.0018899503329007, 0.0003029205145516,
		 -0.0149522583367926, -0.0038087520140601,
		 0.0491371796734768, 0.0272190299168137,
		 -0.0519458381078751, -0.3644418948359564,
		 0.7771857516997478, -0.4813596512592012,
		 -0.0612733590679088, 0.1432942383510542,
		 0.0076074873252848, -0.0316950878103452,
		 -0.0005421323316355, 0.0033824159513594};
	final double[] gla16 = new double[]
		{-0.0033824159513594, -0.0005421323316355,
		 0.0316950878103452, 0.0076074873252848,
		 -0.1432942383510542, -0.0612733590679088,
		 0.4813596512592012, 0.7771857516997478,
		 0.3644418948359564, -0.0519458381078751,
		 -0.0272190299168137, 0.0491371796734768,
		 0.0038087520140601, -0.0149522583367926,
		 -0.0003029205145516, 0.0018899503329007};


	int L;
	double[] h;
	double[] g;

	protected Filter()
		{
		}

	/**
	 * Determines which wavelet Filter, and corresponding scaling Filter,
	 * based on the character string `choice.'
	 * <p/>
	 * Input:
	 * choice = character string (e.g., "haar", "d4, "d6", "d8", "la8", "la16")
	 * <p/>
	 * Output:
	 * fills in h, g, and L for a wavelet Filter structure
	 */
	public Filter(String choice)
		{
		if ("haar".equals(choice))
			{
			L = 2;
			h = hhaar;
			g = ghaar;
			}
		else if ("d4".equals(choice))
			{
			L = 4;
			h = hd4;
			g = gd4;
			}
		else if ("d6".equals(choice))
			{
			L = 6;
			h = hd6;
			g = gd6;
			}
		else if ("d8".equals(choice))
			{
			L = 8;
			h = hd8;
			g = gd8;
			}
		else if ("la8".equals(choice))
			{
			L = 8;
			h = hla8;
			g = gla8;
			}
		else if ("la16".equals(choice))
			{
			L = 16;
			h = hla16;
			g = gla16;
			}
		else
			throw new IllegalArgumentException("...unimplemented wavelet choice...");
		}


	public String toString()
		{
		StringBuffer sb = new StringBuffer();
		sb.append("Wavelet Filter (L = ").append(L).append("):\n");
		sb.append("  h := \n");
		printdvec(sb, h);
		sb.append("  g := \n");
		printdvec(sb, g);
		return sb.toString();
		}


	private static void printdvec(StringBuffer sb, double[] v)
		{
		for (int i = 0; i < v.length; i++)
			{
			sb.append(v[i]).append(' ');
			}
		sb.append('\n');
		}


	/**
	 * This function converts the basic Haar wavelet Filter
	 * (h0 = -0.7017, h1 = 0.7017) into a Filter of length `translate.'
	 * The Filter is re-normalized so that it's sum of squares equals 1.
	 * <p/>
	 * Input:
	 * f     = haar wavelet Filter structure
	 * translate = integer
	 * <p/>
	 * Output:
	 * fills in h, g, and L for a wavelet Filter structure
	 */
	public Filter convert_haar(int scale)
		{
		Filter out = new Filter();
		double inv_sqrt_scale = 1.0 / Math.sqrt(scale);

		out.h = new double[L * scale + 1];
		out.g = new double[L * scale + 1];

		for (int l = 0; l < scale; l++)
			{
			out.h[l] = h[0] * inv_sqrt_scale;
			out.g[l] = g[0] * inv_sqrt_scale;
			out.h[l + scale] = h[1] * inv_sqrt_scale;
			out.g[l + scale] = g[1] * inv_sqrt_scale;
			}
		return (out);
		}
	}
