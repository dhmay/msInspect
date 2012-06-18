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

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.Arrays;
import java.io.PrintStream;



public class Transform
	{
	final static double inv_sqrt_2 = 1.0 / Math.sqrt(2.0);

	private static int wrap(int i, int length)
		{
		if (i >= length) return i - length;
		if (i < 0) return i + length;
		return i;
		}


	/**
	 * Compute the discrete wavelet Transform (DWT).  This method uses the
	 * pyramid algorithm and was adapted from pseudo-code written by
	 * D. B. Percival.  Periodic boundary conditions are assumed.
	 *
	 * @param Vin  vector of wavelet smooths (data if first iteration)
	 * @param f    wavelet Filter structure (e.g., Haar, D(4), LA(8), ...)
	 * @param Wout OUT vector of wavelet coefficients
	 * @param Vout OUT vector of wavelet smooths
	 */
	public static void dwt(double[] Vin, int M, Filter f, double[] Wout, double[] Vout)
		{
		for (int t = 0; t < M / 2; t++)
			{
			Wout[t] = Vout[t] = 0.0;
			for (int l = 0, k = 2 * t + 1; l < f.L; l++)
				{
				Wout[t] += f.h[l] * Vin[k];
				Vout[t] += f.g[l] * Vin[k];
				k = k == 0 ? M - 1 : k - 1;
				}
			}
		}


	public static void dwt(float[] Vin, int M, Filter f, float[] Wout, float[] Vout)
		{
		for (int t = 0; t < M / 2; t++)
			{
			Wout[t] = Vout[t] = 0.0F;
			for (int l = 0, k = 2 * t + 1; l < f.L; l++)
				{
				Wout[t] += f.h[l] * Vin[k];
				Vout[t] += f.g[l] * Vin[k];
				k = k == 0 ? M - 1 : k - 1;
				}
			}
		}


	/**
	 * Compute the inverse DWT via the pyramid algorithm.  This code was
	 * adapted from pseudo-code written by D. B. Percival.  Periodic
	 * boundary conditions are assumed.
	 *
	 * @param Win vector of wavelet coefficients
	 * @param Vin vector of wavelet smooths
	 * @param M   length of Win, Vin
	 * @param f   wavelet Filter structure (e.g., Haar, D(4), LA(8), ...)
	 * @param Xout OUT vector of reconstructed wavelet smooths (eventually the data)
	 */
	public static void idwt(double[] Win, double[] Vin, int M, Filter f, double[] Xout)
		{
		for (int t = 0, m = 0, n = 1; t < M; t++)
			{
			int u = t;
			int i = 1;
			int j = 0;
			Xout[m] = Xout[n] = 0.0;
			for (int l = 0; l < f.L / 2; l++)
				{
				Xout[m] += f.h[i] * Win[u] + f.g[i] * Vin[u];
				Xout[n] += f.h[j] * Win[u] + f.g[j] * Vin[u];
				u = u + 1 >= M ? 0 : u + 1;
				i += 2;
				j += 2;
				}
			m += 2;
			n += 2;
			}
		}


	public static void idwt(float[] Win, float[] Vin, int M, Filter f, float[] Xout)
		{
		for (int t = 0, m = 0, n = 1; t < M; t++)
			{
			int u = t;
			int i = 1;
			int j = 0;
			Xout[m] = Xout[n] = 0.0F;
			for (int l = 0; l < f.L / 2; l++)
				{
				Xout[m] += f.h[i] * Win[u] + f.g[i] * Vin[u];
				Xout[n] += f.h[j] * Win[u] + f.g[j] * Vin[u];
				u = u + 1 >= M ? 0 : u + 1;
				i += 2;
				j += 2;
				}
			m += 2;
			n += 2;
			}
		}


	/**
	 * Compute the maximal overlap discrete wavelet Transform (MODWT).
	 * This method uses the pyramid algorithm and was adapted from
	 * pseudo-code written by D. B. Percival.
	 *
	 * @param Vin  vector of wavelet smooths (data if k=1)
	 * @param N    length of Vin
	 * @param k    iteration (1, 2, ...)
	 * @param f    wavelet Filter structure (e.g., Haar, D(4), LA(8), ...)
	 * @param Wout OUT vector of wavelet coefficients
	 * @param Vout OUT vector of wavelet smooths
	 */
	public static void modwt(double[] Vin, int N, int k, Filter f, double[] Wout, double[] Vout)
		{
		double[] ht = new double[f.L];
		double[] gt = new double[f.L];
		double inv_sqrt_2 = 1.0 / Math.sqrt(2.0);
		int pow2_k = pow2(k - 1);

		for (int l = 0; l < f.L; l++)
			{
			ht[l] = f.h[l] * inv_sqrt_2;
			gt[l] = f.g[l] * inv_sqrt_2;
			}

		for (int t = 0; t < N; t++)
			{
			int j = t;
			Wout[t] = Vout[t] = 0.0;
			for (int l = 0; l < f.L; l++)
				{
				Wout[t] += ht[l] * Vin[j];
				Vout[t] += gt[l] * Vin[j];
				j = wrap(j - pow2_k, N);
				}
			}
		}


	public static void modwt(float[] Vin, int N, int k, Filter f, float[] Wout, float[] Vout)
		{
		double[] ht = new double[f.L];
		double[] gt = new double[f.L];
		double inv_sqrt_2 = 1.0 / Math.sqrt(2.0);
		int pow2_k = pow2(k - 1);

		for (int l = 0; l < f.L; l++)
			{
			ht[l] = f.h[l] * inv_sqrt_2;
			gt[l] = f.g[l] * inv_sqrt_2;
			}

		for (int t = 0; t < N; t++)
			{
			int j = t;
			Wout[t] = Vout[t] = 0.0F;
			for (int l = 0; l < f.L; l++)
				{
				Wout[t] += ht[l] * Vin[j];
				Vout[t] += gt[l] * Vin[j];
				j = wrap(j - pow2_k, N);
				}
			}
		}



	/**
	 * Compute the inverse MODWT via the pyramid algorithm.  Adapted from
	 * pseudo-code written by D. B. Percival.
	 * @param Win  vector of wavelet coefficients
	 * @param Vin  vector of wavelet smooths
	 * @param N    length of Win, Vin
	 * @param k    detail number
	 * @param f    wavelet Filter structure
	 * @param Vout OUT vector of wavelet smooths
	 */

	public static void imodwt(double[] Win, double[] Vin, int N, int k, Filter f, double[] Vout)
		{
		double[] ht = new double[f.L];
		double[] gt = new double[f.L];
		int pow2_k = pow2(k - 1);

		for (int l = 0; l < f.L; l++)
			{
			ht[l] = f.h[l] * inv_sqrt_2;
			gt[l] = f.g[l] * inv_sqrt_2;
			}

		for (int t = 0; t < N; t++)
			{
			int j = t;
			Vout[t] = 0.0;
			for (int l = 0; l < f.L; l++)
				{
				Vout[t] += (ht[l] * Win[j]) + (gt[l] * Vin[j]); /* GMA */
				j = wrap(j + pow2_k, N);                                        /* GMA */
				}
			}
		}


	public static void imodwt(float[] Win, float[] Vin, int N, int k, Filter f, float[] Vout)
		{
		double[] ht = new double[f.L];
		double[] gt = new double[f.L];
		int pow2_k = pow2(k - 1);

		for (int l = 0; l < f.L; l++)
			{
			ht[l] = f.h[l] * inv_sqrt_2;
			gt[l] = f.g[l] * inv_sqrt_2;
			}

		for (int t = 0; t < N; t++)
			{
			int j = t;
			Vout[t] = 0.0F;
			for (int l = 0; l < f.L; l++)
				{
				Vout[t] += (ht[l] * Win[j]) + (gt[l] * Vin[j]); /* GMA */
				j = wrap(j + pow2_k, N);                                        /* GMA */
				}
			}
		}



	/**
	 * The functions for computing wavelet transforms assume periodic
	 * boundary conditions, regardless of the data's true nature.  By
	 * adding a `backwards' version of the data to the end of the current
	 * data vector, we are essentially reflecting the data.  This allows
	 * the periodic methods to work properly.
	 */

	public static void reflect_vector(double[] Xin, int N, double[] Xout)
		{
		for (int t = 0; t < N; t++)
			Xout[t] = Xin[t];
		for (int t = 0; t < N; t++)
			Xout[N + t] = Xin[N - 1 - t];
		}


	public static void reflect_vector(float[] Xin, int N, float[] Xout)
		{
		for (int t = 0; t < N; t++)
			Xout[t] = Xin[t];
		for (int t = 0; t < N; t++)
			Xout[N + t] = Xin[N - 1 - t];
		}


	/**
	 * Peform the discrete wavelet Transform (DWT) or maximal overlap
	 * discrete wavelet Transform (MODWT) to a time-series and obtain a
	 * specified number (K) of wavelet coefficients and subsequent
	 * wavelet smooth.  Reflection is used as the default boundary
	 * condition.
	 * @param X        time-series (vector of data)
	 * @param K        number of details desired
	 * @param f        wavelet Filter structure
	 * @param method   character string (either "dwt" or "modwt")
	 * @param boundary boundary condition (either "period" or "reflect")
	 * @param Xout     vectors to be used to return results, may be null
	 * @return         matrix of wavelet coefficients and smooth
	 */

	public static double[][] decompose(double[] X, int N, int K, Filter f, String method, String boundary, double[][] Xout)
		{
		int scale, length;
		double[] Vin, Vout;

		if (null == Xout)
			Xout = new double[K+1][];

		if (!"dwt".equals(method) && !"modwt".equals(method))
			throw new IllegalArgumentException("...must choose DWT or MODWT...");
		if (null != boundary && !"periodic".equals(boundary) && !"reflection".equals(boundary))
			throw new IllegalArgumentException("...boundary must be 'periodic' or 'reflection'...");
		if ("dwt".equals(method))
			{
			int t = 1;
			while (t < N) t *= 2;
			if (t != N)
				throw new IllegalArgumentException("...data must have dyadic length for DWT...");
			}

		/* The choice of boundary methods affect things... */
		if ("reflection".equals(boundary))
			{
			length = 2 * N;
			scale = length;
			Vin = new double[length];
			reflect_vector(X, N, Vin);
			}
		else
			{
			length = N;
			scale = N;
			Vin = X;
			}

		Vout = new double[length];
		for (int k = 1; k <= K; k++)
			{
			if ("dwt".equals(method))
				{
				Xout[k-1] = realloc(Xout[k-1], scale/2);
				if (k == K)
					Vout = Xout[k] = realloc(Xout[k], scale/2);
				else
					Vout = realloc(Vout, scale/2);
				dwt(Vin, scale, f, Xout[k-1], Vout);
				}
			else
				{
				Xout[k-1] = realloc(Xout[k-1], length);
				if (k == K)
					Vout = Xout[k] = realloc(Xout[k], length);
				else
					Vout = realloc(Vout, length);
				modwt(Vin, length, k, f, Xout[k-1], Vout);
				}

			scale /= 2;

			// swap input and output vectors (don't reuse X)
			double[] t = Vin == X ? null : Vin; Vin = Vout; Vout = t;
			}
		return Xout;
		}


	public static float[][] decompose(float[] X, int N, int K, Filter f, String method, String boundary, float[][] Xout)
		{
		int scale, length;
		float[] Vin, Vout;

		if (null == Xout)
			Xout = new float[K+1][];

		if (!"dwt".equals(method) && !"modwt".equals(method))
			throw new IllegalArgumentException("...must choose DWT or MODWT...");
		if (null != boundary && !"periodic".equals(boundary) && !"reflection".equals(boundary))
			throw new IllegalArgumentException("...boundary must be 'periodic' or 'reflection'...");
		if ("dwt".equals(method))
			{
			int t = 1;
			while (t < N) t *= 2;
			if (t != N)
				throw new IllegalArgumentException("...data must have dyadic length for DWT...");
			}

		/* The choice of boundary methods affect things... */
		if ("reflection".equals(boundary))
			{
			length = 2 * N;
			scale = length;
			Vin = new float[length];
			reflect_vector(X, N, Vin);
			}
		else
			{
			length = N;
			scale = N;
			Vin = X;
			}

		Vout = new float[length];
		for (int k = 1; k <= K; k++)
			{
			if ("dwt".equals(method))
				{
				Xout[k-1] = realloc(Xout[k-1], scale/2);
				if (k == K)
					Vout = Xout[k] = realloc(Xout[k], scale/2);
				else
					Vout = realloc(Vout, scale/2);
				dwt(Vin, scale, f, Xout[k-1], Vout);
				}
			else
				{
				Xout[k-1] = realloc(Xout[k-1], length);
				if (k == K)
					Vout = Xout[k] = realloc(Xout[k], length);
				else
					Vout = realloc(Vout, length);
				modwt(Vin, length, k, f, Xout[k-1], Vout);
				}

			scale /= 2;

			// swap input and output vectors (don't reuse X)
			float[] t = Vin == X ? null : Vin; Vin = Vout; Vout = t;
			}
		return Xout;
		}


	/**
	 * Peform a multirsolution analysis using the DWT or MODWT matrix
	 * obtained from `decompose.'  The inverse Transform will be applied
	 * to selected wavelet detail coefficients.  The wavelet smooth
	 * coefficients from the original Transform are added to the K+1
	 * column in order to preserve the additive decomposition.
	 * Reflection is used as the default boundary condition.
	 *
	 * @param Xin    matrix from `decompose'
	 * @param N      number of rows in Xin
	 * @param K      number of details in Xin
	 * @param f      wavelet Filter structure
	 * @param method character string (either "dwt" or "modwt")
	 * @param Xmra   vectors to be used to return results, may be null
	 * @return matrix containg K wavelet details and 1 wavelet smooth
	 */
	public static double[][] multiresolution(double[][] Xin, int N, int K, Filter f,
	                                   String method, String boundary, double[][] Xmra)
		{
		if (null == Xmra)
			Xmra = new double[K+1][];
		int length;

		if (!"dwt".equals(method) && !"modwt".equals(method))
			throw new IllegalArgumentException("...must choose DWT or MODWT...");
		if (null != boundary && !"periodic".equals(boundary) && !"reflection".equals(boundary))
			throw new IllegalArgumentException("...boundary must be 'periodic' or 'reflection'...");

		if ("reflection".equals(boundary))
			length = 2 * N;
		else
			length = N;

		double[] zero = new double[length];
		double[] Xout = new double[length];
		double[] Win = new double[length];

		for (int k = 1; k <= K; k++)
			{
			if ("dwt".equals(method))
				idwt(Xin[k-1], zero, N / pow2(k), f, Xout);
			else
				imodwt(Xin[k-1], zero, length, k, f, Xout);

			for (int i = k - 1; i >= 1; i--)
				{
				// swap arrays
				double[] t = Win; Win = Xout; Xout = t;
				if ("dwt".equals(method))
					idwt(zero, Win, N / pow2(i), f, Xout);
				else
					imodwt(zero, Win, length, i, f, Xout);
				}

			Xmra[k-1] = realloc(Xmra[k-1], length);
			arraycopy(Xout, Xmra[k-1]); // CONSIDER: avoid this copy?
			}

		/* One more iteration is required on the wavelet smooth coefficients
		   to complete the additive decomposition. */

		if ("dwt".equals(method))
			idwt(zero, Xin[K], N / pow2(K), f, Xout);
		else
			imodwt(zero, Xin[K], length, K, f, Xout);

		for (int i = K - 1; i >= 1; i--)
			{
			// swap arrays
			double[] t = Win; Win = Xout; Xout = t;
			if ("dwt".equals(method))
				idwt(zero, Win, N / pow2(i), f, Xout);
			else
				imodwt(zero, Win, length, i, f, Xout);
			}
		Xmra[K] = realloc(Xmra[K], Xout.length); // too big if "dwt"?
		arraycopy(Xout, Xmra[K]);
		return Xmra;
		}


	public static float[][] multiresolution(float[][] Xin, int N, int K, Filter f,
	                                   String method, String boundary, float[][] Xmra)
		{
		if (null == Xmra)
			Xmra = new float[K+1][];
		int length;

		if (!"dwt".equals(method) && !"modwt".equals(method))
			throw new IllegalArgumentException("...must choose DWT or MODWT...");
		if (null != boundary && !"periodic".equals(boundary) && !"reflection".equals(boundary))
			throw new IllegalArgumentException("...boundary must be 'periodic' or 'reflection'...");

		if ("reflection".equals(boundary))
			length = 2 * N;
		else
			length = N;

		float[] zero = new float[length];
		float[] Xout = new float[length];
		float[] Win = new float[length];

		for (int k = 1; k <= K; k++)
			{
			if ("dwt".equals(method))
				idwt(Xin[k-1], zero, N / pow2(k), f, Xout);
			else
				imodwt(Xin[k-1], zero, length, k, f, Xout);

			for (int i = k - 1; i >= 1; i--)
				{
				// swap arrays
				float[] t = Win; Win = Xout; Xout = t;
				if ("dwt".equals(method))
					idwt(zero, Win, N / pow2(i), f, Xout);
				else
					imodwt(zero, Win, length, i, f, Xout);
				}

			Xmra[k-1] = realloc(Xmra[k-1], length);
			arraycopy(Xout, Xmra[k-1]); // CONSIDER: avoid this copy?
			}

		/* One more iteration is required on the wavelet smooth coefficients
		   to complete the additive decomposition. */

		if ("dwt".equals(method))
			idwt(zero, Xin[K], N / pow2(K), f, Xout);
		else
			imodwt(zero, Xin[K], length, K, f, Xout);

		for (int i = K - 1; i >= 1; i--)
			{
			// swap arrays
			float[] t = Win; Win = Xout; Xout = t;
			if ("dwt".equals(method))
				idwt(zero, Win, N / pow2(i), f, Xout);
			else
				imodwt(zero, Win, length, i, f, Xout);
			}
		Xmra[K] = realloc(Xmra[K], Xout.length); // too big if "dwt"?
		arraycopy(Xout, Xmra[K]);
		return Xmra;
		}


	public static void thresholdHard(double[][] Xin, double[] threshold)
		{
		int K = Xin.length - 1;

		for (int k = 0; k <= K; k++)
			{
			double[] s = Xin[k];
			int N = s.length;
			double t = threshold[k];
			if (0.0 == t)
				continue;
			if (Double.MAX_VALUE == t)
				{
				arrayzero(s);
				continue;
				}
			for (int i = 0; i < N; i++)
				{
				double d = s[i];
				double dm = Math.abs(d);
				if (dm <= t)
					s[i] = 0.0;
				}
			}
		}


	// threshold using different threshold value for each value of s
	public static void thresholdSoft(double[] S, double[] T)
		{
		int N = S.length;
		for (int i = 0; i < N; i++)
			{
			double t = T[i];
			double s = S[i];
			S[i] = Math.abs(s) <= t ? 0.0 : s < 0 ? s + t : s - t;
			}
		}


	public static void thresholdSoft(float[] S, float[] T)
		{
		int N = S.length;
		for (int i = 0; i < N; i++)
			{
			float t = T[i];
			float s = S[i];
			S[i] = Math.abs(s) <= t ? 0.0F : s < 0 ? s + t : s - t;
			}
		}


	public static void thresholdSoft(double[] S, double t)
		{
		if (0.0 == t)
			return;
		int N = S.length;
		if (Double.MAX_VALUE == t)
			{
			arrayzero(S);
			return;
			}
		for (int i = 0; i < N; i++)
			{
			double s = S[i];
			S[i] = Math.abs(s) <= t ? 0.0 : s < 0 ? s + t : s - t;
			}
		}


	public static void thresholdSoft(double[][] Xin, double[] threshold)
		{
		int K = Xin.length - 1;

		for (int k = 0; k <= K; k++)
			{
			thresholdSoft(Xin[k], threshold[k]);
			}
		}


	static void arraycopy(double[] src, double[] dst)
		{
		assert src.length == dst.length;
		System.arraycopy(src, 0, dst, 0, src.length);
		}


	static void arraycopy(float[] src, float[] dst)
		{
		assert src.length == dst.length;
		System.arraycopy(src, 0, dst, 0, src.length);
		}


	public static void arrayzero(double[] zero)
		{
		Arrays.fill(zero, 0, zero.length, 0.0);
		}


	public static void arrayzero(float[] zero)
		{
		Arrays.fill(zero, 0, zero.length, 0.0F);
		}



	/**
	 * Helper to facilitate reusing arrays and hopefully reduce allocations
	 *
	 * @param array
	 * @param length
	 * @return
	 */
	public static double[] realloc(double[] array, int length)
		{
		if (null == array || array.length < length)
			array = new double[length];
		boolean debug = false;
		assert true == (debug = true);
		if (debug)
			Arrays.fill(array, Double.NaN);
		return array;
		}


	public static float[] realloc(float[] array, int length)
		{
		if (null == array || array.length < length)
			array = new float[length];
		boolean debug = false;
		assert true == (debug = true);
		if (debug)
			Arrays.fill(array, Float.NaN);
		return array;
		}


	public static void PrintMatrix(PrintStream out, double[] signal, double[][] mra, int n, int levels, int start, int end)
		{
		StringBuffer buf = new StringBuffer();
		start = Math.max(1, start);
		end = Math.min(n, end);
		for (int i = start; i <= end; i++)
			{
			double s = 0.0;
			if (null != signal)
				out.print(format(signal[i], buf));
			for (int j = 1; j <= levels + 1; j++)
				{
				s += mra[j][i];
				out.print("  " + format(mra[j][i], buf));
				}
			out.println("  " + format(s, buf));
			}
		}


	final static NumberFormat _format = NumberFormat.getNumberInstance();
	final static FieldPosition _pos = new FieldPosition(NumberFormat.INTEGER_FIELD);

	static
	{
	_format.setMaximumFractionDigits(4);
	_format.setMinimumFractionDigits(4);
	}

	private static StringBuffer format(double d, StringBuffer buf)
		{
		//	str = " " + str;
		buf.setLength(0);
		buf.append(d < 0 ? '-' : ' ');
		_format.format(Math.abs(d), buf, _pos);
		return buf;
		}


	static final int pow2(int k)
		{
		return 1 << k;
		}
	}