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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import modwt.Filter;
import modwt.Transform;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.datastructure.IntegerArray;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;

import javax.media.jai.JAI;
import javax.media.jai.PlanarImage;
import javax.media.jai.RasterFactory;
import javax.media.jai.RenderedOp;
import javax.media.jai.operator.DFTDescriptor;

import java.awt.Point;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferFloat;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.awt.image.renderable.ParameterBlock;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

/**
 * User: mbellew
 * Date: Jun 2, 2004
 * Time: 1:57:21 PM
 */
public class Spectrum
{
    //private static Logger _log = Logger.getLogger(Spectrum.class);

    public static final float HYDROGEN_ION_MASS = (float) (1.0078250 - 5.485e-4);  //h - electron
    public static final double LN2 = Math.log(2.0);


    /**
     *for micromass, doesn't work for bruker
     public static float[][] ExpandZeroes(float[][] spectrum)
     {
     final float ZERO = 0F; 			// set to -2 for debugging
     final float ZEROFILL = 0F;     	// set to -1 for debugging
     final int notZero = 0;
     final int oneZero = 1;
     final int twoZero = 2;

     float[] mzIN = spectrum[0];
     float[] intIN = spectrum[1];
     FloatArray mzOUT = new FloatArray();
     FloatArray intOUT = new FloatArray();

     int state = notZero;
     float interval = Math.min(0.5F, mzIN[1] - mzIN[0]);
     for (int i=0 ; i<mzIN.length ; i++)
     {
     float mz = mzIN[i];
     float f = intIN[i];
     if (f != 0)
     {
     state = notZero;
     mzOUT.add(mz);
     intOUT.add(f);
     }
     else if (state == notZero)
     {
     state = oneZero;
     mzOUT.add(mz);
     intOUT.add(ZERO);
     }
     else
     {
     state = twoZero;
     if (i >= 2)
     interval = mzIN[i-1] - mzIN[i-2];
     float mzFrom = mzIN[i-1];
     float mzTo = mzIN[i] - interval/2;
     for (float m = mzFrom+interval ; m < mzTo ; m+=interval)
     {
     mzOUT.add(m);
     intOUT.add(ZEROFILL);
     }
     mzOUT.add(mz);
     intOUT.add(ZERO);
     }
     }
     return new float[][] {mzOUT.toArray(null), intOUT.toArray(null)};
     }
     */



    /**
     * For both micromass and bruker, intervals between readings are proportional to sqrt(mz).
     * This function infers the sampling frequency from a spectrum.
     *
     * r.min should be an actual start value from a real spectrum (e.g. not just 400)
     *
     * This is expensive so it's done once per MSRun, see MSRun.getTemplateSpectrum().
     * Since we do it once, don't try to optimize too much.
     */
    public static float[] GenerateSpectrumTemplate(float[] mzArray, FloatRange r)
    {
        // convert to sqrt(mz) domain
        double start_sqrt = Math.sqrt(r.min);
        double end_sqrt = Math.sqrt(r.max);
        double[] mzArray_sqrt = new double[mzArray.length];
        for (int i=0 ; i<mzArray.length ; i++)
            mzArray_sqrt[i] = Math.sqrt((double)mzArray[i]);

/*
		double gaps[] = new double[mzArray.length-1];
		int lenGaps = 0;
		for (int i=1 ; i<mzArray.length ; i++)
			{
			double d = mzArray_sqrt[i] - mzArray_sqrt[i-1];
			if (d > 0 && d < .025)
				gaps[lenGaps++] = d;
			}
		Arrays.sort(gaps, 0, lenGaps);
*/
        // find smallest gap
        double min = 0.001;
        for (int i=1 ; i<mzArray_sqrt.length ; i++)
        {
            double d = mzArray_sqrt[i] - mzArray_sqrt[i-1];
            if (d < min)
                min = d;
        }

        // find average of small gaps
        double sum = 0.0;
        int count = 0;
        double cutoff = min * 1.5;
        for (int i=1 ; i<mzArray_sqrt.length ; i++)
        {
            double d = mzArray_sqrt[i] - mzArray_sqrt[i-1];
            if (d < cutoff)
            {
                sum += d;
                count++;
            }
        }

        double interval_sqrt = sum/count;
        //double interval_upper = interval_sqrt * 1.5;
        double interval_lower = interval_sqrt * 0.5;

        // We have a really good interval value, but we'll drift from reported mz values if we use it directly.
        // "merge" input data and calculate missing values
        int lenTemp = (int)Math.floor((end_sqrt-start_sqrt)/interval_sqrt) + 200;
        double[] template_sqrt = new double[lenTemp];

        int dst = 0, src = 0;
        mzArray_sqrt[0] = start_sqrt;
        mzArray_sqrt[mzArray_sqrt.length-1] = end_sqrt;
        double prev = start_sqrt;
        for ( ; src<mzArray_sqrt.length ; src++)
        {
            double next = mzArray_sqrt[src];
            double gap = next - prev;
            int c = (int)((gap + interval_lower) / interval_sqrt);
            assert c<2 || gap/c > interval_sqrt*0.9 && gap/c < interval_sqrt*1.1;
            for (int i=1 ; i<c ; i++)
                template_sqrt[dst++] = prev + i * gap/c;
            template_sqrt[dst++] = prev = next;
        }
        int len = dst;

        // convert back to mz domain
        float[] template = new float[len];
        for (int i=0 ; i<len ; i++)
        {
            double mz = template_sqrt[i] * template_sqrt[i];
            template[i] = (float)(Math.round(10000.0 * mz) / 10000.0);
            //System.out.println(template[i]);
        }

        return template;
    }


    static int find(float[] s, float f)
    {
        int i = Arrays.binarySearch(s, f);
        i = i<0 ? -(i+1) : i;
        i = Math.min(s.length-1, i);
        return i;
    }


    public static float[][] CombineRawSpectra(float[] template, List<float[][]> spectra, FloatRange rangeIN)
    {
        int from = find(template, rangeIN.min);
        if (from > 0 && template[from] > rangeIN.min)
            from--;
        int to = find(template, rangeIN.max);
        if (to == template.length)
            to--;
        float[] mzArray = new float[to-from+1];
        System.arraycopy(template, from, mzArray, 0, mzArray.length);

        // now combine intensity values
        float[] iArray = new float[mzArray.length];
        for (int s=0 ; s<spectra.size() ; s++)
        {
            float[][] spectrum = spectra.get(s);
            from = find(spectrum[0], rangeIN.min);
            to = find(spectrum[0], rangeIN.max);
            if (spectrum[0][to] > rangeIN.max)
                --to;
            for (int i=from, j=0 ; i<to ; i++)
            {
                float f = spectrum[0][i];
                while (mzArray[j] < f)
                    j++;
                // which is closer mzArray[j] or mzArray[j-1];
                if (j == 0 || mzArray[j] - f  < f - mzArray[j-1])
                    iArray[j] += spectrum[1][i];
                else
                    iArray[j-1] += spectrum[1][i];
            }
        }
        int sz = spectra.size();
        for (int i=0 ; i<iArray.length ; i++)
            iArray[i] /= sz;
        return new float[][] {mzArray, iArray};
    }


    /**
     * resample and translate to zero charge domain
     * assumes the input is higher resolution than 1/scale*resolution
     * <p/>
     * CONSIDER: There's probably a better, curve estimating/fitting algorithm
     */
    public static float[][] TranslateZeroCharge(float[][] spectrum, FloatRange r, int charge, int resolution)
    {
        if (spectrum.length != 2)
            throw new IllegalArgumentException();
        if (spectrum[0].length != spectrum[1].length)
            throw new IllegalArgumentException();

        float[][] out = new float[2][resolution * ((int) r.max - (int) r.min) + 1];
        int lenOut = out[0].length;
        for (int i = 0; i < lenOut; i++)
            out[0][i] = r.min + i * 1.0F / resolution;

        int len = spectrum[0].length;
        int count = 0;
        float sum = 0.0F;
        int indexLast = -1;
        for (int i = 0; i < len; i++)
        {
            float m = spectrum[0][i];
            m = (m - HYDROGEN_ION_MASS) * charge;
            if (m < r.min - 1 || m > r.max + 1)
                continue;
            float f = spectrum[1][i];
            int index = Math.round((m - r.min) * resolution);
            if (index < 0 || index >= lenOut)
                continue;
            if (index != indexLast)
            {
                if (count > 0)
                    out[1][indexLast] = sum / count;
                indexLast = index;
                sum = 0.0F;
                count = 0;
            }
            sum += f;
            count++;
        }
        if (count > 0)
            out[1][indexLast] = sum / count;
        return out;
    }


    public static float[][] ResampleSpectrum(float[][] spectrum, FloatRange r, int resolution, boolean zeroCharge)
    {
        float[] signal = Resample(spectrum, r, resolution);
        float[] mz = new float[signal.length];
        double in = 1.0 / resolution;
        double start = r.min;
        if (zeroCharge)
            start -= HYDROGEN_ION_MASS;
        for (int i = 0; i < signal.length; i++)
            mz[i] = (float) (start + i * in);
        return new float[][]{mz, signal};
    }


    static CPUTimer timerResample = new CPUTimer("Spectrum.Resample()");

    public static float[] Resample(float[][] spectrum, FloatRange r, int resolution)
    {
        try
        {
            assert timerResample.start();
            return Resample2(spectrum, r, resolution);
        }
        finally
        {
            assert timerResample.stop();
        }
    }


    /**
     * resample
     * assumes the input is higher resolution than 1/scale*resolution
     * <p/>
     * CONSIDER: There's probably a better, curve estimating/fitting algorithm
     */
    private static float[] Resample1(float[][] spectrum, FloatRange r, int resolution)
    {
        if (spectrum.length != 2)
            throw new IllegalArgumentException();
        if (spectrum[0].length != spectrum[1].length)
            throw new IllegalArgumentException();

        // ExpandZeroes
        // It's probably correct to ExpandZeros here, however
        //  a) it is expensive
        //  b) it only affects small noise and edge of peaks
        // I don't think it actually affects results

        float[] out = new float[resolution * ((int) r.max - (int) r.min) + 1];
        int lenOut = out.length;

        int count = 0;
        float sum = 0.0F;
        int indexLast = -1;

        int start = Arrays.binarySearch(spectrum[0], r.min-1);
        start = start < 0 ? -(start+1) : start;
        int end = Arrays.binarySearch(spectrum[0], r.max+1);
        end = end < 0 ? -(end+1) : end;

        for (int i=start ; i<end ; i++)
        {
            float m = spectrum[0][i];
            float x = spectrum[1][i];
            int index = Math.round((m - r.min) * resolution);
            if (index < 0 || index >= lenOut)
                continue;
            if (index != indexLast)
            {
                if (count > 0)
                    out[indexLast] = sum / count;
                indexLast = index;
                sum = 0.0F;
                count = 0;
            }
            sum += x;
            count++;
        }
        if (count > 0)
            out[indexLast] = sum / count;
        return out;
    }


    private static float[] Resample2(float[][] spectrum, FloatRange r, int resolution)
    {
        if (spectrum.length != 2)
            throw new IllegalArgumentException();
        if (spectrum[0].length != spectrum[1].length)
            throw new IllegalArgumentException();

        // ExpandZeroes
        // It's probably correct to ExpandZeros here, however
        //  a) it is expensive
        //  b) it only affects small noise and edge of peaks
        // I don't think it actually affects results

        float[] out = new float[resolution * ((int) r.max - (int) r.min) + 1];
        float[] weights = new float[out.length];

        int start = Arrays.binarySearch(spectrum[0], r.min-1F/resolution);
        start = start < 0 ? -(start+1) : start;
        int end = Arrays.binarySearch(spectrum[0], r.max+1F/resolution);
        end = end < 0 ? -(end+1) : end;

        for (int i=start ; i<end ; i++)
        {
            float m = spectrum[0][i];
            float x = spectrum[1][i];
            float bucket = (m - r.min) * resolution;
            int index = (int)Math.floor(bucket);
            float frac = bucket-index;

            if (index >= 0 && index < out.length)
            {
                out[index] += x*(1-frac);
                weights[index] += (1-frac);
            }
            if (index+1 >= 0 &&index+1 < out.length)
            {
                out[index+1] += x*frac;
                weights[index+1] += frac;
            }
        }
        for (int i=0 ; i<out.length; i++)
        {
            assert weights[i] != 0 || out[i] == 0;
            out[i] /= Math.max(1.0,weights[i]);
        }

        // this algorithm does tend to smooth nicely
        // however, you sometimes get a largish peak right on
        // a bucket boundary, this won't smooth very well.
        SmoothALittle(out);

        return out;
    }


    /**
     * don't put too much effort in this yet,
     * Are there standard centroiding algorithms?
     * I think Tim Randolph has some ideas
     */

    public static int[] PickPeakIndexes(float[] signal, double minFilter)
    {
        IntegerArray arr = new IntegerArray();
        int len = signal.length;

        float prev = -Float.MIN_VALUE;
        float curr;
        for (int i = 0; i < len;)
        {
            for (; i < len && prev <= (curr = signal[i]); i++)
                prev = curr;
            if (prev > minFilter)
                arr.add(i - 1);
            for (; i < len && prev >= (curr = signal[i]); i++)
                prev = curr;
            assert i==len || !Float.isNaN(signal[i]);
        }
        return arr.toArray(null);
    }


    public static Peak[] PickPeaks(float[][] spectrum, double minFilter)
    {
        List<Peak> arr = new ArrayList<Peak>();
        float[] signal = spectrum[1];
        int len = spectrum[0].length;

        float prev = -Float.MAX_VALUE;
        float curr;
        for (int i = 0; i < len;)
        {
            for (; i < len && prev <= (curr = signal[i]); i++)
                prev = curr;
            if (prev > minFilter)
                arr.add(new CentroidedPeak(spectrum[0][i - 1], spectrum[1][i - 1]));
            for (; i < len && prev >= (curr = signal[i]); i++)
                prev = curr;
        }
        return (Peak[]) arr.toArray(new Peak[arr.size()]);
    }


    public static void SmoothALittle(float[] x)
    {
        float last = x[0];
        float cur = x[0];
        int stop = x.length - 1;
        float next;
        for (int i = 0; i < stop; i++)
        {
            next = x[i + 1];
            x[i] = (last + cur + cur + next) / 4.0F;
            last = cur;
            cur = next;
        }
        next = x[stop];
        x[stop] = (last + cur + cur + next) / 4.0F;
    }


    public static void SmoothALittle(double[] x)
    {
        double last = x[0];
        double cur = x[0];
        int stop = x.length - 1;
        double next;
        for (int i = 0; i < stop; i++)
        {
            next = x[i + 1];
            x[i] = (last + cur + cur + next) / 4.0F;
            last = cur;
            cur = next;
        }
        next = x[stop];
        x[stop] = (last + cur + cur + next) / 4.0F;
    }


    public static Comparator<Peak> comparePeakIntensityDesc = new Comparator<Peak>()
    {
        public int compare(Peak a, Peak b)
        {
            //if (a == null || b == null) return a==null ? 1 : -1;
            return _compareDesc(a.intensity, b.intensity);
        }
    };

    public static Comparator<Peak> comparePeakMzAsc = new Comparator<Peak>()
    {
        public int compare(Peak a, Peak b)
        {
            return _compareAsc(a.mz, b.mz);
        }
    };

    public static Comparator<Feature> comparePeakMassAsc = new Comparator<Feature>()
    {
        public int compare(Feature a, Feature b)
        {
            return _compareAsc(a.mass, b.mass);
        }
    };

    public static Comparator<Peak> comparePeakScanAsc = new Comparator<Peak>()
    {
        public int compare(Peak a, Peak b)
        {
            return _compareAsc(a.scan, b.scan);
        }
    };

    public static Comparator<Feature> compareFeatureLengthDesc = new Comparator<Feature>()
    {
        public int compare(Feature a, Feature b)
        {
            return _compareAsc(a.getScanCount(), b.getScanCount());
        }
    };



    static class FloatImage extends PlanarImage
    {
        WritableRaster raster;


        FloatImage(int w, int h, DataBuffer buf)
        {
            tileWidth = width = w;
            tileHeight = height = h;
            tileGridXOffset = tileGridYOffset = minX = minY = 0;
            assert buf.getSize() == w * h;
            sampleModel = RasterFactory.createBandedSampleModel(DataBuffer.TYPE_FLOAT, w, h, 1);
            raster = RasterFactory.createWritableRaster(sampleModel, buf, new Point(0, 0));
        }


        public Raster getData()
        {
            return raster;
        }


        public Raster getTile(int tileX, int tileY)
        {
            return raster;
        }
    }


    public static float[] FFTsmooth(float[] x, double smoothfactor, boolean cliff)
    {
        try
        {
            DataBuffer buf = new DataBufferFloat(x, x.length);
            FloatImage img = new FloatImage(x.length, 1, buf);
            RenderedOp opDFT = JAI.create("DFT", new ParameterBlock()
                    .add(DFTDescriptor.SCALING_UNITARY)
                    .add(DFTDescriptor.REAL_TO_COMPLEX)
                    .addSource(img));

            float[] fft = ((DataBufferFloat) opDFT.getData().getDataBuffer()).getData();
            int complexLen = fft.length / 2;

            // UNDONE cache, this seems overexpensive
            /*
               int len = complexLen / 2;
               float[] g = new float[len];
               double s = len / smoothfactor;
               double s2 = -1.0 / (2.0 * s * s);
               int z = cliff ? (int) s : 0;
               for (int i = 0; i < z; i++)
                   g[i] = 1F;
               for (int i = 0; i < len - z; i++)
                   {
                   double d = Math.exp((double)i * i * s2);
                   if (d < Float.MIN_VALUE) // arbitrary cut-off for small..
                       break;
                   g[z + i] = (float) d;
                   }
               */
            int len = complexLen / 2;
            float[] g = NormalP(len, 1/smoothfactor, false, false, null);

            int off = 2;
            for (int i = 0; i < len; i++)
            {
                fft[off++] *= g[i];	// Real
                fft[off++] *= g[i]; // Complex
            }
            for (int i = len - 2; i >= 0; i--)
            {
                fft[off++] *= g[i];	// Real
                fft[off++] *= g[i]; // Complex
            }
            assert off == fft.length;

            RenderedOp opIDFT = JAI.create("IDFT", new ParameterBlock()
                    .add(DFTDescriptor.SCALING_UNITARY)
                    .add(DFTDescriptor.COMPLEX_TO_REAL)
                    .addSource(opDFT));
            DataBufferFloat bufReturn = (DataBufferFloat) opIDFT.getData().getDataBuffer();
            return bufReturn.getData();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
            return null;
        }
    }

    /**
     * Computes points in Gaussian distribution function, over range [0,1]
     * with mean=0
     *
     * @param len       number of points to compute in the distrubition (min 2)
     * @param sigma     sqrt(variance)
     * @param g         optional array to copy result into
     * @return
     */
    private static float[] NormalP(int len, double sigma, boolean even, boolean normalize, float[] g)
    {
        g = realloc(g, len);
        double variance = sigma * sigma;

        double delta;// = 1.0/(len-1);
        double x; // = 0;

        if (even)
        {
            delta = 1.0/(len-0.5);
            x = delta/2.0;
        }
        else
        {
            delta = 1.0/(len-1.0);
            x = 0.0;
        }

        double s2 = -1.0 / (2.0 * variance);
        double scale = normalize ? 1/(sigma * Math.sqrt(2 * Math.PI)) : 1;
        for ( int i=0 ; i<len ; i++, x+= delta)
        {
            double d = scale * Math.exp((double)x * x * s2);
            if (d < Float.MIN_VALUE) // arbitrary cut-off for small, avoid denormals
                break;
            g[i] = (float) d;
        }
        return g;
    }


    /*
      * Computes points in Gaussian distribution function, over range [0,1]
      * with mean=0.5
      *
      * @param len       number of points to compute in the distrubition
      * @param sigma     sqrt(variance)
      * @param g         optional array to copy result into
      * @return
     private static float[] NormalD(int len, double sigma, float[] g)
         {
         assert null == "NYI";
         return null;
         }
      */


    public static double HellingerDistance(float[] p, float[] q)
    {
        assert p.length == q.length;
        double d = 0.0;
        for (int i = 0; i < p.length; i++)
        {
            double x = Math.sqrt(p[i]) - Math.sqrt(q[i]);
            d += x * x;
        }
        return d;
    }


    public static float KLPoissonDistance(float mass, float[] q)
    {
        float[] p = Poisson(mass);
        assert q.length == p.length;

        double diff = 0.0;
        double sumP = 0.0;
        double sumQ = 0.0;
        int pqMinLength = Math.min(p.length, q.length);
        for (int k = 0; k < pqMinLength ; k++)
        {
            diff += p[k] * Math.log((double)p[k] / q[k]);
            sumP += p[k];
            sumQ += q[k];
        }
        double kl = diff / LN2;

        assert kl > -0.0001;
        assert Math.abs(sumP-1.0) < 0.001 && Math.abs(sumQ-1.0) < 0.001;

        kl = Math.max(kl,0);
        return (float)(diff / Math.log(2.0));
    }


    public static float KLPoissonDistanceSymmetric(float mass, float[] signal)
    {
        float[] p = Poisson(mass);
        float diff = 0.0F;
        for (int k = 0; k < 4; k++)
        {
            float f0 = p[k];
            float f1 = signal[k];
            diff += f0 * Math.log(f0 / f1) + f1 * Math.log(f1 / f0);
        }
        return (float)(diff / Math.log(2.0));
    }


    public static double KLDistanceSymmetric(float[] signal1, int off1, float[] signal2, int off2, int len)
    {
        double sum1 = 0.0, sum2 = 0.0;
        for (int i = 0; i < len; i++)
        {
            sum1 += signal1[off1 + i];
            sum2 += signal2[off2 + i];
        }

        double diff = 0.0;
        double min = 0.001 / len;
        for (int i = 0; i < len; i++)
        {
            double s1 = Math.max(signal1[off1 + i] / sum1, min);
            double s2 = Math.max(signal2[off2 + i] / sum2, min);
            diff += s1 * Math.log(s1 / s2) + s2 * Math.log(s2 / s1);
        }
        return diff / Math.log(2.0);
    }


    public static float KLGayDistance(float m, float[] signal)
    {
        float[] p = _gay(m);
        float diff = 0.0F;
        for (int k = 0; k < 4; k++)
        {
            float f0 = p[k];
            float f1 = signal[k];
            diff += f0 * Math.log(f0 / f1) + f1 * Math.log(f1 / f0);
        }
        return diff;
    }


    /* CONSIDER: Better mathematically justified ways to compute noise? */
    public static float Noise(float[] signal, int start, int end)
    {
        int len = end - start;
        float[] copy = new float[len];
        System.arraycopy(signal, start, copy, 0, len);
        Arrays.sort(copy);
        int quartile = len * 3 / 4;
        return copy[quartile];
    }



    static float[][] _poisson = new float[6400 / 10][6];


    static
    {
        for (int i = 0; i < _poisson.length; i++)
        {
            int mass = i * 10 + 5;
            float mu = mass * 0.0005556F; // martin's empircal lambda = 1/1800 (was mass * 0.000489F)
            double muExp = Math.exp(-mu);
            _poisson[i][0] = (float) (1 * muExp / 1);
            _poisson[i][1] = (float) (mu * muExp / 1);
            _poisson[i][2] = (float) (Math.pow(mu, 2) * muExp / 2);
            _poisson[i][3] = (float) (Math.pow(mu, 3) * muExp / 6);
            _poisson[i][4] = (float) (Math.pow(mu, 4) * muExp / 24);
            _poisson[i][5] = (float) (Math.pow(mu, 5) * muExp / 120);
            float s = _poisson[i][0] + _poisson[i][1] + _poisson[i][2] + _poisson[i][3] + _poisson[i][4] + _poisson[i][5];
            _poisson[i][0] /= s;
            _poisson[i][1] /= s;
            _poisson[i][2] /= s;
            _poisson[i][3] /= s;
            _poisson[i][4] /= s;
            _poisson[i][5] /= s;
        }
    }


    public static float[] Poisson(float m)
    {
        int i = ((int) m - 5) / 10;
        i = Math.max(0, Math.min(_poisson.length - 1, i));
//System.err.println("**Poisson: m=" + m + ", i=" + i + ", p1=" + _poisson[i][0] + ", p2=" + _poisson[i][1]);        
        return _poisson[i];
    }

    /**
     * Return the index of the peak that, according to the Poisson distribution for this mass, should be most intense
     * @param mass
     * @return
     */
    public static int calcMaxIdealPeakIndex(float mass)
    {
        //Compare intensities of the peak that has theoretical max intensity
        int theoreticalMaxIndex = 0;
        float[] p = Poisson(mass);
        float maxTheoreticalInt = p[0];
        for (int i=1; i<p.length; i++)
            if (p[i] > maxTheoreticalInt)
            {
                maxTheoreticalInt = p[i];
                theoreticalMaxIndex = i;
            }
        return theoreticalMaxIndex;
    }

    /**
     * Return the indexes of the peaks in descending order of theoretical intensity.
     * Relies on no exact ties from Poisson();
     * @param mass
     * @return
     */
    public static int[] calcIdealPeakIntensityOrderDesc(float mass)
    {
        //Compare intensities of the peak that has theoretical max intensity
        float[] p = Poisson(mass);
        float[] pSorted = new float[p.length];
        System.arraycopy(p, 0, pSorted, 0, p.length);
        Arrays.sort(pSorted);
        //sorted in ASCENDING ORDER

        int[] result = new int[p.length];
        for (int i=0; i<p.length; i++)
        {
            float match = pSorted[pSorted.length-i-1];
            for (int j=0; j<p.length; j++)
                if (p[j] == match)
                {
                    result[i] = j;
                    break;
                }
        }
        return result;
    }


    // from article in Electrophoresis 1999 by S. Gay et al
    static double[] m0 = new double[]{0.0576, -0.2553, 0.5827, -0.8436, 0.6182};
    static double[] m1 = new double[]{-0.1207, 0.4688, -0.7373, 0.3845, 0.2701};
    static double[] m2 = new double[]{0.0569, -0.1173, -0.0838, 0.3068, 0.0865};
    static double[] m3 = new double[]{0.0296, -0.1293, 0.1330, 0.1151, 0.0204};
    static double[] m4 = new double[]{-0.0101, -0.0118, 0.0789, 0.0292, 0.0040};
    static double[] m5 = new double[]{-0.0132, 0.0448, 0.0256, 0.0079, 0.0008};
    static float[][] _gay = new float[6400 / 10][6];


    static
    {
        for (int i = 0; i < _gay.length; i++)
        {
            int m = i * 10 + 5;
            double x = ((double) m - 800.0) / 2200.0;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            _gay[i][0] = (float) (m0[0] * x4 + m0[1] * x3 + m0[2] * x2 + m0[3] * x + m0[4]);
            _gay[i][1] = (float) (m1[0] * x4 + m1[1] * x3 + m1[2] * x2 + m1[3] * x + m1[4]);
            _gay[i][2] = (float) (m2[0] * x4 + m2[1] * x3 + m2[2] * x2 + m2[3] * x + m2[4]);
            _gay[i][3] = (float) (m3[0] * x4 + m3[1] * x3 + m3[2] * x2 + m3[3] * x + m3[4]);
            _gay[i][4] = (float) (m4[0] * x4 + m4[1] * x3 + m4[2] * x2 + m4[3] * x + m4[4]);
            _gay[i][5] = (float) (m5[0] * x4 + m5[1] * x3 + m5[2] * x2 + m5[3] * x + m5[4]);
            float s = _gay[i][0] + _gay[i][1] + _gay[i][2] + _gay[i][3] + _gay[i][4] + _gay[i][5];
            _gay[i][0] /= s;
            _gay[i][1] /= s;
            _gay[i][2] /= s;
            _gay[i][3] /= s;
            _gay[i][4] /= s;
            _gay[i][5] /= s;
        }
    }


    static float[] _gay(float m)
    {
        int i = ((int) m - 5) / 10;
        i = Math.max(0, Math.min(_gay.length - 1, i));
        return _gay[i];
    }


    public static float[][] Centroid(float[][] spectrum)
    {
        // this isn't exactly well defined, since spectrum is not evenly sampled
        float n = Spectrum.Noise(spectrum[1], 0, spectrum[1].length);
        Spectrum.SmoothALittle(spectrum[1]);
        Spectrum.SmoothALittle(spectrum[1]);
        Peak[] peaks = Spectrum.PickPeaks(spectrum, n);
        float[][] spectrumPeaks = new float[2][peaks.length];
        for (int i = 0; i < peaks.length; i++)
        {
            Peak p = peaks[i];
            spectrumPeaks[0][i] = p.mz;
            spectrumPeaks[1][i] = p.intensity;
        }
        return spectrumPeaks;
    }


    public static void CopyToTSV(float[][] spectrum, Writer out, boolean useHeader) throws IOException
    {
        if (useHeader)
            out.write("mz\tintensity\n");
        for (int i = 0; i < spectrum[0].length; i++)
        {
            out.write(String.valueOf(spectrum[0][i]));
            out.write('\t');
            out.write(String.valueOf(spectrum[1][i]));
            out.write('\n');
        }
    }


    public static class Peak
    {
        public int scan = -1;
        public float mz;
        public float intensity;
        public float background = 0F;
        public float median = 0F;
        public boolean excluded = false;

        public Peak()
        {
        }

        public Peak(int scan, float mz, float intensity)
        {
            this.scan = scan;
            this.mz = mz;
            this.intensity = intensity;
        }

        public Peak(float mz, float intensity)
        {
            this.mz = mz;
            this.intensity = intensity;
        }

        public Peak(Peak p)
        {
            scan = p.scan;
            mz = p.mz;
            intensity = p.intensity;
            background = p.background;
            median = p.median;
        }

        public String toString()
        {
            return scan + "\t" + mz + "\t" + intensity;
        }

        public boolean equals(Object o)
        {
            if (null == o || o.getClass() != this.getClass())
                return false;

            Peak p = (Peak) o;
            return p.scan == scan && p.mz == mz && p.intensity == intensity;
        }

        public int hashCode()
        {
            return scan ^ Float.floatToIntBits(mz) ^ Float.floatToIntBits(intensity);
        }

        public int getScan()
        {
            return scan;
        }

        public void setScan(int scan)
        {
            this.scan = scan;
        }

        public float getMz()
        {
            return mz;
        }

        public void setMz(float mz)
        {
            this.mz = mz;
        }

        public float getIntensity()
        {
            return intensity;
        }

        public void setIntensity(float intensity)
        {
            this.intensity = intensity;
        }

        public float getBackground()
        {
            return background;
        }

        public void setBackground(float background)
        {
            this.background = background;
        }

        public float getMedian()
        {
            return median;
        }

        public void setMedian(float median)
        {
            this.median = median;
        }
    }

    /*
     * this is just a marker class to distinguish between raw and centroided peaks
    */
    public static class CentroidedPeak extends Peak
    {
        public CentroidedPeak(int scan, float mz, float intensity)
        {
            super(scan, mz, intensity);
        }

        public CentroidedPeak(float mz, float intensity)
        {
            super(mz, intensity);
        }
    }


    static final int _compareAsc(float a, float b)
    {
        return a == b ? 0 : a < b ? -1 : 1;
    }


    static final int _compareDesc(float a, float b)
    {
        return a == b ? 0 : a < b ? 1 : -1;
    }


    /*
     public static boolean IsLockSpray(float[][] fullSpectrum)
         {
         float[][] spectrum = ResampleSpectrum(fullSpectrum, new FloatRange(780.0F, 7900.F), 36, false);
         spectrum[1] = Spectrum.FFTsmooth(spectrum[1], 8, false);
         Peak[] peaks = Spectrum.PickPeaks(spectrum, 20.0F);
         Arrays.sort(peaks, Spectrum.comparePeakIntensityDesc);
         if (peaks.length < 3)
             return false;
         if (Math.abs(peaks[0].mz - 785.8333) > .06)
             return false;
         if (Math.abs(peaks[1].mz - 786.3333) > .06)
             return false;
         if (Math.abs(peaks[2].mz - 786.8333) > .06)
             return false;
         return true;
         }
         */


    public static void NormalizeSum(float[] x)
    {
        double sum = 0.0;
        for (int i=0 ; i<x.length ; i++)
            sum += x[i];
        for (int i=0 ; i<x.length ; i++)
            x[i] /= sum;
    }


    public static void NormalizeSquares(float[] x)
    {
        double sum = 0.0;
        for (int i=0 ; i<x.length ; i++)
            sum += x[i] * x[i];
        double sqr = Math.sqrt(sum);
        for (int i=0 ; i<x.length ; i++)
            x[i] /= sqr;
    }


    public static double Correlation(float[] X, float[] Y)
    {
        assert X.length == Y.length;
        if (X.length != Y.length)
            throw new IllegalArgumentException("Spectrum.Correlation(), arrays need be equal length");
        int n = X.length;
        double Sx=0.0, Sx2=0.0, Sy=0.0, Sy2=0.0, Sxy=0.0;
        for (int i=0 ; i<n ; i++)
        {
            double x = X[i], y = Y[i];
            Sx += x;
            Sx2 += x*x;
            Sy += y;
            Sy2 += y*y;
            Sxy += x*y;
        }
        return (n * Sxy - Sx*Sy) / Math.sqrt((n*Sx2 - Sx*Sx) * (n*Sy2 - Sy*Sy));
    }


    /*
     static class array
         {
         float[] v;

         array(float[] v)
             {
             this.v = v;
             }

         float get(int i)
             {
             //return (i<0) ? v[0] : i>=v.length ? v[v.length-1] : v[i];
             return (i < 0 || i >= v.length) ? 0 : v[i];
             }
         }


     public static float[] average(float[] signal, int scale)
         {
         int len = scale * 2 + 1;
         array s = new array(signal);
         float[] x = new float[signal.length];

         double sum = 0.0;
         for (int i = -scale; i < signal.length; i++)
             {
             sum += s.get(i + scale) - s.get(i - scale);
             if (i > 0)
                 x[i] = (float) (sum / len);
             }
         return x;
         }
     */


    /*
      * assumes input has been preprocessed
      *   lightly smoothed
      *   removed background
      *
      * UNDONE: does not filter out small peaks in low signal area (median is zero or close to zero)
      *
      * UNDONE: CombinedStrategy uses WaveletD3(), NOT the same as WaveletPeaks()
      *
      * CONSIDER: rewrite as HaarPeaks() e.g. remove extraneous wavelet trappings. Has this already been done?
      */
    public static Peak[] WaveletPeaks(float[][] spectrum)
    {
        float[] x = spectrum[1];
        int K = 4;
        int N = spectrum[0].length;

        //
        // D1
        //
        float[] mm0 = MedianWindow(x, N, 36 * 5, false);
//		Transform.thresholdSoft(x, mm0);
        float[][] modwt = Transform.decompose(x, N, K, new Filter("haar"), "modwt", "periodic", null);
        float[] D1 = modwt[2];
        modwt[2] = new float[D1.length]; // so we can reuse modwt

        //
        // D2
        //
        //double m = MedianSampled(D1, true);
        //Transform.thresholdSoft(D1, m);
        Transform.decompose(D1, N, K, new Filter("haar"), "modwt", "periodic", modwt);
        float[] D2 = modwt[2];

        // threshold, smooth
        //double[] mm2 = MedianWindow(D2, N, 36 * 5, true);
        //Transform.thresholdSoft(D2, m);
        SmoothALittle(D2);

        // peaks are at D2 minima
        Rotate(D2, -7);
        int[] minima = FindMinimaIndexes(D2);
        List<Peak> list = new ArrayList<Peak>();
        for (int i=0 ; i<minima.length ; i++)
        {
            int curr = minima[i];
            if (x[curr] < 2*mm0[curr])
                continue;
            float mz = spectrum[0][curr];
            float in = (float)x[curr];
            Peak p = new Peak(mz, in);
            list.add(p);
        }
        return (Peak[])list.toArray(new Peak[list.size()]);
    }


    public static Peak[] WaveletPeaksD3(float[][] spectrum)
    {
        float[] m = MedianWindow(spectrum[1], spectrum[0].length, 72, false);
        float[] d3 = Spectrum.WaveletD3(spectrum[1], null);
        int[] indexes = PickPeakIndexes(d3, .5);
        List<Peak> l = new ArrayList<Peak>();
        for (int i=0 ; i<indexes.length; i++)
        {
            int p = indexes[i];
            if (spectrum[1][p] < 3 * Math.max(1,m[p]))
                continue;
            Peak peak = new Peak(0, spectrum[0][p], spectrum[1][p]);
            l.add(peak);
        }
        return (Peak[])l.toArray(new Peak[0]);
    }


    /** find smallest (most negative) values between zero crossings */
    public static int[] FindMinimaIndexes(float[] signal)
    {
        IntegerArray arr = new IntegerArray();
        int len = signal.length;

        double min;
        int pos;
        for (int i = 0; i < len;)
        {
            // find falling zero crossing
            for (; i < len && signal[i] >= 0; i++)
                ;
            if (i == len)
                break;
            pos = i;
            min = signal[pos];
            // find rising zero crossing
            for (; i < len && signal[i] <= 0; i++)
            {
                if (signal[i] < min)
                {
                    pos = i;
                    min = signal[pos];
                }
            }
            arr.add(pos);
        }
        return arr.toArray(null);
    }


    public static void Rotate(float[] x, int d)
    {
        if (d < 0)
            d += x.length;
        // kinda slow, but don't need another array
        Reverse(x, 0, x.length);
        Reverse(x, 0, d);
        Reverse(x, d, x.length-d);
    }


    public static void Reverse(float[] x, int start, int len)
    {
        float t;
        for (int i=start, j=start+len-1 ; i<j ; i++, j--)
        {
            t = x[i]; x[i] = x[j]; x[j] = t;
        }
    }


    /** returns reconstruction of last L levels of wavelet decomposition
     * not including the S 'remainder'
     *
     * this is not a real analysis method, just for visual examination
     */
    public static float[] WaveletFilter(float[] signalF, int K, int L, float[] bg)
    {
        double[] x = PadToDouble(signalF, pow2(K));
        int N = x.length;

        double[][] t = Transform.decompose(x, N, K, new Filter("haar"), "modwt", "periodic", null);
        double[][] m = Transform.multiresolution(t, N, K, new Filter("haar"), "modwt", "periodic", null);

        Arrays.fill(x, 0F);
        for (int d = K - L; d < K; d++)
        {
            double[] D = m[d];
            for (int i = 0; i < N; i++)
                x[i] += D[i];
        }
        if (bg != null)
            UnpadToFloat(m[K], pow2(K), bg);
        return UnpadToFloat(x, pow2(K), null);
    }

    /**
     * Hardcode K=3.  Uses tmp argument because of backward compatibility
     * @param X spectra
     * @param tmp Note!  This is not passed in, as you might expect.  It's ignored entirely.  This is to preserve
     * the legacy behavior of WaveletD3, which accepted but ignored the tmp argument.
     * @return transformed spectra
     */
    public static float[] WaveletD3(float[] X, Pair<float[][], float[][]> tmp)
    {
        return WaveletDX(X, null, 3);
    }

    /**
     * Hardcode K=4.  Doesn't use tmp argument because of backward compatibility
     * @param X spectra
     * @return transformed spectra
     */
    public static float[] WaveletD4(float[] X)
    {
        return WaveletDX(X, null, 4);
    }

    /**
     * Call WaveletDX and ignore the result of Transform.decompose()
     * @param X spectra
     * @param K
     * @return transformed spectra
     */
    public static float[] WaveletDX(float[] X, int K)
    {
        return WaveletDX(X, null, K);
    }

    /**
     * Perform wavelet transformation.  Return the specified level of the result of Transform.multiresolution()
     * @param X spectra
     * @param tmp  If this is supplied, populate it with the results of Transform.decompose() and
     * Transform.multiresolution()
     * @param K
     * @return transformed spectra
     */
    public static float[] WaveletDX(float[] X, Pair<float[][], float[][]> tmp, int K)
    {
        if (null == tmp)
            tmp = new Pair<float[][], float[][]>(null, null);
        int N = X.length;
        tmp.first = Transform.decompose(X, N, K, new Filter("haar"), "modwt",
                "periodic", tmp.first);
        tmp.second = Transform.multiresolution(tmp.first, N, K, new Filter("haar"), "modwt",
                "periodic", tmp.second);
        return (tmp.second)[K-1];
    }

    public static double ChooseThreshold(double[] x, double f)
    {
        if (true)
        {
            double m = MedianSampled(x, true);
            return m * f;
        }
        else
        {
            // mean
            double Sx = 0.0;
            for (int i=0 ; i<x.length ; i++)
                Sx += Math.abs(x[i]);
            return (Sx / x.length) * f;
        }
    }


    public static double Median(double[] x, int start, int len, boolean fABS, double[] t)
    {
        assert null == t || t.length >= len;
        t = realloc(t, len);
        if (!fABS)
            System.arraycopy(x, start, t, 0, len);
        else
        {
            for (int i = 0 ; i<len ; i++)
                t[i] = Math.abs(x[i+start]);
        }
        Arrays.sort(t);
        return t[len/2];
    }


    public static float Median(float[] x, int start, int len, boolean fABS, float[] t)
    {
        assert null == t || t.length >= len;
        t = realloc(t, len);
        if (!fABS)
            System.arraycopy(x, start, t, 0, len);
        else
        {
            for (int i = 0 ; i<len ; i++)
                t[i] = Math.abs(x[i+start]);
        }
        Arrays.sort(t);
        return t[len/2];
    }


    public static final float Median(float a, float b, float c)
    {
        if (a < b)
            return c < a ? a : b < c ? b : c;
        else // b < a
            return a < c ? a : c < b ? b : c;
    }


    public static float[] MedianSmooth(float[] x)
    {
        return MedianSmooth(x, x.length, null);
    }

    /** @deprecated */
    public static float[] MedianSmooth(float[] x, int len)
    {
        return MedianSmooth(x, len, null);
    }

    public static float[] MedianSmooth(float[] x, int len, float[] in)
    {
        float[] y = Spectrum.realloc(in, len);
        if (len < 3)
            return (float[])x.clone();
        for (int i=1 ; i<len-1 ; i++)
            y[i] = Median(x[i-1],x[i],x[i+1]);
        y[0] = y[1];
        y[len-1] = y[len-2];
        return y;
    }


    public static float[] MinimaWindow(float[] x, int len, int windowSize, float[] result)
    {
        result = realloc(result, len);

        float[] buckets = new float[((x.length-1) / windowSize) + 1];
        System.arraycopy(x, 0, result, 0, len);
        for (int b = 0 ; b<buckets.length ; b++)
        {
            int start = b*windowSize;
            int end = Math.min(len, (b+1)*windowSize);
            float min = Float.MAX_VALUE;
            for (int i=start ; i<end && min > 0 ; i++)
                if (result[i] < min) min = result[i];
            buckets[b] = min;
        }
        for (int b=0 ; b<buckets.length-1 ; b++)
            buckets[b] = Math.min(buckets[b], buckets[b+1]);
        for (int b=buckets.length-1 ; b>0 ; b--)
            buckets[b] = Math.min(buckets[b-1], buckets[b]);
        SmoothALittle(buckets);
        result = Interpolate(buckets, len, windowSize, result);
        return result;
    }


    public static float[] MedianWindow(float[] x, int len, int windowSize, boolean fABS)
    {
        int windows = ((len-1)/windowSize)+1;
        float[] buckets = new float[windows];
        for (int b=0 ; b<windows ; b++)
        {
            int windowStart = b * windowSize;
            int windowEnd = Math.min(x.length, (b+1) * windowSize);
            int windowLen = windowEnd - windowStart;
            buckets[b] = Median(x, windowStart, windowLen, fABS, null);
        }
        for (int b=0 ; b<buckets.length-1 ; b++)
            buckets[b] = Math.min(buckets[b], buckets[b+1]);
        for (int b=buckets.length-1 ; b>0 ; b--)
            buckets[b] = Math.min(buckets[b-1], buckets[b]);
        SmoothALittle(buckets);
        float[] result = Interpolate(buckets, len, windowSize, null);
        return result;
    }


    private static float[] Interpolate(float[] m, int len, int windowSize, float[] result)
    {
        windowSize = Math.min(len, windowSize);
        try
        {
            result = realloc(result, len);
            // fill non-interpolated part
            Arrays.fill(result, 0, windowSize, m[0]);
            Arrays.fill(result, result.length-windowSize, result.length, m[m.length-1]);
            // fill interpolated part
            for (int w=0 ; w<m.length-1 ; w++)
            {
                int windowStart = len * w/m.length + windowSize/2;
                int windowEnd = len * (w+1)/m.length + windowSize/2;
                int windowLen = windowEnd - windowStart;
                double v = m[w];
                double d = (m[w+1]-m[w]) / windowLen;
                for (int i=0; i<windowLen ; i++)
                    result[windowStart+i] = (float)(v + d*i);
            }
            return result;
        }
        catch (RuntimeException x)
        {
            x.printStackTrace();
            throw x; // just a place to put a breakpoint
        }
    }


    public static double MedianSampled(double[] x, boolean fABS)
    {
        double p = 0.5;
        int len = x.length;
        int sample = 100 + (int)Math.ceil(Math.sqrt(len));

        double[] copy = null;
        if (sample < len / 4)
        {
            copy = new double[sample];
            for (int i=0 ; i<sample ; i++)
            {
                int r = (int)Math.floor(Math.random() * len);
                copy[i] = fABS ? Math.abs(x[r]) : x[r];
            }
        }
        else
        {
            copy = new double[len];
            for (int i=0 ; i<len ; i++)
                copy[i] = fABS ? Math.abs(x[i]) : x[i];
        }
        Arrays.sort(copy);
        double l = (copy.length * p) - 1;
        int i = Math.max(0, (int)Math.floor(l));
        double m =  copy[i];
        return m;
    }


    public static float MedianSampled(float[] x, boolean fABS)
    {
        double p = 0.5;
        int len = x.length;
        int sample = 100 + (int)Math.ceil(Math.sqrt(len));

        float[] copy = null;
        if (sample < len / 4)
        {
            copy = new float[sample];
            for (int i=0 ; i<sample ; i++)
            {
                int r = (int)Math.floor(Math.random() * len);
                copy[i] = fABS ? Math.abs(x[r]) : x[r];
            }
        }
        else
        {
            copy = new float[len];
            for (int i=0 ; i<len ; i++)
                copy[i] = fABS ? Math.abs(x[i]) : x[i];
        }
        Arrays.sort(copy);
        double l = (copy.length * p) - 1;
        int i = Math.max(0, (int)Math.floor(l));
        float m =  copy[i];
        return m;
    }


    public static float[][] RemoveBackground(float[][] spectra)
    {
        int imzMax = spectra[0].length;
        float[][] background = new float[spectra.length][];
        for (int s = 0; s < spectra.length; s++)
        {
            float[] S = spectra[s];
            background[s] = Spectrum.MinimaWindow(spectra[s], spectra[s].length, 72, null);
            //background[s] = Spectrum.MedianWindow(S, imzMax, 72, false);
            for (int i=0 ; i<imzMax ; i++)
                S[i] = Math.max(0, S[i]-background[s][i]);
        }
        if (spectra.length == 1)
            return background;

        float[] bg=null, row=null;
        for (int i=0 ; i<imzMax ; i++)
        {
            row = Spectrum.getRow(spectra,i,row);
            bg = Spectrum.MinimaWindow(row, row.length, 15, bg);
            //bg = Spectrum.MedianWindow(row, row.length, 10, false);
            for (int s=0 ; s<spectra.length ; s++)
            {
                float b = bg[s];
                background[s][i] += b;
                spectra[s][i] = Math.max(0, spectra[s][i]-b);
            }
        }
        return background;
    }



    // not in Spectrum.java since this is experimental
    public static void threshold(double[][] Xin, double[] threshold)
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
                Arrays.fill(s, 0, s.length, 0.0);
                continue;
            }
            for (int i = 0; i < N; i++)
            {
                double d = s[i];
                double dm = Math.abs(d);
                s[i] = dm <= t ? 0.0 : d * Math.max(0, 1 - Math.exp(1 - dm / t));
            }
        }
    }


    public static double[] PadToDouble(float[] x, int pad)
    {
        return PadToDouble(x, x.length,  pad, null);
    }


    // since we're copying anyway, we might as well pad with something
    public static double[] PadToDouble(float[] x, int length, int pad, double[] y)
    {
        int N = length + 2 * pad;
        y = Transform.realloc(y, N);

        double first = x[0];
        double last = x[length-1];
        double middle = (first+last)/2;
        double d = (last-first)/(pad*2);

        // pad
        for (int i = 0; i < pad; i++)
            y[i] = middle - i * d;
        d = x[length - 1];
        for (int i = 0; i < pad; i++)
            y[y.length - i - 1] = middle + i * d;

        // copy/convert
        for (int i = 0; i < length; i++)
            y[i + pad] = x[i];
        return y;
    }


    public static float[] UnpadToFloat(double[] y, int pad, float[] x)
    {
        int N = y.length - 2 * pad;
        if (null == x || x.length != N)
            x = new float[N];
        for (int i = 0; i < N; i++)
            x[i] = (float) y[i + pad];
        return x;
    }


    static int pow2(int k)
    {
        return 1 << k;
    }


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


    public static float[] getRow(float[][] m, int r, float[] out)
    {
        out = realloc(out, m.length);
        for (int s = 0; s < m.length; s++)
            out[s] = m[s][r];
        return out;
    }


    public static void setRow(float[][] m, int r, float[] in)
    {
        for (int s = 0; s < m.length; s++)
            m[s][r] = in[s];
    }


    public static void main(String[] args)
    {
        float[] g10 = NormalP(11, .01, false, false, null);
        float[] g20 = NormalP(21, .01, false, false, null);
        for (int i=0 ; i<g10.length ; i++)
            System.out.println("" + (((double)i)/(g10.length-1)) + "\t" + g10[i]);
        System.out.println();
        for (int i=0 ; i<g20.length ; i++)
            System.out.println("" + (((double)i)/(g20.length-1)) + "\t" + g20[i]);

        assert 2F == Median(1, 2, 3);
        assert 2F == Median(1, 3, 2);
        assert 2F == Median(2, 1, 2);
        assert 2F == Median(2, 1, 3);
        assert 2F == Median(2, 3, 1);
        assert 2F == Median(3, 2, 1);
    }
}

