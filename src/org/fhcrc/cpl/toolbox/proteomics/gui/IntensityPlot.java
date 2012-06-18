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
package org.fhcrc.cpl.toolbox.proteomics.gui;

import org.fhcrc.cpl.toolbox.datastructure.FloatArray;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.*;
import java.util.Arrays;

/**
 */

public class IntensityPlot implements java.io.Serializable
{
    public static final String[] COLOR_SCHEMES =
            {"Cyan", "Grayscale", "Fire", "Rainbow", "Fancy"};

    FloatArray x = new FloatArray();
    FloatArray y = new FloatArray();
    FloatArray z = new FloatArray();




    public void setData(FloatArray x, FloatArray y, FloatArray z)
    {
        if (x.size() != y.size() || x.size() != z.size())
            throw new IllegalArgumentException();
        this.x = x;
        this.y = y;
        this.z = z;
    }


    int parse(Reader in)
            throws IOException
    {
        StreamTokenizer tokenizer = new StreamTokenizer(in);
        tokenizer.parseNumbers();
        int tt;
        int c = 0;
        while (StreamTokenizer.TT_EOF != (tt = tokenizer.nextToken()))
        {
            if (StreamTokenizer.TT_NUMBER != tt)
                continue;
            int m = c % 3;
            c++;
            switch (m)
            {
                case 0:
                    x.add((float)tokenizer.nval);
                    break;
                case 1:
                    y.add((float)tokenizer.nval);
                    break;
                case 2:
                    z.add((float)tokenizer.nval);
                    break;
            }
        }
        x.setSize(z.size());
        y.setSize(z.size());
        return x.size();
    }


    private final int bound(int a, int n, int x)
    {
        return a > x ? x : a < n ? n : a;
    }


    private final float bound(float x)
    {
        return x > 1.0 ? 1.0F : x < 0.0 ? 0.0F : x;
    }


    public void plotLog(BufferedImage image, float threshold, int yscale, String scheme)
    {
        ColorMap colorMap = mapForScheme(scheme);
        int type = image.getType();

//dhmay removing these checks along with adding support for blingy main image colors
//		if (type != BufferedImage.TYPE_USHORT_GRAY &&
//		        type != BufferedImage.TYPE_BYTE_GRAY &&
//		        type != BufferedImage.TYPE_3BYTE_BGR)
//			throw new IllegalArgumentException();

        int width = image.getWidth();
        int height = image.getHeight();
        WritableRaster r = image.getRaster();
        int maxColor = image.getType() == BufferedImage.TYPE_USHORT_GRAY ? (1 << 16) - 1 : (1 << 8) - 1;
        float nullColor = maxColor + 1.0F;

        FloatRange rangeZ = z.getRange();
        threshold = Math.max(Math.max(1, threshold), rangeZ.min);
        rangeZ.max = Math.max(threshold + 1, rangeZ.max);

        double logMaxIntensity = Math.log(rangeZ.max);
        double logMinIntensity = Math.log(threshold);
        double scale1 = 1.0 / (logMaxIntensity - logMinIntensity);
        double scale2 = 1.0 / Math.log(1 + rangeZ.max - threshold);

        float gray;
        int length = x.size();
        int xcurrent = width;
        float[] samples = new float[height];
        for (int i = 0; i < length; i++)
        {
            float intensity = z.get(i);
            int xval = (int)x.get(i);
            int yval = (int)Math.round(y.get(i) * yscale);
            if (yval < 0 || yval >= height)
                continue;

            if (intensity <= threshold)
                gray = 0.0F;
            else
            {
                if (false)
                    // too linear as threshold increases
                    gray = (float)((Math.log((double)intensity) - logMinIntensity) * scale1);
                else
                    // grainer, but more detail
                    gray = (float)(Math.log(1 + intensity - threshold) * scale2);
                gray = bound(gray);
            }

            if (xval != xcurrent)
            {
                if (xcurrent < width)
                {
                    if (yscale > 1)
                        interpolate(samples, maxColor);
                    plotColumn(r, xcurrent, height, samples, colorMap, maxColor);
                }
                Arrays.fill(samples, yscale > 1 ? nullColor : maxColor);
                xcurrent = xval;
            }
            int y = height - 1 - yval;
            samples[y] = Math.min(samples[y], maxColor - gray * maxColor);
        }
        if (xcurrent < width)
        {
            if (yscale > 1)
                interpolate(samples, maxColor);
            plotColumn(r, xcurrent, height, samples, colorMap, maxColor);
        }
    }


    private void plotColumn(java.awt.image.WritableRaster r,
                            int xcurrent, int height,
                            float[] samples, ColorMap scheme,
                            int maxColor)
    {
        int b = r.getNumBands();
        float colors[][] = scheme.remap(samples, maxColor);

        for (int i=0 ; i<b && i<colors.length ; i++)
            r.setSamples(xcurrent, 0, 1, height, i, colors[i]);
    }


    public void plotSqrt(BufferedImage image, float threshold, int yscale, String scheme)
    {
        ColorMap colorMap = mapForScheme(scheme);
        int type = image.getType();
        if (type != BufferedImage.TYPE_USHORT_GRAY &&
                type != BufferedImage.TYPE_BYTE_GRAY &&
                type != BufferedImage.TYPE_3BYTE_BGR)
            throw new IllegalArgumentException();

        int width = image.getWidth();
        int height = image.getHeight();
        WritableRaster r = image.getRaster();
        int maxColor = image.getType() == BufferedImage.TYPE_USHORT_GRAY ? (1 << 16) - 1 : (1 << 8) - 1;
        float nullColor = maxColor + 1.0F;

        FloatRange rangeZ = z.getRange();
        threshold = Math.max(Math.max(1, threshold), rangeZ.min);
        rangeZ.max = Math.max(threshold + 1, rangeZ.max);

//		double logMaxIntensity = Math.sqrt(rangeZ.max);
//		double logMinIntensity = Math.sqrt(threshold);
        double scale2 = 1.0 / Math.sqrt(rangeZ.max);

        float gray;
        int length = x.size();
        int xcurrent = width;
        float[] samples = new float[height];
        for (int i = 0; i < length; i++)
        {
            float intensity = z.get(i);
            int xval = (int)x.get(i);
            int yval = (int)Math.round(y.get(i) * yscale);
            if (yval < 0 || yval >= height)
                continue;

            if (intensity <= threshold)
                gray = 0.0F;
            else
            {
                gray = (float)(Math.sqrt(intensity) * scale2);
                gray = bound(gray);
            }

            if (xval != xcurrent)
            {
                if (xcurrent < width)
                {
                    if (yscale > 1)
                        interpolate(samples, maxColor);
                    plotColumn(r, xcurrent, height, samples, colorMap, maxColor);
                }
                Arrays.fill(samples, yscale > 1 ? nullColor : maxColor);
                xcurrent = xval;
            }
            int y = height - 1 - yval;
            samples[y] = Math.min(samples[y], maxColor - gray * maxColor);
        }
        if (xcurrent < width)
        {
            if (yscale > 1)
                interpolate(samples, maxColor);
            plotColumn(r, xcurrent, height, samples, colorMap, maxColor);
        }
    }


    private void interpolate(float[] samples, int maxColor)
    {
        float nullColor = maxColor + 1;
        int height = samples.length;
        float p = nullColor;
        float c = samples[0];
        float n;
        for (int j = 0; j < height; j++)
        {
            n = j < height - 1 ? samples[j + 1] : nullColor;
            if (nullColor == c)
            {
                if (nullColor != p && nullColor != n)
                    samples[j] = (p + n) / 2.0F;
                else
                    samples[j] = maxColor;
            }
            p = c;
            c = n;
        }
    }


    public void plot(BufferedImage image, float threshold)
    {
        int width = image.getWidth();
        int height = image.getHeight();
        WritableRaster r = image.getRaster();
        int maxColor =
                image.getType() == BufferedImage.TYPE_BYTE_GRAY ? 0x000000ff : 0x0000ffff;

        FloatRange rangeZ = z.getRange();

        double scaleIntensity = maxColor / (rangeZ.max - threshold);

        int gray;
        int length = x.size();
        for (int i = 0; i < length; i++)
        {
            float intensity = z.get(i);
            if (intensity < threshold)
                continue;
            gray = (int)((intensity - threshold) * scaleIntensity);
            gray = bound(gray, 0, maxColor);

            int xval = (int)x.get(i);
            int yval = (int)y.get(i);

            if (xval < width && yval < height)
            {
//                r.setSample(xval, height - yval - 1, 0, maxColor - gray);

                r.setSample(xval, height - yval - 1, 0,  0);

            }
        }

        Graphics gfx = image.getGraphics();
        gfx.setColor(Color.BLACK);
    }


    public BufferedImage plot(float threshold, boolean log, String colorScheme)
    {
        FloatRange rangeX = x.getRange();
        FloatRange rangeY = y.getRange();
        int width = (int)rangeX.max + 1;
        int height = (int)rangeY.max + 1;

        // TYPE_USHORT_GRAY does not seem to save reliably!
        int type = BufferedImage.TYPE_INT_RGB;
        BufferedImage image = new BufferedImage(width, height, type);
        Graphics gfx = image.getGraphics();
        gfx.setColor(Color.WHITE);
        gfx.fillRect(0, 0, width, height);

        if (log)
        {
            plotLog(image, threshold, 1, colorScheme);
        }
        else
            plot(image, threshold);
        return image;
    }


    public void plot(OutputStream out, String format, boolean log, String colorScheme)
            throws IOException
    {
        BufferedImage image = plot(100, log, colorScheme);
        writePlot(out, image, format);
    }


    public static void writePlot(OutputStream out, Image image, String format) throws IOException
    {
        if (!(image instanceof BufferedImage))
            throw new java.lang.UnsupportedOperationException();

        ImageIO.write((BufferedImage)image, format, out);
    }


    public static void writePlot(File f, Image image, String format)
            throws IOException
    {
//		if (!f.exists())
//			f.createNewFile();
//		else if (!f.isFile())
//			throw new IOException("Not a file: " + f.getPath());

        if (!(image instanceof BufferedImage))
            throw new java.lang.UnsupportedOperationException();

        ImageIO.write((BufferedImage)image, format, f);
    }


//    public static void main(String[] args)
//            throws java.io.FileNotFoundException, java.io.IOException
//    {
//        String pathIn = null, pathOut = null;
//        boolean log = false;
//
//        if (args.length < 2 || args.length > 3)
//        {
//            usage();
//        }
//        for (int i = 0; i < args.length; i++)
//        {
//            if ("-log".equals(args[i]))
//                log = true;
//            else if (pathIn == null)
//                pathIn = args[i];
//            else
//                pathOut = args[i];
//        }
//        if (pathIn == null || pathOut == null)
//            usage();
//
//        run(pathIn, pathOut, log);
//    }


    private static void usage()
    {
        System.err.println("usage: java org.fhcrc.cpl.toolbox.proteomics.gui.IntensityPlot [-log] input.tsv output.[tiff|png]");
        System.exit(1);
    }


    private static void run(String pathIn, String pathOut, boolean log, String colorScheme)
            throws IOException
    {
        long start = System.currentTimeMillis();
        IntensityPlot plot = new IntensityPlot();

        Reader in = new FileReader(pathIn);
        int rows = plot.parse(new FileReader(pathIn));
        System.out.println("read " + rows + " rows.");
        in.close();

        OutputStream out = new FileOutputStream(pathOut);
        String ext = pathOut.substring(pathOut.indexOf('.') + 1);
        plot.plot(out, ext, log, colorScheme);
        out.close();

        File f = new File(pathOut);
        System.out.println("wrote " + f.getAbsolutePath() + ", " + f.length() + " bytes.");
        System.out.println("" + ((System.currentTimeMillis() - start) / 1000.0) + " seconds.");
    }


    /**
     * for performance we use instances of color maps (rather than singletons)
     * we cache the output array for reuse
     */
    public static abstract class ColorMap
    {
        private float[] remap(float[] x, float[] out, float[] y, float range)
        {
            out = Spectrum.realloc(out, x.length);
            float z, dz;
            for (int i = 0; i < x.length; i++)
            {
                z = x[i] / range * ((float)(y.length - 1));
                int b = (int)(z);
                if (b > y.length - 2) b = y.length - 2;
                if (b < 0) b = 0;
                dz = Math.max(Math.min(z - (float)b, 1), 0);
                out[i] = range * ((1 - dz) * y[b] + dz * y[b + 1]);
            }
            return out;
        }


        public float[][] remap(float[] x, float max)
        {
            if (null == out)
                out = new float[3][x.length];
            out[0] = remap(x, out[0], _rmap, max);
            out[1] = remap(x, out[1], _gmap, max);
            out[2] = remap(x, out[2], _bmap, max);
            return out;
        }

        float[] _rmap;
        float[] _gmap;
        float[] _bmap;

        float[][] out;
    }


    public static class RainbowMap extends ColorMap
    {
        public RainbowMap()
        {
            _rmap   = new float[] {1, 1, 0, 0, 0, 1, 1}; //red
            _gmap = new float[] {0, 1, 1, 1, 0, 0, 1}; //green
            _bmap  = new float[] {0, 0, 0, 1, 1, 1, 1}; //blue
        }
    }

    public static class FireMap extends ColorMap
    {
        public FireMap()
        {
            _rmap   = new float[] {1,    1, 1, 1}; //red
            _gmap = new float[] {0, 0.5f, 1, 1}; //green
            _bmap  = new float[] {0,    0, 0, 1}; //blue
        }
    }

    public static class FancyMap extends ColorMap
    {
        public FancyMap()
        {
            _rmap   = new float[] {.5f,    .2f, .625f, .93f}; //red
            _gmap = new float[] {.1f, 0.2f, .8f, .95f}; //green
            _bmap  = new float[] {0.2f,    .65f, .7f, 1}; //blue
        }
    }

    public static class InvFireMap extends ColorMap
    {
        public InvFireMap()
        {
            _rmap = new float[] {1, 1, 1,     1,    1,     1}; //red
            _gmap = new float[] {1, 1, 0.75f, 0.5f, 0.25f, 0}; //green
            _bmap = new float[] {1, 0, 0,     0,    0,     0}; //blue
        }
    }

    public static class GrayMap extends ColorMap
    {
        public float[][] remap(float[] x, float max)
        {
            return new float[][] {x, x, x};
        }
    }


    public static class CyanMap extends ColorMap
    {
        public float[][] remap(float[] x, float max)
        {
            return new float[][] {x};
        }
    }

    public static ColorMap mapForScheme(String name)
    {
        if (name.equals("Rainbow"))
            return new RainbowMap();
        else if (name.equals("Fire"))
            return new FireMap();
        else if (name.equals("Cyan"))
            return new CyanMap();
        else if (name.equals("Heat"))
            return new InvFireMap();
        else if (name.equals("Fancy"))
            return new FancyMap();
        else
            return new GrayMap();
    }
}
