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
import org.fhcrc.cpl.toolbox.proteomics.gui.IntensityPlot;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.viewer.feature.*;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.BaseFeatureStrategy;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.toolbox.gui.ScrollableImage;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.gui.SurfaceFrame;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.datastructure.*;
import org.fhcrc.cpl.toolbox.proteomics.Scan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.swixml.SwingEngine;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;


/**
 * User: mbellew
 * Date: May 24, 2004
 * Time: 4:30:30 PM
 */
public class MSDetailPanel extends ScrollableImage
{
    Logger _log = Logger.getLogger(MSDetailPanel.class);

    static final String LOW_THRESHOLD = "MSDetailPanel.lowThreshold";
    static final String SCALE = "MSDetailPanel.scale";
    static final String COLORSCHEME = "MSDetailPanel.colorscheme";
    static final String PREPROCESS = "MSDetailPanel.preprocess";
    static final String FINDFEATURES = "MSDetailPanel.findfeatures";
    static final String SHOWCHARTLINE = "MSDetailPanel.showChartLine";

    static Font _smallFont = Font.decode("Verdana PLAIN 9");

    //hardcoded initial dimensions
    protected Dimension mPanelDimension = new Dimension(43,2048);

    Tree2D _featuresVisible = null;
    java.util.List _extractedFeatures = null;



    class ImageParameters implements Cloneable
    {
        FloatRange _mzRange = new FloatRange();
        //		float _mzMin;
        int _yScale;
        int _scanMin;
        int _scanCount;
        Rectangle _rectZoom = new Rectangle();

        protected Object clone()
        {
            try
            {
                ImageParameters c = (ImageParameters)super.clone();
                c._rectZoom = (Rectangle)_rectZoom.clone();
                return c;
            }
            catch (CloneNotSupportedException x)
            {
                return null;
            }
        }
    };
    ImageParameters _curr = new ImageParameters();

    Thread _renderThread = null;

    ListenerHelper helper = new ListenerHelper(this);


    public MSDetailPanel()
    {
        super("MSDetailPanel");
        setToolTipText("I'm a tool tip");

        Application app = Application.getInstance();
        helper.addListener(app, "repaint_propertyChange", MSDetailPanel.SHOWCHARTLINE);
        helper.addListener(app, "repaint_propertyChange", SharedProperties.CHART_RANGE);
        helper.addListener(app, "repaint_propertyChange", SharedProperties.FEATURE_RANGES);
        helper.addListener(app, "repaint_propertyChange", "featureSelector");

        helper.addListener(app, "update_propertyChange", MSDetailPanel.FINDFEATURES);
        helper.addListener(app, "update_propertyChange", SharedProperties.MS_RUN);
        helper.addListener(app, "update_propertyChange", SharedProperties.SELECTED_POINT);
        helper.addListener(app, "update_propertyChange", FeatureExtractor.DEFAULT_EXTRACTOR_PROPERTYNAME);
        helper.addListener(app, "update_propertyChange", LOW_THRESHOLD);
        helper.addListener(app, "update_propertyChange", SCALE);
        helper.addListener(app, "update_propertyChange", COLORSCHEME);
        helper.addListener(app, "update_propertyChange", PREPROCESS);

        helper.addListener(this, "_componentResized");
        helper.addListener(this, "_mouseClicked");
    }


    static boolean scheduledUpdate = false;

    public void update_propertyChange(PropertyChangeEvent event)
    {
        if (scheduledUpdate)
            return;
        scheduledUpdate = true;

        EventQueue.invokeLater(new Runnable()
        {
            public void run()
            {
                scheduledUpdate = false;
                updateImage();
            }
        });
    }


    public void repaint_propertyChange(PropertyChangeEvent event)
    {
        if (getFindFeatures() && SharedProperties.CHART_RANGE.equals(event.getPropertyName()))
        {
            Spectrum.Peak p = (Spectrum.Peak)ApplicationContext.getProperty(SharedProperties.SELECTED_POINT);
            if (null != getImage() && null == _renderThread && null != p)
            {
                // Is it scan at a time??, then update features
                // this is hacky
                Class c = (Class) ApplicationContext.getProperty(FeatureExtractor.DEFAULT_EXTRACTOR_PROPERTYNAME);
                if (c == FeatureStrategyWavelet.class || c == FeatureStrategyCentroided.class)
                {
                    MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
                    Pair chartRange = (Pair)event.getNewValue();
                    if (null != chartRange && ((Spectrum.Peak)chartRange.first).scan == ((Spectrum.Peak)chartRange.second).scan)
                    {
                        int scan = ((Spectrum.Peak)chartRange.first).scan;
                        int scanIndex = run.getIndexForScanNum(scan);
                        if (scanIndex < 0) scanIndex = -(scanIndex+1);
                        ImageGenerator gen = new ImageGenerator(run, scanIndex, p.mz);
                        gen._findDetailFeatures(_curr);
                    }
                }
            }
        }

        repaint();
    }


    public void _componentResized(ComponentEvent event)
    {
        updateImage();
    }


    protected MSRun getRun()
    {
        return (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
    }

    public void _mouseClicked(MouseEvent event)
    {
        if (null == _curr || 0 == _curr._yScale)
            return;

        MSRun run = getRun();

        Feature fBest = hitTest(event);
        if (fBest != null)
        {
            ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, fBest);
            ApplicationContext.setProperty(SharedProperties.SELECTED, fBest);
            if (null != run)
            {
                int n = run.getIndexForScanNum(fBest.scan, true);
                ApplicationContext.setProperty(SharedProperties.MS_SCAN, run.getScan(n));
            }
        }
        else
        {
            if (null == run)
                return;

            int x = event.getX();
            int y = getHeight() - event.getY() - 1;

            if (x + _curr._scanMin >= 1 && x + _curr._scanMin <= run.getScanCount())
            {
                MSRun.MSScan scan = run.getScan(x + _curr._scanMin);
                ApplicationContext.setProperty(SharedProperties.MS_SCAN, scan);
                float mz = y / (float)_curr._yScale + _curr._mzRange.min;
//				Spectrum.Peak selectedPoint = (Spectrum.Peak)app.getProperty(SharedProperties.SELECTED_POINT);
//				if (null != selectedPoint)
//					app.setProperty(SharedProperties.SELECTED_POINT, new Spectrum.Peak(x + _scanMin, selectedPoint.mz, 0F));
//				else
                ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, new Spectrum.Peak(scan.getNum(), mz, 0F));
            }
        }
    }


    public String getToolTipText(MouseEvent event)
    {
        if (null == _curr || 0 == _curr._yScale)
            return "";

        int height = getHeight();
        int y = height - event.getY() - 1;
        int x = event.getX();

        MSRun run = getRun();
        if (null == run)
            return "";

        //check the scan under the mouse event to make sure it's within range for the run
        int scanIndex = x + _curr._scanMin;
        if (scanIndex >= run.getScanCount())
            return "";

        MSRun.MSScan scan = run.getScan(scanIndex);
        float mz = (float)y / _curr._yScale + _curr._mzRange.min;
        String text = "" + scan.getNum() + ", " + mz;

        Feature fBest = hitTest(event);
        if (null != fBest)
        {
            StringBuffer sb = new StringBuffer("<html>" +
                    fBest.scan + "<br>" + fBest.mz + "mz (" + fBest.getMass() + ")<br>" + fBest.charge + "+<br>" +
                    ((int)fBest.intensity) + "i" +
                    (0 == fBest.background ? "" : " / " + Math.round(fBest.background) + "bg") +
                    (0 == fBest.getMedian() ? "" : " / " + Math.round(fBest.getMedian()) + "md"));
            if (fBest instanceof Feature)
            {
                Feature fr = (Feature)fBest;
                if (fr.getScanFirst() != fr.getScanLast())
                    sb.append("<br>" + fr.getScanFirst() + "-" + fr.getScanLast());
            }
            if (null != fBest.getDescription())
                sb.append("<br>" + fBest.getDescription());
            sb.append("</html>");
            text = sb.toString();
        }
        return text;
    }


    private Feature hitTest(MouseEvent event)
    {
        if (null == _featuresVisible)
            return null;

        int height = getHeight();
        int y = height - event.getY() - 1;
        int x = event.getX();
        int scanClicked = _curr._scanMin + x;
        MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);

        float mz = (float)y / _curr._yScale + _curr._mzRange.min;

        Feature fBest = null;
        double dist = 0.2;
        ArrayList near = _featuresVisible.getPoints(scanClicked - 3, mz - 1, scanClicked + 3, mz + 1);
        for (Iterator it = near.iterator(); it.hasNext();)
        {
            Feature f = (Feature)it.next();
            double d = Math.abs(mz - f.mz);
            if (d < dist)
            {
                int scanIndex = run.getIndexForScanNum(f.scan, true);
                if (Math.abs(scanIndex - scanClicked) <= 3)
                {
                    dist = d;
                    fBest = f;
                }
            }
        }
        return fBest;
    }


    public Dimension getPreferredSize()
    {
        return super.getPreferredSize();
    }


    protected synchronized void updateImage()
    {
        this._featuresVisible = null;
        this.setImage(null, null);

        Spectrum.Peak p = (Spectrum.Peak)ApplicationContext.getProperty(SharedProperties.SELECTED_POINT);
        if (null == p)
            return;
        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        if (null == run)
            return;

        ApplicationContext.setProperty(SharedProperties.ZOOM_REGION, null);

        if (null != _renderThread && _renderThread.isAlive())
        {
            _renderThread.interrupt();
            _renderThread = null;
        }
        int index = run.getIndexForScanNum(p.scan, true);
        _renderThread = new Thread(new ImageGenerator(run, index, (double)p.mz));
        _renderThread.setPriority(3);
        _renderThread.start();
    }


    public void setImage(Image image, ImageParameters p)
    {
        _curr = p;
        super.setImage(image);
    }


    static Object _imageGeneratorRunLock = new Object();

    private class ImageGenerator implements Runnable
    {
        MSRun _run;
        int _scanIndex;
        double _mz;
        int _lowThreshold;
        boolean _process = false;
        boolean _findFeatures = false;
        String _scale = "log";


        public ImageGenerator(MSRun run, int scanIndex, double mz)
        {
            _run = run;
            _scanIndex = scanIndex;
            _mz = mz;
            _lowThreshold = getLowThreshold();
            _scale = getScale();
            _process = getPreprocess();
            _findFeatures = getFindFeatures();
        }


        public void run()
        {
            try
            {
                final ImageParameters newImage = new ImageParameters();

                // need to synchronize this or we will trample member variables
                synchronized (_imageGeneratorRunLock)
                {
                    _generateDetailImage(newImage);
                    if (_findFeatures)
                        _findDetailFeatures(newImage);
                }
            }
            finally
            {
                synchronized(MSDetailPanel.this)
                {
                    if (_renderThread == Thread.currentThread())
                        _renderThread = null;
                }
            }
        }


        private void _generateDetailImage(ImageParameters newImage) // Rectangle rectZoom /* OUT */)
        {
            _log.debug("_generateDetailImage() thread " + Thread.currentThread());
            Thread me = Thread.currentThread();

            boolean centroided = _run.getHeaderInfo().getDataProcessing().getCentroided() == 1;
            int width = getWidth();
            newImage._scanCount = width - 3;
            int height = getHeight();

            // want range of about +/-10 but keep yscale < 30
            float mzRange = 20;
            newImage._yScale = (int)Math.min(height / mzRange, 36);
            mzRange = height / newImage._yScale;
            newImage._scanMin = Math.max(0, _scanIndex - newImage._scanCount / 2);
            int scanMax = Math.min(newImage._scanMin + newImage._scanCount, _run.getScanCount());
            newImage._scanCount = scanMax - newImage._scanMin;
            newImage._mzRange.min = (int)Math.floor(_mz - (mzRange / 2));
            newImage._mzRange.max = newImage._mzRange.min + mzRange + 1;

            Rectangle rect = new Rectangle(Math.max(10,width), Math.max(10,height));
            BufferedImage img = new BufferedImage((int)rect.getWidth(), (int)rect.getHeight(), BufferedImage.TYPE_3BYTE_BGR);
            Graphics2D g = (Graphics2D)img.getGraphics();
            g.setFont(_smallFont);
            g.setBackground(Color.WHITE);
            g.clearRect((int)rect.getX(), (int)rect.getY(), (int)rect.getWidth(), (int)rect.getHeight());


            // collect all the points in +/- 10 scans and within +/- 10D
            FloatArray x = new FloatArray();
            FloatArray y = new FloatArray();
            FloatArray z = new FloatArray();

            newImage._rectZoom.x = newImage._scanMin;
            newImage._rectZoom.y = (int)Math.ceil(newImage._mzRange.min);
            newImage._rectZoom.width = newImage._scanCount;
            newImage._rectZoom.height = (int)Math.round(newImage._mzRange.max - newImage._mzRange.min);

            FloatRange mzWindow = new FloatRange(newImage._mzRange.min - 2, newImage._mzRange.max + 2);

            MSRun.MSScan[] scans = new MSRun.MSScan[scanMax - newImage._scanMin];
            for (int s = 0; s < scanMax - newImage._scanMin; s++)
                scans[s] = _run.getScan(newImage._scanMin + s);

            float[][] spectra = new float[scans.length][];
            for (int s = 0; s < scans.length; s++)
            {
                if (me.isInterrupted())
                    return;
                float[][] spectrum = scans[s].getSpectrum();
                if (_process && centroided)
                    spectrum = FeatureStrategyCentroided.cleanSpectrum(spectrum);
                float[] signal = Spectrum.Resample(spectrum, mzWindow, 36);
                spectra[s] = signal;
            }

            if (_process)
            {
                Spectrum.RemoveBackground(spectra);
                if (me.isInterrupted())
                    return;
            }

            for (int s = 0; s < scans.length; s++)
            {
                float[] signal = spectra[s];
                float base = mzWindow.min - newImage._mzRange.min;
                for (int i = 0; i < signal.length; i++)
                {
                    float mzOff = base + i / 36F;
                    float in = signal[i];
                    if (mzOff < 0 || in == 0) continue;
                    x.add((float)s);
                    y.add(mzOff);
                    z.add(in);
                }
            }


            IntensityPlot plot = new IntensityPlot();
            String scheme = MSDetailPanel.getColorScheme();
            plot.setData(x, y, z);
            if ("sqrt".equals(_scale))
                plot.plotSqrt(img, _lowThreshold, newImage._yScale, scheme);
            else
                plot.plotLog(img, _lowThreshold, newImage._yScale, scheme);

            g.setColor(Color.BLACK);
            for (int w = (int)rect.getWidth(), m = (int)newImage._mzRange.min; m <= newImage._mzRange.max; m++)
            {
                int t = height - 1 - (m - (int)newImage._mzRange.min) * newImage._yScale;
                if (0 == (m % 10))
                {
                    //g.drawLine(w - 2, t + 1, w, t + 1);
                    //g.drawLine(w - 6, t, w, t);
                    //g.drawLine(w - 2, t - 1, w, t - 1);
                    int stringWidth = g.getFontMetrics().stringWidth(String.valueOf(m));
                    g.drawLine(w - 2, t, w, t);
                    g.drawString(String.valueOf(m), w - stringWidth - 2, t + 3);
                }
                else if (0 == (m % 5))
                {
                    g.drawLine(w - 2, t + 1, w, t + 1);
                    g.drawLine(w - 5, t, w, t);
                    g.drawLine(w - 2, t - 1, w, t - 1);
                }
                else
                    g.drawLine(w - 3, t, w, t);
            }
            if (me.isInterrupted())
                return;

            final ImageParameters p = (ImageParameters)newImage.clone();
            final BufferedImage i = img;
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    _extractedFeatures = null;
                    _featuresVisible = null;
                    ApplicationContext.setProperty(SharedProperties.ZOOM_REGION, p._rectZoom);
                    setImage(i, p);
                }
            });
        }


        private void _findDetailFeatures(ImageParameters newImage)
        {
            Thread me = Thread.currentThread();
            MSRun.MSScan[] scans;
            java.util.List featureSets = new java.util.LinkedList();
            FeatureSet extractFeatures = null; // the _main_ features

            int scanMax = Math.min(newImage._scanMin + newImage._scanCount, _run.getScanCount());

            // Show information related to the feature extraction algorithm
            long msStart = 0; // don't want to time finding background features
            try
            {
                ApplicationContext.setMessage("Looking for features. . .");
                int t = scanMax - newImage._scanMin + 32;
                int len = 128;
                for (; len < t; len = len * 2)
                    ;
                int start = Math.max(0, _scanIndex - len / 2);
                int end = Math.min(_run.getScanCount(), start + len);
                len = end - start;
                scans = new MSRun.MSScan[len];
                for (int s = 0; s < len; s++)
                    scans[s] = _run.getScan(start + s);

                Class c = (Class)ApplicationContext.getProperty(FeatureExtractor.DEFAULT_EXTRACTOR_PROPERTYNAME);
                if (-1 != c.getName().toLowerCase().indexOf("wavelet2d"))
                {
                    FeatureSet s;
                    Spectrum.Peak[] peaks = ExtractMaxima2D.analyze(scans,
                            new FloatRange(newImage._mzRange.min - 0.1F, newImage._mzRange.max + 0.1F),
                            new FeatureStrategyWavelet2D.SmoothWavelet(),
                            -Float.MAX_VALUE);
                    if (me.isInterrupted())
                        return;
                    s = new FeatureSet(peaks, new Color(0, 0, 0, 0x7f));
                    s.setStyle(2); // pixel
                    featureSets.add(s);
                }
                else if (-1 != c.getName().toLowerCase().indexOf("combined")) //  || -1 != c.getTextCode().toLowerCase().indexOf("clusters"))
                {
                    try
                    {
                        FeatureSet s = new FeatureStrategyPeaks(_run, newImage._scanMin, newImage._scanCount, 6, new FloatRange(newImage._mzRange.min - 10F, newImage._mzRange.max + 10F), 4.0).analyze();
                        s.setColor(new Color(0, 0, 0, 0x40));
                        s.setStyle(2); // pixel
                        featureSets.add(s);
                    }
                    catch (InterruptedException xx)
                    {
                        return;
                    }
                }
                else if (-1 != c.getName().toLowerCase().indexOf("gross"))
                {
                    msStart = System.currentTimeMillis(); // HACK: we're finding background and foreground features in parallel
                    FeatureSet[] s = ExtractEdgeFeatures.EdgeFeatures(scans, new FloatRange(newImage._mzRange.min - 10F, newImage._mzRange.max + 10F), 0.0F);
                    if (me.isInterrupted())
                        return;
                    s[0].setStyle(2); // pixel
                    s[1].setStyle(2); // pixel
                    featureSets.add(s[0]);
                    featureSets.add(s[1]);
                    if (-1 != c.getName().toLowerCase().indexOf("gross"))
                    {
                        // don't waste time using FeatureExtractor, we've done the work already
                        featureSets.add(s[2]);
                        extractFeatures = s[2];
                        extractFeatures.setStyle(1);
                        extractFeatures.setColor(Color.BLACK);
                    }
                }
                else if (-1 != c.getName().toLowerCase().indexOf("smooth"))
                {
                    Spectrum.Peak[] peaks = ExtractMaxima2D.analyze(scans, new FloatRange(newImage._mzRange.min - 10F, newImage._mzRange.max + 10F));
                    if (me.isInterrupted())
                        return;
                    FeatureSet s = new FeatureSet(peaks, Color.BLACK);
                    s.setStyle(2); // pixel
                    featureSets.add(s);
                }

                if (me.isInterrupted())
                    return;

                // features on selected scan(s)
                if (null == extractFeatures)
                {
                    if (0 == msStart)
                        msStart = System.currentTimeMillis();
                    FloatRange mzRangeExtract = new FloatRange(newImage._mzRange.min - 10.0F, newImage._mzRange.max + 10.0F);

                    Class featureStrategyClass = FeatureExtractor.getDefaultClass();
                    String featureStrategyPackageName = null;
                    //TODO:  this should NOT be necessary.  But something about the way one-jar loads
                    //the new-school classes is proving a problem.
                    try
                    {
                        featureStrategyPackageName = featureStrategyClass.getPackage().getName();
                    }
                    catch (NullPointerException e)
                    {
                        String className = featureStrategyClass.getName();
                        featureStrategyPackageName = className.substring(0, className.lastIndexOf("."));
                    }
                    
                    //determine whether this is an old-school feature strategy and invoke appropriately
                    if ("org.fhcrc.cpl.viewer.feature".equals(featureStrategyPackageName))
                    {
                        FeatureExtractor alg = FeatureExtractor.getDefault(_run, newImage._scanMin, newImage._scanCount, 6, mzRangeExtract, 2.0);
                        // UNDONE: need dialog to set maxCharge
                        if (FeatureExtractor.TYPE_1D == alg.getType())
                            alg = FeatureExtractor.getDefault(_run, _scanIndex, 1, 6, mzRangeExtract, 2.0);

                        extractFeatures = alg.analyze();

                    }
                    else
                    {

                        FeatureFinder featureFinder =
                                new FeatureFinder(_run, newImage._scanMin, newImage._scanCount,
                                        PeakCombiner.DEFAULT_MAX_CHARGE,
                                        mzRangeExtract,
                                        featureStrategyClass, false);

                        featureFinder.setStatusListener(new BaseFeatureStrategy.StatusListener()
                        {
                            public void progress(float percent)
                            {
                                float p = Math.round(percent * 10) / 10F;
                                ApplicationContext.setMessage(TextProvider.getText("FINDING_ALL_FEATURES_IN_FILE","FILEPATH",_run.getFile().getName()) +
                                        TextProvider.getText("PERCENT_PERCENT_COMPLETE_DOT","PERCENT",String.valueOf(p)));
                            }
                        });

                        extractFeatures = featureFinder.findPeptides();
                    }


                }
                if (extractFeatures != null)
                {
                    extractFeatures.setColor(Color.BLACK);
                    if (me.isInterrupted())
                        return;
                    featureSets.add(extractFeatures);
//							if (-1 != alg.getClass().getTextCode().indexOf("FeatureStrategyPeaks"))
//								extractFeatures.setStyle(2); // pixels
//							else
                    extractFeatures.setStyle(4); // 1=X 4=box
                }

                if (me.isInterrupted())
                    return;

                final java.util.List sets = featureSets;
                EventQueue.invokeLater(new Runnable()
                {
                    public void run()
                    {
                        _extractedFeatures = sets;
                        repaint();
                    }
                });

            }
            catch (InterruptedException ex)
            {
                return;
            }
            finally
            {
                long msEnd = System.currentTimeMillis();
                if (null == extractFeatures) // interrupted probably
                    ApplicationContext.setMessage(" ");
                else
                {
                    FeatureSet.FeatureSelector sel = getSelector(_run, newImage);
                    int c = extractFeatures.getFeatures(sel).length;
                    String time = NumberFormat.getInstance().format((msEnd - msStart) / 1000.0);
                    ApplicationContext.setMessage("Found " + c + " feature" + (c == 1 ? "" : "s") + " in " + time + "s");
                }
            }
        }
    }


    static Color TRANSLUCENT_RED = new Color(1.0F, 0.0F, 0.0F, 0.4F);
    static Color TRANSLUCENT_BLUE = new Color(0.0F, 0.0F, 1.0F, 0.6F);


    public void paint(Graphics g)
    {
        _log.debug("MSDetailPanel.paint()");
        super.paint(g);

        if (null == this._curr)
            return;

        Graphics2D graphics = (Graphics2D)g;

        //
        // draw SpectrumComponent range
        //

        Pair chartRange = (Pair)ApplicationContext.getProperty(SharedProperties.CHART_RANGE);
        if (getShowChartRange() && null != chartRange)
        {
            int height = getHeight();
            Spectrum.Peak start = (Spectrum.Peak)chartRange.first;
            Spectrum.Peak end = (Spectrum.Peak)chartRange.second;

            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null == run || null == _curr)
                return;
            int xStart = run.getIndexForScanNum(start.scan) - _curr._scanMin;
            int xEnd = run.getIndexForScanNum(end.scan) - _curr._scanMin;
            int yStart = (int)(height - 1 - (start.mz - _curr._mzRange.min) * _curr._yScale);
            int yEnd = (int)(height - 1 - (end.mz - _curr._mzRange.min) * _curr._yScale);
            if ("fire".equals(getColorScheme().toLowerCase()) || "rainbow".equals(getColorScheme().toLowerCase()))
                graphics.setColor(TRANSLUCENT_BLUE);
            else
                graphics.setColor(TRANSLUCENT_RED);
            graphics.drawLine(xStart, yStart, xEnd, yEnd);
        }


        //
        // draw features
        //

        // for scale get TIC of selected scan,
        // CONSIDER use median TIC of whole file so it doesn't change
        float ticScale = 1;
        try
        {
            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
            int scanMin = _curr._scanMin;
            float[] tic = new float[_curr._scanCount];
            int s = 0;
            for (s = 0; s < _curr._scanCount && scanMin + s < run.getScanCount(); s++)
                tic[s] = run.getScan(scanMin + s).getTotIonCurrent();
            float med = Spectrum.Median(tic, 0, s, false, null);
            ticScale = Math.max(ticScale, med / 200000F);
        }
        catch (Throwable t)
        {
        }


        // only need to recompute if features change, but this is convienent
        _featuresVisible = new Tree2D();

        MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
        if (null == run)
            return;
        ArrayList featureSetsApp = (ArrayList)ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
        ArrayList featureSets = new ArrayList();
        if (null != featureSetsApp)
            featureSets.addAll(featureSetsApp);
        if (null != _extractedFeatures)
            featureSets.addAll(_extractedFeatures);

        for (int i = 0; i < featureSets.size(); i++)
        {
            FeatureSet fs = (FeatureSet)featureSets.get(i);
            if (!fs.isDisplayed())
                continue;
            Color color = fs.getColor();
            Color light = new Color(color.getRed(), color.getGreen(), color.getBlue(), 0x40);

            int height = getHeight();
            float mzMin = _curr._mzRange.min;

            FeatureSet.FeatureSelector selector = getSelector(run, _curr);
            Feature[] features = fs.getFeatures(selector);
            graphics.setColor(color);
            for (int j = 0; j < features.length; j++)
            {
                Feature f = features[j];

                // UNDONE: convert to retention time
                int index = run.getIndexForScanNum(f.scan, true);

                if (fs.getStyle() != 2) // don't hit test the "pixel" features
                    _featuresVisible.add(index, f.mz, f);

                int yDraw = (int)(height - 1 - (f.mz - mzMin) * _curr._yScale);
                float mzTop = f.mz + Math.max(0, f.peaks - 1) / Math.max(1F, f.charge);
                int yTop = (int)(height - 1 - (mzTop - mzMin) * _curr._yScale);
                int xDraw = index - _curr._scanMin;

                int leftMargin = 2;
                int rightMargin = 2;
                if (f instanceof Feature)
                {
                    Feature fr = (Feature)f;
                    int indexStart = run.getIndexForScanNum(fr.getScanFirst(), true);
                    int indexLast = run.getIndexForScanNum(fr.getScanLast(), true);

                    leftMargin = fr.getScanFirst() == 0 ? 2 : Math.min(40, Math.max(2, index - indexStart));
                    rightMargin = fr.getScanLast() == 0 ? 2 : Math.min(40, Math.max(2, indexLast - index));
                }

                if (fs.getStyle() == 1 || fs.getStyle() == 4) // X
                {
                    int l = (int)Math.log(f.intensity / ticScale + 6);
                    g.drawLine(xDraw - l, yDraw - l, xDraw + l, yDraw + l);
                    g.drawLine(xDraw - l, yDraw + l, xDraw + l, yDraw - l);
                }
                else if (fs.getStyle() == 2) // pixel
                {
                    g.drawLine(xDraw, yDraw, xDraw, yDraw);
                }
                else // +
                {
                    if (null != IsotopicLabelExtraInfoDef.getLabel(f) &&
                            0 < IsotopicLabelExtraInfoDef.getLabelCount(f))
                    {
                        float mzHeavy =
                                f.mz +
                                        (IsotopicLabelExtraInfoDef.getLabelCount(f) *
                                                IsotopicLabelExtraInfoDef.getLabel(f).getHeavy()
                                                / f.charge);
                        int y2 = (int)(height - 1 - (mzHeavy - mzMin) * _curr._yScale);
                        graphics.drawLine(xDraw - 2, yDraw, xDraw + 2, yDraw);
                        graphics.drawLine(xDraw - 1, yDraw - 1, xDraw + 1, yDraw - 1);
                        graphics.drawLine(xDraw - 2, y2, xDraw + 2, y2);
                        graphics.drawLine(xDraw - 1, y2 + 1, xDraw + 1, y2 + 1);
                        graphics.drawLine(xDraw, yDraw, xDraw, y2);
                    }
                    else
                    {
                        g.drawLine(xDraw - 2, yDraw, xDraw + 2, yDraw);
                        g.drawLine(xDraw, yDraw - 1, xDraw, yDraw + 1);
                    }
                }
                if (fs.getStyle() == 4 && f instanceof Feature)
                {
                    Feature fr = (Feature)f;
                    if (fr.peaks < 2 && fr.scanFirst == fr.scanLast)
                        continue;

                    int bottom = yDraw + 2;
                    int top = yTop - 2;
                    if (fr.peaks < 2)
                        g.drawLine(xDraw - leftMargin, yDraw, xDraw + rightMargin, yDraw);
                    else if (fr.scanFirst == fr.scanLast)
                        g.drawLine(xDraw, top - 1, xDraw, bottom + 1);
                    else
                    {
                        g.setColor(light);
                        g.drawOval(xDraw - leftMargin - 1, top - 1, leftMargin + rightMargin + 2, bottom - top + 2);
                        g.setColor(color);
                    }
                }
            }
        }
    }


    private FeatureSet.FeatureSelector getSelector(MSRun run, ImageParameters img)
    {
        float mzMin = img._mzRange.min;
        float mzMax = mzMin + img._rectZoom.height;

        FeatureSet.FeatureSelector selector = (FeatureSet.FeatureSelector)ApplicationContext.getProperty("featureSelector");
        if (null == selector)
            selector = new FeatureSet.FeatureSelector();
        else
            selector = (FeatureSet.FeatureSelector)selector.clone();

        selector.setMinMz(mzMin - 30); // need a pretty wide range here because of icat
        selector.setMaxMz(mzMax + 1);
        // map scanIndex to scanNum
        selector.setScanFirst(run.getScan(Math.max(0, img._scanMin - 10)).getNum());
        selector.setScanLast(run.getScan(Math.min(run.getScanCount() - 1, img._scanMin + img._scanCount + 10)).getNum());
        return selector;
    }


    // CONSIDER: just use statics instead of application properties
    static int getLowThreshold()
    {
        try
        {
            String t = (String)ApplicationContext.getProperty(LOW_THRESHOLD);
            return Integer.parseInt(t);
        }
        catch (Exception x)
        {
            return 0;
        }
    }


    static void setLowThreshold(int low)
    {
        ApplicationContext.setProperty(LOW_THRESHOLD, "" + low);
    }


    static void setScale(String scale)
    {
        ApplicationContext.setProperty(SCALE, scale);
    }

    static void setColorScheme(String scheme)
    {
        ApplicationContext.setProperty(COLORSCHEME, scheme);
    }


    static void setPreprocess(boolean r)
    {
        ApplicationContext.setProperty(PREPROCESS, r ? Boolean.TRUE : Boolean.FALSE);
    }


    static boolean getPreprocess()
    {
        try
        {
            Boolean bool = (Boolean)ApplicationContext.getProperty(MSDetailPanel.PREPROCESS);
            return null == bool ? false : bool.booleanValue();
        }
        catch (Exception x)
        {
        }
        return false;
    }


    static void setFindFeatures(boolean r)
    {
        ApplicationContext.setProperty(FINDFEATURES, r ? Boolean.TRUE : Boolean.FALSE);
    }


    static boolean getFindFeatures()
    {
        try
        {
            Boolean bool = (Boolean)ApplicationContext.getProperty(MSDetailPanel.FINDFEATURES);
            return null == bool ? true : bool.booleanValue();
        }
        catch (Exception x)
        {
        }
        return true;
    }


    static void setShowChartRange(boolean r)
    {
        ApplicationContext.setProperty(SHOWCHARTLINE, r ? Boolean.TRUE : Boolean.FALSE);
    }


    static boolean getShowChartRange()
    {
        try
        {
            Boolean bool = (Boolean)ApplicationContext.getProperty(MSDetailPanel.SHOWCHARTLINE);
            return null == bool ? true : bool.booleanValue();
        }
        catch (Exception x)
        {
        }
        return true;
    }


    static String getScale()
    {
        String scale = null;
        try
        {
            scale = (String)ApplicationContext.getProperty(SCALE);
        }
        catch (Exception x)
        {
        }
        return null == scale ? "log" : scale;
    }

    static String getColorScheme()
    {
        String s = null;
        try
        {
            s = (String)ApplicationContext.getProperty(COLORSCHEME);
        }
        catch (Exception x)
        {
        }
        return null == s ? "Cyan" : s;
    }


    public static class ShowDialogAction extends AbstractAction
    {
        ShowDialogAction()
        {
            super("Detail Pane Settings");
        }

        public void actionPerformed(ActionEvent e)
        {
            new DetailDialog().setVisible(true);
        }
    }


    public static class Show3DAction extends AbstractAction
    {
        Show3DAction()
        {
            super("3D view of Detail Pane");
        }

        public void actionPerformed(ActionEvent e)
        {
            Rectangle rect = (Rectangle)ApplicationContext.getProperty(SharedProperties.ZOOM_REGION);
            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null == rect || null == run)
                return;

            ApplicationContext.setMessage(rect.toString());

            float[][] matrix = new float[rect.width][];
            FloatRange r = new FloatRange(rect.y, rect.y + rect.height);
            int start = Math.max(0, rect.x);
            int end = Math.min(run.getScanCount() - 1, rect.x + rect.width);
            for (int s = start; s < end; s++)
            {
                Scan scan = run.getScan(s);
                matrix[s - start] = Spectrum.Resample(scan.getSpectrum(), r, 36);
            }
            JFrame frame = SurfaceFrame.ShowSurfaceFrame(matrix);
            frame.setVisible(true);
        }
    }


    public static class DetailDialog extends JDialog
    {
        public JSpinner lowSpinner;
        public JButton okButton;
        public JButton cancelButton;
        public JButton applyButton;
        public JPanel contentPanel;
        public JRadioButton radioLogScale;
        public JRadioButton radioSqrtScale;
        public JCheckBox checkPreprocess;
        public JCheckBox checkDetectFeatures;
        public JComboBox heatmapColorScheme;
        public JCheckBox checkShowChartRange;

        DetailDialog()
        {
            super(ApplicationContext.getFrame(), "Detail Pane Settings", false);

            try
            {
                new SwingEngine(this).render("org/fhcrc/cpl/viewer/gui/DetailDialog.xml");
                assert null != contentPanel;
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage("error creating dialog", x);
                throw new RuntimeException(x);
            }

            this.setContentPane(contentPanel);
            ButtonGroup group = new ButtonGroup();
            group.add(radioSqrtScale);
            group.add(radioLogScale);

            String[] s = IntensityPlot.COLOR_SCHEMES;
            for (int j = 0; j < s.length; j++) heatmapColorScheme.addItem(s[j]);

            lowSpinner.setValue(new Integer(getLowThreshold()));
            ((SpinnerNumberModel)lowSpinner.getModel()).setMinimum(new Integer(0));

            if ("sqrt".equals(getScale()))
                radioSqrtScale.setSelected(true);
            else
                radioLogScale.setSelected(true);
            heatmapColorScheme.setSelectedItem(getColorScheme());

            checkPreprocess.setSelected(getPreprocess());
            checkDetectFeatures.setSelected(getFindFeatures());
            checkShowChartRange.setSelected(getShowChartRange());

            ListenerHelper helper = new ListenerHelper(this);
            helper.addListener(cancelButton,"cancel_actionPerformed");
            helper.addListener(okButton,"ok_actionPerformed");
            helper.addListener(applyButton,"apply_actionPerformed");

            this.setSize(300, 220);
            doLayout();
        }


        public void ok_actionPerformed(ActionEvent e)
        {
            apply_actionPerformed(null);
            this.dispose();
        }

        public void cancel_actionPerformed(ActionEvent e)
        {
            this.dispose();
        }

        public void apply_actionPerformed(ActionEvent e)
        {
            if ((lowSpinner.getValue() instanceof Integer))
            {
                Integer I = (Integer)lowSpinner.getValue();
                int low = Math.max(0, I.intValue());
                MSDetailPanel.setLowThreshold(low);
            }
            MSDetailPanel.setScale(radioSqrtScale.isSelected() ? "sqrt" : "log");
            MSDetailPanel.setColorScheme((String)heatmapColorScheme.getSelectedItem());
            MSDetailPanel.setPreprocess(checkPreprocess.isSelected());
            MSDetailPanel.setFindFeatures(checkDetectFeatures.isSelected());
            MSDetailPanel.setShowChartRange(checkShowChartRange.isSelected());
        }

    }
}
