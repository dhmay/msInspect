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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.gui.IntensityPlot;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.FloatArray;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;

import javax.swing.*;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;
import java.awt.event.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.Comparator;
import java.util.Arrays;

/**
 * HeatMapPanel - Main JPanel, contains and controls images for
 * each of the individual charges
 */
public class HeatMapPanel extends JPanel implements PropertyChangeListener, ComponentListener
{
    private static Logger _log = Logger.getLogger(HeatMapPanel.class);
    private static final int MAX_CHARGE = 6;
    private static final int RESAMPLE_FREQ = 36;

    private final JFrame frame;

    private FeatureSet featureSet = null;
    private Feature[] featureClone = null;
    private int currentFeatureIndex = -1;

    // image panels register themselves here; add one so we can index directly by charge
    private ChargeImagePanel[] chargeImagePanels = new ChargeImagePanel[MAX_CHARGE + 1];
    private int[] chargeOffset = new int[MAX_CHARGE + 2]; // Add an extra 1 so that we can store charge 6 end offset

    private Comparator sortOrder = featureChargeMassAscComparator;

    private int panelWidth = 600;
    private int panelHeight = 600;

    /**
     *
     */
    public HeatMapPanel(JFrame frame)
    {
        this.frame = frame;

        setPreferredSize(new Dimension(panelWidth, panelHeight));

        // Listen for changes to the run, feature ranges, or selected point
        ApplicationContext.addPropertyChangeListener(SharedProperties.MS_RUN, this);
        ApplicationContext.addPropertyChangeListener(SharedProperties.FEATURE_RANGES, this);
        ApplicationContext.addPropertyChangeListener(SharedProperties.SELECTED_POINT, this);
        ApplicationContext.addPropertyChangeListener(SharedProperties.SELECTED, this);

        setLayout(new GridLayout(0,2,5,5));
        for (int charge = 1; charge <= MAX_CHARGE; charge++)
            add(new ChargePanel(charge));

        getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, 0), "previousFeature");
        getInputMap().put(KeyStroke.getKeyStroke('-'), "previousFeature");
        getActionMap().put("previousFeature", previousFeatureAction);

        getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, 0), "nextFeature");
        getInputMap().put(KeyStroke.getKeyStroke('+'), "nextFeature");
        getInputMap().put(KeyStroke.getKeyStroke('='), "nextFeature"); // no need to hit shift
        getActionMap().put("nextFeature", nextFeatureAction);

        addComponentListener(this);
    }

    /**
     * Check that the currently loaded feature set matches the one currently displayed
     */
    boolean updateFeatureSet(boolean displayWhenDone)
    {
        FeatureSet fs = HeatMapFrame.getDisplayedFeatureSet();

        if ( fs == featureSet && null != featureSet )
        {
            if (displayWhenDone)
                frame.setVisible(true);
            return true;
        }

        IntensityWindowScanner.load(this, fs, displayWhenDone);
        return false;
    }

    /**
     * Update the heat map panel and the individual charge images
     */
    void updatePanel(FeatureSet fs, boolean forceDisplay)
    {
        featureSet = fs;
        if ( null == fs )
            clearChargeImages();
        else
            updateChargeImages();
        // ????? Set highlight

        if ( forceDisplay && null != fs )
            frame.setVisible(true);
    }

    /**
     *
     */
    boolean setHighlight(int index)
    {
        int charge = -1;
        int newFeatureIndex = -1;

        if (index >= 0)
            charge = featureClone[index].charge;

        for (int i = 1; i <= MAX_CHARGE; i++)
            if ( i == charge )
                newFeatureIndex = chargeImagePanels[i].setHighlight(index);
            else
                chargeImagePanels[i].clearHighlight();

        // Highlighted feature did not change
        if (newFeatureIndex == currentFeatureIndex)
            return false;

        currentFeatureIndex = newFeatureIndex;
        return true;
    }

    private boolean isHitCloseEnough(Feature feature, int index)
    {
        if (index < 0 || index >= featureClone.length)
            return false;
        Feature f = featureClone[index];
        return (feature.charge == f.charge &&
                Math.abs(feature.mz - f.mz) <= 0.01  &&
                Math.abs(feature.scan - f.scan) <= 1);
    }

    boolean setHighlight(Feature feature)
    {
        int index = Arrays.binarySearch(featureClone, feature, sortOrder);

        if ( index < 0 )
        {
            // Be a little flexible on mz comparisons
            int i = - index - 1;
            if (isHitCloseEnough(feature, i))
                index = i;
            else if (isHitCloseEnough(feature, i - 1))
                index = i - 1;
            else if (isHitCloseEnough(feature, i + 1))
                index = i + 1;
        }
        return setHighlight(index);
    }

    void broadcastHighlight()
    {
        MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
        Feature f = currentFeatureIndex < 0 ? null : featureClone[currentFeatureIndex];
        ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, f);
        ApplicationContext.setProperty(SharedProperties.SELECTED, f);
        if (null != run)
        {
            int i = getIndexForScanNum(run, f.scan);
            ApplicationContext.setProperty(SharedProperties.MS_SCAN, run.getScan(i));
        }
    }

    void setHighlightAndBroadcast(int index)
    {
        if (setHighlight(index))
            broadcastHighlight();
    }

    /**
     *
     */
    void clearHighlights()
    {
        for (int charge = 1; charge <= MAX_CHARGE; charge++)
            chargeImagePanels[charge].clearHighlight();
        currentFeatureIndex = -1;
    }

    /**
     *
     */
    void clearChargeImages()
    {
        clearHighlights();
        frame.setTitle(getTitle(null));
        for(int charge = 0; charge <= MAX_CHARGE + 1; charge++)
            chargeOffset[charge] = 0;
        for (int charge = 1; charge <= MAX_CHARGE; charge++)
            chargeImagePanels[charge].clearImage();
        repaint();
    }

    /**
     *
     */
    void updateChargeImages()
    {
        clearHighlights();
        frame.setTitle(getTitle(featureSet));

        Feature[] f = featureSet.getFeatures();
        featureClone = new Feature[f.length];
        System.arraycopy(f, 0, featureClone, 0, f.length);
        Arrays.sort(featureClone, sortOrder);

        int charge = 0;
        for (int i = 0; i < f.length; i++)
            while (featureClone[i].charge > charge)
                chargeOffset[charge++] = i;

        for( ; charge <= MAX_CHARGE + 1; charge++)
            chargeOffset[charge] = f.length;

        for (charge = 1; charge <= MAX_CHARGE; charge++)
            chargeImagePanels[charge].update();
        repaint();
    }

    /**
     * Build a frame title from the given feature set
     */
    static String getTitle(FeatureSet fs)
    {
        StringBuffer sb = new StringBuffer(WorkbenchFrame.getAppName());
        sb.append(" - FeatureHeatMap");
        if ( null != fs )
        {
            sb.append(" -- ");
            sb.append(fs.getSourceFile().getName());
        }
        return sb.toString();
    }

    /**
     * Get the index corresponding to the given scan number. Since we may be handed
     * ms2 or ms3 scans, pick the previous scan if not found
     */
    static int getIndexForScanNum(MSRun run, int n)
    {
        int i = run.getIndexForScanNum(n);
        if ( i == -1 )
            return 0;
        if ( i < 0 )
            return -2 - i;
        return i;
    }

    /**
     *
     */
    public void propertyChange(PropertyChangeEvent event)
    {
        _log.debug("Heat map panel saw an event " + event.getPropertyName() + " " + event);

        // Loading a new run or selecting a new set of features invalidates current
        // images even if we are not displayed.
        if ( SharedProperties.MS_RUN.equals(event.getPropertyName()) )
        {
            frame.setVisible(false); // Disable the heat map; no longer valid
            return;
        }

        // Need to know about new feature sets if we are displayed *or*
        // processing another in the background.
        if (!frame.isVisible() && !IntensityWindowScanner.isRunning())
            return;

        // New Feature Set loaded, created, or selected
        if ( SharedProperties.FEATURE_RANGES.equals(event.getPropertyName()) )
        {
            SwingUtilities.invokeLater(new Runnable()
            {
                public void run()
                {
                    updateFeatureSet(false);
                }
            });
        }

        // Don't mark the current feature unless we are visible
        if (!frame.isVisible())
            return;

        // New point or Feature. If Feature, highlight it
        if ( SharedProperties.SELECTED.equals(event.getPropertyName()) ||
             SharedProperties.SELECTED_POINT.equals(event.getPropertyName()) )
        {
            if (null != event.getNewValue() && event.getNewValue() instanceof Feature)
            {
                Feature f = (Feature)event.getNewValue();
                // _log.debug("Heat map saw feature selection (" + f.charge + ", " + f.mz + ", " + f.scan + ")");
                setHighlight(f);
            } else {
                clearHighlights();
            }
        }
    }

    Action previousFeatureAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent e)
        {
            if (currentFeatureIndex <= 0 || currentFeatureIndex >= featureClone.length )
                return;
            int charge = featureClone[currentFeatureIndex].charge;
            int index = currentFeatureIndex - 1;
            if ( charge != featureClone[index].charge )
                return;
            setHighlightAndBroadcast(index);
        }
    };

    Action nextFeatureAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent e)
        {
            if (currentFeatureIndex < 0 || currentFeatureIndex >= featureClone.length - 1 )
                return;
            int charge = featureClone[currentFeatureIndex].charge;
            int index = currentFeatureIndex + 1;
            if ( charge != featureClone[index].charge )
                return;
            setHighlightAndBroadcast(index);
        }
    };

    /**
     *
     */
    public Comparator getSortOrder()
    {
        return sortOrder;
    }

    /**
     *
     */
    public void setSortOrder(Comparator sortOrder)
    {
        this.sortOrder = sortOrder;
        updateChargeImages();
    }

    /**
     *
     */
    public void setSortOrderMass()
    {
        setSortOrder(featureChargeMassAscComparator);
    }

    /**
     *
     */
    public void setSortOrderKl()
    {
        setSortOrder(featureChargeKlDscComparator);
    }

    /**
     *
     */
    public void setSortOrderScan()
    {
        setSortOrder(featureChargeScanAscComparator);
    }

    /**
     *
     */
    public void setSortOrderIntensity()
    {
        setSortOrder(featureChargeIntensityAscComparator);
    }

    /**
     *
     */
    public void setSortOrderTotalIntensity()
    {
        setSortOrder(featureChargeTotalIntensityAscComparator);
    }

    /*----------------------------------------------------------------------*
     * Comparators
     *----------------------------------------------------------------------*/

    /**
     *
     */
    static Comparator featureScanAscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            return (f1.scan < f2.scan ? -1 : (f1.scan == f2.scan ? 0 : 1));
        }
    };

    /**
     *
     */
    static Comparator featureChargeMassAscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.charge != f2.charge )
                return f1.charge < f2.charge ? -1 : 1;
            if ( f1.mz != f2.mz )
                return f1.mz < f2.mz ? -1 : 1;
            return f1.scan < f2.scan ? -1 : (f1.scan > f2.scan ? 1 : 0);
        }
    };

    /**
     *
     */
    static Comparator featureChargeScanAscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.charge != f2.charge )
                return f1.charge < f2.charge ? -1 : 1;
            if ( f1.scan != f2.scan )
                return f1.scan < f2.scan ? -1 : 1;
            return f1.mz < f2.mz ? -1 : (f1.mz > f2.mz ? 1 : 0);
        }
    };

    /**
     *
     */
    static Comparator featureChargeKlDscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.charge != f2.charge )
                return f1.charge < f2.charge ? -1 : 1;
            if ( f1.kl != f2.kl )
                return f1.kl < f2.kl ? 1 : -1;
            if ( f1.scan != f2.scan )
                return f1.scan < f2.scan ? -1 : 1;
            return f1.mz < f2.mz ? -1 : (f1.mz > f2.mz ? 1 : 0);
        }
    };

    /**
     *
     */
    static Comparator featureChargeIntensityAscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.charge != f2.charge )
                return f1.charge < f2.charge ? -1 : 1;
            if ( f1.intensity != f2.intensity )
                return f1.intensity < f2.intensity ? -1 : 1;
            if ( f1.scan != f2.scan )
                return f1.scan < f2.scan ? -1 : 1;
            return f1.mz < f2.mz ? -1 : (f1.mz > f2.mz ? 1 : 0);
        }
    };

    /**
     *
     */
    static Comparator featureChargeTotalIntensityAscComparator = new Comparator()
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.charge != f2.charge )
                return f1.charge < f2.charge ? -1 : 1;
            if ( f1.totalIntensity != f2.totalIntensity )
                return f1.totalIntensity < f2.totalIntensity ? -1 : 1;
            if ( f1.scan != f2.scan )
                return f1.scan < f2.scan ? -1 : 1;
            return f1.mz < f2.mz ? -1 : (f1.mz > f2.mz ? 1 : 0);
        }
    };


    /*----------------------------------------------------------------------*
     * Component actions, currently unused
     *----------------------------------------------------------------------*/
    public void componentResized(ComponentEvent e)
    {
        panelWidth = getSize().width;
        panelHeight = getSize().height;

        for (int i = 1; i <= MAX_CHARGE; i++)
            chargeImagePanels[i].setSize();

        revalidate();
    }

    public void componentHidden(ComponentEvent e)
    {
    }

    public void componentMoved(ComponentEvent e)
    {
    }

    public void componentShown(ComponentEvent e)
    {
    }

    /**
     * Panel displaying the heat map for a single charge
     */
    class ChargePanel extends JPanel
    {
        /**
         *
         */
        public ChargePanel(int charge)
        {
            JLabel chargeLabel = new JLabel("Charge " + charge);
            if (false)
            {
                setLayout(new GridBagLayout());
                GridBagConstraints c = new GridBagConstraints();
                c.anchor = GridBagConstraints.CENTER;
                c.weightx = 0.0;
                c.weighty = 0.0;
                c.gridx = 0;
                c.gridy = 0;
                c.insets = new Insets(4, 4, 4, 4);
                add(new ChargeImagePanel(charge), c);
                c.gridx = 0;
                c.gridy = 1;
                c.anchor = GridBagConstraints.PAGE_END;
                c.insets = new Insets(0, 0, 0, 0);
                add(chargeLabel, c);
            }
            else
            {
                setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
                ChargeImagePanel cip = new ChargeImagePanel(charge);
                cip.setAlignmentX(Component.CENTER_ALIGNMENT);
                add(cip);
                add(chargeLabel);
            }
        }
    }

    /**
     * Image for a single charge
     */
    class ChargeImagePanel extends JPanel implements MouseListener
    {
        private int _charge = 0;
        BufferedImage img;
        int imgWidth, imgHeight;
        int highlightX = -1;
        float xscale;

        /**
         *
         */
        public ChargeImagePanel(int charge)
        {
            _charge = charge;
            chargeImagePanels[_charge] = this;
            setBorder(BorderFactory.createLineBorder(Color.BLACK));
            setSize();
            addMouseListener(this);
        }

        public void setSize()
        {
            int w = panelWidth/2 - 10;
            int h = panelHeight/3 - 20;
            setPreferredSize(new Dimension(imgWidth, imgHeight));
//            setMinimumSize(new Dimension(imgWidth, imgHeight));
            revalidate();
            update();
            repaint();
        }

        public void update()
        {
            imgWidth = getWidth();
            imgHeight = getHeight();

            // Can't do anything yet...
            if (imgWidth == 0 || imgHeight == 0)
                return;

            img = new BufferedImage(imgWidth, imgHeight, BufferedImage.TYPE_3BYTE_BGR);
            Graphics2D g = (Graphics2D)img.getGraphics();
            g.setBackground(Color.GRAY);
            g.clearRect(0, 0, imgWidth, imgHeight);

            float maxColor = (1 << 8) - 1;

            int start = chargeOffset[_charge - 1];
            int end = chargeOffset[_charge];

            if ( end - start <= 0 )
                return;

            xscale = (float)imgWidth / (end - start);
            if (xscale > 20.f)
                xscale = 20.f;

            FloatArray x = new FloatArray();
            FloatArray y = new FloatArray();
            FloatArray z = new FloatArray();

            // TODO: Minor speedup: if i is the same on consecutive iterations, just copy image column
            for (int xi = 0; xi < imgWidth; xi++)
            {
                int i = (int)(xi / xscale) + start;

                if ( i >= end )
                    break;

                // Can happen if feature set is out of range for this run (user error)...
                if ( null == featureClone[i].intensityWindow )
                    continue;

                float yscale = imgHeight / (featureClone[i].intensityWindow.length - 1.f);

                // Scale each window so that monoisotopic peak is maximum intensity,
                // then clamp so IntensityPlot does not rescale
                float intensityScale = 0.f;
                int center = featureClone[i].intensityLeadingPeaks * RESAMPLE_FREQ;
                for (int j = center - 10; j < center + 10; j++)
                    if (featureClone[i].intensityWindow[j] > intensityScale)
                        intensityScale = featureClone[i].intensityWindow[j];

                if (intensityScale == 0.f)
                {
                    // no signal near center peak; consider entire window
                    for (int j = 0; j < featureClone[i].intensityWindow.length; j++)
                        if (featureClone[i].intensityWindow[j] > intensityScale)
                            intensityScale = featureClone[i].intensityWindow[j];
                }

                if ( intensityScale != 0.f )
                    intensityScale = maxColor / intensityScale;
                else
                    intensityScale = 1.f; // punt

                int maxy = featureClone[i].intensityWindow.length;
                for (int yi = 0; yi < imgHeight; yi++)
                {
                    x.add((float)xi);
                    y.add((float)yi);
                    int j = (int) (yi / yscale);
                    float c = featureClone[i].intensityWindow[j] * intensityScale;
                    z.add( c > maxColor ? maxColor : (c < 0 ? 0 : c));
                }
            }

            IntensityPlot plot = new IntensityPlot();
            String scheme = "Heat";
            plot.setData(x, y, z);
            plot.plotSqrt(img, -1.f, 1, scheme);
        }

        public int setHighlight(int featureIndex)
        {
            // _log.debug("  Highlight charge " +_charge + ": " + chargeOffset[_charge - 1] + " <= " + featureIndex + " < " + chargeOffset[_charge]);

            if ( featureIndex < chargeOffset[_charge - 1] ||
                 featureIndex >= chargeOffset[_charge] )
            {
                clearHighlight();
                return -1;
            }
            highlightX = (int)((featureIndex - chargeOffset[_charge - 1] + .5) * xscale);
            repaint();
            return featureIndex;
        }

        public void clearImage()
        {
            int width = img.getWidth();
            int height = img.getHeight();

            Graphics2D g = (Graphics2D)img.getGraphics();
            g.setBackground(Color.GRAY);
            g.clearRect(0, 0, width, height);
        }


        public void clearHighlight()
        {
            if ( highlightX >= 0 )
            {
                highlightX = -1;
                repaint();
            }
        }

        /**
         *
         */
        public int getCharge()
        {
            return _charge;
        }

        /**
         *
         */
        public void paintComponent(Graphics g)
        {
            if ( img != null )
            {
                g.drawImage(img, 0, 0, this);
                if ( highlightX >= 0 && highlightX < imgWidth)
                {
                    g.setColor(Color.BLUE);
                    // line
                    g.drawLine(highlightX, 0, highlightX, imgHeight - 1);
                    // top mark
                    g.fillPolygon(new int[]{highlightX-3, highlightX+3, highlightX},
                                  new int[]{0, 0, 7},
                                  3);
                    // bottom mark
                    g.fillPolygon(new int[]{highlightX-3, highlightX+3, highlightX},
                                  new int[]{imgHeight, imgHeight, imgHeight-7},
                                  3);
                }
            }
        }

        /*----------------------------------------------------------------------*
         * Mouse events
         *----------------------------------------------------------------------*/

        /**
         *
         */
        public void mouseClicked(MouseEvent e)
        {
            ChargeImagePanel chargeImagePanel = (ChargeImagePanel)e.getComponent();

            clearHighlights();

            int featureIndex = (int)(e.getX() / xscale) + chargeOffset[_charge - 1];

            if ( featureIndex < chargeOffset[_charge] &&
                 featureIndex < featureClone.length )
                HeatMapPanel.this.setHighlightAndBroadcast(featureIndex);
        }

        public void mousePressed(MouseEvent e)
        {
        }

        public void mouseReleased(MouseEvent e)
        {
        }

        public void mouseEntered(MouseEvent e)
        {
        }

        public void mouseExited(MouseEvent e)
        {
        }

    }

    /**
     * Runnable to gather intensity windows for features that do not yet
     * have them
     */
    static class IntensityWindowScanner implements Runnable
    {
        static IntensityWindowScanner instance = null;
        static Object lock = new Object();
        Thread worker = null;

        HeatMapPanel panel = null;
        FeatureSet featureSet = null;
        boolean displayWhenDone;

        private IntensityWindowScanner(HeatMapPanel panel, FeatureSet featureSet, boolean displayWhenDone)
        {
            this.panel = panel;
            this.featureSet = featureSet;
            this.displayWhenDone = displayWhenDone;
        }

        static void load(HeatMapPanel panel, FeatureSet featureSet, boolean displayWhenDone)
        {
            synchronized(lock)
            {
                // Abort current load, if any
                if ( null != instance && null != instance.worker )
                {
                    instance.worker.interrupt();
                    instance.worker = null;
                    instance = null;
                }
                instance = new IntensityWindowScanner(panel, featureSet, displayWhenDone);

                Thread w = new Thread(instance);
                w.setPriority(Thread.MIN_PRIORITY);
                w.start();

                instance.worker = w;
            }
        }

        public void run()
        {
            Thread currentThread = Thread.currentThread();

            try
            {
                ApplicationContext.setMessage("Started extracting intensity windows...");

                if (null == featureSet)
                {
                    ApplicationContext.setMessage("HeatMap: Could not find a displayed feature set.");
                    update();
                    return;
                }

                MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
                Feature[] features = featureSet.getFeatures();

                if (null == run || null == features || 0 >= features.length)
                {
                    ApplicationContext.setMessage("HeatMap: Feature set contained nothing to display.");
                    featureSet = null;
                    update();
                    return;
                }

                Feature[] f = new Feature[features.length];
                System.arraycopy(features, 0, f, 0, f.length);
                Arrays.sort(f, featureScanAscComparator);

                MSRun.MSScan scan = null;
                float[][] spectrum = null;
                for (int i = 0; i < f.length; i++)
                {

                    if ( 0 == (i % 100) )
                    {
                        float p = Math.round(1000.f * i / f.length) / 10.f;
                        ApplicationContext.setMessage("Extracting intensity windows, " + String.valueOf(p) + "% complete...");
                    }

                    // If feature intensity window was already extracted, just skip.
                    if ( f[i].intensityLeadingPeaks == 3 && f[i].intensityTrailingPeaks == 3 )
                        continue;

                    int n = getIndexForScanNum(run, f[i].scan);
                    // This can happen if user, in error, applies feature set from one run to another.
                    if ( n >= run.getScanCount() )
                        continue;

                    scan = run.getScan(n);
                    spectrum = scan.getSpectrum();
                    if (null == spectrum)
                    {
                        _log.error("Failed to get spectrum for scan " + f[i].scan);
                        ApplicationContext.setMessage("Failed to get spectrum for scan " + f[i].scan);
                        featureSet = null;
                        update();
                        return;
                    }

                    if (currentThread.isInterrupted())
                        throw new InterruptedException();

                    f[i].intensityWindow = Spectrum.Resample(spectrum, new FloatRange(f[i].mz - 3, f[i].mz + 3), RESAMPLE_FREQ);
                    f[i].intensityLeadingPeaks = 3;
                    f[i].intensityTrailingPeaks = 3;
                }
                ApplicationContext.setMessage("");
                update();
            }
            catch (InterruptedException e)
            {
                _log.debug("loader interrupted");
                // just ignore; we're aborting the load if we get here.
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Error extracting intensity windows", e);
                ApplicationContext.setMessage("Error extracting intensity windows");
            }
            finally
            {
                synchronized(lock)
                {
                    worker = null;
                }
            }
        }

        static synchronized boolean isRunning()
        {
            synchronized(lock)
            {
                return null != instance && null != instance.worker;
            }
        }

        /**
         *
         */
        private void update()
        {
            final HeatMapPanel hmp = panel;
            final FeatureSet fs = featureSet;
            final boolean dwd = displayWhenDone;
            SwingUtilities.invokeLater(new Runnable()
            {
                public void run()
                {
                    hmp.updatePanel(fs, dwd);
                }
            });
        }

    }

}
