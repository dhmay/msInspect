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

import org.fhcrc.cpl.toolbox.proteomics.gui.IntensityPlot;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.gui.SavePartialMzxmlDialog;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.fhcrc.cpl.toolbox.gui.ScrollableImage;
import org.fhcrc.cpl.toolbox.gui.ImagePanelLayer;
import org.fhcrc.cpl.toolbox.gui.BaseImagePanelLayer;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.apache.log4j.Logger;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.PopupMenuListener;
import javax.swing.event.PopupMenuEvent;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.beans.PropertyChangeEvent;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.prefs.Preferences;


/**
 * User: mbellew
 * Date: May 26, 2004
 * Time: 10:06:01 AM
 */
public class MSImageComponent
{
    protected static Logger _log = Logger.getLogger(MSImageComponent.class);

    public static String getPrefColorScheme()
    {
        String colorScheme =
                (String) ApplicationContext.getProperty(MSImageComponent.COLORSCHEME_PROPNAME);
        if (colorScheme == null)
        {
            Preferences prefs = Preferences.userNodeForPackage(Application.class);
            colorScheme = prefs.get(MSImageComponent.COLORSCHEME_PROPNAME,
                    MSImageComponent.DEFAULT_COLORSCHEME);
        }
        return colorScheme;
    }

    public static void initializeColorMenu(JMenu colorMenu)
    {
        //if only one language supported, disable the menu and hide it
        for (String colorScheme : IntensityPlot.COLOR_SCHEMES)
        {
            JMenuItem menuItem = new JMenuItem(colorScheme);
            ColorActionListener listener = new ColorActionListener(colorScheme);
            menuItem.addActionListener(listener);
            colorMenu.add(menuItem);
        }
    }

    /**
     * actionlistener class for color menu items
     */
    protected static class ColorActionListener implements ActionListener {
        String _colorScheme;

        public ColorActionListener(String colorScheme)
        {
            _colorScheme = colorScheme;
        }

        public void actionPerformed(ActionEvent event)
        {
            ApplicationContext.setProperty(MSImageComponent.COLORSCHEME_PROPNAME, _colorScheme);
        }
    }

    public static final int DEFAULT_FEATURE_SCAN_TOLERANCE = 3;
    public static final float DEFAULT_FEATURE_MZ_TOLERANCE = 3;

    public static final String COLORSCHEME_PROPNAME = "MSImageComponent.colorScheme";
    public static final String DEFAULT_COLORSCHEME = "Fancy";

    static Font _smallFont = Font.decode("Verdana PLAIN 9");
    static Color _ruleColor = new Color(255, 255, 192); // plaf??

    private JScrollPane scrollPane;
    private MSImagePanel imagePanel;
    private ColumnPanel columnPanel;
    private RowPanel rowPanel;

    private Rectangle highlightRect; // in mz/scan coordinates
    private int highlightScan = -1;
    private double[] highlightMasses = null;
    private java.util.List<FeatureSet> featureSets;

    private int ROWHEADER_WIDTH = 30;
    private int COLHEADER_HEIGHT = 25;

    private MSRun _run = null;

    private ListenerHelper helper = new ListenerHelper(this);

    //the current zoom level
    protected static double _zoomFactor = 1;
    //right-click context menu
    protected RightMousePopupMenu _rightMousePopupMenu = null;

    protected boolean _inZoomProcess = false;

    protected Point _currentMousePosition;

    //the effective size of the image panel, taking zoom into account
    Dimension _imagePanelSize;
    Point _lastMouseRightClick;

    java.util.List<ImagePanelLayer> imagePanelLayers = new ArrayList<ImagePanelLayer>();

    //This single tree contains all features, from all _displayed_ featuresets.  The two dimensions
    //are scan _index_ and m/z.  This structure is ideal for quickly finding the closest
    //single feature to a given point (e.g., the mouse position).
    //For MS/MS features, the actual index is not stored here.  Instead, the closest
    //lower single MS scan index is stored.
    protected Tree2D _allDisplayedFeaturesScanIndexMzTree;

    public MSImageComponent()
    {
        imagePanel = new MSImagePanel();
        columnPanel = new ColumnPanel();
        rowPanel = new RowPanel();
//		columnPanel.setPreferredSize(new Dimension(640, COLHEADER_HEIGHT));
//		rowPanel.setPreferredSize(new Dimension(ROWHEADER_WIDTH, 400));

        helper.addListener(Application.getInstance(), "app_propertyChange");

        EventQueue.invokeLater(new Runnable()
        {
            public void run()
            {
                setRun((MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN));
            }
        });

        _allDisplayedFeaturesScanIndexMzTree = new Tree2D();
    }

    public MSImageComponent(Image image)
    {
        this();
        setImage(image);
    }


    public void setRun(MSRun run)
    {
        if (_run == run)
            return;
        _run = run;
        if (null == run)
            setImage(null);
        else
            setImage(run.getImage(getPrefColorScheme()));
    }

    /**
     * Flashy transition to new zoom level.  Swoopy
     * @param newZoomFactor
     */
    public void transitionToZoomFactor(double newZoomFactor)
    {
        if (newZoomFactor == _zoomFactor)
            return;

        //more steps for more change in zoom level
        int numSteps=(int)(Math.abs(newZoomFactor - _zoomFactor) * 15);

        double startingZoomFactor = _zoomFactor;

        //indicate that this zoom step is an intermediate step
        _inZoomProcess=true;

        for (double i=1; i<= numSteps-1; i++)
        {
            double currentStep = startingZoomFactor +
                    (i * ((newZoomFactor - startingZoomFactor) / numSteps));
            setZoomFactor(currentStep);
            try
            {
                Thread.sleep(5);
            }
            catch (Exception e) {}
        }
        _inZoomProcess=false;
        setZoomFactor(newZoomFactor);

        _rightMousePopupMenu.updateZoomList();
    }



    /**
     * Recenter on a feature, and populate the detail panes appropriately
     * @param feature
     */
    public void selectFeature(Feature feature)
    {
        int scanNum = feature.getScan();
        int scanIndex = scanNum;
        if (null != _run && scanNum <= _run.getScanCount())
        {

            scanIndex = _run.getIndexForScanNum(scanNum);
        }

        float mz = feature.getMz();

        ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, feature);
        ApplicationContext.setProperty(SharedProperties.SELECTED, feature);
        if (null != _run)
        {
            int n = _run.getIndexForScanNum(feature.scan, true);
            ApplicationContext.setProperty(SharedProperties.MS_SCAN, _run.getScan(n));
        }

        Point mousePoint = convertScanMzToMousePoint(scanIndex, mz);
        recenter((int) mousePoint.getX(), (int) mousePoint.getY());
        imagePanel.selectPoint(scanIndex, mz);
    }

    /**
     * Tree2D doesn't allow node removal, so have to rebuild from the ground up.
     * MS2 scans are assigned the scan index of their _precursor scan_
     */
    protected void rebuildScanIndexMzTree()
    {
        MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);

        //if no run is loaded, nothing to do here
        if (run == null)
            return;

        _allDisplayedFeaturesScanIndexMzTree = new Tree2D();


        for (FeatureSet featureSet : featureSets)
        {
            if (featureSet.isDisplayed())
            {

                for (Feature feature : featureSet.getFeatures())
                {
                    int scanIndex = run.getIndexForScanNum(feature.getScan());

                    if (scanIndex < 0)
                    {
                        int ms2ScanIndex = run.getIndexForMS2ScanNum(feature.getScan());

                        if (ms2ScanIndex > 0)
                        {
                            MSRun.MSScan ms2Scan = run.getMS2Scan(ms2ScanIndex);
                            int precursorScanNum = ms2Scan.getPrecursorScanNum();
                            if (precursorScanNum < 0)
                                precursorScanNum = -precursorScanNum;
                            int precursorScanIndex = run.getIndexForScanNum(precursorScanNum);
                            if (precursorScanIndex < 0)
                                precursorScanIndex = -precursorScanIndex;

                            if (precursorScanIndex > run.getScanCount() - 1)
                                continue;

                            MSRun.MSScan precursorScan =
                                    run.getScan(precursorScanIndex);
                            scanIndex = precursorScan.getIndex();
                        }
                        else
                        {
                            scanIndex = -scanIndex;
                        }
                    }
                    _allDisplayedFeaturesScanIndexMzTree.add(scanIndex, feature.getMz(),
                                                             feature);
                }
            }
        }
    }


    /**
     * convert a point given by a MouseEvent into a scan and mz value,
     * returned as x and y of a new Point, respectively
     * @param mousePoint
     * @return
     */
    public Point convertMousePointToScanIndexMz(Point mousePoint)
    {
        int x = (int) mousePoint.getX();
        int y = (int) mousePoint.getY();
        //scale to account for zoom factor
        x = (int) ( x / _zoomFactor);
        y = (int) ( y / _zoomFactor);

        y = imagePanel.getImageHeight() - y - 1;

        return new Point(x,y);
    }

    public Point convertMousePointToScanMz(Point mousePoint)
    {
        Point scanIndexMzPoint = convertMousePointToScanIndexMz(mousePoint);
        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        MSRun.MSScan scan = run.getScan((int) scanIndexMzPoint.getX());
        int scanNum = (int) scanIndexMzPoint.getX();
        if (null != scan)
        {
            scanNum = scan.getNum();
        }
        return new Point(scanNum, (int) scanIndexMzPoint.getY());
    }

    /**
     * Convert a scan(index), mz pair to a mouse point
     * @param scanIndex
     * @param mz
     * @return
     */
    public Point convertScanMzToMousePoint(int scanIndex, float mz)
    {
        int x = (int) (_zoomFactor * scanIndex);
        int y = (int) ((_zoomFactor * imagePanel.getImageHeight()) -
                (_zoomFactor * mz));
        return new Point(x,y);
    }

    /**
     * Change the zoom level to a new value and redraw everything immediately.
     * Set the scrollbars appropriately
     * @param zoomFactor
     */
    public void setZoomFactor(double zoomFactor)
    {
        if (zoomFactor == _zoomFactor)
            return;
        _zoomFactor = zoomFactor;

        _imagePanelSize = new Dimension((int)(_zoomFactor * imagePanel.getImageWidth()),
                (int)(_zoomFactor * imagePanel.getImageHeight()));
        imagePanel.invalidate();
//        imagePanel.revalidate();
        rowPanel.invalidate();
//        rowPanel.revalidate();
        columnPanel.invalidate();
//        columnPanel.revalidate();
        getScrollPane().invalidate();
        getScrollPane().revalidate();

        //calling paintImmediately() to facilitate swoopiness
        getScrollPane().paintImmediately(0,0,getScrollPane().getWidth(),getScrollPane().getHeight());

        int horizScrollCenter = (int) (_zoomFactor * _lastMouseRightClick.getX());
        int vertScrollCenter = (int) (_zoomFactor * _lastMouseRightClick.getY());

        recenter(horizScrollCenter, vertScrollCenter);
    }

    /**
     * Recenter on the specified point, with x and y in mouse coordinates
     * @param x
     * @param y
     */
    protected void recenter(int x, int y)
    {
        Rectangle viewportRect = getScrollPane().getViewportBorderBounds();

        int horizScrollBarValue = x - (int) ((viewportRect.getWidth()) / 2);
        int vertScrollBarValue = y - (int) ((viewportRect.getHeight()) / 2);

        //for some reason it's necessary to call scrollRectToVisible() before setting
        //the scrollbar values.  Otherwise, the scrollbar values get capped
        getScrollPane().getViewport().scrollRectToVisible(new Rectangle(horizScrollBarValue, vertScrollBarValue, 1,1));
        getScrollPane().getHorizontalScrollBar().setValue(horizScrollBarValue);
        getScrollPane().getVerticalScrollBar().setValue(vertScrollBarValue);
    }

    //return the current zoom factor
    public static double getZoomFactor()
    {
        return _zoomFactor;
    }

    /**
     * Handling for property changes we care about: run, scan, zoom region, feature
     * ranges, highlighted features
     * @param event
     */
    public void app_propertyChange(PropertyChangeEvent event)
    {
        String propName = event.getPropertyName();
        if (SharedProperties.MS_RUN.equals(propName))
        {
            MSRun run = (MSRun)event.getNewValue();
            setRun(run);
        }
        else if (SharedProperties.MS_SCAN.equals(propName))
        {
            if (null == event.getNewValue())
                highlightScan = -1;
            else
            {
                MSRun.MSScan scan = (MSRun.MSScan)event.getNewValue();
                highlightScan = scan.getIndex();
            }
            imagePanel.repaint();
        }
        else if (SharedProperties.ZOOM_REGION.equals(propName))
        {
            if (null == event.getNewValue() || event.getNewValue() instanceof Rectangle)
                highlightRect = (Rectangle)event.getNewValue();
            imagePanel.repaint();
        }
        else if (SharedProperties.FEATURE_RANGES.equals(propName))
        {
            featureSets = (java.util.List<FeatureSet>)event.getNewValue();
            _log.debug("Got a FEATURE_RANGES app property event.  New number of sets: " + featureSets.size());
            int newSetsIndex = 0;
            for (int i=0; i<imagePanelLayers.size(); i++)
            {
                ImagePanelLayer oldLayer = imagePanelLayers.get(i);

                if (oldLayer instanceof MSImagePanel.FeatureSetImagePanelLayer)
                {
                    if (i < featureSets.size() - 1)
                    {
                        int displayStyle =
                                MSImagePanel.FeatureSetImagePanelLayer.DISPLAY_STYLE_X;
                        if (0 == (i & 1))
                            displayStyle =
                                    MSImagePanel.FeatureSetImagePanelLayer.DISPLAY_STYLE_PLUS;
                        MSImagePanel.FeatureSetImagePanelLayer newSetLayer =
                                imagePanel.createFeatureSetLayer(featureSets.get(newSetsIndex++),
                                        "Feature Set " + newSetsIndex, displayStyle);

                        imagePanelLayers.set(i, newSetLayer);
                    }
                    else
                    {
                        removeImagePanelLayer(i--);
                        //if we just removed the last layer, we're done
                        if (i >= imagePanelLayers.size() - 1)
                            break;
                    }
                }
            }
            for (int i=newSetsIndex; i<featureSets.size(); i++)
            {
                int displayStyle = MSImagePanel.FeatureSetImagePanelLayer.DISPLAY_STYLE_X;
                if (0 == (i & 1))
                    displayStyle = MSImagePanel.FeatureSetImagePanelLayer.DISPLAY_STYLE_PLUS;
                FeatureSet newFeatureSet = (FeatureSet) featureSets.get(newSetsIndex++);
                _log.debug("Displaying new feature set " + (i+1) + ", source file: " +
                           newFeatureSet.getSourceFile().getAbsolutePath() +
                           ", features: " + newFeatureSet.getFeatures().length);
                MSImagePanel.FeatureSetImagePanelLayer newSetLayer =
                        imagePanel.createFeatureSetLayer(newFeatureSet,
                                "Feature Set " + newSetsIndex,
                                displayStyle);
                addImagePanelLayer(newSetLayer);
            }
            //rebuild the tree that holds all features' scan indices and m/z info
            rebuildScanIndexMzTree();
            imagePanel.repaint();
        }
        else if (SharedProperties.HIGHLIGHT_FEATURES.equals(propName) ||
                SharedProperties.HIGHLIGHT_FEATURES2.equals(propName) ||
                SharedProperties.HIGHLIGHT_FEATURES3.equals(propName))
        {
            imagePanel.repaint();
        }
        else if (MSImageComponent.COLORSCHEME_PROPNAME.equals(propName))
        {
            Preferences prefs = Preferences.userNodeForPackage(Application.class);
            prefs.put(MSImageComponent.COLORSCHEME_PROPNAME,
                    (String) ApplicationContext.getProperty(MSImageComponent.COLORSCHEME_PROPNAME));

            if (_run != null)
            {
                _run.invalidateImage();
                setImage(_run.getImage(getPrefColorScheme()));
                imagePanel.repaint();              
            }
        }
    }

    /**
     * Listens for events from a transparency slider and adjusts the transparency
     * of a layer accordingly
     */
    protected class LayerAdjustmentListener implements AdjustmentListener
    {
        protected ImagePanelLayer layer = null;

        public LayerAdjustmentListener(ImagePanelLayer layer)
        {
            this.layer = layer;
        }
        public void adjustmentValueChanged(AdjustmentEvent e)
        {
            layer.setTranslucence(e.getValue());
            imagePanel.repaint();
        }
    }

    /**
     * Menu action.  Shows sliders for adjusting the transparency of all
     * adjustable layers
     *
     *
     */
    public void showLayerTransparencyDialog()
    {
        JDialog dialog = new JDialog(ApplicationContext.getFrame(),
                TextProvider.getText("LAYER_TRANSPARENCY"));
        dialog.setLocation(ApplicationContext.getFrame().getX() + 50,
                ApplicationContext.getFrame().getY() + 50);
        dialog.setLayout(new FlowLayout());

        dialog.add(new JPanel());
        for (ImagePanelLayer layer : imagePanelLayers)
        {
            JPanel scrollBarPanel = new JPanel();
            scrollBarPanel.setLayout(new GridBagLayout());

            scrollBarPanel.setSize(new Dimension(80,200));

            JScrollBar layerScrollBar = new JScrollBar();
            layerScrollBar.setMinimum(0);
            layerScrollBar.setMaximum(255);
            layerScrollBar.setPreferredSize(new Dimension(20,150));
            layerScrollBar.setSize(new Dimension(20,150));
            layerScrollBar.setMaximumSize(new Dimension(20,150));            
            layerScrollBar.setValue(layer.getTranslucence());

            GridBagConstraints scrollBarGridBagConstraints =
                    new GridBagConstraints();
            scrollBarGridBagConstraints.gridwidth =
                    GridBagConstraints.REMAINDER;
            scrollBarPanel.add(layerScrollBar, scrollBarGridBagConstraints);
            scrollBarPanel.add(new JLabel(layer.getName()));

            layerScrollBar.addAdjustmentListener(new LayerAdjustmentListener(layer));
            dialog.add(scrollBarPanel);
        }
        dialog.add(new JPanel());


        dialog.setSize(90 * imagePanelLayers.size() + 20,200);
        dialog.setVisible(true);
    }

    public void addImagePanelLayer(ImagePanelLayer newLayer)
    {
        imagePanelLayers.add(newLayer);
        ((WorkbenchFrame) ApplicationContext.getFrame()).showLayerTransparencyAction.setEnabled(true);
    }

    public void removeImagePanelLayer(int i)
    {
        imagePanelLayers.remove(i);
        if (imagePanelLayers.size() == 0)
            ((WorkbenchFrame) ApplicationContext.getFrame()).showLayerTransparencyAction.setEnabled(false);                
    }

    public void setImage(Image image)
    {
        imagePanel.setImage(image);

        if (null != image)
        {
            _imagePanelSize = new Dimension((int)(_zoomFactor * image.getWidth(imagePanel)),
                    (int)(_zoomFactor * image.getHeight(imagePanel)));
        }
        if (null != scrollPane)
        {
            rowPanel.revalidate();
            columnPanel.revalidate();
            imagePanel.revalidate(); //Make sure the scroll pane knows our size has changed
            scrollPane.revalidate();
            scrollPane.paintImmediately(0,0,getScrollPane().getWidth(),getScrollPane().getHeight());
        }
    }


    public MSImagePanel getImagePanel()
    {
        return imagePanel;
    }


    public JComponent getColumnPanel()
    {
        return columnPanel;
    }


    public JComponent getRowPanel()
    {
        return rowPanel;
    }


    public JScrollPane getScrollPane()
    {
        if (null == scrollPane)
        {
            setScrollPane(new JScrollPane());


        }
        return scrollPane;
    }


    public void setScrollPane(JScrollPane scrollPane)
    {
        scrollPane.setViewportView(imagePanel);
        scrollPane.setColumnHeaderView(columnPanel);
        scrollPane.setRowHeaderView(rowPanel);

        // multiple monitors don't work with 'smart' scrolling
        scrollPane.getRowHeader().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
        scrollPane.getColumnHeader().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
        scrollPane.getViewport().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);

        this.scrollPane = scrollPane;
    }

    /**
     * horizontal scale bar
     */
    class ColumnPanel extends JComponent
    {
        ColumnPanel()
        {
            setBackground(new Color(0x00d0ffef));
        }

        public Dimension getPreferredSize() {
            if (_imagePanelSize == null)
                return super.getPreferredSize();
            return new Dimension(imagePanel.getWidth(),COLHEADER_HEIGHT);
        }

        public void setWidth(int width)
        {
            setPreferredSize(new Dimension(width, 20));
            setSize(width, 20);
            invalidate();
        }

        NumberFormat onedecimal = NumberFormat.getInstance();
        {
            onedecimal.setMinimumFractionDigits(1);
            onedecimal.setMaximumFractionDigits(1);
            onedecimal.setGroupingUsed(false);
        }

        public void paint(Graphics graphics)
        {
            int height = COLHEADER_HEIGHT;
            int width = imagePanel.getImageWidth();
            int viewWidth = (int) (_zoomFactor * width);
            graphics.setColor(_ruleColor); // taupe
            graphics.fillRect(0, 0, viewWidth, height);
            graphics.setColor(Color.BLACK);
            graphics.setFont(_smallFont);
            int distanceBetweenTicks = Math.max(1,(int) (10 / _zoomFactor));
            for (int x = distanceBetweenTicks; x < width; x += distanceBetweenTicks)
            {
                int viewX = (int)(_zoomFactor * x);
                int index = x -1 ;

                if (0 == (x % (10 * distanceBetweenTicks)))
                {
                    if (null != _run && x <= _run.getScanCount())
                    {
                        MSRun.MSScan scan = _run.getScan(index);
                        double t = scan.getDoubleRetentionTime();
                        String time = null;
                        if (t <= 0.0)
                            time = "" + scan.getNum();
                        else
                            time = "" + onedecimal.format(new Double(t));
                        graphics.drawString(time, viewX - 7, height - 10);
                    }
                }
                int viewIndex = (int) (_zoomFactor * index);
                if (0 == (x % (5 * distanceBetweenTicks)))
                    graphics.drawLine(viewIndex, height - 7, viewIndex, height - 1);
                else
                    graphics.drawLine(viewIndex, height - 4, viewIndex, height - 1);
            }
        }
    }

    /**
     * vertical scale bar
     */
    class RowPanel extends JComponent
    {
        RowPanel()
        {
            // taupe more or less
            setBackground(new Color(0x00d0ffef));
        }


        public void setHeight(int height)
        {
            setPreferredSize(new Dimension(ROWHEADER_WIDTH, height));
            setSize(ROWHEADER_WIDTH, height);
            invalidate();
        }

        public Dimension getPreferredSize() {
            if (_imagePanelSize == null)
                return super.getPreferredSize();
            return new Dimension(ROWHEADER_WIDTH, imagePanel.getHeight());
        }

        public void paint(Graphics graphics)
        {
            int height = imagePanel.getImageHeight();
            int viewHeight = (int) (_zoomFactor * imagePanel.getImageHeight());
            int width = ROWHEADER_WIDTH;
            graphics.setColor(_ruleColor); // taupe
            graphics.fillRect(0, 0, width, viewHeight);
            graphics.setColor(Color.BLACK);
            graphics.setFont(_smallFont);
            int distanceBetweenTicks = Math.max(1,(int) (50 / _zoomFactor));
            for (int y = 0; y < height; y += distanceBetweenTicks)
            {
                int viewY = (int) (y * _zoomFactor);

                graphics.drawString("" + y, 0, viewHeight - viewY + 3);
                graphics.drawLine(width - 4, viewHeight - viewY - 1, width - 1, viewHeight - viewY - 1);

            }
        }
    }




    public class MSImagePanel extends ScrollableImage
    {
        public FeatureSetImagePanelLayer createFeatureSetLayer(FeatureSet featureSet,
                                                               String name,
                                                               int displayStyle)
        {
            return new MSImagePanel.FeatureSetImagePanelLayer(featureSet, name,
                                                              displayStyle);
        }

        public class FeatureSetImagePanelLayer extends BaseImagePanelLayer
        {
            protected FeatureSet featureSet;
            protected int displayStyle = DISPLAY_STYLE_X;

            public static final int DISPLAY_STYLE_X = 0;
            public static final int DISPLAY_STYLE_PLUS = 1;

            public FeatureSetImagePanelLayer(FeatureSet newFeatureSet,
                                             String name,
                                             int newDisplayStyle)
            {
                _log.debug("FeatureSetImagePanelLayer: constructor, featureset name = " + name);
                featureSet = newFeatureSet;
                displayStyle = newDisplayStyle;
                setName(name);             
            }


            public void draw(Graphics graphics, int imageWidth, int imageHeight,
                             int minVisibleScan, int maxVisibleScan,
                             double minVisibleMz, double maxVisibleMz)
            {
                //don't display feature sets during transition -- too slow
                _log.debug("FeatureSetImagePanelLayer.draw 1, name=" + getName());

                if (_inZoomProcess || !featureSet.isDisplayed())
                    return;
                Feature[] featureRanges = featureSet.getFeatures();
                if (featureRanges == null)
                    return;

                graphics.setColor(adjustColor(featureSet.getColor()));

                _log.debug("FeatureSetImagePanelLayer.draw, beginning to draw " + featureRanges.length +
                           " features.");
                for (Feature f : featureRanges)
                {

                    if (f.mz < minVisibleMz || f.mz > maxVisibleMz)
                        continue;
                    // UNDONE: use retention time!!!!
                    int scanIndex = _run.getIndexForScanNum(f.scan);
                    scanIndex = scanIndex < 0 ? -(scanIndex+1) : scanIndex;

                    if (scanIndex < minVisibleScan || scanIndex > maxVisibleScan ||
                            f.mz < minVisibleMz || f.mz > maxVisibleMz)
                        continue;
                    int y = imageHeight - (int)Math.floor(f.mz) - 1;

                    if (IsotopicLabelExtraInfoDef.hasLabel(f) &&
                            IsotopicLabelExtraInfoDef.getLabelCount(f) > 0)
                    {
                        float mzHeavy = f.mz +
                                (IsotopicLabelExtraInfoDef.getLabelCount(f) *
                                        IsotopicLabelExtraInfoDef.getLabel(f).getHeavy() / f.charge);
                        int y2 = imageHeight - (int)Math.floor(mzHeavy) - 1;
                        graphics.drawLine(scanIndex - 2, y, scanIndex + 2, y);
                        graphics.drawLine(scanIndex - 1, y-1, scanIndex + 1, y-1);

                        graphics.drawLine(scanIndex-2, y2, scanIndex+2, y2);
                        graphics.drawLine(scanIndex-1, y2+1, scanIndex+1, y2+1);

                        graphics.drawLine(scanIndex, y, scanIndex, y2);
                    }
                    else
                    {
                        switch (displayStyle)
                        {
                            case DISPLAY_STYLE_PLUS:
                                // Draw "odd" sets as a '+'
                                graphics.drawLine(scanIndex - 2, y, scanIndex + 2, y);
                                graphics.drawLine(scanIndex, y+2, scanIndex, y-2);
                                break;
                            case DISPLAY_STYLE_X:
                                // Draw "even" sets as a 'X'
                                graphics.drawLine(scanIndex - 2, y-2, scanIndex + 2, y+2);
                                graphics.drawLine(scanIndex - 2, y+2, scanIndex + 2, y-2);
                                break;
                        }
                    }
                }
            }
        }



        /**
         * The topmost layer:  highlight rectangle, highlighted scan,
         * mouse-selection rectangle.
         */
        protected class TopPanelLayer extends BaseImagePanelLayer
        {
            protected final Color HIGHLIGHT_BASE_COLOR =
                    new Color(Color.YELLOW.getRed(), Color.YELLOW.getGreen(),
                            Color.YELLOW.getBlue(), 192);

            public TopPanelLayer()
            {
                setName("Top");
            }

            public void draw(Graphics graphics,
                             int imageWidth, int imageHeight,
                             int minVisibleScan, int maxVisibleScan,
                             double minVisibleMz, double maxVisibleMz)
            {
                if (null != highlightRect)
                {
                    Rectangle viewRect = new Rectangle(highlightRect);

                    viewRect.y = getImageHeight() - (highlightRect.y + highlightRect.height) - 1;
                    //graphics.setColor(Color.RED);
                    graphics.setColor(adjustColor(TRANSLUCENT_CYAN));

                    int x0 = viewRect.x, x1 = viewRect.x + viewRect.width - 1;
                    int y0 = viewRect.y, y1 = viewRect.y + viewRect.height - 1;

                    graphics.drawLine(x0, y0, x1, y0);
                    graphics.drawLine(x0, y0+1, x0, y0+3);
                    graphics.drawLine(x1, y0+1, x1, y0+3);
                    graphics.drawLine(x0, y1, x1, y1);
                    graphics.drawLine(x0, y1-1, x0, y1-3);
                    graphics.drawLine(x1, y1-1, x1, y1-3);

                    if (null != highlightMasses)
                    {
                        graphics.setColor(adjustColor(HIGHLIGHT_BASE_COLOR));
                        Rectangle bounds = graphics.getClipBounds();
                        for (int i = 0; i < highlightMasses.length; i++)
                            graphics.drawLine((int)bounds.getMinX(),
                                    getImageHeight() - (int)highlightMasses[i] - 1,
                                    (int)bounds.getMaxX(),
                                    getImageHeight() - (int)highlightMasses[i] - 1);
                    }
                }


                if (-1 != highlightScan)
                {
                    int x = highlightScan;
                    if (x >= 0 && x < getImageWidth())
                    {
                        graphics.setColor(adjustColor(TRANSLUCENT_RED));
                        graphics.drawLine(x, 0, x, getImageHeight() - 1);
                    }
                }

                //Draw the mouse selection rectangle
                if (mouseSelectionRect != null)
                {
                    graphics.setColor(SELECTIONRECT_COLOR);
                    graphics.drawRect((int) mouseSelectionRect.getX(),
                            (int) (getImageHeight() - (mouseSelectionRect.getY() + mouseSelectionRect.getHeight())),
                            (int) Math.abs(mouseSelectionRect.getWidth()),
                            (int) Math.abs(mouseSelectionRect.getHeight()));
                    graphics.setColor(SELECTIONRECT_FILLCOLOR);
                    graphics.fillRect((int) mouseSelectionRect.getX(),
                            (int) (getImageHeight() - (mouseSelectionRect.getY() + mouseSelectionRect.getHeight())),
                            (int) Math.abs(mouseSelectionRect.getWidth()),
                            (int) Math.abs(mouseSelectionRect.getHeight()));
                }

                //todo: this really should be tied to FeatureSets or done away with
                Feature[] highlightedFeatures = (Feature[]) ApplicationContext.getProperty(SharedProperties.HIGHLIGHT_FEATURES);
                highlightFeatures(graphics, highlightedFeatures, Color.ORANGE,
                        minVisibleScan, maxVisibleScan, minVisibleMz, maxVisibleMz);
                Feature[] highlightedFeatures2 = (Feature[]) ApplicationContext.getProperty(SharedProperties.HIGHLIGHT_FEATURES2);
                highlightFeatures(graphics, highlightedFeatures2, Color.RED,
                        minVisibleScan, maxVisibleScan, minVisibleMz, maxVisibleMz);
                Feature[] highlightedFeatures3 = (Feature[]) ApplicationContext.getProperty(SharedProperties.HIGHLIGHT_FEATURES3);
                highlightFeatures(graphics, highlightedFeatures3, Color.BLUE,
                        minVisibleScan, maxVisibleScan, minVisibleMz, maxVisibleMz);
            }

        }

        protected ImagePanelLayer topPanelLayer = new TopPanelLayer();
        boolean mouseOver = false;

        //These rectangles will both be in correct scanindex, mz coordinates
        Point mouseSelectionStartingPoint = null;
        public Rectangle mouseSelectionRect = null;

        protected boolean button1IsPressed = false;

        /**
         * Return the height of the image.  This is important because it _doesn't_
         * change with the zoom factor
         * @return
         */
        public int getImageHeight()
        {
            if (getImage() == null)
                return super.getHeight();
            return getImage().getHeight(this);
        }

        /**
         * Return the width of the image.  This is important because it _doesn't_
         * change with the zoom factor
         * @return
         */
        public int getImageWidth()
        {
            if (getImage() == null)
                return super.getWidth();
            return getImage().getWidth(this);
        }

        /**
         * This DOES scale with the zoom factor
         * @return
         */
        public int getHeight()
        {
            return (int) getPreferredSize().getHeight();
        }

        /**
         * This DOES scale with the zoom factor
         * @return
         */
        public int getWidth()
        {
            return (int) getPreferredSize().getWidth();
        }

        /**
         * This DOES scale with the zoom factor, which is necessary in order
         * to tell the scrollpane how big it is
         * @return
         */
        public Dimension getPreferredSize() {
            if (_imagePanelSize == null)
                return super.getPreferredSize();
            return _imagePanelSize;
        }

        /**
         * If there's a current mouse-selected rectangle, zoom to it
         */
        public void zoomToSelection()
        {
            if (mouseSelectionRect == null)
                return;
            recenter((int) mouseSelectionRect.getCenterX(), (int) mouseSelectionRect.getCenterY());

            Rectangle viewportRect = getScrollPane().getViewportBorderBounds();

            double xZoomFactor = ((viewportRect.getWidth()-ROWHEADER_WIDTH) / getImageWidth()) / (Math.abs(mouseSelectionRect.getWidth() / getImageWidth()));
            double yZoomFactor = ((viewportRect.getHeight()-COLHEADER_HEIGHT) / getImageHeight()) / (Math.abs(mouseSelectionRect.getHeight() / getImageHeight()));
            double newZoomFactor = Math.min(xZoomFactor, yZoomFactor);
            transitionToZoomFactor(newZoomFactor);
        }

        public MSImagePanel()
        {
            super("MSImagePanel");

            _rightMousePopupMenu = new RightMousePopupMenu();
            _rightMousePopupMenu.setImagePanel(this);

            setComponentPopupMenu(_rightMousePopupMenu);

            //dragging listener for selection box creation
            this.addMouseMotionListener(new java.awt.event.MouseMotionAdapter()
            {
                public void mouseMoved(MouseEvent event)
                {
                    _currentMousePosition = event.getPoint();
                }

                //drag selection is done with left button only
                public void mouseDragged(MouseEvent event)
                {
                    if (button1IsPressed)
                    {
                        Point convertedEventPoint = convertMousePointToScanIndexMz(event.getPoint());
                        int topX = Math.min((int) convertedEventPoint.getX(),
                                (int) mouseSelectionStartingPoint.getX());
                        int topY = Math.min((int) convertedEventPoint.getY(),
                                (int) mouseSelectionStartingPoint.getY());
                        int width = (int) Math.abs(convertedEventPoint.getX() -
                                mouseSelectionStartingPoint.getX());
                        int height = (int) Math.abs(convertedEventPoint.getY() -
                                mouseSelectionStartingPoint.getY());

                        mouseSelectionRect = new Rectangle(topX, topY, width, height);

                        repaint();

                        _rightMousePopupMenu.enableSelectionItems();
                    }
                }
            }
            );

            this.addMouseListener(new java.awt.event.MouseAdapter()
            {
                public void mouseClicked(MouseEvent event)
                {

                    switch (event.getButton())
                    {
                        case MouseEvent.BUTTON1:
                        {
                            Point scanIndexMzPoint =
                                    convertMousePointToScanIndexMz(event.getPoint());
                            selectPoint((int) scanIndexMzPoint.getX(),
                                        (float) scanIndexMzPoint.getY());

                            mouseSelectionRect = null;
                            _rightMousePopupMenu.disableSelectionItems();
                        }
                        default:
                        {
                            //popup menu is taken care of separately.
                            //Only need to repaint to get rid of selection rectangle
                            repaint();
                        }
                    }
                }

                public void mouseEntered(MouseEvent event)
                {
                    mouseOver = true;
                }


                public void mouseExited(MouseEvent event)
                {
                    mouseOver = false;
                }

                public void mousePressed(MouseEvent event)
                {
                    if (event.getButton() == MouseEvent.BUTTON1)
                    {
                        button1IsPressed = true;
                        mouseSelectionStartingPoint = convertMousePointToScanIndexMz(event.getPoint());
//                            mouseSelectionRect = null;
//                            _rightMousePopupMenu.disableSelectionItems();
                    }
                }

                public void mouseReleased(MouseEvent event)
                {
                    if (event.getButton() == MouseEvent.BUTTON1)
                    {
                        button1IsPressed = false;
                    }
                }
            });

            // won't call getToolTipText() if we never setToolTipText().
            setToolTipText("");
        }

        public void saveMzXmlRegion()
        {
            if (mouseSelectionRect != null)
            {
                SavePartialMzxmlDialog saveDialog = new SavePartialMzxmlDialog();
                //have to convert scan indices to scan numbers
                int lowScanIndex = (int) mouseSelectionRect.getX();
                int highScanIndex = (int) (mouseSelectionRect.getX() + mouseSelectionRect.getWidth());
                MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);

                //trim the scan window, just in case selection is out of whack, because
                //otherwise we get nullpointerexceptions
                lowScanIndex = Math.max(lowScanIndex, 0);
                highScanIndex = Math.min(highScanIndex, run.getScanCount() - 1);

                int lowScanNumber = run.getScan(lowScanIndex).getNum();
                int highScanNumber = run.getScan(highScanIndex).getNum();
                saveDialog.setRegionToSave(lowScanNumber, highScanNumber,
                        (int) mouseSelectionRect.getY(),
                        (int) (mouseSelectionRect.getY() + mouseSelectionRect.getHeight()));
                saveDialog.setVisible(true);
            }
        }

        /**
         * Select a point in the run -- populate the detail and spectrum windows appropriately
         * @param scanNumber
         * @param mz
         */
        public void selectPoint(int scanNumber, float mz)
        {
            MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null == run || scanNumber < 0 || scanNumber >= run.getScanCount())
                return;
            MSRun.MSScan scan = run.getScan(scanNumber);
            ApplicationContext.setProperty(SharedProperties.MS_SCAN, scan);
            java.util.List<FeatureSet> featureSets =
                    (java.util.List<FeatureSet>)
                            ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
//            Feature feat = FeatureSet.hitTest(featureSets, scanNumber, mz, 3, 3);
            Feature feat = getNearestFeature(scanNumber, mz,
                    DEFAULT_FEATURE_SCAN_TOLERANCE,
                    DEFAULT_FEATURE_MZ_TOLERANCE);
            if (null != feat)
            {
                ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, feat);
                ApplicationContext.setProperty(SharedProperties.SELECTED, feat);
            }
            else
                ApplicationContext.setProperty(SharedProperties.SELECTED_POINT,
                                               new Spectrum.Peak(scan.getNum(), mz, 0.0F));
        }

        /**
         *  Figure out the right text to display, first accounting for the zoom
         * factor to figure out the _effective_ position of the mouse
         * @param event
         * @return
         */
        public String getToolTipText(MouseEvent event)
        {
            Point scanMzPoint = convertMousePointToScanIndexMz(event.getPoint());
            int scanIndex = (int) scanMzPoint.getX();
            int mz = (int) scanMzPoint.getY();

            Feature nearestFeature =
                    getNearestFeature(scanIndex, mz, DEFAULT_FEATURE_SCAN_TOLERANCE,
                                      DEFAULT_FEATURE_MZ_TOLERANCE);
            if (null != nearestFeature)
            {
                String message = "(" + nearestFeature.scan + "," + nearestFeature.mz + ")" +
                        " " + nearestFeature.getScanCount() + " scans " + nearestFeature.getScanFirst() + "-" + nearestFeature.getScanLast() + ", " + nearestFeature.intensity + "i " + nearestFeature.charge + "+";
                return message;
            }

            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);

            if (null != run && scanIndex >= 0 && scanIndex < run.getScanCount())
            {
                MSRun.MSScan scan = run.getScan(scanIndex);
                if (null != scan)
                {
                    String message = "(" + scan.getNum() + "," + mz + ") ";
                    return message;
                }
            }
            return null;
        }

        /**
         * Find the nearest feature to the (scan, mz) position
         *
         * @param scanIndex
         * @param mz
         * @return
         */
        public Feature getNearestFeature(int scanIndex, float mz,
                                         int scanTolerance, float mzTolerance)
        {
            Feature nearestFeature = null;
            if (null != featureSets)
            {
                java.util.List<Feature> nearbyFeatures =
                        (java.util.List<Feature>)
                                _allDisplayedFeaturesScanIndexMzTree.getPoints(
                                        scanIndex-scanTolerance, mz-mzTolerance,
                                        scanIndex+scanTolerance,  mz+mzTolerance);
                if (nearbyFeatures.size() == 0)
                    return null;
                double minDistance = Double.MAX_VALUE;
                Feature result = null;
                for (Feature feature : nearbyFeatures)
                {
                    int featureScanIndex = _run.getIndexForScanNum(feature.scan);
                    scanIndex = scanIndex < 0 ? -(scanIndex+1) : scanIndex;
                    double distance =
                            Math.sqrt(Math.pow(featureScanIndex - scanIndex, 2) *
                                    Math.pow(feature.mz - mz, 2));
                    if (distance < minDistance)
                    {
                        result = feature;
                        minDistance = distance;
                    }
                }
                return result;
            }
            return nearestFeature;
        }


        Color TRANSLUCENT_RED = new Color(1F, 0F, 0F, 0.5F);
        Color TRANSLUCENT_CYAN = new Color(0F, 1F, 1F, 0.5F);
        Color SELECTIONRECT_COLOR = new Color(0.4F, 1F, 0.4F, 0.8F);
        Color SELECTIONRECT_FILLCOLOR = new Color(0.3F, .9F, 0.3F, 0.1F);

        /**
         * Paint, scaling appropriately according to zoom
         * @param graphics
         */
        public void paint(Graphics graphics)
        {
            paint(graphics, true);
        }

        /**
         * Paint the image window, being awfully darn careful about things
         * that need to get adjusted for zoom level
         * @param graphics
         * @param scaleWithZoom If true, scale according to zoom. Set to false
         * for saving an image to a file
         */
        public void paint(Graphics graphics, boolean scaleWithZoom)
        {
            if (scaleWithZoom)
                ((Graphics2D) (graphics)).scale(_zoomFactor, _zoomFactor);

            super.paint(graphics);

            if (null == _run)
                return;



            Rectangle bounds = graphics.getClipBounds();
            int minScan = null == bounds ? 0 : (int)bounds.getX();
            int maxScan = null == bounds ? Integer.MAX_VALUE : (int)(minScan + bounds.getWidth());
            double minMz = null == bounds ? 0 : getImageHeight() - bounds.getY() - bounds.getHeight();
            double maxMz = null == bounds ? Integer.MAX_VALUE : getImageHeight() - bounds.getY();

            for (ImagePanelLayer imagePanelLayer : imagePanelLayers)
            {
                imagePanelLayer.draw(graphics, getImageWidth(), getImageHeight(), minScan, maxScan, minMz, maxMz);
            }


            topPanelLayer.draw(graphics, getImageWidth(), getImageHeight(), minScan, maxScan, minMz, maxMz);


        }

        /**
         * draw circles around highlighted features in a given color
         * @param graphics
         * @param highlightedFeatures
         * @param color
         * @param minScan
         * @param maxScan
         * @param minMz
         * @param maxMz
         */
        protected void highlightFeatures(Graphics graphics, Feature[] highlightedFeatures, Color color,
                                         int minScan, int maxScan, double minMz, double maxMz
        )
        {
            if (null != highlightedFeatures)
            {
                //Clone so we don't have to reset line width
                Graphics2D g = (Graphics2D)graphics.create();
                g.setColor(color);
                g.setStroke(new BasicStroke(2.0f));
                for (int i = 0; i < highlightedFeatures.length; i++)
                {
                    Feature feature = highlightedFeatures[i];
                    int scanIndex = _run.getIndexForScanNum(feature.scan);
                    scanIndex = scanIndex < 0 ? -(scanIndex+1) : scanIndex;
                    if (scanIndex < minScan || scanIndex > maxScan || feature.mz < minMz || feature.mz > maxMz)
                        continue;

                    g.drawOval(scanIndex - 8, getImageHeight() - (int)feature.mz - 9, 16, 16);
                }
            }
        }
    }


    public static class SaveImageAction extends AbstractAction
    {
        public SaveImageAction()
        {
            super("Save Image...");
        }

        public void actionPerformed(ActionEvent e)
        {
            MSImageComponent comp = (MSImageComponent)ApplicationContext.getProperty("MSImageComponent");
            if (null == comp)
                return;

            String name = "";
            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null != run)
            {
                name = run.getFile().getPath();
                if (name.toLowerCase().endsWith(".mzxml"))
                    name = name.substring(0, name.length() - ".mzxml".length());
                name += ".png";
            }

            WorkbenchFileChooser chooser = new WorkbenchFileChooser();
            if (0 < name.length())
                chooser.setSelectedFile(new File(name));
            int res = chooser.showSaveDialog(ApplicationContext.getFrame());
            if (res != JFileChooser.APPROVE_OPTION)
                return;
            File f = chooser.getSelectedFile();
            if (null == f)
                return;

            comp.saveImage(f);
        }
    }

    /**
     * Save the image to a file, no rescaling
     * @param f
     */
    public void saveImage(File f)
    {
        saveImage(f, Integer.MAX_VALUE, Integer.MAX_VALUE, true);
    }

    /**
     * Saves the image to a file, rescaling the image if necessary, to make it fit both the maxWidth and
     * maxHeight constraints.
     *
     * DOES NOT distort the image
     * @param f
     * @param maxWidth
     * @param maxHeight
     */
    public void saveImage(File f, int maxWidth, int maxHeight, boolean includeTIC)
    {
        BufferedImage checkedImage = (BufferedImage) imagePanel.getImage();
        String ext = f.getName().substring(f.getName().lastIndexOf('.')+1).toLowerCase();

        String formats[] = ImageIO.getWriterFormatNames();
        setImage(checkedImage);

        String format = null;

        for (int i = 0; i < formats.length && null == format; i++)
        {
            if (formats[i].toLowerCase().equals(ext))
                format = formats[i];
        }
        if (null == format)
        {
            f = new File(f.getPath() + ".png");
            format = "PNG";
        }


        BufferedImage imageBW = (BufferedImage)this.imagePanel.getImage();
        int width = imageBW.getWidth();
        int height = imageBW.getHeight();
        Image copy = null;

        if (null == featureSets)
        {
            copy = imageBW;
//			new BufferedImage(w, h, imageBW.getType());
//			imageBW.copyData(copy.getRaster());
        }
        else
        {
            // need to translate to color
            copy = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);

            Raster rBW = imageBW.getRaster();
            WritableRaster rCopy = ((BufferedImage) copy).getRaster();
            // could try to do this all at once, but it requires very large sample array
            // do row by row
            int[] samples = new int[width];
            for (int y = 0 ; y<height ; y++)
            {
                rBW.getSamples(0, y, width, 1, 0, samples);
                rCopy.setSamples(0, y, width, 1, 0, samples);
                rCopy.setSamples(0, y, width, 1, 1, samples);
                rCopy.setSamples(0, y, width, 1, 2, samples);
            }

            // we don't want to draw highlights, save and restore these variables
            Rectangle r = highlightRect;
            int s = highlightScan;
            highlightRect = null;
            highlightScan = -1;

            // draw features
            imagePanel.paint(copy.getGraphics(), false);
            highlightRect = r;
            highlightScan = s;
        }

        //cut out the Total Ion Chromatogram
        if (!includeTIC)
        {
            float minMz = _run.getMzRange().min;

            CropImageFilter cif = new CropImageFilter(0, 0, width, (int) (height-minMz));
            FilteredImageSource fis = new FilteredImageSource(copy.getSource(), cif);
            Image croppedImage = imagePanel.createImage(fis);

            int w = croppedImage.getWidth(null), h = croppedImage.getHeight(null);
            int type = BufferedImage.TYPE_INT_RGB;  // many options
            copy = new BufferedImage(w, h, type);
            Graphics2D g2 = ((BufferedImage) copy).createGraphics();
            g2.drawImage(croppedImage, 0, 0, null);
            g2.dispose();
        }

        if (maxWidth < width || maxHeight < height)
        {
            if (maxWidth > width) maxWidth = width;
            if (maxHeight > height) maxHeight = height;
            ApplicationContext.setMessage("Rescaling image...");
            float imageWidthHeightRatio = (float) width / (float) height;
            float specifiedWidthHeightRatio = (float) maxWidth / (float) maxHeight;
            int newWidth, newHeight;
            if (imageWidthHeightRatio > specifiedWidthHeightRatio)
            {
                newWidth = maxWidth;
                newHeight = (int) (newWidth / imageWidthHeightRatio);
            }
            else
            {
                newHeight = maxHeight;
                newWidth = (int) (imageWidthHeightRatio * (float) newHeight);
            }
//            imageToWrite = copy.getScaledInstance(newWidth, newHeight, Image.SCALE_DEFAULT);
//            copy = new BufferedImage(imageToWrite);

            BufferedImage scaledImage = new BufferedImage(newWidth,
                    newHeight, BufferedImage.TYPE_INT_RGB);
            Graphics2D graphics2D = scaledImage.createGraphics();
            graphics2D.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            graphics2D.drawImage(copy, 0, 0, newWidth, newHeight, null);

            copy = scaledImage;
        }

        try
        {
            ApplicationContext.setMessage("Writing image " + f.getPath() + "...");
            IntensityPlot.writePlot(f, copy, format);
            ApplicationContext.setMessage("Successfully saved file " + f.getPath());
        }
        catch (IOException x)
        {
            ApplicationContext.errorMessage(x.getMessage(), null);
        }
    }

    /**
     * Our context menu.  There is a set of standard items that always comes up, related
     * to zoom and mzXML-saving.  It is also possible to add feature-specific items that
     * come up when the user right-clicks close to a feature.
     *
     * This is done by analyzing the Extra Information Types that a feature has.  Each
     * Extra Information Type may optionally provide MenuItems for a feature that will
     * show up here.  E.g., an item that kicks off a Google search for the peptide.
     */
    protected class RightMousePopupMenu extends JPopupMenu
    {
        //zoom level items
        JMenuItem item100, item25, item50, item150, item200, item300, item400;
        //mzxml-saving item
        JMenuItem itemSaveMzXml;
        //selection-zooming item
        JMenuItem itemZoomToSelection;

        Feature selectedFeature = null;
        java.util.List<Component> featureMenuItems = new ArrayList<Component>();

        protected MSImagePanel _imagePanel = null;

        //for easy reference
        protected ArrayList<JMenuItem> _menuItems;

        public RightMousePopupMenu()
        {
            super();

            item100 = new JMenuItem("100%");
            item25 = new JMenuItem("25%");
            item50 = new JMenuItem("50%");
            item150 = new JMenuItem("150%");
            item200 = new JMenuItem("200%");
            item300 = new JMenuItem("300%");
            item400 = new JMenuItem("400%");
            itemSaveMzXml = new JMenuItem(TextProvider.getText("SAVE_MZXML_REGION"));
            itemZoomToSelection = new JMenuItem(TextProvider.getText("ZOOM_TO_SELECTION"));

            _menuItems = new ArrayList<JMenuItem>();
            _menuItems.add(item100);
            _menuItems.add(item25);
            _menuItems.add(item50);
            _menuItems.add(item150);
            _menuItems.add(item200);
            _menuItems.add(item300);
            _menuItems.add(item400);
            _menuItems.add(itemSaveMzXml);
            _menuItems.add(itemZoomToSelection);

            item100.setActionCommand("100");
            item25.setActionCommand("25");
            item50.setActionCommand("50");
            item150.setActionCommand("150");
            item200.setActionCommand("200");
            item300.setActionCommand("300");
            item400.setActionCommand("400");
            itemSaveMzXml.setActionCommand("save_mzxml");
            itemZoomToSelection.setActionCommand("zoom_to_selection");

            ListenerHelper helper = new ListenerHelper(this);
            for (int i=0; i<_menuItems.size(); i++)
            {
                helper.addListener(_menuItems.get(i),"menuItem_actionPerformed");
            }

            //layout
            add(item100);
            add(new JPopupMenu.Separator());
            add(item150);
            add(item200);
            add(item300);
            add(item400);
            add(new JPopupMenu.Separator());
            add(item25);
            add(item50);
            add(itemZoomToSelection);
            add(new JPopupMenu.Separator());
            add(itemSaveMzXml);
            

            //initially disable 100% zoom, because that's what we're at
            item100.setEnabled(false);

            //initially disable saving mzXml and zoom to region, since no region selected
            itemSaveMzXml.setEnabled(false);
            itemZoomToSelection.setEnabled(false);

            //for figuring out mouse location
            this.addPopupMenuListener(new RightPopupListener());
        }

        public void setSelectedFeature(Feature feature)
        {
            selectedFeature = feature;
        }

        public void enableSelectionItems()
        {
            itemSaveMzXml.setEnabled(true);
            itemZoomToSelection.setEnabled(true);
        }

        public void disableSelectionItems()
        {
            itemSaveMzXml.setEnabled(false);
            itemZoomToSelection.setEnabled(false);
        }

        public void setImagePanel(MSImagePanel imagePanel)
        {
            _imagePanel = imagePanel;
        }

        /**
         * Remove all special items related to selected feature
         */
        public void removeFeatureSpecificItems()
        {
            _menuItems.removeAll(featureMenuItems);
            for (Component component : featureMenuItems)
                remove(component);
            featureMenuItems = new ArrayList<Component>();
        }

        /**
         *  Add items related to the selected feature, by querying all extra info
         * types associated with the feature
         */
        public void addFeatureSpecificItems()
        {

            for (FeatureExtraInformationDef extraInfoDef :
                    selectedFeature.determineExtraInformationTypes())
            {
                java.util.List<JMenuItem> thisDefMenuItems =
                        extraInfoDef.createPopupMenuItems(selectedFeature);
                if (thisDefMenuItems != null && thisDefMenuItems.size() > 0)
                {
                    Component separator = new JPopupMenu.Separator();
                    add(separator);
                    featureMenuItems.add(separator);
                    for (JMenuItem menuItem : thisDefMenuItems)
                    {
                        featureMenuItems.add(menuItem);
                        add(menuItem);
                    }
                }
            }
        }

        public void menuItem_actionPerformed(ActionEvent event)
        {
            if ("save_mzxml".equals(event.getActionCommand()))
            {
                _imagePanel.saveMzXmlRegion();
            }
            else if  ("zoom_to_selection".equals(event.getActionCommand()))
            {
                _imagePanel.zoomToSelection();
            }
            else
            {
                int zoomFactorPercent = Integer.parseInt(event.getActionCommand());
                transitionToZoomFactor(zoomFactorPercent * .01);
            }
        }

        /**
         * Set the correct disabled item representing the current zoom level
         */
        public void updateZoomList()
        {
            String zoomFactorPercent = "" + (int)(_zoomFactor * 100);
            for (int i=0; i<_menuItems.size(); i++)
            {
                JMenuItem currentItem = _menuItems.get(i);

                if (zoomFactorPercent.equals(currentItem.getActionCommand()))
                    currentItem.setEnabled(false);
                else
                    currentItem.setEnabled(true);
            }
        }
    }


    /**
     * This is necessary because the automatic event listening supplied by
     * setComponentPopupMenu() cuts in front of the regular mouseclick listening
     * we're doing.  Meaning that the only way to determine the location of
     * the mouseclick is by asking the popup menu itself.
     *
     * Also, here, we handle creation of special menu items related to the selected feature
     */
    protected class RightPopupListener implements PopupMenuListener
    {
        /**
         * Can't ask for the location here because the menu isn't visible yet
         * @param e
         */
        public void popupMenuWillBecomeVisible(PopupMenuEvent e)
        {
            Point convertedMousePoint = convertMousePointToScanIndexMz(_currentMousePosition);
            Feature nearestFeature =
                    imagePanel.getNearestFeature((int) convertedMousePoint.getX(),
                            (float) convertedMousePoint.getY(),
                            DEFAULT_FEATURE_SCAN_TOLERANCE,
                            DEFAULT_FEATURE_MZ_TOLERANCE);
            _rightMousePopupMenu.setSelectedFeature(nearestFeature);
            _rightMousePopupMenu.removeFeatureSpecificItems();
            if (nearestFeature != null)
            {
                _rightMousePopupMenu.addFeatureSpecificItems();
            }
        }

        /**
         * You'd think this would get kicked off _after_ the event processing
         * that we use to handle zoom, and therefore it'd be useless for storing
         * the mouseclick location. You'd be wrong.  For some reason this happens
         * first.  And a good thing, too!
         * @param e
         */
        public void popupMenuWillBecomeInvisible(PopupMenuEvent e)
        {
            //get the absolute screen location
            int x = (int) _rightMousePopupMenu.getLocationOnScreen().getX();
            int y = (int) _rightMousePopupMenu.getLocationOnScreen().getY();

            //subtract out the upper-left-corner of the image panel
            x = x - (int) imagePanel.getLocationOnScreen().getX();
            y = y - (int) imagePanel.getLocationOnScreen().getY();

            //account for zoom
            x = (int) (x / _zoomFactor);
            y = (int) (y / _zoomFactor);

            _lastMouseRightClick = new Point(x,y);
        }

        /**
         * Nothin'
         * @param e
         */
        public void popupMenuCanceled(PopupMenuEvent e)
        {
        }
    }
}

