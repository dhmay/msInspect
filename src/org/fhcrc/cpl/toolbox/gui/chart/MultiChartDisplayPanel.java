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

package org.fhcrc.cpl.toolbox.gui.chart;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.imageio.ImageIO;
import java.util.List;
import java.util.ArrayList;
import java.util.Date;
import java.awt.event.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;

/**
 * display multiple charts in a panel.  Supports saving of all charts at once, to a directory
 *
 * Only a tabbed display is supported at the moment
 */
public abstract class MultiChartDisplayPanel extends JPanel
{
    protected static Logger _log =
            Logger.getLogger(MultiChartDisplayPanel.class);


    public static final int DEFAULT_WIDTH = 800;
    public static final int DEFAULT_HEIGHT = 800;

    public static MultiChartDisplayPanel _singletonInstance = null;

    public static boolean _hiddenMode = false;

    protected final List<PanelWithChart> chartPanels = new ArrayList<PanelWithChart>();

    protected long lastResizeTime = 0;

    protected JDialog dialog = null;

    //delay in resizing chart, after component resize.  For performance
    protected int resizeDelayMS = 200;


    public MultiChartDisplayPanel()
    {
        super();
        //setSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
        addComponentListener(new ResizeListener());
    }

    public static void destroySingleton()
    {
        _singletonInstance = null;
    }

    /**
     * Adds a new menu item to the PanelWithChart for saving of all charts
     * @param newPanel
     */
    public void addChartPanel(PanelWithChart newPanel)
    {
        synchronized (chartPanels)
        {
            chartPanels.add(newPanel);
        }
        resizeChartPanel(newPanel);
        addChartPanelToGUI(newPanel);
        newPanel.addSeparatorToPopupMenu();
        newPanel.addItemToPopupMenu(createSaveAllMenuItem());
    }

    protected abstract void addChartPanelToGUI(PanelWithChart newPanel);

    protected abstract void removeChartPanelFromGUI(PanelWithChart panel);



    public void removeAllCharts()
    {
        for (PanelWithChart panel : chartPanels)
        {
            removeChartPanelFromGUI(panel);
        }
        synchronized(chartPanels)
        {
            chartPanels.removeAll(chartPanels);                
        }
    }

    public int getNumCharts()
    {
        return chartPanels.size();
    }

    public void removeChartPanel(PanelWithChart panel)
    {
        synchronized (chartPanels)
        {
            chartPanels.remove(panel);
        }
        removeChartPanelFromGUI(panel);
    }

    /**
     * Creates a "save all charts" menu item
     * @return
     */
    protected JMenuItem createSaveAllMenuItem()
    {
        JMenuItem saveAllMenuItem = new JMenuItem("Save All Charts");
        saveAllMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JFileChooser fc = new JFileChooser();
                    fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                    int chooserStatus = fc.showOpenDialog(MultiChartDisplayPanel.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outDir = fc.getSelectedFile();
                    try
                    {
                        saveAllChartsToFiles(outDir);
                        ApplicationContext.infoMessage("All charts saved to directory " + outDir.getAbsolutePath());
                    }
                    catch (IOException e)
                    {
                        ApplicationContext.errorMessage("Failed to save charts to directory " +
                                outDir.getAbsolutePath(),e);
                    }
                }
            }
        );
        return saveAllMenuItem;
    }

    /**
     * Resize a chart panel appropriately
     * @param pwc
     */
    protected void resizeChartPanel(PanelWithChart pwc)
    {
    }

    protected abstract void resizeChartInFocus();

    protected class TabFocusListener implements ChangeListener
    {
        public void stateChanged(ChangeEvent e)
        {
            resizeChartInFocus();
        }
    }

    /**
     * Resize all chart panels when this panel is resized
     */
    protected class ResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
              resizeChartInFocus();
        }

        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }

    public static MultiChartDisplayPanel getSingletonInstance()
    {
        if (_singletonInstance == null)
            _singletonInstance = new TabbedMultiChartDisplayPanel();
        return _singletonInstance;
    }

    /**
     * Returns the old singleton instance, if any;
     * @return
     */
    public static MultiChartDisplayPanel createNewSingletonInstance()
    {
        MultiChartDisplayPanel result = _singletonInstance;
        _singletonInstance = new TabbedMultiChartDisplayPanel();
        return result;
    }

    public JDialog displayInDialog()
    {
        return displayInDialog("Charts", _hiddenMode);
    }

    
    public JDialog displayInDialog(String dialogTitle, boolean hidden)
    {
        setSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
        
        dialog = new JDialog();
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setTitle(dialogTitle);
        dialog.setSize(new Dimension(getWidth() + 10, getHeight() + 25));
        dialog.add(this);
        if (!hidden)
            dialog.setVisible(true);
                                           
        return dialog;
    }

    public JDialog getDialog()
    {
        return dialog;
    }

    public static void addAndDisplayChartOnSingleton(PanelWithChart pwc)
    {
        if (_singletonInstance == null || getSingletonInstance().getDialog() == null || !getSingletonInstance().getDialog().isDisplayable())
        {
            getSingletonInstance().displayInDialog();
        }
        getSingletonInstance().addChartPanel(pwc);
    }

    /**
     * Saves all charts to files in a given directory
     * @param outputDir
     * @throws IOException
     */
    public void saveAllChartsToFiles(File outputDir)
            throws IOException
    {
        if (!outputDir.isDirectory())
            throw new IOException("saveAllChartsToFiles: File " + outputDir.getAbsolutePath() + " is not a directory");        
        if (!outputDir.exists() || !outputDir.canWrite())
            throw new IOException("saveAllChartsToFiles: nonexistent or unwriteable directory " + outputDir.getAbsolutePath());
        for (PanelWithChart pwc : getChartPanels())
            pwc.saveChartToImageFile(new File(outputDir, createChartFileName(pwc.getName())));
    }

    /**
     * vertical orientation
     * @param charts
     * @return
     */
    public static BufferedImage createImageForAllCharts(List<PanelWithChart> charts)
    {
        int maxWidth = 0;
        int totalHeight = 0;
        List<Integer> heights = new ArrayList<Integer>();

        for (PanelWithChart pwc : charts)
        {
            int currentWidth = pwc.getWidth();
            if (currentWidth == 0)
                currentWidth = PanelWithChart.DEFAULT_WIDTH_FOR_IMAGE_FILE;
            if (pwc.getWidth() > maxWidth)
                maxWidth = pwc.getWidth();

            int currentHeight = pwc.getHeight();
            if (currentHeight == 0)
                currentHeight = PanelWithChart.DEFAULT_HEIGHT_FOR_IMAGE_FILE;
            heights.add(currentHeight);
            totalHeight += currentHeight;
        }

        BufferedImage bigImage = new BufferedImage(maxWidth, totalHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = bigImage.createGraphics();

        int currentHeightCum = 0;
        for (int i=0; i<charts.size(); i++)
        {
            PanelWithChart pwc = charts.get(i);
            Image thisChartImage = pwc.createImage();
            g.drawImage(thisChartImage, 0, currentHeightCum, null);
            currentHeightCum += heights.get(i);
        }

        g.dispose();

        return bigImage;
    }

    /**
     * Dump all charts to one file, vertically
     * @param outFile
     */
    public void saveAllChartsToFile(File outFile)
            throws IOException
    {
        BufferedImage bigImage = createImageForAllCharts(getChartPanels());
        ImageIO.write(bigImage,"png",outFile);
    }

    /**
     * Create a safe name for a chart
     * @param chartTitle
     * @return
     */
    protected String createChartFileName(String chartTitle)
    {
        String result = chartTitle + ".png";
        try
        {
            result = URLEncoder.encode(chartTitle, "US-ASCII") + ".png";
        }
        catch (UnsupportedEncodingException e)
        {
            ApplicationContext.infoMessage("Error encoding filename for chart named " + chartTitle);
        }
        return result;
    }


    public List<PanelWithChart> getChartPanels()
    {
        return chartPanels;
    }

    public static boolean isHiddenMode()
    {
        return _hiddenMode;
    }

    public static void setHiddenMode(boolean _hiddenMode)
    {
        MultiChartDisplayPanel._hiddenMode = _hiddenMode;

        if (_singletonInstance != null && _singletonInstance.getDialog() != null)
        {
            _singletonInstance.getDialog().setVisible(!_hiddenMode);
        }
    }

    public int getResizeDelayMS()
    {
        return resizeDelayMS;
    }

    public void setResizeDelayMS(int resizeDelayMS)
    {
        this.resizeDelayMS = resizeDelayMS;
    }
}
