/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
import java.util.List;
import java.util.ArrayList;
import java.util.Date;
import java.awt.event.ComponentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;

/**
 * display multiple charts in a panel.  Supports saving of all charts at once, to a directory
 *
 * Only a tabbed display is supported at the moment
 */
public class MultiChartDisplayPanel extends JPanel
{
    protected static Logger _log =
            Logger.getLogger(MultiChartDisplayPanel.class);


    public static final int DEFAULT_WIDTH = 800;
    public static final int DEFAULT_HEIGHT = 800;

    public static MultiChartDisplayPanel _singletonInstance = null;

    public static boolean _hiddenMode = false;

    protected List<PanelWithChart> chartPanels = null;

    public static final int DISPLAY_STYLE_TABBED = 0;

    protected int displayStyle = DISPLAY_STYLE_TABBED;

    protected JComponent componentHoldingChildren = null;

    protected long lastResizeTime = 0;

    protected JDialog dialog = null;


    public MultiChartDisplayPanel()
    {
        super();
        chartPanels = new ArrayList<PanelWithChart>();
        setSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
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
        getComponentHoldingChildren().add(newPanel);
        newPanel.addSeparatorToPopupMenu();
        newPanel.addItemToPopupMenu(createSaveAllMenuItem());
    }

    public void removeAllCharts()
    {
        for (PanelWithChart panel : chartPanels)
        {
            getComponentHoldingChildren().remove(panel);
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
        getComponentHoldingChildren().remove(panel);
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
     * Abstraction for JTabbedPane and whatever else we implement
     * @return
     */
    protected JComponent getComponentHoldingChildren()
    {
        if (componentHoldingChildren == null)
        {
            switch (displayStyle)
            {
                case DISPLAY_STYLE_TABBED:
                    JTabbedPane tabbedPane = new JTabbedPane();
                    componentHoldingChildren = tabbedPane;
                    tabbedPane.addChangeListener(new TabFocusListener());
                    break;
            }
            add(componentHoldingChildren);
        }
        return componentHoldingChildren;
    }


    /**
     * Resize a chart panel appropriately
     * @param pwc
     */
    protected void resizeChartPanel(PanelWithChart pwc)
    {



//
//        int newPanelWidth = componentForSizing.getSize().width - 20;
//        int newPanelHeight = componentForSizing.getSize().height - 50;
//
//        pwc.setPreferredSize(new Dimension(newPanelWidth, newPanelHeight));
//        pwc.updateUI();
    }

    protected void resizeChartInFocus()
    {
        long currentTime = new Date().getTime();
        long timeDiff = currentTime - lastResizeTime;

        if (timeDiff > 200)
        {
            switch (displayStyle)
            {
                case DISPLAY_STYLE_TABBED:
                    JTabbedPane tabbedPane = (JTabbedPane) componentHoldingChildren;
                    if (tabbedPane != null)
                    {
                        PanelWithChart selectedChart = (PanelWithChart) tabbedPane.getSelectedComponent();
                        if (selectedChart != null)
                        {
                            Point chartLocation = selectedChart.getLocation();
                            int newChartWidth = getWidth() - 2 * (5 + chartLocation.x - getLocation().x);
                            int newChartHeight = getHeight() - 10 - (chartLocation.y - getLocation().y);
                            selectedChart.setPreferredSize(new Dimension(newChartWidth, newChartHeight));

                            selectedChart.updateUI();
                        }
                    }
            }
        }
        lastResizeTime = currentTime;
    }

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
//            synchronized (chartPanels)
//            {
//                for (PanelWithChart pwc : chartPanels)
//                {
//                    resizeChartPanel(pwc);
//                }
//            }

              resizeChartInFocus();
        }

        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }

    public static MultiChartDisplayPanel getSingletonInstance()
    {
        if (_singletonInstance == null)
            _singletonInstance = new MultiChartDisplayPanel();
        return _singletonInstance;
    }

    public JDialog displayInDialog()
    {
        return displayInDialog("Charts");
    }

    public JDialog displayInDialog(String dialogTitle)
    {
        dialog = new JDialog();
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setTitle(dialogTitle);
        dialog.setSize(new Dimension(getWidth() + 10, getHeight() + 25));
        dialog.add(this);
        if (!_hiddenMode)
            dialog.setVisible(true);
                                           
        return dialog;
    }

    public JDialog getDialog()
    {
        return dialog;
    }

    public static void addAndDisplayChartOnSingleton(PanelWithChart pwc)
    {
        if (_singletonInstance == null || getSingletonInstance().getDialog() == null)
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
}
