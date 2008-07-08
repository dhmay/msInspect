/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

package org.fhcrc.cpl.viewer.gui.util;

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItemSource;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.LegendItem;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.*;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.ms2.Fractionation2DUtilities;
import org.fhcrc.cpl.viewer.gui.ImagePanel;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.labkey.common.tools.ApplicationContext;
import org.labkey.common.tools.TextProvider;
import org.labkey.common.tools.Rounder;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseAdapter;
import java.awt.image.BufferedImage;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * PanelWithChart implementation to make it easy to put out Line Charts
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithBlindImageChart extends PanelWithChart
{
    static Logger _log = Logger.getLogger(PanelWithBlindImageChart.class);

    protected BufferedImage image = null;
    protected ImagePanel imagePanel = null;
    protected JPopupMenu popupMenu = null;

    protected JScrollPane scrollPane = null;


    public PanelWithBlindImageChart()
    {
        super();
        scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        add(scrollPane);
    }

    public PanelWithBlindImageChart(String dataSetName)
    {
        this();
        setName(dataSetName);
    }

    public void setPreferredSize(Dimension dimension)
    {
        scrollPane.setPreferredSize(dimension);       
    }


    public PanelWithBlindImageChart(File imageFile, String name)
    {
        this(name);
        try
        {
            setImage(imageFile);
        }
        catch (IOException e)
        {
            ApplicationContext.errorMessage("Error while displaying chart image " +
                                            imageFile.getAbsolutePath(),e);
        }
    }

    public PanelWithBlindImageChart(BufferedImage image, String name)
    {
        this(name);
        setImage(image);
    }

    public void setImage(BufferedImage image)
    {
        setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));

        imagePanel = new ImagePanel(image);
        imagePanel.setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));
        scrollPane.setViewportView(imagePanel);

        popupMenu = new JPopupMenu();
        JMenuItem saveImageMenuItem = new JMenuItem("Save Image");
        saveImageMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    WorkbenchFileChooser wfc = new WorkbenchFileChooser();
                    int chooserStatus = wfc.showOpenDialog(PanelWithBlindImageChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = wfc.getSelectedFile();
                    try
                    {
                        saveChartToImageFile(outFile);
                    }
                    catch (IOException e)
                    {
                        ApplicationContext.errorMessage("Failed to save image to file " + outFile.getAbsolutePath(), e);
                    }
                }
            }
        );
        popupMenu.add(saveImageMenuItem);

        PopupListener myPopupListener = new PopupListener(popupMenu);
        scrollPane.addMouseListener(myPopupListener);
        imagePanel.addMouseListener(myPopupListener);
        this.addMouseListener(myPopupListener);
    }

    protected class PopupListener extends MouseAdapter
    {
        protected JPopupMenu popup = null;

        public PopupListener(JPopupMenu popup)
        {
            this.popup = popup;
        }

        public void mousePressed(MouseEvent e) {
            maybeShowPopup(e);
        }

        public void mouseReleased(MouseEvent e) {
            maybeShowPopup(e);
        }

        private void maybeShowPopup(MouseEvent e) {
            if (e.isPopupTrigger()) {
                popup.show(e.getComponent(),
                        e.getX(), e.getY());
            }
        }
    }




    public void setImage(File imageFile)
            throws IOException
    {
        image = ImageIO.read(imageFile);
        _log.debug("Loaded image " + imageFile.getAbsolutePath() + ", size: " +
                image.getWidth() + ", " + image.getHeight());
        setImage(image);
    }


    public void addItemToPopupMenu(JMenuItem item)
    {
    }

    public void addSeparatorToPopupMenu()
    {
    }



    /**
     * Don't add the usual "save data" items to the popup
     */
    protected void initPopupMenu()
    {
        JMenuItem saveImageMenuItem = new JMenuItem("Save Chart");
        saveImageMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    WorkbenchFileChooser wfc = new WorkbenchFileChooser();
                    int chooserStatus = wfc.showOpenDialog(PanelWithBlindImageChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = wfc.getSelectedFile();
                    try
                    {
                        saveChartToImageFile(outFile);
                    }
                    catch (Exception e)
                    {
                        ApplicationContext.errorMessage("Error saving image file " +
                                outFile.getAbsolutePath(),e);
                    }
                }
            }
        );
        addItemToPopupMenu(saveImageMenuItem);
    }

    public void saveChartToImageFile(File outFile) throws IOException
    {
        ImageIO.write(image,"png",outFile);
    }

    /**
     * Do nothing
     * @param outFile
     * @param delimiter
     */
    protected void saveChartDataToFile(File outFile, String delimiter)
    {

    }
}
