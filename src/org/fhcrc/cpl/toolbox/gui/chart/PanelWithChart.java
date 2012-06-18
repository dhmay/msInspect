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

import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.RectangleInsets;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.imageio.ImageIO;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * Generic render-a-chart-in-a-panel class.  Can export the chart data to
 * .csv, .tsv files
 */
public class PanelWithChart extends JPanel
{
    protected static Logger _log =
            Logger.getLogger(PanelWithChart.class);

    protected Plot _plot = null;
    protected JFreeChart _chart = null;
    protected ChartPanel _chartPanel = null;
    boolean showLegend = false;


    public static final int DEFAULT_WIDTH_FOR_IMAGE_FILE = 800;
    public static final int DEFAULT_HEIGHT_FOR_IMAGE_FILE = 600;


    public PanelWithChart()
    {
        super();
    }

    public PanelWithChart(String name)
    {
        this();
        setName(name);
    }

    public PanelWithChart(Plot plot)
    {
        this();
        init(plot);
    }

    public PanelWithChart(JFreeChart chart)
    {
        init(chart);
    }

    protected void init(JFreeChart chart)
    {
        if (_chartPanel != null)
            remove(_chartPanel);
        _chart = chart;
        _plot = chart.getPlot();
        //dhmay changing the useBuffer arg to true, 20090915, with jfree 1.0.13.  Much prettier, better performance
        _chartPanel = new ChartPanel(_chart, true);
        _chartPanel.setDisplayToolTips(true);
        add(_chartPanel);

        if (_plot instanceof XYPlot)
        {
            //only add .tsv and .csv save options if this is an XYPlot.
            //Otherwise, no way to get at the data generically
            initPopupMenu();
            
            //dhmay adding 2009/09/15.  As of jfree 1.0.13, several defaults changed annoyingly
            ((XYPlot)_plot).setDomainGridlinePaint(Color.LIGHT_GRAY);
            ((XYPlot)_plot).setRangeGridlinePaint(Color.LIGHT_GRAY);
            ((XYPlot)_plot).setAxisOffset(new RectangleInsets(0,0,0,0));
        }

        //dhmay adding 2009/09/15.  As of jfree 1.0.13, several defaults changed annoyingly
        _chart.setBackgroundPaint(new Color(210, 210, 210));
        _plot.setBackgroundPaint(Color.white);
    }


    protected void init(Plot plot)
    {
        _plot = plot;
        if (_chart == null)
            init(new JFreeChart(null, null, _plot, showLegend));
        else
            init(_chart);
    }

    public JFreeChart getChart()
    {
        return _chart;

    }


    public String getToolTipText(MouseEvent e)
    {
        if (_chartPanel != null)
        {
            return _chartPanel.getToolTipText(e);
        }
        return null;
    }

    public Plot getPlot()
    {
        return _plot;
    }

    public ChartPanel getChartPanel()
    {
        return _chartPanel;
    }

    public void addItemToPopupMenu(JMenuItem item)
    {
        JPopupMenu popup = _chartPanel.getPopupMenu();
        popup.add(item);
    }

    public void addSeparatorToPopupMenu()
    {
        _chartPanel.getPopupMenu().add(new JSeparator());
    }

    /**
     * Add two new menu items to the popup menu, for saving to TSV and CSV files
     */
    protected void initPopupMenu()
    {
        addSeparatorToPopupMenu();

        //TSV
        JMenuItem saveTSVMenuItem = new JMenuItem(TextProvider.getText("SAVE_DATA_AS_TSV"));
        saveTSVMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JFileChooser fc = new JFileChooser();
                    int chooserStatus = fc.showOpenDialog(PanelWithChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = fc.getSelectedFile();
                    saveChartDataToTSV(outFile);
                }
            }
        );
        addItemToPopupMenu(saveTSVMenuItem);

        //CSV
        JMenuItem saveCSVMenuItem = new JMenuItem(TextProvider.getText("SAVE_DATA_AS_CSV"));
        saveCSVMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JFileChooser wfc = new JFileChooser();
                    int chooserStatus = wfc.showOpenDialog(PanelWithChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = wfc.getSelectedFile();
                    saveChartDataToCSV(outFile);
                }
            }
        );
        addItemToPopupMenu(saveCSVMenuItem);
    }

    protected void saveChartDataToCSV(File outFile)
    {
        saveChartDataToFile(outFile, ",");
    }

    protected void saveChartDataToTSV(File outFile)
    {
        saveChartDataToFile(outFile, "\t");
    }

    public void saveChartToImageFile(File outFile) throws IOException
    {
        ImageIO.write(createImage(),"png",outFile);
    }

    public void saveChartToImageFile(File outFile, int width, int height) throws IOException
    {
        ImageIO.write(createImage(width, height),"png",outFile);
    }

    public BufferedImage createImage(int width, int height)
    {
        return _chart.createBufferedImage(width,height);
    }

    public BufferedImage createImage()
    {
        int chartWidth = getWidth();
        int chartHeight = getHeight();

        if (chartWidth == 0)
            chartWidth = DEFAULT_WIDTH_FOR_IMAGE_FILE;
        if (chartHeight == 0)
            chartHeight = DEFAULT_HEIGHT_FOR_IMAGE_FILE;
        return createImage(chartWidth,chartHeight);
    }

    protected void saveChartDataToFile(File outFile, String delimiter)
    {
        _log.debug("saveChartDataToFile1, *delimiter*=*" + delimiter + "*");
        XYPlot xyPlot = (XYPlot) _plot;
        XYDataset dataset = xyPlot.getDataset();

        boolean hasZValues = false;
        if (dataset instanceof XYZDataset)
            hasZValues = true;

        PrintWriter pw = null;

        try
        {
            pw = new PrintWriter(outFile);
            int seriesCount = dataset.getSeriesCount();
            String headerLine = null;

            if (seriesCount > 1)
                headerLine = "series" + delimiter + "x" + delimiter + "y";
            else
                headerLine = "x" + delimiter + "y";
            if (hasZValues)
                headerLine = headerLine + delimiter + "z";

            pw.println(headerLine);

            for (int i = 0; i < seriesCount; i++)
            {
                int itemCount = dataset.getItemCount(i);
                for (int j = 0; j < itemCount; j++)
                {
                    String fileLine = null;
                    if (seriesCount > 1)
                        fileLine = i + delimiter + dataset.getX(i, j) + delimiter + dataset.getY(i, j);
                    else
                        fileLine = dataset.getX(i, j) + delimiter + dataset.getY(i, j);

                    if (hasZValues)
                        fileLine = fileLine + delimiter + ((XYZDataset) dataset).getZ(i,j);
                    pw.println(fileLine);
                }
                pw.flush();
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_SAVING_CHART_DATA"), e);
        }
        finally
        {
            if (pw != null)
                pw.close();
        }
    }

    public void setPreferredSize(Dimension newSize)
    {
        super.setPreferredSize(newSize);
        if (_chartPanel != null)
        {
            _chartPanel.setPreferredSize(newSize);
        }
    }


    public boolean isShowLegend()
    {
        return showLegend;
    }

    public void setShowLegend(boolean showLegend)
    {
        this.showLegend = showLegend;
    }

    public ChartDialog displayDialog(String title)
    {
        ChartDialog result = new ChartDialog(this);
        result.setTitle(title);
        result.setVisible(true);
        return result;
    }
    
    public MultiChartDisplayPanel displayInTab()
    {
        MultiChartDisplayPanel.addAndDisplayChartOnSingleton(this);
        return MultiChartDisplayPanel.getSingletonInstance();
    }

    public MultiChartDisplayPanel displayInTab(String newName)
    {
        setName(newName);
        return displayInTab();
    }
}
