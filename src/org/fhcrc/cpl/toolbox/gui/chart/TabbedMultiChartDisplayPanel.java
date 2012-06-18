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

import javax.swing.*;
import java.util.Date;
import java.awt.*;

/**
 * display multiple charts in a panel.  Supports saving of all charts at once, to a directory
 *
 * Only a tabbed display is supported at the moment
 */
public class TabbedMultiChartDisplayPanel extends MultiChartDisplayPanel
{
    protected static Logger _log =
            Logger.getLogger(TabbedMultiChartDisplayPanel.class);

    protected JTabbedPane tabbedPane = null;

    public static final int TABBEDPANE_VERTICAL_SLOP = 5;

    protected GridBagConstraints gbc = null;

    //dhmay adding 20100120 to resize things differently on mac
    protected boolean isMacOS = false;

    public TabbedMultiChartDisplayPanel()
    {
        super();
        init();
    }

    protected void init()
    {
        setLayout(new FlowLayout());
//        gbc = new GridBagConstraints();
//        gbc.fill = GridBagConstraints.BOTH;
//        gbc.gridwidth = GridBagConstraints.REMAINDER;
//        gbc.insets = new Insets(0, 0, 0, 0);
        if (tabbedPane == null)
        {
            tabbedPane = new JTabbedPane();
            tabbedPane.addChangeListener(new TabFocusListener());
//            tabbedPane.setLayout(new GridBagLayout());
            add(tabbedPane);//, gbc);
            tabbedPane.updateUI();
        }
        String vers = System.getProperty("os.name").toLowerCase();
        if (vers.contains("mac"))
        {
            isMacOS = true;
        }

    }


    protected void addChartPanelToGUI(PanelWithChart newPanel)
    {
        tabbedPane.add(newPanel);
        if (chartPanels.size() == 1)
            resizeChartInFocus();
    }

//    public void setSize(int width, int height)
//    {
//        adjustTabbedPaneSize(width, height);
//        super.setSize(width, height);
//    }
//
    public void setPreferredSize(Dimension dim)
    {
        super.setPreferredSize(dim);
        tabbedPane.setPreferredSize(dim);
        tabbedPane.updateUI();
        resizeChartInFocus();

    }



    protected void removeChartPanelFromGUI(PanelWithChart panel)
    {
        tabbedPane.remove(panel);
        if (chartPanels.isEmpty())
            tabbedPane.setPreferredSize(new Dimension(0,0));
    }            

    public void resizeChartInFocus()
    {
        long currentTime = new Date().getTime();
        long timeDiff = currentTime - lastResizeTime;
//System.err.println("Tabbed resizeC 1, diff=" + timeDiff + " vs. " + resizeDelayMS + ", tP null? " + (tabbedPane == null));
        if (timeDiff > resizeDelayMS)
        {
            if (tabbedPane != null)
            {
                PanelWithChart selectedChart = (PanelWithChart) tabbedPane.getSelectedComponent();
                //System.err.println("\tsC null? " + (selectedChart == null));

                if (selectedChart != null)
                {                         
                    Point chartLocation = selectedChart.getLocation();
                    int widthSlop = 5;
                    int heightSlop=10;
                    if (isMacOS)
                    {
                        widthSlop=15;
                        heightSlop=30;
                    }
                    int newChartWidth = (int) getVisibleRect().getWidth() - 2 * (widthSlop + chartLocation.x - getLocation().x);
                    int newChartHeight = (int) getVisibleRect().getHeight() - heightSlop - (chartLocation.y - getLocation().y);
//                    tabbedPane.setPreferredSize(new Dimension((int) getVisibleRect().getWidth(), (int)getVisibleRect().getHeight()));
//System.err.println("\t\tsizing to " + newChartWidth + ", " + newChartHeight);
                    selectedChart.setPreferredSize(new Dimension(newChartWidth, newChartHeight));
                    selectedChart.updateUI();
//throw new RuntimeException();

                }
            }
            lastResizeTime = currentTime;
            
        }


//        if (tabbedPane == null)
//            init();
////System.err.println("Explicitly adjusting tabbedpane to " + width + ", " + height);
//        tabbedPane.setPreferredSize(new Dimension(getWidth(), getHeight()));
    }

    public int getSelectedIndex()
    {
        if (tabbedPane == null)
            return 0;
        return tabbedPane.getSelectedIndex();
    }

    public void setSelectedIndex(int index)
    {
        if (tabbedPane == null)
            return;
        try
        {
            tabbedPane.setSelectedIndex(index);
        }
        catch (Exception e) {}
    }

}
