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


    protected void addChartPanelToGUI(PanelWithChart newPanel)
    {
        if (tabbedPane == null)
        {
            tabbedPane = new JTabbedPane();
            tabbedPane.addChangeListener(new TabFocusListener());
            add(tabbedPane);
        }
        tabbedPane.add(newPanel);        
    }

    protected void removeChartPanelFromGUI(PanelWithChart panel)
    {
        tabbedPane.remove(panel);
    }

    protected void resizeChartInFocus()
    {
        long currentTime = new Date().getTime();
        long timeDiff = currentTime - lastResizeTime;

        if (timeDiff > 200)
        {
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
        lastResizeTime = currentTime;
    }

}
