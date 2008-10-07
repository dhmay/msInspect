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
import org.jfree.chart.ChartPanel;

import javax.swing.*;
import java.util.Date;
import java.util.Map;
import java.util.HashMap;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

/**
 * This implementation uses a dropdown box to control display
 */
public class DropdownMultiChartDisplayPanel extends MultiChartDisplayPanel
{
    protected static Logger _log =
            Logger.getLogger(DropdownMultiChartDisplayPanel.class);

    protected JComboBox comboBox = null;
    protected JPanel panelHoldingChart = null;

    protected Map<String, PanelWithChart> nameChartMap = null;
    protected PanelWithChart selectedChart = null;

    protected boolean initialized = false;

    protected class ComboListener implements ActionListener
    {
        public void actionPerformed(ActionEvent event)
        {
            if (selectedChart != null)
            {
                panelHoldingChart.remove(selectedChart);
            }
            String name = (String) comboBox.getSelectedItem();
            selectedChart = nameChartMap.get(name);
            panelHoldingChart.add(selectedChart);
            resizeChartInFocus();
            panelHoldingChart.updateUI();
        }
    }

    protected void initialize()
    {
        comboBox = new JComboBox();
        comboBox.addActionListener(new ComboListener());
        panelHoldingChart = new JPanel();

        add(comboBox);
        add(panelHoldingChart);

        nameChartMap = new HashMap<String, PanelWithChart>();

        initialized = true;
    }

    protected void removeChartPanelFromGUI(PanelWithChart newPanel)
    {
        nameChartMap.remove(newPanel.getName());           
        comboBox.removeItem(newPanel.getName());
    }


    protected void addChartPanelToGUI(PanelWithChart newPanel)
    {
        if (!initialized)
            initialize();

        //prevent name collisions by renaming the chart if necessary
        int appendedNum = 1;
        while (nameChartMap.containsKey(newPanel.getName()))
            appendedNum++;
        if (appendedNum > 1)
            newPanel.setName(newPanel.getName() + appendedNum);

        nameChartMap.put(newPanel.getName(), newPanel);
        comboBox.addItem(newPanel.getName());
    }

    protected void resizeChartInFocus()
    {
        long currentTime = new Date().getTime();
        long timeDiff = currentTime - lastResizeTime;

        if (timeDiff > 200)
        {
            if (selectedChart != null)
            {
                if (selectedChart != null)
                {
                    int newChartWidth = getWidth() - 10;
                    int newChartHeight = getHeight() - comboBox.getHeight() - 20;
                    selectedChart.setPreferredSize(new Dimension(newChartWidth, newChartHeight));

                    selectedChart.updateUI();
                }
            }
        }
        lastResizeTime = currentTime;
    }
}
