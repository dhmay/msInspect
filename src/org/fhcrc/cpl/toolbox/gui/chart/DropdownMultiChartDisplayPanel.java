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
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;

import javax.swing.*;
import java.util.*;
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

    protected boolean displaySlideshowButton = false;

    protected JButton slideshowButton = null;

    protected int slideshowDelayMillis = 500;

    /**
     * Control which chart is displayed
     */
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

        slideshowButton = new JButton("Slideshow");
        ListenerHelper helper = new ListenerHelper(this);
        helper.addListener(slideshowButton,"slideshowButton_actionPerformed");
        slideshowButton.setVisible(displaySlideshowButton);
        

        add(slideshowButton);
        add(comboBox);
        add(panelHoldingChart);

        nameChartMap = new HashMap<String, PanelWithChart>();

        initialized = true;
    }

    public void slideshowButton_actionPerformed(ActionEvent event)
    {
        Thread slideshowThread = new Thread(new SlideshowManager(this));
        slideshowThread.start();
    }

    /**
     * Separate thread to manage the slideshow, because if we put the main thread to sleep, UI doesn't update
     */
    protected class SlideshowManager implements Runnable
    {
        protected DropdownMultiChartDisplayPanel dropdownPanel;

        public SlideshowManager(DropdownMultiChartDisplayPanel dropdownPanel)
        {
            this.dropdownPanel = dropdownPanel;
        }

        public void run()
        {
            dropdownPanel.resizeAllCharts();

            try
            {
                Thread.sleep(200);
            }
            catch (InterruptedException e)
            {

            }

            for (int i=0; i<dropdownPanel.comboBox.getItemCount(); i++)
            {
                dropdownPanel.updateUI();
                dropdownPanel.comboBox.setSelectedIndex(i);
                try
                {
                    Thread.sleep(slideshowDelayMillis);
                }
                catch (InterruptedException e)
                {

                }
            }
        }
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

    protected void resizeAllCharts()
    {
        if (getChartPanels() == null)
            return;
        for (PanelWithChart chart : getChartPanels())
        {
            resizeChart(chart);
        }
    }

    protected void resizeChartInFocus()
    {
        resizeChart(selectedChart);
    }

    protected void resizeChart(PanelWithChart chart)
    {
        long currentTime = new Date().getTime();
        long timeDiff = currentTime - lastResizeTime;

        if (timeDiff > 200)
        {
            if (chart != null)
            {
                int newChartWidth = getWidth() - 10;
                int newChartHeight = getHeight() - 20;
                if (comboBox != null)
                    newChartHeight -= comboBox.getHeight();
                chart.setPreferredSize(new Dimension(newChartWidth, newChartHeight));

                chart.updateUI();
            }
        }
        lastResizeTime = currentTime;
    }

    public boolean isDisplaySlideshowButton()
    {
        return displaySlideshowButton;
    }

    public void setDisplaySlideshowButton(boolean displaySlideshowButton)
    {
        this.displaySlideshowButton = displaySlideshowButton;
        if (slideshowButton != null)
            slideshowButton.setVisible(displaySlideshowButton);
        updateUI();
    }

    public int getSlideshowDelayMillis()
    {
        return slideshowDelayMillis;
    }

    public void setSlideshowDelayMillis(int slideshowDelayMillis)
    {
        this.slideshowDelayMillis = slideshowDelayMillis;
    }
}
