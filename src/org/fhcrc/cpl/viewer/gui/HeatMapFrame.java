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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.gui.HeatMapPanel;
import org.fhcrc.cpl.viewer.util.SharedProperties;

import java.awt.Container;
import java.awt.event.ActionEvent;
import java.util.List;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import org.swixml.SwingEngine;

/**
 * HeatMapFrame - Main frame that holds the heat map images
 */
public class HeatMapFrame extends JFrame
{
    private static HeatMapFrame instance = null;
    private HeatMapPanel heatMapPanel;
    
    private HeatMapFrame()
    {
        super(WorkbenchFrame.getAppName() + " - Feature Heat Map");
        try
        {
            JMenuBar menuBar = (JMenuBar)new SwingEngine(this).render("org/fhcrc/cpl/viewer/gui/HeatMapFrameMenu.xml");
            setJMenuBar(menuBar);
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Error loading Heat Map menus", x);
            throw new RuntimeException(x);
        }
        Container contentPane = getContentPane();
        heatMapPanel = new HeatMapPanel(this);
        contentPane.add(heatMapPanel);
        setResizable(true);
        setIconImage(ApplicationContext.getFrame().getIconImage());
        setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        pack();
    }

    /**
     *
     */
    public static synchronized HeatMapFrame getInstance()
    {
        if ( null == instance )
            instance = new HeatMapFrame();
        return instance;
    }

    /**
     * Return the *last* displayed feature set (the one drawn on top).
     * Returns null if none displayed
     */
    public static FeatureSet getDisplayedFeatureSet()
    {
        List featureSets = (List)ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);

        if ( featureSets == null )
            return null;

        for (int i = featureSets.size() - 1; i >= 0; i--)
        {
            FeatureSet fs = (FeatureSet)featureSets.get(i);
            if ( fs.isDisplayed() )
                return fs;
        }

        return null;
    }

    /**
     *
     */
    public void display()
    {
        // Make the HeatMapFrame visible iff the feature windows are ready
        heatMapPanel.updateFeatureSet(true);
    }

    /*----------------------------------------------------------------------*
     * Menu actions
     *----------------------------------------------------------------------*/
    public Action closeWindowAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            setVisible(false);
        }
    };

    public Action sortMassAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            heatMapPanel.setSortOrderMass();
        }
        };
    
    public Action sortKlAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            heatMapPanel.setSortOrderKl();
        }
    };
    
    public Action sortIntensityAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            heatMapPanel.setSortOrderIntensity();
        }
    };

    public Action sortTotalIntensityAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            heatMapPanel.setSortOrderTotalIntensity();
        }
    };
    
    public Action sortScanAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent event)
        {
            heatMapPanel.setSortOrderScan();
        }
    };
    
}

