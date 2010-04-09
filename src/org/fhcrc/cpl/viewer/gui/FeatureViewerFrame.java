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
package org.fhcrc.cpl.viewer.gui;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.gui.chart.TabbedMultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.gui.chart.DropdownMultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ComponentListener;
import java.awt.event.ComponentEvent;
import java.util.*;

/**
 * GUI
 */
public class FeatureViewerFrame extends JDialog
{
    protected static Logger _log = Logger.getLogger(FeatureViewerFrame.class);

    //everything GUI
    protected FeatureViewer featureViewer;

    public FeatureViewerFrame(MSRun run)
    {
        super();

        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setSize(new Dimension(FeatureViewer.DEFAULT_WIDTH+10, FeatureViewer.DEFAULT_HEIGHT+10));
        setTitle("Feature Viewer");

        featureViewer = new FeatureViewer(run);
        add(featureViewer);

        
    }

    public void displayFeature(Feature feature)
    {
        featureViewer.displayFeature(feature);
    }

    public FeatureViewer getFeatureViewer() {
        return featureViewer;
    }

    public void setFeatureViewer(FeatureViewer featureViewer) {
        this.featureViewer = featureViewer;
    }
}
