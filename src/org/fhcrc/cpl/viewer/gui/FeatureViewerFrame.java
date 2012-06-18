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
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;

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

        if (ApplicationContext.getFrame() != null)
            setIconImage(ApplicationContext.getFrame().getIconImage());        

        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setSize(new Dimension(FeatureViewer.DEFAULT_WIDTH+10, FeatureViewer.DEFAULT_HEIGHT+10));
        setTitle("Feature Viewer");

        featureViewer = new FeatureViewer(run);
        add(featureViewer);
    }

    /**
     * Must call this before calling setVisible() or you get random NPEs from AWT
     * @param feature
     */
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
