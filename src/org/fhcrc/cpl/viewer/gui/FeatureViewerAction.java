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

import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

/**
 * User: mbellew
 * Date: May 26, 2004
 * Time: 2:58:14 PM
 *
 * dhmay changing 2009/01/02: adding default filtering of minimum 2 peaks and maximum KL 3.0
 *
 */
public class FeatureViewerAction extends AbstractAction implements PropertyChangeListener
{

    protected FeatureViewerFrame featureViewerFrame;

    public FeatureViewerAction()
    {
        super("Open Feature Viewer");
        ApplicationContext.addPropertyChangeListener(SharedProperties.MS_RUN, this);
    }


    public void actionPerformed(ActionEvent event)
    {
        MSRun run;
        run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        if (null == run)
        {
            ApplicationContext.setMessage("No run loaded.");
            return;
        }

        featureViewerFrame = new FeatureViewerFrame(run);
        featureViewerFrame.setVisible(true); 
        ApplicationContext.setProperty(SharedProperties.FEATURE_VIEWER, featureViewerFrame);

        showCurrentSelectedIfFeature();

        ApplicationContext.addPropertyChangeListener(SharedProperties.SELECTED, new PropertyChangeListener()
        {
            public void propertyChange(PropertyChangeEvent evt)
            {
                showCurrentSelectedIfFeature();
            }
        });

    }

    protected void showCurrentSelectedIfFeature()
    {
        Object currentSelected = ApplicationContext.getProperty(SharedProperties.SELECTED);
        if (currentSelected != null && currentSelected instanceof Feature)
            featureViewerFrame.getFeatureViewer().displayFeature((Feature) currentSelected);
    }

    public void propertyChange(PropertyChangeEvent event)
    {
        setEnabled(null != ApplicationContext.getProperty(SharedProperties.MS_RUN));
    }

    public FeatureViewerFrame getFeatureViewerFrame() {
        return featureViewerFrame;
    }
}
