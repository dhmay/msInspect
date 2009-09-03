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
package org.fhcrc.cpl.viewer.quant.gui;

import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.ChartMouseAndMotionListener;
import org.fhcrc.cpl.toolbox.Rounder;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;


public class ProteinRatioHistMouseListener extends ChartMouseAndMotionListener
{
    protected static Logger _log = Logger.getLogger(ProteinRatioHistMouseListener.class);    

    protected float selectedXMinValue;
    protected float selectedXMaxValue;

    protected JLabel maxLowRatioLabel;
    protected JLabel minHighRatioLabel;

    protected ProteinSummarySelectorFrame proteinSummaryFrame;


    public ProteinRatioHistMouseListener(PanelWithChart panelWithChart,
                                         ProteinSummarySelectorFrame proteinSummaryFrame,
                                         JLabel maxLowRatioLabel, JLabel minHighRatioLabel)
    {
        super(panelWithChart);
        this.maxLowRatioLabel = maxLowRatioLabel;
        this.minHighRatioLabel = minHighRatioLabel;
        this.proteinSummaryFrame = proteinSummaryFrame;
    }

    protected Rectangle2D selectedRegion;
    private Point selectedRegionStart;

    public void mouseMoved(MouseEvent e)
    {
    }

    public void mousePressed(MouseEvent e)
    {
        Rectangle2D screenDataArea = _chartPanel.getScreenDataArea(e.getX(), e.getY());
        if (screenDataArea != null)
        {
            this.selectedRegionStart = getPointInRectangle(e.getX(), e.getY(),
                    screenDataArea);
        }
        else
        {
            this.selectedRegionStart = null;
        }
    }

    public void mouseReleased(MouseEvent e)
    {
        try
        {
            transformAndSaveSelectedRegion();
            selectedXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
            selectedXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());

            drawAllButSelectedRegionHoriz(selectedRegion, new BasicStroke(2.0f), new Color(30,10,30,5));

            this.selectedRegion = null;
            this.selectedRegionStart = null;

            proteinSummaryFrame.minHighRatio = (float) Math.exp(selectedXMaxValue);
            proteinSummaryFrame.maxLowRatio = (float) Math.exp(selectedXMinValue);
            proteinSummaryFrame.updateExtremeRatioGUI();
        }
        catch (Exception ee) {}
    }

    protected void transformAndSaveSelectedRegion()
    {
        selectedXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
        selectedXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());
        maxLowRatioLabel.setText("Max Low Ratio: " + Rounder.round(Math.exp(selectedXMinValue), 2));
        minHighRatioLabel.setText("Min High Ratio: " + Rounder.round(Math.exp(selectedXMaxValue), 2));

//_log.debug("  Transformed: " + selectedXMinValue + ", " + selectedXMaxValue);        
    }


    public void mouseDragged(MouseEvent e)
    {
        if (this.selectedRegionStart == null || e.getX() < this.selectedRegionStart.getX())
        {
            return;
        }

        if(this.selectedRegion != null)
            drawAllButSelectedRegionHoriz(selectedRegion, new BasicStroke(2.0f), new Color(30,10,30,5));

        // Erase the previous zoom rectangle (if any)...
        Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea(
                (int) this.selectedRegionStart.getX(), (int) this.selectedRegionStart.getY());


        this.selectedRegion = new Rectangle2D.Double(
                this.selectedRegionStart.getX(), scaledDataArea.getMinY(),
                Math.abs(e.getX()-selectedRegionStart.getX()), scaledDataArea.getHeight());

        transformAndSaveSelectedRegion();


        // Draw the new zoom rectangle...
        drawAllButSelectedRegionHoriz(selectedRegion, new BasicStroke(2.0f), new Color(30,10,30,5));
    }
}
