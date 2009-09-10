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
import org.apache.log4j.Logger;
import org.jfree.chart.event.ChartChangeListener;
import org.jfree.chart.event.ChartChangeEvent;

import java.util.List;
import java.util.ArrayList;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Rectangle2D;

//todo: fix behavior when user navigates to another window and back.  That's why the ChartChangeListener

public class LogRatioHistMouseListener extends ChartMouseAndMotionListener implements ChartChangeListener
{
    protected static Logger _log = Logger.getLogger(LogRatioHistMouseListener.class);

    protected float selectedXMinValue;
    protected float selectedXMaxValue;

    protected List<ActionListener> rangeUpdateListeners;

    protected Color fillColor = new Color(30,10,30,5);
    protected Stroke stroke = new BasicStroke(2.0f);

    protected boolean regionIsDrawn = false;
    protected boolean shouldRedrawOldBeforeDrawingNew = true;

    public LogRatioHistMouseListener(PanelWithChart panelWithChart)
    {
        super(panelWithChart);
        rangeUpdateListeners = new ArrayList<ActionListener>();
    }

    public LogRatioHistMouseListener(PanelWithChart panelWithChart, ActionListener actionListener)
    {
        this(panelWithChart);
        rangeUpdateListeners.add(actionListener);
    }

    protected Rectangle2D selectedRegion;
    private Point selectedRegionStart;



    public void chartChanged(ChartChangeEvent e)
    {
//System.err.println("CHANGE!!!! " + e.getType());

        regionIsDrawn = false;
    }

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

    public void addRangeUpdateListener(ActionListener actionListener)
    {
        rangeUpdateListeners.add(actionListener);
    }

    public void mouseReleased(MouseEvent e)
    {
        try
        {
            if(this.selectedRegion != null && regionIsDrawn)
            {

                drawOrUndrawRegion();
            }

            transformAndSaveSelectedRegion();
            selectedXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
            selectedXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());

            drawOrUndrawRegion();

            //this.selectedRegion = null;
            this.selectedRegionStart = null;

            for (ActionListener listener : rangeUpdateListeners)
            {
                listener.actionPerformed(null);
            }
        }
        catch (Exception ee) {}
    }

    protected void transformAndSaveSelectedRegion()
    {
        selectedXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
        selectedXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());


        for (ActionListener listener : rangeUpdateListeners)
        {
            listener.actionPerformed(null);
        }
//_log.debug("  Transformed: " + selectedXMinValue + ", " + selectedXMaxValue);        
    }

    public void setSelectedRegionWithChartValues(float minValue, float maxValue)
    {
        Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea();
//System.err.println("Scaled: " + scaledDataArea);
        this.selectedRegion = new Rectangle2D.Double(
                transformXValueToMouse(minValue), scaledDataArea.getMinY(),
                transformXValueToMouse(maxValue), scaledDataArea.getHeight());
//System.err.println(selectedRegion);
    }



    public void mouseDragged(MouseEvent e)
    {

        if (this.selectedRegionStart == null || e.getX() < this.selectedRegionStart.getX())
        {
            return;
        }
//Graphics2D g2 = getChartPanelGraphics();
//Rectangle2D rect = _chartPanel.getScreenDataArea();
//g2.clearRect(0,0, (int) rect.getWidth(), (int) rect.getHeight());
//System.err.println("DRAGGING, drawn? " + regionIsDrawn);
        if(this.selectedRegion != null && regionIsDrawn)
            drawOrUndrawRegion();


        // Erase the previous zoom rectangle (if any)...
        Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea(
                (int) this.selectedRegionStart.getX(), (int) this.selectedRegionStart.getY());


        this.selectedRegion = new Rectangle2D.Double(
                this.selectedRegionStart.getX(), scaledDataArea.getMinY(),
                Math.abs(e.getX()-selectedRegionStart.getX()), scaledDataArea.getHeight());

        transformAndSaveSelectedRegion();


        // Draw the new zoom rectangle...
        drawOrUndrawRegion();
    }

    /**
     * Since the region is XORed, we have to draw it again to undraw it
     */
    protected void drawOrUndrawRegion()
    {
//System.err.println("Drawing! " + selectedRegion);        
        drawAllButSelectedRegionHoriz(selectedRegion, stroke, fillColor, true);
        regionIsDrawn = !regionIsDrawn;
    }

    public float getSelectedXMinValue()
    {
        return selectedXMinValue;
    }

    public void setSelectedXMinValue(float selectedXMinValue)
    {
        this.selectedXMinValue = selectedXMinValue;
    }

    public float getSelectedXMaxValue()
    {
        return selectedXMaxValue;
    }

    public void setSelectedXMaxValue(float selectedXMaxValue)
    {
        this.selectedXMaxValue = selectedXMaxValue;
    }
}
