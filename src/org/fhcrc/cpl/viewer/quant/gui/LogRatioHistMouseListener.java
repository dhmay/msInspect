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
import org.jfree.chart.panel.Overlay;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.event.OverlayChangeListener;

import java.util.List;
import java.util.ArrayList;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Rectangle2D;

//todo: fix behavior when user navigates to another window and back.  That's why the ChartChangeListener

public class LogRatioHistMouseListener extends ChartMouseAndMotionListener
{
    protected static Logger _log = Logger.getLogger(LogRatioHistMouseListener.class);

    //These are in chart-axis scale.  Initially set to sentinels
    protected float selectedRealXMinValue = Float.POSITIVE_INFINITY;
    protected float selectedRealXMaxValue = Float.POSITIVE_INFINITY;

    //listeners to be updated when the selected region changes
    protected List<ActionListener> rangeUpdateListeners;

    protected Color fillColor = new Color(31, 47, 31, 5);//30,10,30,5);   old, from the days of XOR
    protected Stroke stroke = new BasicStroke(2.0f);

    protected boolean shouldRedrawOldBeforeDrawingNew = true;

    //overlay for drawing the selection
    protected SelectedRegionOverlay selectionOverlay;

    public LogRatioHistMouseListener(PanelWithChart panelWithChart)
    {
        super(panelWithChart);
        rangeUpdateListeners = new ArrayList<ActionListener>();
        selectionOverlay = new SelectedRegionOverlay();
        _chartPanel.addOverlay(selectionOverlay);
    }

    public LogRatioHistMouseListener(PanelWithChart panelWithChart, ActionListener actionListener)
    {
        this(panelWithChart);
        rangeUpdateListeners.add(actionListener);
    }

    protected Rectangle2D selectedRegion;
    private Point selectedRegionStart;


    /**
     * When mouse moved, draw the ratio under the mouse pointer
     * @param e
     */
    public void mouseMoved(MouseEvent e)
    {
        double ratio = Rounder.round(Math.exp(transformMouseXValue(e.getX())), 2);
        drawRatioInBox(ratio);
    }

    /**
     * When mouse leaves chart, paint over the ratio
     * @param e
     */
    public void mouseExited(MouseEvent e)
    {
        drawBoxForRatio();
    }

    /**
     * Draw the box to contain the ratio
     */
    protected void drawBoxForRatio()
    {
        Graphics2D g = getChartPanelGraphics();
        g.setPaint(Color.LIGHT_GRAY);
        g.fillRect(15, 15, 40, 12);
    }

    /**
     * Draw the ratio in its box.  Separated from drawBoxForRatio so the box can be drawn empty
     * @param ratio
     */
    protected void drawRatioInBox(double ratio)
    {
        drawBoxForRatio();
        Graphics2D g = getChartPanelGraphics();

        g.setFont(new Font("Verdana", Font.PLAIN, 10));
        g.setColor(Color.BLACK);
        g.setPaint(Color.BLACK);
        g.drawString("" + ratio, 16, 24);
    }

    /**
     * Save the selected area
     * @param e
     */
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
            if(this.selectedRegion != null)
                drawOrUndrawRegion();

            transformAndSaveSelectedRegion();
            selectedRealXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
            selectedRealXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());

            drawOrUndrawRegion();

            this.selectedRegionStart = null;

            for (ActionListener listener : rangeUpdateListeners)
            {
                listener.actionPerformed(null);
            }
        }
        catch (Exception ee) {}
    }

    /**
     * Transform the region into values in units of the chart X axis
     */
    protected void transformAndSaveSelectedRegion()
    {
        selectedRealXMinValue = (float) super.transformMouseXValue(selectedRegion.getX());
        selectedRealXMaxValue = (float) super.transformMouseXValue(selectedRegion.getX()+selectedRegion.getWidth());

        for (ActionListener listener : rangeUpdateListeners)
        {
            listener.actionPerformed(null);
        }
    }

    /**
     * Set the selected region in chartpanel coordinates, based on axis-scale coordinates given 
     * @param minValue
     * @param maxValue
     */
    public void setSelectedRegionWithChartValues(float minValue, float maxValue)
    {
        selectedRealXMinValue = minValue;
        selectedRealXMaxValue = maxValue;

//System.err.println("***DRAWING, " + minValue + ", " + maxValue + ", region: " + selectedRegion);
//        drawOrUndrawRegion();
    }

    /**
     * Undraw the previous selected region (if it was drawn), calculate the new regions, draw again, save
     * the points, and draw the numeric ratio in its little box
     * @param e
     */
    public void mouseDragged(MouseEvent e)
    {

        if (this.selectedRegionStart == null || e.getX() < this.selectedRegionStart.getX())
        {
            return;
        }

        if(this.selectedRegion != null)
            drawOrUndrawRegion();

        // Erase the previous zoom rectangle (if any)...
        Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea();

        this.selectedRegion = new Rectangle2D.Double(
                this.selectedRegionStart.getX(), scaledDataArea.getMinY(),
                Math.min(Math.abs(e.getX() - selectedRegionStart.getX()),
                        _chartPanel.getWidth() - this.selectedRegionStart.getX()),
                scaledDataArea.getHeight());
        transformAndSaveSelectedRegion();

        // Draw the new zoom rectangle...
        drawOrUndrawRegion();

        double ratio = Rounder.round(Math.exp(transformMouseXValue(e.getX())), 2);
        drawRatioInBox(ratio);
    }

    /**
     * Draw or undraw the selected region.  Since the region is XORed with the chart image, have to undraw the
     * old one before drawing the new one.
     */
    protected void drawOrUndrawRegion()
    {
        selectionOverlay.paintOverlay((Graphics2D)_chartPanel.getGraphics(), _chartPanel);
    }

    public float getSelectedRealXMinValue()
    {
        return selectedRealXMinValue;
    }

    public void setSelectedRealXMinValue(float selectedRealXMinValue)
    {
        this.selectedRealXMinValue = selectedRealXMinValue;
    }

    public float getSelectedRealXMaxValue()
    {
        return selectedRealXMaxValue;
    }

    public void setSelectedRealXMaxValue(float selectedRealXMaxValue)
    {
        this.selectedRealXMaxValue = selectedRealXMaxValue;
    }

    /**
     * A ChartPanel overlay that draws in the selected region
     */
    protected class SelectedRegionOverlay implements Overlay
    {
        public void paintOverlay(Graphics2D graphics2D, ChartPanel chartPanel)
        {
            //infinity is a sentinel value to tell us the real x and y values haven't been set yet.
            //If they have been set but selectedRegion is null, transform them into selectedRegion
            if (selectedRegion == null &&
                !Float.isInfinite(selectedRealXMinValue) && !Float.isInfinite(selectedRealXMaxValue))
            {
                Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea();
                double transformedMin = transformXValueToMouse(selectedRealXMinValue);
                double transformedMax = transformXValueToMouse(selectedRealXMaxValue);

                selectedRegion = new Rectangle2D.Double(
                        transformedMin, scaledDataArea.getMinY(),
                        transformedMax - transformedMin, scaledDataArea.getHeight());
            }
            if (selectedRegion != null)
            {
                drawAllButSelectedRegionHoriz(selectedRegion, stroke, fillColor, true, graphics2D);
            }
        }

        public void addChangeListener(OverlayChangeListener overlayChangeListener)
        {
        }

        public void removeChangeListener(OverlayChangeListener overlayChangeListener)
        {
        }
    }    
}
