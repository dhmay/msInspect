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

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.apache.log4j.Logger;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;


/**
 * Superclass for MouseListener/MouseMotionListener combos to hang on charts, with some utility methods
 */
public abstract class ChartMouseAndMotionListener implements MouseListener, MouseMotionListener
{
    protected static Logger _log = Logger.getLogger(ChartMouseAndMotionListener.class);    

    protected ChartPanel _chartPanel;

    protected NumberAxis domainAxis;
    protected NumberAxis rangeAxis;

    public ChartMouseAndMotionListener()
    {

    }


    public ChartMouseAndMotionListener(PanelWithChart panelWithChart)
    {
        _chartPanel = panelWithChart.getChartPanel();
        XYPlot p = _chartPanel.getChart().getXYPlot();
        domainAxis = (NumberAxis) p.getDomainAxis();
        rangeAxis = (NumberAxis) p.getRangeAxis();
    }

    //Copied verbatim from ChartPanel
    /**
     * Returns a point based on (x, y) but constrained to be within the bounds
     * of the given rectangle.  This method could be moved to JCommon.
     *
     * @param x  the x-coordinate.
     * @param y  the y-coordinate.
     * @param area  the rectangle (<code>null</code> not permitted).
     *
     * @return A point within the rectangle.
     */
    protected Point getPointInRectangle(int x, int y, Rectangle2D area)
    {
        x = (int) Math.max(Math.ceil(area.getMinX()), Math.min(x,
                Math.floor(area.getMaxX())));
        y = (int) Math.max(Math.ceil(area.getMinY()), Math.min(y,
                Math.floor(area.getMaxY())));
        return new Point(x, y);
    }




    public void mouseClicked(MouseEvent e)
    {
    }

    public void mousePressed(MouseEvent e)
    {
    }

    public void mouseReleased(MouseEvent e)
    {
    }


    public void mouseEntered(MouseEvent e)
    {
    }

    public void mouseExited(MouseEvent e)
    {
    }

    public void mouseDragged(MouseEvent e)
    {
    }

    public void mouseMoved(MouseEvent e)
    {
    }

    /**
     * Utility method for getting graphics
     * @return
     */
    protected Graphics2D getChartPanelGraphics()
    {
        return (Graphics2D) _chartPanel.getGraphics();
    }



    /**
     * Transform a mouse X value into a value in the units of the X axis of the chart.
     * Note: if there were multiple subplots, this would need to take a MouseEvent to determine which one
     * @param rawValue
     * @return
     */
    protected double transformMouseXValue(double rawValue)
    {
        Rectangle2D screenDataArea = _chartPanel.getScreenDataArea();

        double leftmostOnAxis = domainAxis.getLowerBound();
        double rightmostOnAxis = domainAxis.getUpperBound();
        double leftmostonscreen = screenDataArea.getX();
        double rightmostonscreen = leftmostonscreen+screenDataArea.getWidth();
        double slope = (rightmostOnAxis-leftmostOnAxis)/(rightmostonscreen-leftmostonscreen);
        return ((slope*(rawValue-leftmostonscreen))+leftmostOnAxis);
    }

    /**
     * Transform a value in the units of the X axis of the chart into a mouse x value
     * Note: if there were multiple subplots, this would need to take a MouseEvent to determine which one
     * @param xValue
     * @return
     */
    protected double transformXValueToMouse(double xValue)
    {
        Rectangle2D screenDataArea = _chartPanel.getScreenDataArea();

        double leftmostOnAxis = domainAxis.getLowerBound();
        double rightmostOnAxis = domainAxis.getUpperBound();
        double leftmostonscreen = screenDataArea.getX();
        double rightmostonscreen = leftmostonscreen+screenDataArea.getWidth();
        double slope = (rightmostOnAxis-leftmostOnAxis)/(rightmostonscreen-leftmostonscreen);
        return ((xValue - leftmostOnAxis) / slope) + leftmostonscreen;
    }


    //Lifted from ChartPanel because in JFree it is private
    /**
     * Draws zoom rectangle (if present).
     * The drawing is performed in XOR mode, therefore
     * when this method is called twice in a row,
     * the second call will completely restore the state
     * of the canvas.
     *
     * @param selectedRegion
     * @param stroke
     * @param fillColor
     */
    protected void drawSelectedRegion(Rectangle2D selectedRegion, Stroke stroke,
                                      Color fillColor, boolean xorMode)
    {
        Graphics2D g2 = getChartPanelGraphics();
        // Set XOR mode to draw the zoom rectangle
        if(g2 == null) return;
        Paint origColor = g2.getPaint();

        if (selectedRegion != null)
        {
            if (xorMode)
            {
                g2.setXORMode(fillColor);
            }
            else
                g2.setPaint(fillColor);
            g2.fill(selectedRegion);
            g2.setPaint(Color.black);
            g2.setStroke(stroke);
            if (xorMode)
                g2.setPaint(Color.white);
            else
                g2.setPaint(Color.black);

            g2.draw(selectedRegion);
            g2.setPaintMode();
            g2.setPaint(origColor);
        }
    }

    //Lifted from ChartPanel because in JFree it is private
    /**
     * Draw a box to the left of the selected region, and one to the right, that together
     * encompass everything but the selection
     *
     *
     * @param selectedRegion
     * @param stroke
     * @param color
     */
    protected void drawAllButSelectedRegionHoriz(Rectangle2D selectedRegion, Stroke stroke, Color color,
                                                 boolean xorMode)
    {
        Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea();
        Rectangle2D firstBox = new Rectangle((int) scaledDataArea.getMinX(),(int) scaledDataArea.getMinY(),
                (int) (selectedRegion.getMinX() - scaledDataArea.getMinX()),
                (int)scaledDataArea.getMaxY());
        Rectangle2D secondBox = new Rectangle(
                (int) selectedRegion.getMaxX(),
                (int) scaledDataArea.getMinY(),
                (int) (scaledDataArea.getMaxX() - scaledDataArea.getMinX()),
                (int)scaledDataArea.getMaxY());
        drawSelectedRegion(firstBox, stroke, color, xorMode);
        drawSelectedRegion(secondBox, stroke, color, xorMode);
    }


    public ChartPanel getChartPanel() {
        return _chartPanel;
    }

    public void setChartPanel(ChartPanel _cp) {
        this._chartPanel = _cp;
    }
}
