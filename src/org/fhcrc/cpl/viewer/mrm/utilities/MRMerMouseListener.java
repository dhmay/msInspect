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
package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.fhcrc.cpl.viewer.mrm.MRMTransition;
import org.fhcrc.cpl.viewer.mrm.MRMDaughter;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.InputEvent;
import java.awt.geom.Rectangle2D;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Oct 15, 2007
 * Time: 2:29:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class MRMerMouseListener implements MouseListener, MouseMotionListener {

    private boolean shifted = false;

    public ChartPanel getChartPanel() {
        return _cp;
    }

    public void setChartPanel(ChartPanel _cp) {
        this._cp = _cp;
    }

    protected ChartPanel _cp;

    public MRMerMouseListener(ChartPanel cp, PeaksTableModel ptm) {
        _cp = cp;
        _ptm = ptm;
    }

    public PeaksTableModel getPtm() {
        return _ptm;
    }

    public void setPtm(PeaksTableModel _ptm) {
        this._ptm = _ptm;
    }

    PeaksTableModel _ptm;


    public void mouseClicked(MouseEvent e) {
        if(e.getSource() instanceof ChartPanel &&
            (
             (e.getButton() == MouseEvent.BUTTON2)
             ||
             (e.getButton() == MouseEvent.BUTTON1) && ((e.getModifiers() & InputEvent.CTRL_MASK) != 0)
            )          
           ){
           CenterZoomNumberAxis czna = (CenterZoomNumberAxis) _cp.getChart().getXYPlot().getDomainAxis();
           NumberAxis range = (NumberAxis)_cp.getChart().getXYPlot().getRangeAxis();
           Rectangle2D screenDataArea = _cp.getScreenDataArea(e.getX(), e.getY());
           double y1 = czna.getLowerBound();
           double y2 = czna.getUpperBound();
           double x1 = screenDataArea.getX();
           double x2 = x1 + screenDataArea.getWidth();
           double transformedx = (((y2-y1)/(x2-x1))*(e.getX() - x1))+y1;
           MRMDialog mrmd = (MRMDialog) MRMAncestor();
           PeaksTableModel model = (PeaksTableModel) mrmd.peaksTable.getModel();
           MRMTransition mrt = mrmd.transitionOnPlot;
           mrt.setCalcXatMaxYAllDaughters(transformedx);
           mrt.setCalcMaxYAllDaughters(range.getLowerBound()+ 0.95*(range.getUpperBound()-range.getLowerBound()));
           model.setValueAt(new Float(mrt.getCalcXatMaxYAllDaughters()),mrt.getTableRow(), MRMDialog.peaksData.MidTime.colno);
           for(MRMDaughter d: mrt.getDaughters().values()){
              model.setValueAt(new Float(mrt.getCalcXatMaxYAllDaughters()),d.getElutionDataTableRow(), MRMDialog.peaksData.MidTime.colno);
           }
           mrmd.updateChartsAndFields(false);
        }
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {
            _cp.mouseClicked(e);
        } else {
            _cp.mouseClicked(e);
        }
    }

    protected Rectangle2D coElutionRegion;
    private Point coElutionStart;

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
    private Point getPointInRectangle(int x, int y, Rectangle2D area) {
        x = (int) Math.max(Math.ceil(area.getMinX()), Math.min(x,
                Math.floor(area.getMaxX())));
        y = (int) Math.max(Math.ceil(area.getMinY()), Math.min(y,
                Math.floor(area.getMaxY())));
        return new Point(x, y);
    }

    public void mousePressed(MouseEvent e) {
       if((e.isShiftDown() || e.getButton() == MouseEvent.BUTTON3) || shifted ) {
            if (this.coElutionRegion == null) {
                shifted = true;
                Rectangle2D screenDataArea = _cp.getScreenDataArea(e.getX(), e.getY());
                if (screenDataArea != null) {
                    this.coElutionStart = getPointInRectangle(e.getX(), e.getY(),
                            screenDataArea);
                }
                else {
                    this.coElutionStart = null;
                }
            }
        } else {
           _cp.mousePressed(e);
       }
    }

    Component MRMAncestor() {
        Component c = null;
        for(c = _cp; c != null && c.getClass() != MRMDialog.class; c=c.getParent());
        return (c == null) ? null : (MRMDialog) c;
    }


    public void mouseReleased(MouseEvent e) {
     try {
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {
//            Rectangle2D scaledDataArea = _chartPanel.getScreenDataArea(
//                    (int) this.coElutionStart.getX(), (int) this.coElutionStart.getY());
            JFreeChart jfc = _cp.getChart();
            XYPlot p = jfc.getXYPlot();
            CenterZoomNumberAxis czna = (CenterZoomNumberAxis)p.getDomainAxis();
            Rectangle2D screenDataArea = _cp.getScreenDataArea(e.getX(), e.getY());
            Rectangle2D plotboundaries = _cp.getChartRenderingInfo().getPlotInfo().getPlotArea();

            double leftmostOnAxis = czna.getLowerBound();
            double rightmostOnAxis = czna.getUpperBound();
            double leftmostOnRange = this.coElutionRegion.getX();
            double rightmostOnRange = this.coElutionRegion.getX()+this.coElutionRegion.getWidth();
            double leftmostonscreen = screenDataArea.getX();
            double rightmostonscreen = leftmostonscreen+screenDataArea.getWidth();
            double slope = (rightmostOnAxis-leftmostOnAxis)/(rightmostonscreen-leftmostonscreen);
            double transformedLeftRange = (slope*(leftmostOnRange-leftmostonscreen))+leftmostOnAxis;
            double transformedRightRange = (slope*(rightmostOnRange-leftmostonscreen))+leftmostOnAxis;
            shifted = false;
            MRMDialog ultimateParent = (MRMDialog)MRMAncestor();
            if(ultimateParent != null) {
                MRMTransition transition = ultimateParent.transitionOnPlot;
                MRMTransition mrmt = transition;
                if(mrmt != null) {
                    int row = mrmt.getTableRow();
                    _ptm.data[row][MRMDialog.peaksData.CoStart.colno] = new Float(0f);
                    _ptm.data[row][MRMDialog.peaksData.CoEnd.colno] = new Float(10000000f);
                    _ptm.setValueAt(new Float(transformedRightRange),row, MRMDialog.peaksData.CoEnd.colno);
                    _ptm.setValueAt(new Float(transformedLeftRange),row,MRMDialog.peaksData.CoStart.colno);
                }
            }
            Graphics2D g2 = (Graphics2D) _cp.getGraphics();
            if(this.coElutionRegion != null) drawCoElutionRegion(g2);
            this.coElutionRegion = null;
            this.coElutionStart = null;
        } else {
            _cp.mouseReleased(e);
        }
     } catch (Exception ee) {}  
    }

    public void mouseEntered(MouseEvent e) {
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {

        } else {
            _cp.mouseEntered(e);
        }
    }

    public void mouseExited(MouseEvent e) {
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {
         
        } else {
            _cp.mouseExited(e);
        }
    }

    //Lifted from ChartPanel because in JFree it is private
    /**
     * Draws zoom rectangle (if present).
     * The drawing is performed in XOR mode, therefore
     * when this method is called twice in a row,
     * the second call will completely restore the state
     * of the canvas.
     *
     * @param g2 the graphics device.
     */
    private void drawCoElutionRegion(Graphics2D g2) {
        // Set XOR mode to draw the zoom rectangle
        if(g2 == null) return;
        Paint origColor = g2.getPaint();
//        g2.setXORMode(Color.gray);
        g2.setXORMode(new Color(30,10,30,5));
        if (this.coElutionRegion != null) {
                g2.fill(this.coElutionRegion);
                g2.setPaint(Color.white);
                g2.setStroke(new BasicStroke(3.0f));
                g2.draw(this.coElutionRegion);
        }
        g2.setPaintMode();
        g2.setPaint(origColor);
    }

    public void mouseDragged(MouseEvent e) {
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {
            if (this.coElutionStart == null) {
                return;
            }
            Graphics2D g2 = (Graphics2D) _cp.getGraphics();
            if(this.coElutionRegion != null) drawCoElutionRegion(g2);
            if(e.getX() < this.coElutionStart.getX()) return;
            // Erase the previous zoom rectangle (if any)...
            Rectangle2D scaledDataArea = _cp.getScreenDataArea(
                    (int) this.coElutionStart.getX(), (int) this.coElutionStart.getY());

                // selected rectangle shouldn't extend outside the data area...
            double xmax = Math.min(e.getX(), scaledDataArea.getMaxX());
            double ymax = Math.min(e.getY(), scaledDataArea.getMaxY());
/*
            this.coElutionRegion = new Rectangle2D.Double(
                        this.coElutionStart.getX(), this.coElutionStart.getY(),
                        xmax - this.coElutionStart.getX(), ymax - this.coElutionStart.getY());
*/
            this.coElutionRegion = new Rectangle2D.Double(
                        this.coElutionStart.getX(), scaledDataArea.getMinY(),
                        Math.abs(e.getX()-coElutionStart.getX()), scaledDataArea.getHeight());

            // Draw the new zoom rectangle...
            drawCoElutionRegion(g2);

            g2.dispose();

        } else {
            _cp.mouseDragged(e);
        }
    }

    public void mouseMoved(MouseEvent e) {
        if((e.isShiftDown() || e.getButton()==MouseEvent.BUTTON3) || shifted) {
            _cp.mouseMoved(e);
        } else {
            _cp.mouseMoved(e);
        }
    }
}
