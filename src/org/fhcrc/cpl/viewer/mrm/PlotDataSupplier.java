package org.fhcrc.cpl.viewer.mrm;

import org.jfree.data.xy.XYSeries;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 19, 2007
 * Time: 12:05:16 PM
 * To change this template use File | Settings | File Templates.
 */
public interface PlotDataSupplier {
    public XYSeries getGraphData();
    public void setGraphData(XYSeries xys);
    public Paint getGraphColor();
    public void setGraphColor(Paint p);
    public boolean isContinuous();
}
