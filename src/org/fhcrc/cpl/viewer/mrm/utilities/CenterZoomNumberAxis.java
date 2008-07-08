package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.event.AxisChangeListener;
import org.jfree.data.Range;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 21, 2007
 * Time: 12:32:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class CenterZoomNumberAxis extends NumberAxis {
     public CenterZoomNumberAxis() {
         super();
     } 

     public CenterZoomNumberAxis(String label) {
         super(label);
     }

    public MRMDialog.domainAxisZoomCoordinator getTheListener() {
        return theListener;
    }

    public void setTheListener(MRMDialog.domainAxisZoomCoordinator theListener) {
        this.theListener = theListener;
    }

    private MRMDialog.domainAxisZoomCoordinator theListener;

    public void addChangeListener(AxisChangeListener listener) {
        super.addChangeListener(listener);
        if(listener.getClass() == MRMDialog.domainAxisZoomCoordinator.class) {
            setTheListener((MRMDialog.domainAxisZoomCoordinator)listener);
        }
    }


    public void zoomRange(double lowerPercent, double upperPercent) {
        double start = this.getRange().getLowerBound()+this.getRange().getLength()*lowerPercent;
        double end =   this.getRange().getLowerBound()+this.getRange().getLength()*upperPercent;
        Range adjusted = null;
        if (isInverted()) {
            adjusted = new Range(1-end,1-start);
        }
        else {
            adjusted = new Range(start,end);
       }
        setRange(adjusted);
        configure();
    }

}
