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
