package org.fhcrc.cpl.viewer.util;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Jan 16, 2008
 * Time: 12:33:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class ElutionDataPoint {
    public enum timeUnits {Seconds,Minutes,Hours,Percent,Undetermined};
    public enum intensityUnits {Counts, Percent, Undetermined};

    protected Double time;
    protected Double intensity;
    protected timeUnits tunits;
    protected intensityUnits iunits;

    public Double getTime() {
        return time;
    }

    public void setTime(Double time) {
        this.time = time;
    }

    public Double getIntensity() {
        return intensity;
    }

    public void setIntensity(Double intensity) {
        this.intensity = intensity;
    }

    public timeUnits getTunits() {
        return tunits;
    }

    public void setTunits(timeUnits tunits) {
        this.tunits = tunits;
    }

    public intensityUnits getIunits() {
        return iunits;
    }

    public void setIunits(intensityUnits iunits) {
        this.iunits = iunits;
    }

    public ElutionDataPoint(Double t, Double i) {
        this.setTime(t);
        this.setIntensity(i);
        this.setTunits(timeUnits.Seconds);
        this.setIunits(intensityUnits.Counts);
    }

}
