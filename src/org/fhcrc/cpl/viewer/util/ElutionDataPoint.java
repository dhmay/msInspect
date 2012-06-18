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
