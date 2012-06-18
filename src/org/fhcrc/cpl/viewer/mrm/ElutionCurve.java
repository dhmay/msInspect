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
package org.fhcrc.cpl.viewer.mrm;

import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;

import java.awt.geom.Line2D;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 19, 2007
 * Time: 11:47:13 AM
 * To change this template use File | Settings | File Templates.
 */
public class ElutionCurve {

    long minScan;
    public long getMinScan() {
        return minScan;
    }

    public void setMinScan(long minScan) {
        this.minScan = minScan;
    }

    long maxScan;
    public long getMaxScan() {
        return maxScan;
    }

    public void setMaxScan(long maxScan) {
        this.maxScan = maxScan;
    }

    public int getMinArrayIndex() {
        return minArrayIndex;
    }

    public void setMinArrayIndex(int minArrayIndex) {
        this.minArrayIndex = minArrayIndex;
    }

    public int getMaxArrayIndex() {
        return maxArrayIndex;
    }

    public void setMaxArrayIndex(int maxArrayIndex) {
        this.maxArrayIndex = maxArrayIndex;
    }

    int minArrayIndex;
    int maxArrayIndex;
    int minOrigArrayIndex;

    public int getMinOrigArrayIndex() {
        return minOrigArrayIndex;
    }

    public void setMinOrigArrayIndex(int minOrigArrayIndex) {
        this.minOrigArrayIndex = minOrigArrayIndex;
    }

    public int getMaxOrigArrayIndex() {
        return maxOrigArrayIndex;
    }

    public void setMaxOrigArrayIndex(int maxOrigArrayIndex) {
        this.maxOrigArrayIndex = maxOrigArrayIndex;
    }

    int maxOrigArrayIndex;
    double minElutionTimeSecs;

    public double getMinElutionTimeSecs() {
        return minElutionTimeSecs;
    }

    public void setMinElutionTimeSecs(double minElutionTimeSecs) {
        this.minElutionTimeSecs = minElutionTimeSecs;
    }

    double maxElutionTimeSecs;
    public double getMaxElutionTimeSecs() {
        return maxElutionTimeSecs;
    }

    public void setMaxElutionTimeSecs(double maxElutionTimeSecs) {
        this.maxElutionTimeSecs = maxElutionTimeSecs;
    }

    double AUC;
    public double getAUC() {
        return AUC;
    }

    public void setAUC(double AUC) {
        this.AUC = AUC;
    }

    protected double highestPointX;
    public double getHighestPointX() {
        return this.highestPointX;
    }

    public void setHighestPointX(double highestPointX) {
        this.highestPointX = highestPointX;
    }

    double highestPointY;
    public double getHighestPointY() {
        return highestPointY;
    }

    public void setHighestPointY(double highestPointY) {
        this.highestPointY = highestPointY;
    }

    double centerOfMassX;

    public double getCenterOfMassX() {
        return centerOfMassX;
    }

    public void setCenterOfMassX(double centerOfMassX) {
        this.centerOfMassX = centerOfMassX;
    }

    double backgroundLevel;
    public double getBackgroundLevel() {
        return backgroundLevel;
    }

    public void setBackgroundLevel(double backgroundLevel) {
        this.backgroundLevel = backgroundLevel;
    }

    private float backgroundLevels[][];
    public float[][] getBackgroundLevels() {
        return backgroundLevels;
    }

    public void setBackgroundLevels(float backgroundLevels[][]) {
        this.backgroundLevels = backgroundLevels;
    }
    
    PlotDataSupplier parent;

    public PlotDataSupplier getParent() {
        return parent;
    }

    public void setParent(PlotDataSupplier parent) {
        this.parent = parent;
    }

    // constructors
    public ElutionCurve() {
        parent = null;
    }

    public ElutionCurve(PlotDataSupplier p) {
        setParent(p);
    }

    public ElutionCurve(ElutionCurveStrategy s) {
        setParent(null);
        setStrategy(s);
    }

    public ElutionCurve(PlotDataSupplier pds, ElutionCurveStrategy s) {
        setParent(pds);
        setStrategy(s);
    }

    public ArrayList<Line2D.Double> getSegments() {
        return segments;
    }

    public void setSegments(ArrayList<Line2D.Double> segs) {
        this.segments = segs;
    }

    ArrayList<Line2D.Double> segments;

    public ElutionCurveStrategy getStrategy() {
        return strategy;
    }

    public void setStrategy(ElutionCurveStrategy strategy) {
        this.strategy = strategy;
    }

    ElutionCurveStrategy strategy;

    public XYSeries getGraphRegion() {
        return graphRegion;
    }

    public void setGraphRegion(XYSeries graphRegion) {
        this.graphRegion = graphRegion;
    }

    protected XYSeries graphRegion;

    public void calculateHighPoints() {
        //this is truly annoying:  XYSeries don't have a size or isEmpty
        if(getGraphRegion() == null || getSegments().isEmpty()) {
            return;
        }
        double highx=-1,highy=-1;
        for(Object xydo: getGraphRegion().getItems()){
            XYDataItem xydi = (XYDataItem)xydo;
            double curY = xydi.getY().doubleValue();
            if(curY > highy){
                highy = curY;
                highx = xydi.getX().doubleValue();
            }
        }
        setHighestPointX(highx);
        setHighestPointY(highy);
    }

    public void calculateCMPoints(){
        setCenterOfMassX(-1.0);
        if(getGraphRegion() == null || getSegments().isEmpty()) {
            return;
        }
        double sumX=0.0;
        double sumY=0.0;
        for(Object xydo: getGraphRegion().getItems()){
            XYDataItem xydi = (XYDataItem)xydo;
            sumX += xydi.getX().doubleValue()*xydi.getY().doubleValue();
            sumY += xydi.getY().doubleValue();
        }
        setCenterOfMassX(sumX/sumY);
    }


}
