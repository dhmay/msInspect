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

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;

import java.awt.*;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
     * Model a single transition
 */
public class MRMTransition implements PlotDataSupplier{
    public static final Paint COLOR_SERIES[] =
    /*   DefaultDrawingSupplier.DEFAULT_PAINT_SEQUENCE; */
               {
          Color.BLUE, Color.GREEN, Color.MAGENTA,
          Color.RED,  Color.ORANGE,Color.PINK
       };

    protected float precursorMz;
    protected String name;
    protected MSRun run;
    protected XYSeries parentData;
    protected int currentDaughterIndex = 0;
    protected TreeMap<Float,MRMDaughter> Daughters;

    public boolean isContinuous() {
        return true;
    }

    public HashMap<MRMDaughter, ElutionCurveStrategy> getElutionCurves() {
        return elutionCurves;
    }

    public void setElutionCurves(HashMap<MRMDaughter, ElutionCurveStrategy> elutionCurves) {
        this.elutionCurves = elutionCurves;
    }

    protected HashMap<MRMDaughter,ElutionCurveStrategy> elutionCurves = new HashMap<MRMDaughter,ElutionCurveStrategy>();

    public XYSeries getParentData() {
        return parentData;
    }

    public void setParentData(XYSeries parentData) {
        this.parentData = parentData;
    }

    public Paint getGraphColor() {
        return graphColor;
    }

    public void setGraphColor(Paint graphColor) {
        this.graphColor = graphColor;
    }

    protected Paint graphColor;

    public MRMTransition(float precursorMz, MSRun run)
    {
        init();
        setPrecursorMz(precursorMz);
        setRun(run);
    }

    public String toString()
    {
        String retVal = "MRMTransition: precursorMz " + getPrecursorMz() + ", " +
                "n of daughters: ";
        if(getDaughters() == null) {
            retVal += " 0";
        } else {
            retVal += getDaughters().size() +":";
            for(MRMDaughter mrmd: getDaughters().values()){
               retVal += " " + mrmd.getMeanMz();
            }
        }
        return retVal;

   }

    protected float quality;
    public float getQuality() {
        return quality;
    }

    public void setQuality(float quality) {
        this.quality = quality;
    }

    public void calculateQuality(){
        setQuality(-1.0f);
    }


    protected void init()
    {
        Daughters = new TreeMap<Float,MRMDaughter>();
        graphColor = Color.BLACK;
        quality = -1.0f;
    }

    public float getPrecursorMz()
    {
        return precursorMz;
    }

    public void setPrecursorMz(float precursorMz)
    {
        this.precursorMz = precursorMz;
    }

    public MSRun getRun()
    {
        return run;
    }

    public void setRun(MSRun run)
    {
        this.run = run;
    }

    public void setGraphData(XYSeries xys)
    {
        this.parentData = xys;
    }

    public XYSeries getGraphData()
    {
        return this.parentData;
    }

    public Map<Float, MRMDaughter> getDaughters() {
        return Daughters;
    }

    public void setDaughters(TreeMap<Float, MRMDaughter> daughters) {
        Daughters = daughters;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getCurrentDaughterIndex() {
        return currentDaughterIndex;
    }

    public boolean setCurrentDaughterIndex(int n) {
        if(n < 0 || n >= getDaughters().size()) return false;
        currentDaughterIndex = n;
        return true;
    }

    public MRMDaughter getCurrentDaughter() {
        if(getCurrentDaughterIndex() < 0) return null;
        return getDaughters().values().toArray(new MRMDaughter[0])[getCurrentDaughterIndex()];
    }

    public int getMinScanOfDaughters() {
       int retVal = 10000000;
       for(MRMDaughter d: getDaughters().values()){
          if(retVal > d.getStartScan()) retVal = d.getStartScan();
       }
       return retVal;
    }

    public double getMinTimeOfAllDaughters() {
        double retVal = 1000000;
        for(MRMDaughter d: getDaughters().values()){
           if(retVal > d.continDaughterData.getDataItem(0).getX().doubleValue()) retVal = d.continDaughterData.getDataItem(0).getX().doubleValue();
        }
        return retVal;
    }

    public double getMaxTimeOfAllDaughters() {
         double retVal = -1;
         for(MRMDaughter d: getDaughters().values()){
            if(retVal < d.continDaughterData.getDataItem(d.continDaughterData.getItemCount()-1).getX().doubleValue()) retVal = d.continDaughterData.getDataItem(d.continDaughterData.getItemCount()-1).getX().doubleValue();
         }
         return retVal;
    }

    public double getCalcMaxYAllDaughters() {
        return calcMaxYAllDaughters;
    }

    public void setCalcMaxYAllDaughters(double calcMaxYAllDaughters) {
        this.calcMaxYAllDaughters = calcMaxYAllDaughters;
    }

    public double getCalcXatMaxYAllDaughters() {
        return calcXatMaxYAllDaughters;
    }

    public void setCalcXatMaxYAllDaughters(double calcXatMaxYAllDaughters) {
        this.calcXatMaxYAllDaughters = calcXatMaxYAllDaughters;
    }

    protected double calcMaxYAllDaughters;
    protected double calcXatMaxYAllDaughters;

    public void calcMaxYofAllDaughters(){
        double maxYval = -1;
        double maxXval = -1;
        setCalcMaxYAllDaughters(maxYval);
        setCalcXatMaxYAllDaughters(maxXval);
        for(MRMDaughter d: getDaughters().values()){
           ElutionCurve ec = d.getBestElutionCurve();
           if(ec == null || ec.getHighestPointY() <= 0) continue;
           if(ec.getHighestPointY() > maxYval) {
              maxYval = ec.getHighestPointY();
              maxXval = ec.getHighestPointX();
           }
        }
        setCalcXatMaxYAllDaughters(maxXval);
        setCalcMaxYAllDaughters(maxYval);
    }

    public double getMaxYofAllDaughters() {
        double retVal = -1;
        if(getDaughters() == null || getDaughters().size() == 0) return retVal;
        for(MRMDaughter d: getDaughters().values()){
           XYSeries xys = d.getContinDaughterData();
           if(xys == null  || xys.getItemCount() == 0) {
               continue;
           }
           for(Object o: xys.getItems()) {
               XYDataItem xydi = (XYDataItem)o;
               double curY = xydi.getY().doubleValue();
               if(curY>retVal){
                   retVal = curY;
               }
           }
        }
        return retVal;
    }

    public int getMaxScanOfDaughters(){
        int retVal = 0;
        for(MRMDaughter d: getDaughters().values()){
           if(retVal < d.getEndScan()) retVal = d.getEndScan();
        }
        return retVal;
    }

    protected double elutionRegionStart=-1;
    protected double elutionRegionEnd=-1;

    public int getTableRow() {
        return tableRow;
    }

    public void setTableRow(int tableRow) {
        this.tableRow = tableRow;
    }

    protected int tableRow;


    public double getElutionRegionStart() {
        return elutionRegionStart;
    }

    public void setElutionRegionStart(double elutionRegionStart) {
        this.elutionRegionStart = elutionRegionStart;
    }

    public double getElutionRegionEnd() {
        return elutionRegionEnd;
    }

    public void setElutionRegionEnd(double elutionRegionEnd) {
        this.elutionRegionEnd = elutionRegionEnd;
    }

    public static class PrecursorMzComparator implements Comparator<MRMTransition>
    {
        public int compare(MRMTransition o1, MRMTransition o2)
        {
            float mzDiff = o1.getPrecursorMz() - o2.getPrecursorMz();
            return mzDiff == 0 ? 0 : mzDiff > 0 ? 1 : -1;
        }
    }

    private static final double UNSET_CODE = 100000000000.0;
    public double calculateMinOfAllBestDaughterCurves() {
       double retVal = UNSET_CODE;
        for(MRMDaughter d: getDaughters().values()){
           if(d != null) {
              ElutionCurve ec = d.getBestElutionCurve();
              if(ec != null) {
                 double curLeft = ec.getMinElutionTimeSecs();
                 if(curLeft < retVal) retVal = curLeft;
              }
           }
        }
       if(retVal == UNSET_CODE) retVal = -1; 
       return retVal;
    }

    public double calculateMaxOfAllBestDaughterCurves() {
        double retVal = -1;
        for(MRMDaughter d: getDaughters().values()){
           if(d != null) {
              ElutionCurve ec = d.getBestElutionCurve();
              if(ec != null) {
                 double curRight = ec.getMaxElutionTimeSecs();
                 if(curRight > retVal) retVal = curRight;
              }
           }
        }
        return retVal;
    }
}
