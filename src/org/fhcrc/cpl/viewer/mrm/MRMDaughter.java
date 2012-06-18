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

import org.fhcrc.cpl.viewer.util.ElutionDataPoint;
import org.jfree.data.xy.XYSeries;

import java.awt.*;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.TreeMap;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Dec 5, 2006
 * Time: 11:35:35 AM
 * To change this template use File | Settings | File Templates.
 */
public class MRMDaughter implements PlotDataSupplier {

    protected float meanMz;
    protected float lowMz;
    protected float highMz;

    protected MRMTransition precursor;
    protected boolean continuous;
    protected Paint daughterPlotColor;


    protected TreeMap<Integer, ElutionDataPoint> scanVals;
    protected int startScan;
    protected int endScan=0;
    protected String name;
    protected XYSeries daughterData;
    protected XYSeries continDaughterData;
    protected int elutionDataTableRow;
    protected ElutionCurve bestElutionCurve;
    protected float quality;

    public void setQuality(float q) {
        this.quality = q;
    }

    public float getQuality() {
        return this.quality;
    }

    public void calculateQuality() {
        this.setQuality(-1.0f);
    }

    public Paint getDaughterPlotColor() {
        return daughterPlotColor;
    }

    public void setDaughterPlotColor(Paint daughterPlotColor) {
        this.daughterPlotColor = daughterPlotColor;
    }

    public XYSeries getContinDaughterData() {
        return continDaughterData;
    }

    public void setContinDaughterData(XYSeries continDaughterData) {
        this.continDaughterData = continDaughterData;
    }

    public int getElutionDataTableRow() {
        return elutionDataTableRow;
    }

    public void setElutionDataTableRow(int elutionDataTableRow) {
        this.elutionDataTableRow = elutionDataTableRow;
    }

    public ElutionCurve getBestElutionCurve() {
        return bestElutionCurve;
    }

    public void setBestElutionCurve(ElutionCurve bestElutionCurve) {
        this.bestElutionCurve = bestElutionCurve;
    }

    public TreeMap<Integer, ElutionDataPoint> getScanVals() {
        return scanVals;
    }

    public void setScanVals(TreeMap<Integer, ElutionDataPoint> scanVals) {
        this.scanVals = scanVals;
    }


    public MRMDaughter(float meanMz, float lowMz, float highMz,
                         int startScan, int endScan, MRMTransition precursor)
    {
        init();
        setMeanMz(meanMz);
        setLowMz(lowMz);
        setHighMz(highMz);
        setPrecursor(precursor);
        setStartScan(startScan);
        setEndScan(endScan);
        setPrecursor(precursor);
        setContinuous(false);
        setBestElutionCurve(null);
        setQuality(-1.0f);
    }

    protected void init() {
       scanVals = new TreeMap<Integer,ElutionDataPoint>();
       continuous = false;
       elutionCurves = new LinkedList<ElutionCurve>();
    }

    public int getStartScan()
    {
        return startScan;
    }

    public void setStartScan(int startScan)
    {
        this.startScan = startScan;
    }

    public int getEndScan()
    {
        return endScan;
    }

    public void setEndScan(int endScan)
    {
        this.endScan = endScan;
    }

    public float getMeanMz()
    {
        return meanMz;
    }

    public void setMeanMz(float meanMz)
    {
        this.meanMz = meanMz;
    }

    public float getLowMz()
    {
        return lowMz;
    }

    public void setLowMz(float lowMz)
    {
        this.lowMz = lowMz;
    }

    public float getHighMz()
    {
        return highMz;
    }

    public void setHighMz(float highMz)
    {
        this.highMz = highMz;
    }

    public void setGraphData(XYSeries xys)
    {
        this.daughterData = xys;
    }

    public String getName()
    {
        return this.name;
    }

    public void setName(String n)
    {
        this.name = n;
    }
    
    public MRMTransition getPrecursor() {
        return precursor;
    }

    public void setPrecursor(MRMTransition precursor) {
        this.precursor = precursor;
    }

    public String toString() {
        return "MRMDaughter: meanMz=" + meanMz + ", precursorMz=" + precursor.getPrecursorMz() +
                ", lowMz=" + lowMz + ", highMz=" + highMz +
                ", startScan=" + startScan + ", endScan=" + endScan;
    }

    public void addScanVal(int scanNumber,ElutionDataPoint edp)
    {
        if (!getScanVals().containsKey(scanNumber))
        {
            scanVals.put(scanNumber,edp);
            if (scanNumber > endScan)
            {
                setEndScan(scanNumber);
            }
        }
    }
   

    /**
     * Can't use toArray() because we're converting Integers to ints
     * @return
     */
    public Integer[] getScanNumbers()
    {
        return (Integer []) getScanVals().keySet().toArray();

    }

    public int numScans()
    {
        return getScanVals().size();
    }

    public boolean isContinuous()
    {
        return this.continuous;
    }

    public void setContinuous(boolean contin)
    {
        this.continuous = contin;
    }

    public Paint getGraphColor()
    {
        return this.daughterPlotColor;
    }

    public void setGraphColor(Paint c)
    {
        this.daughterPlotColor = c;
    }

    public XYSeries getGraphData()
    {
        return this.daughterData;
    }

    public LinkedList<ElutionCurve> getElutionCurves() {
        return elutionCurves;
    }

    public void setElutionCurves(LinkedList<ElutionCurve> elutionCurves) {
        this.elutionCurves = elutionCurves;
    }

    public XYSeries getDaughterData() {
        return daughterData;
    }

    public void setDaughterData(XYSeries daughterData) {
        this.daughterData = daughterData;
    }

    public XYSeries makeDaughterSeries(MRMDaughter s, boolean contin)
    {
       XYSeries result = new XYSeries(s.getName());

       for (int scanNumber : s.getScanVals().keySet())
       {
          ElutionDataPoint edp = s.getScanVals().get(scanNumber);
          if(contin)
          {
             result.add(edp.getTime(),edp.getIntensity());
          }
          else
          {
             result.add((Number)edp.getTime(),0d);
             result.add(edp.getTime(),edp.getIntensity());
             result.add((Number)edp.getTime(),0d);
          }
       }
       return result;
    }

    public XYSeries makeDaughterSeries(MRMDaughter s) {
        return(makeDaughterSeries(s,s.isContinuous()));
    }

    public XYSeries makeDaughterSeries() {
        return(makeDaughterSeries(this,this.isContinuous()));
    }


    LinkedList<ElutionCurve> elutionCurves;

    public static class MeanMzAscComparator implements Comparator<MRMDaughter> {

        public int compare(MRMDaughter o1, MRMDaughter o2)
        {
            float mzDiff = o1.getMeanMz() - o2.getMeanMz();
            return mzDiff == 0 ? 0 : mzDiff > 0 ? 1 : -1;
        }
    }
        public static class PrecursorMzComparator implements Comparator<MRMDaughter>
    {
        public int compare(MRMDaughter o1, MRMDaughter o2)
        {
            float mzDiff = o1.getPrecursor().getPrecursorMz() - o2.getPrecursor().getPrecursorMz();
            return mzDiff == 0 ? 0 : mzDiff > 0 ? 1 : -1;
        }
    }

}
