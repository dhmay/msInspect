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
import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.fhcrc.cpl.viewer.mrm.utilities.PeaksTableModel;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.DrawingSupplier;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.*;
import java.awt.geom.Line2D;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: May 10, 2007
 * Time: 11:42:06 AM
 * To change this template use File | Settings | File Templates.
 */
public class Utils {
    /**
    * For a given scan, within an MZ range, feed back the highest peak found.
    *
    * TODO: This could be faster, if necessary: we could do a binary search and then expand up and down.
    *  @param scan
    * @param minMZ
    * @param maxMZ
    * @return
    */
   public static double getMaxIntensityForScan(MSRun.MSScan scan,
                                               double minMZ,
                                               double maxMZ)
   {
       double maxIntensityThisScan = 0;

       float[][] data = scan.getSpectrum();

       for (int j = 0; j < data[0].length; j++)
       {
           //since MZ values are monotonically increasing with j, when we
           // hit the max we can stop.
           if (data[0][j] > maxMZ)
               break;
           if (data[0][j] >= minMZ)
           {
               maxIntensityThisScan =
                       Math.max(maxIntensityThisScan, data[1][j]);
           }
       }
//_log.debug("       max intensity: scanIndex " + scanIndex + ", scanNumber " + run.getScan(scanIndex).getNum() + ", maxInt=" + maxIntensityThisScan);

       return maxIntensityThisScan;
      }

    private static final int PALE_CONSTANT=100;
    public static Color paleColor(Color origColor) {

        int origRed = origColor.getRed();
        int origGreen = origColor.getGreen();
        int origBlue = origColor.getBlue();

        int maxColorVal = Math.max(Math.max(origRed,origGreen),origBlue);
        if(maxColorVal+PALE_CONSTANT > 255) {
            int adjust = (maxColorVal+PALE_CONSTANT) - 255;
            origRed = Math.max(0,origRed-adjust);
            origGreen = Math.max(0,origGreen-adjust);
            origBlue = Math.max(0,origBlue-adjust);
        }

        Color paleColor = new Color(origRed+PALE_CONSTANT,origGreen+PALE_CONSTANT,origBlue+PALE_CONSTANT);
        return paleColor;
    }

     public static DrawingSupplier plainDrawingSupplier(Color c)
    {
        Paint black[] = new Paint[1];
        black[0] = c;
        return new DefaultDrawingSupplier(
              black,
              DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
              DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
              DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
              DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE
        );
     }

    public static boolean allYsAreLEVal(PlotDataSupplier pds, double val) {
         XYSeries xys = pds.getGraphData();
         int i = 0;
         for(Object oitem: xys.getItems()) {
             i++;
             if(((XYDataItem)oitem).getY().doubleValue() > val){
                 return false;
             }
         }
         return true;
     }

     public static boolean allYsAre0(PlotDataSupplier pds) {
         return allYsAreLEVal(pds,0.0);
     }

     public static int allYsGEVal(XYSeries xys, double val) {
        int retVal = 0;
        for(Object oitem: xys.getItems()) {
           if(((XYDataItem)oitem).getY().doubleValue() > val){
              retVal++;
           }
        }
        return retVal;
     }

    public static int allYsGEArr(PlotDataSupplier pds, float arr[][]) {
        int retVal = 0;
        XYSeries xys = pds.getGraphData();
        float xysarr[][] = PDStoArray(pds);
        for(int i=0; i<xysarr[0].length; i++) {
           if(xysarr[1][i] > arr[1][i]){
              retVal++;
           }
        }
        return retVal;
     }

    public static int allYsGEArr(XYSeries xys, float arr[][]){
       int retVal = 0;
       for(int i = 0; i<xys.getItemCount(); i++) {
           if(xys.getY(i).floatValue() > arr[1][i]) {
               retVal++;
           }
       }
       return retVal;
    }

    public static Polygon legendThing(int width, int height)
    {
       int[] xes = new int[100];
       int[] yes = new int[100];
       int leftX = -1*(width/2);
       int bottomY = -1*(height/2);
       int pointCounter = 0;
       xes[0] = leftX;
       yes[0] = bottomY;
       int curX = leftX;
       int curY = bottomY;
       while (curY <= height/2) {
           xes[pointCounter] = curX;
           yes[pointCounter] = curY;
           pointCounter++;
           curX=width/2;
           xes[pointCounter] = curX;
           yes[pointCounter] = curY;
           pointCounter++;
           curY++;
           xes[pointCounter] = curX;
           yes[pointCounter] = curY;
           pointCounter++;
           curX=leftX;
       }
       return new Polygon(xes,yes,pointCounter-1);
    }

    public static float[][] PDStoArray(PlotDataSupplier pds)
    {
       float retval[][] = null;
       if(pds == null) return retval;
       XYSeries xys = pds.getGraphData();
       int itemCount = xys.getItemCount();
       retval = new float[2][itemCount];
       for (int i = 0; i < itemCount; i++) {
          retval[0][i] = xys.getX(i).floatValue();
          Number y = xys.getY(i);
          if (y != null) {
             retval[1][i] = y.floatValue();
          }  else {
             retval[1][i] = Float.NaN;
          }
       }
       return retval;
    }

    public static float[][] XYSeriesToArray(XYSeries xys) {
/*
       float retval[][] = null;
       if(xys == null) return retval;
       int itemCount = xys.getItemCount();
       retval = new float[2][itemCount];
       for (int i = 0; i < itemCount; i++) {
          retval[0][i] = xys.getX(i).floatValue();
          Number y = xys.getY(i);
          if (y != null) {
             retval[1][i] = y.floatValue();
          }  else {
             retval[1][i] = Float.NaN;
          }
       }
*/     double tmp[][] = xys.toArray();
       float retval[][] = new float[2][tmp[0].length];
       for(int i = 0; i<2; i++)
          for(int j= 0; j<tmp[0].length; j++)
             retval[i][j] = (float)tmp[i][j];
       return retval; 
    }

    public static XYSeries ArrayToXYSeries(float farr[][], Comparable name) {
        XYSeries retval = new XYSeries(name);
        if(farr == null || farr[0].length == 0 || farr.length != 2 ) return retval;
        for(int i = 0; i<farr[0].length; i++) {
            retval.add(farr[0][i],farr[1][i],false);
        }
        return retval;
    }


    public static void ArrayRefillPDS(PlotDataSupplier pds, float newVals[][])
    {
       if(pds == null || newVals == null) return;
       XYSeries xys = pds.getGraphData();
       xys.clear();
       for(int i = 0; i<newVals[0].length; i++)
          xys.add(newVals[0][i],newVals[1][i],false);
    }

    public static XYLineAnnotation line2Annotation(Line2D.Double l2dd, Stroke s, Paint p){
        return new XYLineAnnotation(l2dd.getX1(), l2dd.getY1(),  l2dd.getX2(), l2dd.getY2(),s,p);
    }

    public static int getMinScanForPlot(MRMTransition plot)
    {
       if(plot == null) return 1;
       int min = 10000000;
       for(MRMDaughter d :  plot.getDaughters().values()) {
          if(d.getStartScan() < min) min = d.getStartScan();
       }
       return min;
    }

    public static int getMaxScanForPlot(MRMTransition plot)
    {
       if(plot == null)return 100;
       int max = 0;
       for(MRMDaughter d: plot.getDaughters().values()){
          if(d.getEndScan() > max) max = d.getEndScan();
       }

       return max;
    }

    public static double getWeightedAverageDaughtersTime(XYSeriesCollection daughterSet) {
       double retVal = 0.0;
       double yTot=0.0;
       double xWeight=0.0;
       for(Object xyso: daughterSet.getSeries()){
          XYSeries xys = (XYSeries)xyso;
          for(Object xydio: xys.getItems()){
              XYDataItem xydi = (XYDataItem)xydio;
              xWeight += (xydi.getX().doubleValue()*xydi.getY().doubleValue());
              yTot += xydi.getY().doubleValue();
          }
       }
       if(yTot > 0) retVal = xWeight/yTot;
       return retVal;
    }

    public static int findIndexLEXvalue(PlotDataSupplier pds,float xvalue) {
        XYSeries xys = pds.getGraphData();
        if(xys.getItems().size() == 0) return -1;
        float curx=-1;
        int i;
        for(i = 0; i<xys.getItemCount(); i++){
            curx = xys.getDataItem(i).getX().floatValue();
            if(curx >= xvalue) break;
        }
        if(i == 0 && curx > xvalue) return 0;
        if(curx == xvalue) return i;
        return i-1;
    }

    public static boolean qualColIsEmpty() {
       boolean retVal = true;
       for(int i=0; i< ((PeaksTableModel)(MRMDialog.peaksTable.getModel())).data.length; i++){
          Object curQual = ((PeaksTableModel)(MRMDialog.peaksTable.getModel())).data[i][MRMDialog.peaksData.Quality.colno];
          if(curQual != null && (curQual instanceof Float && ((Float)curQual).floatValue() != -1.0F)) {
             return false;
          }
       }
       return retVal;
    }


}