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

import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;

import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 29, 2007
 * Time: 4:32:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class NoiseSubtractionElutionCurveStrategy extends ElutionCurveStrategy {

    private int noiseSearchWindowWidth = 3;
    public int getNoiseSearchWindowWidth() {
        return noiseSearchWindowWidth;
    }
    public void setNoiseSearchWindowWidth(int w){
        this.noiseSearchWindowWidth = w;
    }
    public ElutionCurve getBestParentCurve() {
        return this.bestParentCurve;
    }

    public ElutionCurve getBestDaughterCurve() {
        return this.bestDaughterCurve;
    }
     
    public void setBestParentCurve(ElutionCurve ec) {
        this.bestParentCurve = ec;
    }

    public void setBestDaughterCurve(ElutionCurve ec) {
        this.bestDaughterCurve = ec;
    }

    public boolean isBestParentCurve(ElutionCurve ec) {
        return ec == getBestParentCurve();
    }

    public boolean isBestDaughterCurve(ElutionCurve ec) {
        return ec == getBestDaughterCurve();
    }

    public MRMTransition getParent() {
        return this.parent;
    }

    public void setParent(MRMTransition p) {
        this.parent = p;
    }

    public MRMDaughter getDaughter() {
        return this.daughter;
    }

    public void setDaughter(MRMDaughter d) {
        this.daughter = d;
    }

    public List<ElutionCurve> getParentCurves() {
        return this.parentCurves;
    }

    public void setParentCurves(List<ElutionCurve> pc) {
        this.parentCurves = pc;
    }

    public List<ElutionCurve> getDaughterCurves() {
         return this.daughterCurves;
     }

     public void setDaughterCurves(List<ElutionCurve> dc) {
         this.daughterCurves = dc;
     }

    public NoiseSubtractionElutionCurveStrategy(MRMTransition p, MRMDaughter d) {
        super(p,d);
    }

    public NoiseSubtractionElutionCurveStrategy() {};


    public double getNoisePctOfHighestPeak() {
        return noisePctOfHighestPeak;
    }

    public void setNoisePctOfHighestPeak(double noisePctOfHighestPeak) {
        this.noisePctOfHighestPeak = noisePctOfHighestPeak;
    }

    protected double noisePctOfHighestPeak = 1.0;

    //very primitive.
    public double noiseLevel(PlotDataSupplier pds) {
        double retVal = 0.0;
        if(pds.getGraphData() == null || pds.getGraphData().getItemCount() == 0) return retVal;
        //return highestPeak(pds)*getNoisePctOfHighestPeak()/100.0;
        return getParent().getMaxYofAllDaughters()*getNoisePctOfHighestPeak()/100.0;
    }

    public double highestPeak(PlotDataSupplier pds) {
        double retVal = 0.0;
        if(pds.getGraphData() == null || pds.getGraphData().getItemCount() == 0) return retVal;
        for(Object oitem: pds.getGraphData().getItems()) {
            double curY =  ((XYDataItem)oitem).getY().doubleValue();
            if(curY > retVal){
                retVal = curY;
            }
        }
        return retVal;
    }

    public ElutionCurve bestDaughterCurve(double noiseFloor,PlotDataSupplier pds){
        return null;
    }

    public double meanWindow(ArrayList<XYDataItem> data, int startAt) {
        double retVal = 0.0;
        int windowWidth = noiseSearchWindowWidth;
        int i = 0;
        int counter = 0;
        int iters = 0;
        for(i = Math.max((int)(startAt-((windowWidth+0.5)/2)),0); i<(data.size()) && (iters<=windowWidth); i++) {
            double tmp = data.get(i).getY().doubleValue();
            iters++;
            if(tmp > 0.0) {
               counter++;
               retVal += data.get(i).getY().doubleValue();
            }
        }
        retVal = retVal/(double)(counter);
        return retVal;
    }


    ElutionCurve findCurveStartingAt(PlotDataSupplier pds, int startIndex) {
       ElutionCurve retVal = null;
       if(pds.getGraphData() == null || pds.getGraphData().getItemCount() == 0) return retVal;
       if(pds.getGraphData().getItemCount() <= startIndex) return retVal; 
       XYSeries workingSeries = pds.getGraphData();
       ArrayList<XYDataItem> allPts = new ArrayList(workingSeries.getItems());
       float noise[][];
       int tmpStartIndex = startIndex;
        MRMTransition mrmt = null;
        if(pds instanceof MRMTransition) {
           mrmt = (MRMTransition)pds;
        } else {
            if(pds instanceof MRMDaughter) {
                mrmt = ((MRMDaughter)pds).getPrecursor();
            }
        }
       if(!pds.isContinuous() && pds instanceof MRMDaughter) {
          MRMDaughter origD= (MRMDaughter)pds;
          MRMDaughter tmpDaughter = new MRMDaughter(origD.getMeanMz(),origD.getLowMz(),origD.getHighMz(),origD.getStartScan(),origD.getEndScan(),origD.getPrecursor());
          workingSeries = origD.getContinDaughterData();
          allPts=new ArrayList(workingSeries.getItems());
          tmpDaughter.setGraphData(workingSeries);
          tmpStartIndex=Utils.findIndexLEXvalue(tmpDaughter,origD.getGraphData().getDataItem(startIndex).getX().floatValue());
       }

        noise = Utils.XYSeriesToArray(workingSeries);
        float tmp[][]=Utils.XYSeriesToArray(workingSeries);
        Spectrum.RemoveBackground(noise);
        for(int jj=0; jj<noise[0].length;jj++){
            noise[0][jj] =  tmp[0][jj];
            noise[1][jj] = Math.max(0.0f,tmp[1][jj]-noise[1][jj]);
        };

       //At least 4 points to make a curve
       if(Utils.allYsGEArr(workingSeries,noise)<4) return retVal;
       double terminus = mrmt.getElutionRegionEnd();

       //find curve starting point
       int curveCounter = 0;
       int curveStart = -1;
       for(int i = tmpStartIndex; i<allPts.size()  ; i++) {
           if(terminus != -1.0d && allPts.get(i).getX().doubleValue()>=terminus)break;
           double curMean = meanWindow(allPts,i);
           if(curMean > noise[1][i]) {
               if(curveStart == -1) curveStart = i;
               curveCounter++;
               //very primitive
               if(curveCounter == 3)break;
           } else {
               curveStart = -1;
           }
       }
       if(curveStart == -1) return null;
       retVal = new ElutionCurve(pds,this);
       retVal.setMinArrayIndex(curveStart);
       retVal.setStrategy(this);
       retVal.setMinElutionTimeSecs(allPts.get(curveStart).getX().doubleValue());
       int startPointInOrigArray = Utils.findIndexLEXvalue(pds,workingSeries.getDataItem(curveStart).getX().floatValue());
       retVal.setMinOrigArrayIndex(startPointInOrigArray);

       //find curve ending point
       curveCounter = 0;
       int curveEnd = -1;

       for(int i = curveStart+1; i<allPts.size();i++){
          if(terminus != -1.0d && allPts.get(i).getX().doubleValue()>=terminus){
              curveEnd = i;
              break;
          }
          double curMean = meanWindow(allPts,i);
          if(curMean <= noise[1][i]) {
              if(curveEnd == -1) curveEnd = i;
              curveCounter++;
              if(curveCounter == 3)break;
          } else {
              curveEnd = -1;
          }
       }
       if(curveEnd == -1) curveEnd = allPts.size()-1;
       int endPointInOrigArray = Utils.findIndexLEXvalue(pds,workingSeries.getDataItem(curveEnd).getX().floatValue());
       retVal.setMaxOrigArrayIndex(endPointInOrigArray);

       retVal.setMaxArrayIndex(curveEnd);
       retVal.setMaxElutionTimeSecs(allPts.get(curveEnd).getX().doubleValue());
       ArrayList<Line2D.Double> segments = new ArrayList<Line2D.Double>();
       XYSeries graphRegion = new XYSeries(getDaughter().getName()+"_ecurve");

       /* data tracing */

       graphRegion.add(allPts.get(curveStart).getX(),0.0);
       segments.add(
          new Line2D.Double(
             allPts.get(curveStart).getX().doubleValue(),
             0.0D,
             allPts.get(curveStart).getX().doubleValue(),
             allPts.get(curveStart).getY().doubleValue()
          )
       );
       for(int i = curveStart; i<curveEnd; i++) {
           if(allPts.get(i).getY().doubleValue() <= 0.0D || ((i<(curveEnd-1))) && allPts.get(i+1).getY().doubleValue()<=0.0D) continue;
           segments.add(
               new Line2D.Double(
                  allPts.get(i).getX().doubleValue(),
                  allPts.get(i).getY().doubleValue(),
                  allPts.get(i+1).getX().doubleValue(),
                  allPts.get(i+1).getY().doubleValue()
               )
           );
           graphRegion.add(allPts.get(i).getX(),allPts.get(i).getY());
       }
       segments.add(
           new Line2D.Double(
              allPts.get(curveEnd).getX().doubleValue(),
              0.0D,
              allPts.get(curveEnd).getX().doubleValue(),
              allPts.get(curveEnd).getY().doubleValue()
           )
        );

       float meanNoise = 0;
       for(int i = curveStart;i<=curveEnd;i++)meanNoise+= noise[1][i];
       meanNoise = meanNoise/(curveEnd-curveStart+1);

        for(int i = curveStart; i<curveEnd; i++) {
            if(allPts.get(i).getY().doubleValue() < meanNoise && allPts.get(i+1).getY().doubleValue() < meanNoise) continue;
            segments.add(
                new Line2D.Double(
                   allPts.get(i).getX().doubleValue(),
                   meanNoise,
                   allPts.get(i+1).getX().doubleValue(),
                   meanNoise
                )
            );
        }

       retVal.setSegments(segments);
       retVal.setGraphRegion(graphRegion);

       float bglevels[][] = new float[2][curveEnd-curveStart+1];
       for(int i = curveStart; i<= curveEnd; i++) {
           bglevels[0][i-curveStart] = noise[0][i];
           bglevels[1][i-curveStart] = noise[1][i];
       }
       retVal.setBackgroundLevels(bglevels);
       retVal.setBackgroundLevel(meanNoise);
       retVal.calculateHighPoints();
       retVal.calculateCMPoints();
       calculateAUC(retVal);
       return retVal;
    }

    public List<ElutionCurve> calculateElutionCurves(PlotDataSupplier pds) {
        Vector<ElutionCurve> retVal = new Vector<ElutionCurve>();
        if(pds.getGraphData() == null || pds.getGraphData().getItemCount() == 0) return retVal;
        //Scan pds left to right by noisesearchwindow means
        MRMTransition mrmt = null;
        if(pds instanceof MRMTransition) {
           mrmt = (MRMTransition)pds;
        } else {
            if(pds instanceof MRMDaughter) {
                mrmt = ((MRMDaughter)pds).getPrecursor();
            }
        }
        int curCurveIndex = 0;
        if(mrmt != null && mrmt.getElutionRegionStart() != -1)
            curCurveIndex = Utils.findIndexLEXvalue(pds,(float)mrmt.getElutionRegionStart());
        if(curCurveIndex == -1) return retVal;
        ElutionCurve curElute = null;
        int maxIndex = 1000000000;
        if(mrmt != null && mrmt.getElutionRegionEnd() != -1) {
            maxIndex = Utils.findIndexLEXvalue(pds,(float)mrmt.getElutionRegionEnd());
        }
        for(;;) {
           curElute = findCurveStartingAt(pds,curCurveIndex);
           if(curElute != null){
               retVal.add(curElute); 
           } else {
               break;
           }
           curCurveIndex = curElute.getMaxOrigArrayIndex() + noiseSearchWindowWidth;
           //curCurveIndex = curElute.getMaxArrayIndex() + noiseSearchWindowWidth;
           if(curCurveIndex >= maxIndex) break;
        }
        return retVal;
    }

    public List<ElutionCurve> calculateParentElutionCurves(PlotDataSupplier pds) {
       PlotDataSupplier defaultPds = (pds == null) ? getParent() : pds;
       List<ElutionCurve> retVal = calculateElutionCurves(defaultPds);
       setParentCurves(retVal);
       return retVal;
    }

    public List<ElutionCurve> calculateDaughterElutionCurves(PlotDataSupplier pds) {
        PlotDataSupplier defaultPds = (pds == null) ? getDaughter() : pds;
        List<ElutionCurve> retVal = calculateElutionCurves(defaultPds);
        setDaughterCurves(retVal);
        return retVal;
     }
          
    public ElutionCurve bestParentCurve(double noiseFloor,PlotDataSupplier pds){
        return null;
    }

    public ElutionCurve maxAUC(List<ElutionCurve> curves) {
        ElutionCurve retVal = null;
        double maxAUC = -1.0;
        if(curves == null || curves.size()==0) return retVal;
        for(ElutionCurve ec: curves) {
            double curAUC = ec.getAUC();
            if(curAUC > maxAUC) {
                maxAUC = curAUC;
                retVal = ec;
            }
        }
        return retVal;
    }

    public void calculateBestCurves() {
       //If there is no precursor scan, choose curve with the highest AUC
       //temporary
       //if(getParentCurves() == null || getParentCurves().isEmpty()) {
          setBestDaughterCurve(maxAUC(getDaughterCurves()));
       //} else {
          //temporary
          setBestParentCurve(maxAUC(getParentCurves()));
       //}
    }

    public void calculateAUC(ElutionCurve ec){
        ec.setAUC(-1.0);
        double sum = 0;
        if(ec.getGraphRegion() == null || ec.getGraphRegion().getItemCount()==0) {
            return;
        }
        float noise[][] = ec.getBackgroundLevels();
        double meanNoise = ec.getBackgroundLevel();
        for(int i = 0; i<ec.getGraphRegion().getItemCount(); i++){
            double bg = meanNoise;
            if(noise!=null)bg = noise[1][i];
            double curY = ec.getGraphRegion().getY(i).doubleValue();
            sum += Math.max((curY-bg),0);
        }
        ec.setAUC(sum);
    }
}
