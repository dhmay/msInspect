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

import java.lang.reflect.Constructor;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 29, 2007
 * Time: 2:12:21 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ElutionCurveStrategy {
    public ElutionCurveStrategy(MRMTransition p, MRMDaughter d) {
        setParent(p);
        setDaughter(d);
    }

    public ElutionCurveStrategy(){};

    public static ElutionCurveStrategy getInstance(MRMTransition p, MRMDaughter d, Class ecurvestrat)
    {
       try  {
           Constructor cons = ecurvestrat.getConstructor();
           ElutionCurveStrategy retVal =
              (ElutionCurveStrategy)cons.newInstance();
           retVal.setParent(p);
           retVal.setDaughter(d);
           return retVal;
       }   catch (Exception x)
       {
                   x.printStackTrace();
                   System.exit(1);
                   return null; // make compiler happy
       }
   }
 
    protected MRMTransition parent;

    public abstract MRMTransition getParent();

    public abstract void setParent(MRMTransition _parent);

    public abstract MRMDaughter getDaughter();

    public abstract void setDaughter(MRMDaughter _daughter);

    protected MRMDaughter daughter;

    public abstract int getNoiseSearchWindowWidth();

    public abstract void setNoiseSearchWindowWidth(int noiseSearchWindowWidth);

    protected int noiseSearchWindowWidth = 3;

    public abstract double noiseLevel(PlotDataSupplier pds);

    public abstract double highestPeak(PlotDataSupplier pds);

    public abstract List<ElutionCurve> getParentCurves();

    public abstract void setParentCurves(List<ElutionCurve> parentCurves);

    public abstract List<ElutionCurve> getDaughterCurves();

    public abstract void setDaughterCurves(List<ElutionCurve> daughterCurves);

    protected List<ElutionCurve> parentCurves = null;
    protected List<ElutionCurve> daughterCurves = null;

    public abstract List<ElutionCurve> calculateParentElutionCurves(PlotDataSupplier pds);
    
    public abstract List<ElutionCurve> calculateDaughterElutionCurves(PlotDataSupplier pds);

    protected ElutionCurve bestDaughterCurve;
    public abstract ElutionCurve getBestDaughterCurve();
    public abstract void setBestDaughterCurve(ElutionCurve ec);

    protected ElutionCurve bestParentCurve;
    public abstract ElutionCurve getBestParentCurve();
    public abstract void setBestParentCurve(ElutionCurve ec);

    public abstract void calculateBestCurves();
    public boolean isBestParentCurve(ElutionCurve ec) {
         return ec == getBestParentCurve();
     }

     public boolean isBestDaughterCurve(ElutionCurve ec) {
         return ec == getBestDaughterCurve();
     }

     public abstract void calculateAUC(ElutionCurve ec);

}
