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
package org.fhcrc.cpl.viewer.feature.extraction;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategy;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.BaseFeatureStrategy;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategyWindow;

import java.util.*;


/**
 * This class exists to bridge the gap between old-school and new-school feature finders.
 * Eventually I'd like to destroy it.  Meantime, any code that needs to invoke feature-
 * finding and might have to use old-school or new-school structures should use this.
 */
public class FeatureFindingBroker
{
    static Logger _log = Logger.getLogger(FeatureFindingBroker.class);

    public static final String OLD_FEATURE_FINDING_PACKAGE_NAME =
            "org.fhcrc.cpl.viewer.feature";
    public static final String NEW_FEATURE_FINDING_PACKAGE_NAME =
            "org.fhcrc.cpl.viewer.feature.extraction.strategy";

    public static FeatureSet findPeptides(MSRun run, int startScan, int scanCount,
                                          int maxCharge, FloatRange mzRange,
                                          int dumpWindowSize,
                                          int accurateMassAdjustmentScans,
                                          Class featureStrategyClass,
                                          boolean writeStatus,
                                          boolean peakRidgeWalkSmoothed,
                                          boolean plotStatistics)
            throws InterruptedException
    {
        return findPeptides(run, startScan, scanCount, maxCharge, mzRange, dumpWindowSize, accurateMassAdjustmentScans,
                featureStrategyClass, writeStatus, peakRidgeWalkSmoothed, plotStatistics, FeatureStrategyWindow.DEFAULT_WINDOW_WIDTH);
    }

    public static FeatureSet findPeptides(MSRun run, int startScan, int scanCount,
                                          int maxCharge, FloatRange mzRange,
                                          int dumpWindowSize,
                                          int accurateMassAdjustmentScans,
                                          Class featureStrategyClass,
                                          boolean writeStatus,
                                          boolean peakRidgeWalkSmoothed,
                                          boolean plotStatistics, int scanWindowSize)
            throws InterruptedException
    {
        float minMz = mzRange.min;
        float maxMz = mzRange.max;
        FeatureSet featureSet = null;

        boolean oldSchoolStrategy = isOldSchoolStrategy(featureStrategyClass);
        if (oldSchoolStrategy)
        {
            FeatureExtractor.setDefault(featureStrategyClass);
            FeatureExtractor find =
                    FeatureExtractor.getDefault(run, startScan, scanCount,
                            maxCharge,
                            new FloatRange(minMz, maxMz), 2.0);
            if (dumpWindowSize > 0)
                find.setDumpWindowSize(dumpWindowSize);
            if (accurateMassAdjustmentScans > 0)
                find.setAccurateMassAdjustmentScans(accurateMassAdjustmentScans);

            if (writeStatus)
            {
                find.setStatusListener(new FeatureExtractor.StatusListener()
                {
                    public void progress(float percent)
                    {
                        ApplicationContext.setMessage(String.valueOf(percent) +
                                "% complete.");
                    }
                });
            }
            featureSet = find.analyze();
            Arrays.sort(featureSet.getFeatures(), new Feature.IntensityDescComparator());
        }
        else
        {
            FeatureFinder featureFinder =
                    new FeatureFinder(run, startScan, scanCount,
                            maxCharge,
                            new FloatRange(minMz, maxMz),
                            featureStrategyClass, plotStatistics, scanWindowSize);
            if (dumpWindowSize > 0)
                featureFinder.setDumpWindowSize(dumpWindowSize);
            if (accurateMassAdjustmentScans > 0)
                featureFinder.setAccurateMassAdjustmentScans(accurateMassAdjustmentScans);
            featureFinder.setPeakRidgeWalkSmoothed(peakRidgeWalkSmoothed);
            if (writeStatus)
            {
                featureFinder.setStatusListener(new BaseFeatureStrategy.StatusListener()
                {
                    public void progress(float percent)
                    {
                        ApplicationContext.setMessage(String.valueOf(percent) +
                                "% complete.");
                    }
                });
            }

            featureSet = featureFinder.findPeptides();
            if (plotStatistics)
                featureFinder.plotStatistics();
        }

        return featureSet;
    }

    public static Class getFeatureStrategyClass(String strategyClassName)
            throws ClassNotFoundException
    {
        String newSchoolStrategyClassName = strategyClassName;
        Class result = null;

        if (!newSchoolStrategyClassName.contains("."))
            newSchoolStrategyClassName = NEW_FEATURE_FINDING_PACKAGE_NAME
                    + "." + strategyClassName;

        try
        {
            result = Class.forName(newSchoolStrategyClassName);
        }
        catch (ClassNotFoundException e)
        {
            {
                if (!strategyClassName.contains("."))
                {

                    String oldSchoolStrategyClassName =
                            OLD_FEATURE_FINDING_PACKAGE_NAME + "." +
                                    strategyClassName;
                    try
                    {
                        result = Class.forName(oldSchoolStrategyClassName);
                    }
                    catch (ClassNotFoundException x)
                    {
                        throw new ClassNotFoundException("Could not load class: " +
                                newSchoolStrategyClassName +
                                " or " + oldSchoolStrategyClassName);
                    }
                }
                else
                    throw new ClassNotFoundException("Could not load class: " +
                            strategyClassName);
            }
        }
        return result;
    }

    public static boolean isOldSchoolStrategy(Class strategyClass)
    {
        if (FeatureStrategy.class.isAssignableFrom(strategyClass))
            return false;
        return true;
    }

}
