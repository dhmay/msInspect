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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHeatMap;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.io.File;


/**
 * Command linemodule for adduct detection. Quick and dirty.  Lots of hacks.
 */
public class DetectAdductsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(DetectAdductsCLM.class);

    protected File[] ms1FeatureFiles;

    protected double massWavelength =
            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH;
    protected int scanWindowSize = 1;
    protected int maxRelativeDaltons = 100;
    protected int minRelativeDaltons = 0;

    protected float matcherMassTolerance = 10;
    protected int matcherMassToleranceType = FeatureSetMatcher.DELTA_MASS_TYPE_PPM;
    protected int matcherScanTolerance=100;

    protected File outZeroBucketFeatureFileDir = null;

    protected int minRelativeSeconds = -60;
    protected int maxRelativeSeconds = 60;
    protected int secondsIncrement = 1;




    public DetectAdductsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "detectadducts";
        mShortDescription = "Display a chart showing MS1 feature masses at different intervals above MS2 feature masses";
        mHelpMessage = "Display a chart showing MS1 feature masses at different intervals above MS2 feature masses.  Highly experimental";

        CommandLineArgumentDefinition[] argDefs =
            {
                    createUnnamedSeriesFileArgumentDefinition(true,
                            "Feature File of features to interrogate"),
                    new DecimalArgumentDefinition("masswavelength", false,
                            "Mass wavelength", massWavelength),
                    new IntegerArgumentDefinition("scanwindowsize", false,
                            "Size of the scan window (including identity scan)", scanWindowSize),
                    new IntegerArgumentDefinition("maxrelativedaltons", false,
                            "Maximum Daltons, relative to original mass", maxRelativeDaltons),
                    new IntegerArgumentDefinition("minrelativedaltons", false,
                            "Minimum Daltons, relative to original mass", minRelativeDaltons),
                    new DeltaMassArgumentDefinition("deltamass", false,
                            "Mass tolerance around each Dalton increment",
                            new DeltaMassArgumentDefinition.DeltaMassWithType(matcherMassTolerance, matcherMassToleranceType)),
                    new IntegerArgumentDefinition("minrelativeseconds", false,
                            "", minRelativeSeconds),
                    new IntegerArgumentDefinition("maxrelativeseconds", false,
                            "", maxRelativeSeconds),
                    new IntegerArgumentDefinition("secondsincrement", false,
                            "", secondsIncrement),
            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        ms1FeatureFiles =
                this.getUnnamedSeriesFileArgumentValues();


        massWavelength = getDoubleArgumentValue("masswavelength");
        scanWindowSize = getIntegerArgumentValue("scanwindowsize");

        maxRelativeDaltons = getIntegerArgumentValue("maxrelativedaltons");
        minRelativeDaltons = getIntegerArgumentValue("minrelativedaltons");

        if (minRelativeDaltons < 0)
            throw new ArgumentValidationException("What's the point of minrelativedaltons<0?");

        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassAndType =
                getDeltaMassArgumentValue("deltamass");
        matcherMassTolerance = deltaMassAndType.getDeltaMass();
        matcherMassToleranceType = deltaMassAndType.getDeltaMassType();
        minRelativeSeconds = getIntegerArgumentValue("minrelativeseconds");
        maxRelativeSeconds = getIntegerArgumentValue("maxrelativeseconds");
        secondsIncrement = getIntegerArgumentValue("secondsincrement");



    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        int numSecondsIncrements = (maxRelativeSeconds - minRelativeSeconds) / secondsIncrement;
        int numMassIncrements = (maxRelativeDaltons - minRelativeDaltons);

        double minRelativeMass = minRelativeDaltons * massWavelength;
        double maxRelativeMass = maxRelativeDaltons * massWavelength;


        float[] massBuckets = new float[numMassIncrements];
        float[] timeBuckets = new float[numSecondsIncrements];

        for (int i=0; i<numMassIncrements; i++)
            massBuckets[i] = minRelativeDaltons + i;
        for (int i=0; i<numSecondsIncrements; i++)
            timeBuckets[i] = minRelativeSeconds + (i * secondsIncrement);

        try
        {
            for (int i=0; i<ms1FeatureFiles.length; i++)
            {
                File ms1File = ms1FeatureFiles[i];
                ApplicationContext.infoMessage("Processing MS1 file " + ms1File.getName());
                
                FeatureSet ms1FeatureSet = new FeatureSet(ms1File);

                Feature[] ms1FeaturesMassOrdered = ms1FeatureSet.getFeatures();
                Arrays.sort(ms1FeaturesMassOrdered, new Feature.MassAscComparator());

                ApplicationContext.infoMessage("Number of features to interrogate: " + ms1FeaturesMassOrdered.length);

                float[][] heatMapData = new float[numSecondsIncrements][numMassIncrements];

                int minMassIndex = 0;
                for (int j=0; j< ms1FeaturesMassOrdered.length; j++)
                {
                    Feature baseFeature = ms1FeaturesMassOrdered[j];
                    if (j % (ms1FeaturesMassOrdered.length / 20) == 0)
                        ApplicationContext.infoMessage((j * 100 / ms1FeaturesMassOrdered.length) + "% complete");
                    double minRelativeFeatureMass = baseFeature.getMass() + minRelativeMass -
                            MassUtilities.calculateAbsoluteDeltaMass(baseFeature.getMass(), matcherMassTolerance, matcherMassToleranceType);
                    double maxRelativeFeatureMass = baseFeature.getMass() + maxRelativeMass +
                            MassUtilities.calculateAbsoluteDeltaMass(baseFeature.getMass(), matcherMassTolerance, matcherMassToleranceType);

                    for (int k=minMassIndex; k<ms1FeaturesMassOrdered.length; k++)
                    {
                        //skip identity comparison
                        if (j==k)
                            continue;

                        Feature relativeFeature = ms1FeaturesMassOrdered[k];
                        if (relativeFeature.getMass() < minRelativeFeatureMass)
                        {
                            minMassIndex++;
                            continue;
                        }
                        else if (relativeFeature.getMass() > maxRelativeFeatureMass)
                            break;
                        float massDifference = relativeFeature.getMass() - baseFeature.getMass();
                        //if same mass, forget it
                        if (Math.abs(massDifference) < 0.5)
                            continue;
                        int massBucketIndex = (int) (massDifference / massWavelength) - minRelativeDaltons;
                        float bucketCenter = (float) (baseFeature.getMass() + minRelativeDaltons + (massBucketIndex * massWavelength));
                        //too far out of center of bucket
                        if (Math.abs(relativeFeature.getMass() - bucketCenter) >
                            MassUtilities.calculateAbsoluteDeltaMass(bucketCenter, matcherMassTolerance, matcherMassToleranceType))
                        {
//System.err.println("In mass range, mass=" + relativeFeature.getMass() + ", bucket=" + bucketCenter + ", dist to bucket center is " + Math.abs(relativeFeature.getMass() - bucketCenter) + ", more than " + AmtUtilities.calculateAbsoluteDeltaMass(bucketCenter, matcherMassTolerance, matcherMassToleranceType));
                                continue;
                        }

                        //in mass range.  Check seconds range
                        float timeDifference = relativeFeature.getTime() - baseFeature.getTime();
                        if (timeDifference < minRelativeSeconds || timeDifference > maxRelativeSeconds)
                            continue;
                        int secondsBucketIndex = (int) ((timeDifference - minRelativeSeconds) / secondsIncrement);
//System.err.println(massDifference + ", " + timeDifference + "... " + massBucketIndex + ", " + secondsBucketIndex);
                        //in range.  Add to correct bucket
//TODO: this sucks, stop doing it                        
if (secondsBucketIndex > numSecondsIncrements-1 || massBucketIndex > numMassIncrements-1) continue;
                        heatMapData[secondsBucketIndex][massBucketIndex]++;
                    }
                }
                PanelWithHeatMap pwhm = new PanelWithHeatMap(timeBuckets, massBuckets, heatMapData, "features at relative masses, times");
                pwhm.setAxisLabels("Relative time","Relative Mass (Da)");
                pwhm.displayInTab();
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

    }



}

