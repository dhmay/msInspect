/* 
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.commandline.arguments.DeltaMassArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.MassCalibrationUtilities;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.viewer.amt.AmtFeatureSetMatcher;
import org.fhcrc.cpl.viewer.amt.BaseAmtFeatureSetMatcherImpl;
import org.fhcrc.cpl.viewer.amt.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.plot.PlotOrientation;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.io.File;


/**
 * Command linemodule for adduct detection. Quick and dirty.  Lots of hacks.
 */
public class DetectAdductsCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(DetectAdductsCLM.class);

    protected File[] ms1FeatureFiles;

    protected File baseFeatureSetDir;

    protected double massWavelength =
            MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH;
    protected int scanWindowSize = 1;
    protected double minRelativeMass = 0;
    protected double maxRelativeMass = 100;

    protected float matcherMassTolerance = 10;
    protected int matcherMassToleranceType = AmtFeatureSetMatcher.DELTA_MASS_TYPE_PPM;
    protected int matcherScanTolerance=100;

    protected boolean onlyBaseFeaturesWithSTY=false;

    protected int currentMinMassIndexHACK = 0;

    protected File outZeroBucketFeatureFileDir = null;


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
                    createDirectoryToReadArgumentDefinition(
                            "basefeaturesdir",true,
                            "Directory of features files of known features"),
                    createUnnamedSeriesArgumentDefinition(ArgumentDefinitionFactory.FILE_TO_READ, true,
                            "Feature File of features to interrogate"),
                    createDecimalArgumentDefinition("masswavelength", false,
                            "Size of the mass bucket, in Daltons", massWavelength),
                    createIntegerArgumentDefinition("scanwindowsize", false,
                            "Size of the scan window (including identity scan)", scanWindowSize),
                    createIntegerArgumentDefinition("minrelativemass", false,
                            "Minimum relative mass", (int) minRelativeMass),
                    createIntegerArgumentDefinition("maxrelativemass", false,
                            "Maximum relative mass", (int) minRelativeMass),
                    createBooleanArgumentDefinition("onlybasefeatureswithsty", false,
                            "Only use base features with S, T, or Y in peptide ID", onlyBaseFeaturesWithSTY),
                    createDeltaMassArgumentDefinition("basefeaturematchdeltamass", false,
                            "Mass tolerance for matching features from the base set, to filter out identity features",
                            new DeltaMassArgumentDefinition.DeltaMassWithType(matcherMassTolerance, matcherMassToleranceType)),
                    createIntegerArgumentDefinition("basefeaturematchdeltascan", false,
                            "Scan tolerance for matching features from the base set, to filter out identity features", matcherScanTolerance),
                    createDirectoryToReadArgumentDefinition("outzerobucketfeaturesdir", false,
                            "Directory for outputting feature files containing all features that occur in the bucket that contains 0")
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

        minRelativeMass = (double) getIntegerArgumentValue("minrelativemass");

        minRelativeMass =  minRelativeMass -
                (minRelativeMass% massWavelength) -
                (massWavelength /2);
        maxRelativeMass = (double) getIntegerArgumentValue("maxrelativemass");
        maxRelativeMass = maxRelativeMass -
                ((maxRelativeMass% massWavelength) / 2);

        baseFeatureSetDir =
                getFileArgumentValue("basefeaturesdir");
        onlyBaseFeaturesWithSTY = getBooleanArgumentValue("onlybasefeatureswithsty");

        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassAndType =
                getDeltaMassArgumentValue("basefeaturematchdeltamass");
        matcherMassTolerance = deltaMassAndType.getDeltaMass();
        matcherMassToleranceType = deltaMassAndType.getDeltaMassType();
        matcherScanTolerance = getIntegerArgumentValue("basefeaturematchdeltascan");

        outZeroBucketFeatureFileDir = getFileArgumentValue("outzerobucketfeaturesdir");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        XYSeries[] allSeries = new XYSeries[ms1FeatureFiles.length];
        try
        {
            for (int i=0; i<ms1FeatureFiles.length; i++)
            {

                File ms1File = ms1FeatureFiles[i];
                ApplicationContext.infoMessage("Processing MS1 file " + ms1File.getName());
                
                FeatureSet ms1FeatureSet = new FeatureSet(ms1File);

                FeatureSet baseFeatureSet = null;

                String pepXmlFilename =
                        (ms1File.getName().substring(0,
                                ms1File.getName().indexOf(".")) + ".pep.xml");
                String tsvFilename =
                        (ms1File.getName().substring(0,
                                ms1File.getName().indexOf(".")) + ".tsv");
                String filteredTsvFilename =
                        (ms1File.getName().substring(0,
                                ms1File.getName().indexOf(".")) + ".filtered.tsv");
                boolean foundIt = false;
                for (String potentialMs2Filename : baseFeatureSetDir.list())
                {
                    if (potentialMs2Filename.equalsIgnoreCase(pepXmlFilename) ||
                        potentialMs2Filename.equalsIgnoreCase(tsvFilename) ||
                         potentialMs2Filename.equalsIgnoreCase(filteredTsvFilename)   )
                    {
                        File baseFeatureFile = new File(baseFeatureSetDir.getAbsolutePath() +
                                File.separatorChar + potentialMs2Filename);
                        baseFeatureSet = new FeatureSet(baseFeatureFile);
                        ApplicationContext.setMessage("Located MS2 file " + baseFeatureFile.getName() +
                                " with " + baseFeatureSet.getFeatures().length + " features");
                        foundIt=true;
                        break;
                    }
                 }
                if (!foundIt)
                    throw new CommandLineModuleExecutionException(
                            "Failed to find embedded MS2 or mzxml file for ms1 file " + ms1File.getName());

                allSeries[i] = processFeatureSet(ms1FeatureSet, baseFeatureSet);
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        if (allSeries.length == 1)
        {
            DefaultCategoryDataset dataset = new DefaultCategoryDataset();

            for (XYDataItem dataItem : (List<XYDataItem>) allSeries[0].getItems())
            {
                dataset.addValue(dataItem.getY(), "", dataItem.getX().doubleValue());
            }


            JFreeChart chart = ChartFactory.createBarChart(null, null, null, dataset,
                    PlotOrientation.VERTICAL, false, true, false);
            ChartDialog chartDialog = new ChartDialog(chart);
            chart.removeLegend();

            chartDialog.setSize(800,600);
            chartDialog.setVisible(true);
     

        }
        else
        {
            PanelWithLineChart panelWithLineChart = new PanelWithLineChart();
            for (XYSeries series : allSeries)
                panelWithLineChart.addSeries(series);
            ChartDialog chartDialog = new ChartDialog(panelWithLineChart);
            
            chartDialog.setSize(800,600);
            chartDialog.setVisible(true);
        }
    }

    protected void stripNonSTYFromFeatureSet(FeatureSet featureSet)
    {
        List<Feature> keptFeatureList = new ArrayList<Feature>();
        for (Feature feature : featureSet.getFeatures())
        {
            String peptideSequence = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptideSequence.contains("S") ||
                peptideSequence.contains("T") ||
                peptideSequence.contains("Y"))
                keptFeatureList.add(feature);
        }

        _log.debug("stripNonSTY: Removed " + keptFeatureList.size() +
                   " non-STY-bearing features out of " + featureSet.getFeatures().length +
                   " (" + (keptFeatureList.size() * 100 / featureSet.getFeatures().length) + "%)");

        featureSet.setFeatures(keptFeatureList.toArray(new Feature[keptFeatureList.size()]));
    }

    protected XYSeries processFeatureSet(FeatureSet ms1FeatureSet, FeatureSet baseFeatureSet)
    {
        _log.debug("Base features: " + baseFeatureSet.getFeatures().length);

        removeBaseSetMatches(ms1FeatureSet, baseFeatureSet);

        MS2ExtraInfoDef.removeAllButFirstFeatureForEachPeptide(baseFeatureSet);
        if (onlyBaseFeaturesWithSTY)
            stripNonSTYFromFeatureSet(baseFeatureSet);


        Feature[] ms1FeaturesMassOrdered = ms1FeatureSet.getFeatures();
        Arrays.sort(ms1FeaturesMassOrdered, new Feature.MassAscComparator());

        Feature[] baseFeaturesMassOrdered =  baseFeatureSet.getFeatures();
        Arrays.sort(baseFeaturesMassOrdered, new Feature.MassAscComparator());


        double maxFeatureMass = ms1FeaturesMassOrdered[ms1FeaturesMassOrdered.length - 1].getMass();
//        for (Feature feature : ms1FeaturesMassOrdered)
//            if (feature.getMass() > maxFeatureMass)
//                maxFeatureMass = feature.getMass();
//System.err.println("Min relative mass: " + minRelativeMass + ", max relative mass: " + maxRelativeMass);
        double rangeSize = maxRelativeMass - minRelativeMass;
        int numBuckets = (int)(rangeSize/ massWavelength);
//System.err.println("Number of buckets: " + numBuckets);
        ApplicationContext.infoMessage("Number of features to interrogate: " + ms1FeaturesMassOrdered.length);
        int fullDataIndex = 0;

        double[] massBuckets = new double[numBuckets];

        int localCurrentMinMassIndex = 0;

        List<Feature> zeroBucketFeatures = new ArrayList<Feature>();
        boolean shouldFillZeroBucket = (outZeroBucketFeatureFileDir != null);
        for (int i=0; i< baseFeaturesMassOrdered.length; i++)
        {
            if (i % (baseFeaturesMassOrdered.length / 20) == 0)
                ApplicationContext.infoMessage((i * 100 / baseFeaturesMassOrdered.length) + "% complete");
            Feature currentFeature = baseFeaturesMassOrdered[i];
            int minBucketScan = currentFeature.getScanFirst() - (scanWindowSize / 2);
            int maxBucketScan = currentFeature.getScanLast() + (scanWindowSize / 2);
            if (scanWindowSize % 2 == 0)
                maxBucketScan--;
//System.err.println("" + (i+1) + "               , size=" + scanWindowSize + ", min=" + minScan + ", max=" + maxScan);
//System.err.println("***Feature mass: " + currentFeature.getMass());

            for (int j = 0; j < massBuckets.length; j++)
            {
                double windowMinMass = currentFeature.getMass() + minRelativeMass +
                                       (j * massWavelength);
                if (windowMinMass > maxFeatureMass)
                    break;
                double windowMaxMass = windowMinMass + massWavelength;
//System.err.println("Feature mass: " + currentFeature.getMass() + ", min bucket: " + windowMinMass + ", max bucket: " + windowMaxMass);

                double relativeMinMass = windowMinMass - currentFeature.getMass();
                double relativeMaxMass = windowMaxMass - currentFeature.getMass();

                if (shouldFillZeroBucket && relativeMinMass < 0 && relativeMaxMass > 0)
                {
                    List<Feature> zeroBucketFeaturesThisFeature = findFeaturesInBox(windowMinMass,
                            windowMaxMass, minBucketScan, maxBucketScan, ms1FeatureSet.getFeatures(),
                            localCurrentMinMassIndex);
//System.err.println("Added " + zeroBucketFeaturesThisFeature.size() + " zero-bucket features");
                    zeroBucketFeatures.addAll(zeroBucketFeaturesThisFeature);
                    massBuckets[j] += zeroBucketFeaturesThisFeature.size();
                }
                else
                    massBuckets[j] += findNumFeaturesInBox(windowMinMass,
                            windowMaxMass, minBucketScan, maxBucketScan, ms1FeatureSet.getFeatures(),
                            localCurrentMinMassIndex);
                if (j == 0)
                    localCurrentMinMassIndex =  currentMinMassIndexHACK;                    
            }
            
        }
        XYSeries series = new XYSeries(ms1FeatureSet.getSourceFile().getName());
        for (int i=0; i<massBuckets.length; i++)
        {
            series.add((i * massWavelength) + minRelativeMass + (massWavelength / 2), massBuckets[i]);
//System.err.println("Bucket center: " + (i * massWavelength) + minRelativeMass + (massWavelength / 2) + ", min: " + (i * massWavelength) + minRelativeMass + ", max: " + (i * massWavelength) + minRelativeMass + (massWavelength) + ", center mass cluster deviation: " + (((i * massWavelength) + minRelativeMass + (massWavelength / 2)) % 1.000476));
        }

        if (shouldFillZeroBucket)
        {
            FeatureSet zeroBucketFeatureSet =
                    new FeatureSet(zeroBucketFeatures.toArray(new Feature[zeroBucketFeatures.size()]));
            try
            {
                File outZeroBucketFile = new File(outZeroBucketFeatureFileDir, ms1FeatureSet.getSourceFile().getName());

                zeroBucketFeatureSet.save(outZeroBucketFile);
                ApplicationContext.infoMessage("Saved zero-bucket file " + outZeroBucketFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Error saving zero-bucket features",e);
            }
        }

        return series;
    }

    protected int findNumFeaturesInBox(double minMass, double maxMass,
                                          double minScan, double maxScan,
                                          Feature[] ms1FeaturesMassOrdered,
                                          int localCurrentMinMassIndex)
    {
        int numFeaturesInBox = 0;
        int index = localCurrentMinMassIndex;
        while (ms1FeaturesMassOrdered[index].getMass() < minMass)
            index++;
        currentMinMassIndexHACK = index;
        while (index < ms1FeaturesMassOrdered.length &&
                ms1FeaturesMassOrdered[index].getMass() <= maxMass)
        {
            if (ms1FeaturesMassOrdered[index].getScanLast() >= minScan &&
                    ms1FeaturesMassOrdered[index].getScanFirst() <= maxScan)
                numFeaturesInBox++;
            index++;
        }
        return numFeaturesInBox;
    }

    protected List<Feature> findFeaturesInBox(double minMass, double maxMass,
                                          double minScan, double maxScan,
                                          Feature[] ms1FeaturesMassOrdered,
                                          int localCurrentMinMassIndex)
    {
        List<Feature> resultList = new ArrayList<Feature>();
        int index = localCurrentMinMassIndex;
        while (ms1FeaturesMassOrdered[index].getMass() < minMass)
            index++;
        currentMinMassIndexHACK = index;
        while (index < ms1FeaturesMassOrdered.length &&
                ms1FeaturesMassOrdered[index].getMass() <= maxMass)
        {
            if (ms1FeaturesMassOrdered[index].getScanLast() >= minScan &&
                    ms1FeaturesMassOrdered[index].getScanFirst() <= maxScan)
                resultList.add(ms1FeaturesMassOrdered[index]);
            index++;
        }
        return resultList;
    }

    protected void removeBaseSetMatches(FeatureSet ms1FeatureSet, FeatureSet baseFeatureSet)
    {
System.err.println("Interrogation features before removing base set matches: " + ms1FeatureSet.getFeatures().length);
        ClusteringFeatureSetMatcher fsm =
                new ClusteringFeatureSetMatcher(matcherMassTolerance, matcherMassToleranceType, matcherScanTolerance);
        fsm.setElutionMode(BaseAmtFeatureSetMatcherImpl.ELUTION_MODE_SCAN);
        AmtFeatureSetMatcher.FeatureMatchingResult matchingResult =
                fsm.matchFeatures(baseFeatureSet, ms1FeatureSet);
        Set<Feature> ms1MatchingFeatures =
                matchingResult.getSlaveSetFeatures();

        List<Feature> ms1UnmatchedFeatures = new ArrayList<Feature>();
        for (Feature feature : ms1FeatureSet.getFeatures())
        {
            if (!ms1MatchingFeatures.contains(feature))
                ms1UnmatchedFeatures.add(feature);
        }
        

        ms1FeatureSet.setFeatures(ms1UnmatchedFeatures.toArray(new Feature[ms1UnmatchedFeatures.size()]));
System.err.println("Interrogation features after removing base set matches: " + ms1FeatureSet.getFeatures().length);

    }

}

