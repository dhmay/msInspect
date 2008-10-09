/* 
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.gui.util.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.awt.*;


/**
 * test
 */
public class PeptideQuantVisualizationCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptideQuantVisualizationCLM.class);


    protected int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;


    int numHeavyPeaksToPlot = 4;

    protected File outDir = null;

    protected int scanImageHeight = 100;
    protected int maxScansImageHeight = 4000;
    protected int spectrumImageHeight = 1000;
    protected int imageWidth = 1000;

    protected File pepXmlFile;
    protected File mzXmlDir;

    protected Set<String> peptidesToExamine;
    protected Set<String> proteinsToExamine;


    protected float minPeptideProphet = 0;

    protected File mzXmlFile = null;

    protected int numPaddingScans = 2;
    protected float mzPadding = 1;

    Set<String> peptidesFound = new HashSet<String>();

    public PeptideQuantVisualizationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "quantcharts";

        mHelpMessage ="quantcharts";
        mShortDescription = "PeptideSpectrumChart";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createIntegerArgumentDefinition("spectrumimageheight", false, "Spectrum image height",
                                spectrumImageHeight),
                        createIntegerArgumentDefinition("resolution", false, "resolution (number of breaks per Thompson",
                                resolution),
                        createIntegerArgumentDefinition("width", false, "Image width (used for both)", imageWidth),
                        createIntegerArgumentDefinition("scanimageheight", false,
                                "Height of EACH per-scan image, in the output file",
                                scanImageHeight),
                        createDirectoryToReadArgumentDefinition("outdir", true, "Output directory"),
                        createIntegerArgumentDefinition("maxscansimageheight", false,
                                "Maximum overall height for the all-scans line plot image (overrides scansfileimageheight)",
                                maxScansImageHeight),
                        createFileToReadArgumentDefinition("pepxml", true, "pepXML file"),
                        createFileToReadArgumentDefinition("mzxml", false, "mzXML file"),

                        createDirectoryToReadArgumentDefinition("mzxmldir", false, "mzXML directory"),
                        createStringArgumentDefinition("peptides", false, "comma-separated list of peptides to examine"),
                        createStringArgumentDefinition("proteins", false, "comma-separated list of proteins to examine"),
                        createDecimalArgumentDefinition("minpprophet", false, "minimum PeptideProphet value",
                                minPeptideProphet),
                        createIntegerArgumentDefinition("paddingscans", false,
                                "number of scans before and after quant envelope to display", numPaddingScans),
                        createDecimalArgumentDefinition("mzpadding", false,
                                "amount of m/z space to display around quant", mzPadding),
                        createIntegerArgumentDefinition("numpeaksaboveheavy", false,
                                "number of peaks above the heavy-ion monoisotope to display", numHeavyPeaksToPlot)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        resolution = getIntegerArgumentValue("resolution");


        spectrumImageHeight = getIntegerArgumentValue("spectrumimageheight");
        imageWidth = getIntegerArgumentValue("width");
        scanImageHeight = getIntegerArgumentValue("scanimageheight");
        maxScansImageHeight = getIntegerArgumentValue("maxscansimageheight");

        pepXmlFile = getFileArgumentValue("pepxml");

        mzXmlFile = getFileArgumentValue("mzxml");
        mzXmlDir = getFileArgumentValue("mzxmldir");

        if (mzXmlFile == null && mzXmlDir == null)
            throw new ArgumentValidationException("Must supply mzxml or mzxmldir argument");

        outDir = getFileArgumentValue("outdir");

        minPeptideProphet = getFloatArgumentValue("minpprophet");

        if ((hasArgumentValue("peptides") && hasArgumentValue("proteins")) ||
                (!hasArgumentValue("peptides") && !hasArgumentValue("proteins")))
            throw new ArgumentValidationException("Please supply either peptides or proteins (but not both");

        if (hasArgumentValue("peptides"))
        {
            String peptidesString = getStringArgumentValue("peptides");
            String[] peptidesArray = peptidesString.split(",");
            peptidesToExamine = new HashSet<String>();
            for (String peptide : peptidesArray)
                peptidesToExamine.add(peptide);
        }

        if (hasArgumentValue("proteins"))
        {
            String proteinsString = getStringArgumentValue("proteins");
            String[] proteinsArray = proteinsString.split(",");
            proteinsToExamine = new HashSet<String>();
            for (String protein : proteinsArray)
                proteinsToExamine.add(protein);
        }

        numPaddingScans = getIntegerArgumentValue("paddingscans");
        mzPadding = getFloatArgumentValue("mzpadding");
        numHeavyPeaksToPlot = getIntegerArgumentValue("numpeaksaboveheavy");
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Iterator<FeatureSet> fsi = null;
        try
        {
            fsi = new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
        }
        catch (IOException e)
        {
            try
            {
                FeatureSet featureSet = new FeatureSet(pepXmlFile);
                List<FeatureSet> dummyList = new ArrayList<FeatureSet>();
                dummyList.add(featureSet);
                fsi = dummyList.iterator();
            }
            catch (Exception e2)
            {
                throw new CommandLineModuleExecutionException("Failed to open file ",e2);
            }
        }

        while (fsi.hasNext())
        {
            FeatureSet fraction = fsi.next();
            ApplicationContext.infoMessage("Handling fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(fraction));
            handleFraction(fraction);
        }
    }

    protected void handleFraction(FeatureSet featureSet)
            throws CommandLineModuleExecutionException
    {
        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet(minPeptideProphet);
        featureSet = featureSet.filter(sel);

        Arrays.sort(featureSet.getFeatures(), new Feature.ScanAscComparator());

        MSRun run = null;
        File fileToLoad = mzXmlFile;
        if (fileToLoad == null)
        {
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            try
            {
                fileToLoad = CommandLineModuleUtilities.findFileWithPrefix(fractionName, mzXmlDir, "mzXML");
                ApplicationContext.infoMessage("Found mzXML file " + fileToLoad.getAbsolutePath());
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Couldn't find mzXML file for fraction " +
                        fractionName);
            }
        }

        try
        {
            run = MSRun.load(fileToLoad.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to load run",e);
        }

        if (proteinsToExamine != null)
            handleProteinsInRun(featureSet, run, proteinsToExamine);
        else
            handlePeptidesInRun(featureSet, run, peptidesToExamine, null, outDir);

        if (peptidesFound.isEmpty())
            ApplicationContext.infoMessage("No peptides found with passing quantitation in the files specified!");
    }

    protected void handleProteinsInRun(FeatureSet featureSet, MSRun run, Set<String> proteinsToHandle)
            throws CommandLineModuleExecutionException
    {
        for (String protein : proteinsToHandle)
        {
            Set<String> peptidesThisProtein = new HashSet<String>();
            for (Feature feature : featureSet.getFeatures())
            {
                List<String> proteinsThisFeature = MS2ExtraInfoDef.getProteinList(feature);
                if (proteinsThisFeature != null && proteinsThisFeature.contains(protein))
                    peptidesThisProtein.add(MS2ExtraInfoDef.getFirstPeptide(feature));
            }

            if (!peptidesThisProtein.isEmpty())
            {
                ApplicationContext.infoMessage("Protein " + protein + ": finding " + peptidesThisProtein.size() +
                        " peptides in run " + MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                File proteinOutDir = new File(outDir, protein);
                proteinOutDir.mkdir();
                handlePeptidesInRun(featureSet, run, peptidesThisProtein, protein, proteinOutDir);
            }
        }
    }

    protected void handlePeptidesInRun(FeatureSet featureSet, MSRun run, Set<String> peptidesToHandle,
                                       String proteinName, File outputDir)
            throws CommandLineModuleExecutionException
    {
        String fraction = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);

        Map<String, List<Feature>> peptideFeaturesMap = new HashMap<String, List<Feature>>();
        for (Feature feature : featureSet.getFeatures())
        {
            String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (featurePeptide == null)
                continue;
            if (featurePeptide != null && peptidesToHandle.contains(featurePeptide) &&
                    IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
                if (!peptidesFound.contains(featurePeptide))
                    ApplicationContext.infoMessage("\tFound peptide " + featurePeptide);
                peptidesFound.add(featurePeptide);
                List<Feature> featuresThisPeptide = peptideFeaturesMap.get(featurePeptide);
                if (featuresThisPeptide == null)
                {
                    featuresThisPeptide = new ArrayList<Feature>();
                    peptideFeaturesMap.put(featurePeptide, featuresThisPeptide);
                }
                featuresThisPeptide.add(feature);
            }
        }

        for (String peptide : peptideFeaturesMap.keySet())
        {
            int numTotalEventsThisPeptide = peptideFeaturesMap.get(peptide).size();
            ApplicationContext.infoMessage("\tHandling peptide " + peptide + " with " + numTotalEventsThisPeptide + " events");           
            List<Feature> nonOverlappingFeatures = findNonOverlappingEvents(peptideFeaturesMap.get(peptide));

            for (Feature feature : nonOverlappingFeatures)
            {
                int firstQuantScan = IsotopicLabelExtraInfoDef.getLightFirstScan(feature);
                int lastQuantScan = IsotopicLabelExtraInfoDef.getLightLastScan(feature);

                if (firstQuantScan <= 0)
                    firstQuantScan = feature.getScan();
                if (lastQuantScan <= 0)
                    lastQuantScan = feature.getScan();

                int minScanIndex = Math.max(run.getIndexForScanNum(firstQuantScan) - numPaddingScans, 0);
                int maxScanIndex = Math.min(run.getIndexForScanNum(lastQuantScan) + numPaddingScans,
                        run.getScanCount()-1);

                int minScan = run.getScanNumForIndex(minScanIndex);
                int maxScan = run.getScanNumForIndex(maxScanIndex);
                int minMz = (int) (IsotopicLabelExtraInfoDef.getLightMass(feature) / feature.getCharge() - mzPadding);
                int maxMz = (int) (IsotopicLabelExtraInfoDef.getHeavyMass(feature) / feature.getCharge() +
                        numHeavyPeaksToPlot / feature.getCharge() + mzPadding);
                createChartsForEvent(run, feature.getScan(), minScan, maxScan, minMz, maxMz, firstQuantScan, lastQuantScan,
                        outputDir, peptide, fraction, proteinName);
            }
        }
    }

    protected List<Feature> findNonOverlappingEvents(List<Feature> features)
    {
        List<Feature> result = new ArrayList<Feature>();
        Map<Integer, List<Feature>> chargeFeatureMap = new HashMap<Integer, List<Feature>>();
        for (Feature feature : features)
        {
            List<Feature> featuresThisCharge = chargeFeatureMap.get(feature.getCharge());
            if (featuresThisCharge == null)
            {
                featuresThisCharge = new ArrayList<Feature>();
                chargeFeatureMap.put(feature.getCharge(), featuresThisCharge);
            }
            featuresThisCharge.add(feature);
        }
        for (List<Feature> featuresWithinCharge : chargeFeatureMap.values())
        {
            while (!featuresWithinCharge.isEmpty())
                result.add(findFeatureRepresentingFirstFeatureOverlap(featuresWithinCharge));
        }
        return result;
    }

    /**
     * Find all features overlapping in scans with the first feature.  Choose on of them to represent them all
     * ("median" by scan).
     * SIDE EFFECT: Remove them all from featuresWithinCharge
     * @param featuresWithinCharge
     * @return
     */
    protected Feature findFeatureRepresentingFirstFeatureOverlap(List<Feature> featuresWithinCharge)
    {
        List<Feature> featuresOverlappingFirst = new ArrayList<Feature>();
        featuresOverlappingFirst.add(featuresWithinCharge.get(0));

        Feature result = featuresOverlappingFirst.get(0);

        if (featuresWithinCharge.size() > 1)
        {
            int firstFeatureStart = IsotopicLabelExtraInfoDef.getLightFirstScan(result);
            int firstFeatureEnd = IsotopicLabelExtraInfoDef.getLightLastScan(result);
            if (firstFeatureStart <= 0)
                firstFeatureStart = result.getScan();
            if (firstFeatureEnd <= 0)
                firstFeatureEnd = result.getScan();

            for (int i=1; i<featuresWithinCharge.size(); i++)
            {
                Feature compareFeature = featuresWithinCharge.get(i);
                int compareFeatureStart = IsotopicLabelExtraInfoDef.getLightFirstScan(compareFeature);
                int compareFeatureEnd = IsotopicLabelExtraInfoDef.getLightLastScan(compareFeature);
                if (compareFeatureStart <= 0)
                    compareFeatureStart = compareFeature.getScan();
                if (compareFeatureEnd <= 0)
                    compareFeatureEnd = compareFeature.getScan();
                if (Math.max(firstFeatureStart, compareFeatureStart) < Math.min(compareFeatureEnd, compareFeatureEnd))
                    featuresOverlappingFirst.add(compareFeature);
            }


            //find the "median" feature by scan (round down if even number)
            Collections.sort(featuresOverlappingFirst, new Feature.ScanAscComparator());
            int numOverlapping = featuresOverlappingFirst.size();
            if (numOverlapping % 2 == 1)
                result = featuresOverlappingFirst.get((numOverlapping-1)/2);
            else
                result = featuresOverlappingFirst.get((numOverlapping-2)/2);
            ApplicationContext.infoMessage("\t\tCollapsing " + featuresOverlappingFirst.size() +
                    " features into one event, charge " + result.getCharge());

        }

        featuresWithinCharge.removeAll(featuresOverlappingFirst);

        return result;
    }

    protected void createChartsForEvent(MSRun run, int scan, int minScan, int maxScan, int minMz, int maxMz,
                                        int scanLine1, int scanLine2, File outputDir, String peptide, String fraction,
                                        String proteinName)
            throws CommandLineModuleExecutionException
    {
        PanelWithSpectrumChart spectrumPanel =
                new PanelWithSpectrumChart(run, minScan, maxScan, minMz, maxMz, resolution, scanLine1, scanLine2,
                        true);
        spectrumPanel.setName("Spectrum");
        spectrumPanel.setVisible(true);

        DropdownMultiChartDisplayPanel multiChartPanelForScans = new DropdownMultiChartDisplayPanel();
        multiChartPanelForScans.setDisplaySlideshowButton(true);

        spectrumPanel.setSize(new Dimension(imageWidth, spectrumImageHeight));
        spectrumPanel.setMinimumSize(new Dimension(imageWidth, spectrumImageHeight));

        PanelWithChart spectrumChartToDisplay = spectrumPanel;

        Map<Integer, PanelWithLineChart> scanChartMap = spectrumPanel.getScanLineChartMap();
        List<Integer> allScans = new ArrayList<Integer>(scanChartMap.keySet());
        Collections.sort(allScans);

        for (int scanForChart : allScans)
            multiChartPanelForScans.addChartPanel(scanChartMap.get(scanForChart));

        String filePrefix = peptide + "_" + fraction + "_" + scan;
        if (proteinName != null)
            filePrefix = proteinName + "_" + filePrefix;

        File outSpectrumFile = new File(outputDir, filePrefix + "_spectrum.png");
        File outScansFile = new File(outputDir, filePrefix + "_scans.png");


        try
        {
            spectrumChartToDisplay.saveChartToImageFile(outSpectrumFile);
            ApplicationContext.infoMessage("Wrote spectrum to image " + outSpectrumFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to save image file",e);
        }

        try
        {
            spectrumPanel.savePerScanSpectraImage(imageWidth, scanImageHeight,
                    maxScansImageHeight, outScansFile);
            ApplicationContext.infoMessage("Wrote scans to image " + outScansFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to write image file " +
                    outScansFile.getAbsolutePath(),e);
        }

    }

}
