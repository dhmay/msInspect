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
package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.gui.util.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;


/**
 * Visualize quantitation events.  This class figures out what events to visualize and gathers all
 * the data from each event, but the actual creation of charts is done in PanelWithSpectrumChart.
 *
 * Summaries of the quantitation event information are written in TSV and HTML format.  Qurate
 * can load the TSV version
 */
public class QuantitationVisualizer
{
    protected static Logger _log = Logger.getLogger(QuantitationVisualizer.class);


    protected int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;
    public static final int DEFAULT_IMAGE_HEIGHT_3D = 900;
    public static final int DEFAULT_IMAGE_WIDTH_3D = 900;

    public static final int DEFAULT_SINGLE_SCAN_IMAGE_HEIGHT = 100;
    public static final int DEFAULT_MAX_SINGLE_SCANS_TOTAL_IMAGE_HEIGHT = 4000;
    public static final int DEFAULT_SPECTRUM_IMAGE_HEIGHT = 700;
    public static final int DEFAULT_SPECTRUM_IMAGE_WIDTH = 900;


    protected Iterator<FeatureSet> featureSetIterator;


    int numHeavyPeaksToPlot = 4;

    protected File outDir = null;
    protected File outTsvFile = null;
    protected File outHtmlFile = null;

    protected float peakSeparationMass = PanelWithSpectrumChart.DEFAULT_PEAK_SEPARATION_MASS;
    protected float peakTolerancePPM = PanelWithSpectrumChart.DEFAULT_PEAK_TOLERANCE_PPM;

    protected int scanImageHeight = DEFAULT_SINGLE_SCAN_IMAGE_HEIGHT;
    protected int maxScansImageHeight = DEFAULT_MAX_SINGLE_SCANS_TOTAL_IMAGE_HEIGHT;
    protected int spectrumImageHeight = DEFAULT_SPECTRUM_IMAGE_HEIGHT;
    protected int imageWidth = DEFAULT_SPECTRUM_IMAGE_WIDTH;

    protected File mzXmlDir;

    //The amount by which two features can differ and still be combined into the same single representation.
    //TODO: convert this to log space, parameterize
    protected float maxCombineFeatureRatioDiff = 0.2f;

    protected boolean writeHTMLAndText = true;

    protected float minPeptideProphet = 0;

    protected File mzXmlFile = null;

    protected int numPaddingScans = 3;
    protected float mzPadding = 1;

    protected Set<String> peptidesFound = new HashSet<String>();
    protected int sidebarWidth = 180;

    protected Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>> proteinPeptideFractionChargeFilesMap =
            new HashMap<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>>();


    //For QuantEventInfo objects with no associated protein
    protected static final String DUMMY_PROTEIN_NAME = "DUMMY_PROTEIN";

    protected boolean show3DPlots = true;
    protected int rotationAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_ROTATION;
    protected int tiltAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_TILT;
    protected int imageHeight3D = DEFAULT_IMAGE_HEIGHT_3D;
    protected int imageWidth3D = DEFAULT_IMAGE_WIDTH_3D;
    protected boolean show3DAxes = true;
    //Controls whether we write event info directly to the chart images
    protected boolean writeInfoOnCharts = false;

    int scan = 0;


    protected Set<String> peptidesToExamine;
    protected Set<String> proteinsToExamine;
    protected Set<String> fractionsToExamine;


    protected PrintWriter outHtmlPW;
    protected PrintWriter outTsvPW;

    protected boolean showProteinColumn = true;


    public QuantitationVisualizer()
    {
    }

    /**
     * Visualize just the events specified
     * @param quantEvents
     * @throws IOException
     */
    public void visualizeQuantEvents(List<QuantEventInfo> quantEvents) throws IOException
    {
        if (outHtmlFile == null)
            outHtmlFile = new File(outDir,"quantitation.html");
        if (outTsvFile == null)
            outTsvFile = new File(outDir,"quantitation.tsv");
        if (writeHTMLAndText)
        {
                outHtmlPW = new PrintWriter(outHtmlFile);
                outTsvPW = new PrintWriter(outTsvFile);


            QuantEventInfo.writeHeader(outHtmlPW, outTsvPW, showProteinColumn, show3DPlots);
        }

        Map<String, List<QuantEventInfo>> fractionEventMap = new HashMap<String, List<QuantEventInfo>>();

        for (QuantEventInfo quantEvent : quantEvents)
        {
            List<QuantEventInfo> eventList = fractionEventMap.get(quantEvent.getFraction());
            if (eventList == null)
            {
                eventList = new ArrayList<QuantEventInfo>();
                fractionEventMap.put(quantEvent.getFraction(), eventList);
            }
            eventList.add(quantEvent);
        }
        for (String fraction : fractionEventMap.keySet())
        {
            File mzXmlFile = CommandLineModuleUtilities.findFileWithPrefix(fraction, mzXmlDir, "mzXML");
            MSRun run = MSRun.load(mzXmlFile.getAbsolutePath());

            for (QuantEventInfo quantEvent : fractionEventMap.get(fraction))
            {
                createChartsForEvent(run, outDir, quantEvent.getProtein(), fraction, quantEvent);                
            }
        }

        if (writeHTMLAndText)
        {
            QuantEventInfo.writeFooterAndClose(outHtmlPW, outTsvPW);
            ApplicationContext.infoMessage("Saved HTML file " + outHtmlFile.getAbsolutePath());
            ApplicationContext.infoMessage("Saved TSV file " + outTsvFile.getAbsolutePath());
        }        
    }

    /**
     * Iterate through all fractions, finding and visualizing the selected events
     */
    public void visualizeQuantEvents() throws IOException
    {
        //trying to prevent memory leaks
        ImageIO.setUseCache(false);
        
        if (outHtmlFile == null)
            outHtmlFile = new File(outDir,"quantitation.html");
        if (outTsvFile == null)
            outTsvFile = new File(outDir,"quantitation.tsv");


        if (writeHTMLAndText)
        {
                outHtmlPW = new PrintWriter(outHtmlFile);
                outTsvPW = new PrintWriter(outTsvFile);


            QuantEventInfo.writeHeader(outHtmlPW, outTsvPW, showProteinColumn, show3DPlots);
        }

        boolean processedAFraction = false;
        while (featureSetIterator.hasNext())
        {
            FeatureSet fraction = featureSetIterator.next();
            ApplicationContext.infoMessage("Handling fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(fraction));
            if (fractionsToExamine == null || fractionsToExamine.contains(MS2ExtraInfoDef.getFeatureSetBaseName(fraction)))
            {
                handleFraction(fraction);
                processedAFraction = true;
            }
            //trying to prevent memory overflow
            System.gc();
        }
        if (!processedAFraction)
            ApplicationContext.infoMessage("WARNING: no fractions processed");
        if (writeHTMLAndText)
        {
            QuantEventInfo.writeFooterAndClose(outHtmlPW, outTsvPW);
            ApplicationContext.infoMessage("Saved HTML file " + outHtmlFile.getAbsolutePath());
            ApplicationContext.infoMessage("Saved TSV file " + outTsvFile.getAbsolutePath());
        }

    }

    /**
     * Visualize selected events from a fraction
     * @param featureSet
     * @throws CommandLineModuleExecutionException
     */
    protected void handleFraction(FeatureSet featureSet) throws IOException
    {
        String errorMessage = "";

        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinPProphet(minPeptideProphet);
        featureSet = featureSet.filter(sel);

        Arrays.sort(featureSet.getFeatures(), new Feature.ScanAscComparator());

        MSRun run = null;
        File fileToLoad = mzXmlFile;
        if (fileToLoad == null)
        {
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            errorMessage = "Error finding mzXML files";
            fileToLoad = CommandLineModuleUtilities.findFileWithPrefix(fractionName, mzXmlDir, "mzXML");
            ApplicationContext.infoMessage("Found mzXML file " + fileToLoad.getAbsolutePath());

        }

        run = MSRun.load(fileToLoad.getAbsolutePath());


        if (proteinsToExamine != null)
            handleProteinsInRun(featureSet, run, proteinsToExamine);
        else if (peptidesToExamine != null)
            handlePeptidesInRun(featureSet, run, peptidesToExamine, DUMMY_PROTEIN_NAME, outDir);
        else
        {
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            if (fractionsToExamine == null || fractionsToExamine.contains(fractionName))
                handleScanInRun(featureSet, fractionName, run);
        }
    }

    /**
     * Visualize all charts for a particular scan number.  Awfully inefficient, but
     * it doesn't matter because this is only invoked for single-scan runs
     * @param featureSet
     * @param fraction
     * @param run
     * @throws CommandLineModuleExecutionException
     */
    protected void handleScanInRun(FeatureSet featureSet, String fraction, MSRun run)
    {
        for (Feature feature : featureSet.getFeatures())
        {
            if (scan == feature.getScan())
                createChartsForEvent(run, outDir, DUMMY_PROTEIN_NAME, fraction, feature);
        }
    }

    /**
     * Display all passing events for a particular set of proteins
     * @param featureSet
     * @param run
     * @param proteinsToHandle
     * @throws CommandLineModuleExecutionException
     */
    protected void handleProteinsInRun(FeatureSet featureSet, MSRun run, Set<String> proteinsToHandle)
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

    /**
     * Display all passing events for a particular set of peptides
     * @param featureSet
     * @param run
     * @param peptidesToHandle
     * @param proteinName
     * @param outputDir
     * @throws CommandLineModuleExecutionException
     */
    protected void handlePeptidesInRun(FeatureSet featureSet, MSRun run, Set<String> peptidesToHandle,
                                       String proteinName, File outputDir)
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
                createChartsForEvent(run, outputDir, proteinName, fraction, feature);
            }
        }
    }

    /**
     * Given a list of features (presumably with the same peptide ID, in the same charge and
     * modification state), collapse them down to a single feature representing each set of
     * overlapping events
     * @param features
     * @return
     */
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
     * Find all features overlapping in scans with the first feature.  Choose one of them to represent them all
     * ("median" by scan).  Store the others in Feature.comprised
     * SIDE EFFECT: Remove them all from featuresWithinCharge
     * @param featuresWithinCharge
     * @return
     */
    protected Feature findFeatureRepresentingFirstFeatureOverlap(List<Feature> featuresWithinCharge)
    {
        List<Feature> featuresOverlappingFirst = new ArrayList<Feature>();
        List<Spectrum.Peak> otherFeaturesAsPeaks = new ArrayList<Spectrum.Peak>();

        featuresOverlappingFirst.add(featuresWithinCharge.get(0));

        Feature firstFeature = featuresOverlappingFirst.get(0);
        Feature result = firstFeature;

        if (featuresWithinCharge.size() > 1)
        {
            int firstFeatureStart = IsotopicLabelExtraInfoDef.getLightFirstScan(firstFeature);
            int firstFeatureEnd = IsotopicLabelExtraInfoDef.getLightLastScan(firstFeature);
            if (firstFeatureStart <= 0)
                firstFeatureStart = firstFeature.getScan();
            if (firstFeatureEnd <= 0)
                firstFeatureEnd = firstFeature.getScan();
            double firstFeatureRatio = IsotopicLabelExtraInfoDef.getRatio(firstFeature);
            for (int i=1; i<featuresWithinCharge.size(); i++)
            {
                Feature compareFeature = featuresWithinCharge.get(i);
                int compareFeatureStart = IsotopicLabelExtraInfoDef.getLightFirstScan(compareFeature);
                int compareFeatureEnd = IsotopicLabelExtraInfoDef.getLightLastScan(compareFeature);
                if (compareFeatureStart <= 0)
                    compareFeatureStart = compareFeature.getScan();
                if (compareFeatureEnd <= 0)
                    compareFeatureEnd = compareFeature.getScan();
                if (Math.max(firstFeatureStart, compareFeatureStart) < Math.min(firstFeatureEnd, compareFeatureEnd) &&
                        Math.abs(IsotopicLabelExtraInfoDef.getRatio(compareFeature) - firstFeatureRatio) <
                                maxCombineFeatureRatioDiff)
                {
                    featuresOverlappingFirst.add(compareFeature);
                    otherFeaturesAsPeaks.add(compareFeature);
                }
            }


            //find the "median" feature by scan (round down if even number)
            Collections.sort(featuresOverlappingFirst, new Feature.ScanAscComparator());
            int numOverlapping = featuresOverlappingFirst.size();
            if (numOverlapping % 2 == 1)
                result = featuresOverlappingFirst.get((numOverlapping-1)/2);
            else
                result = featuresOverlappingFirst.get((numOverlapping-2)/2);
        }
        ApplicationContext.infoMessage("\t\tcharge " + result.getCharge() + ", ratio " +
                IsotopicLabelExtraInfoDef.getRatio(result) + ", scans " +
                IsotopicLabelExtraInfoDef.getLightFirstScan(result) + "-" +
                IsotopicLabelExtraInfoDef.getLightLastScan(result) + ", represents " +
                featuresOverlappingFirst.size() + " event(s)");


        featuresWithinCharge.removeAll(featuresOverlappingFirst);
        result.comprised = otherFeaturesAsPeaks.toArray(new Spectrum.Peak[otherFeaturesAsPeaks.size()]);
        return result;
    }

    /**
     * Call PanelWithSpectrumChart to create all the charts for a quantitative event. Save
     * all the charts.
     * Writes a row to the output files
     * @param run
     * @param outputDir
     * @param protein
     * @param fraction
     * @param feature
     * @throws CommandLineModuleExecutionException
     */
    protected void createChartsForEvent(MSRun run,File outputDir, String protein, String fraction,
                                        Feature feature)
    {
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        QuantEventInfo quantEvent =
                new QuantEventInfo(feature, fraction);
    }

    protected void createChartsForEvent(MSRun run,File outputDir, String protein, String fraction,
                                        QuantEventInfo quantEvent)
    {
        String filePrefix = quantEvent.getPeptide() + "_" + fraction + "_" + quantEvent.getCharge() + "_" + quantEvent.getScan();


        //chart output files
        if (quantEvent.getSpectrumFile() == null)
            quantEvent.setSpectrumFile(new File(outputDir, filePrefix + "_spectrum.png"));
        if (quantEvent.getScansFile() == null)
            quantEvent.setScansFile(new File(outputDir, filePrefix + "_scans.png"));
        if (quantEvent.getIntensitySumFile() == null)
            quantEvent.setIntensitySumFile(new File(outputDir, filePrefix + "_intensitysum.png"));
        if (quantEvent.getFile3D() == null)
            quantEvent.setFile3D(new File(outputDir, filePrefix + "_3D.png"));



        int firstLightQuantScan = quantEvent.getFirstLightQuantScan();
        int lastLightQuantScan = quantEvent.getLastLightQuantScan();
        int firstHeavyQuantScan = quantEvent.getFirstLightQuantScan();
        int lastHeavyQuantScan = quantEvent.getLastHeavyQuantScan();

        //Add padding around quant event limits
        int minScanIndex = Math.max(Math.abs(run.getIndexForScanNum(Math.min(firstLightQuantScan, firstHeavyQuantScan))) -
                numPaddingScans, 0);
        int maxScanIndex = Math.min(Math.abs(run.getIndexForScanNum(Math.min(lastLightQuantScan, lastHeavyQuantScan))) +
                numPaddingScans, run.getScanCount()-1);

        int minScan = run.getScanNumForIndex(minScanIndex);
        int maxScan = run.getScanNumForIndex(maxScanIndex);

        float minMz = quantEvent.getLightMz() - mzPadding;
        float maxMz = quantEvent.getHeavyMz() + numHeavyPeaksToPlot / quantEvent.getCharge() + mzPadding;
        _log.debug("Building chart for event:\n\t" + quantEvent);
        _log.debug("Scan=" + quantEvent.getScan() + ", ratio=" + quantEvent.getRatio() +
                ", lightInt=" + quantEvent.getLightIntensity() + ", heavyInt=" +
                quantEvent.getHeavyIntensity() + ", minScanIndex=" + minScanIndex +
                ", maxScanIndex=" + maxScanIndex + ", minMz=" +
                minMz + ", maxMz=" + maxMz);


        if (protein != null && !protein.equals(DUMMY_PROTEIN_NAME))
            filePrefix = protein + "_" + filePrefix;

        //create the PanelWithSpectrumChart and make it generate everything
        PanelWithSpectrumChart spectrumPanel =
                new PanelWithSpectrumChart(run, minScan, maxScan, minMz, maxMz, 
                        firstLightQuantScan, lastLightQuantScan, firstHeavyQuantScan, lastHeavyQuantScan,
                        quantEvent.getLightMz(), quantEvent.getHeavyMz(), quantEvent.getCharge());
        spectrumPanel.setResolution(resolution);
        spectrumPanel.setGenerateLineCharts(true);
        spectrumPanel.setGenerate3DChart(show3DPlots);
        spectrumPanel.setIdEventScan(quantEvent.getScan());
        spectrumPanel.setIdEventMz(quantEvent.getMz());
        spectrumPanel.setOtherEventScans(quantEvent.getOtherEventScans());
        spectrumPanel.setOtherEventMZs(quantEvent.getOtherEventMZs());
        spectrumPanel.setName("Spectrum");
        spectrumPanel.setContourPlotRotationAngle(rotationAngle3D);
        spectrumPanel.setContourPlotTiltAngle(tiltAngle3D);
        spectrumPanel.setContourPlotWidth(imageWidth3D);
        spectrumPanel.setContourPlotHeight(imageHeight3D);
        spectrumPanel.setContourPlotShowAxes(show3DAxes);
        spectrumPanel.setSize(new Dimension(imageWidth, spectrumImageHeight));

        spectrumPanel.setPeakSeparationMass(peakSeparationMass);
        spectrumPanel.setPeakTolerancePPM(peakTolerancePPM);

        spectrumPanel.generateCharts();
        spectrumPanel.setVisible(true);
        spectrumPanel.setMinimumSize(new Dimension(imageWidth, spectrumImageHeight));



        Map<Integer, PanelWithLineChart> scanChartMap = spectrumPanel.getScanLineChartMap();
        List<Integer> allScans = new ArrayList<Integer>(scanChartMap.keySet());
        Collections.sort(allScans);

        //Save the chart images, with sidebar data in case we need it.  Clunky.
        try
        {
            //the call to spectrumPanel.isSpecifiedScanFoundMS1() is a hack to find out the scan level of the ID scan
            saveChartToImageFile(spectrumPanel, quantEvent.getSpectrumFile(), sidebarWidth,
                    quantEvent.getPeptide(), quantEvent.getCharge(), quantEvent.getLightMz(), quantEvent.getHeavyMz(),
                    quantEvent.getLightIntensity(), quantEvent.getHeavyIntensity(), quantEvent.getRatio(),
                    firstLightQuantScan, lastLightQuantScan,
                    firstHeavyQuantScan, lastHeavyQuantScan, quantEvent.calcLightNeutralMass(),
                    quantEvent.getScan(), spectrumPanel.isSpecifiedScanFoundMS1() ? 1 : 2);
            ApplicationContext.infoMessage("Wrote spectrum to image " + quantEvent.getSpectrumFile().getAbsolutePath());
            saveChartToImageFile(spectrumPanel.getIntensitySumChart(), quantEvent.getIntensitySumFile(),
                    sidebarWidth,
                    quantEvent.getPeptide(), quantEvent.getCharge(), quantEvent.getLightMz(), quantEvent.getHeavyMz(),
                    quantEvent.getLightIntensity(), quantEvent.getHeavyIntensity(), quantEvent.getRatio(),
                    firstLightQuantScan, lastLightQuantScan,
                    firstHeavyQuantScan, lastHeavyQuantScan, quantEvent.calcLightNeutralMass(),
                    quantEvent.getScan(), spectrumPanel.isSpecifiedScanFoundMS1() ? 1 : 2);
            ApplicationContext.infoMessage("Wrote intensity sum to image " +
                    quantEvent.getIntensitySumFile().getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new RuntimeException("Failed to save image file",e);
        }

        if (show3DPlots)
        {
            try
            {
                //the call to spectrumPanel.isSpecifiedScanFoundMS1() is a hack to find out the scan level of the ID scan
                saveChartToImageFile(spectrumPanel.getContourPlot(), quantEvent.getFile3D(), sidebarWidth,
                        quantEvent.getPeptide(), quantEvent.getCharge(), quantEvent.getLightMz(), quantEvent.getHeavyMz(),
                        quantEvent.getLightIntensity(), quantEvent.getHeavyIntensity(), quantEvent.getRatio(),
                        firstLightQuantScan, lastLightQuantScan,
                        firstHeavyQuantScan, lastHeavyQuantScan, quantEvent.calcLightNeutralMass(),
                        quantEvent.getScan(), spectrumPanel.isSpecifiedScanFoundMS1() ? 1 : 2);
                ApplicationContext.infoMessage("Wrote 3D plot to image " + quantEvent.getFile3D().getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new RuntimeException("Failed to save 3D image file",e);
            }
        }

        try
        {
            spectrumPanel.savePerScanSpectraImage(imageWidth, scanImageHeight,
                    maxScansImageHeight, quantEvent.getScansFile());
            ApplicationContext.infoMessage("Wrote scans to image " + quantEvent.getScansFile().getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new RuntimeException("Failed to write image file " +
                    quantEvent.getScansFile().getAbsolutePath(),e);
        }

        String outChartsRelativeDirName = "";
        if (proteinsToExamine != null)
            outChartsRelativeDirName = protein + File.separatorChar;
        if (writeHTMLAndText)
        {
            outHtmlPW.println(quantEvent.createOutputRowHtml(outChartsRelativeDirName, showProteinColumn, show3DPlots));
            outTsvPW.println(quantEvent.createOutputRowTsv(showProteinColumn, show3DPlots));
        }

        
        //record event
        Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>> peptideFractionChargeFilesMap =
                proteinPeptideFractionChargeFilesMap.get(protein);
        if (peptideFractionChargeFilesMap == null)
        {
            peptideFractionChargeFilesMap = new HashMap<String, Map<String, Map<Integer, List<Pair<File, File>>>>>();
            proteinPeptideFractionChargeFilesMap.put(protein, peptideFractionChargeFilesMap);
        }
        Map<String, Map<Integer, List<Pair<File, File>>>> fractionChargeFilesMap =
                peptideFractionChargeFilesMap.get(quantEvent.getPeptide());
        if (fractionChargeFilesMap == null)
        {
            fractionChargeFilesMap = new HashMap<String, Map<Integer, List<Pair<File, File>>>>();
            peptideFractionChargeFilesMap.put(quantEvent.getPeptide(), fractionChargeFilesMap);
        }
        Map<Integer, List<Pair<File, File>>> chargeFilesMap = fractionChargeFilesMap.get(fraction);
        if (chargeFilesMap == null)
        {
            chargeFilesMap = new HashMap<Integer, List<Pair<File, File>>>();
            fractionChargeFilesMap.put(fraction, chargeFilesMap);
        }
        List<Pair<File, File>> filesList = chargeFilesMap.get(quantEvent.getCharge());
        if (filesList == null)
        {
            filesList = new ArrayList<Pair<File, File>>();
            chargeFilesMap.put(quantEvent.getCharge(), filesList);
        }
        filesList.add(new Pair<File, File>(quantEvent.getSpectrumFile(), quantEvent.getScansFile()));

        //Now is a VERY good time to do GC
        System.gc();
    }

    /**
     * Save chart to an image file, with or without the sidebar information
     * @param chartPanel
     * @param outFile
     * @param sidebarWidth
     * @param peptide
     * @param charge
     * @param lightMz
     * @param heavyMz
     * @param lightIntensity
     * @param heavyIntensity
     * @param ratio
     * @param lightMinQuantScan
     * @param lightMaxQuantScan
     * @param heavyMinQuantScan
     * @param heavyMaxQuantScan
     * @param lightMass
     * @param idScan
     * @param idScanLevel
     * @throws IOException
     */
    public void saveChartToImageFile(PanelWithChart chartPanel,
                                     File outFile, int sidebarWidth,
                                     String peptide, int charge, float lightMz, float heavyMz,
                                     float lightIntensity, float heavyIntensity, float ratio,
                                     int lightMinQuantScan, int lightMaxQuantScan,
                                     int heavyMinQuantScan, int heavyMaxQuantScan,
                                     float lightMass,
                                     int idScan, int idScanLevel) throws IOException
    {
        BufferedImage spectrumImage = chartPanel.createImage();
        BufferedImage imageToWrite = spectrumImage;

        if (writeInfoOnCharts)
        {
            int fullImageWidth = spectrumImage.getWidth() + sidebarWidth;

            imageToWrite = new BufferedImage(fullImageWidth,
                    spectrumImage.getHeight(), BufferedImage.TYPE_INT_RGB);

            Graphics2D g = imageToWrite.createGraphics();
            g.drawImage(spectrumImage, sidebarWidth, 0, null);

            //write in sidebar
            int lineHeight = 20;
            int lineNum = 1;
            int indent = 5;
            g.setPaint(Color.RED);
            g.drawString(peptide, indent, lineNum++ * lineHeight);
            g.drawString("Charge=" + charge, indent, lineNum++ * lineHeight);
            g.drawString("Light mass=" + lightMass, indent, lineNum++ * lineHeight);
            g.drawString("Light m/z=" + lightMz, indent, lineNum++ * lineHeight);
            g.drawString("Heavy m/z=" + heavyMz, indent, lineNum++ * lineHeight);
            g.drawString("Light int=" + lightIntensity, indent, lineNum++ * lineHeight);
            g.drawString("Heavy int=" + heavyIntensity, indent, lineNum++ * lineHeight);
            g.drawString("Ratio=" + ratio, indent, lineNum++ * lineHeight);
            g.drawString("MinscanLt=" + lightMinQuantScan, indent, lineNum++ * lineHeight);
            g.drawString("MaxscanLt=" + lightMaxQuantScan, indent, lineNum++ * lineHeight);
            g.drawString("MinScanHv=" + heavyMinQuantScan, indent, lineNum++ * lineHeight);
            g.drawString("MaxScanHv=" + heavyMaxQuantScan, indent, lineNum++ * lineHeight);
            g.drawString("ID scan=" + idScan, indent, lineNum++ * lineHeight);
            g.drawString("IDscan level=" + idScanLevel, indent, lineNum++ * lineHeight);

            g.dispose();
        }
        ImageIO.write(imageToWrite,"png",outFile);
    }


    public int getResolution()
    {
        return resolution;
    }

    public void setResolution(int resolution)
    {
        this.resolution = resolution;
    }

    public int getNumHeavyPeaksToPlot()
    {
        return numHeavyPeaksToPlot;
    }

    public void setNumHeavyPeaksToPlot(int numHeavyPeaksToPlot)
    {
        this.numHeavyPeaksToPlot = numHeavyPeaksToPlot;
    }

    public File getOutDir()
    {
        return outDir;
    }

    public void setOutDir(File outDir)
    {
        this.outDir = outDir;
    }

    public int getScanImageHeight()
    {
        return scanImageHeight;
    }

    public void setScanImageHeight(int scanImageHeight)
    {
        this.scanImageHeight = scanImageHeight;
    }

    public int getMaxScansImageHeight()
    {
        return maxScansImageHeight;
    }

    public void setMaxScansImageHeight(int maxScansImageHeight)
    {
        this.maxScansImageHeight = maxScansImageHeight;
    }

    public int getSpectrumImageHeight()
    {
        return spectrumImageHeight;
    }

    public void setSpectrumImageHeight(int spectrumImageHeight)
    {
        this.spectrumImageHeight = spectrumImageHeight;
    }

    public int getImageWidth()
    {
        return imageWidth;
    }

    public void setImageWidth(int imageWidth)
    {
        this.imageWidth = imageWidth;
    }

    public File getMzXmlDir()
    {
        return mzXmlDir;
    }

    public void setMzXmlDir(File mzXmlDir)
    {
        this.mzXmlDir = mzXmlDir;
    }

    public float getMaxCombineFeatureRatioDiff()
    {
        return maxCombineFeatureRatioDiff;
    }

    public void setMaxCombineFeatureRatioDiff(float maxCombineFeatureRatioDiff)
    {
        this.maxCombineFeatureRatioDiff = maxCombineFeatureRatioDiff;
    }

    public boolean isWriteHTMLAndText()
    {
        return writeHTMLAndText;
    }

    public void setWriteHTMLAndText(boolean writeHTMLAndText)
    {
        this.writeHTMLAndText = writeHTMLAndText;
    }

    public float getMinPeptideProphet()
    {
        return minPeptideProphet;
    }

    public void setMinPeptideProphet(float minPeptideProphet)
    {
        this.minPeptideProphet = minPeptideProphet;
    }

    public File getMzXmlFile()
    {
        return mzXmlFile;
    }

    public void setMzXmlFile(File mzXmlFile)
    {
        this.mzXmlFile = mzXmlFile;
    }

    public int getNumPaddingScans()
    {
        return numPaddingScans;
    }

    public void setNumPaddingScans(int numPaddingScans)
    {
        this.numPaddingScans = numPaddingScans;
    }

    public float getMzPadding()
    {
        return mzPadding;
    }

    public void setMzPadding(float mzPadding)
    {
        this.mzPadding = mzPadding;
    }

    public Set<String> getPeptidesFound()
    {
        return peptidesFound;
    }

    public void setPeptidesFound(Set<String> peptidesFound)
    {
        this.peptidesFound = peptidesFound;
    }

    public int getSidebarWidth()
    {
        return sidebarWidth;
    }

    public void setSidebarWidth(int sidebarWidth)
    {
        this.sidebarWidth = sidebarWidth;
    }

    public Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>> getProteinPeptideFractionChargeFilesMap()
    {
        return proteinPeptideFractionChargeFilesMap;
    }

    public void setProteinPeptideFractionChargeFilesMap(Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>> proteinPeptideFractionChargeFilesMap)
    {
        this.proteinPeptideFractionChargeFilesMap = proteinPeptideFractionChargeFilesMap;
    }

    public PrintWriter getOutHtmlPW()
    {
        return outHtmlPW;
    }

    public void setOutHtmlPW(PrintWriter outHtmlPW)
    {
        this.outHtmlPW = outHtmlPW;
    }

    public PrintWriter getOutTsvPW()
    {
        return outTsvPW;
    }

    public void setOutTsvPW(PrintWriter outTsvPW)
    {
        this.outTsvPW = outTsvPW;
    }

    public boolean isShow3DPlots()
    {
        return show3DPlots;
    }

    public void setShow3DPlots(boolean show3DPlots)
    {
        this.show3DPlots = show3DPlots;
    }

    public int getRotationAngle3D()
    {
        return rotationAngle3D;
    }

    public void setRotationAngle3D(int rotationAngle3D)
    {
        this.rotationAngle3D = rotationAngle3D;
    }

    public int getTiltAngle3D()
    {
        return tiltAngle3D;
    }

    public void setTiltAngle3D(int tiltAngle3D)
    {
        this.tiltAngle3D = tiltAngle3D;
    }

    public int getImageHeight3D()
    {
        return imageHeight3D;
    }

    public void setImageHeight3D(int imageHeight3D)
    {
        this.imageHeight3D = imageHeight3D;
    }

    public int getImageWidth3D()
    {
        return imageWidth3D;
    }

    public void setImageWidth3D(int imageWidth3D)
    {
        this.imageWidth3D = imageWidth3D;
    }

    public boolean isShow3DAxes()
    {
        return show3DAxes;
    }

    public void setShow3DAxes(boolean show3DAxes)
    {
        this.show3DAxes = show3DAxes;
    }

    public int getScan()
    {
        return scan;
    }

    public void setScan(int scan)
    {
        this.scan = scan;
    }

    public Set<String> getPeptidesToExamine()
    {
        return peptidesToExamine;
    }

    public void setPeptidesToExamine(Set<String> peptidesToExamine)
    {
        this.peptidesToExamine = peptidesToExamine;
    }

    public Set<String> getProteinsToExamine()
    {
        return proteinsToExamine;
    }

    public void setProteinsToExamine(Set<String> proteinsToExamine)
    {
        this.proteinsToExamine = proteinsToExamine;
    }

    public Set<String> getFractionsToExamine()
    {
        return fractionsToExamine;
    }

    public void setFractionsToExamine(Set<String> fractionsToExamine)
    {
        this.fractionsToExamine = fractionsToExamine;
    }

    public Iterator<FeatureSet> getFeatureSetIterator()
    {
        return featureSetIterator;
    }

    public void setFeatureSetIterator(Iterator<FeatureSet> featureSetIterator)
    {
        this.featureSetIterator = featureSetIterator;
    }

    public boolean isWriteInfoOnCharts()
    {
        return writeInfoOnCharts;
    }

    public void setWriteInfoOnCharts(boolean writeInfoOnCharts)
    {
        this.writeInfoOnCharts = writeInfoOnCharts;
    }

    public File getOutHtmlFile()
    {
        return outHtmlFile;
    }

    public void setOutHtmlFile(File outHtmlFile)
    {
        this.outHtmlFile = outHtmlFile;
    }

    public File getOutTsvFile()
    {
        return outTsvFile;
    }

    public void setOutTsvFile(File outTsvFile)
    {
        this.outTsvFile = outTsvFile;
    }

    public float getPeakSeparationMass()
    {
        return peakSeparationMass;
    }

    public void setPeakSeparationMass(float peakSeparationMass)
    {
        this.peakSeparationMass = peakSeparationMass;
    }

    public float getPeakTolerancePPM()
    {
        return peakTolerancePPM;
    }

    public void setPeakTolerancePPM(float peakTolerancePPM)
    {
        this.peakTolerancePPM = peakTolerancePPM;
    }

}
