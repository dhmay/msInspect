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
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
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
    protected int spectrumImageHeight = 800;
    protected int imageWidth = 1200;

    protected File pepXmlFile;
    protected File mzXmlDir;

    //The amount by which two features can differ and still be combined into the same single representation.
    //TODO: convert this to log space, parameterize
    protected float maxCombineFeatureRatioDiff = 0.2f;

    protected Set<String> peptidesToExamine;
    protected Set<String> proteinsToExamine;
    protected Set<String> fractionsToExamine;



    protected float minPeptideProphet = 0;

    protected File mzXmlFile = null;

    protected int numPaddingScans = 3;
    protected float mzPadding = 1;

    protected Set<String> peptidesFound = new HashSet<String>();
    protected int sidebarWidth = 180;

    protected Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>> proteinPeptideFractionChargeFilesMap =
            new HashMap<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>>();

    protected static final String DUMMY_PROTEIN_NAME = "DUMMY_PROTEIN";

    protected PrintWriter outHtmlPW;

    protected boolean show3DPlots = true;
    protected int rotationAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_ROTATION;
    protected int tiltAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_TILT;

    int scan = 0;


    public PeptideQuantVisualizationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "quantcharts";

        mHelpMessage ="Create charts showing the MS1 spectra near quantitative events for specified peptides or " +
                "proteins.  Includes an HTML index to all of the charts, broken down by protein, peptide, " +
                "fraction, charge, etc.";
        mShortDescription = "Create charts showing the MS1 spectra near quantitative events for specified peptides " +
                "or proteins";

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
                        createIntegerArgumentDefinition("scan", false, "Scan number of desired quantitation event",scan),
                        createStringArgumentDefinition("fractions", false, "Fraction containing desired quantitation event"),

                        createDecimalArgumentDefinition("minpprophet", false, "minimum PeptideProphet value",
                                minPeptideProphet),
                        createIntegerArgumentDefinition("paddingscans", false,
                                "number of scans before and after quant envelope to display", numPaddingScans),
                        createDecimalArgumentDefinition("mzpadding", false,
                                "amount of m/z space to display around quant", mzPadding),
                        createIntegerArgumentDefinition("numpeaksaboveheavy", false,
                                "number of peaks above the heavy-ion monoisotope to display", numHeavyPeaksToPlot),
                        createBooleanArgumentDefinition("show3dplots", false,
                                "Show 3D plot? (takes more time)", show3DPlots),
                        createIntegerArgumentDefinition("3drotation", false,
                                "Rotation angle for 3D plot", rotationAngle3D),
                        createIntegerArgumentDefinition("3dtilt", false,
                                "Tilt angle for 3D plot", tiltAngle3D),
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

        int numModeArgs = 0;
        if (hasArgumentValue("peptides")) numModeArgs++;
        if (hasArgumentValue("proteins")) numModeArgs++;
        if (hasArgumentValue("scan")) numModeArgs++;

        if (numModeArgs != 1)
            throw new ArgumentValidationException(
                    "Please supply either peptides, or proteins, or scan (but no more than one");

        if (hasArgumentValue("peptides"))
        {
            String peptidesString = getStringArgumentValue("peptides");
            String[] peptidesArray = peptidesString.split(",");
            peptidesToExamine = new HashSet<String>();
                peptidesToExamine.addAll(Arrays.asList(peptidesArray));
        }

        if (hasArgumentValue("proteins"))
        {
            String proteinsString = getStringArgumentValue("proteins");
            String[] proteinsArray = proteinsString.split(",");
            proteinsToExamine = new HashSet<String>();
            proteinsToExamine.addAll(Arrays.asList(proteinsArray));
        }

        scan = getIntegerArgumentValue("scan");
        if (hasArgumentValue("fractions"))
        {
            String fractionsString = getStringArgumentValue("fractions");
            String[] fractionsArray = fractionsString.split(",");
            fractionsToExamine = new HashSet<String>();
            fractionsToExamine.addAll(Arrays.asList(fractionsArray));
for (String fraction : fractionsToExamine) System.err.println(fraction);            
        }

        numPaddingScans = getIntegerArgumentValue("paddingscans");
        mzPadding = getFloatArgumentValue("mzpadding");
        numHeavyPeaksToPlot = getIntegerArgumentValue("numpeaksaboveheavy");

        show3DPlots = getBooleanArgumentValue("show3dplots");
        rotationAngle3D = getIntegerArgumentValue("3drotation");
        tiltAngle3D = getIntegerArgumentValue("3dtilt");
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        File outHtmlFile = new File(outDir,"quantitation.html");

        try
        {
            outHtmlPW = new PrintWriter(outHtmlFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("ERROR: Failed to open output HTML document");
        }

        outHtmlPW.println(HtmlGenerator.createDocumentBeginning("Quantitative Events"));
        String contourHeaderText = "";
        if (show3DPlots)
             contourHeaderText = "<th>3D</th>";
        outHtmlPW.print("<table border=\"1\"><tr><th>Protein</th><th>Peptide</th><th>Fraction</th><th>Charge</th>" +
                        "<th>Scan</th><th>Spectrum</th><th>Scans</th>" + contourHeaderText +
                "<th>Ratio</th><th>Light</th><th>Heavy</th>" +
                "<th>FirstScan</th><th>LastScan</th></tr>");

        outHtmlPW.flush();



        Iterator<FeatureSet> fsi = null;
        try
        {
            fsi = new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
        }
        catch (IOException e)
        {
            try
            {
                if (fractionsToExamine != null)
                    throw new CommandLineModuleExecutionException("'fractions' argument provided on a non-pepxml file.  Quitting");
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

        boolean processedAFraction = false;
        while (fsi.hasNext())
        {
            FeatureSet fraction = fsi.next();
            ApplicationContext.infoMessage("Handling fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(fraction));
            if (fractionsToExamine == null || fractionsToExamine.contains(MS2ExtraInfoDef.getFeatureSetBaseName(fraction)))
            {
                handleFraction(fraction);
                processedAFraction = true;
            }
        }
        if (!processedAFraction)
            ApplicationContext.infoMessage("WARNING: no fractions processed");

        outHtmlPW.println("</table>");
        outHtmlPW.println(HtmlGenerator.createDocumentEnd());
        outHtmlPW.flush();
        outHtmlPW.close();
        ApplicationContext.infoMessage("Saved HTML file " + outHtmlFile.getAbsolutePath());

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
        else if (peptidesToExamine != null)
            handlePeptidesInRun(featureSet, run, peptidesToExamine, DUMMY_PROTEIN_NAME, outDir);
        else
        {
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            if (fractionsToExamine == null || fractionsToExamine.contains(fractionName))
                handleScanInRun(featureSet, fractionName, run);
        }
    }

    protected void handleScanInRun(FeatureSet featureSet, String fraction, MSRun run)
            throws CommandLineModuleExecutionException
    {
        for (Feature feature : featureSet.getFeatures())
        {
            if (scan == feature.getScan())
                createChartsForEvent(run, outDir, DUMMY_PROTEIN_NAME, fraction, feature);
        }
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
                createChartsForEvent(run, outputDir, proteinName, fraction, feature);
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
                    featuresOverlappingFirst.add(compareFeature);
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

        return result;
    }

    protected void createChartsForEvent(MSRun run,File outputDir, String protein, String fraction, Feature feature)
            throws CommandLineModuleExecutionException
    {
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature);
        float heavyIntensity =                 (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
        float ratio =                 (float) IsotopicLabelExtraInfoDef.getRatio(feature);
        int charge = feature.getCharge();
        int firstQuantScan = IsotopicLabelExtraInfoDef.getLightFirstScan(feature);
        int lastQuantScan = IsotopicLabelExtraInfoDef.getLightLastScan(feature);

        if (firstQuantScan <= 0)
            firstQuantScan = feature.getScan();
        if (lastQuantScan <= 0)
            lastQuantScan = feature.getScan();

        int minScanIndex = Math.max(Math.abs(run.getIndexForScanNum(firstQuantScan)) - numPaddingScans, 0);
        int maxScanIndex = Math.min(Math.abs(run.getIndexForScanNum(lastQuantScan)) + numPaddingScans,
                run.getScanCount()-1);

        int minScan = run.getScanNumForIndex(minScanIndex);
        int maxScan = run.getScanNumForIndex(maxScanIndex);

        float lightNeutralMass = (float) IsotopicLabelExtraInfoDef.getLightMass(feature) - Spectrum.HYDROGEN_ION_MASS;
        float heavyNeutralMass = (float) IsotopicLabelExtraInfoDef.getHeavyMass(feature) - Spectrum.HYDROGEN_ION_MASS;

        float lightMz = lightNeutralMass / feature.getCharge() + Spectrum.HYDROGEN_ION_MASS;
        float heavyMz = heavyNeutralMass / feature.getCharge() + Spectrum.HYDROGEN_ION_MASS;
        float minMz = lightMz - mzPadding;
        float maxMz = heavyMz + numHeavyPeaksToPlot / feature.getCharge() + mzPadding;
        _log.debug("Building chart for feature:\n\t" + feature);

        PanelWithSpectrumChart spectrumPanel =
                new PanelWithSpectrumChart(run, minScan, maxScan, minMz, maxMz, 
                        firstQuantScan, lastQuantScan,
                        lightMz, heavyMz);
        spectrumPanel.setResolution(resolution);
        spectrumPanel.setGenerateLineCharts(true);
        spectrumPanel.setGenerate3DChart(show3DPlots);
        spectrumPanel.setScanToCheckLevel(feature.getScan());
        spectrumPanel.setName("Spectrum");
        spectrumPanel.setContourPlotRotationAngle(rotationAngle3D);
        spectrumPanel.setContourPlotTiltAngle(tiltAngle3D);


        spectrumPanel.generateChart();
        spectrumPanel.setVisible(true);

        DropdownMultiChartDisplayPanel multiChartPanelForScans = new DropdownMultiChartDisplayPanel();
        multiChartPanelForScans.setDisplaySlideshowButton(true);

        spectrumPanel.setSize(new Dimension(imageWidth, spectrumImageHeight));
        spectrumPanel.setMinimumSize(new Dimension(imageWidth, spectrumImageHeight));


        Map<Integer, PanelWithLineChart> scanChartMap = spectrumPanel.getScanLineChartMap();
        List<Integer> allScans = new ArrayList<Integer>(scanChartMap.keySet());
        Collections.sort(allScans);

        for (int scanForChart : allScans)
            multiChartPanelForScans.addChartPanel(scanChartMap.get(scanForChart));

        String filePrefix = peptide + "_" + fraction + "_" + charge + "_" + feature.getScan();
        if (protein != null && !protein.equals(DUMMY_PROTEIN_NAME))
            filePrefix = protein + "_" + filePrefix;

        File outSpectrumFile = new File(outputDir, filePrefix + "_spectrum.png");
        File outScansFile = new File(outputDir, filePrefix + "_scans.png");
        File out3DFile = new File(outputDir, filePrefix + "_3D.png");



        try
        {
            //the call to spectrumPanel.isSpecifiedScanFoundMS1() is a hack to find out the scan level of the ID scan
            saveChartToImageFileWithSidebar(spectrumPanel, outSpectrumFile, sidebarWidth,
                    peptide, charge, lightMz, heavyMz, lightIntensity, heavyIntensity, ratio,
                    firstQuantScan, lastQuantScan, lightNeutralMass,
                    feature.getScan(), spectrumPanel.isSpecifiedScanFoundMS1() ? 1 : 2);
            ApplicationContext.infoMessage("Wrote spectrum to image " + outSpectrumFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to save image file",e);
        }

        if (show3DPlots)
        {
            try
            {
                //the call to spectrumPanel.isSpecifiedScanFoundMS1() is a hack to find out the scan level of the ID scan
                saveChartToImageFileWithSidebar(spectrumPanel.getContourPlot(), out3DFile, sidebarWidth,
                        peptide, charge, lightMz, heavyMz, lightIntensity, heavyIntensity, ratio,
                        firstQuantScan, lastQuantScan, lightNeutralMass,
                        feature.getScan(), spectrumPanel.isSpecifiedScanFoundMS1() ? 1 : 2);
                ApplicationContext.infoMessage("Wrote 3D plot to image " + out3DFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to save 3D image file",e);
            }
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

        String outRelativeDirName = "";
        if (proteinsToExamine != null)
            outRelativeDirName = protein + File.separatorChar;

        String spectrumLink = HtmlGenerator.createLink(outRelativeDirName + outSpectrumFile.getName(),"Spectrum");
        String scansLink = HtmlGenerator.createLink(outRelativeDirName + outScansFile.getName(),"Scans");
        String contourLinkText = "";
        if (show3DPlots)
        {
           contourLinkText = "</td><td>" + HtmlGenerator.createLink(outRelativeDirName + out3DFile.getName(),"3D");
        }

        String proteinColumn = (proteinsToExamine != null) ? "<td>" + protein + "</td>" : "";
        outHtmlPW.print("<tr>" + proteinColumn + "<td>" + peptide + "</td><td>" +
                fraction + "</td><td>" + charge + "</td><td>" + feature.getScan() + "</td><td>" + spectrumLink + "</td><td>" +
                scansLink + contourLinkText +  "</td><td>" + ratio + "</td><td>" + lightIntensity + "</td><td>" + heavyIntensity + "</td><td>" + firstQuantScan + "</td><td>" + lastQuantScan+ "</td></tr>\n");


        //record event
        Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>> peptideFractionChargeFilesMap =
                proteinPeptideFractionChargeFilesMap.get(protein);
        if (peptideFractionChargeFilesMap == null)
        {
            peptideFractionChargeFilesMap = new HashMap<String, Map<String, Map<Integer, List<Pair<File, File>>>>>();
            proteinPeptideFractionChargeFilesMap.put(protein, peptideFractionChargeFilesMap);
        }
        Map<String, Map<Integer, List<Pair<File, File>>>> fractionChargeFilesMap =
                peptideFractionChargeFilesMap.get(peptide);
        if (fractionChargeFilesMap == null)
        {
            fractionChargeFilesMap = new HashMap<String, Map<Integer, List<Pair<File, File>>>>();
            peptideFractionChargeFilesMap.put(peptide, fractionChargeFilesMap);
        }
        Map<Integer, List<Pair<File, File>>> chargeFilesMap = fractionChargeFilesMap.get(fraction);
        if (chargeFilesMap == null)
        {
            chargeFilesMap = new HashMap<Integer, List<Pair<File, File>>>();
            fractionChargeFilesMap.put(fraction, chargeFilesMap);
        }
        List<Pair<File, File>> filesList = chargeFilesMap.get(charge);
        if (filesList == null)
        {
            filesList = new ArrayList<Pair<File, File>>();
            chargeFilesMap.put(charge, filesList);
        }
        filesList.add(new Pair<File, File>(outSpectrumFile, outScansFile));
    }

    public void saveChartToImageFileWithSidebar(PanelWithChart chartPanel,
                                                File outFile, int sidebarWidth,
                                                String peptide, int charge, float lightMz, float heavyMz,
                                                float lightIntensity, float heavyIntensity, float ratio,
                                                int minQuantScan, int maxQuantScan, float lightMass,
                                                int idScan, int idScanLevel) throws IOException
    {
        BufferedImage spectrumImage = chartPanel.createImage();

        int fullImageWidth = spectrumImage.getWidth() + sidebarWidth;

        BufferedImage imageToWrite = new BufferedImage(fullImageWidth,
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
        g.drawString("Min scan=" + minQuantScan, indent, lineNum++ * lineHeight);
        g.drawString("Max scan=" + maxQuantScan, indent, lineNum++ * lineHeight);
        g.drawString("ID scan=" + idScan, indent, lineNum++ * lineHeight);
        g.drawString("IDscan level=" + idScanLevel, indent, lineNum++ * lineHeight);



        g.dispose();

        ImageIO.write(imageToWrite,"png",outFile);
    }

}
