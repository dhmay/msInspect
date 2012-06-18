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
package org.fhcrc.cpl.viewer.quant.gui;

import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.fhcrc.cpl.viewer.quant.QuantEventAssessor;
import org.fhcrc.cpl.viewer.quant.turk.TurkUtilities;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
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


    public static final float AMINOACID_MODIFICATION_EQUALITY_MASS_TOLERANCE = 0.25f;
    public static final float AMINOACID_MODIFICATION_EQUALITY_MZ_TOLERANCE = 0.25f;


    protected Iterator<FeatureSet> featureSetIterator;

    //a non-displayed button that's used to add listeners who care when we complete each event.
    //TODO: do this in some less silly way
    protected JButton dummyProgressButton = new JButton();

    //Have to stop somewhere.  Chart will extend this number of peaks beyond the heavy monoisotope, plus slop
    int numHeavyPeaksToPlot = 4;

    protected File outDir = null;
    protected File outTsvFile = null;
    //for Mechanical Turk HITs
    protected File outTurkFile = null;
    protected File outHtmlFile = null;
    //Should we append TSV output to the existing file, if one exists?
    protected boolean appendTsvOutput = false;
    protected boolean tsvFileAlreadyExists = false;

    //Do we actually create the charts, or just create the output files?
    protected boolean shouldCreateCharts = true;
    //Should all events that are written be marked Bad?  If not, mark Unknown
    protected boolean markAllEventsBad = false;
    //should we do an automated assessment of the event?
    protected boolean shouldAssessEvents = true;

    //The assumed difference in mass between isotopic peaks
    protected float peakSeparationMass = PanelWithSpectrumChart.DEFAULT_PEAK_SEPARATION_MASS;
    //Mass tolerance around each peak, for summary chart
    protected float peakTolerancePPM = PanelWithSpectrumChart.DEFAULT_PEAK_TOLERANCE_PPM;

    //height and width of images
    protected int scanImageHeight = DEFAULT_SINGLE_SCAN_IMAGE_HEIGHT;
    protected int maxScansImageHeight = DEFAULT_MAX_SINGLE_SCANS_TOTAL_IMAGE_HEIGHT;
    protected int spectrumImageHeight = DEFAULT_SPECTRUM_IMAGE_HEIGHT;
    protected int imageWidth = DEFAULT_SPECTRUM_IMAGE_WIDTH;

    //Directory of mzXML files
    protected File mzXmlDir;
    //Single mzXML file
    protected File mzXmlFile = null;

    //The amount by which two features can differ and still be combined into the same single representation.
    //TODO: convert this to log space, parameterize
    protected float maxCombineFeatureRatioDiff = 0.2f;

    //Should we write out HTML and tsv text for our events as we build them?
    protected boolean writeHTMLAndText = true;

    //PeptideProphet minimum
    protected float minPeptideProphet = 0;

    //Scans to display around event
    protected int numPaddingScans = 3;
    //m/z padding to display around event
    protected float mzPadding = 1.5f;

    //keeps track of the unique peptides of the events plotted
    protected Set<String> peptidesFound = new HashSet<String>();

    //A ghastly structure to hold references to all of the tsv and html files created
    //In order: protein, peptide, fraction, charge, pairs of <tsv,html> files
    protected Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>>
            proteinPeptideFractionChargeFilesMap =
            new HashMap<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>>();

    //For QuantEvent objects with no associated protein
    protected static final String DUMMY_PROTEIN_NAME = "DUMMY_PROTEIN";

    //should we show 3D plots?
    protected boolean show3DPlots = true;
    //Full control of 3D plot parameters
    protected int rotationAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_ROTATION;
    protected int tiltAngle3D = PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_TILT;
    protected int imageHeight3D = DEFAULT_IMAGE_HEIGHT_3D;
    protected int imageWidth3D = DEFAULT_IMAGE_WIDTH_3D;
    protected boolean show3DAxes = true;
    //Controls whether we write event info directly to the chart images
    protected boolean writeInfoOnCharts = false;
    protected boolean writeTheoreticalPeaksOnCharts = false;
    //If we're adding sidebar information to the images, width of the sidebar
    protected int sidebarWidth = 180;

    //a single scan that we want to visualize.  This should probably be controlled somewhere else
    int scan = 0;

    //Sets of things we want to visualize.  This should probably be controlled somewhere else
    protected Set<String> peptidesToExamine;
    protected Set<String> proteinsToExamine;
    protected Set<String> fractionsToExamine;

    //PrintWriters for the output files
    protected PrintWriter outHtmlPW;
    protected PrintWriter outTsvPW;
    protected PrintWriter outTurkPW;

    //Unique identifier (per file) for Turk HITs
    protected int currentTurkID = 0;
    //URL prefix to add before all image names in Turk output file.  Should end with "/"
    protected String turkImageURLPrefix = "";

    protected int turkChartWidth = 690;
    protected int turkChartHeight = 460;


    //Should we include the Protein column in the output files?
    //TODO: get rid of this, always show protein column
    protected boolean showProteinColumn = true;

    public QuantitationVisualizer()
    {
    }

    /**
     * Visualize just the events specified
     * TODO: fold parts of this in with the no-arg visualizeQuantEvents
     * @param quantEvents
     * @throws IOException
     * @return the list of QuantEvents in the order in which they were written to the file(s)
     */
    public List<QuantEvent> visualizeQuantEvents(List<QuantEvent> quantEvents, boolean saveInProteinDirs)
            throws IOException
    {
        if (outHtmlFile == null)
            outHtmlFile = new File(outDir,"quantitation.html");
        if (outTsvFile == null)
            outTsvFile = new File(outDir,"quantitation.tsv");
        tsvFileAlreadyExists = outTsvFile.exists();
        if (writeHTMLAndText)
        {
            outHtmlPW = new PrintWriter(outHtmlFile);
            outTsvPW = new PrintWriter(new FileOutputStream(outTsvFile, appendTsvOutput));             

            //if we're appending, don't write header.  Cheating by writing it to a fake file
            PrintWriter conditionalTsvPW = outTsvPW;
            if (appendTsvOutput && tsvFileAlreadyExists)
                conditionalTsvPW = new PrintWriter(TempFileManager.createTempFile("fake_file",
                        "fake_file_for_quantvisualizer"));

            QuantEvent.writeHeader(outHtmlPW, conditionalTsvPW, showProteinColumn, show3DPlots);
            _log.debug("Wrote header to HTML file " + outHtmlFile.getAbsolutePath() + " and tsv file " +
                    outTsvFile.getAbsolutePath());
            TempFileManager.deleteTempFiles("fake_file_for_quantvisualizer");
        }
        if (outTurkFile != null)
        {
            outTurkPW = new PrintWriter(outTurkFile);
            outTurkPW.println(TurkUtilities.createTurkHITFileHeaderLine());
            outTurkPW.flush();
        }

        //map from fraction name to list of events in that fraction
        Map<String, List<QuantEvent>> fractionEventMap = new HashMap<String, List<QuantEvent>>();

        for (QuantEvent quantEvent : quantEvents)
        {
            List<QuantEvent> eventList = fractionEventMap.get(quantEvent.getFraction());
            if (eventList == null)
            {
                eventList = new ArrayList<QuantEvent>();
                fractionEventMap.put(quantEvent.getFraction(), eventList);
            }
            eventList.add(quantEvent);

        }
        //sort events by scan within fractions, to keep recently-used scans in cache
        Comparator<QuantEvent> scanAscComp = new QuantEvent.ScanAscComparator();
        for (List<QuantEvent> eventList : fractionEventMap.values())
            Collections.sort(eventList, scanAscComp);
        int numEventsProcessed = 0;

        //listeners that want to be updated when we finish an event
        ActionListener[] progressListeners = dummyProgressButton.getActionListeners();

        List<QuantEvent> resortedEvents = new ArrayList<QuantEvent>();
        for (String fraction : fractionEventMap.keySet())
        {
            File mzXmlFile = CommandLineModuleUtilities.findFileWithPrefix(fraction + ".", mzXmlDir, "mzXML");
            MSRun run = MSRun.load(mzXmlFile.getAbsolutePath());

            for (QuantEvent quantEvent : fractionEventMap.get(fraction))
            {
                resortedEvents.add(quantEvent);
                File outDirThisEvent = outDir;
                if (saveInProteinDirs && shouldCreateCharts)
                {
                    String protein = quantEvent.getProtein();
                    outDirThisEvent = new File(outDir, protein);
                    outDirThisEvent.mkdir();
                }
                handleEvent(run, outDirThisEvent, quantEvent.getProtein(), fraction, quantEvent);
                numEventsProcessed++;

                if (progressListeners != null)
                {
                    ActionEvent event = new ActionEvent(dummyProgressButton, 0,"" + numEventsProcessed);
                    for (ActionListener listener : progressListeners)
                        listener.actionPerformed(event);
                }
            }
        }

        if (writeHTMLAndText)
        {
            QuantEvent.writeFooterAndClose(outHtmlPW, outTsvPW);
            ApplicationContext.infoMessage("Saved HTML file " + outHtmlFile.getAbsolutePath());
            ApplicationContext.infoMessage("Saved TSV file " + outTsvFile.getAbsolutePath());
        }
        if (outTurkFile != null)
        {
            try
            {
                outTurkPW.close();
            }
            catch (Exception e) {}
        }

        return resortedEvents;
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

        _log.debug("visualizeQuantEvents begin, write HTML and text? " + writeHTMLAndText);
        if (writeHTMLAndText)
        {
            outHtmlPW = new PrintWriter(outHtmlFile);
            outTsvPW = new PrintWriter(new FileOutputStream(outTsvFile, appendTsvOutput));
            _log.debug("opened HTML file " + outHtmlFile.getAbsolutePath() + " and tsv file " +
                    outTsvFile.getAbsolutePath() + " for writing.");

            //if we're appending, don't write header.  Cheating by writing it to a fake file
            PrintWriter conditionalTsvPW = outTsvPW;
            if (appendTsvOutput && tsvFileAlreadyExists)
                conditionalTsvPW = new PrintWriter(TempFileManager.createTempFile("fake_file",
                        "fake_file_for_quantvisualizer"));

            QuantEvent.writeHeader(outHtmlPW, conditionalTsvPW, showProteinColumn, show3DPlots);
            _log.debug("Wrote HTML and TSV header");
            TempFileManager.deleteTempFiles("fake_file_for_quantvisualizer");
        }

        boolean processedAFraction = false;
        while (featureSetIterator.hasNext())
        {
            FeatureSet fraction = featureSetIterator.next();
            ApplicationContext.infoMessage("Evaluating fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(fraction));
            if (fractionsToExamine == null || fractionsToExamine.contains(MS2ExtraInfoDef.getFeatureSetBaseName(fraction)))
            {
                ApplicationContext.infoMessage("Handling fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(fraction));
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
            QuantEvent.writeFooterAndClose(outHtmlPW, outTsvPW);
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
        }

        run = MSRun.load(fileToLoad.getAbsolutePath());
        _log.debug("\tLoaded mzXML file " + fileToLoad.getAbsolutePath());
        if (proteinsToExamine != null)
            handleProteinsInRun(featureSet, run, proteinsToExamine);
        else if (peptidesToExamine != null)
            handlePeptidesInRun(featureSet, run, peptidesToExamine, DUMMY_PROTEIN_NAME, outDir);
        else if (scan != 0)
        {
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            if (fractionsToExamine == null || fractionsToExamine.contains(fractionName))
                handleScanInRun(featureSet, fractionName, run, scan);
        }
        else
        {
            _log.debug("\t processing all scans");
            //handle all scans
            String fractionName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            List<QuantEvent> allQuantEvents = new ArrayList<QuantEvent>();
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    allQuantEvents.add(new QuantEvent(feature, fractionName));
            }
            if (allQuantEvents.isEmpty())
                ApplicationContext.infoMessage("\tSkipping empty fraction " + fractionName);
            else
            {
                List<QuantEvent> nonOverlappingQuantEvents = findNonOverlappingEvents(allQuantEvents);
                for (QuantEvent quantEvent : nonOverlappingQuantEvents)
                {
                    handleEvent(run, outDir, DUMMY_PROTEIN_NAME, fractionName, quantEvent);
                }
                ApplicationContext.infoMessage("\tProcessed " + nonOverlappingQuantEvents.size() + " non-overlapping events for " +
                        featureSet.getFeatures().length + " features");
            }
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
    protected void handleScanInRun(FeatureSet featureSet, String fraction, MSRun run, int scanToHandle)
    {
        for (Feature feature : featureSet.getFeatures())
        {
            if (scanToHandle == feature.getScan())
                handleFeature(run, outDir, DUMMY_PROTEIN_NAME, fraction, feature);
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
                {
                    String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                    {
                      
                        ApplicationContext.infoMessage("Found peptide " + peptide + ", protein " + protein);
                        peptidesThisProtein.add(peptide);
                    }
                }

            }

            if (!peptidesThisProtein.isEmpty())
            {
                ApplicationContext.infoMessage("Protein " + protein + ": Finding " + peptidesThisProtein.size() +
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

        List<QuantEvent> allQuantEventsAllPeptides = new ArrayList<QuantEvent>();
//        String labeledResidue = null;
//        float labelMassDiff = 0;
        for (Feature feature : featureSet.getFeatures())
        {
            if (peptidesToHandle.contains(MS2ExtraInfoDef.getFirstPeptide(feature)) &&
                    MS2ExtraInfoDef.getPeptideProphet(feature) >= this.minPeptideProphet &&
                    IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
//                if (labeledResidue == null)
//                {
//                    AnalyzeICAT.IsotopicLabel label = IsotopicLabelExtraInfoDef.getLabel(feature);
//                    if (label != null)
//                    {
//                        labeledResidue = "" + label.getResidue();
//                        labelMassDiff = label.getHeavy() - label.getLight();
//                        _log.debug("Found label: " + labeledResidue + ", " + labelMassDiff);
//                    }
//                }
                allQuantEventsAllPeptides.add(new QuantEvent(feature, fraction));
            }
        }
//        if (labeledResidue == null)
//            ApplicationContext.infoMessage("WARNING: unable to determine modification used for quantitation.  " +
//                    "Cannot collapse light and heavy states.");

        List<QuantEvent> nonOverlappingEventsAllPeptides =
                findNonOverlappingQuantEventsAllPeptides(allQuantEventsAllPeptides);

        for (QuantEvent quantEvent : nonOverlappingEventsAllPeptides)
        {
            int numTotalEvents = 1;
            if (quantEvent.getOtherEvents() != null)
                numTotalEvents += quantEvent.getOtherEvents().size();
            ApplicationContext.infoMessage("\tHandling peptide " +  quantEvent.getPeptide() +
                    ", charge " + quantEvent.getCharge() + " with " + numTotalEvents + " events");
            handleEvent(run, outputDir, proteinName, fraction, quantEvent);
        }
    }

    /**
     * Find overlapping events in a list of events that may come from different peptides, fractions, etc.
     * Return a list with one member of each of those overlapping event lists
     * @param quantEvents
     * @return
     */
    public List<QuantEvent> findNonOverlappingQuantEventsAllPeptides(
            List<QuantEvent> quantEvents)
    {
        List<QuantEvent> result = new ArrayList<QuantEvent>();
        Map<String, Map<String, Map<Integer, List<QuantEvent>>>>
                peptideFractionChargeQuantEventsMap =
                new HashMap<String, Map<String, Map<Integer, List<QuantEvent>>>>();
        _log.debug("findNonOverlappingQuantEventsAllPeptides 1");
        for (QuantEvent quantEvent : quantEvents)
        {
            String peptide = quantEvent.getPeptide();

            Map<String, Map<Integer, List<QuantEvent>>> fractionChargesMap =
                    peptideFractionChargeQuantEventsMap.get(peptide);
            if (fractionChargesMap == null)
            {
                fractionChargesMap = new HashMap<String, Map<Integer, List<QuantEvent>>>();
                peptideFractionChargeQuantEventsMap.put(peptide, fractionChargesMap);
            }

            String fraction = quantEvent.getFraction();

            Map<Integer, List<QuantEvent>> chargeEventsMap =
                    fractionChargesMap.get(fraction);
            if (chargeEventsMap == null)
            {
                chargeEventsMap = new HashMap<Integer, List<QuantEvent>>();
                fractionChargesMap.put(fraction, chargeEventsMap);
            }

            List<QuantEvent> eventsToEvaluate = chargeEventsMap.get(quantEvent.getCharge());

            if (eventsToEvaluate == null)
            {
                eventsToEvaluate = new ArrayList<QuantEvent>();
                chargeEventsMap.put(quantEvent.getCharge(), eventsToEvaluate);
            }
            eventsToEvaluate.add(quantEvent);
        }
        _log.debug("Built map with " + peptideFractionChargeQuantEventsMap.size() + " peptides");

        //combine modification states that represent light and heavy versions of same species.
        //To do this, first break up everything by peptide, fraction and charge.  Then associate everything
        //within a charge by light mz and find overlapping events within that list.  Don't bother looking at
        //heavy -- we're not dealing with multiple-label scenarios.
        for (String peptide : peptideFractionChargeQuantEventsMap.keySet())
        {
            _log.debug("Map peptide " + peptide);
            Map<String, Map<Integer, List<QuantEvent>>> fractionChargesMap =
                    peptideFractionChargeQuantEventsMap.get(peptide);
            for (String fraction : fractionChargesMap.keySet())
            {
                _log.debug("  Map fraction " + fraction);
                Map<Integer, List<QuantEvent>> chargeEventsMap =
                        fractionChargesMap.get(fraction);
                for (int charge : chargeEventsMap.keySet())
                {
                    _log.debug("  Map charge " + charge + ", " + chargeEventsMap.get(charge) + " events");
                    Map<Float, List<QuantEvent>> lightMzEventsMap = new HashMap<Float, List<QuantEvent>>();
                    for (QuantEvent event : chargeEventsMap.get(charge))
                    {
                        boolean foundIt = false;
                        for (Float lightMz : lightMzEventsMap.keySet())
                        {
                            if (Math.abs(event.getLightMz() - lightMz) < AMINOACID_MODIFICATION_EQUALITY_MZ_TOLERANCE)
                            {
                                lightMzEventsMap.get(lightMz).add(event);
                                foundIt = true;
                                break;
                            }
                        }
                        if (!foundIt)
                        {
                            List<QuantEvent> eventList = new ArrayList<QuantEvent>();
                            eventList.add(event);
                            lightMzEventsMap.put(event.getLightMz(), eventList);
                        }
                    }
                    _log.debug("Collapsing events from " + lightMzEventsMap.size() + " distinct mzs");
                    for (List<QuantEvent> eventList : lightMzEventsMap.values())
                    {
                        List<QuantEvent> nonOverlappingList = findNonOverlappingEvents(eventList);
                        if (nonOverlappingList.size() < eventList.size())
                            _log.debug("    Reduced " + eventList.size() + " events to " + nonOverlappingList.size() +
                                    " non-overlapping events.");
                        result.addAll(nonOverlappingList);
                    }
                }
            }
        }



//        for (String peptide : peptideFractionChargeModsQuantEventsMap.keySet())
//        {
//            Map<String, Map<Integer, Map<String, List<QuantEvent>>>> fractionChargesMap = peptideFractionChargeModsQuantEventsMap.get(peptide);
//            _log.debug("processing peptide " + peptide + " with " + fractionChargesMap.size() + " fractions");
//
//            for (String fraction : fractionChargesMap.keySet())
//            {
//                Map<Integer, Map<String, List<QuantEvent>>> chargeModificationsMap =
//                        fractionChargesMap.get(fraction);
//                for (int charge : chargeModificationsMap.keySet())
//                {
//                    Map<String, List<QuantEvent>> modificationsEventsMap = chargeModificationsMap.get(charge);
//                    for (String modifications : modificationsEventsMap.keySet())
//                    {
//                        List<QuantEvent> eventList = modificationsEventsMap.get(modifications);
//                        _log.debug("\tCharge " + charge + ", mods " + modifications + ": " + eventList.size() + " events");
//                        result.addAll(findNonOverlappingEvents(eventList));
//                    }
//                }
//            }
//        }
        return result;
    }

    /**
     * Find modification states that are equivalent except for modification masses representing the difference
     * between light and heavy labels.
     */
//    protected List<Pair<String,String>> pairLightAndHeavyModificationStates(Collection<String> modificationStates,
//                                                          String labeledResidue, float labelMassDiff)
//    {
//        ModResidueMassAscComparator modResidueMassAscComparator = new ModResidueMassAscComparator();
//
//        List<Pair<String,String>> lightHeavyModStatePairs = new ArrayList<Pair<String,String>>();
//        Set<String> alreadyAssignedCompareModStates = new HashSet<String>();
//        for (String baseModificationState : modificationStates)
//        {
//            _log.debug("      pairLightAndHeavy, checking " + baseModificationState);
//            Map<Integer, List<ModifiedAminoAcid>> baseModStateMap =
//                    MS2ExtraInfoDef.parsePositionModifiedAminoAcidListMapString(baseModificationState);
//            for (String compareModificationState : modificationStates)
//            {
//                if (baseModificationState.equals(compareModificationState) ||
//                        alreadyAssignedCompareModStates.contains(compareModificationState))
//                    continue;
//                boolean foundDifference = false;
//                Map<Integer, List<ModifiedAminoAcid>> compareModStateMap =
//                        MS2ExtraInfoDef.parsePositionModifiedAminoAcidListMapString(
//                                compareModificationState);
//                if (baseModStateMap.size() != compareModStateMap.size())
//                {
//                    continue;
//                }
//                for (int position : baseModStateMap.keySet())
//                {
//                    if (!compareModStateMap.containsKey(position))
//                    {
//                        foundDifference = true;
//                        break;
//                    }
//                    List<ModifiedAminoAcid> baseModsThisPosition = baseModStateMap.get(position);
//                    List<ModifiedAminoAcid> compareModsThisPosition = compareModStateMap.get(position);
//
//                    if (baseModsThisPosition.size() != compareModsThisPosition.size())
//                    {
//                        foundDifference = true;
//                        break;
//                    }
//                    Collections.sort(baseModsThisPosition, modResidueMassAscComparator);
//                    Collections.sort(compareModsThisPosition, modResidueMassAscComparator);
//
//                    for (int i=0; i<baseModsThisPosition.size(); i++)
//                    {
//                        ModifiedAminoAcid baseMod = baseModsThisPosition.get(i);
//                        ModifiedAminoAcid compareMod = compareModsThisPosition.get(i);
//                        if (!baseMod.getAminoAcidAsString().equals(compareMod.getAminoAcidAsString()))
//                        {
//                            foundDifference = true;
//                            break;
//                        }
//                        float absMassDiff = (float) Math.abs(baseMod.getMass() - compareMod.getMass());
//                        if (absMassDiff > AMINOACID_MODIFICATION_EQUALITY_MASS_TOLERANCE)
//                        {
//                            if ((!baseMod.getAminoAcidAsString().equals(labeledResidue)) ||
//                                    (absMassDiff - labelMassDiff > AMINOACID_MODIFICATION_EQUALITY_MASS_TOLERANCE))
//                            {
//                                foundDifference = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//                if (foundDifference)
//                    continue;
//                alreadyAssignedCompareModStates.add(baseModificationState);
//                alreadyAssignedCompareModStates.add(compareModificationState);
//                lightHeavyModStatePairs.add(
//                        new Pair<String,String>(baseModificationState, compareModificationState));
//                break;
//            }
//        }
//        return lightHeavyModStatePairs;
//    }

    protected class ModResidueMassAscComparator implements Comparator<ModifiedAminoAcid>
    {
         public int compare (ModifiedAminoAcid o1, ModifiedAminoAcid o2)
         {
             int result = (o1.getAminoAcidAsString().compareTo(o2.getAminoAcidAsString()));
             if (result == 0)
             {
                 if (o1.getMass() > o2.getMass()) return 1;
                 if (o1.getMass() < o2.getMass()) return -1;
                 return 0;
             }
             return result;
         }
    }

    /**
     * cover method.  Turn the features into events, generate the nonoverlapping events, and
     * return the feature results
     * @param features
     * @return
     */
    public List<Feature> findNonOverlappingFeatures(List<Feature> features)
    {
        List<QuantEvent> quantEvents = new ArrayList<QuantEvent>();
        for (Feature feature : features)
            quantEvents.add(new QuantEvent(feature, ""));
        List<QuantEvent> nonOverlappingEvents = findNonOverlappingEvents(quantEvents);
        List<Feature> result = new ArrayList<Feature>();
        for (QuantEvent quantEvent : nonOverlappingEvents)
        {
            Feature representativeFeature = quantEvent.getSourceFeature();
            if (quantEvent.getOtherEvents() != null && !quantEvent.getOtherEvents().isEmpty())
            {
                List<Spectrum.Peak> otherFeaturesAsPeaks = new ArrayList<Spectrum.Peak>();
                for (QuantEvent otherEvent : quantEvent.getOtherEvents())
                    otherFeaturesAsPeaks.add(otherEvent.getSourceFeature());
                representativeFeature.comprised =
                        otherFeaturesAsPeaks.toArray(new Spectrum.Peak[otherFeaturesAsPeaks.size()]);
            }
            result.add(representativeFeature);
        }
        return result;
    }

    /**
     * Given a list of features (presumably with the same peptide ID, in the same charge and
     * modification state), collapse them down to a single feature representing each set of
     * overlapping events
     * @param quantEvents
     * @return
     */
    public List<QuantEvent> findNonOverlappingEvents(List<QuantEvent> quantEvents)
    {
        List<QuantEvent> result = new ArrayList<QuantEvent>();

        while (!quantEvents.isEmpty())
        {
            List<QuantEvent> eventsOverlappingFirst = findEventsOverlappingFirst(quantEvents);
            QuantEvent firstRepresentative = eventsOverlappingFirst.get(0);

            firstRepresentative.setOtherEvents(new ArrayList<QuantEvent>());
            for (int i=1; i<eventsOverlappingFirst.size(); i++)
                firstRepresentative.getOtherEvents().add(eventsOverlappingFirst.get(i));
            result.add(firstRepresentative);
        }
        return result;
    }

    /**
     * Find all events overlapping in scans with the first event
     * SIDE EFFECT: Remove them all from eventsWithinCharge
     * @param eventsWithinCharge
     * @return
     */
    public List<QuantEvent> findEventsOverlappingFirst(List<QuantEvent> eventsWithinCharge)
    {
        if (eventsWithinCharge.size() == 1)
        {
            List<QuantEvent> result = new ArrayList<QuantEvent>(eventsWithinCharge);
            eventsWithinCharge.remove(0);
            return result;
        }

        QuantEvent firstEvent = eventsWithinCharge.get(0);
        eventsWithinCharge.remove(firstEvent);
        List<QuantEvent> eventsOverlappingFirst = new ArrayList<QuantEvent>();
        eventsOverlappingFirst.add(firstEvent);
        eventsOverlappingFirst.addAll(findEventsOverlappingEvent(firstEvent, eventsWithinCharge));

        //find the "median" feature by scan (round down if even number)
        Collections.sort(eventsOverlappingFirst, new QuantEvent.ScanAscComparator());
        int numOverlapping = eventsOverlappingFirst.size();
        int medianEventIndex = (numOverlapping)/2;
        if (numOverlapping % 2 == 1)
            medianEventIndex = (numOverlapping-1)/2;

        List<QuantEvent> sortedForResult = new ArrayList<QuantEvent>();
        sortedForResult.add(eventsOverlappingFirst.get(medianEventIndex));
        for (int i=0; i<eventsOverlappingFirst.size(); i++)
            if (i != medianEventIndex)
                sortedForResult.add(eventsOverlappingFirst.get(i));

        QuantEvent representativeEvent = sortedForResult.get(0);
        _log.debug("\t" + representativeEvent.getPeptide() +
                ", fraction " + representativeEvent.getFraction() + ", charge " +
                representativeEvent.getCharge() + ", ratio " +
                representativeEvent.getRatio() + ", scans " +
                representativeEvent.getFirstLightQuantScan() + "-" +
                representativeEvent.getLastLightQuantScan() + ", represents " +
                sortedForResult.size() + " event(s)");
        if (sortedForResult.size() > 1)
            for (int i=1; i<sortedForResult.size(); i++)
                _log.debug("\t\tother event " + i +
                        ", ratio=" + sortedForResult.get(i).getRatio());
        return sortedForResult;
    }

    /**
     * Find all events overlapping in scans with the base event
     * SIDE EFFECT: Remove them all from otherEvents
     * @param baseEvent
     * @param otherEvents
     * @return
     */
    public List<QuantEvent> findEventsOverlappingEvent(QuantEvent baseEvent,
                                                           List<QuantEvent> otherEvents)
    {
        List<QuantEvent> eventsOverlappingBaseEvent = new ArrayList<QuantEvent>();
        if (otherEvents.size() > 0)
        {
            //take a broad interpretation of overlap -- overlapping either light or heavy
            int firstEventStart =
                    Math.min(baseEvent.getFirstLightQuantScan(), baseEvent.getFirstHeavyQuantScan());
            int firstEventEnd =
                    Math.max(baseEvent.getLastLightQuantScan(), baseEvent.getLastHeavyQuantScan());
            if (firstEventStart <= 0)
                firstEventStart = baseEvent.getScan();
            if (firstEventEnd <= 0)
                firstEventEnd = baseEvent.getScan();
            double firstFeatureRatio = baseEvent.getRatio();
            _log.debug("Looking for events overlapping " + firstEventStart + "-" + firstEventEnd + ", ratio " + firstFeatureRatio + ", out of " + otherEvents.size());


            for (QuantEvent compareEvent : otherEvents)
            {
                int compareFeatureStart =
                        Math.min(compareEvent.getFirstLightQuantScan(), compareEvent.getFirstHeavyQuantScan());
                int compareFeatureEnd =
                        Math.max(compareEvent.getLastLightQuantScan(), compareEvent.getLastHeavyQuantScan());
                if (compareFeatureStart <= 0)
                    compareFeatureStart = compareEvent.getScan();
                if (compareFeatureEnd <= 0)
                    compareFeatureEnd = compareEvent.getScan();
                _log.debug("\tChecking " + compareFeatureStart + "-" + compareFeatureEnd + ", ratio " + compareEvent.getRatio());

                //first check ratio similarity, because that's simplest
                if (Math.abs(compareEvent.getRatio() - firstFeatureRatio) >
                        maxCombineFeatureRatioDiff)
                    continue;


                if (Math.max(firstEventStart, compareFeatureStart) < Math.min(firstEventEnd, compareFeatureEnd))
                {
                    _log.debug("\t\tMatch!");
                    eventsOverlappingBaseEvent.add(compareEvent);
                }
            }
            otherEvents.removeAll(eventsOverlappingBaseEvent);
        }

//        result.comprised = otherFeaturesAsPeaks.toArray(new Spectrum.Peak[otherFeaturesAsPeaks.size()]);
        return eventsOverlappingBaseEvent;
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
    protected void handleFeature(MSRun run,File outputDir, String protein, String fraction,
                                        Feature feature)
    {
        QuantEvent quantEvent =
                new QuantEvent(feature, fraction);
        handleEvent(run, outputDir, protein, fraction, quantEvent);
    }

    /**
     * Given a QuantEvent and a run, create the PanelWithSpectrumChart object that will create all the charts,
     * and generate them
     * @param run
     * @param quantEvent
     * @return
     */
    public PanelWithSpectrumChart createPanelWithSpectrumChart(MSRun run, QuantEvent quantEvent)
    {
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
        List<Integer> otherEventScans = new ArrayList<Integer>();
        List<Float> otherEventMzs = new ArrayList<Float>();
        if (quantEvent.getOtherEvents() != null)
            for (QuantEvent otherEvent : quantEvent.getOtherEvents())
            {
                otherEventScans.add(otherEvent.getScan());
                otherEventMzs.add(otherEvent.getMz());
            }
        spectrumPanel.setOtherEventScans(otherEventScans);
        spectrumPanel.setOtherEventMZs(otherEventMzs);
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

        return spectrumPanel;
    }

    /**
     * Side effect!  Sets QuantEvent.ratioOnePeak, but only if charts are created
     * @param run
     * @param outputDir
     * @param protein
     * @param fraction
     * @param quantEvent
     */
    protected void handleEvent(MSRun run,File outputDir, String protein, String fraction,
                               QuantEvent quantEvent)
    {
        if (markAllEventsBad)
        {
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_BAD);
        }
        if (shouldAssessEvents)
        {
            //only assess if not already assessed
            if (quantEvent.getAlgorithmicAssessment() == null)
            {
                QuantEventAssessor eventAssessor = new QuantEventAssessor();
//                eventAssessor.setLabelType(labelType);
                eventAssessor.assessQuantEvent(quantEvent, run);
//System.err.println("*******ASSESSED!!!!! " + quantEvent.getAlgorithmicAssessment());
            }
//else System.err.println("NOT ASSESSING!");            
        }
        if (shouldCreateCharts || outTurkPW != null)
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

            PanelWithSpectrumChart spectrumPanel = createPanelWithSpectrumChart(run, quantEvent);
            //The single-peak ratio calculated by QuantEventAssessor is better than the slapdash one I calculate
            //in spectrumPanel, so if we've got that, use it.
            if (quantEvent.getAlgorithmicAssessment() != null)
                quantEvent.setRatioOnePeak(quantEvent.getAlgorithmicAssessment().getSinglePeakRatio());
            else
                quantEvent.setRatioOnePeak(spectrumPanel.getRatioOnePeak());

            Map<Integer, PanelWithLineChart> scanChartMap = spectrumPanel.getScanLineChartMap();
            List<Integer> allScans = new ArrayList<Integer>(scanChartMap.keySet());
            Collections.sort(allScans);

            if (outTurkPW != null)
            {
                try
                {

                    outTurkPW.println(saveTurkImage(quantEvent, spectrumPanel, outDir, currentTurkID));
                    outTurkPW.flush();
                    currentTurkID++;
                }
                catch (Exception e)
                {
                    throw new RuntimeException("Failed to save image file",e);
                }
            }

            if (shouldCreateCharts)
            {
                //Save the chart images, with sidebar data in case we need it.  Clunky.
                File currentImageFile = new File("");
                try
                {
                    currentImageFile = quantEvent.getSpectrumFile();
                    saveChartToImageFile(quantEvent, spectrumPanel, currentImageFile,
                            writeInfoOnCharts, 0, 0, false);
                    ApplicationContext.infoMessage("Wrote spectrum to image " +
                            quantEvent.getSpectrumFile().getAbsolutePath());
                    currentImageFile = quantEvent.getIntensitySumFile();
                    saveChartToImageFile(quantEvent, spectrumPanel.getIntensitySumChart(),
                            currentImageFile,
                            writeInfoOnCharts, 0, 0, false);
                    ApplicationContext.infoMessage("Wrote intensity sum to image " +
                            quantEvent.getIntensitySumFile().getAbsolutePath());
                    currentImageFile = quantEvent.getScansFile();
                    spectrumPanel.savePerScanSpectraImage(imageWidth, scanImageHeight,
                            maxScansImageHeight, currentImageFile);
                    ApplicationContext.infoMessage("Wrote scans to image " + quantEvent.getScansFile().getAbsolutePath());
                }
                catch (Exception e)
                {
                    throw new RuntimeException("Failed to save image file " + currentImageFile.getAbsolutePath(),e);
                }
            }

            if (show3DPlots)
            {
                try
                {
                    saveChartToImageFile(quantEvent, spectrumPanel.getContourPlot(),
                            quantEvent.getFile3D(), writeInfoOnCharts, 0, 0, false);
                    ApplicationContext.infoMessage("Wrote 3D plot to image " + quantEvent.getFile3D().getAbsolutePath());
                }
                catch (Exception e)
                {
                    throw new RuntimeException("Failed to save 3D image file",e);
                }
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
        } //end if shouldCreateCharts

        String outChartsRelativeDirName = "";
        if (proteinsToExamine != null)
            outChartsRelativeDirName = protein + File.separatorChar;
        if (writeHTMLAndText)
        {
            outHtmlPW.println(quantEvent.createOutputRowHtml(outChartsRelativeDirName, showProteinColumn, show3DPlots));
            outHtmlPW.flush();
            outTsvPW.println(quantEvent.createOutputRowTsv(showProteinColumn, show3DPlots));
            outTsvPW.flush();
        }

    }

    /**
     * Save the intensity-sum image to a file, appropriately sized for the Mechanical Turk.
     * Return a String that's appropriate to be one line of a turk HIT file
     * @param quantEvent
     * @param spectrumPanel
     * @param outDir
     * @param turkId
     * @return
     * @throws IOException
     */
    public String saveTurkImage(QuantEvent quantEvent, PanelWithSpectrumChart spectrumPanel, File outDir, int turkId)
            throws IOException
    {
        String imageFileName = TurkUtilities.createTurkImageFileName(turkId);
        saveChartToImageFile(quantEvent, spectrumPanel.getIntensitySumChart(),
                new File(outDir, imageFileName), true, turkChartWidth, turkChartHeight, true);
        return TurkUtilities.createTurkHITFileLine(quantEvent, turkId, turkImageURLPrefix);
    }



    public static PanelWithPeakChart buildTheoreticalPeakChart(QuantEvent quantEvent, int chartWidth, int chartHeight)
    {
        return buildTheoreticalPeakChart(quantEvent.getLightMz(), quantEvent.getHeavyMz(), quantEvent.getCharge(),
                quantEvent.getRatio(), chartWidth, chartHeight);
    }

    /**
     * Build a peaks chart showing the theoretical distribution, assuming the ratio is correct
     * @param lightMz
     * @param heavyMz
     * @param charge
     * @param ratio
     * @param chartWidth
     * @param chartHeight
     * @return
     */
    public static PanelWithPeakChart buildTheoreticalPeakChart(float lightMz, float heavyMz, int charge, float ratio,
                                                        int chartWidth, int chartHeight)
    {
        float lightNeutralMass = (lightMz - Spectrum.HYDROGEN_ION_MASS) *charge;
        //Hardcoded 6 is number of peakd returned by Poisson()
        //Can't just use the result of Poisson(), as that's a static array and we're gonna mess with it
        float[] lightTheoreticalPeaks = new float[6];
        System.arraycopy(Spectrum.Poisson(lightNeutralMass), 0, lightTheoreticalPeaks, 0,
                lightTheoreticalPeaks.length);
//System.err.println("***LIGHT, mz=" + lightMz + ", mass=" + lightNeutralMass + ", charge=" + charge);
//for (float peak : lightTheoreticalPeaks)
//    System.err.println("\t" + peak);
        float[] lightPeakMzs = new float[6];
        for (int i=0; i<6; i++)
            lightPeakMzs[i] = lightMz + (Spectrum.HYDROGEN_ION_MASS * i / charge);
        PanelWithPeakChart theoreticalPeaksChart = new PanelWithPeakChart(lightPeakMzs, lightTheoreticalPeaks,
                "Theoretical Peaks");

        float heavyNeutralMass = (heavyMz - Spectrum.HYDROGEN_ION_MASS) * charge;

        float[] heavyTheoreticalPeaks = new float[6];
        System.arraycopy(Spectrum.Poisson(heavyNeutralMass), 0, heavyTheoreticalPeaks, 0, heavyTheoreticalPeaks.length);
        //fix 0 ratios at 0.001 to avoid divide by 0
        for (int i=0; i<heavyTheoreticalPeaks.length; i++)
            heavyTheoreticalPeaks[i] *= 1 / Math.max(ratio, 0.001f);

        float[] heavyPeakMzs = new float[6];
        for (int i=0; i<6; i++)
            heavyPeakMzs[i] = heavyMz +
                    (Spectrum.HYDROGEN_ION_MASS * i  / charge);
        //Adjust heavy peaks if light peaks intrude.  Light appears in front of heavy
        for (int i=0; i<heavyPeakMzs.length; i++)
        {
            for (int j=0; j<lightPeakMzs.length; j++)
                if (heavyPeakMzs[i] - lightPeakMzs[j] < 0.1)
                    heavyTheoreticalPeaks[i] += lightTheoreticalPeaks[j];
        }
//System.err.println("***HEAVY, mz=" + heavyMz + ", mass=" + heavyNeutralMass);
//for (float peak : heavyTheoreticalPeaks)
//    System.err.println("\t" + peak);
        theoreticalPeaksChart.addData(heavyPeakMzs, heavyTheoreticalPeaks, "heavy");

        theoreticalPeaksChart.setPreferredSize(new Dimension(chartWidth, chartHeight));
        theoreticalPeaksChart.setSize(new Dimension(chartWidth, chartHeight));
        theoreticalPeaksChart.getChart().removeLegend();

        //remove axes from chart
        ((XYPlot)theoreticalPeaksChart.getPlot()).getDomainAxis().setVisible(false);
        ((XYPlot)theoreticalPeaksChart.getPlot()).getRangeAxis().setVisible(false);

        return theoreticalPeaksChart;
    }



    /**
     * Cover method, harvests everything from quantEvent.  Still silly, because quantEvent is also passed
     * @param quantEvent
     * @param chartPanel
     * @throws IOException
     */
    public void saveChartToImageFile(QuantEvent quantEvent, PanelWithChart chartPanel,
                                     File outputFile, boolean writeInfo, int width, int height, boolean overrideSize) throws IOException
    {
        saveChartToImageFile(chartPanel, outputFile, sidebarWidth,
                quantEvent.getCharge(), quantEvent.getLightMz(), quantEvent.getHeavyMz(),
                quantEvent.getRatio(),
                writeInfo,
                width, height, overrideSize);
    }

    /**
     * Save chart to an image file, with or without the sidebar information and/or theoretical peaks
     * @param chartPanel
     * @param outFile
     * @param sidebarWidth
     * @param charge
     * @param lightMz
     * @param heavyMz
     * @param ratio
     * @throws IOException
     */
    public void saveChartToImageFile(PanelWithChart chartPanel,
                                     File outFile, int sidebarWidth,
                                     int charge, float lightMz, float heavyMz,
                                     float ratio,
                                     boolean writeChartInfo,
                                     int width, int height, boolean overrideSize) throws IOException
    {
        BufferedImage spectrumImage = null;
        if (overrideSize)
            spectrumImage = chartPanel.createImage(width, height);
        else
            spectrumImage = chartPanel.createImage();
        BufferedImage imageToWrite = spectrumImage;

        if (writeChartInfo)
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


            if (writeChartInfo)
            {
                g.setPaint(Color.WHITE);
//                g.drawString(peptide, indent, lineNum++ * lineHeight);
//                g.drawString("Charge=" + charge, indent, lineNum++ * lineHeight);
//                g.drawString("Light mass=" + lightMass, indent, lineNum++ * lineHeight);
//                g.drawString("Light m/z=" + lightMz, indent, lineNum++ * lineHeight);
//                g.drawString("Heavy m/z=" + heavyMz, indent, lineNum++ * lineHeight);
//                g.drawString("Light int=" + lightIntensity, indent, lineNum++ * lineHeight);
//                g.drawString("Heavy int=" + heavyIntensity, indent, lineNum++ * lineHeight);
                g.drawString("Ratio=" + ratio, indent, lineNum++ * lineHeight);


//                g.drawString("MinscanLt=" + lightMinQuantScan, indent, lineNum++ * lineHeight);
//                g.drawString("MaxscanLt=" + lightMaxQuantScan, indent, lineNum++ * lineHeight);
//                g.drawString("MinScanHv=" + heavyMinQuantScan, indent, lineNum++ * lineHeight);
//                g.drawString("MaxScanHv=" + heavyMaxQuantScan, indent, lineNum++ * lineHeight);
//                g.drawString("ID scan=" + idScan, indent, lineNum++ * lineHeight);
//                g.drawString("IDscan level=" + idScanLevel, indent, lineNum++ * lineHeight);

                //theoretical peaks in bottom left
                int theoreticalPeaksHeight = (int) (sidebarWidth * 2.0 / 3.0);
                int theoreticalPeaksTop = spectrumImage.getHeight() - theoreticalPeaksHeight - 10;
                g.drawString("Ideal Peaks", indent, theoreticalPeaksTop - 20);
                //combined theoretical peak distribution chart
                PanelWithPeakChart theoreticalPeakChart = buildTheoreticalPeakChart(lightMz, heavyMz, charge, ratio,
                        sidebarWidth, (int) (sidebarWidth * 2.0 / 3.0));
                g.drawImage(theoreticalPeakChart.createImage(sidebarWidth,
                        (int) (sidebarWidth * 2.0 / 3.0)), 0, theoreticalPeaksTop, null);
            }
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

    public boolean isWriteTheoreticalPeaksOnCharts()
    {
        return writeTheoreticalPeaksOnCharts;
    }

    public void setWriteTheoreticalPeaksOnCharts(boolean writeTheoreticalPeaksOnCharts)
    {
        this.writeTheoreticalPeaksOnCharts = writeTheoreticalPeaksOnCharts;
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

    public File getOutTurkFile()
    {
        return outTurkFile;
    }

    public void setOutTurkFile(File outTurkFile)
    {
        this.outTurkFile = outTurkFile;
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

    public boolean isAppendTsvOutput()
    {
        return appendTsvOutput;
    }

    public void setAppendTsvOutput(boolean appendTsvOutput)
    {
        this.appendTsvOutput = appendTsvOutput;
    }

    public void addProgressListener(ActionListener listener)
    {
        dummyProgressButton.addActionListener(listener);
    }

    public boolean isShouldCreateCharts()
    {
        return shouldCreateCharts;
    }

    public void setShouldCreateCharts(boolean shouldCreateCharts)
    {
        this.shouldCreateCharts = shouldCreateCharts;
    }

    public boolean isMarkAllEventsBad()
    {
        return markAllEventsBad;
    }

    public void setMarkAllEventsBad(boolean markAllEventsBad)
    {
        this.markAllEventsBad = markAllEventsBad;
    }

    public String getTurkImageURLPrefix()
    {
        return turkImageURLPrefix;
    }

    public void setTurkImageURLPrefix(String turkImageURLPrefix)
    {
        this.turkImageURLPrefix = turkImageURLPrefix;
    }
}
