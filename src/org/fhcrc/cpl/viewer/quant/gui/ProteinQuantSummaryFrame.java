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
package org.fhcrc.cpl.viewer.quant.gui;

import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.QuantitationUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.widget.SwingWorkerWithProgressBarDialog;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.XYPlot;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.imageio.ImageIO;
import java.io.File;
import java.io.StringWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.awt.*;
import java.awt.event.*;


/**
 * A dialog box with a table of quantitative events for a single protein, or a few.  Controls for viewing
 * event details and for generating charts from a list of events.  The latter action closes the dialog
 * and adds the events to the events in the main Qurate window
 */
public class ProteinQuantSummaryFrame extends JDialog
{
    protected static Logger _log = Logger.getLogger(ProteinQuantSummaryFrame.class);

    protected int width = 800;
    protected int height = 900;

    protected float maxChartDisplayLogRatio = (float) Math.log(20f);
    protected float minChartDisplayLogRatio = (float) Math.log(1f/20f);

    protected int LOGRATIO_HISTOGRAM_PANEL_HEIGHT = 150;

    //min and max ratio for display in table.  Ratio must be > min OR < max, or both
    protected float minHighRatio = 0f;
    protected float maxLowRatio = 999f;


    //protein-level info
    protected List<String> proteinNames;
    protected List<Float> proteinRatio;

    //Quantitative events
    protected List<QuantEvent> quantEvents;
    protected List<QuantEvent> selectedQuantEvents;

    protected String labeledResidue = null;
    protected float labelMassDiff = 0f;

    protected int labelType = -1;

    //needed for chart generation
    protected File mzXmlDir;

    //Output location
    protected File outDir;
    protected File outFile;
    //If the output file exists already, do we append to it?
    protected boolean appendOutput = true;

    //dimensions
    protected int fullWidth = 1000;
    protected int fullHeight = 500;
    protected int summaryPanelHeight = 100;
    protected int propertiesWidth = 400;
    protected int propertiesHeight = 500;

    //This is hacky.  It's for sizing the window appropriately when we know the number of table rows.  There
    //is probably a much more reasonable way to do this in AWT.
    //TODO: do this with less tacky
    protected final int TITLEBAR_HEIGHT = 55;
    protected final int STATUSPANEL_HEIGHT = 25;
    protected final int SUMMARYPANEL_HEIGHT = 81;
    protected final int TABLEROW_HEIGHT = 17;

    //container declarations
    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public JScrollPane eventsScrollPane;
    public JPanel eventsPanel;
    public JPanel logRatioHistogramPanel;

    protected PanelWithHistogram logRatioHistogram;

    protected Map<String, List<String>> proteinGenesMap;

    //single event details components
    protected QuantEvent.QuantEventPropertiesTable eventPropertiesTable;
    protected JDialog eventPropertiesDialog;

    //Should we roll in events that overlap the selected events?
    protected boolean shouldAddOverlappingEvents = true;

    JTextArea proteinNameTextArea;
    JTextArea proteinRatioTextArea;
    JButton buildChartsForSelectedButton = new JButton("Build Selected Charts");
    JButton showPropertiesButton = new JButton("Show Event Properties");

    //an array of quantitation events that are already built.  These will be shown as already selected,
    //and un-unselectable, and will not be rebuilt.
    protected List<QuantEvent> existingQuantEvents;


    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    // ratio filter control
    public JLabel maxLowRatioLabel;
    public JLabel minHighRatioLabel;
    public JLabel numPassingEventsLabel;



    //event properties
    protected QuantEventsSummaryTable eventsTable;

    public ProteinQuantSummaryFrame()
    {
        initGUI();
    }

    public ProteinQuantSummaryFrame(File outDir, File mzXmlDir, File outFile, boolean appendOutput)
            throws IllegalArgumentException
    {
        this();
        this.outDir = outDir;
        this.outFile = outFile;
        this.mzXmlDir = mzXmlDir;
        this.appendOutput = appendOutput;
    }

    /**
     * Find all the peptides contributing to the ratio for the FIRST OCCURRENCE of a
     * protein in the pepXML file, find all quantitative events for those peptides in the
     * pepXML file, and show them
     * @param pepXmlFile
     * @param proteins
     */
    public void displayData(File pepXmlFile, List<ProtXmlReader.Protein> proteins)
    {
        Collections.sort(proteins, new Comparator<ProtXmlReader.Protein>()
        {
            public int compare(ProtXmlReader.Protein o1, ProtXmlReader.Protein o2)
            {             
                return o1.getProteinName().compareTo(o2.getProteinName());
            }
        });

        StringBuffer proteinNameLabelTextBuf = new StringBuffer();
        StringBuffer ratioLabelTextBuf = new StringBuffer();

        this.proteinNames = new ArrayList<String>();
        List<ProtXmlReader.QuantitationRatio> quantRatios = new ArrayList<ProtXmlReader.QuantitationRatio>();
        for (int i=0; i<proteins.size(); i++)
        {
            String proteinName = proteins.get(i).getProteinName();
            if (i>0)
            {
                proteinNameLabelTextBuf.append(",");
                ratioLabelTextBuf.append(",");
                //5 proteins per line
                if (i%5==0)
                {
                    proteinNameLabelTextBuf.append("\n");
                    ratioLabelTextBuf.append("\n");
                }
            }
            proteinNameLabelTextBuf.append(proteinName);
            ProtXmlReader.QuantitationRatio quantRatio = proteins.get(i).getQuantitationRatio();
            quantRatios.add(quantRatio);
            ratioLabelTextBuf.append(quantRatio.getRatioMean());
            proteinNames.add(proteinName);
        }
        if (proteinNames.size() == 1)
            eventsTable.hideProteinColumn();
        if (proteinGenesMap == null)
            eventsTable.hideGeneColumn();
        else
        {
            eventsTable.setProteinGenesMap(proteinGenesMap);
System.err.println("*****SET!!! " + proteinGenesMap.size());
        }
System.err.println("@@@Should see set");

        proteinNameTextArea.setText(proteinNameLabelTextBuf.toString());
        proteinRatioTextArea.setText(ratioLabelTextBuf.toString());

        contentPanel.updateUI();


        quantEvents = new ArrayList<QuantEvent>();

        Map<String, Set<String>> peptideProteinsQuantMap = new HashMap<String, Set<String>>();
        for (int i=0; i<proteins.size(); i++)
        {
            for (String peptide : quantRatios.get(i).getPeptides())
            {
                Set<String> proteinsThisPep = peptideProteinsQuantMap.get(peptide);
                if (proteinsThisPep == null)
                {
                    proteinsThisPep = new HashSet<String>();
                    peptideProteinsQuantMap.put(peptide, proteinsThisPep);
                }
                proteinsThisPep.add(proteins.get(i).getProteinName());
            }
        }


        try
        {
            PepXMLFeatureFileHandler.PepXMLFeatureSetIterator fsi =
                    new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
            int numFractions = 0;
            while (fsi.hasNext())
            {
                boolean thisFracHasEvents = false;
                FeatureSet featureSet = fsi.next();
                setMessage("Checking fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                _log.debug("Checking fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                //check all features to see if they're in our list of peptides.  If so, add to quantEvents
                for (Feature feature : featureSet.getFeatures())
                {
                    String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                    if (peptideProteinsQuantMap.containsKey(peptide) &&
                            IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        thisFracHasEvents = true;
                        //pick up the labeled residue from the first feature
                        if (labeledResidue == null)
                        {
                            AnalyzeICAT.IsotopicLabel label = IsotopicLabelExtraInfoDef.getLabel(feature);
                            if (label != null)
                            {
                                labeledResidue = "" + label.getResidue();
                                labelMassDiff = label.getHeavy() - label.getLight();
                                _log.debug("Found label: " + labeledResidue + ", " + labelMassDiff);
                            }
                        }
                        QuantEvent quantEvent =
                                new QuantEvent(feature, MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                        quantEvent.setProtein(new ArrayList<String>(peptideProteinsQuantMap.get(peptide)).get(0));
                        quantEvents.add(quantEvent);
                    }
                }
                if (thisFracHasEvents)
                    numFractions++;
            }
            if (numFractions < 2)
                setMessage("Loaded all quantitation events from 1 fraction");
            else
                setMessage("Loaded all quantitation events from " + numFractions + " separate fractions");
            if (labeledResidue == null)
                infoMessage("WARNING: unable to determine modification used for quantitation.  " +
                        "Cannot collapse light and heavy states or perform assessment.");
            else
            {
                labelType = QuantitationUtilities.inferLabelType(labeledResidue, labelMassDiff);
            }
        }
        catch (Exception e)
        {
            errorMessage("Failed to load features from pepXML file",e);
            return;
        }

        //sort by peptide, then fraction, then charge, then modifications
        Collections.sort(quantEvents,
                new QuantEvent.ProteinPeptideFractionChargeModificationsRatioAscComparator());
        displayEvents();
        if (quantRatios.size() == 1)
        {
            eventsTable.setLogRatioHeaderRatio(quantRatios.get(0).getRatioMean());
        }

        updateExtremeRatioGUI();
        
    }


    /**
     * Initialize the GUI components
     */
    protected void initGUI()
    {
        //Global stuff
        setSize(fullWidth, fullHeight);

        eventPropertiesTable = new QuantEvent.QuantEventPropertiesTable();
        eventPropertiesTable.setVisible(true);
        JScrollPane eventPropsScrollPane = new JScrollPane();
        eventPropsScrollPane.setViewportView(eventPropertiesTable);
        eventPropsScrollPane.setSize(propertiesWidth, propertiesHeight);
        eventPropertiesDialog = new JDialog(this, "Event Properties");
        eventPropertiesDialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        eventPropertiesDialog.setSize(propertiesWidth, propertiesHeight);
        eventPropertiesDialog.setContentPane(eventPropsScrollPane);

        ListenerHelper helper = new ListenerHelper(this);
        setTitle("Protein Summary");
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(5,5,5,5);
        gbc.weighty = 1;
        gbc.weightx = 1;

        try
        {
            (getOwner()).setIconImage(ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif")));
        }
        catch (Exception e)
        {
        }

        try
        {
            Localizer.renderSwixml("org/fhcrc/cpl/viewer/quant/gui/ProteinQuantSummaryFrame.xml",this);
            assert null != contentPanel;
            setContentPane(contentPanel);            
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("error creating dialog", x);
            throw new RuntimeException(x);
        }

        //status message
        messageLabel.setBackground(Color.WHITE);
        messageLabel.setFont(Font.decode("verdana plain 12"));
        messageLabel.setText(" ");

        buildChartsForSelectedButton.setEnabled(false);
        helper.addListener(buildChartsForSelectedButton, "buttonBuildCharts_actionPerformed");
        showPropertiesButton.setEnabled(false);
        helper.addListener(showPropertiesButton, "buttonShowProperties_actionPerformed");

        //summary panel
        summaryPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
        summaryPanel.setPreferredSize(new Dimension(fullWidth, summaryPanelHeight));
        summaryPanel.setMinimumSize(new Dimension(200, summaryPanelHeight));
        gbc.fill = GridBagConstraints.NONE;

        JLabel proteinNameLabel = new JLabel("Protein(s):");
        gbc.gridwidth = 1;
        summaryPanel.add(proteinNameLabel, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        proteinNameTextArea = new JTextArea("");
        proteinNameTextArea.setEditable(false);
        JScrollPane proteinNamesScrollPane = new JScrollPane();
        proteinNamesScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        proteinNamesScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        proteinNamesScrollPane.setViewportView(proteinNameTextArea);
        proteinNamesScrollPane.setMinimumSize(new Dimension(500, 100));
        summaryPanel.add(proteinNamesScrollPane, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(buildChartsForSelectedButton, gbc);
        gbc.gridwidth = 1;
        JLabel proteinRatioLabel = new JLabel("Ratio(s):");
        summaryPanel.add(proteinRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        proteinRatioTextArea = new JTextArea("");
        proteinRatioTextArea.setEditable(false);
        JScrollPane proteinRatiosScrollPane = new JScrollPane();
        proteinRatiosScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        proteinRatiosScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        proteinRatiosScrollPane.setViewportView(proteinRatioTextArea);
        proteinRatiosScrollPane.setMinimumSize(new Dimension(500, 100));

        summaryPanel.add(proteinRatiosScrollPane, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(showPropertiesButton, gbc);

        gbc.fill = GridBagConstraints.BOTH;

        eventsScrollPane = new JScrollPane();
        eventsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        eventsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        eventsPanel = new JPanel();
        eventsPanel.setLayout(new GridBagLayout());

        eventsTable = new QuantEventsSummaryTable();
        ListSelectionModel tableSelectionModel = eventsTable.getSelectionModel();
        tableSelectionModel.addListSelectionListener(new EventsTableListSelectionHandler());
        eventsScrollPane.setViewportView(eventsTable);
        eventsScrollPane.setMinimumSize(new Dimension(400, 400));
                                                                    

        gbc.insets = new Insets(0,0,0,0);             
        mainPanel.add(eventsScrollPane, gbc);

        logRatioHistogramPanel.setBorder(BorderFactory.createTitledBorder("Log Ratios"));
        maxLowRatioLabel = new JLabel("Max Low Ratio: ");
        minHighRatioLabel = new JLabel("Min High Ratio: ");
        numPassingEventsLabel = new JLabel("Events Retained: ");
        gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth = 1;
        logRatioHistogramPanel.add(maxLowRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        logRatioHistogramPanel.add(minHighRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        logRatioHistogramPanel.add(numPassingEventsLabel, gbc);
    }



    /**
     * Chart-building is long-running and needs to provide user feedback on status, so it
     * runs in a SwingWorker that displays a progress bar.
     *
     * This disposes the parent (ProteinQuantSummaryFrame) at the end, which is a bit goofy,
     * but that's the signal for the charts that we build here to be added to the QuantitationReviewer.
     * Would probably be nice to send that signal some other way.
     */
    protected class ChartBuilderWorker extends
            SwingWorkerWithProgressBarDialog<Throwable, String>
    {
        protected QuantitationVisualizer quantVisualizer;
        protected static final String expressionForLabel = "Processed " +
                SwingWorkerWithProgressBarDialog.CURRENT_VALUE_TOKEN + " of " +
                SwingWorkerWithProgressBarDialog.MAX_VALUE_TOKEN + " events";

        ChartBuilderWorker(JDialog parent, QuantitationVisualizer quantVisualizer)
        {
            super(parent, 0, selectedQuantEvents.size(), 0, expressionForLabel, "Building Charts...");
            this.quantVisualizer = quantVisualizer;
            this.quantVisualizer.setLabelType(labelType);
            quantVisualizer.addProgressListener(new ProgressBarUpdater(progressBar));
        }

        protected class ProgressBarUpdater implements ActionListener
        {
            protected JProgressBar progressBar;

            public ProgressBarUpdater(JProgressBar progressBar)
            {
                this.progressBar = progressBar;
            }

            public void actionPerformed(ActionEvent event)
            {
                int numProcessed = Integer.parseInt(event.getActionCommand());
                updateLabelText(numProcessed);
                progressBar.setValue(numProcessed);
            }
        }

        public Throwable doInBackground()
        {
            try
            {
                quantVisualizer.visualizeQuantEvents(selectedQuantEvents, true);
            }
            catch (IOException e)
            {
                ApplicationContext.errorMessage(e.getMessage(),e);
                return e;
            }
            finally
            {
                if (progressDialog != null) progressDialog.dispose();
            }
            return null;
        }

        protected void done()
        {
            try
            {
                Throwable throwable = get();
                if (throwable == null)
                {
                    infoMessage("Saved chart summary to file " + quantVisualizer.getOutTsvFile().getAbsolutePath());
                    parent.dispose();
                }
                else
                {
                    errorMessage("Error building charts", throwable);
                }
            }
            catch (ExecutionException e)
            {
                errorMessage("Error building charts",e);
            }
            catch (InterruptedException e)
            {

            }
        }
    }

    public void buttonShowProperties_actionPerformed(ActionEvent event)
    {
        eventPropertiesDialog.setVisible(true);
    }

    /**
     * Build charts (in a separate worker thread) and dispose()
     * @param event
     */
    public void buttonBuildCharts_actionPerformed(ActionEvent event)
    {
        selectedQuantEvents = eventsTable.getSelectedEvents();
        if (selectedQuantEvents.isEmpty())
        {
            infoMessage("No new events selected");
            return;
        }

        if (shouldAddOverlappingEvents)
        {
            QuantitationVisualizer quantVisualizer = new QuantitationVisualizer();
            ApplicationContext.infoMessage("Finding all overlapping events, starting with " +
                    selectedQuantEvents.size());

            //This is not as efficient as it could be.  Rather than looking for overlapping events
            //just for our selected events, I find all sets of overlapping events and then look
            //for selected events in them.  Doing it in a more targeted way would get pretty complicated,
            //though, and it's not a big burden to check them all.
            List<QuantEvent> allOverlappingEvents =
                    quantVisualizer.findNonOverlappingQuantEventsAllPeptides(quantEvents,
                            labeledResidue, labelMassDiff);
            _log.debug("Got overlapping events, " + allOverlappingEvents.size());
            List<QuantEvent> eventsRepresentingSelectedAndOverlap =
                    new ArrayList<QuantEvent>();
            for (QuantEvent quantEvent : allOverlappingEvents)
            {
                if (selectedQuantEvents.contains(quantEvent))
                {
                    eventsRepresentingSelectedAndOverlap.add(quantEvent);
                    continue;
                }
                for (QuantEvent otherEvent : quantEvent.getOtherEvents())
                {
                    if (selectedQuantEvents.contains(otherEvent))
                    {
                        eventsRepresentingSelectedAndOverlap.add(quantEvent);
                    }
                }
            }
            selectedQuantEvents = eventsRepresentingSelectedAndOverlap;
            ApplicationContext.infoMessage("Including overlapping events, selected events: " +
                    selectedQuantEvents.size());
        }

        ApplicationContext.infoMessage(selectedQuantEvents.size() + " events selected for charts");



        setMessage("Building charts for " + selectedQuantEvents.size() + " events...");
        QuantitationVisualizer quantVisualizer = new QuantitationVisualizer();
        quantVisualizer.setMzXmlDir(mzXmlDir);
        String allProteinNamesUnderscores = concatProteinNames("_");
        File proteinOutDir = new File(outDir, allProteinNamesUnderscores);
        proteinOutDir.mkdir();
        quantVisualizer.setOutDir(outDir);
        File outputFile = outFile;
        if (outputFile == null)
            outputFile = (new File(outDir, "quantitation_" + allProteinNamesUnderscores + ".tsv"));
        quantVisualizer.setOutTsvFile(outputFile);
        quantVisualizer.setOutHtmlFile(new File(outDir, "quantitation_" +allProteinNamesUnderscores + ".html"));
        quantVisualizer.setAppendTsvOutput(appendOutput);

        ChartBuilderWorker swingWorker = new ChartBuilderWorker(this, quantVisualizer);
        swingWorker.execute();
    }

    /**
     * Shut down any stray children and die
     */
    public void dispose()
    {
        if (eventPropertiesDialog != null)
            eventPropertiesDialog.dispose();
        super.dispose();
    }

    /**
     * Concatenate all values in proteinNames, separated by separatorString
     * @param separatorString
     * @return
     */
    protected String concatProteinNames(String separatorString)
    {
        StringBuffer allProteinNamesStringBuf = new StringBuffer();

        for (int i=0; i<proteinNames.size(); i++)
        {
            if (i > 0)
                allProteinNamesStringBuf.append(separatorString);
            allProteinNamesStringBuf.append(proteinNames.get(i));
        }
        return allProteinNamesStringBuf.toString();
    }

    /**
     * Populate the table with the current quantEvents
     */
    protected void displayEvents()
    {
        String proteinOrProteins = "Protein";
        if (proteinNames.size() > 1)
            proteinOrProteins = "Proteins";
        setTitle("Event Summary for " + proteinOrProteins + " " + concatProteinNames(","));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 1;
        gbc.weightx = 1;

        List<Integer> alreadySelectedEventIndices = new ArrayList<Integer>();
        if (existingQuantEvents != null)
        {
            for (int i=0; i<quantEvents.size(); i++)
            {
                for (QuantEvent existingEvent : existingQuantEvents)
                {

                    if (quantEvents.get(i).isSameEvent(existingEvent))
                    {
                        alreadySelectedEventIndices.add(i);
                        quantEvents.set(i, existingEvent);

//                        quantEvents.get(i).setQuantCurationStatus(existingEvent.getQuantCurationStatus());
//                        quantEvents.get(i).setIdCurationStatus(existingEvent.getIdCurationStatus());

                        break;
                    }
                }
            }
        }

        List<Float> eventLogRatios = new ArrayList<Float>();

        for (QuantEvent event : quantEvents)
        {
            float boundedProteinLogRatio = (float) Math.min(maxChartDisplayLogRatio,
                    Math.max(minChartDisplayLogRatio, Math.log(event.getRatio())));
            eventLogRatios.add(boundedProteinLogRatio);
        }

        eventsTable.displayEvents(quantEvents, alreadySelectedEventIndices);

        buildChartsForSelectedButton.setEnabled(true);
        showPropertiesButton.setEnabled(true);


        logRatioHistogram = new PanelWithHistogram(eventLogRatios, "Protein Log Ratios", 200);

        LogRatioHistMouseListener histMouseListener =
                new LogRatioHistMouseListener(logRatioHistogram);
        histMouseListener.addRangeUpdateListener(
                new LogRatioHistogramListener(this, histMouseListener));
        ChartPanel histChartPanel = logRatioHistogram.getChartPanel();
        histChartPanel.removeMouseListener(histChartPanel);
        histChartPanel.removeMouseMotionListener(histChartPanel);
        logRatioHistogram.getChartPanel().addMouseListener(histMouseListener);
        logRatioHistogram.getChartPanel().addMouseMotionListener(histMouseListener);
        mainPanel.updateUI();

        int chartWidth = width - 25;
        int chartHeight = LOGRATIO_HISTOGRAM_PANEL_HEIGHT-70;
        Dimension histDimension =  new Dimension(chartWidth, chartHeight);
        logRatioHistogram.setPreferredSize(histDimension);
        logRatioHistogram.setSize(histDimension);
        logRatioHistogram.getChart().removeLegend();

        //remove axes from chart
        ((XYPlot)logRatioHistogram.getPlot()).getDomainAxis().setVisible(false);
        ((XYPlot)logRatioHistogram.getPlot()).getRangeAxis().setVisible(false);
//System.err.println("**" + chartWidth + ","+chartHeight + "... " + logRatioHistogramPanel.getWidth() + ","+logRatioHistogramPanel.getHeight());
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        logRatioHistogramPanel.add(logRatioHistogram, gbc);

        contentPanel.updateUI();

        fullHeight = Math.min(800, (quantEvents.size() + 1) * TABLEROW_HEIGHT + SUMMARYPANEL_HEIGHT +
                LOGRATIO_HISTOGRAM_PANEL_HEIGHT + STATUSPANEL_HEIGHT + TITLEBAR_HEIGHT);                                
        setSize(fullWidth, fullHeight);        
    }

    protected class LogRatioHistogramListener implements ActionListener
    {
        protected ProteinQuantSummaryFrame proteinSummaryFrame;
        protected LogRatioHistMouseListener logRatioHistMouseListener;

        public LogRatioHistogramListener(ProteinQuantSummaryFrame proteinSummaryFrame,
                                         LogRatioHistMouseListener logRatioHistMouseListener)
        {
            this.proteinSummaryFrame = proteinSummaryFrame;
            this.logRatioHistMouseListener = logRatioHistMouseListener;
        }
        public void actionPerformed(ActionEvent e)
        {
            proteinSummaryFrame.minHighRatio = (float) Math.exp(logRatioHistMouseListener.getSelectedXMaxValue());
            proteinSummaryFrame.maxLowRatio = (float) Math.exp(logRatioHistMouseListener.getSelectedXMinValue());
            proteinSummaryFrame.updateExtremeRatioGUI();
        }
    }

    protected void updateExtremeRatioGUI()
    {
        eventsTable.showOnlyExtremeRatios(maxLowRatio, minHighRatio);
        maxLowRatioLabel.setText("Max Low Ratio: " + Rounder.round(maxLowRatio, 2));
        minHighRatioLabel.setText("Min High Ratio: " + Rounder.round(minHighRatio, 2));
        numPassingEventsLabel.setText("Events Retained: " + eventsTable.getRowCount() + " / " +
                quantEvents.size() + " (" + Rounder.round(100f * eventsTable.getRowCount() / quantEvents.size(),1) + "%)");
    }

    /**
     * Display a dialog box with info message
     * @param message
     */
    protected void infoMessage(String message)
    {
        if (message.length() > 2000)
            message = message.substring(0,1997) + "...";
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Display a dialog box with info message and stack trace
     * @param message
     * @param t
     */
    protected void errorMessage(String message, Throwable t)
    {
        ApplicationContext.errorMessage(message, t);
        if (null != t)
        {
            message = message + "\n" + t.getMessage() + "\n";

            StringWriter sw = new StringWriter();
            PrintWriter w = new PrintWriter(sw);
            t.printStackTrace(w);
            w.flush();
            message += "\n";
            message += sw.toString();
        }
        infoMessage(message);
    }


    /**
     * Set status message.  Separate thread necessary or UI hangs
     * @param message
     */
    public void setMessage(String message)
    {
        if (EventQueue.isDispatchThread())
        {
            if (null == message || 0 == message.length())
                message = " ";
            if (message.length() > 500)
                message = message.substring(0,500);
            messageLabel.setText(message);
        }
        else
        {
            final String msg = message;
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    setMessage(msg);
                }
            });
        }
        statusPanel.updateUI();
    }

    public List<QuantEvent> getSelectedQuantEvents()
    {
        return selectedQuantEvents;
    }

    public File getOutFile()
    {
        return outFile;
    }

    public void setOutFile(File outFile)
    {
        this.outFile = outFile;
    }

    /**
     * display the properties for the selected event, if only one's selected
     */
    public class EventsTableListSelectionHandler implements ListSelectionListener
    {
        public void valueChanged(ListSelectionEvent e)
        {
            if (!e.getValueIsAdjusting())
            {
                int selectedIndex = eventsTable.getSelectedIndex();
                if (selectedIndex >= 0)
                    eventPropertiesTable.displayQuantEvent(quantEvents.get(selectedIndex));
                else
                    eventPropertiesTable.clearProperties();
            }
        }
    }

    public List<QuantEvent> getExistingQuantEvents()
    {
        return existingQuantEvents;
    }

    public void setExistingQuantEvents(List<QuantEvent> existingQuantEvents)
    {
        this.existingQuantEvents = existingQuantEvents;
    }

    public Map<String, List<String>> getProteinGeneMap()
    {
        return proteinGenesMap;
    }

    public void setProteinGeneMap(Map<String, List<String>> proteinGeneMap)
    {
        this.proteinGenesMap = proteinGeneMap;
    }
}
