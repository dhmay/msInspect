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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.widget.SwingWorkerWithProgressBarDialog;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
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
import java.awt.geom.Ellipse2D;
import java.awt.event.*;


/**
 * A dialog box with a table of quantitative events for a single protein, or a few.  Controls for viewing
 * event details and for generating charts from a list of events.  The latter action closes the dialog
 * and adds the events to the events in the main Qurate window
 */
public class ProteinQuantSummaryFrame extends JDialog
{
    protected static Logger _log = Logger.getLogger(ProteinQuantSummaryFrame.class);

    protected int width = 900;
    protected int height = 900;


    protected int LOGRATIO_HISTOGRAM_PANEL_HEIGHT = 150;

    //min and max ratio for display in table.  Ratio must be > min OR < max, or both
    protected float minHighRatio = 0f;
    protected float maxLowRatio = 999f;

    protected int proteinDialogWidth = 500;
    protected int proteinDialogHeight = 800;

    Map<String, List<QuantEvent>> proteinEventsMap;
    protected String proteinTableSelectedProtein;


    //protein-level info
    protected List<String> proteinNames;
    protected List<Float> proteinRatio;

    //Quantitative events
    protected List<QuantEvent> quantEvents;
    protected List<QuantEvent> selectedQuantEvents;

    //SILAC or Acrylamide
    protected int labelType = -1;
    protected String labeledResidue = null;
    protected float labelMassDiff = 0f;

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
    protected final int SUMMARYPANEL_HEIGHT = 41;
    protected final int TABLEROW_HEIGHT = 17;

    protected final int PROTEINTABLE_HISTPANEL_HEIGHT = 150;
    protected final int PROTEINTABLE_SCATTERPLOTPANEL_HEIGHT = 100;


    //container declarations
    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public JScrollPane eventsScrollPane;
    public JPanel eventsPanel;

    public PanelWithLogRatioHistAndFields logRatioHistogramPanel;

    public PanelWithLogRatioHistAndFields perProteinLogRatioHistogramPanel;

    public JPanel perProteinPeptideLogRatioPanel;



    protected Map<String, List<String>> proteinGenesMap;

    //single event details components
    protected QuantEvent.QuantEventPropertiesTable eventPropertiesTable;
    protected JDialog eventPropertiesDialog;

    //Should we roll in events that overlap the selected events?
    protected boolean shouldAddOverlappingEvents = true;

    JTextArea proteinNameTextArea;
    JTextArea proteinRatioTextArea;
    JButton showProteinRatiosButton = new JButton("Show Protein Table");

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

    protected JTable proteinRatiosTable;
    protected JDialog proteinRatiosDialog;


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



        buildChartsForSelectedButton.setEnabled(false);
        helper.addListener(buildChartsForSelectedButton, "buttonBuildCharts_actionPerformed");
        showPropertiesButton.setEnabled(false);
        helper.addListener(showPropertiesButton, "buttonShowProperties_actionPerformed");
        showProteinRatiosButton.setEnabled(false);
        helper.addListener(showProteinRatiosButton, "buttonShowProteinRatios_actionPerformed");


        //summary panel
        summaryPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
        summaryPanel.setPreferredSize(new Dimension(fullWidth, summaryPanelHeight));
        summaryPanel.setMinimumSize(new Dimension(200, summaryPanelHeight));
        gbc.fill = GridBagConstraints.NONE;

        gbc.gridwidth = 1;
        summaryPanel.add(showProteinRatiosButton, gbc);

        gbc.gridwidth = GridBagConstraints.RELATIVE;
        summaryPanel.add(buildChartsForSelectedButton, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(showPropertiesButton, gbc);

        gbc.fill = GridBagConstraints.BOTH;

        eventsScrollPane = new JScrollPane();
        eventsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        eventsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        eventsPanel = new JPanel();
        eventsPanel.setLayout(new GridBagLayout());

        eventsTable = new QuantEventsSummaryTable();
        eventsTable.getSelectionModel().addListSelectionListener(new EventsTableListSelectionHandler());
        eventsScrollPane.setViewportView(eventsTable);
        eventsScrollPane.setMinimumSize(new Dimension(400, 400));

        gbc.insets = new Insets(0,0,0,0);
        mainPanel.add(eventsScrollPane, gbc);

        logRatioHistogramPanel = new PanelWithLogRatioHistAndFields();
        logRatioHistogramPanel.setBorder(BorderFactory.createTitledBorder("Log Ratios"));
        logRatioHistogramPanel.setPreferredSize(new Dimension(width-10, 300));
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty=100;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        add(logRatioHistogramPanel, gbc);

        //status message
        messageLabel = new JLabel();
        messageLabel.setBackground(Color.WHITE);
        messageLabel.setFont(Font.decode("verdana plain 12"));
        messageLabel.setText(" ");
        statusPanel = new JPanel();
        gbc.weighty=1;
        statusPanel.setPreferredSize(new Dimension(width-10, 50));
        statusPanel.add(messageLabel, gbc);

        //per-protein event summary table; disembodied
        //todo: move this into its own class? it's getting kind of complicated
        proteinRatiosTable = new JTable();
        proteinRatiosTable.setVisible(true);
        ListSelectionModel proteinTableSelectionModel = proteinRatiosTable.getSelectionModel();
        proteinTableSelectionModel.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        proteinTableSelectionModel.addListSelectionListener(new ProteinTableListSelectionHandler());
        JScrollPane proteinRatiosScrollPane = new JScrollPane();
        proteinRatiosScrollPane.setViewportView(proteinRatiosTable);
        proteinRatiosScrollPane.setPreferredSize(new Dimension(proteinDialogWidth,
                proteinDialogHeight - PROTEINTABLE_HISTPANEL_HEIGHT - PROTEINTABLE_SCATTERPLOTPANEL_HEIGHT - 70));
        proteinRatiosDialog = new JDialog(this, "Protein Ratios");
        proteinRatiosDialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        proteinRatiosDialog.setSize(proteinDialogWidth, proteinDialogHeight);
        JPanel proteinRatiosContentPanel = new JPanel();
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.fill = GridBagConstraints.BOTH;        
        proteinRatiosContentPanel.add(proteinRatiosScrollPane, gbc);
        proteinRatiosDialog.setContentPane(proteinRatiosContentPanel);
        perProteinLogRatioHistogramPanel = new PanelWithLogRatioHistAndFields();
        perProteinLogRatioHistogramPanel.addRangeUpdateListener(new ProteinTableLogRatioHistogramListener());

        perProteinLogRatioHistogramPanel.setBorder(BorderFactory.createTitledBorder("Log Ratios"));
        perProteinLogRatioHistogramPanel.setPreferredSize(new Dimension(proteinDialogWidth-10,
                PROTEINTABLE_HISTPANEL_HEIGHT));
        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        proteinRatiosDialog.add(perProteinLogRatioHistogramPanel, gbc);

        perProteinPeptideLogRatioPanel = new JPanel();
        perProteinPeptideLogRatioPanel.setBorder(BorderFactory.createTitledBorder("By Peptide"));
        perProteinPeptideLogRatioPanel.setPreferredSize(new Dimension(proteinDialogWidth-10,
                PROTEINTABLE_SCATTERPLOTPANEL_HEIGHT));
        proteinRatiosDialog.add(perProteinPeptideLogRatioPanel, gbc);
    }    

    public void buttonShowProteinRatios_actionPerformed(ActionEvent event)
    {
        proteinRatiosDialog.setVisible(true);
    }

    /**
     * display the ratios for the selected event
     */
    public class ProteinTableListSelectionHandler implements ListSelectionListener
    {
        public void valueChanged(ListSelectionEvent e)
        {
            if (!e.getValueIsAdjusting())
            {
                ListSelectionModel lsm = proteinRatiosTable.getSelectionModel();
                if (lsm.isSelectionEmpty() || (lsm.getMinSelectionIndex() != lsm.getMaxSelectionIndex()) ||
                    lsm.getMinSelectionIndex() < 0)
                    return;
                proteinTableSelectedProtein = (String) proteinRatiosTable.getValueAt(lsm.getMinSelectionIndex(), 0);

                List<QuantEvent> proteinEvents = proteinEventsMap.get(proteinTableSelectedProtein);
                List<Float> eventLogRatios = new ArrayList<Float>();

                Map<String, List<Float>> perPeptideLogRatios = new HashMap<String, List<Float>>();
                for (QuantEvent event : proteinEvents)
                {
                    float logRatio = (float) Math.log(event.getRatio());
                    eventLogRatios.add(logRatio);
                    List<Float> peptideLogRatios = perPeptideLogRatios.get(event.getPeptide());
                    if (peptideLogRatios == null)
                    {
                        peptideLogRatios = new ArrayList<Float>();
                        perPeptideLogRatios.put(event.getPeptide(), peptideLogRatios);
                    }
                    peptideLogRatios.add(logRatio);
                }
                //Add a crosshair at the value of the protein log ratio
                perProteinLogRatioHistogramPanel.setDomainCrosshairValue((float) Math.log(
                        (Double) proteinRatiosTable.getValueAt(lsm.getMinSelectionIndex(), 1)));
                perProteinLogRatioHistogramPanel.setLogRatios(eventLogRatios);
                perProteinLogRatioHistogramPanel.setSize(proteinDialogWidth-10, PROTEINTABLE_HISTPANEL_HEIGHT);
                perProteinLogRatioHistogramPanel.setMaxLowRatio(1000);
                perProteinLogRatioHistogramPanel.setMinHighRatio(0.0001f);

                //Scatterplot of per-peptide log-ratios
                perProteinPeptideLogRatioPanel.removeAll();
                PanelWithScatterPlot peptideProteinLogRatioScatterPlot = new PanelWithScatterPlot();
                int peptideIndex = 0;
                for (List<Float> peptideLogRatios : perPeptideLogRatios.values())
                {
                    List<Float> dummyYValues = new ArrayList<Float>(peptideLogRatios.size());
                    for (int i=0; i<peptideLogRatios.size(); i++)
                        dummyYValues.add((float)peptideIndex);
                    peptideProteinLogRatioScatterPlot.addData(peptideLogRatios, dummyYValues, "asdf");
                    peptideIndex++;
                }
                //remove axes from chart
                ((XYPlot)peptideProteinLogRatioScatterPlot.getPlot()).getRangeAxis().setVisible(false);
                ((XYPlot)peptideProteinLogRatioScatterPlot.getPlot()).getDomainAxis().setVisible(false);

                peptideProteinLogRatioScatterPlot.setPointSize(6);
                peptideProteinLogRatioScatterPlot.getChart().removeLegend();
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().setDomainGridlinesVisible(false);
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().setRangeGridlinesVisible(false);
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().getDomainAxis().setLowerBound(
                        Math.log(perProteinLogRatioHistogramPanel.getMinRatioBound()));
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().getDomainAxis().setUpperBound(
                        Math.log(perProteinLogRatioHistogramPanel.getMaxRatioBound()));
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().getRangeAxis().setLowerBound(-1);
                peptideProteinLogRatioScatterPlot.getChart().getXYPlot().getRangeAxis().setUpperBound(
                        perPeptideLogRatios.size()+1);

                peptideProteinLogRatioScatterPlot.setPreferredSize(
                        new Dimension(proteinDialogWidth-20, PROTEINTABLE_SCATTERPLOTPANEL_HEIGHT-30));
                peptideProteinLogRatioScatterPlot.updateUI();
                GridBagConstraints gbc = new GridBagConstraints();
                gbc.fill = GridBagConstraints.BOTH;
                gbc.anchor = GridBagConstraints.PAGE_START;
                gbc.gridwidth = GridBagConstraints.REMAINDER;
                gbc.insets = new Insets(5,5,5,5);
                gbc.weighty = 1;
                gbc.weightx = 1;
                perProteinPeptideLogRatioPanel.add(peptideProteinLogRatioScatterPlot, gbc);

            }
        }
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

        DefaultTableModel proteinRatiosTableModel = new DefaultTableModel(0, 4)
        {
            //all cells uneditable
            public boolean isCellEditable(int row, int column)
            {
                return false;
            }
            public Class getColumnClass(int columnIndex)
            {
                switch (columnIndex)
                {
                    case 0:return String.class;
                    case 1: return Float.class;
                    case 2: case 3: return Integer.class;
                }
                return String.class;
            }
        };
        proteinRatiosTable.setModel(proteinRatiosTableModel);
        proteinRatiosTable.getColumnModel().getColumn(0).setHeaderValue("Protein");
        proteinRatiosTable.getColumnModel().getColumn(1).setHeaderValue("Ratio");
        proteinRatiosTable.getColumnModel().getColumn(2).setHeaderValue("Quant Peptides");
        proteinRatiosTable.getColumnModel().getColumn(3).setHeaderValue("Events");


        this.proteinNames = new ArrayList<String>();
        List<ProtXmlReader.QuantitationRatio> quantRatios = new ArrayList<ProtXmlReader.QuantitationRatio>();
        proteinRatiosTableModel.setRowCount(proteins.size());            

        for (int i=0; i<proteins.size(); i++)
        {
            String proteinName = proteins.get(i).getProteinName();
            ProtXmlReader.QuantitationRatio quantRatio = proteins.get(i).getQuantitationRatio();
            quantRatios.add(quantRatio);
            proteinNames.add(proteinName);

            //careful -- the 3rd column values are populated below
            proteinRatiosTableModel.setValueAt(proteinName, i, 0);
            proteinRatiosTableModel.setValueAt(Rounder.round(quantRatio.getRatioMean(),2), i, 1);
        }
        if (proteinNames.size() == 1)
            eventsTable.hideProteinColumn();
        if (proteinGenesMap == null)
            eventsTable.hideGeneColumn();
        else
        {
            eventsTable.setProteinGenesMap(proteinGenesMap);
        }

        TableRowSorter<TableModel> sorter = new TableRowSorter<TableModel>(proteinRatiosTableModel);
        proteinRatiosTable.setRowSorter(sorter);        

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

        proteinEventsMap = new HashMap<String, List<QuantEvent>>();
        Map<String, Set<String>> proteinPeptidesMap = new HashMap<String, Set<String>>();

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
                        List<QuantEvent> eventsThisProtein = proteinEventsMap.get(quantEvent.getProtein());
                        if (eventsThisProtein == null)
                        {
                                eventsThisProtein = new ArrayList<QuantEvent>();
                            proteinEventsMap.put(quantEvent.getProtein(), eventsThisProtein);
                        }
                        eventsThisProtein.add(quantEvent);

                        Set<String> peptidesThisProtein = proteinPeptidesMap.get(quantEvent.getProtein());
                        if (peptidesThisProtein == null)
                        {
                            peptidesThisProtein = new HashSet<String>();
                            proteinPeptidesMap.put(quantEvent.getProtein(), peptidesThisProtein);
                        }
                        peptidesThisProtein.add(quantEvent.getPeptide());
                    }
                }
                if (thisFracHasEvents)
                    numFractions++;
            }

            for (int i=0; i<proteins.size(); i++)
            {
                proteinRatiosTableModel.setValueAt(proteinPeptidesMap.get(proteinNames.get(i)).size(), i, 2);
                proteinRatiosTableModel.setValueAt(proteinEventsMap.get(proteinNames.get(i)).size(), i, 3);
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
     * Chart-building is long-running and needs to provide user feedback on status, so it
     * runs in a SwingWorker that displays a progress bar.
     *
     * This disposes the parent (ProteinQuantSummaryFrame) at the end, which is a bit goofy,
     * but that's the signal for the charts that we build here to be added to the QuantitationReviewer.
     * Would probably be cleaner to send that signal some other way.
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
        String allProteinNamesUnderscores = concatProteinNamesMax3("_");
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
    protected String concatProteinNamesMax3(String separatorString)
    {
        StringBuffer allProteinNamesStringBuf = new StringBuffer();

        for (int i=0; i<proteinNames.size(); i++)
        {

            if (i > 0)
                allProteinNamesStringBuf.append(separatorString);
            if (i>2)
            {
                allProteinNamesStringBuf.append("etc");
                break;
            }
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
        if (proteinNames.size() > 4)
            setTitle("Event Summary for " + proteinNames.size() + " Proteins");
        else
            setTitle("Event Summary for " + proteinOrProteins + " " + concatProteinNamesMax3(","));

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
                    //check the events for identity with already-loaded events.  If we get it, grey out and
                    //replace the new one with the old one (to pick up curation statuses
                    if (quantEvents.get(i).isSameEvent(existingEvent))
                    {
                        alreadySelectedEventIndices.add(i);
                        quantEvents.set(i, existingEvent);

                        break;
                    }
                }
            }
        }

        List<Float> eventLogRatios = new ArrayList<Float>();

        for (QuantEvent event : quantEvents)
            eventLogRatios.add((float) Math.log(event.getRatio()));
        eventsTable.displayEvents(quantEvents, alreadySelectedEventIndices);

        buildChartsForSelectedButton.setEnabled(true);
        showPropertiesButton.setEnabled(true);
        showProteinRatiosButton.setEnabled(true);


        logRatioHistogramPanel.setMaxLowRatio(maxLowRatio);
        logRatioHistogramPanel.setMinHighRatio(minHighRatio);
        logRatioHistogramPanel.setLogRatios(eventLogRatios);
        logRatioHistogramPanel.setSize(width-5, LOGRATIO_HISTOGRAM_PANEL_HEIGHT-20);

        logRatioHistogramPanel.addRangeUpdateListener(new LogRatioHistogramListener());

        logRatioHistogramPanel.updateUI();

        contentPanel.updateUI();

        fullHeight = Math.min(800, (quantEvents.size() + 1) * TABLEROW_HEIGHT + SUMMARYPANEL_HEIGHT +
                LOGRATIO_HISTOGRAM_PANEL_HEIGHT + STATUSPANEL_HEIGHT + TITLEBAR_HEIGHT);                                
        setSize(fullWidth, fullHeight);        
    }

    /**
     * A chart listener that picks up events indicating changes to the selected area
     */
    protected class LogRatioHistogramListener implements ActionListener
    {
        public LogRatioHistogramListener()
        {
        }
        public void actionPerformed(ActionEvent e)
        {
            PanelWithLogRatioHistAndFields logRatioHistAndFields = (PanelWithLogRatioHistAndFields) e.getSource();
            minHighRatio = logRatioHistAndFields.getMinHighRatio();
            maxLowRatio = logRatioHistAndFields.getMaxLowRatio();
            updateExtremeRatioGUI();
        }
    }

        /**
     * A chart listener that picks up events indicating changes to the selected area on the protein table
     */
    protected class ProteinTableLogRatioHistogramListener implements ActionListener
    {

        public ProteinTableLogRatioHistogramListener()
        {
        }
        public void actionPerformed(ActionEvent e)
        {
            PanelWithLogRatioHistAndFields logRatioHistAndFields = (PanelWithLogRatioHistAndFields) e.getSource();
            List<QuantEvent> thisProteinEvents = proteinEventsMap.get(proteinTableSelectedProtein);
            float minRatio = logRatioHistAndFields.getMinHighRatio();
            float maxRatio = logRatioHistAndFields.getMaxLowRatio();
            for (int i=0; i<quantEvents.size(); i++)
            {
                QuantEvent quantEvent = quantEvents.get(i);
                if (thisProteinEvents.contains(quantEvent))
                {
                    boolean shouldSelect = false;
                    if (quantEvent.getRatio() >= minRatio || quantEvent.getRatio() <= maxRatio)
                        shouldSelect = true;
                    eventsTable.getModel().setValueAt(shouldSelect, i, 0);
                }
            }
        }
    }



    protected void updateExtremeRatioGUI()
    {
        eventsTable.showOnlyExtremeRatios(maxLowRatio, minHighRatio);
//        numPassingEventsLabel.setText("Events Retained: " + eventsTable.getRowCount() + " / " +
//                quantEvents.size() + " (" + Rounder.round(100f * eventsTable.getRowCount() / quantEvents.size(),1) +
//                "%)");
    }

    /**
     * Display a dialog box with info message
     * @param message
     */
    protected void infoMessage(String message)
    {
        if (message.length() > 2000)
            message = message.substring(0,1997) + "...";
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information",
                JOptionPane.INFORMATION_MESSAGE);
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
