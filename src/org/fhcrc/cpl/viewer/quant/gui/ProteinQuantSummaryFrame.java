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

import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.AnalyzeICAT;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.quant.QuantEventInfo;
import org.fhcrc.cpl.viewer.quant.QuantitationVisualizer;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.SwingUtils;
import org.fhcrc.cpl.toolbox.gui.widget.SwingWorkerWithProgressBarDialog;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.imageio.ImageIO;
import java.io.File;
import java.io.StringWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ExecutionException;
import java.awt.*;
import java.awt.event.*;


/**
 * A dialog box with a table of quantitative events for a single protein.  Controls for viewing
 * event details and for generating charts from a list of events.  The latter action closes the dialog
 * and adds the events to the events in the main Qurate window
 */
public class ProteinQuantSummaryFrame extends JDialog
{
    protected static Logger _log = Logger.getLogger(ProteinQuantSummaryFrame.class);


    //protein-level info
    protected String proteinName;
    protected float proteinRatio;

    //Quantitative events
    protected List<QuantEventInfo> quantEvents;
    protected List<QuantEventInfo> selectedQuantEvents;

    protected String labeledResidue = null;
    protected float labelMassDiff = 0f;

    //needed for chart generation
    protected File mzXmlDir;

    //Output location
    protected File outDir;
    protected File outFile;
    //If the output file exists already, do we append to it?
    protected boolean appendOutput = true;

    //List of row indexes in the table that are shaded.  This is for shading different peptides in
    //alternating colors.  It's distracting if the rows are resorted. Remove?
    protected List<Integer> shadedTableRows = new ArrayList<Integer>();

    //dimensions
    protected int fullWidth = 900;
    protected int fullHeight = 600;
    protected int propertiesWidth = 400;
    protected int propertiesHeight = 500;

    //container declarations
    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public JScrollPane eventsScrollPane;
    public JPanel eventsPanel;

    //single event details components
    protected QuantEventInfo.QuantEventPropertiesTable eventPropertiesTable;
    protected JDialog eventPropertiesDialog;

    //Should we roll in events that overlap the selected events?
    protected boolean shouldAddOverlappingEvents = true;

    JLabel proteinNameLabel = new JLabel("Protein: ");
    JLabel proteinRatioLabel  = new JLabel("Ratio: ");
    JButton buildChartsForSelectedButton = new JButton("Build Selected Charts");
    JButton showPropertiesButton = new JButton("Show Event Properties");

    //an array of quantitation events that are already built.  These will be shown as already selected,
    //and un-unselectable, and will not be rebuilt.
    protected List<QuantEventInfo> existingQuantEvents;


    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    //event properties
    protected QuantEventInfo.QuantEventsSummaryTable eventsTable;

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
     * @param protXmlFile
     * @param pepXmlFile
     * @param proteinName
     */
    public void displayData(File protXmlFile, File pepXmlFile, String proteinName)
    {
        this.proteinName = proteinName;

        proteinNameLabel.setText("Protein: " + proteinName);

        setMessage("Finding protein in file " + protXmlFile.getAbsolutePath());
        ProtXmlReader.Protein protein = null;
        try
        {
            protein = ProteinUtilities.loadFirstProteinOccurrence(protXmlFile, proteinName);
        }
        catch (Exception e)
        {
            errorMessage("Failure reading file " + protXmlFile.getAbsolutePath(),e);
            return;
        }

        if (protein == null || protein.getQuantitationRatio() == null)
        {
            throw new IllegalArgumentException("Protein " + proteinName + " does not occur in the file" +
                    " or is not quantitated");
        }

        ProtXmlReader.QuantitationRatio quantRatio = protein.getQuantitationRatio();

        proteinRatioLabel.setText("Ratio: " + quantRatio.getRatioMean());
        contentPanel.updateUI();

        //Now we've got the peptides that contributed to this ratio
        List<String> proteinPeptidesRatiosUsed = quantRatio.getPeptides();

        quantEvents = new ArrayList<QuantEventInfo>();
        try
        {
            PepXMLFeatureFileHandler.PepXMLFeatureSetIterator fsi =
                    new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
            while (fsi.hasNext())
            {
                FeatureSet featureSet = fsi.next();
                setMessage("Checking fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                _log.debug("Checking fraction " + MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                //check all features to see if they're in our list of peptides.  If so, add to quantEvents
                for (Feature feature : featureSet.getFeatures())
                {
                    if (proteinPeptidesRatiosUsed.contains(MS2ExtraInfoDef.getFirstPeptide(feature)) &&
                            IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
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
                        QuantEventInfo quantEvent =
                                new QuantEventInfo(feature, MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                        quantEvent.setProtein(this.proteinName);
                        quantEvents.add(quantEvent);
                    }
                }
            }
            setMessage("Loaded all quantitation events.");
            if (labeledResidue == null)
                infoMessage("WARNING: unable to determine modification used for quantitation.  " +
                        "Cannot collapse light and heavy states.");
        }
        catch (Exception e)
        {
            errorMessage("Failed to load features from pepXML file",e);
            return;
        }
        //sort by peptide, then fraction, then charge, then modifications
        Collections.sort(quantEvents,
                new QuantEventInfo.PeptideSequenceAscFractionAscChargeModificationsAscRatioAscComparator());
        displayEvents();
    }


    /**
     * Initialize the GUI components
     */
    protected void initGUI()
    {
        //Global stuff
        setSize(fullWidth, fullHeight);

        eventPropertiesTable = new QuantEventInfo.QuantEventPropertiesTable();
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
        gbc.fill = GridBagConstraints.NONE;

        gbc.gridwidth = GridBagConstraints.RELATIVE;
        summaryPanel.add(proteinNameLabel, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(buildChartsForSelectedButton, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;                   
        summaryPanel.add(proteinRatioLabel, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(showPropertiesButton, gbc);

        gbc.fill = GridBagConstraints.BOTH;

        eventsScrollPane = new JScrollPane();
        eventsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        eventsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        eventsPanel = new JPanel();
        eventsPanel.setLayout(new GridBagLayout());

        eventsTable = new QuantEventInfo.QuantEventsSummaryTable();
        ListSelectionModel tableSelectionModel = eventsTable.getSelectionModel();
        tableSelectionModel.addListSelectionListener(new EventsTableListSelectionHandler());
        eventsScrollPane.setViewportView(eventsTable);
        eventsScrollPane.setMinimumSize(new Dimension(400, 400));
        eventsTable.hideProteinColumn();

        gbc.insets = new Insets(0,0,0,0);             
        mainPanel.add(eventsScrollPane, gbc);  
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
            super(parent, 0, selectedQuantEvents.size(), 0, expressionForLabel);
            this.quantVisualizer = quantVisualizer;
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
                quantVisualizer.visualizeQuantEvents(selectedQuantEvents);
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
            List<QuantEventInfo> allOverlappingEvents =
                    quantVisualizer.findNonOverlappingQuantEventsAllPeptides(quantEvents,
                            labeledResidue, labelMassDiff);
            _log.debug("Got overlapping events, " + allOverlappingEvents.size());
            List<QuantEventInfo> eventsRepresentingSelectedAndOverlap =
                    new ArrayList<QuantEventInfo>();
            for (QuantEventInfo quantEvent : allOverlappingEvents)
            {
                if (selectedQuantEvents.contains(quantEvent))
                {
                    eventsRepresentingSelectedAndOverlap.add(quantEvent);
                    continue;
                }
                for (QuantEventInfo otherEvent : quantEvent.getOtherEvents())
                {
                    if (selectedQuantEvents.contains(otherEvent))
                    {
                        eventsRepresentingSelectedAndOverlap.add(quantEvent);
                        continue;
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
        File proteinOutDir = new File(outDir, proteinName);
        proteinOutDir.mkdir();
        quantVisualizer.setOutDir(proteinOutDir);
        File outputFile = outFile;
        if (outputFile == null)
            outputFile = (new File(outDir, "quantitation_" + proteinName + ".tsv"));
        quantVisualizer.setOutTsvFile(outputFile);
        quantVisualizer.setOutHtmlFile(new File(outDir, "quantitation_" + proteinName + ".html"));
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
     * Populate the table with the current quantEvents
     */
    protected void displayEvents()
    {
        setTitle("Protein Summary for " + proteinName);

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
                for (QuantEventInfo existingEvent : existingQuantEvents)
                {
                    if (quantEvents.get(i).isSameEvent(existingEvent))
                    {
                        alreadySelectedEventIndices.add(i);
                        break;
                    }
                }
            }
        }
        eventsTable.displayEvents(quantEvents, alreadySelectedEventIndices);

        buildChartsForSelectedButton.setEnabled(true);
        showPropertiesButton.setEnabled(true);


        contentPanel.updateUI();
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

    public List<QuantEventInfo> getSelectedQuantEvents()
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

    public List<QuantEventInfo> getExistingQuantEvents()
    {
        return existingQuantEvents;
    }

    public void setExistingQuantEvents(List<QuantEventInfo> existingQuantEvents)
    {
        this.existingQuantEvents = existingQuantEvents;
    }
}
