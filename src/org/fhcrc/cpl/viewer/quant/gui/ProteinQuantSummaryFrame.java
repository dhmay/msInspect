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
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.quant.QuantEventInfo;
import org.fhcrc.cpl.viewer.quant.QuantitationVisualizer;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.table.*;
import javax.imageio.ImageIO;
import java.io.File;
import java.io.StringWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.awt.*;
import java.awt.event.*;


/**
 * test
 */
public class ProteinQuantSummaryFrame extends JDialog
{
    protected static Logger _log = Logger.getLogger(ProteinQuantSummaryFrame.class);

    protected boolean done = false;

    protected String proteinName;

    protected float proteinRatio;
    protected List<QuantEventInfo> quantEvents;
    protected List<QuantEventInfo> selectedQuantEvents;
    protected File mzXmlDir;

    protected File outDir;
    protected File outFile;
    protected boolean appendOutput = true;

    protected List<Integer> shadedTableRows = new ArrayList<Integer>();

    protected int fullWidth = 900;
    protected int fullHeight = 600;

    protected int propertiesWidth = 400;
    protected int propertiesHeight = 500;

    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public JScrollPane eventsScrollPane;
    public JPanel eventsPanel;

    protected QuantEventInfo.QuantEventPropertiesTable eventPropertiesTable;
    protected JDialog eventPropertiesDialog;

    protected boolean shouldAddOverlappingEvents = true;

    JLabel proteinNameLabel = new JLabel("Protein: ");
    JLabel proteinRatioLabel  = new JLabel("Ratio: ");
    JButton buildChartsForSelectedButton = new JButton("Build Selected Charts");
    JButton showPropertiesButton = new JButton("Show Event Properties");


    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    //event properties
    protected QuantEventInfo.QuantEventsSummaryTable eventsTable;

    public ProteinQuantSummaryFrame()
    {
        initGUI();
    }

    public ProteinQuantSummaryFrame(File protXmlFile, File pepXmlFile, String proteinName,
                                    File outDir, File mzXmlDir, File outFile, boolean appendOutput)
            throws IllegalArgumentException
    {
        this();
        this.outDir = outDir;
        this.outFile = outFile;
        this.mzXmlDir = mzXmlDir;
        this.appendOutput = appendOutput;
        displayData(protXmlFile, pepXmlFile, proteinName);
    }

    protected void displayData(File protXmlFile, File pepXmlFile, String proteinName)
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
                for (Feature feature : featureSet.getFeatures())
                {
                    if (proteinPeptidesRatiosUsed.contains(MS2ExtraInfoDef.getFirstPeptide(feature)) &&
                            IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        QuantEventInfo quantEvent =
                                new QuantEventInfo(feature, MS2ExtraInfoDef.getFeatureSetBaseName(featureSet));
                        quantEvents.add(quantEvent);
                    }
                }
            }
            setMessage("Loaded all quantitation events.");
        }
        catch (Exception e)
        {
            errorMessage("Failed to load features from pepXML file",e);
            return;
        }
        Collections.sort(quantEvents,
                new QuantEventInfo.PeptideSequenceAscFractionAscChargeModificationsAscRatioAscComparator());
        displayEvents();
    }


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
            ((java.awt.Frame)getOwner()).setIconImage(ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif")));            
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

    public void buttonShowProperties_actionPerformed(ActionEvent event)
    {
        eventPropertiesDialog.setVisible(true);
    }

    public void buttonBuildCharts_actionPerformed(ActionEvent event)
    {
        selectedQuantEvents = eventsTable.getSelectedEvents();
        if (selectedQuantEvents.isEmpty())
        {
            infoMessage("No events selected");
            return;
        }

        if (shouldAddOverlappingEvents)
        {
            QuantitationVisualizer quantVisualizer = new QuantitationVisualizer();
            ApplicationContext.infoMessage("Finding all overlapping events, starting with " +
                    selectedQuantEvents.size());

            List<QuantEventInfo> allOverlappingEvents =
                    quantVisualizer.findNonOverlappingQuantEventsAllPeptides(quantEvents);
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

        QuantitationVisualizer quantVisualizer = new QuantitationVisualizer();
        quantVisualizer.setMzXmlDir(mzXmlDir);
        quantVisualizer.setOutDir(outDir);
        File outputFile = outFile;
        if (outputFile == null)
            outputFile = (new File(outDir, "quantitation_" + proteinName + ".tsv"));
        quantVisualizer.setOutTsvFile(outputFile);
        quantVisualizer.setOutHtmlFile(new File(outDir, "quantitation_" + proteinName + ".htnl"));
        quantVisualizer.setAppendTsvOutput(appendOutput);

        try
        {
            quantVisualizer.visualizeQuantEvents(selectedQuantEvents);
            infoMessage("Saved quantitation summary to " + quantVisualizer.getOutTsvFile().getAbsolutePath());            
        }
        catch (IOException e)
        {
            errorMessage("Error creating charts",e);
        }
        eventPropertiesDialog.dispose();
        done = true;
        setVisible(false);        
    }



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

        eventsTable.displayEvents(quantEvents);

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
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Display a dialog box with info message and stack trace
     * @param message
     * @param t
     */
    protected void errorMessage(String message, Throwable t)
    {
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
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
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
    }

    public boolean isDone()
    {
        return done;
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

//                ListSelectionModel lsm = (ListSelectionModel) e.getSource();
//                if (lsm.isSelectionEmpty()) {
//                    eventPropertiesTable.clearProperties();
//                } else
//                {
//                    // Find out which indexes are selected.
//                    int minIndex = lsm.getMinSelectionIndex();
//                    int maxIndex = lsm.getMaxSelectionIndex();
//                    if (minIndex == maxIndex)
//                        eventPropertiesTable.displayQuantEvent(quantEvents.asdfget(minIndex));
//                    else
//                        eventPropertiesTable.clearProperties();
//                }
            }
        }
    }
}
