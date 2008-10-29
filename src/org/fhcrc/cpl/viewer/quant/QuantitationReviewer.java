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

import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.gui.chart.TabbedMultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.SimpleXMLEventRewriter;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.fhcrc.cpl.viewer.Localizer;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.xml.stream.events.XMLEvent;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.XMLStreamException;
import javax.xml.namespace.QName;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;
import java.awt.event.ComponentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ActionEvent;
import java.awt.*;
import java.io.*;

public class QuantitationReviewer extends JDialog
{
    List<QuantitationVisualizer.QuantEventInfo> quantEvents;
    protected File quantFile;

    protected TabbedMultiChartDisplayPanel multiChartDisplay;

    protected int displayedEventIndex = 0;

    public JPanel contentPanel;
    public JPanel leftPanel;
    public JPanel rightPanel;



    DefaultTableModel propertiesTableModel;
    JTable propertiesTable;
    JScrollPane propertiesScrollPane;

    JButton backButton;
    JButton forwardButton;
    protected JLabel displayStatusLabel;
    protected ButtonGroup curationButtonGroup;
    ButtonModel unknownRadioButtonModel;
    ButtonModel goodRadioButtonModel;
    ButtonModel badRadioButtonModel;
    protected JButton saveChangesButton;
    protected JButton filterPepXMLButton;
    public JSplitPane splitPane;



    protected int propertiesWidth = 180;
    protected int imagePanelWidth = 780;
    protected int fullWidth = 1000;
    protected int fullHeight = 1000;
    protected int propertiesHeight = 250;
    protected int chartPaneHeight = 950;

    protected int leftPanelWidth = 190;
    protected int rightPanelWidth = 790;

    protected static Logger _log = Logger.getLogger(QuantitationReviewer.class);

    public QuantitationReviewer()
    {
        initGUI();
    }

    public QuantitationReviewer(List<QuantitationVisualizer.QuantEventInfo> quantEvents)
    {
        this();

        this.quantEvents = quantEvents;
        displayedEventIndex = 0;
        displayCurrentQuantEvent();
    }

    public QuantitationReviewer(File quantFile)
            throws IOException
    {
        this(QuantitationVisualizer.QuantEventInfo.loadQuantEvents(quantFile));
        this.quantFile = quantFile;
    }

    protected void initGUI()
    {
        try
        {
            Localizer.renderSwixml("org/fhcrc/cpl/viewer/quant/QuantitationReviewer.xml",this);
            assert null != contentPanel;
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("error creating dialog", x);
            throw new RuntimeException(x);
        }

        setSize(fullWidth, fullHeight);

        setContentPane(contentPanel);

//        setLayout(new GridBagLayout());
//        setSize(fullWidth, fullHeight);
//        setMinimumSize(new Dimension(fullWidth, fullHeight));

        ListenerHelper helper = new ListenerHelper(this);

        GridBagConstraints secondToLastGBC = new GridBagConstraints();
        secondToLastGBC.fill = GridBagConstraints.BOTH;
        secondToLastGBC.gridwidth = GridBagConstraints.RELATIVE;
        secondToLastGBC.insets = new Insets(0, 0, 0, 0);
        secondToLastGBC.anchor = GridBagConstraints.PAGE_START;
        GridBagConstraints lastGBC = new GridBagConstraints();
        lastGBC.fill = GridBagConstraints.BOTH;
        lastGBC.gridwidth = GridBagConstraints.REMAINDER;
        lastGBC.insets = new Insets(0, 0, 0, 0);
        lastGBC.anchor = GridBagConstraints.PAGE_START;
        GridBagConstraints lastGBCNoFill = new GridBagConstraints();
        lastGBCNoFill.gridwidth = GridBagConstraints.REMAINDER;
        lastGBCNoFill.insets = new Insets(0, 0, 0, 0);
        lastGBCNoFill.anchor = GridBagConstraints.PAGE_START;




        leftPanel.setLayout(new GridBagLayout());
//        leftPanel.setPreferredSize(new Dimension(1000, fullHeight));
        leftPanel.setMinimumSize(new Dimension(leftPanelWidth+20, 500));


//        JPanel rightPanel = new JPanel();
//        rightPanel.setLayout(new GridBagLayout());
//        rightPanel.setSize(rightPanelWidth, fullHeight);

        //all cells uneditable
        propertiesTableModel = new DefaultTableModel(0, 2)
        {
            public boolean isCellEditable(int row, int column)
            {
                return false;
            }

            public Class getColumnClass(int columnIndex)
            {
                return String.class;
            }
        };
        propertiesTable = new JTable(propertiesTableModel)
        {
            //show tooltip with contents of cells
            public Component prepareRenderer(TableCellRenderer renderer,
                                             int rowIndex, int vColIndex)
            {
                Component c = super.prepareRenderer(renderer, rowIndex, vColIndex);
                if (c instanceof JComponent)
                {
                    JComponent jc = (JComponent)c;
                    jc.setToolTipText((String)getValueAt(rowIndex, vColIndex));
                }
                return c;
            }
        };
 
//        propertiesTable.setLayout(new GridBagLayout());
//        propertiesTable.setMinimumSize(new Dimension(propertiesWidth, propertiesHeight));
        propertiesTable.getColumnModel().getColumn(0).setHeaderValue(TextProvider.getText("PROPERTY_LOWERCASE"));
        propertiesTable.getColumnModel().getColumn(1).setHeaderValue(TextProvider.getText("VALUE_LOWERCASE"));
        propertiesScrollPane = new JScrollPane();
        propertiesScrollPane.setViewportView(propertiesTable);
        propertiesScrollPane.setPreferredSize(new Dimension(propertiesWidth, propertiesHeight));
        propertiesScrollPane.setMinimumSize(new Dimension(propertiesWidth, propertiesHeight));

        leftPanel.add(propertiesScrollPane, lastGBC);

        //fields related to curation and control
        JPanel etcPanel = new JPanel();
        etcPanel.setLayout(new GridBagLayout());
        backButton = new JButton("<");
        forwardButton = new JButton(">");
        helper.addListener(backButton, "buttonBack_actionPerformed");
        helper.addListener(forwardButton, "buttonForward_actionPerformed");
        etcPanel.add(backButton, secondToLastGBC);
        etcPanel.add(forwardButton, lastGBCNoFill);

        displayStatusLabel = new JLabel("Event");
        etcPanel.add(displayStatusLabel, lastGBC);

        JPanel curationPanel = new JPanel();
        curationPanel.setLayout(new GridBagLayout());
        JLabel quantAssessmentLabel = new JLabel("Assessment:");
        curationPanel.add(quantAssessmentLabel, lastGBC);
        curationButtonGroup = new ButtonGroup();
        JRadioButton unknownRadioButton = new JRadioButton("Unknown");
        JRadioButton goodRadioButton = new JRadioButton("Good");
        JRadioButton badRadioButton = new JRadioButton("Bad");
        curationButtonGroup.add(unknownRadioButton);
        curationButtonGroup.add(goodRadioButton);
        curationButtonGroup.add(badRadioButton);
        unknownRadioButtonModel = unknownRadioButton.getModel();
        goodRadioButtonModel = goodRadioButton.getModel();
        badRadioButtonModel = badRadioButton.getModel();

        helper.addListener(unknownRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(goodRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(badRadioButton, "buttonCuration_actionPerformed");

        curationPanel.add(unknownRadioButton, lastGBC);
        curationPanel.add(badRadioButton, lastGBC);
        curationPanel.add(goodRadioButton, lastGBC);

        saveChangesButton = new JButton("Save Changes");
        helper.addListener(saveChangesButton, "buttonSaveChanges_actionPerformed");
        curationPanel.add(saveChangesButton, lastGBC);

        filterPepXMLButton = new JButton("Filter PepXML...");
        helper.addListener(filterPepXMLButton, "buttonFilterPepXML_actionPerformed");
        curationPanel.add(filterPepXMLButton, lastGBC);

        etcPanel.add(curationPanel, lastGBC);

        leftPanel.add(etcPanel, lastGBC);
        leftPanel.addComponentListener(new LeftPanelResizeListener());


        //display charts
        multiChartDisplay = new TabbedMultiChartDisplayPanel();
        multiChartDisplay.setResizeDelayMS(0);


        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        rightPanel.addComponentListener(new RightPanelResizeListener());
        rightPanel.add(multiChartDisplay, lastGBC);
    }

    public void buttonBack_actionPerformed(ActionEvent event)
    {
        if (displayedEventIndex > 0)
        {
            displayedEventIndex--;
            displayCurrentQuantEvent();
        }
    }

    public void buttonForward_actionPerformed(ActionEvent event)
    {
        if (displayedEventIndex < quantEvents.size()-1)
        {
            displayedEventIndex++;
            displayCurrentQuantEvent();
        }
    }

    public void buttonSaveChanges_actionPerformed(ActionEvent event)
    {
        try
        {
            QuantitationVisualizer.QuantEventInfo.saveQuantEventsToTSV(quantEvents, quantFile, true, true);
            infoMessage("Saved changes to file " + quantFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            errorMessage("ERROR: failed to save file " + quantFile.getAbsolutePath(), e);
        }
    }

    public void buttonFilterPepXML_actionPerformed(ActionEvent event)
    {
        WorkbenchFileChooser wfc = new WorkbenchFileChooser();
        wfc.setDialogTitle("Choose PepXML File to Filter");
        int chooserStatus = wfc.showOpenDialog(this);
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return;
        File file = wfc.getSelectedFile();
        if (null == file)
            return;

        WorkbenchFileChooser wfc2 = new WorkbenchFileChooser();
        wfc2.setSelectedFile(file);
        wfc2.setDialogTitle("Choose Output File");

        chooserStatus = wfc2.showOpenDialog(this);
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return;
        File outFile = wfc2.getSelectedFile();
        if (null == outFile)
            return;

        try
        {
            filterBadEventsFromFile(quantEvents, file, outFile);
            infoMessage("Filtered bad quantitation events from file " + file.getAbsolutePath());
        }
        catch (Exception e)
        {
            errorMessage("Error filtering bad events from file " + file.getAbsolutePath() + " into file " +
                    outFile.getAbsolutePath(), e);
        }



    }

    public void buttonCuration_actionPerformed(ActionEvent event)
    {
        QuantitationVisualizer.QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

        ButtonModel selectedButtonModel = curationButtonGroup.getSelection();
        if (selectedButtonModel == goodRadioButtonModel)
            quantEvent.setCurationStatus(QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_GOOD);
        else if (selectedButtonModel == badRadioButtonModel)
            quantEvent.setCurationStatus(QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_BAD);
        else
            quantEvent.setCurationStatus(QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_UNKNOWN);
    }

    protected void updateUIAfterChange()
    {
        if (displayedEventIndex > 0)
            backButton.setEnabled(true);
        else
            backButton.setEnabled(false);
        if (displayedEventIndex < quantEvents.size()-1)
            forwardButton.setEnabled(true);
        else
            forwardButton.setEnabled(false);

        displayStatusLabel.setText("Event " + (displayedEventIndex+1) + " / " + quantEvents.size());

        ButtonModel buttonModelToSelect = null;
        switch (quantEvents.get(displayedEventIndex).getCurationStatus())
        {
            case QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_UNKNOWN:
                buttonModelToSelect = unknownRadioButtonModel;
                break;
            case QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_GOOD:
                buttonModelToSelect = goodRadioButtonModel;
                break;
            case QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_BAD:
                buttonModelToSelect = badRadioButtonModel;
                break;
        }
        curationButtonGroup.setSelected(buttonModelToSelect, true);
        multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()));
        multiChartDisplay.updateUI();
    }





    protected void displayCurrentQuantEvent()
    {
        QuantitationVisualizer.QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

        PanelWithBlindImageChart spectrumChart =
                new PanelWithBlindImageChart(quantEvent.getSpectrumFile(), "Spectrum");
        PanelWithBlindImageChart scansChart =
                new PanelWithBlindImageChart(quantEvent.getScansFile(), "Scans");

        rightPanel.remove(multiChartDisplay);
        multiChartDisplay = new TabbedMultiChartDisplayPanel();
        multiChartDisplay.setResizeDelayMS(0);
        GridBagConstraints lastGBC = new GridBagConstraints();
        lastGBC.fill = GridBagConstraints.BOTH;
        lastGBC.gridwidth = GridBagConstraints.REMAINDER;
        lastGBC.insets = new Insets(0, 0, 0, 0);
        lastGBC.anchor = GridBagConstraints.PAGE_START;        
        rightPanel.add(multiChartDisplay, lastGBC);

        multiChartDisplay.addChartPanel(spectrumChart);
        multiChartDisplay.addChartPanel(scansChart);

        if (quantEvent.getFile3D() != null)
        {
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart(quantEvent.getFile3D(), "3D"));
        }

        clearProperties();
        Map<String, String> propMap = quantEvent.getNameValueMapNoCharts();
        for (String propName : QuantitationVisualizer.QuantEventInfo.dataColumnNames)
        {
            if (propMap.containsKey(propName))
                addPropertyToModel(propName, propMap.get(propName));
        }
        updateUIAfterChange();
    }

    protected void clearProperties()
    {
        for (int i=propertiesTableModel.getRowCount()-1; i>=0; i--)
        {
            propertiesTableModel.setValueAt(null, i, 0);
            propertiesTableModel.setValueAt(null, i, 1);
        }
        propertiesTableModel.setRowCount(0);
    }

    protected void addPropertyToModel(Object propertyName, Object propertyValue)
    {
        int numRows = propertiesTableModel.getRowCount();
        propertiesTableModel.setRowCount(numRows + 1);
        propertiesTableModel.setValueAt(propertyName, numRows, 0);
        propertiesTableModel.setValueAt(propertyValue, numRows, 1);
    }

    public void setSize(int width, int height)
    {
        super.setSize(width, height);
    }

    public int getDisplayedEventIndex()
    {
        return displayedEventIndex;
    }

    public void setDisplayedEventIndex(int displayedEventIndex)
    {
        this.displayedEventIndex = displayedEventIndex;
    }

    protected class ResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
//System.err.println("***" + rightPanel.getWidth() + ", " +  rightPanel.getHeight());
//              multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()));
        }
        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }


    protected class RightPanelResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
              multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()));
        }
        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }

    protected class LeftPanelResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
            propertiesScrollPane.setPreferredSize(new Dimension(leftPanel.getWidth(), propertiesScrollPane.getHeight()));
            propertiesScrollPane.setSize(new Dimension(leftPanel.getWidth(), propertiesScrollPane.getHeight()));
//            propertiesScrollPane.updateUI();
        }
        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }


    protected void infoMessage(String message)
    {
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
    }

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

    public static void filterBadEventsFromFile(List<QuantitationVisualizer.QuantEventInfo> quantEvents,
                                               File pepXmlFile, File outFile)
            throws IOException, XMLStreamException
    {
        Map<String, List<Integer>> fractionBadQuantScanListMap = new HashMap<String, List<Integer>>();
        for (QuantitationVisualizer.QuantEventInfo quantEvent : quantEvents)
        {
            if (quantEvent.getCurationStatus() == QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_BAD)
            {
                String fraction = quantEvent.getFraction();
                List<Integer> thisFractionList =
                        fractionBadQuantScanListMap.get(fraction);
                if (thisFractionList == null)
                {
                    thisFractionList = new ArrayList<Integer>();
                    fractionBadQuantScanListMap.put(fraction, thisFractionList);
                }
                thisFractionList.add(quantEvent.getScan());
                for (int scan : quantEvent.getOtherEventScans())
                    if (!thisFractionList.contains(scan))
                        thisFractionList.add(scan);
                ApplicationContext.infoMessage("Stripping " + thisFractionList.size() + " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
        }
        StripQuantPepXmlRewriter quantStripper = new StripQuantPepXmlRewriter(pepXmlFile, outFile,
                fractionBadQuantScanListMap);
        quantStripper.rewrite();
        quantStripper.close();


        /*
        Map<String, Map<Integer, List<QuantitationVisualizer.QuantEventInfo>>> fractionScanQuantInfoListMap =
                new HashMap<String, Map<Integer, List<QuantitationVisualizer.QuantEventInfo>>>();

        for (QuantitationVisualizer.QuantEventInfo quantEvent : quantEvents)
        {
            String fraction = quantEvent.getFraction();
            Map<Integer, List<QuantitationVisualizer.QuantEventInfo>> thisFractionMap =
                    fractionScanQuantInfoListMap.get(fraction);
            if (thisFractionMap == null)
            {
                thisFractionMap = new HashMap<Integer, List<QuantitationVisualizer.QuantEventInfo>>();
                fractionScanQuantInfoListMap.put(fraction, thisFractionMap);
            }

            for (int scan : quantEvent.getAllScans())
            {
                List<QuantitationVisualizer.QuantEventInfo> thisScanList = thisFractionMap.get(scan);
                if (thisScanList == null)
                {
                    thisScanList = new ArrayList<QuantitationVisualizer.QuantEventInfo>();
                    thisFractionMap.put(scan, thisScanList);
                }
                thisScanList.add(quantEvent);
            }
        }

        PepXMLFeatureFileHandler.PepXMLFeatureSetIterator featureSetIterator =
                 new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
        List<File> tempFeatureFiles = new ArrayList<File>();
        int numSetsProcessed = 0;
        String tempFileDummyString = "DUMMY_TEMPFILE_STRING_QUANTREVIEWER";
        while (featureSetIterator.hasNext())
        {
            ApplicationContext.infoMessage("\tProcessing fraction " + (numSetsProcessed+1) + "...");

            FeatureSet featureSet = featureSetIterator.next();
            String fraction = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            Map<Integer, List<QuantitationVisualizer.QuantEventInfo>> thisFractionMap =
                    fractionScanQuantInfoListMap.get(fraction);
            if (thisFractionMap != null)
            {
                List<Feature> featuresToKeep = new ArrayList<Feature>();
                for (Feature feature : featureSet.getFeatures())
                {
                    int scan = feature.getScan();
                    List<QuantitationVisualizer.QuantEventInfo> thisScanList = thisFractionMap.get(scan);
                    if (thisScanList == null)
                        continue;
                    boolean keepThisFeature = true;
                    for (QuantitationVisualizer.QuantEventInfo quantInfo : thisScanList)
                    {
                        if (quantInfo.getCurationStatus() == QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_BAD &&
                                quantInfo.getPeptide().equals(MS2ExtraInfoDef.getFirstPeptide(feature)))
                        {
                            keepThisFeature = false;
                            _log.debug("Filtering feature with scan " + scan);
                        }

                    }
                    if (keepThisFeature)
                        featuresToKeep.add(feature);
                }
                featureSet.setFeatures(featuresToKeep.toArray(new Feature[featuresToKeep.size()]));
            }
            String baseName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
            if (baseName == null)
            {
                baseName = pepXmlFile.getName();
                if (numSetsProcessed > 0 || featureSetIterator.hasNext())
                    baseName = baseName + "_" + numSetsProcessed;
            }
            File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", tempFileDummyString);
            featureSet.savePepXml(thisFractionFile);
            tempFeatureFiles.add(thisFractionFile);
            numSetsProcessed++;
        }

        ApplicationContext.infoMessage("Saving output file " + outFile.getAbsolutePath() + "...");

        if (numSetsProcessed == 1)
        {
            FileReader in = new FileReader(tempFeatureFiles.get(0));
            FileWriter out = new FileWriter(outFile);
            int c;

            while ((c = in.read()) != -1)
                out.write(c);

            in.close();
            out.close();
        }
        else
        {
            ApplicationContext.infoMessage("\tCombining individual fraction files... " +
                    outFile.getAbsolutePath() + "...");
            new PepXMLFeatureFileHandler().combinePepXmlFiles(tempFeatureFiles, outFile);
        }
        ApplicationContext.infoMessage("Done.");
        TempFileManager.deleteTempFiles(tempFileDummyString);
        */
    }

       /**
     *
     */
    static class StripQuantPepXmlRewriter extends SimpleXMLEventRewriter
    {
        Map<String, List<Integer>> fractionBadQuantScansMap;
        protected boolean insideSkippedQuantEvent = false;
        protected boolean insideScanWithSkippedEvent = false;

        List<Integer> currentFractionBadQuantScans;

        public StripQuantPepXmlRewriter(File inputFile, File outputFile,
                                        Map<String, List<Integer>> fractionBadQuantScansMap)
        {
            super(inputFile.getAbsolutePath(), outputFile.getAbsolutePath());
            this.fractionBadQuantScansMap = fractionBadQuantScansMap;
        }

        public void add(XMLEvent event)
                throws XMLStreamException
        {
            if (!insideSkippedQuantEvent)
            {
                super.add(event);
            }
        }


        public void handleStartElement(StartElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();

            if ("msms_run_summary".equals(qname.getLocalPart()))
            {
                Attribute baseNameAttr = event.getAttributeByName(new QName("base_name"));
                String baseName = baseNameAttr.getValue();
                currentFractionBadQuantScans = fractionBadQuantScansMap.get(baseName);
            }
            else if ("spectrum_query".equals(qname.getLocalPart()))
            {
                int scan = Integer.parseInt(event.getAttributeByName(new QName("start_scan")).getValue());
                if (currentFractionBadQuantScans != null && currentFractionBadQuantScans.contains(scan))
                {
                    insideScanWithSkippedEvent = true;
                    _log.debug("Skipping scan " + scan);
                }
            }
            else if ("analysis_result".equals(qname.getLocalPart()))
            {
                if (insideScanWithSkippedEvent)
                {
                    String analysisType = event.getAttributeByName(new QName("analysis")).getValue();
                    if ("q3".equals(analysisType) || "xpress".equals(analysisType))
                        insideSkippedQuantEvent = true;
                }
            }

            add(event);
        }

        public void handleEndElement(EndElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();

            add(event);
            if (insideSkippedQuantEvent && "analysis_result".equals(qname.getLocalPart()))
                insideSkippedQuantEvent = false;
            else if ("spectrum_query".equals(qname.getLocalPart()))
                insideScanWithSkippedEvent = false;            
        }
    }

    public File getQuantFile()
    {
        return quantFile;
    }

    public void setQuantFile(File quantFile)
    {
        this.quantFile = quantFile;
    }
}
