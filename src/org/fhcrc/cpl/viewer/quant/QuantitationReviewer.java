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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.SimpleXMLEventRewriter;
import org.fhcrc.cpl.toolbox.TempFileManager;
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
import javax.imageio.ImageIO;
import java.util.*;
import java.util.List;
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
    public JPanel navigationPanel;
    public JPanel curationPanel;
    public JPanel globalOpsPanel;
    



    DefaultTableModel propertiesTableModel;
    JTable propertiesTable;
    JScrollPane propertiesScrollPane;
//    JPanel propertiesPanel;

    JButton backButton;
    JButton forwardButton;
    protected ButtonGroup curationButtonGroup;
    ButtonModel unknownRadioButtonModel;
    ButtonModel goodRadioButtonModel;
    ButtonModel badRadioButtonModel;
    protected JButton saveChangesButton;
    protected JButton filterPepXMLButton;
    protected JButton helpButton;

    public JSplitPane splitPane;


    protected int leftPanelWidth = 250;
    protected int rightPanelWidth = 790;    

    protected int imagePanelWidth = 780;
    protected int fullWidth = 1000;
    protected int fullHeight = 1000;
    protected int propertiesWidth = leftPanelWidth - 20;

    protected int propertiesHeight = 250;
    protected int chartPaneHeight = 950;



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

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(5,5,5,5);
        gbc.weighty = 1;
        gbc.weightx = 1;


        leftPanel.setLayout(new GridBagLayout());
        leftPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
//        leftPanel.setMinimumSize(new Dimension(leftPanelWidth, 500));

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
//        propertiesScrollPane.setPreferredSize(new Dimension(leftPanelWidth, propertiesHeight));
        propertiesScrollPane.setMinimumSize(new Dimension(propertiesWidth, propertiesHeight));
//        leftPanel.add(propertiesScrollPane, gbc);

        //fields related to navigation
        navigationPanel = new JPanel();
        backButton = new JButton("<");
        backButton.setToolTipText("Previous Event");
        backButton.setMaximumSize(new Dimension(50, 30));
        forwardButton = new JButton(">");
        forwardButton.setToolTipText("Next Event");
        forwardButton.setMaximumSize(new Dimension(50, 30));

        helper.addListener(backButton, "buttonBack_actionPerformed");
        helper.addListener(forwardButton, "buttonForward_actionPerformed");
        gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        gbc.anchor = GridBagConstraints.WEST;
        navigationPanel.add(backButton, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        navigationPanel.add(forwardButton, gbc);
        gbc.fill = GridBagConstraints.BOTH;
        navigationPanel.setBorder(BorderFactory.createTitledBorder("Event"));
        gbc.anchor = GridBagConstraints.PAGE_START;
        

        curationPanel = new JPanel();
        curationPanel.setLayout(new GridBagLayout());
        curationPanel.setBorder(BorderFactory.createTitledBorder("Assessment"));

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

        gbc.anchor = GridBagConstraints.WEST;

        curationPanel.add(unknownRadioButton, gbc);
        curationPanel.add(badRadioButton, gbc);
        curationPanel.add(goodRadioButton, gbc);
        gbc.anchor = GridBagConstraints.PAGE_START;       



        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.CENTER;
        globalOpsPanel = new JPanel();
        globalOpsPanel.setBorder(BorderFactory.createTitledBorder("Actions"));          
        globalOpsPanel.setLayout(new GridBagLayout());
        saveChangesButton = new JButton("Save Changes");
        saveChangesButton.setToolTipText("Save your curation changes to the same file you loaded");

        helper.addListener(saveChangesButton, "buttonSaveChanges_actionPerformed");
        globalOpsPanel.add(saveChangesButton, gbc);
        filterPepXMLButton = new JButton("Filter PepXML...");
        filterPepXMLButton.setToolTipText("Filter curated 'Bad' events from a pepXML file, " +
                "saving the results to an output file");
        helper.addListener(filterPepXMLButton, "buttonFilterPepXML_actionPerformed");
        globalOpsPanel.add(filterPepXMLButton, gbc);

        helpButton = new JButton("Help");
        helpButton.setToolTipText("Get help with this interface");
        helper.addListener(helpButton, "buttonHelp_actionPerformed");
        globalOpsPanel.add(helpButton, gbc);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;



        leftPanel.addComponentListener(new LeftPanelResizeListener());
        gbc.weighty = 2;
        leftPanel.add(propertiesScrollPane, gbc);
        gbc.weighty = 1;
        gbc.anchor = GridBagConstraints.PAGE_END;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        
        leftPanel.add(curationPanel, gbc);
        //gbc.weighty = 2;
        leftPanel.add(navigationPanel, gbc);
        //gbc.weighty = 3;
        leftPanel.add(globalOpsPanel, gbc);
        //gbc.weighty = 1;
        gbc.anchor = GridBagConstraints.PAGE_START;




        //display charts
        multiChartDisplay = new TabbedMultiChartDisplayPanel();
        multiChartDisplay.setResizeDelayMS(0);




        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        rightPanel.addComponentListener(new RightPanelResizeListener());
        rightPanel.add(multiChartDisplay, gbc);
//        rightPanel.setBorder(BorderFactory.createLineBorder(Color.black));

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

    public void buttonHelp_actionPerformed(ActionEvent event)
    {
        String tempFileDummyString = "TEMP_FILE_CONTROLLER_QREVIEWER_HELP";
        File helpFile =
                TempFileManager.createTempFile(tempFileDummyString + ".html", tempFileDummyString);
        StringBuffer helpBuf = new StringBuffer();
        helpBuf.append("<H2>Qurate Help</H2>");
        helpBuf.append("<b>Qurate</b> allows you to view quantitative events in several different ways, " +
                "and to judge them as 'good' or 'bad' events.  'Bad' events can be filtered out of the " +
                "pepXML file that they came from.");
        helpBuf.append("<H3>Properties</H3>The Properties table shows you information about the " +
                "quantitative event. Multiple overlapping events (same peptide, same charge) are rolled " +
                "into a single event if they have similar ratios.  Any events other than the one described " +
                "in these properties will appear in 'OtherEventScans' and 'OtherEventMZs'.");
        helpBuf.append("<H3>Assessment</H3>Provide your assessment of the quantitative event here.  " +
                "Your assessments will not be saved unless you use the 'Save Changes' button.");
        helpBuf.append("<H3>Actions</H3>");
        String[] actionsHelpStrings = new String[]
                {
                    "<b>Save Changes:</b> Save your assessments back to the same file that you loaded.",
                    "<b>Filter PepXML...:</b> Filter all of your 'Bad' assessments out of the " +
                                "pepXML file that they came from.  You must specify both the source " +
                                "pepXML file and the output pepXML file (which may be the same)",
                };
        helpBuf.append(HtmlGenerator.createList(false, Arrays.asList(actionsHelpStrings)));
        helpBuf.append("<H3>Charts</H3>");
        helpBuf.append("The charts help you assess the quality of quantitative events.");

        String[] chartsHelpStrings = new String[]
                {
                    "<b>Spectrum:</b> A heatmap showing the intensity of the region of spectra around the" +
                            " quantitative event.  White means intensity 0, and higher intensities are " +
                            "darker blue.  Horizontal axis shows increasing scan number, Vertical axis " +
                            "shows increasing m/z.  Vertical red lines descending from the top indicate " +
                            "the scan range used by the quantitation software to determine the heavy " +
                            "intensity.  Red lines from the bottom indicate the light intensity range.  " +
                            "Red tick marks on the left and right indicate the light and heavy monoisotopic " +
                            "m/z values.  Red 'X' marks indicate identification events.",
                    "<b>Scans:</b> A separate chart for each scan, showing peak intensities.  All intensities " +
                            "are normalized to the highest-intensity value in any scan (which may or may not " +
                            "be included in the peptide quantitation).  Scan numbers are given at the left -- " +
                            "numbers in green were used in quantitation, numbers in red were not.",
                    "<b>3D:</b> A 3D representation of the same data in the Spectrum chart.  Perspective is that " +
                            "of someone standing at a scan higher than the last scan shown on the chart, near " +
                            "the m/z of the light monoisotopic peak.  Z axis is intensity.  'Veritical' red lines " +
                            "are the same as in the spectrum chart.  'Horizontal' red lines indicate the light " +
                            "and heavy monoisotopic peaks, just like the tick marks in the spectrum chart.  Blue " +
                            "'X' marks indicate the event used in this quantitative event.  Yellow 'X' marks indicate" +
                            "other identification events folded into this one representation."
                };
        helpBuf.append(HtmlGenerator.createList(false, Arrays.asList(chartsHelpStrings)));

        String helpDocString = HtmlGenerator.createHtmlDocument(helpBuf.toString(), "Qurate Help");
        PrintWriter helpPW = null;

        try
        {
            helpPW = new PrintWriter(helpFile);
            helpPW.println(helpDocString);
            helpPW.flush();
            HtmlViewerPanel.showFileInDialog(helpFile, "Qurate Help");
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("Error displaying help", e);
        }
        finally
        {
            if (helpPW != null)
                helpPW.close();
        }
        TempFileManager.deleteTempFiles(tempFileDummyString);
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
        navigationPanel.setBorder(BorderFactory.createTitledBorder(
                "Event " + (displayedEventIndex+1) + " / " + quantEvents.size()));

//        displayStatusLabel.setText();

        ButtonModel buttonModelToSelect = null;
//System.err.println("status " + quantEvents.get(displayedEventIndex).getCurationStatus());
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

//        PanelWithBlindImageChart spectrumChart =
//                new PanelWithBlindImageChart(quantEvent.getSpectrumFile(), "Spectrum");
//        PanelWithBlindImageChart scansChart =
//                new PanelWithBlindImageChart(quantEvent.getScansFile(), "Scans");

        List<PanelWithChart> multiChartPanels = multiChartDisplay.getChartPanels();

        if (multiChartPanels == null || multiChartPanels.isEmpty())
        {
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("Spectrum"));
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("Scans"));
            if (quantEvent.getFile3D() != null)
                multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("3D"));    
        }

        try
        {
            PanelWithBlindImageChart spectrumChart = (PanelWithBlindImageChart) multiChartPanels.get(0);
            PanelWithBlindImageChart scansChart = (PanelWithBlindImageChart) multiChartPanels.get(1);

            spectrumChart.setImage(ImageIO.read(quantEvent.getSpectrumFile()));
            scansChart.setImage(ImageIO.read(quantEvent.getScansFile()));

            if (quantEvent.getFile3D() != null)
            {
                PanelWithBlindImageChart chart3D = (PanelWithBlindImageChart) multiChartPanels.get(2);
                chart3D.setImage(ImageIO.read(quantEvent.getFile3D()));
            }
        }
        catch (IOException e)
        {
            ApplicationContext.errorMessage("Failure displaying charts",e);
        }
//        rightPanel.remove(multiChartDisplay);
//        multiChartDisplay = new TabbedMultiChartDisplayPanel();
//        multiChartDisplay.setResizeDelayMS(0);
//        GridBagConstraints lastGBC = new GridBagConstraints();
//        lastGBC.fill = GridBagConstraints.BOTH;
//        lastGBC.gridwidth = GridBagConstraints.REMAINDER;
//        lastGBC.insets = new Insets(0, 0, 0, 0);
//        lastGBC.anchor = GridBagConstraints.NORTH;
//        rightPanel.add(multiChartDisplay, lastGBC);

//        multiChartDisplay.addChartPanel(spectrumChart);
//        multiChartDisplay.addChartPanel(scansChart);
//
//        if (quantEvent.getFile3D() != null)
//        {
//            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart(quantEvent.getFile3D(), "3D"));
//        }

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
              multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()-5));
        }
        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }

    protected class LeftPanelResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
            propertiesScrollPane.setPreferredSize(new Dimension(leftPanel.getWidth()-15,
                    propertiesScrollPane.getHeight()));
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
