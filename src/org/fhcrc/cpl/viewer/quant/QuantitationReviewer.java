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
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.SimpleXMLEventRewriter;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;

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
import java.awt.event.*;
import java.awt.*;
import java.io.*;

public class QuantitationReviewer extends JFrame
{
    List<QuantEventInfo> quantEvents;
    protected File quantFile;

    protected TabbedMultiChartDisplayPanel multiChartDisplay;

    protected int displayedEventIndex = 0;

    public JPanel contentPanel;
    public JPanel leftPanel;
    public JPanel rightPanel;
    public JPanel navigationPanel;
    public JPanel curationPanel;
    public JPanel theoreticalPeaksPanel;

    protected PanelWithPeakChart theoreticalPeaksChart;
    
    public Action helpAction = new HelpAction();
    public Action openFileAction;
    public Action saveAction = new SaveAction();
    public Action filterPepXMLAction;





    DefaultTableModel propertiesTableModel;
    JTable propertiesTable;
    JScrollPane propertiesScrollPane;
//    JPanel propertiesPanel;

    JButton backButton;
    JButton forwardButton;
    protected ButtonGroup quantCurationButtonGroup;
    ButtonModel unknownRadioButtonModel;
    ButtonModel goodRadioButtonModel;
    ButtonModel badRadioButtonModel;
    protected ButtonGroup idCurationButtonGroup;
    ButtonModel idUnknownRadioButtonModel;
    ButtonModel idGoodRadioButtonModel;
    ButtonModel idBadRadioButtonModel;
    protected JButton saveChangesButton;
    protected JButton filterPepXMLButton;
    protected JTextField commentTextField;

    public JSplitPane splitPane;

    public JPanel statusPanel;
    public JLabel messageLabel;



    protected int leftPanelWidth = 250;
    protected int rightPanelWidth = 790;    

    protected int imagePanelWidth = 780;
    protected int fullWidth = 1000;
    protected int fullHeight = 1000;
    protected int propertiesWidth = leftPanelWidth - 20;

    protected int propertiesHeight = 250;
    protected int chartPaneHeight = 950;
    protected int theoreticalPeaksPanelHeight = 150;



    protected static Logger _log = Logger.getLogger(QuantitationReviewer.class);

    public QuantitationReviewer()
    {
        initGUI();
    }

    public QuantitationReviewer(List<QuantEventInfo> quantEvents)
    {
        this();
        displayQuantEvents(quantEvents);
    }

    public QuantitationReviewer(File quantFile)
            throws IOException
    {
        this();
        displayQuantFile(quantFile);
    }

    public void displayQuantEvents(List<QuantEventInfo> quantEvents)
    {
        this.quantEvents = quantEvents;
        displayedEventIndex = 0;
        displayCurrentQuantEvent();
    }

    public void displayQuantFile(File quantFile)
            throws IOException
    {
        this.quantFile = quantFile;
        displayQuantEvents(QuantEventInfo.loadQuantEvents(quantFile));
        setMessage("Loaded quantitation events from file " + quantFile.getAbsolutePath());
    }

    protected void initGUI()
    {
        try
        {
            setIconImage(ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif")));
        }
        catch (Exception e)
        {
        }

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

        openFileAction = new OpenFileAction(this);
        filterPepXMLAction = new FilterPepXMLAction(this);


        try
        {
            JMenuBar jmenu = (JMenuBar)Localizer.getSwingEngine(this).render(
                            "org/fhcrc/cpl/viewer/quant/QuantitationReviewerMenu.xml");
            for (int i=0 ; i<jmenu.getMenuCount() ; i++)
                jmenu.getMenu(i).getPopupMenu().setLightWeightPopupEnabled(false);
            this.setJMenuBar(jmenu);
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_LOADING_MENUS"), x);
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

        JPanel quantCurationPanel = new JPanel();
        quantCurationPanel.setLayout(new GridBagLayout());
        quantCurationPanel.setBorder(BorderFactory.createTitledBorder("Qurate"));

        quantCurationButtonGroup = new ButtonGroup();
        JRadioButton unknownRadioButton = new JRadioButton("?");
        JRadioButton goodRadioButton = new JRadioButton("Good");
        JRadioButton badRadioButton = new JRadioButton("Bad");
        quantCurationButtonGroup.add(unknownRadioButton);
        quantCurationButtonGroup.add(goodRadioButton);
        quantCurationButtonGroup.add(badRadioButton);
        unknownRadioButtonModel = unknownRadioButton.getModel();
        goodRadioButtonModel = goodRadioButton.getModel();
        badRadioButtonModel = badRadioButton.getModel();

        helper.addListener(unknownRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(goodRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(badRadioButton, "buttonCuration_actionPerformed");

        gbc.anchor = GridBagConstraints.WEST;

        quantCurationPanel.add(unknownRadioButton, gbc);
        quantCurationPanel.add(badRadioButton, gbc);
        quantCurationPanel.add(goodRadioButton, gbc);
        gbc.anchor = GridBagConstraints.PAGE_START;       



        JPanel idCurationPanel = new JPanel();
        idCurationPanel.setLayout(new GridBagLayout());
        idCurationPanel.setBorder(BorderFactory.createTitledBorder("ID"));

        idCurationButtonGroup = new ButtonGroup();
        JRadioButton idUnknownRadioButton = new JRadioButton("?");
        JRadioButton idGoodRadioButton = new JRadioButton("Good");
        JRadioButton idBadRadioButton = new JRadioButton("Bad");
        idCurationButtonGroup.add(idUnknownRadioButton);
        idCurationButtonGroup.add(idGoodRadioButton);
        idCurationButtonGroup.add(idBadRadioButton);
        idUnknownRadioButtonModel = idUnknownRadioButton.getModel();
        idGoodRadioButtonModel = idGoodRadioButton.getModel();
        idBadRadioButtonModel = idBadRadioButton.getModel();

        helper.addListener(idUnknownRadioButton, "buttonIDCuration_actionPerformed");
        helper.addListener(idGoodRadioButton, "buttonIDCuration_actionPerformed");
        helper.addListener(idBadRadioButton, "buttonIDCuration_actionPerformed");

        gbc.anchor = GridBagConstraints.WEST;

        idCurationPanel.add(idUnknownRadioButton, gbc);
        idCurationPanel.add(idBadRadioButton, gbc);
        idCurationPanel.add(idGoodRadioButton, gbc);


        gbc.gridwidth = GridBagConstraints.RELATIVE;
        curationPanel.add(quantCurationPanel,gbc);

        gbc.gridwidth = GridBagConstraints.REMAINDER;
        curationPanel.add(idCurationPanel,gbc);

        commentTextField = new JTextField();
        commentTextField.addKeyListener(
                new KeyAdapter() {
                    public void keyReleased(KeyEvent e) {
                        if (quantEvents == null)
                            return;
                        QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);
                        //save the comment, being careful about tabs and new lines
                        quantEvent.setComment(commentTextField.getText().replace("\t"," ").replace("\n"," "));
                    }
                    public void keyTyped(KeyEvent e) {
                    }
                    public void keyPressed(KeyEvent e) {
                    }
                }
        );
        curationPanel.add(commentTextField,gbc);

        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.CENTER;
        theoreticalPeaksPanel = new JPanel();
        theoreticalPeaksPanel.setBorder(BorderFactory.createTitledBorder("Theoretical Peaks"));
        theoreticalPeaksPanel.setLayout(new GridBagLayout());
        theoreticalPeaksPanel.setMinimumSize(new Dimension(leftPanelWidth-10, theoreticalPeaksPanelHeight));
        theoreticalPeaksPanel.setMaximumSize(new Dimension(1200, theoreticalPeaksPanelHeight));

        showTheoreticalPeaks();


        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;



        leftPanel.addComponentListener(new LeftPanelResizeListener());
        gbc.weighty = 2;
        leftPanel.add(propertiesScrollPane, gbc);
        gbc.weighty = 1;
        gbc.anchor = GridBagConstraints.PAGE_END;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        
        leftPanel.add(curationPanel, gbc);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        leftPanel.add(theoreticalPeaksPanel, gbc);
        leftPanel.add(navigationPanel, gbc);

        gbc.fill = GridBagConstraints.BOTH;         

        gbc.weighty = 1;
        gbc.anchor = GridBagConstraints.PAGE_START;




        //display charts
        multiChartDisplay = new TabbedMultiChartDisplayPanel();
        multiChartDisplay.setResizeDelayMS(0);




        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        rightPanel.addComponentListener(new RightPanelResizeListener());
        rightPanel.add(multiChartDisplay, gbc);
//        rightPanel.setBorder(BorderFactory.createLineBorder(Color.black));

        messageLabel.setBackground(Color.WHITE);
        messageLabel.setFont(Font.decode("verdana plain 12"));
        messageLabel.setText(" ");            
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



    public void buttonCuration_actionPerformed(ActionEvent event)
    {
        QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

        ButtonModel selectedButtonModel = quantCurationButtonGroup.getSelection();
        if (selectedButtonModel == goodRadioButtonModel)
            quantEvent.setQuantCurationStatus(QuantEventInfo.CURATION_STATUS_GOOD);
        else if (selectedButtonModel == badRadioButtonModel)
            quantEvent.setQuantCurationStatus(QuantEventInfo.CURATION_STATUS_BAD);
        else
            quantEvent.setQuantCurationStatus(QuantEventInfo.CURATION_STATUS_UNKNOWN);
    }

    public void buttonIDCuration_actionPerformed(ActionEvent event)
    {
        QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

        ButtonModel selectedButtonModel = idCurationButtonGroup.getSelection();
        if (selectedButtonModel == idGoodRadioButtonModel)
            quantEvent.setIdCurationStatus(QuantEventInfo.CURATION_STATUS_GOOD);
        else if (selectedButtonModel == idBadRadioButtonModel)
        {
            //setting ID to bad is special -- also sets quantitation to bad
            quantEvent.setIdCurationStatus(QuantEventInfo.CURATION_STATUS_BAD);
            quantEvent.setQuantCurationStatus(QuantEventInfo.CURATION_STATUS_BAD);
            quantCurationButtonGroup.setSelected(badRadioButtonModel, true);
        }
        else
            quantEvent.setIdCurationStatus(QuantEventInfo.CURATION_STATUS_UNKNOWN);
    }

    protected void updateUIAfterChange()
    {
        QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

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
//System.err.println("status " + quantEvents.get(displayedEventIndex).getQuantCurationStatus());
        switch (quantEvent.getQuantCurationStatus())
        {
            case QuantEventInfo.CURATION_STATUS_UNKNOWN:
                buttonModelToSelect = unknownRadioButtonModel;
                break;
            case QuantEventInfo.CURATION_STATUS_GOOD:
                buttonModelToSelect = goodRadioButtonModel;
                break;
            case QuantEventInfo.CURATION_STATUS_BAD:
                buttonModelToSelect = badRadioButtonModel;
                break;
        }
        quantCurationButtonGroup.setSelected(buttonModelToSelect, true);

        switch (quantEvent.getIdCurationStatus())
        {
            case QuantEventInfo.CURATION_STATUS_UNKNOWN:
                buttonModelToSelect = idUnknownRadioButtonModel;
                break;
            case QuantEventInfo.CURATION_STATUS_GOOD:
                buttonModelToSelect = idGoodRadioButtonModel;
                break;
            case QuantEventInfo.CURATION_STATUS_BAD:
                buttonModelToSelect = idBadRadioButtonModel;
                break;
        }
        idCurationButtonGroup.setSelected(buttonModelToSelect, true);        

        multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()));
        multiChartDisplay.updateUI();

        commentTextField.setText(quantEvent.getComment() != null ? quantEvent.getComment() : "");

        showTheoreticalPeaks();
        
    }


    protected void showTheoreticalPeaks()
    {
        int chartWidth = Math.max(200, leftPanelWidth-30);
        int chartHeight = Math.max(theoreticalPeaksPanel.getHeight()-35, theoreticalPeaksPanelHeight-35);
        if (theoreticalPeaksChart != null)
        {
//            chartWidth = theoreticalPeaksChart.getWidth();
//            chartHeight = theoreticalPeaksChart.getHeight();

            theoreticalPeaksPanel.remove(0);
        }

        QuantEventInfo quantEvent = null;
        if (quantEvents != null)
        {
            quantEvent = quantEvents.get(displayedEventIndex);

            float lightNeutralMass = (quantEvent.getLightMz() - Spectrum.HYDROGEN_ION_MASS) *
                    quantEvent.getCharge();

            float[] lightTheoreticalPeaks = Spectrum.Poisson(lightNeutralMass);
            float[] lightPeakMzs = new float[6];
            for (int i=0; i<6; i++)
                lightPeakMzs[i] = quantEvent.getLightMz() + (Spectrum.HYDROGEN_ION_MASS * i / quantEvent.getCharge());
            theoreticalPeaksChart = new PanelWithPeakChart(lightPeakMzs, lightTheoreticalPeaks,
                    "Theoretical Peaks");

            float heavyNeutralMass = (quantEvent.getHeavyMz() - Spectrum.HYDROGEN_ION_MASS) *
                    quantEvent.getCharge();

            float[] heavyTheoreticalPeaks = Spectrum.Poisson(heavyNeutralMass);
            for (int i=0; i<heavyTheoreticalPeaks.length; i++)
                heavyTheoreticalPeaks[i] *= 1 / quantEvent.getRatio();

            float[] heavyPeakMzs = new float[6];
            for (int i=0; i<6; i++)
                heavyPeakMzs[i] = quantEvent.getHeavyMz() + (Spectrum.HYDROGEN_ION_MASS * i  / quantEvent.getCharge());
            for (int i=0; i<heavyPeakMzs.length; i++)
            {
                for (int j=0; j<lightPeakMzs.length; j++)
                    if (heavyPeakMzs[i] - lightPeakMzs[j] < 0.1)
                        heavyTheoreticalPeaks[i] += lightTheoreticalPeaks[j];
            }
            theoreticalPeaksChart.addData(heavyPeakMzs, heavyTheoreticalPeaks, "heavy");

            theoreticalPeaksChart.setPreferredSize(new Dimension(chartWidth, chartHeight));
            theoreticalPeaksChart.setSize(new Dimension(chartWidth, chartHeight));
            theoreticalPeaksChart.getChart().removeLegend();

            ((XYPlot)theoreticalPeaksChart.getPlot()).getDomainAxis().setVisible(false);
            ((XYPlot)theoreticalPeaksChart.getPlot()).getRangeAxis().setVisible(false);


            GridBagConstraints gbc = new GridBagConstraints();
            gbc.fill = GridBagConstraints.BOTH;
            gbc.anchor = GridBagConstraints.PAGE_START;
            gbc.gridwidth = GridBagConstraints.REMAINDER;
            theoreticalPeaksPanel.add(theoreticalPeaksChart, gbc);

            theoreticalPeaksPanel.setToolTipText("LightMass=" + lightNeutralMass + ", HeavyMass=" + heavyNeutralMass + 
                    ", Ratio=" + quantEvent.getRatio());
//            theoreticalPeaksChart.setToolTipText(theoreticalPeaksPanel.getToolTipText());
            theoreticalPeaksChart.updateUI();
        }
        else
        {
//            theoreticalPeaksChart = new PanelWithPeakChart();
        }


//        theoreticalPeaksChart.setPreferredSize(theoreticalPeaksChart.getSize());


    }


    protected void displayCurrentQuantEvent()
    {
        QuantEventInfo quantEvent = quantEvents.get(displayedEventIndex);

//        PanelWithBlindImageChart spectrumChart =
//                new PanelWithBlindImageChart(quantEvent.getSpectrumFile(), "Spectrum");
//        PanelWithBlindImageChart scansChart =
//                new PanelWithBlindImageChart(quantEvent.getScansFile(), "Scans");

        List<PanelWithChart> multiChartPanels = multiChartDisplay.getChartPanels();

        if (multiChartPanels == null || multiChartPanels.isEmpty())
        {
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("Intensity Sum"));
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("Spectrum"));
            multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("Scans"));
            if (quantEvent.getFile3D() != null)
                multiChartDisplay.addChartPanel(new PanelWithBlindImageChart("3D"));    
        }

        try
        {
            PanelWithBlindImageChart intensitySumChart = (PanelWithBlindImageChart) multiChartPanels.get(0);            
            PanelWithBlindImageChart spectrumChart = (PanelWithBlindImageChart) multiChartPanels.get(1);
            PanelWithBlindImageChart scansChart = (PanelWithBlindImageChart) multiChartPanels.get(2);


            spectrumChart.setImage(ImageIO.read(quantEvent.getSpectrumFile()));
            scansChart.setImage(ImageIO.read(quantEvent.getScansFile()));
            if (quantEvent.getIntensitySumFile() != null)
                intensitySumChart.setImage(ImageIO.read(quantEvent.getIntensitySumFile()));

            if (quantEvent.getFile3D() != null)
            {
                PanelWithBlindImageChart chart3D = (PanelWithBlindImageChart) multiChartPanels.get(3);
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
        for (String propName : QuantEventInfo.dataColumnNames)
        {
            if (propMap.containsKey(propName))
                addPropertyToModel(propName, propMap.get(propName));
        }

        String lightOrHeavyID = "Light";
        if (Math.abs(quantEvent.getHeavyMz() - quantEvent.getMz()) < 0.25)
            lightOrHeavyID = "Heavy";
        //this is to accommodate legacy files without MZ
        if (quantEvent.getMz() == 0)
            lightOrHeavyID = "Unknown";
        addPropertyToModel("ID Light/Heavy", lightOrHeavyID);                   


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

    public static void filterBadEventsFromFile(List<QuantEventInfo> quantEvents,
                                               File pepXmlFile, File outFile)
            throws IOException, XMLStreamException
    {
        Map<String, List<Integer>> fractionBadQuantScanListMap = new HashMap<String, List<Integer>>();
        Map<String, List<Integer>> fractionBadIDScanListMap = new HashMap<String, List<Integer>>();

        for (QuantEventInfo quantEvent : quantEvents)
        {
            if (quantEvent.getIdCurationStatus() == QuantEventInfo.CURATION_STATUS_BAD)
            {
                String fraction = quantEvent.getFraction();
                List<Integer> thisFractionList =
                        fractionBadIDScanListMap.get(fraction);
                if (thisFractionList == null)
                {
                    thisFractionList = new ArrayList<Integer>();
                    fractionBadIDScanListMap.put(fraction, thisFractionList);
                }
                thisFractionList.add(quantEvent.getScan());
                for (int scan : quantEvent.getOtherEventScans())
                    if (!thisFractionList.contains(scan))
                        thisFractionList.add(scan);
                ApplicationContext.infoMessage("Stripping ID for " + thisFractionList.size() +
                        " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
            //only if the ID was unknown or good do we check quant -- quant is automatically
            //filtered for bad IDs
            else if (quantEvent.getQuantCurationStatus() == QuantEventInfo.CURATION_STATUS_BAD)
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
                ApplicationContext.infoMessage("Stripping Quantitation for " + thisFractionList.size() +
                        " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
        }
        StripQuantPepXmlRewriter quantStripper = new StripQuantPepXmlRewriter(pepXmlFile, outFile,
                fractionBadIDScanListMap, fractionBadQuantScanListMap);
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
                        if (quantInfo.getQuantCurationStatus() == QuantitationVisualizer.QuantEventInfo.CURATION_STATUS_BAD &&
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


    /**
     *
     */
    static class StripQuantPepXmlRewriter extends SimpleXMLEventRewriter
    {
        Map<String, List<Integer>> fractionBadQuantScansMap;
        Map<String, List<Integer>> fractionBadIDScansMap;

        protected boolean insideSkippedQuantEvent = false;
        protected boolean insideSkippedSpectrumQuery = false;

        protected boolean insideScanWithSkippedEvent = false;

        List<Integer> currentFractionBadQuantScans;
        List<Integer> currentFractionBadIDScans;


        public StripQuantPepXmlRewriter(File inputFile, File outputFile,
                                        Map<String, List<Integer>> fractionBadIDScansMap,
                                        Map<String, List<Integer>> fractionBadQuantScansMap)
        {
            super(inputFile.getAbsolutePath(), outputFile.getAbsolutePath());
            this.fractionBadQuantScansMap = fractionBadQuantScansMap;
            this.fractionBadIDScansMap = fractionBadQuantScansMap;
        }

        public void add(XMLEvent event)
                throws XMLStreamException
        {
            if (!insideSkippedQuantEvent && !insideSkippedSpectrumQuery)
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
                currentFractionBadIDScans = fractionBadIDScansMap.get(baseName);
            }
            else if ("spectrum_query".equals(qname.getLocalPart()))
            {
                int scan = Integer.parseInt(event.getAttributeByName(new QName("start_scan")).getValue());
                if (currentFractionBadIDScans != null && currentFractionBadIDScans.contains(scan))
                {
                    insideSkippedSpectrumQuery = true;
                    _log.debug("Skipping ID for scan " + scan);
                }
                else if (currentFractionBadQuantScans != null && currentFractionBadQuantScans.contains(scan))
                {
                    insideScanWithSkippedEvent = true;
                    _log.debug("Skipping quantitation for scan " + scan);
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
            {
                insideScanWithSkippedEvent = false;
                insideSkippedSpectrumQuery = false;
            }
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

    public static class HelpAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
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
    }

    protected class OpenFileAction extends AbstractAction
    {
        protected Component parentComponent;

        public OpenFileAction(Component parentComponent)
        {
            this.parentComponent = parentComponent;
        }

        public void actionPerformed(ActionEvent event)
        {
            WorkbenchFileChooser wfc = new WorkbenchFileChooser();
            wfc.setDialogTitle("Open Quantitation Event File");
            int chooserStatus = wfc.showOpenDialog(parentComponent);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            quantFile = wfc.getSelectedFile();
            try
            {
                displayQuantFile(quantFile);
            }
            catch (IOException e)
            {
                ApplicationContext.errorMessage("Failed to open quantitation file " +
                        quantFile.getAbsolutePath(),e);
            }
        }
    }

    protected class SaveAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            try
            {
                QuantEventInfo.saveQuantEventsToTSV(quantEvents, quantFile, true, true);
                setMessage("Saved changes to file " + quantFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                errorMessage("ERROR: failed to save file " + quantFile.getAbsolutePath(), e);
            }
        }
    }

    protected class FilterPepXMLAction extends AbstractAction
    {
        protected Component parentComponent;

        public FilterPepXMLAction(Component parentComponent)
        {
            this.parentComponent = parentComponent;
        }

        public void actionPerformed(ActionEvent event)
        {
            WorkbenchFileChooser wfc = new WorkbenchFileChooser();
            wfc.setDialogTitle("Choose PepXML File to Filter");
            int chooserStatus = wfc.showOpenDialog(parentComponent);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            final File file = wfc.getSelectedFile();
            if (null == file)
                return;

            WorkbenchFileChooser wfc2 = new WorkbenchFileChooser();
            wfc2.setSelectedFile(file);
            wfc2.setDialogTitle("Choose Output File");

            chooserStatus = wfc2.showOpenDialog(parentComponent);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            final File outFile = wfc2.getSelectedFile();
            if (null == outFile)
                return;

            try
            {
                setMessage("Filtering bad events from file...");
                //do this in a separate thread so the UI doesn't freeze
                Thread filterThread = new Thread(new Runnable()
                {
                    public void run()
                    {
                        try
                        {
                            filterBadEventsFromFile(quantEvents, file, outFile);
                            setMessage("Filtered events saved to file " + outFile.getAbsolutePath());
                        }
                        catch (Exception e)
                        {
                            errorMessage("Failed to filter file " + file.getAbsolutePath(),e);
                        }
                    }
                });

                filterThread.start();
            }
            catch (Exception e)
            {
                errorMessage("Error filtering bad events from file " + file.getAbsolutePath() + " into file " +
                        outFile.getAbsolutePath(), e);
            }
        }
    }
}
