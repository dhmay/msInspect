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

import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.gui.chart.TabbedMultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.toolbox.gui.widget.SplashFrame;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLEventRewriter;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.gui.ViewerInteractiveModuleFrame;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.ViewerUserManualGenerator;
import org.fhcrc.cpl.viewer.quant.commandline.PeptideQuantVisualizationCLM;
import org.fhcrc.cpl.viewer.quant.commandline.ProteinQuantChartsCLM;
import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.XYPlot;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
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
import java.net.URL;

/**
 * This is the GUI for Qurate.  It uses SwiXML for the menu and for the broad outlines, but most of it
 * is done right here.
 */
public class QuantitationReviewer extends JDialog
{
    //Quantitation events
    List<QuantEvent> quantEvents;
    //Loaded quantitation summary file
    protected File quantFile;

    //The multi-chart display panel
    protected TabbedMultiChartDisplayPanel multiChartDisplay;

    //This is the only way we keep track of the currently-displayed quantitation event.  An index into quantEvents
    protected int displayedEventIndex = 0;

    //everything
    public JPanel contentPanel;
    public JSplitPane splitPane;
    //left of splitpane
    public JPanel leftPanel;
    //right of splitpane
    public JPanel rightPanel;
    //for navigating between events
    public JPanel navigationPanel;
    JButton backButton;
    JButton forwardButton;
    JButton showEventSummaryButton;

    public ProteinQuantChartsCLM settingsCLM;

    protected SplashFrame splashFrame;
    protected final URL splashImageURL = QuantitationReviewer.class.getResource("qurate_splash.gif");

    //For assigning event statuses
    public JPanel curationPanel;
    protected ButtonGroup quantCurationButtonGroup;
    //need a reference to this in order to change title
    protected JRadioButton onePeakRatioRadioButton;
    ButtonModel unknownRadioButtonModel;
    ButtonModel goodRadioButtonModel;
    ButtonModel badRadioButtonModel;
    ButtonModel onePeakRadioButtonModel;

    protected ButtonGroup idCurationButtonGroup;
    ButtonModel idUnknownRadioButtonModel;
    ButtonModel idGoodRadioButtonModel;
    ButtonModel idBadRadioButtonModel;
    protected JButton saveChangesButton;
    protected JButton filterPepXMLButton;
    protected JTextField commentTextField;

    protected ProteinQuantSummaryFrame quantSummaryFrame;

    //theoretical peak distribution
    public JPanel theoreticalPeaksPanel;
    protected PanelWithPeakChart theoreticalPeaksChart;


    protected ProteinSummarySelectorFrame proteinSummarySelector;

    //menu actions
    public Action helpAction = new HelpAction();
    public Action exitAction = new ExitAction();    
    public Action openFileAction;
    public Action createChartsAction;
    public Action saveAction = new SaveAction(this);
    public Action filterPepXMLAction;
    public Action proteinSummaryAction;
    public Action aboutAction = new AbstractAction("About")
    {
        public void actionPerformed(ActionEvent e)
        {
            showSplashScreen();
        }
    };


    protected QuantEventsSummaryTable eventSummaryTable;
    protected JDialog eventSummaryDialog;


    //event properties
//    DefaultTableModel propertiesTableModel;
    protected QuantEvent.QuantEventPropertiesTable propertiesTable;
    protected JScrollPane propertiesScrollPane;

    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    //Sizes of things
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

    public QuantitationReviewer(boolean showSplash, boolean showFileOpen)
    {
        if (showSplash)
            showSplashScreen();

        initGUI();
        if (showFileOpen)
            openFileAction.actionPerformed(null);

        if (showSplash)
            splashFrame.dispose();
    }

    /**
     * No-arg constructor doesn't pop up a file chooser, but it does show splash screen
     */
    public QuantitationReviewer()
    {
        this(true, false);
    }

    public QuantitationReviewer(List<QuantEvent> quantEvents)
    {
        showSplashScreen();

        initGUI();
        displayQuantEvents(quantEvents);
        splashFrame.dispose();
    }

    public QuantitationReviewer(File quantFile)
            throws IOException
    {
        showSplashScreen();


        initGUI();
        displayQuantFile(quantFile);
        splashFrame.dispose();
    }

    public void displayQuantEvents(List<QuantEvent> quantEvents)
    {
//        this.quantEvents = new ArrayList<QuantEvent>(quantEvents);
        this.quantEvents = quantEvents;
//        Collections.sort(quantEvents, new QuantEvent.PeptideSequenceAscFractionAscChargeModificationsAscRatioAscComparator());
        displayedEventIndex = 0;
        eventSummaryTable.displayEvents(quantEvents);       
        displayCurrentQuantEvent();
    }

    public void displayQuantFile(File quantFile)
            throws IOException
    {
        this.quantFile = quantFile;
        quantEvents = QuantEvent.loadQuantEvents(quantFile);
        //handling for empty file
        if (quantEvents != null && !quantEvents.isEmpty())
        {
            displayQuantEvents(quantEvents);
            setMessage("Loaded quantitation events from file " + quantFile.getAbsolutePath());            
        }
        else
            setMessage("No events found in file " + quantFile.getAbsolutePath());
    }



    /**
     * Initialize all GUI components and display the UI
     */
    protected void initGUI()
    {

        settingsCLM = new ProteinQuantChartsCLM();

        setTitle("Qurate");
        try
        {
            setIconImage(ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif")));
        }
        catch (Exception e)
        {
        }

        try
        {
            Localizer.renderSwixml("org/fhcrc/cpl/viewer/quant/gui/QuantitationReviewer.xml",this);
            assert null != contentPanel;
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("error creating dialog", x);
            throw new RuntimeException(x);
        }

        //Menu
        openFileAction = new OpenFileAction(this);
        createChartsAction = new CreateChartsAction();
        filterPepXMLAction = new FilterPepXMLAction(this);
        proteinSummaryAction = new ProteinSummaryAction(this);

        try
        {
            JMenuBar jmenu = (JMenuBar)Localizer.getSwingEngine(this).render(
                    "org/fhcrc/cpl/viewer/quant/gui/QuantitationReviewerMenu.xml");
            for (int i=0 ; i<jmenu.getMenuCount() ; i++)
                jmenu.getMenu(i).getPopupMenu().setLightWeightPopupEnabled(false);
            this.setJMenuBar(jmenu);
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_LOADING_MENUS"), x);
            throw new RuntimeException(x);
        }

        //Global stuff
        setSize(fullWidth, fullHeight);
        setContentPane(contentPanel);
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

        //Properties panel stuff
        propertiesTable = new QuantEvent.QuantEventPropertiesTable();
        propertiesScrollPane = new JScrollPane();
        propertiesScrollPane.setViewportView(propertiesTable);
        propertiesScrollPane.setMinimumSize(new Dimension(propertiesWidth, propertiesHeight));

        //event summary table; disembodied
        eventSummaryTable = new QuantEventsSummaryTable();
        eventSummaryTable.setVisible(true);
        eventSummaryTable.hideSelectionColumn();
        ListSelectionModel tableSelectionModel = eventSummaryTable.getSelectionModel();
        tableSelectionModel.addListSelectionListener(new EventSummaryTableListSelectionHandler());        
        JScrollPane eventSummaryScrollPane = new JScrollPane();
        eventSummaryScrollPane.setViewportView(eventSummaryTable);
        eventSummaryScrollPane.setSize(propertiesWidth, propertiesHeight);
        eventSummaryDialog = new JDialog(this, "All Events");
        eventSummaryDialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        eventSummaryDialog.setSize(950, 450);
        eventSummaryDialog.setContentPane(eventSummaryScrollPane);

        //fields related to navigation
        navigationPanel = new JPanel();
        backButton = new JButton("<");
        backButton.setToolTipText("Previous Event");
        backButton.setMaximumSize(new Dimension(50, 30));
        backButton.setEnabled(false);
        forwardButton = new JButton(">");
        forwardButton.setToolTipText("Next Event");
        forwardButton.setMaximumSize(new Dimension(50, 30));
        forwardButton.setEnabled(false);
        showEventSummaryButton = new JButton("Show All");
        showEventSummaryButton.setToolTipText("Show all events in a table");
        showEventSummaryButton.setEnabled(false);

        helper.addListener(backButton, "buttonBack_actionPerformed");
        helper.addListener(forwardButton, "buttonForward_actionPerformed");
        helper.addListener(showEventSummaryButton, "buttonShowEventSummary_actionPerformed");

        gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        gbc.anchor = GridBagConstraints.WEST;
        navigationPanel.add(backButton, gbc);
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        navigationPanel.add(forwardButton, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        navigationPanel.add(showEventSummaryButton, gbc);        
        gbc.fill = GridBagConstraints.BOTH;
        navigationPanel.setBorder(BorderFactory.createTitledBorder("Event"));
        gbc.anchor = GridBagConstraints.PAGE_START;
        
        //Fields related to curation of events
        curationPanel = new JPanel();
        curationPanel.setLayout(new GridBagLayout());
        curationPanel.setBorder(BorderFactory.createTitledBorder("Assessment"));
        //Quantitation curation
        JPanel quantCurationPanel = new JPanel();
        quantCurationPanel.setLayout(new GridBagLayout());
        quantCurationPanel.setBorder(BorderFactory.createTitledBorder("Quantitation"));
        quantCurationButtonGroup = new ButtonGroup();
        JRadioButton unknownRadioButton = new JRadioButton("?");
        JRadioButton goodRadioButton = new JRadioButton("Good");
        JRadioButton badRadioButton = new JRadioButton("Bad");
        onePeakRatioRadioButton = new JRadioButton("1-Peak");

        unknownRadioButton.setEnabled(false);
        goodRadioButton.setEnabled(false);
        badRadioButton.setEnabled(false);
        onePeakRatioRadioButton.setEnabled(false);

        quantCurationButtonGroup.add(unknownRadioButton);
        quantCurationButtonGroup.add(goodRadioButton);
        quantCurationButtonGroup.add(badRadioButton);
        quantCurationButtonGroup.add(onePeakRatioRadioButton);

        unknownRadioButtonModel = unknownRadioButton.getModel();
        goodRadioButtonModel = goodRadioButton.getModel();
        badRadioButtonModel = badRadioButton.getModel();
        onePeakRadioButtonModel = onePeakRatioRadioButton.getModel();

        helper.addListener(unknownRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(goodRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(badRadioButton, "buttonCuration_actionPerformed");
        helper.addListener(onePeakRadioButtonModel, "buttonCuration_actionPerformed");

        gbc.anchor = GridBagConstraints.WEST;
        quantCurationPanel.add(unknownRadioButton, gbc);
        quantCurationPanel.add(badRadioButton, gbc);
        quantCurationPanel.add(goodRadioButton, gbc);
        quantCurationPanel.add(onePeakRatioRadioButton, gbc);

        gbc.anchor = GridBagConstraints.PAGE_START;
        //ID curation
        JPanel idCurationPanel = new JPanel();
        idCurationPanel.setLayout(new GridBagLayout());
        idCurationPanel.setBorder(BorderFactory.createTitledBorder("ID"));
        idCurationButtonGroup = new ButtonGroup();
        JRadioButton idUnknownRadioButton = new JRadioButton("?");
        JRadioButton idGoodRadioButton = new JRadioButton("Good");
        JRadioButton idBadRadioButton = new JRadioButton("Bad");
        idUnknownRadioButton.setEnabled(false);
        idGoodRadioButton.setEnabled(false);
        idBadRadioButton.setEnabled(false);        

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

        //curation comment
        commentTextField = new JTextField();
        commentTextField.setToolTipText("Comment on this event");
        //saves after every keypress.  Would be more efficient to save when navigating away or saving to file
        commentTextField.addKeyListener(
                new KeyAdapter() {
                    public void keyReleased(KeyEvent e) {
                        if (quantEvents == null)
                            return;
                        QuantEvent quantEvent = quantEvents.get(displayedEventIndex);
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

        //Theoretical peak distribution
        gbc.fill = GridBagConstraints.NONE;
        gbc.anchor = GridBagConstraints.CENTER;
        theoreticalPeaksPanel = new JPanel();
        theoreticalPeaksPanel.setBorder(BorderFactory.createTitledBorder("Theoretical Peaks"));
        theoreticalPeaksPanel.setLayout(new GridBagLayout());
        theoreticalPeaksPanel.setMinimumSize(new Dimension(leftPanelWidth-10, theoreticalPeaksPanelHeight));
        theoreticalPeaksPanel.setMaximumSize(new Dimension(1200, theoreticalPeaksPanelHeight));
        showTheoreticalPeaks();


        //Add everything to the left panel
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

        //Chart display
        multiChartDisplay = new TabbedMultiChartDisplayPanel();
        multiChartDisplay.setResizeDelayMS(0);

        rightPanel.addComponentListener(new RightPanelResizeListener());
        rightPanel.add(multiChartDisplay, gbc);

        //status message
        messageLabel.setBackground(Color.WHITE);
        messageLabel.setFont(Font.decode("verdana plain 12"));
        messageLabel.setText(" ");

        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        //paranoia.  Sometimes it seems Qurate doesn't exit when you close every window
        addWindowStateListener(new WindowStateListener()
        {
            public void windowStateChanged(WindowEvent e)
            {
                if (e.getNewState() == WindowEvent.WINDOW_CLOSED)
                {
                    dispose();
                    System.exit(0);
                }
            }
        });

    }

    //button actions

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

    public void buttonShowEventSummary_actionPerformed(ActionEvent event)
    {
        eventSummaryDialog.setVisible(true);
    }

    public void buttonCuration_actionPerformed(ActionEvent event)
    {
        QuantEvent quantEvent = quantEvents.get(displayedEventIndex);

        ButtonModel selectedButtonModel = quantCurationButtonGroup.getSelection();
        if (selectedButtonModel == goodRadioButtonModel)
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_GOOD);
        else if (selectedButtonModel == badRadioButtonModel)
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_BAD);
        else if (selectedButtonModel == onePeakRadioButtonModel)
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_RATIO_ONEPEAK);        
        else
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_UNKNOWN);
    }

    public void buttonIDCuration_actionPerformed(ActionEvent event)
    {
        QuantEvent quantEvent = quantEvents.get(displayedEventIndex);

        ButtonModel selectedButtonModel = idCurationButtonGroup.getSelection();
        if (selectedButtonModel == idGoodRadioButtonModel)
            quantEvent.setIdCurationStatus(QuantEvent.CURATION_STATUS_GOOD);
        else if (selectedButtonModel == idBadRadioButtonModel)
        {
            //setting ID to bad is special -- also sets quantitation to bad
            quantEvent.setIdCurationStatus(QuantEvent.CURATION_STATUS_BAD);
            quantEvent.setQuantCurationStatus(QuantEvent.CURATION_STATUS_BAD);
            quantCurationButtonGroup.setSelected(badRadioButtonModel, true);
        }
        else
            quantEvent.setIdCurationStatus(QuantEvent.CURATION_STATUS_UNKNOWN);
    }

    /**
     * Update lots of UI components after a change of quantitation event
     */
    protected void updateUIAfterChange()
    {
        QuantEvent quantEvent = quantEvents.get(displayedEventIndex);

        if (displayedEventIndex > 0)
            backButton.setEnabled(true);
        else
            backButton.setEnabled(false);
        if (displayedEventIndex < quantEvents.size()-1)
            forwardButton.setEnabled(true);
        else
            forwardButton.setEnabled(false);
        showEventSummaryButton.setEnabled(true);

        onePeakRatioRadioButton.setText("" + Rounder.round(quantEvent.getRatioOnePeak(), 2));
        Enumeration<AbstractButton> quantButtons = quantCurationButtonGroup.getElements();
        while (quantButtons.hasMoreElements())
            quantButtons.nextElement().setEnabled(true);
        Enumeration<AbstractButton> idButtons = idCurationButtonGroup.getElements();
        while (idButtons.hasMoreElements())
            idButtons.nextElement().setEnabled(true);        

        navigationPanel.setBorder(BorderFactory.createTitledBorder(
                "Event " + (displayedEventIndex+1) + " / " + quantEvents.size()));

        ButtonModel buttonModelToSelect = null;
        switch (quantEvent.getQuantCurationStatus())
        {
            case QuantEvent.CURATION_STATUS_UNKNOWN:
                buttonModelToSelect = unknownRadioButtonModel;
                break;
            case QuantEvent.CURATION_STATUS_GOOD:
                buttonModelToSelect = goodRadioButtonModel;
                break;
            case QuantEvent.CURATION_STATUS_BAD:
                buttonModelToSelect = badRadioButtonModel;
                break;
            case QuantEvent.CURATION_STATUS_RATIO_ONEPEAK:
                buttonModelToSelect = onePeakRadioButtonModel;
                break;
        }
        quantCurationButtonGroup.setSelected(buttonModelToSelect, true);

        switch (quantEvent.getIdCurationStatus())
        {
            case QuantEvent.CURATION_STATUS_UNKNOWN:
                buttonModelToSelect = idUnknownRadioButtonModel;
                break;
            case QuantEvent.CURATION_STATUS_GOOD:
                buttonModelToSelect = idGoodRadioButtonModel;
                break;
            case QuantEvent.CURATION_STATUS_BAD:
                buttonModelToSelect = idBadRadioButtonModel;
                break;
        }
        idCurationButtonGroup.setSelected(buttonModelToSelect, true);        

        multiChartDisplay.setPreferredSize(new Dimension(rightPanel.getWidth(), rightPanel.getHeight()));
        multiChartDisplay.updateUI();

        commentTextField.setText(quantEvent.getComment() != null ? quantEvent.getComment() : "");
        eventSummaryTable.getSelectionModel().setSelectionInterval(displayedEventIndex, displayedEventIndex);
        
        showTheoreticalPeaks();
    }


    /**
     * Calculate and show theoretical isotopic distribution peaks, with light encroaching on heavy
     * if necessary
     */
    protected void showTheoreticalPeaks()
    {
        int chartWidth = Math.max(200, leftPanelWidth-30);
        int chartHeight = Math.max(theoreticalPeaksPanel.getHeight()-35, theoreticalPeaksPanelHeight-35);
        if (theoreticalPeaksChart != null && theoreticalPeaksChart.getComponentCount() > 0)
            theoreticalPeaksPanel.remove(0);

        QuantEvent quantEvent = null;
        if (quantEvents != null)
        {
            quantEvent = quantEvents.get(displayedEventIndex);

            float lightNeutralMass = (quantEvent.getLightMz() - Spectrum.HYDROGEN_ION_MASS) *
                    quantEvent.getCharge();
            //Hardcoded 6 is number of peakd returned by Poisson()
            //Can't just use the result of Poisson(), as that's a static array and we're gonna mess with it
            float[] lightTheoreticalPeaks = new float[6];
            System.arraycopy(Spectrum.Poisson(lightNeutralMass), 0, lightTheoreticalPeaks, 0,
                    lightTheoreticalPeaks.length);

            float[] lightPeakMzs = new float[6];
            for (int i=0; i<6; i++)
                lightPeakMzs[i] = quantEvent.getLightMz() + (Spectrum.HYDROGEN_ION_MASS * i / quantEvent.getCharge());
            theoreticalPeaksChart = new PanelWithPeakChart(lightPeakMzs, lightTheoreticalPeaks,
                    "Theoretical Peaks");

            float heavyNeutralMass = (quantEvent.getHeavyMz() - Spectrum.HYDROGEN_ION_MASS) *
                    quantEvent.getCharge();

            float[] heavyTheoreticalPeaks = new float[6];
            System.arraycopy(Spectrum.Poisson(heavyNeutralMass), 0, heavyTheoreticalPeaks, 0, heavyTheoreticalPeaks.length);
            //fix 0 ratios at 0.001 to avoid divide by 0
            for (int i=0; i<heavyTheoreticalPeaks.length; i++)
                heavyTheoreticalPeaks[i] *= 1 / Math.max(quantEvent.getRatio(), 0.001f);

            float[] heavyPeakMzs = new float[6];
            for (int i=0; i<6; i++)
                heavyPeakMzs[i] = quantEvent.getHeavyMz() +
                        (Spectrum.HYDROGEN_ION_MASS * i  / quantEvent.getCharge());
            //Adjust heavy peaks if light peaks intrude.  Light appears in front of heavy
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

            //remove axes from chart
            ((XYPlot)theoreticalPeaksChart.getPlot()).getDomainAxis().setVisible(false);
            ((XYPlot)theoreticalPeaksChart.getPlot()).getRangeAxis().setVisible(false);

            GridBagConstraints gbc = new GridBagConstraints();
            gbc.fill = GridBagConstraints.BOTH;
            gbc.anchor = GridBagConstraints.PAGE_START;
            gbc.gridwidth = GridBagConstraints.REMAINDER;
            theoreticalPeaksPanel.add(theoreticalPeaksChart, gbc);

            theoreticalPeaksPanel.setToolTipText("LightMass=" + lightNeutralMass + ", HeavyMass=" + heavyNeutralMass + 
                    ", Ratio=" + quantEvent.getRatio());
            theoreticalPeaksChart.updateUI();
        }
    }

    /**
     * Take care of the charts and the properties panel
     */
    protected void displayCurrentQuantEvent()
    {
        QuantEvent quantEvent = quantEvents.get(displayedEventIndex);

        List<PanelWithChart> multiChartPanels = multiChartDisplay.getChartPanels();

        //first-time initialization
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

        propertiesTable.displayQuantEvent(quantEvent);

        updateUIAfterChange();
    }

    public int getDisplayedEventIndex()
    {
        return displayedEventIndex;
    }

    public void setDisplayedEventIndex(int displayedEventIndex)
    {
        this.displayedEventIndex = displayedEventIndex;
    }

    /**
     * Manually manage the size of the multi-chart panel
     */
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

    /**
     * Manually manage the size of the properties table
     */
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
        ApplicationContext.errorMessage(message, t);
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information",
                                      JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Remove all the events the user has designated as 'bad' from the pepXML file they came from
     * TODO: report how many events weren't found
     * @param quantEvents
     * @param pepXmlFile
     * @param outFile
     * @throws IOException
     * @throws XMLStreamException
     */
    public static void filterBadEventsFromFile(List<QuantEvent> quantEvents,
                                               File pepXmlFile, File outFile)
            throws IOException, XMLStreamException
    {
        Map<String, List<Integer>> fractionBadQuantScanListMap = new HashMap<String, List<Integer>>();
        Map<String, List<Integer>> fractionBadIDScanListMap = new HashMap<String, List<Integer>>();
        Map<String, Map<Integer, Float>> fractionScanNewRatioMap = new HashMap<String, Map<Integer, Float>>();


        for (QuantEvent quantEvent : quantEvents)
        {
            String fraction = quantEvent.getFraction();

            if (quantEvent.getIdCurationStatus() == QuantEvent.CURATION_STATUS_BAD)
            {
                List<Integer> thisFractionList =
                        fractionBadIDScanListMap.get(fraction);
                if (thisFractionList == null)
                {
                    thisFractionList = new ArrayList<Integer>();
                    fractionBadIDScanListMap.put(fraction, thisFractionList);
                }
                thisFractionList.add(quantEvent.getScan());
                for (QuantEvent otherEvent : quantEvent.getOtherEvents())
                    if (!thisFractionList.contains(otherEvent.getScan()))
                        thisFractionList.add(otherEvent.getScan());
                ApplicationContext.infoMessage("Stripping ID for " + (quantEvent.getOtherEvents().size() + 1) +
                        " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
            //only if the ID was unknown or good do we check quant -- quant is automatically
            //filtered for bad IDs
            else if (quantEvent.getQuantCurationStatus() == QuantEvent.CURATION_STATUS_BAD)
            {
                List<Integer> thisFractionList =
                        fractionBadQuantScanListMap.get(fraction);
                if (thisFractionList == null)
                {
                    thisFractionList = new ArrayList<Integer>();
                    fractionBadQuantScanListMap.put(fraction, thisFractionList);
                }
                thisFractionList.add(quantEvent.getScan());
                for (QuantEvent otherEvent : quantEvent.getOtherEvents())
                    if (!thisFractionList.contains(otherEvent.getScan()))
                        thisFractionList.add(otherEvent.getScan());
                ApplicationContext.infoMessage("Stripping Quantitation for " + (quantEvent.getOtherEvents().size() + 1) +
                        " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
            //dhmay adding 20090723
            else if (quantEvent.getQuantCurationStatus() == QuantEvent.CURATION_STATUS_RATIO_ONEPEAK)
            {
                Map<Integer, Float> thisFractionMap = fractionScanNewRatioMap.get(fraction);
                if (thisFractionMap == null)
                {
                    thisFractionMap = new HashMap<Integer, Float>();
                    fractionScanNewRatioMap.put(fraction, thisFractionMap);
                }
                float newRatio = quantEvent.getRatioOnePeak();
                thisFractionMap.put(quantEvent.getScan(), newRatio);
                for (QuantEvent otherEvent : quantEvent.getOtherEvents())
                    if (!thisFractionMap.containsKey(otherEvent.getScan()))
                        thisFractionMap.put(otherEvent.getScan(), newRatio);
                ApplicationContext.infoMessage("Changing ratios to " + newRatio + " for " +
                        (quantEvent.getOtherEvents().size() + 1) +
                        " events for peptide " +
                        quantEvent.getPeptide() + " from fraction " + fraction);
            }
        }
        StripQuantOrIDPepXmlRewriter quantStripper = new StripQuantOrIDPepXmlRewriter(pepXmlFile, outFile,
                fractionBadIDScanListMap, fractionBadQuantScanListMap, fractionScanNewRatioMap);
        quantStripper.rewrite();
        quantStripper.close();
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

    /**
     *  pepXML rewriter that can strip out quantitation events or entire spectrum_query tags.
     * dhmay 20090723: adding ability to alter ratios of events
     */
    static class StripQuantOrIDPepXmlRewriter extends SimpleXMLEventRewriter
    {
        //Track the bad stuff by fraction name and scan number
        Map<String, List<Integer>> fractionBadQuantScansMap;
        Map<String, List<Integer>> fractionBadIDScansMap;
        Map<String, Map<Integer, Float>> fractionScanNewRatioMap;

        protected boolean insideSkippedQuantEvent = false;
        protected boolean insideSkippedSpectrumQuery = false;

        protected boolean insideScanWithSkippedEvent = false;
        protected boolean insideScanWithRatioChange = false;

        List<Integer> currentFractionBadQuantScans;
        List<Integer> currentFractionBadIDScans;
        Map<Integer, Float> currentFractionScanNewRatioMap;

        float currentScanNewRatio;


        public StripQuantOrIDPepXmlRewriter(File inputFile, File outputFile,
                                        Map<String, List<Integer>> fractionBadIDScansMap,
                                        Map<String, List<Integer>> fractionBadQuantScansMap,
                                        Map<String, Map<Integer, Float>> fractionScanNewRatioMap)
        {
            super(inputFile.getAbsolutePath(), outputFile.getAbsolutePath());
            this.fractionBadQuantScansMap = fractionBadQuantScansMap;
            this.fractionBadIDScansMap = fractionBadIDScansMap;
            this.fractionScanNewRatioMap = fractionScanNewRatioMap;

        }

        public void add(XMLEvent event)
                throws XMLStreamException
        {
            if (!insideSkippedQuantEvent && !insideSkippedSpectrumQuery)
            {
                super.add(event);
            }
        }

        /**
         * special handling for keeping track of fraction and dealing with the skipped stuff
         * @param event
         * @throws XMLStreamException
         */
        public void handleStartElement(StartElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();
            String elementName = qname.getLocalPart();
            XMLEvent eventToAdd = event;

            if ("msms_run_summary".equals(elementName))
            {
                Attribute baseNameAttr = event.getAttributeByName(new QName("base_name"));
                String baseName = baseNameAttr.getValue();
                currentFractionBadQuantScans = fractionBadQuantScansMap.get(baseName);
                currentFractionBadIDScans = fractionBadIDScansMap.get(baseName);
                currentFractionScanNewRatioMap = fractionScanNewRatioMap.get(baseName);
            }
            else if ("spectrum_query".equals(elementName))
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
                else if (currentFractionScanNewRatioMap != null && currentFractionScanNewRatioMap.containsKey(scan))
                {
                    insideScanWithRatioChange = true;
                    currentScanNewRatio = currentFractionScanNewRatioMap.get(scan);
                }
            }
            else if ("analysis_result".equals(elementName))
            {
                if (insideScanWithSkippedEvent)
                {
                    String analysisType = event.getAttributeByName(new QName("analysis")).getValue();
                    if ("q3".equals(analysisType) || "xpress".equals(analysisType))
                        insideSkippedQuantEvent = true;
                }
            }
            //Adjust Q3 or XPress ratios.  Note: doesn't affect Q3's Q2 ratio or KL scores
            else if ("xpressratio_result".equals(elementName) || "q3ratio_result".equals(elementName))
            {
                if (insideScanWithRatioChange)
                {
                    float lightArea = (float) Rounder.round(Float.parseFloat(event.getAttributeByName(new QName("light_area")).getValue()), 4);

                    //guard against division by 0.  If we do this, light area 100 is arbitrary
                    boolean lightWas0 = false;
                    if (lightArea == 0)
                    {
                        lightWas0 = true;
                        lightArea = 100f;
                    }
                    float newHeavyArea = (float) Rounder.round(lightArea / currentScanNewRatio, 4);

                    _log.debug("Changing ratio to " + currentScanNewRatio + ".  Light=" + lightArea + ", heavy=" + newHeavyArea + ".  Light was 0? " + lightWas0);


                    SimpleStartElement newEvent = new SimpleStartElement(qname.getLocalPart());
                    Iterator<Attribute> attIter = event.getAttributes();
                    while (attIter.hasNext())
                    {
                        Attribute oldAttr = attIter.next();
                        String attrName = oldAttr.getName().getLocalPart();
                        if (attrName.equals("decimal_ratio"))
                        {
                            newEvent.addAttribute("decimal_ratio", currentScanNewRatio);
                        }
                        else if (attrName.equals("heavy_area"))
                        {
                            newEvent.addAttribute("heavy_area", newHeavyArea);                            
                        }
                        else if (attrName.equals("light_area"))
                        {
                            newEvent.addAttribute("light_area", lightArea);
                        }
                        else if (attrName.equals("heavy2light_ratio"))
                        {
                            newEvent.addAttribute("heavy2light_ratio", 1.0f / currentScanNewRatio);                            
                        }
                        else
                            newEvent.addAttribute(attrName, oldAttr.getValue());
                    }
                    eventToAdd = newEvent.getEvent();
                }
            }

            add(eventToAdd);
        }

        /**
         * Pop out of skipping mode
         * @param event
         * @throws XMLStreamException
         */
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
                insideScanWithRatioChange = false;
            }
        }
    }


    /**
     * Display help from static help file
     */
    public static class HelpAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            try
            {
                JDialog dialog = HtmlViewerPanel.showResourceInDialog(
                                "org/fhcrc/cpl/viewer/quant/gui/qurate_help.html", "Qurate Help");
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Error displaying help", e);
            }
        }
    }

    /**
     * "Create Charts" in the File menu. Create charts for running this on, by invoking the
     * 'quantchart' command programmatically, letting user specify args
     */
    protected class CreateChartsAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            Thread createChartsThread = new Thread(new Runnable()
            {
                public void run()
                {
                    PeptideQuantVisualizationCLM createChartsModule = new PeptideQuantVisualizationCLM();
                    ViewerInteractiveModuleFrame interactFrame =
                            new ViewerInteractiveModuleFrame(createChartsModule, true, null);
                    interactFrame.setUserManualGenerator(new ViewerUserManualGenerator());
                    boolean shouldExecute = interactFrame.collectArguments();

                    if (shouldExecute)
                    {

                        try
                        {
                            setMessage("Building charts.  This could take a while.  Details on command line.");
                            contentPanel.updateUI();
                            createChartsModule.execute();                              
                            infoMessage("Saved charts.  Opening summary file " +
                                    createChartsModule.getOutTsvFile().getAbsolutePath());
                            contentPanel.updateUI();
                            try
                            {
                                displayQuantFile(createChartsModule.getOutTsvFile());
                            }
                            catch (IOException e)
                            {
                                ApplicationContext.errorMessage("Failed to open quantitation file " +
                                        quantFile.getAbsolutePath(),e);
                            }                }
                        catch (Exception e)
                        {
                            String message = "Error creating charts: " + e.getMessage();

                            errorMessage(message,e);
                        }


                    }
                }
            });
            createChartsThread.start();
        }
    }

    /**
     * Open a tsv file
     */
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
            if (quantFile != null)
                wfc.setSelectedFile(quantFile);
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
                errorMessage("Failed to open quantitation file " +
                        quantFile.getAbsolutePath(),e);
            }
        }
    }

    /**
     * Save changes back to the file we opened
     * TODO: allow saving to different file?
     */
    protected class SaveAction extends AbstractAction
    {
        protected Component parentComponent;

        public SaveAction(Component parentComponent)
        {
            this.parentComponent = parentComponent;
        }

        public void actionPerformed(ActionEvent event)
        {
            WorkbenchFileChooser wfc = new WorkbenchFileChooser();
            if (quantFile != null)
                wfc.setSelectedFile(quantFile);
            wfc.setDialogTitle("Output File");
            int chooserStatus = wfc.showOpenDialog(parentComponent);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            quantFile = wfc.getSelectedFile();
            try
            {
                QuantEvent.saveQuantEventsToTSV(quantEvents, quantFile, true, true);
                setMessage("Saved changes to file " + quantFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                errorMessage("ERROR: failed to save file " + quantFile.getAbsolutePath(), e);
            }
        }

    }

    public void showProteinQuantSummaryFrame(List<ProtXmlReader.Protein> proteins)
    {
        List<QuantEvent> selectedQuantEvents = null;

        if (quantSummaryFrame != null)
            quantSummaryFrame.dispose();

        File outFile = settingsCLM.outFile;
        if (outFile == null)
        {
            //We may not actually want to keep the output file, in which case we need a temp file            
            outFile = TempFileManager.createTempFile("qurate_ProteinSelectedActionListener.tsv",
                "DUMMY_ProteinSelectedActionListener_CALLER");
        }

        try
        {
            quantSummaryFrame =
                    new ProteinQuantSummaryFrame(settingsCLM.outDir, settingsCLM.mzXmlDir, outFile,
                            settingsCLM.appendOutput);
            quantSummaryFrame.setExistingQuantEvents(quantEvents);
            quantSummaryFrame.displayData(settingsCLM.pepXmlFile, proteins);

            quantSummaryFrame.setModal(true);
            quantSummaryFrame.setVisible(true);

            selectedQuantEvents = quantSummaryFrame.getSelectedQuantEvents();
            quantSummaryFrame.dispose();
        }
        catch (IllegalArgumentException e)
        {
            infoMessage(e.getMessage());
            return;
        }
        finally
        {
            if (quantSummaryFrame != null)
                quantSummaryFrame.dispose();
        }
        //will have no effect if temp file not created
        TempFileManager.deleteTempFiles("DUMMY_ProteinSelectedActionListener_CALLER");

        if (selectedQuantEvents == null ||
                selectedQuantEvents.isEmpty())
            return;

        setMessage(selectedQuantEvents.size() + " events selected for charts");
        if (quantEvents == null)
            quantEvents = new ArrayList<QuantEvent>();
        quantEvents.addAll(selectedQuantEvents);
        eventSummaryTable.setEvents(selectedQuantEvents);
        displayedEventIndex = quantEvents.size() - selectedQuantEvents.size();
        displayCurrentQuantEvent();
    }

    /**
     * Clean up the windows that might be open
     */
    public void dispose()
    {
        _log.debug("QuantitationReviewer.dispose()");
        if (quantSummaryFrame != null)
            quantSummaryFrame.dispose();
        if (proteinSummarySelector != null)
            proteinSummarySelector.dispose();
        super.dispose();
    }


    protected class ProteinSelectedActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
            if (proteinSummarySelector.getSelectedProteins() != null &&
                    !proteinSummarySelector.getSelectedProteins().isEmpty())
            {
                showProteinQuantSummaryFrame(proteinSummarySelector.getSelectedProteins());
            }
        }
    }

    protected class ProteinSummaryAction extends AbstractAction
    {
        List<QuantEvent> selectedQuantEvents = null;

        protected Component parentComponent;

        public ProteinSummaryAction(Component parentComponent)
        {
            this.parentComponent = parentComponent;
        }


        public void actionPerformed(ActionEvent event)
        {
            ViewerInteractiveModuleFrame interactFrame =
                    new ViewerInteractiveModuleFrame(settingsCLM, true, null);
            interactFrame.setModal(true);
            interactFrame.setTitle("Settings");
            interactFrame.setUserManualGenerator(new ViewerUserManualGenerator());
            boolean hasRunSuccessfully = interactFrame.collectArguments();
            interactFrame.dispose();
            if (!hasRunSuccessfully) return;

            if (settingsCLM.proteins != null)
            {
                showProteinQuantSummaryFrame(settingsCLM.proteins);
            }
            else
            {
                if (proteinSummarySelector == null || !proteinSummarySelector.isVisible())
                    try
                    {
                        proteinSummarySelector = new ProteinSummarySelectorFrame();
                        proteinSummarySelector.setMinProteinProphet(settingsCLM.minProteinProphet);
                        proteinSummarySelector.setMinRatio(settingsCLM.minRatio);
                        proteinSummarySelector.setMaxRatio(settingsCLM.maxRatio);
                        proteinSummarySelector.setProteinGeneMap(settingsCLM.proteinGeneListMap);
                        proteinSummarySelector.addSelectionListener(new ProteinSelectedActionListener());
                        proteinSummarySelector.displayProteins(settingsCLM.protXmlFile);
                        proteinSummarySelector.setVisible(true);
                    }
                    catch (Exception e)
                    {
                        errorMessage("Error opening ProtXML file " + settingsCLM.protXmlFile.getAbsolutePath(),e);
                    }
            }
        }
    }


    /**
     * Action to quit
     */
    protected class ExitAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            setVisible(false);
            dispose();
        }
    }

    /**
     * Action to remove bad events and IDs from the file they came from.  Choose source and output pepXML
     * file
     */
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

    public File getQuantFile()
    {
        return quantFile;
    }

    public void setQuantFile(File quantFile)
    {
        this.quantFile = quantFile;
    }

    /**
     * Displays the splash screen.  Side effect: sets the splashFrame variable, so it can be disposed later
     */
    protected void showSplashScreen()
    {
        splashFrame = new SplashFrame(splashImageURL);
        splashFrame.setVisible(true);
    }

     /**
      * display the properties for the selected event, if only one's selected
      */
     public class EventSummaryTableListSelectionHandler implements ListSelectionListener
     {
         public void valueChanged(ListSelectionEvent e)
         {
             if (!e.getValueIsAdjusting())
             {
                 if (eventSummaryTable.getSelectedIndex() >= 0)
                 {
                     displayedEventIndex = eventSummaryTable.getSelectedIndex();
//System.err.println(quantEvents.get(displayedEventIndex));
                     displayCurrentQuantEvent();
                 }
             }
         }
     }

}