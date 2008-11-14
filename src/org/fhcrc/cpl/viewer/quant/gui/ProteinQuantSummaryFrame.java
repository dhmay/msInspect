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
public class ProteinQuantSummaryFrame extends JFrame
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
    DefaultTableModel eventsTableModel;
    JTable eventsTable;

    public ProteinQuantSummaryFrame()
    {
        initGUI();
    }

    public ProteinQuantSummaryFrame(File protXmlFile, File pepXmlFile, String proteinName,
                                    File outDir, File mzXmlDir, File outFile, boolean appendOutput)
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
            infoMessage("Protein " + proteinName + " does not occur in the file" +
                    " or has no quantitation");
            return;
        }

        ProtXmlReader.QuantitationRatio quantRatio = protein.getQuantitationRatio();

        proteinRatioLabel.setText("Ratio: " + quantRatio.getRatioMean());
        contentPanel.updateUI();

        List<String> proteinPeptidesRatiosUsed = quantRatio.getPeptides();

        ApplicationContext.infoMessage("Protein ratio: " + quantRatio.getRatioMean() + ", Std Dev: " +
                quantRatio.getRatioStandardDev());

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
            setIconImage(ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif")));
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

        //summary panel
        summaryPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        gbc.fill = GridBagConstraints.NONE;        
        summaryPanel.add(proteinNameLabel, gbc);
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        summaryPanel.add(proteinRatioLabel, gbc);
        buildChartsForSelectedButton.setEnabled(false);
        helper.addListener(buildChartsForSelectedButton, "buttonBuildCharts_actionPerformed");
        summaryPanel.add(buildChartsForSelectedButton, gbc);
        helper.addListener(showPropertiesButton, "buttonShowProperties_actionPerformed");
        showPropertiesButton.setEnabled(false);
        summaryPanel.add(showPropertiesButton, gbc);

        gbc.fill = GridBagConstraints.BOTH;


        eventsScrollPane = new JScrollPane();
        eventsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        eventsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        eventsPanel = new JPanel();
        eventsPanel.setLayout(new GridBagLayout());

       //Properties panel stuff
        eventsTableModel = new DefaultTableModel(0, 8)
        {
            //all cells uneditable
            public boolean isCellEditable(int row, int column)
            {
                if (column == 0)
                    return true;
                return false;
            }

            public Class getColumnClass(int columnIndex)
            {
                switch (columnIndex)
                {
                    case 0:
                        return Boolean.class;
                    case 6:
                        return JSlider.class;
                    default:
                        return String.class;
                }
            }
        };
        eventsTable = new JTable(eventsTableModel)
        {
            protected Color altRowColor = new Color(235, 235, 235);
            /**
             * Shades alternate peptides in different colors.
             */
            public Component prepareRenderer(TableCellRenderer renderer, int row, int column)
            {
                Component c = super.prepareRenderer(renderer, row, column);
                if (isCellSelected(row, column) == false)
                {
                    Color rowColor = UIManager.getColor("Table.background");
                    if (shadedTableRows.contains(row))
                        rowColor = altRowColor;
                    c.setBackground(rowColor);
                    c.setForeground(UIManager.getColor("Table.foreground"));
                } else
                {
                    c.setBackground(UIManager.getColor("Table.selectionBackground"));
                    c.setForeground(UIManager.getColor("Table.selectionForeground"));
                }
                return c;
            }
        };
        ListSelectionModel tableSelectionModel = eventsTable.getSelectionModel();
        tableSelectionModel.addListSelectionListener(new EventsTableListSelectionHandler());

        TableColumn peptideColumn = eventsTable.getColumnModel().getColumn(1);
        peptideColumn.setPreferredWidth(170);
        peptideColumn.setMinWidth(140);


        TableColumn checkboxColumn = eventsTable.getColumnModel().getColumn(0);
        checkboxColumn.setHeaderRenderer(new CheckBoxHeader(new SelectAllListener()));
//        checkboxColumn.setCellEditor(new DefaultCellEditor(new JCheckBox()));
        checkboxColumn.setPreferredWidth(20);
        checkboxColumn.setMaxWidth(20);

        eventsTable.getColumnModel().getColumn(1).setHeaderValue("Peptide");
        eventsTable.getColumnModel().getColumn(2).setHeaderValue("Charge");
        eventsTable.getColumnModel().getColumn(3).setHeaderValue("Probability");
        eventsTable.getColumnModel().getColumn(4).setHeaderValue("Ratio");
        eventsTable.getColumnModel().getColumn(5).setHeaderValue("Light");
        eventsTable.getColumnModel().getColumn(6).setHeaderValue("Heavy");
        eventsTable.getColumnModel().getColumn(7).setHeaderValue("LogRatio");


        TableColumn logRatioSliderColumn = eventsTable.getColumnModel().getColumn(7);
        JSliderRenderer sliderRenderer = new JSliderRenderer();
        logRatioSliderColumn.setCellRenderer(sliderRenderer);
        logRatioSliderColumn.setPreferredWidth(280);
        logRatioSliderColumn.setMinWidth(100);

        eventsScrollPane.setViewportView(eventsTable);
        eventsScrollPane.setMinimumSize(new Dimension(400, 400));

        gbc.insets = new Insets(0,0,0,0);             
        mainPanel.add(eventsScrollPane, gbc);

        setVisible(true);
    }

    public void buttonShowProperties_actionPerformed(ActionEvent event)
    {
        eventPropertiesDialog.setVisible(true);
    }

    public void buttonBuildCharts_actionPerformed(ActionEvent event)
    {
        selectedQuantEvents = new ArrayList<QuantEventInfo>();
        for (int i=0; i<eventsTableModel.getRowCount(); i++)
        {
            Boolean isSelected = (Boolean) eventsTableModel.getValueAt(i, 0);
            if (isSelected)
                selectedQuantEvents.add(quantEvents.get(i));
        }
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

        done = true;
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

        shadedTableRows = new ArrayList<Integer>();
        boolean shaded = true;
        String previousPeptide = "";
        for (QuantEventInfo quantEvent : quantEvents)
        {
            if (!previousPeptide.equals(quantEvent.getPeptide()))
            {
                shaded = !shaded;
                previousPeptide = quantEvent.getPeptide();
            }

            int numRows = eventsTableModel.getRowCount();

            if (shaded)
                shadedTableRows.add(numRows);

            eventsTableModel.setRowCount(numRows + 1);
//            JCheckBox thisEventCheckBox = new JCheckBox();
            eventsTableModel.setValueAt(new Boolean(false), numRows, 0);
            eventsTableModel.setValueAt(quantEvent.getPeptide(), numRows, 1);
            eventsTableModel.setValueAt("" + quantEvent.getCharge(), numRows, 2);
            eventsTableModel.setValueAt("" + quantEvent.getPeptideProphet(), numRows, 3);
            eventsTableModel.setValueAt("" + quantEvent.getRatio(), numRows, 4);
            eventsTableModel.setValueAt("" + quantEvent.getLightIntensity(), numRows, 5);
            eventsTableModel.setValueAt("" + quantEvent.getHeavyIntensity(), numRows, 6);

            float ratioBound = 10f;
            float logRatioBounded =
                    (float) Math.log(Math.min(ratioBound, Math.max(1.0f / ratioBound, quantEvent.getRatio())));
            int logRatioIntegerizedHundredScale =
                    (int) (logRatioBounded * 100 / (2 * Math.log(ratioBound))) + 50;
            eventsTableModel.setValueAt(logRatioIntegerizedHundredScale, numRows, 7);

        }
        buildChartsForSelectedButton.setEnabled(true);
        showPropertiesButton.setEnabled(true);


        contentPanel.updateUI();
    }


    class CheckBoxHeader extends JCheckBox
            implements TableCellRenderer, MouseListener
    {
        protected CheckBoxHeader rendererComponent;
        protected int column;
        protected boolean mousePressed = false;
        public CheckBoxHeader(ItemListener itemListener) {
            rendererComponent = this;
            rendererComponent.addItemListener(itemListener);
        }
        public Component getTableCellRendererComponent(
                JTable table, Object value,
                boolean isSelected, boolean hasFocus, int row, int column) {
            if (table != null) {
                JTableHeader header = table.getTableHeader();
                if (header != null) {
                    rendererComponent.setForeground(header.getForeground());
                    rendererComponent.setBackground(header.getBackground());
                    rendererComponent.setFont(header.getFont());
                    header.addMouseListener(rendererComponent);
                }
            }
            setColumn(column);
            rendererComponent.setText("Check All");
            setBorder(UIManager.getBorder("TableHeader.cellBorder"));
            return rendererComponent;
        }
        protected void setColumn(int column) {
            this.column = column;
        }
        public int getColumn() {
            return column;
        }
        protected void handleClickEvent(MouseEvent e) {
            if (mousePressed) {
                mousePressed=false;
                JTableHeader header = (JTableHeader)(e.getSource());
                JTable tableView = header.getTable();
                TableColumnModel columnModel = tableView.getColumnModel();
                int viewColumn = columnModel.getColumnIndexAtX(e.getX());
                int column = tableView.convertColumnIndexToModel(viewColumn);

                if (viewColumn == this.column && e.getClickCount() == 1 && column != -1) {
                    doClick();
                }
            }
        }
        public void mouseClicked(MouseEvent e) {
            handleClickEvent(e);
            ((JTableHeader)e.getSource()).repaint();
        }
        public void mousePressed(MouseEvent e) {
            mousePressed = true;
        }
        public void mouseReleased(MouseEvent e) {
        }
        public void mouseEntered(MouseEvent e) {
        }
        public void mouseExited(MouseEvent e) {
        }
    }

    class SelectAllListener implements ItemListener
    {
        public void itemStateChanged(ItemEvent e) {
            Object source = e.getSource();
            if (source instanceof AbstractButton == false) return;
            boolean checked = e.getStateChange() == ItemEvent.SELECTED;
            for(int x = 0, y = eventsTable.getRowCount(); x < y; x++)
            {
                eventsTable.setValueAt(new Boolean(checked),x,0);
            }
        }
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


    public class JSliderRenderer implements TableCellRenderer
    {
        protected JSlider slider = null;

        public JSliderRenderer()
        {
            slider = new JSlider();
            slider.setMinimum(0);
            slider.setMaximum(100);
            slider.setPaintLabels(false);
            slider.setPaintTicks(false);
            slider.setMajorTickSpacing(25);
            slider.setPreferredSize(new Dimension(280, 15));
            slider.setPreferredSize(new Dimension(100, 15));
            slider.setToolTipText("Log ratio, bounded at 0.1 and 10");
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column)
        {
            Integer val = (Integer)value;
            slider.setValue(val.intValue());
            return slider;
        }
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
                ListSelectionModel lsm = (ListSelectionModel) e.getSource();
                if (lsm.isSelectionEmpty()) {
                    eventPropertiesTable.clearProperties();
                } else {
                    // Find out which indexes are selected.
                    int minIndex = lsm.getMinSelectionIndex();
                    int maxIndex = lsm.getMaxSelectionIndex();
                    if (minIndex == maxIndex)
                        eventPropertiesTable.displayQuantEvent(quantEvents.get(minIndex));
                    else
                        eventPropertiesTable.clearProperties();
                }
            }
        }
    }
}
