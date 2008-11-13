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


    protected int fullWidth = 900;
    protected int fullHeight = 600;    

    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public JScrollPane eventsScrollPane;
    public JPanel eventsPanel;

    JLabel proteinNameLabel = new JLabel("Protein: ");
    JLabel proteinRatioLabel  = new JLabel("Ratio: ");
    JButton buildChartsForSelectedButton = new JButton("Build Selected Charts");

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

    public ProteinQuantSummaryFrame(File protXmlFile, File pepXmlFile, String proteinName, File outDir, File mzXmlDir)
    {
        this();
        this.outDir = outDir;
        this.mzXmlDir = mzXmlDir;
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
        Collections.sort(quantEvents, new QuantEventInfo.RatioAscComparator());
        displayEvents();
    }


    protected void initGUI()
    {
        //Global stuff
        setSize(fullWidth, fullHeight);

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
        gbc.fill = GridBagConstraints.BOTH;


        eventsScrollPane = new JScrollPane();
        eventsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        eventsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        eventsPanel = new JPanel();
        eventsPanel.setLayout(new GridBagLayout());

       //Properties panel stuff
        eventsTableModel = new DefaultTableModel(0, 7)
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
        eventsTable = new JTable(eventsTableModel);

        TableColumn peptideColumn = eventsTable.getColumnModel().getColumn(1);
        peptideColumn.setPreferredWidth(170);
        peptideColumn.setMinWidth(140);


        TableColumn checkboxColumn = eventsTable.getColumnModel().getColumn(0);
        checkboxColumn.setHeaderRenderer(new CheckBoxHeader(new SelectAllListener()));
//        checkboxColumn.setCellEditor(new DefaultCellEditor(new JCheckBox()));
        checkboxColumn.setPreferredWidth(20);
        checkboxColumn.setMaxWidth(20);

        eventsTable.getColumnModel().getColumn(1).setHeaderValue("Peptide");
        eventsTable.getColumnModel().getColumn(2).setHeaderValue("Probability");
        eventsTable.getColumnModel().getColumn(3).setHeaderValue("Ratio");
        eventsTable.getColumnModel().getColumn(4).setHeaderValue("Light");
        eventsTable.getColumnModel().getColumn(5).setHeaderValue("Heavy");
        eventsTable.getColumnModel().getColumn(6).setHeaderValue("LogRatio");


        TableColumn logRatioSliderColumn = eventsTable.getColumnModel().getColumn(6);
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
            done = true;
            return;
        }

        ApplicationContext.infoMessage(selectedQuantEvents.size() + " events selected for charts");

        QuantitationVisualizer quantVisualizer = new QuantitationVisualizer();
        quantVisualizer.setMzXmlDir(mzXmlDir);
        quantVisualizer.setOutDir(outDir);
        quantVisualizer.setOutTsvFile(new File(outDir, "quantitation_" + proteinName + ".tsv"));
        quantVisualizer.setOutHtmlFile(new File(outDir, "quantitation_" + proteinName + ".htnl"));

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

        for (QuantEventInfo quantEvent : quantEvents)
        {
//            JPanel eventPanel = createPanelForEvent(quantEvent);
//            eventsPanel.add(eventPanel, gbc);

            int numRows = eventsTableModel.getRowCount();
            eventsTableModel.setRowCount(numRows + 1);
//            JCheckBox thisEventCheckBox = new JCheckBox();
            eventsTableModel.setValueAt(new Boolean(false), numRows, 0);
            eventsTableModel.setValueAt(quantEvent.getPeptide(), numRows, 1);              
            eventsTableModel.setValueAt("" + quantEvent.getPeptideProphet(), numRows, 2);
            eventsTableModel.setValueAt("" + quantEvent.getRatio(), numRows, 3);
            eventsTableModel.setValueAt("" + quantEvent.getLightIntensity(), numRows, 4);
            eventsTableModel.setValueAt("" + quantEvent.getHeavyIntensity(), numRows, 5);

            float ratioBound = 10f;
            float logRatioBounded =
                    (float) Math.log(Math.min(ratioBound, Math.max(1.0f / ratioBound, quantEvent.getRatio())));
            int logRatioIntegerizedHundredScale =
                    (int) (logRatioBounded * 100 / (2 * Math.log(ratioBound))) + 50;
System.err.println("ratio: " + quantEvent.getRatio() + ", log: " + logRatioBounded + ", integerized: " + logRatioIntegerizedHundredScale);
            eventsTableModel.setValueAt(logRatioIntegerizedHundredScale, numRows, 6);

        }
        buildChartsForSelectedButton.setEnabled(true);

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
}
