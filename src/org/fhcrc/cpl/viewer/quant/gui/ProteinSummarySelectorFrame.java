/*
 * Copyright (c) 2003-2012 Fred Hutchinson Cancer Research Center
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

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.SwingUtils;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.xml.stream.XMLStreamException;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 *  A window displaying a table with info on all the quantitated proteins in a protXML file.  A single
 * row, or multiple may be selected and the quantitative events pulled up in a ProteinQuantSummaryFrame.
 *
 * Gene information for proteins can optionally be pulled from a protein-gene mapping file, in which case the
 * Gene column of the table will be populated.
 *
 * The table is sortable on all columns
 */
public class ProteinSummarySelectorFrame extends JFrame
{
    protected static Logger _log = Logger.getLogger(ProteinSummarySelectorFrame.class);

    //dimensions
    protected int width = 750;
    protected int height = 800;

    //This is hacky.  It's for sizing the window appropriately when we know the number of table rows.  There
    //is probably a much more reasonable way to do this in AWT.
    protected final int TITLEBAR_HEIGHT = 55;
    protected final int STATUSPANEL_HEIGHT = 25;
    protected final int SUMMARYPANEL_HEIGHT = 81;
    protected int LOGRATIO_HISTOGRAM_PANEL_HEIGHT = 150;
    protected final int TABLEROW_HEIGHT = 17;

    //the main table
    protected ProteinSummaryTable proteinSummaryTable;
    protected ListSelectionModel tableSelectionModel;

    //Track the selected protein
    protected List<ProtXmlReader.Protein> selectedProteins = null;

    //All the proteins in the table
    protected java.util.List<ProtXmlReader.Protein> proteins;

    //map from proteins to protein group numbers, so that you san sort on group number
    protected Map<ProtXmlReader.Protein, Integer> proteinGroupNumberMap;

    //ProteinProphet cutoff for display in the table
    protected float minProteinProphet = 0.75f;

    //min and max ratio for display in table.  Protein must be > min OR < max, or both
    protected float minHighRatio = 0f;
    protected float maxLowRatio = 999f;

    protected boolean allowMultipleSelection = true;

    //Containers
    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;
    public PanelWithLogRatioHistAndFields logRatioHistogramPanel;
    
    //buttons
    protected JButton buttonShowEvents  = new JButton("Show Events");
    protected JButton buttonSaveTSV  = new JButton("Save Table");
    protected JButton buttonSelectedProtein  = new JButton("DUMMY");
    protected JButton buttonSelectAllVisible = new JButton("Select All Visible");
    protected JButton buttonDeselectAll = new JButton("Clear All");


//    protected PanelWithHistogram logRatioHistogram;

    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    //protein ratio filter control
    public JLabel maxLowRatioLabel;
    public JLabel minHighRatioLabel;
    public JLabel numPassingProteinsLabel;


    //Map from proteins to genes, for displaying a Gene column
    protected Map<String, List<String>> proteinGeneMap;

    public ProteinSummarySelectorFrame()
    {
        super();
    }

    public ProteinSummarySelectorFrame(boolean allowMultipleSelection)
    {                                             
        this.allowMultipleSelection = allowMultipleSelection;
    }

    public ProteinSummarySelectorFrame(File protXmlFile)
            throws XMLStreamException, FileNotFoundException
    {
        this();
        displayProteins(protXmlFile);
    }

    /**
     * Initialize GUI components
     */
    protected void initGUI()
    {
        setTitle("Protein Summary");
        setSize(width, height);

        try
        {
            Localizer.renderSwixml("org/fhcrc/cpl/viewer/quant/gui/ProteinSummarySelectorFrame.xml",this);
            assert null != contentPanel;
            setContentPane(contentPanel);     
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("error creating dialog", x);
            throw new RuntimeException(x);
        }

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 1;
        gbc.weightx = 1;

        ListenerHelper helper = new ListenerHelper(this);

        summaryPanel.setMaximumSize(new Dimension(1200, 80));

        gbc.fill = GridBagConstraints.NONE;
       
        gbc.insets = new Insets(5,5,5,5);
        buttonShowEvents.setEnabled(false);
        helper.addListener(buttonShowEvents, "buttonShowEvents_actionPerformed");
        gbc.gridwidth = 1;
        summaryPanel.add(buttonShowEvents, gbc);

        buttonSelectAllVisible.setEnabled(false);
        helper.addListener(buttonSelectAllVisible, "buttonSelectAllVisible_actionPerformed");
        summaryPanel.add(buttonSelectAllVisible, gbc);

        buttonDeselectAll.setEnabled(false);
        helper.addListener(buttonDeselectAll, "buttonDeselectAll_actionPerformed");
        summaryPanel.add(buttonDeselectAll, gbc);

        buttonSaveTSV.setEnabled(false);
        buttonSaveTSV.setToolTipText("Save the contents of this table to a tab-separated-value file");        
        helper.addListener(buttonSaveTSV, "buttonSaveTSV_actionPerformed");
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        summaryPanel.add(buttonSaveTSV, gbc);

        JButton buttonCancel = new JButton("Cancel");
        helper.addListener(buttonCancel, "buttonCancel_actionPerformed");
        gbc.gridwidth = GridBagConstraints.REMAINDER;        
        summaryPanel.add(buttonCancel, gbc);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(0,0,0,0);

        proteinSummaryTable = new ProteinSummaryTable();
        if (proteinGeneMap != null)
            proteinSummaryTable.proteinGeneMap = proteinGeneMap;
        tableSelectionModel = proteinSummaryTable.getSelectionModel();

        if (allowMultipleSelection)
            tableSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        //when a protein is selected, enable the "show events" button
        tableSelectionModel.addListSelectionListener(new ListSelectionListener()
        {
            public void valueChanged(ListSelectionEvent e)
            {
                if (!e.getValueIsAdjusting())
                {
                    ListSelectionModel lsm = (ListSelectionModel) e.getSource();
                    if (!lsm.isSelectionEmpty())
                    {
                        buttonShowEvents.setEnabled(true);
                    }
                }
            }
        });

        JScrollPane summaryTableScrollPane = new JScrollPane();
        summaryTableScrollPane.setViewportView(proteinSummaryTable);
        mainPanel.add(summaryTableScrollPane, gbc);

        logRatioHistogramPanel = new PanelWithLogRatioHistAndFields();
               logRatioHistogramPanel.setBorder(BorderFactory.createTitledBorder("Log Ratios"));
        logRatioHistogramPanel.setPreferredSize(new Dimension(width-10, 300));
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty=100;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        add(logRatioHistogramPanel, gbc);

        statusPanel = new JPanel();
        gbc.weighty=1;
        statusPanel.setPreferredSize(new Dimension(width-10, 50));
        messageLabel = new JLabel();
        statusPanel.add(messageLabel, gbc);
        add(statusPanel, gbc);        

        contentPanel.updateUI();
    }

    /**
     * Define the mapping from protein names to gene symbols
     * @param proteinGeneMap
     */
    public void setProteinGeneMap(Map<String, List<String>> proteinGeneMap)
    {
        this.proteinGeneMap = proteinGeneMap;
        if (proteinSummaryTable != null)
            proteinSummaryTable.proteinGeneMap = proteinGeneMap;
    }

    /**
     * Display all the proteins in the protXML file that pass the ProteinProphet threshold and have ratios
     * @param protXmlFile
     * @throws XMLStreamException
     * @throws FileNotFoundException
     */
    public void displayProteins(File protXmlFile)
            throws XMLStreamException, FileNotFoundException, IllegalArgumentException
    {
        initGUI();        
        proteins = new ArrayList<ProtXmlReader.Protein>();

        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);
        ProtXmlReader.ProteinGroupIterator groupIterator = protXmlReader.iterator();

        proteinGroupNumberMap = new HashMap<ProtXmlReader.Protein, Integer>();
        List<Float> proteinLogRatios = new ArrayList<Float>();
        while (groupIterator.hasNext())
        {
            ProteinGroup proteinGroup = groupIterator.next();
            for (ProtXmlReader.Protein protein : proteinGroup.getProteins())
            {
                ProtXmlReader.QuantitationRatio quantRatio = protein.getQuantitationRatio();
                if (protein.getProbability() > minProteinProphet && quantRatio != null)
                {
                    proteinLogRatios.add((float) Math.log(quantRatio.getRatioMean()));
                    proteins.add(protein);
                    proteinGroupNumberMap.put(protein, proteinGroup.getGroupNumber());
                }
            }
        }

        if (proteins.isEmpty())
            throw new IllegalArgumentException("No quantified proteins to display");

        Collections.sort(proteins, new ProteinRatioAscComparator());
        for (ProtXmlReader.Protein protein : proteins)
        {
            proteinSummaryTable.addProtein(protein, proteinGroupNumberMap.get(protein));
        }
        proteinSummaryTable.updateUI();
        updateExtremeRatioGUI();
        buttonSaveTSV.setEnabled(true);
        buttonSelectAllVisible.setEnabled(true);
        buttonDeselectAll.setEnabled(true);

        logRatioHistogramPanel.setMaxLowRatio(maxLowRatio);
        logRatioHistogramPanel.setMinHighRatio(minHighRatio);
        logRatioHistogramPanel.setLogRatios(proteinLogRatios);        
        logRatioHistogramPanel.setSize(width-5, LOGRATIO_HISTOGRAM_PANEL_HEIGHT-20);

        logRatioHistogramPanel.addRangeUpdateListener(new LogRatioHistogramListener(this));

        logRatioHistogramPanel.updateUI();
        
        height = Math.min(700, Math.max((proteins.size() + 1) * TABLEROW_HEIGHT,500) + SUMMARYPANEL_HEIGHT +
                STATUSPANEL_HEIGHT + TITLEBAR_HEIGHT + LOGRATIO_HISTOGRAM_PANEL_HEIGHT);
        setSize(width, height);
    }


    /**
     * A chart listener that picks up events indicating changes to the selected area
     */
    protected class LogRatioHistogramListener implements ActionListener
    {
        protected ProteinSummarySelectorFrame proteinSummaryFrame;

        public LogRatioHistogramListener(ProteinSummarySelectorFrame proteinSummaryFrame)
        {
            this.proteinSummaryFrame = proteinSummaryFrame;
        }
        public void actionPerformed(ActionEvent e)
        {
            PanelWithLogRatioHistAndFields logRatioHistAndFields = (PanelWithLogRatioHistAndFields) e.getSource();
            proteinSummaryFrame.minHighRatio = logRatioHistAndFields.getMinHighRatio();
            proteinSummaryFrame.maxLowRatio = logRatioHistAndFields.getMaxLowRatio();
            proteinSummaryFrame.updateExtremeRatioGUI();
        }
    }

    /**
     * In response to a user action restricting ratios, modify what's shown on the table and reflected in the labels.
     * Called from two places
     */
    protected void updateExtremeRatioGUI()
    {
        proteinSummaryTable.showOnlyExtremeRatios(maxLowRatio, minHighRatio);
    }


    /**
     * Add a listener for selecting a row in the table.  This is used for populating the selected protein
     * @param listener
     */
    public void addSelectionListener(ActionListener listener)
    {
        buttonSelectedProtein.addActionListener(listener);
    }

    /**
     * Show the detail screen with a table of the actual events for this/these protein(s)
     * @param event
     */
    public void buttonShowEvents_actionPerformed(ActionEvent event)
    {
        int[] selectedRows = proteinSummaryTable.getSelectedRows();
        if (selectedRows.length == 0)
            return;
        selectedProteins = new ArrayList<ProtXmlReader.Protein>(selectedRows.length);
        for (int selectedIndex : selectedRows)
            selectedProteins.add(proteins.get(selectedIndex));

        ActionListener[] buttonListeners = buttonSelectedProtein.getActionListeners();

        if (buttonListeners != null)
        {
            for (ActionListener listener : buttonListeners)
                listener.actionPerformed(event);
        }
    }

    /**
     * Select all currently-visible (i.e., passes filter) rows
     * @param event
     */
    public void buttonSelectAllVisible_actionPerformed(ActionEvent event)
    {      
        for (int i=0; i<proteinSummaryTable.getRowCount(); i++)
            tableSelectionModel.addSelectionInterval(i, i);
    }

    /**
     * Deselect all currently-selected rows, whether or not they're currently passing the filter
     * @param event
     */
    public void buttonDeselectAll_actionPerformed(ActionEvent event)
    {
        tableSelectionModel.clearSelection();
    }

    /**
     * Save the table contents as a TSV file
     * @param event
     */
    public void buttonSaveTSV_actionPerformed(ActionEvent event)
    {
        WorkbenchFileChooser wfc = new WorkbenchFileChooser();
        int chooserStatus = wfc.showSaveDialog(this);
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return;
        File outTsvFile = wfc.getSelectedFile();
        try
        {
            SwingUtils.SaveTableAsTSV(proteinSummaryTable, outTsvFile);
            setMessage("Saved table to file " + outTsvFile.getAbsolutePath());
        }
        catch (IOException e)
        {
            errorMessage("Failed to save file " + outTsvFile.getAbsolutePath(), e);
        }
    }

    /**
     * Quit
     * @param event
     */
    public void buttonCancel_actionPerformed(ActionEvent event)
    {
        selectedProteins = null;
        setVisible(false);
    }


    /**
     * Sort proteins by ratio, ascending
     */
    public static class ProteinRatioAscComparator implements Comparator<ProtXmlReader.Protein>
    {
        public int compare(ProtXmlReader.Protein o1, ProtXmlReader.Protein o2)
        {
            if (o1.getQuantitationRatio().getRatioMean() > o2.getQuantitationRatio().getRatioMean())
                return 1;
            if (o1.getQuantitationRatio().getRatioMean() < o2.getQuantitationRatio().getRatioMean())
                return -1;
            return 0;
        }
    }

    public void displayProteins()
    {
        proteinSummaryTable.clearRows();
        for (ProtXmlReader.Protein protein : proteins)
            proteinSummaryTable.addProtein(protein, proteinGroupNumberMap.get(protein));
        proteinSummaryTable.updateUI();
        buttonSaveTSV.setEnabled(true);                
    }

    /**
     * Sortable on all columns
     */
    public static final class ProteinSummaryTable extends JTable
    {
        protected Map<String, List<String>> proteinGeneMap;

        protected int ratioColumnIndex = 0;

        protected TableRowSorter<TableModel> sorter;

        protected static final String[] columnTitles = new String[]
                {
                        "Protein", "Group", "Genes", "Prob", "Ratio", "IDPeptides", "QuantEvents"
                };

     protected static final int[] columnWidths = new int[]
                {
                        100, 70, 170, 70, 70, 90, 90
                };

        DefaultTableModel model = new DefaultTableModel(0, columnTitles.length)
            {
                //all cells uneditable
                public boolean isCellEditable(int row, int column)
                {
                    return false;
                }

                public Class getColumnClass(int columnIndex)
                {
                    String columnTitle = columnTitles[columnIndex];
                    if ("Group".equals(columnTitle) || "IDPeptides".equals(columnTitle) ||
                            "QuantEvents".equals(columnTitle))
                        return Integer.class;
                    if ("Probability".equals(columnTitle) || "Ratio".equals(columnTitle))
                        return Float.class;
                    return String.class;
                }
            };

        public ProteinSummaryTable()
        {
            setModel(model);
            for (int i=0; i<columnTitles.length; i++)
            {
                getColumnModel().getColumn(i).setHeaderValue(columnTitles[i]);
                getColumnModel().getColumn(i).setPreferredWidth(columnWidths[i]);
            }

            getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

            sorter = new TableRowSorter<TableModel>(model);
            setRowSorter(sorter);
        }

        /**
         * Filters rows based on having extreme ratios
         */
        protected class RatioRowFilter extends RowFilter<TableModel, Object>
        {
            protected float maxLowRatioValue;
            protected float minHighRatioValue;

            public RatioRowFilter(float maxLowRatioValue, float minHighRatioValue)
            {
                this.maxLowRatioValue = maxLowRatioValue;
                this.minHighRatioValue = minHighRatioValue;
            }
            public boolean include(RowFilter.Entry entry)
            {
                float ratio = (Float) entry.getValue(ratioColumnIndex);
                return isRatioIncluded(ratio);
            }
            public boolean isRatioIncluded(float ratio)
            {
                return (ratio <= maxLowRatioValue || ratio >= minHighRatioValue);
            }
        }

        public void showOnlyExtremeRatios(float maxLowRatioValue, float minHighRatioValue)
        {
            RowFilter<TableModel, Object> rf = new RatioRowFilter(maxLowRatioValue, minHighRatioValue);
            sorter.setRowFilter(rf);
        }

        /**
         * Override parent behavior, which doesn't take sorting into account
         * @return
         */
        public int[] getSelectedRows()
        {
            int[] rawRows = super.getSelectedRows();
            for (int i=0; i<rawRows.length; i++)
            {
                rawRows[i] = convertRowIndexToModel(rawRows[i]);
            }
            return rawRows;
        }

        /**
         * Remove all properties from table
         */
        public void clearRows()
        {
            while (model.getRowCount() > 0)
            {
                model.removeRow(0);
            }
        }

        public void addProtein(ProtXmlReader.Protein protein, int proteinGroupNumber)
        {
            int numRows = model.getRowCount();
            model.setRowCount(numRows + 1);
            int currentColIndex = 0;
            model.setValueAt(protein.getProteinName(), numRows, currentColIndex++);
            model.setValueAt(proteinGroupNumber, numRows, currentColIndex++);
            if (proteinGeneMap != null && proteinGeneMap.containsKey(protein.getProteinName()))
            {
                List<String> genes = proteinGeneMap.get(protein.getProteinName());
                if (genes != null && !genes.isEmpty())
                {
                    StringBuffer genesStringBuf = new StringBuffer(genes.get(0));
                    for (int i=1; i<genes.size(); i++)
                        genesStringBuf.append("," + genes.get(i));
                    model.setValueAt(genesStringBuf.toString(), numRows, currentColIndex++);
                }
            }
            else model.setValueAt("", numRows, currentColIndex++);
            model.setValueAt(protein.getProbability(), numRows, currentColIndex++);
            ratioColumnIndex = currentColIndex;
            model.setValueAt(protein.getQuantitationRatio().getRatioMean(), numRows, currentColIndex++);
            model.setValueAt(protein.getUniquePeptidesCount(), numRows, currentColIndex++);            
            model.setValueAt(protein.getQuantitationRatio().getRatioNumberPeptides(), numRows, currentColIndex++);
        }
    }

    public List<ProtXmlReader.Protein> getSelectedProteins()
    {
        return selectedProteins;
    }

    public void setSelectedProteins(List<ProtXmlReader.Protein> selectedProteins)
    {
        this.selectedProteins = selectedProteins;
    }

    public float getMinProteinProphet()
    {
        return minProteinProphet;
    }

    public void setMinProteinProphet(float minProteinProphet)
    {
        this.minProteinProphet = minProteinProphet;
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

    public float getMinHighRatio()
    {
        return minHighRatio;
    }

    public void setMinHighRatio(float minHighRatio)
    {
        this.minHighRatio = minHighRatio;
    }

    public float getMaxLowRatio()
    {
        return maxLowRatio;
    }

    public void setMaxLowRatio(float maxLowRatio)
    {
        this.maxLowRatio = maxLowRatio;
    }
}
