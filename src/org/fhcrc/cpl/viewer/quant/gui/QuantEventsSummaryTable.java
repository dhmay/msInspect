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

import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.fhcrc.cpl.viewer.quant.QuantEventAssessor;
import org.fhcrc.cpl.toolbox.Rounder;

import javax.swing.*;
import javax.swing.table.*;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.event.*;

/**
 *  Display quantitative event info in a table.
 *
 * Ratio is indicated by a number, and also by a slider indicating the ratio in log space
 */
public class QuantEventsSummaryTable extends JTable
{
    //shading for alternating peptides.  Not a great idea if the table is re-sorted
    protected List<Integer> shadedTableRows = new ArrayList<Integer>();

    //List of events that have already been selected, should not be allowed to be deselected
    protected List<Integer> alreadySelectedRows = new ArrayList<Integer>();

    protected List<QuantEvent> quantEvents = new ArrayList<QuantEvent>();

    protected Map<String, Integer> fractionNameNumberMap = new HashMap<String, Integer>();

    protected TableColumn logRatioSliderColumn;
    protected TableColumn proteinColumn;
    protected TableColumn fractionColumn;
    protected TableColumn scanColumn;
    protected TableColumn assessmentColumn;
    protected TableColumn geneColumn;

    protected Map<String, List<String>> proteinGenesMap;

    protected int quantCurationColumnIndex;
    protected int quantAlgAssessmentColumnIndex;



//    protected QuantEventChangeListener changeListener = new QuantEventChangeListener();

    protected int ratioColumnIndex = 0;

    protected TableRowSorter<TableModel> sorter;




    DefaultTableModel model = new DefaultTableModel(0, 13)
    {
        //all cells uneditable
        public boolean isCellEditable(int row, int column)
        {
            if (column == 0)
            {
                if (alreadySelectedRows == null || !alreadySelectedRows.contains(row))
                    return true;
            }
            return false;
        }

        public Class getColumnClass(int columnIndex)
        {
//            String columnName = getColumnName(columnIndex); //(String) getColumn(columnIndex).getHeaderValue();*///
            if (columnIndex < 0)
                return String.class;
            int viewColumnIndex = convertColumnIndexToView(columnIndex);
            if (viewColumnIndex < 0)
                return String.class;
            TableColumn column =
                    getColumnModel().getColumn(viewColumnIndex);
            String columnName = (String) column.getHeaderValue();
            if ("Ratio".equals(columnName))// || "Light".equals(columnName) || "Heavy".equals(columnName))
                return Float.class;
            else if ("LogRatio".equals(columnName))
                return JSlider.class;
            else if ("Scan".equals(columnName))
                return Integer.class;
            else return super.getColumnClass(columnIndex);

//            else return String.class;
        }
    };


    /**
     * Hide the Protein column.  There's no undoing this
     */
    public void hideProteinColumn()
    {
        this.removeColumn(proteinColumn);
//        proteinColumn.setMaxWidth(0);
    }

    /**
     * Hide the fraction column.  There's no undoing this
     */
    public void hideFractionColumn()
    {
        this.removeColumn(fractionColumn);

    }

    /**
     * Hide the fraction column.  There's no undoing this
     */
    public void hideScanColumn()
    {
        this.removeColumn(scanColumn);
    }

    /**
     * Hide the Assessment column.  There's no undoing this
     */
    public void hideAssessmentColumn()
    {
        this.removeColumn(assessmentColumn);
//        assessmentColumn.setMaxWidth(0);
    }

    /**
     * Hide the Gene column.  There's no undoing this
     */
    public void hideGeneColumn()
    {
        this.removeColumn(geneColumn);
//        geneColumn.setMaxWidth(0);
    }

    public QuantEventsSummaryTable()
    {
        setModel(model);
        sorter = new TableRowSorter<TableModel>(model);


        int columnNum = 0;

        List<String> columnNames = new ArrayList<String>();

        geneColumn = getColumnModel().getColumn(columnNum++);
        geneColumn.setHeaderValue("Gene");
        columnNames.add("Gene");
        geneColumn.setPreferredWidth(90);

        proteinColumn = getColumnModel().getColumn(columnNum++);
        proteinColumn.setHeaderValue("Protein");
        columnNames.add("Protein");
        proteinColumn.setPreferredWidth(90);

        TableColumn peptideColumn = getColumnModel().getColumn(columnNum++);
        peptideColumn.setHeaderValue("Peptide");
        columnNames.add("Peptide");
        peptideColumn.setPreferredWidth(170);
        peptideColumn.setMinWidth(140);

        fractionColumn = getColumnModel().getColumn(columnNum++);
        fractionColumn.setHeaderValue("Fraction");
        columnNames.add("Fraction");

        scanColumn = getColumnModel().getColumn(columnNum++);
        scanColumn.setHeaderValue("Scan");
        columnNames.add("Scan");

        getColumnModel().getColumn(columnNum).setHeaderValue("Charge");
        columnNames.add("Charge");
        getColumnModel().getColumn(columnNum++).setPreferredWidth(45);
        getColumnModel().getColumn(columnNum).setHeaderValue("Prob");
        columnNames.add("Prob");
        
        getColumnModel().getColumn(columnNum++).setPreferredWidth(50);
        getColumnModel().getColumn(columnNum).setHeaderValue("Ratio");
        columnNames.add("Ratio");

        ratioColumnIndex = columnNum;
        getColumnModel().getColumn(columnNum++).setPreferredWidth(50);
        getColumnModel().getColumn(columnNum++).setHeaderValue("Light");
        columnNames.add("Light");

        getColumnModel().getColumn(columnNum++).setHeaderValue("Heavy");
        columnNames.add("Heavy");


        logRatioSliderColumn = getColumnModel().getColumn(columnNum);
        logRatioSliderColumn.setHeaderValue("LogRatio");
        columnNames.add("LogRatio");

        JSliderRenderer sliderRenderer = new JSliderRenderer();
        logRatioSliderColumn.setCellRenderer(sliderRenderer);
        logRatioSliderColumn.setPreferredWidth(280);
        logRatioSliderColumn.setMinWidth(100);
        //special comparator for slider column
        sorter.setComparator(columnNum++, new Comparator<Integer>() {
            public int compare(Integer o1, Integer o2)
            {
                return o1 > o2 ? 1 : o1 < o2 ? -1 : 0;
            }
        });
        quantAlgAssessmentColumnIndex = columnNum;
        assessmentColumn = getColumnModel().getColumn(columnNum++);
        assessmentColumn.setHeaderValue("Assessment");
        assessmentColumn.setCellRenderer(new FlagQuantStatusRenderer());
        columnNames.add("Assessment");

        quantCurationColumnIndex = columnNum;
        TableColumn evaluationColumn =  getColumnModel().getColumn(columnNum++);
        evaluationColumn.setHeaderValue("Evaluation");
        evaluationColumn.setCellRenderer(new CurationStatusRenderer());
        columnNames.add("Evaluation");


        getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);


        setRowSorter(sorter);
//        model.setColumnIdentifiers(columnNames.toArray(new String[columnNames.size()]));
    }

    public void setSelectionMode(int mode)
    {
        getSelectionModel().setSelectionMode(mode);
    }

    /**
     * Returns model, not view, index
     * @return
     */
    public int getSelectedIndex()
    {
        int[] selectedIndices = getSelectedIndices();
        if (selectedIndices == null || selectedIndices.length < 1 || selectedIndices.length > 1)
            return -1;
        return selectedIndices[0];
    }

    public int[] getSelectedIndices()
    {
        int[] rawRows = super.getSelectedRows();
        if (rawRows == null)
            return null;
        for (int i=0; i<rawRows.length; i++)
        {
            rawRows[i] = convertRowIndexToModel(rawRows[i]);
        }
        return rawRows;
    }

    protected Color altRowColor = new Color(235, 235, 235);
    /**
     * Shades alternate peptides in different colors.
     */
    public Component prepareRenderer(TableCellRenderer renderer, int row, int column)
    {
        Component c = super.prepareRenderer(renderer, row, column);
        Color previousForegroundColor = c.getForeground();
        //Need to be careful: don't recolor the foreground if it's a special color
        boolean shouldRecolorForeground = previousForegroundColor != Color.red &&
                previousForegroundColor != Color.green && previousForegroundColor != Color.blue;
        if (isCellSelected(row, column))
        {
            c.setBackground(UIManager.getColor("Table.selectionBackground"));
            if (shouldRecolorForeground)
            {
                Color selectedForegroundColor = UIManager.getColor("Table.selectionForeground");
                if (alreadySelectedRows != null && alreadySelectedRows.contains(row))
                    selectedForegroundColor = Color.GRAY;
                c.setForeground(selectedForegroundColor);
            }
        }
        else
        {
            Color rowColor = UIManager.getColor("Table.background");
            if (shadedTableRows.contains(row))
                rowColor = altRowColor;
            c.setBackground(rowColor);
            if (shouldRecolorForeground)
            {
                Color unselectedForegroundColor = UIManager.getColor("Table.foreground");
                if (alreadySelectedRows != null && alreadySelectedRows.contains(row))
                    unselectedForegroundColor = Color.GRAY;
                c.setForeground(unselectedForegroundColor);
            }
        }
        return c;
    }

    /**
     * Remove all properties from table
     */
    public void clearProperties()
    {
        while (model.getRowCount() > 0)
        {
            model.removeRow(0);
        }
    }

    /**
     * Map fraction names to numbers
     * @param events
     */
    protected void buildFractionNameNumberMap(List<QuantEvent> events)
    {
        Set<String> fractionNames = new HashSet<String>();
        for (QuantEvent quantEvent : events)
        {
            String fractionName = quantEvent.getFraction();
            if (fractionName != null)
            {
                fractionNames.add(fractionName);
            }
        }
        if (fractionNames.isEmpty())
            fractionNameNumberMap = null;
        else
        {
            fractionNameNumberMap = new HashMap<String, Integer>();
            List<String> fractionNamesList = new ArrayList<String>(fractionNames);
            Collections.sort(fractionNamesList);
            for (int i=0; i<fractionNamesList.size(); i++)
                fractionNameNumberMap.put(fractionNamesList.get(i), i+1);
        }
    }

    public void setEvents(List<QuantEvent> events)
    {
        buildFractionNameNumberMap(events);

        for (QuantEvent quantEvent : events)
        {
            addEvent(quantEvent, false);
        }
        quantEvents = new ArrayList<QuantEvent>(events);
    }

    public void addEvent(QuantEvent quantEvent, boolean alreadySelected)
    {
        String previousPeptide = "";
        int numRows = model.getRowCount();
        boolean previousRowShaded = false;

        quantEvent.addQuantCurationStatusListener(new QuantEventChangeListener(numRows));
        quantEvent.addAlgAssessmentStatusListener(new QuantEventAlgAssessmentChangeListener(numRows));
        
        
        if (numRows > 0)
        {
            previousPeptide = model.getValueAt(numRows-1, 3).toString();
            if (shadedTableRows.contains(numRows-1))
                previousRowShaded = true;
        }
        String peptide = quantEvent.getPeptide();
        boolean shaded = ((previousRowShaded && peptide.equals(previousPeptide)) ||
                (!previousRowShaded && !peptide.equals(previousPeptide)));
        if (shaded)
            shadedTableRows.add(numRows);

        model.setRowCount(numRows + 1);

        String fraction = quantEvent.getFraction();
        int fractionNum = 0;
        if (fractionNameNumberMap != null && fraction != null &&
            fractionNameNumberMap.containsKey(fraction))
            fractionNum = fractionNameNumberMap.get(fraction);



        int colNum = 0;
        String geneValue = "";
        String protein = quantEvent.getProtein();
        if (proteinGenesMap != null && proteinGenesMap.containsKey(protein))
        {
            //it would be better to do this once and cache it
            StringBuffer geneValueBuf = new StringBuffer();
            boolean first = true;
            for (String gene : proteinGenesMap.get(protein))
            {
                if (!first)
                    geneValueBuf.append(",");
                geneValueBuf.append(gene);
                first = false;
            }
            geneValue = geneValueBuf.toString();
        }
        model.setValueAt(geneValue, numRows, colNum++);
        model.setValueAt(protein, numRows, colNum++);
        model.setValueAt(peptide, numRows, colNum++);
        model.setValueAt("" + fractionNum, numRows, colNum++);
        model.setValueAt(quantEvent.getScan(), numRows, colNum++);
        model.setValueAt("" + quantEvent.getCharge(), numRows, colNum++);
        model.setValueAt("" + Rounder.round(quantEvent.getPeptideProphet(),3), numRows, colNum++);
        model.setValueAt((float) Rounder.round(quantEvent.getRatio(),3), numRows, colNum++);
        model.setValueAt("" + Rounder.round(quantEvent.getLightIntensity(),1), numRows, colNum++);
        model.setValueAt("" + Rounder.round(quantEvent.getHeavyIntensity(),1), numRows, colNum++);
        model.setValueAt(integerizeRatio(quantEvent.getRatio()), numRows, colNum++);
        String assessmentString = "";
        QuantEventAssessor.QuantEventAssessment assessment = quantEvent.getAlgorithmicAssessment();
        if (assessment != null)
            assessmentString = QuantEventAssessor.flagReasonCodes[assessment.getStatus()];
        model.setValueAt(assessmentString, numRows, colNum++);
        model.setValueAt(QuantEvent.convertCurationStatusToString(quantEvent.getQuantCurationStatus()),
                numRows, colNum++);

        if (alreadySelected)
        {
            alreadySelectedRows.add(numRows);
            model.setValueAt(true, numRows, 0);
        }
    }

    protected int integerizeRatio(float ratio)
    {
        float ratioBound = 10f;
        float logRatioBounded =
                (float) Math.log(Math.min(ratioBound, Math.max(1.0f / ratioBound, ratio)));
        return (int) (logRatioBounded * 100 / (2 * Math.log(ratioBound))) + 50;
    }

    public void displayEvents(List<QuantEvent> quantEvents)
    {
        displayEvents(quantEvents, null);
    }

    /**
     * Display a list of events.  Indicate selection for the events that have already been selected
     * @param quantEvents
     * @param alreadySelectedEventIndices
     */
    public void displayEvents(List<QuantEvent> quantEvents, List<Integer> alreadySelectedEventIndices)
    {
        clearProperties();
        buildFractionNameNumberMap(quantEvents);
        if (fractionNameNumberMap == null || fractionNameNumberMap.isEmpty() || fractionNameNumberMap.size() == 1)
            hideFractionColumn();
        this.quantEvents = quantEvents;
        shadedTableRows = new ArrayList<Integer>();
        boolean anyEventHasAssessment = false;

        for (int i=0; i<quantEvents.size(); i++)
        {
            QuantEvent quantEvent = quantEvents.get(i);
            if (quantEvent.getAlgorithmicAssessment() != null)
                anyEventHasAssessment = true;
            boolean alreadySelected = (alreadySelectedEventIndices != null &&
                alreadySelectedEventIndices.contains(i));
            addEvent(quantEvent, alreadySelected);
        }
//        if (!anyEventHasAssessment)
//            hideAssessmentColumn();
        if (proteinGenesMap == null)
        {
            hideGeneColumn();
        }
    }

    /**
     * Return all checked rows except the alreadySelectedRows rows
     * @return
     */
    public List<QuantEvent> getSelectedEvents()
    {
        List<QuantEvent> selectedQuantEvents = new ArrayList<QuantEvent>();
        int[] selectedIndices = getSelectedIndices();
        if (selectedIndices != null)
        {
            for (int i : selectedIndices)
                selectedQuantEvents.add(quantEvents.get(i));
        }
        return selectedQuantEvents;
    }

    public class JSliderRenderer implements TableCellRenderer
    {
        protected JSlider slider = new JSlider();

        public JSliderRenderer()
        {
            slider.setMinimum(0);
            slider.setMaximum(100);
            slider.setPaintLabels(false);
            slider.setPaintTicks(false);
            slider.setMajorTickSpacing(25); 
            slider.setPreferredSize(new Dimension(280, 15));
            slider.setPreferredSize(new Dimension(100, 15));
            slider.setToolTipText("Log ratio, bounded at 0.1 and 10");
        }

        public JSliderRenderer(float ratio)
        {
            this();
            slider.setValue(integerizeRatio(ratio));
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column)
        {
            //if value is an integer, set the slider.  Otherwise, it's a String for the header, ignore
            if (Integer.class.isAssignableFrom(value.getClass()))
                slider.setValue((Integer)value);
            return slider;
        }
    }

    public class CurationStatusRenderer extends DefaultTableCellRenderer
    {
        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column)
        {
            Component c =
              super.getTableCellRendererComponent(table, value,
                                                  isSelected, hasFocus,
                                                  row, column);
            int status = QuantEvent.parseCurationStatusString((String) value);
            Color color = c.getForeground();

            switch (status)
            {
                case QuantEvent.CURATION_STATUS_GOOD:
                    color = Color.green;
                    break;
                case QuantEvent.CURATION_STATUS_RATIO_ONEPEAK:
                    color = Color.blue;
                    break;
                case QuantEvent.CURATION_STATUS_BAD:
                    color = Color.red;
                    break;
                default:
                    color = Color.black;
                    break;

            }
            c.setForeground(color);
            return c;

        }        
    }

    public class FlagQuantStatusRenderer extends DefaultTableCellRenderer
    {
        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column)
        {
            Component c = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
            QuantEvent event = quantEvents.get(table.convertRowIndexToModel(row));
            int status = QuantEventAssessor.FLAG_REASON_UNEVALUATED;
            if (event.getAlgorithmicAssessment() != null)
                    status = event.getAlgorithmicAssessment().getStatus();//QuantEventAssessor.parseAssessmentCodeString((String) value);
            Color color;
            switch (status)
            {
                case QuantEventAssessor.FLAG_REASON_OK:
                    color = Color.green;
                    break;
                case QuantEventAssessor.FLAG_REASON_UNEVALUATED:
                    color = Color.black;                    
                    break;
                default:
                    color = Color.red;
                    break;                         
            }
            c.setForeground(color);
            return c;
        }
    }

    /**
     * Make the header of the logratio column display a slider with the given ratio
     * @param ratio
     */
    public void setLogRatioHeaderRatio(float ratio)
    {
        JSliderRenderer renderer = new JSliderRenderer(ratio);
        renderer.slider.setToolTipText("Protein log ratio");
        logRatioSliderColumn.setHeaderRenderer(renderer);
    }

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
            return include(ratio);
        }

        public boolean include(float ratio)
        {
            return (ratio <= maxLowRatioValue || ratio >= minHighRatioValue);
        }
    }

    public void showOnlyExtremeRatios(float maxLowRatioValue, float minHighRatioValue)
    {
        RowFilter<TableModel, Object> rf = new RatioRowFilter(maxLowRatioValue, minHighRatioValue);
        sorter.setRowFilter(rf);
    }


    class SelectAllListener implements ItemListener
    {
        public void itemStateChanged(ItemEvent e)
        {
            Object source = e.getSource();
            if (!(source instanceof AbstractButton)) return;
            boolean checked = e.getStateChange() == ItemEvent.SELECTED;
            for(int x = 0, y = getRowCount(); x < y; x++)
            {
                setValueAt(checked, x, 0);
            }
        }
    }

    public Map<String, Integer> getFractionNameNumberMap()
    {
        return fractionNameNumberMap;
    }

    public void setFractionNameNumberMap(Map<String, Integer> fractionNameNumberMap)
    {
        this.fractionNameNumberMap = fractionNameNumberMap;
    }

    public Map<String, List<String>> getProteinGenesMap()
    {
        return proteinGenesMap;
    }

    public void setProteinGenesMap(Map<String, List<String>> proteinGenesMap)
    {
        this.proteinGenesMap = proteinGenesMap;
    }

    protected class QuantEventChangeListener implements ActionListener
    {
        protected int row;
        public QuantEventChangeListener(int row)
        {
            this.row = row;
        }

        public void actionPerformed(ActionEvent event)
        {
            int currentValue =QuantEvent.parseCurationStatusString(
                    (String) model.getValueAt(row, quantCurationColumnIndex));
            int newValue = quantEvents.get(row).getQuantCurationStatus();
            if (currentValue != newValue)
                model.setValueAt(QuantEvent.convertCurationStatusToString(
                        newValue), row, quantCurationColumnIndex);
            updateUI();
        }
    }

    protected class QuantEventAlgAssessmentChangeListener implements ActionListener
    {
        protected int row;
        public QuantEventAlgAssessmentChangeListener(int row)
        {
            this.row = row;
        }

        public void actionPerformed(ActionEvent event)
        {
            //adding evaluation status listener
            int currentAlgValue = QuantEventAssessor.FLAG_REASON_UNEVALUATED;
            String currentValueString = (String) model.getValueAt(row, quantAlgAssessmentColumnIndex);
            if (currentValueString != null && currentValueString.length() > 0)
                currentAlgValue =QuantEventAssessor.parseAssessmentCodeString(currentValueString);
            int newAlgValue = QuantEventAssessor.FLAG_REASON_UNEVALUATED;
            QuantEventAssessor.QuantEventAssessment assessment = quantEvents.get(row).getAlgorithmicAssessment();
            if (assessment != null)
                newAlgValue = assessment.getStatus();
            if (currentAlgValue != newAlgValue)
                model.setValueAt(QuantEventAssessor.getAssessmentCodeDesc(
                        newAlgValue), row, quantAlgAssessmentColumnIndex);
            updateUI();
        }
    }
}
