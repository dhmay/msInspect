package org.fhcrc.cpl.viewer.quant.gui;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.Localizer;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.xml.stream.XMLStreamException;
import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.List;


/**
 *
 */
public class ProteinSummarySelectorFrame extends JFrame
{
    protected static Logger _log = Logger.getLogger(ProteinSummarySelectorFrame.class);

    protected int width = 700;
    protected int height = 800;

    protected ProteinSummaryTable proteinSummaryTable;

    protected ProtXmlReader.Protein selectedProtein = null;

    protected java.util.List<ProtXmlReader.Protein> proteins;

    protected float minProteinProphet = 0.75f;

    public JPanel contentPanel;
    public JPanel summaryPanel;
    public JPanel mainPanel;

    protected JButton buttonOK  = new JButton("OK");

    protected JButton buttonSelectedProtein  = new JButton("DUMMY");


    //Status message
    public JPanel statusPanel;
    public JLabel messageLabel;

    protected Map<String, List<String>> proteinGeneMap;
    

    public ProteinSummarySelectorFrame()
    {
        super();
    }

    public ProteinSummarySelectorFrame(File protXmlFile)
            throws XMLStreamException, FileNotFoundException
    {
        this();
        displayProteins(protXmlFile);
    }

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
//        gbc.anchor=GridBagConstraints.FIRST_LINE_START;
        gbc.insets = new Insets(0,0,0,0);
        gbc.weighty = 1;
        gbc.weightx = 1;

        ListenerHelper helper = new ListenerHelper(this);

        summaryPanel.setMaximumSize(new Dimension(1200, 100));

        gbc.fill = GridBagConstraints.NONE;
       
        gbc.insets = new Insets(5,5,5,5);
        buttonOK.setEnabled(false);
        helper.addListener(buttonOK, "buttonOK_actionPerformed");
        gbc.gridwidth = GridBagConstraints.RELATIVE;
        summaryPanel.add(buttonOK, gbc);

        JButton buttonCancel = new JButton("Cancel");
        helper.addListener(buttonCancel, "buttonCancel_actionPerformed");
        gbc.gridwidth = GridBagConstraints.REMAINDER;        
        summaryPanel.add(buttonCancel, gbc);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.insets = new Insets(0,0,0,0);

        proteinSummaryTable = new ProteinSummaryTable();
        if (proteinGeneMap != null)
            proteinSummaryTable.proteinGeneMap = proteinGeneMap;
        ListSelectionModel tableSelectionModel = proteinSummaryTable.getSelectionModel();
        tableSelectionModel.addListSelectionListener(new MyListSelectionHandler());

        JScrollPane summaryTableScrollPane = new JScrollPane();
        summaryTableScrollPane.setViewportView(proteinSummaryTable);
        mainPanel.add(summaryTableScrollPane, gbc);

        contentPanel.updateUI();
    }

    public void setProteinGeneMap(Map<String, List<String>> proteinGeneMap)
    {
        this.proteinGeneMap = proteinGeneMap;
        if (proteinSummaryTable != null)
            proteinSummaryTable.proteinGeneMap = proteinGeneMap;
    }

    public void displayProteins(File protXmlFile)
            throws XMLStreamException, FileNotFoundException
    {
        initGUI();        
        proteins = new ArrayList<ProtXmlReader.Protein>();

        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);
        ProtXmlReader.ProteinGroupIterator groupIterator = protXmlReader.iterator();

        while (groupIterator.hasNext())
        {
            ProteinGroup proteinGroup = groupIterator.next();
            for (ProtXmlReader.Protein protein : proteinGroup.getProteins())
            {
                if (protein.getProbability() > minProteinProphet && protein.getQuantitationRatio() != null)
                    proteins.add(protein);
            }
        }

        Collections.sort(proteins, new ProteinRatioAscComparator());
        for (ProtXmlReader.Protein protein : proteins)
        {
            proteinSummaryTable.addProtein(protein);
        }
        proteinSummaryTable.updateUI();
    }

    public void addSelectionListener(ActionListener listener)
    {
        buttonSelectedProtein.addActionListener(listener);
    }

    public void buttonOK_actionPerformed(ActionEvent event)
    {
        if (proteinSummaryTable.getSelectedIndex() == -1)
            return;
        selectedProtein = proteins.get(proteinSummaryTable.getSelectedIndex());

        ActionListener[] buttonListeners = buttonSelectedProtein.getActionListeners();

        if (buttonListeners != null)
        {
            for (ActionListener listener : buttonListeners)
                listener.actionPerformed(event);
        }
    }

    public void buttonCancel_actionPerformed(ActionEvent event)
    {
        selectedProtein = null;
        setVisible(false);
    }


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

    public void displayProteins(java.util.List<ProtXmlReader.Protein> proteins)
    {
        this.proteins = proteins;
        displayProteins();
    }

    public void displayProteins()
    {
        proteinSummaryTable.clearRows();
        for (ProtXmlReader.Protein protein : proteins)
            proteinSummaryTable.addProtein(protein);
        proteinSummaryTable.updateUI();
    }

    public static final class ProteinSummaryTable extends JTable
    {
        protected Map<String, List<String>> proteinGeneMap;

        DefaultTableModel model = new DefaultTableModel(0, 5)
            {
                //all cells uneditable
                public boolean isCellEditable(int row, int column)
                {
                    return false;
                }

                public Class getColumnClass(int columnIndex)
                {
                    return String.class;
                }
            };

        public ProteinSummaryTable()
        {
            setModel(model);
            getColumnModel().getColumn(0).setHeaderValue("Protein");
            getColumnModel().getColumn(1).setHeaderValue("Genes");
            getColumnModel().getColumn(2).setHeaderValue("Probability");
            getColumnModel().getColumn(3).setHeaderValue("Ratio");
            getColumnModel().getColumn(4).setHeaderValue("UniquePeptides");

            getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

            TableRowSorter<TableModel> sorter
                    = new TableRowSorter<TableModel>(model);
            setRowSorter(sorter);
        }

        /**
         * Returns model, not view, index
         * @return
         */        
        public int getSelectedIndex()
        {
            ListSelectionModel lsm = this.getSelectionModel();
            if (lsm.isSelectionEmpty())
                return -1;
            // Find out which indexes are selected.
            int minIndex = lsm.getMinSelectionIndex();
            int maxIndex = lsm.getMaxSelectionIndex();
            if (minIndex == maxIndex)
                return convertRowIndexToModel(minIndex);
            else
                return -1;
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

        public void addProtein(ProtXmlReader.Protein protein)
        {
            int numRows = model.getRowCount();
            model.setRowCount(numRows + 1);
            model.setValueAt(protein.getProteinName(), numRows, 0);
            if (proteinGeneMap != null && proteinGeneMap.containsKey(protein.getProteinName()))
            {
                List<String> genes = proteinGeneMap.get(protein.getProteinName());
                if (genes != null && !genes.isEmpty())
                {
                    StringBuffer genesStringBuf = new StringBuffer(genes.get(0));
                    for (int i=1; i<genes.size(); i++)
                        genesStringBuf.append("," + genes.get(i));
                    model.setValueAt(genesStringBuf.toString(), numRows, 1);
                }
            }
            model.setValueAt("" + protein.getProbability(), numRows, 2);            
            if (protein.getQuantitationRatio() != null)
                model.setValueAt("" + protein.getQuantitationRatio().getRatioMean(), numRows, 3);
            model.setValueAt("" + protein.getUniquePeptidesCount(), numRows, 4);

        }
    }

    public ProtXmlReader.Protein getSelectedProtein()
    {
        return selectedProtein;  
    }

    public void setSelectedProtein(ProtXmlReader.Protein selectedProtein)
    {
        this.selectedProtein = selectedProtein;
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
     * display the properties for the selected event, if only one's selected
     */
    public class MyListSelectionHandler implements ListSelectionListener
    {
        public void valueChanged(ListSelectionEvent e)
        {
            if (!e.getValueIsAdjusting())
            {
                ListSelectionModel lsm = (ListSelectionModel) e.getSource();
                if (!lsm.isSelectionEmpty())
                {
                    buttonOK.setEnabled(true);
                }
            }
        }
    }
}
