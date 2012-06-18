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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.FastaLoader;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.MutableTreeNode;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * User: migra
 * Date: Jun 29, 2004
 * Time: 11:48:16 AM
 *
 */
public class OpenFastaDialog extends JDialog implements PropertyChangeListener
{
    private static final double _hMass = PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[PeptideGenerator.H_ION_INDEX];
    private static OpenFastaDialog _instance ;
    static final NumberFormat massFormat = new DecimalFormat("#####.0000");

    private static class ResidueMod
    {
        final double addMass;
        final char residue;
        final String desc;
        final boolean isDefault;
        public ResidueMod(char residue, double addMass, String desc, boolean isDefault)
        {
            this.residue = residue;
            this.addMass = addMass;
            this.desc = desc;
            this.isDefault = isDefault;
        }
        public String toString()
        {
            return desc + " (" + residue + (addMass >= 0.0 ? " +" : " ") + addMass + ")";
        }
    }

    private static final ResidueMod[] residueMods = new ResidueMod[]{
        new ResidueMod('C', 0.0, "Reduced form", false),
        new ResidueMod('C', 58.004, "Iodoacetic Acid", false),
        new ResidueMod('C', 57.0214, "Iodoacetamide", true),     //IAM is our default
        new ResidueMod('K', 6.020124, "SILAC", false),
    };

    private Container contentPanel;
    private JTextField textMass;
    private JTree treePeptide;
    private JLabel labelStatus;
    private File _fastaFile;
    private JTextField textTolerance;
    private JComboBox comboUnits;
    private JComboBox comboMissedCleavages;
    private JComboBox comboCysteine;
    private JButton buttonSearch;
    private JTable tablePeptides;
    private JLabel labelSequence;
    private JButton buttonFindProtein;
    private JTextField textProteinName;
    private JSplitPane splitPane;

    private Protein _protein;

    public static OpenFastaDialog getInstance()
    {
        if (null == _instance)
            _instance = new OpenFastaDialog();

        return _instance;
    }

//TODO: Pre-expand proteins since that's the only interesting thing
    public OpenFastaDialog()
    {
        try
        {
            contentPanel =
                    Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/OpenFastaDialog.xml",
                            this);

        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
            throw new RuntimeException(x);
        }

        setContentPane(contentPanel);

        MutableTreeNode rootNode = new DefaultMutableTreeNode("No Peptides");
        treePeptide.setModel(new DefaultTreeModel(rootNode));
        textTolerance.setText(".1");
        comboUnits.addItem("Daltons");
        comboUnits.addItem("PPM");

        for (int i = 0; i < residueMods.length; i++)
        {
            comboCysteine.addItem(residueMods[i].toString());
            if ( residueMods[i].isDefault )
                comboCysteine.setSelectedIndex(i);
        }

        comboMissedCleavages.addItem("0");
        comboMissedCleavages.addItem("1");
        comboMissedCleavages.addItem("2");
        comboMissedCleavages.setSelectedIndex(0); //1 missed cleavage

        tablePeptides.setModel(new PeptideTableModel(null));

        ActionListener searchListener = new ActionListener(){
            public void actionPerformed(ActionEvent event)
            {
                doSearch();
            }
        };
        textMass.addActionListener(searchListener);
        buttonSearch.addActionListener(searchListener);
        buttonSearch.setDefaultCapable(true);

        tablePeptides.getSelectionModel().addListSelectionListener(new PeptideTableSelectionListener());
        ApplicationContext.addPropertyChangeListener(this);
        this.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        FocusListener topPaneFocusListener = new FocusListener()
        {
            public void focusGained(FocusEvent event)
            {
                getRootPane().setDefaultButton(buttonSearch);
                if (event.getSource() instanceof JTextField)
                {
                    JTextField textField = (JTextField) event.getSource();
                    textField.setSelectionStart(0);
                    textField.setSelectionEnd(textField.getText().length());
                }
            }

            public void focusLost(FocusEvent event)
            {
            }
        };
        comboMissedCleavages.addFocusListener(topPaneFocusListener);
        textMass.addFocusListener(topPaneFocusListener);
        textTolerance.addFocusListener(topPaneFocusListener);
        comboUnits.addFocusListener(topPaneFocusListener);


        textProteinName.addFocusListener(new FocusListener(){
            public void focusGained(FocusEvent event)
            {
                getRootPane().setDefaultButton(buttonFindProtein);
                textProteinName.setSelectionStart(0);
                textProteinName.setSelectionEnd(textProteinName.getText().length());
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void focusLost(FocusEvent event)
            {
            }
        });
        buttonFindProtein.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event)
            {
                _protein = findProtein(textProteinName.getText());
                showProtein(_protein, null);
            }
        });

        ApplicationContext.addPropertyChangeListener(SharedProperties.FEATURE_RANGES, new PropertyChangeListener()
            {
            public void propertyChange(PropertyChangeEvent event)
                {
                    if (null != _protein)
                    {
                        ListSelectionModel lsm = tablePeptides.getSelectionModel();
                        Peptide pep = null;
                        PeptideTableModel tableModel = (PeptideTableModel) tablePeptides.getModel();
                        if (!lsm.isSelectionEmpty())
                            pep = tableModel.getPeptides()[lsm.getMinSelectionIndex()];

                        showProtein(_protein, pep);
                    }
                }
            });
        splitPane.setDividerLocation(200);

        this.setSize(530, 600);
    }


    public File getFastaFile()
    {
        return _fastaFile;
    }

    public void setFastaFile(File fastaFile)
    {
        _fastaFile = fastaFile;
        _protein = null;
        showProtein(null, null);
        treePeptide.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Peptides")));
        this.setTitle(fastaFile.getAbsolutePath());
        ApplicationContext.setProperty("fastaFile", fastaFile);
    }


    //UNDONE: Cache this and listen for change event
    public double getTolerance()
    {
        String toleranceString = textTolerance.getText();
        double tolerance;
        try
        {
            tolerance = Double.parseDouble(toleranceString);
            if (comboUnits.getSelectedIndex() == 1)
                tolerance = (tolerance * getMass()) / 1000000.0; //PPM
        }
        catch (NumberFormatException x)
        {
            return 0;
        }
        return tolerance;
    }

    public double getMass()
    {
        String massString = textMass.getText();
        double mass;
        try
        {
            mass = Double.parseDouble(massString);
        }
        catch (NumberFormatException x)
        {
            return 0;
        }
        return mass;
    }

    private void doSearch()
    {
        double mass = getMass();
        double tolerance = getTolerance();
        PeptideGenerator pepGen = new PeptideGenerator();
        pepGen.setInputFileName(_fastaFile.getAbsolutePath());

        pepGen.setMinMass(mass - tolerance);
        pepGen.setMaxMass(mass + tolerance);
        pepGen.addListener(new PeptideCollector());
        double[] massTab = PeptideGenerator.getMasses(true);
        //Residue modifications
        int i = comboCysteine.getSelectedIndex();
        massTab[residueMods[i].residue] += residueMods[i].addMass;
        pepGen.setMassTable(massTab);
        pepGen.setMaxMissedCleavages(comboMissedCleavages.getSelectedIndex());
        labelStatus.setText("Searching for masses " + String.valueOf(pepGen.getMinMass()) + "-" + String.valueOf(pepGen.getMaxMass()));
        treePeptide.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Searching...")));
        tablePeptides.setModel(new PeptideTableModel(null));
        ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES, null);
        Thread t = new Thread(pepGen);
        t.start();
    }

    //UNDONE: Do this in the background?
    public Protein findProtein(String proteinText)
    {
        proteinText = proteinText.trim().toLowerCase();
        FastaLoader fastaLoader = new FastaLoader(_fastaFile);
        FastaLoader.ProteinIterator iterator = (FastaLoader.ProteinIterator) fastaLoader.iterator();
        List proteins = new ArrayList();
        while (iterator.hasNext())
        {
            Protein protein = (Protein) iterator.next();
            if (protein.getHeader().toLowerCase().indexOf(proteinText) != -1)
            {
                proteins.add(protein);
            }
        }

        if (proteins.size() == 0)
            return null;

        if (proteins.size() == 1)
            return (Protein) proteins.get(0);

        Object[] arr = proteins.toArray();
        return (Protein) JOptionPane.showInputDialog(this, "Multiple proteins found, please pick one", "Select Protein", JOptionPane.PLAIN_MESSAGE, null, arr, arr[0]);
    }


    private void showProtein(Protein protein, Peptide peptide)
    {
        if (null == protein)
        {
            tablePeptides.setModel(new PeptideTableModel(null));
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES, null);
            labelSequence.setText("");
            textProteinName.setText("<no protein>");
            return;
        }

        PeptideGenerator pepGen = new PeptideGenerator();
        double[] massTab = getMassTab();
        pepGen.setMassTable(massTab);
        pepGen.setMaxMissedCleavages(comboMissedCleavages.getSelectedIndex());

        Peptide[] peptides = pepGen.digestProtein(protein);
        PeptideTableModel tableModel =new PeptideTableModel(peptides);
        labelSequence.setText(getProteinHtml(protein, tableModel, peptide));
        tablePeptides.setModel(tableModel);
        //UNDONE: This overwrites search string. Should put protein name elsewhere
        textProteinName.setText(protein.getHeader());
        int index = -1;
        if (null != peptide)
            index = tableModel.findPeptide(peptide);
        if (index != -1)
            tablePeptides.changeSelection(index, 0, true, false);
        //TODO: Put this on another thread
        ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES, tableModel.findHighlightedFeatures());

    }

    public double[] getMassTab()
    {
        double[] massTab = PeptideGenerator.getMasses(true);
        //Residue modifications
        int i = comboCysteine.getSelectedIndex();
        massTab[residueMods[i].residue] += residueMods[i].addMass;
        return massTab;
    }

    public int getMissedCleavages()
    {
        return comboMissedCleavages.getSelectedIndex();
    }
    
    private static final int PLAIN_RESIDUES = 0;
    private static final int PEPTIDE_RESIDUES = 1; //Residues are covered by a peptide
    private static final int FEATURE_PEPTIDE_RESIDUES = 2; //Match a feature within given tolerance
    private static final int SELECTED_PEPTIDE_RESIDUES = 3; //UI selection is over this peptide

    public String getProteinHtml(Protein protein, PeptideTableModel tableModel, Peptide peptide)
    {
        Peptide[] peptides = tableModel.getPeptides();
        String seq = protein.getSequenceAsString();
        StringBuffer sb = new StringBuffer();
        int chIndex = 0;
        sb.append("<html><pre>");
        int pepIndex = 0;
        int residueType;
        while(pepIndex < peptides.length)
        {
            int pepStart = peptides[pepIndex].getStart();
            if (chIndex < pepStart)
            {
                int len  = Math.min(10 - (chIndex % 10), pepStart - chIndex);
                sb.append(seq.substring(chIndex, chIndex + len).toLowerCase());
                chIndex += len;
            }
            else
            {
                residueType = getResidueType(tableModel, pepIndex, peptide);

                int len = Math.min(10 - (chIndex % 10), pepStart + peptides[pepIndex].getLength() - chIndex);
                //If we find an equal or "more selected" peptide in the next spot, render it in preference
                //to this one.
                //TODO: This really should loop for the case where there are 3 peptides that start
                //at same place. First is FEATURE, second is nothing, third is SELECTED
                //Currently first will be rendered FEATURE and rest will be rendered SELECTED
                if (residueType != SELECTED_PEPTIDE_RESIDUES && pepIndex < peptides.length - 1)
                {
                    int resTypeNext = getResidueType(tableModel, pepIndex + 1, peptide);
                    if (resTypeNext >= residueType)
                        len = Math.min(len, peptides[pepIndex+1].getStart() - chIndex);
                }
                if (len <= 0)
                {
                    pepIndex++;
                    continue;
                }

                if (residueType == SELECTED_PEPTIDE_RESIDUES)
                    sb.append("<font color=\"#FF00FF\">");
                else if (residueType == FEATURE_PEPTIDE_RESIDUES)
                    sb.append("<font color=\"#FFA500\">");
                sb.append(seq.substring(chIndex, chIndex + len));
                if (residueType != PEPTIDE_RESIDUES)
                    sb.append("</font>");
                chIndex += len;
            }

            if (chIndex % 60 == 0)
                sb.append("\n");
            else if (chIndex % 10 == 0)
                sb.append(" ");

            while (pepIndex < peptides.length && chIndex >= peptides[pepIndex].getStart() + peptides[pepIndex].getLength())
                pepIndex++;

        }

        while (chIndex < seq.length())
        {
            int len  = Math.min(10 - (chIndex % 10), seq.length() - chIndex);
            sb.append(seq.substring(chIndex, chIndex + len).toLowerCase());
            chIndex += len;
        }

        sb.append("</pre></html>");
        return sb.toString();
    }

    private int getResidueType(PeptideTableModel tableModel, int row, Peptide selected)
    {
        Peptide[] peptides = tableModel.getPeptides();
        Peptide pep = peptides[row];
        if (pep.equals(selected))
            return SELECTED_PEPTIDE_RESIDUES;

        Feature nearestFeature = tableModel.getNearestFeature(row);
        if (null != nearestFeature && Math.abs(nearestFeature.mass - peptides[row].getMass()) <  getTolerance())
            return FEATURE_PEPTIDE_RESIDUES;

        return PEPTIDE_RESIDUES;
    }

    public Feature[] findHighlightedFeatures()
    {
        PeptideTableModel tableModel = (PeptideTableModel) tablePeptides.getModel();
        return tableModel.findHighlightedFeatures();
    }

    public void propertyChange(PropertyChangeEvent event)
    {
        if (!this.isVisible())
            return;

        if (!SharedProperties.SELECTED_POINT.equals(event.getPropertyName()))
	        return;

		if (null != event.getNewValue() && event.getNewValue() instanceof Feature)
		{
			Feature f = (Feature) event.getNewValue();
			if (f.mass > 0)
			{
				textMass.setText("" + f.mass);
				doSearch();
			}
		}
    }

    public class PeptideTableSelectionListener implements ListSelectionListener
    {
        public void valueChanged(ListSelectionEvent e)
        {
            if (e.getValueIsAdjusting())
                return;

            ListSelectionModel lsm = (ListSelectionModel) e.getSource();
            if (!lsm.isSelectionEmpty())
            {
                PeptideTableModel tableModel = (PeptideTableModel) tablePeptides.getModel();
                Peptide pep = tableModel.getPeptides()[lsm.getMinSelectionIndex()];
                labelSequence.setText(getProteinHtml(pep.getProtein(), tableModel, pep));
            }
        }
    }


    private class PeptideTableModel extends AbstractTableModel
    {
        private final String[] colNamesBase = new String[] {"Peptide", "M (monoisotopic)", "[M + H]+", "[M + H]2+"};
        private final String[] colNamesFeatures = new String[] {"Peptide", "M (monoisotopic)", "[M + H]+", "[M + H]2+", "Nearest Feature", "Mass", "Scan"};
        private Peptide[] _peptides;
        private FeatureSet.FeatureSelector _featureSelector = null;
        private double _tolerance;
        private String[] _colNames = colNamesBase;
        private Feature[] _nearestNeighbors = null;
        private Feature[] _highlightedFeatures;


        public PeptideTableModel(Peptide[] peptides)
        {
            _peptides = peptides;
            _featureSelector = (FeatureSet.FeatureSelector) ApplicationContext.getProperty("featureSelector");
            if (null != _featureSelector)
            {
                _colNames = colNamesFeatures;
                if (null != peptides)
                    _nearestNeighbors = new Feature[peptides.length];
            }
            _tolerance = getTolerance();
        }

        public int getRowCount()
        {
            return _peptides == null ? 0 : _peptides.length;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getColumnCount()
        {

            return _colNames.length;
        }

        public Peptide[] getPeptides()
        {
            return _peptides;
        }

        public int findPeptide(Peptide peptide)
        {
            for (int i = 0; i < _peptides.length; i++)
                if (_peptides[i].equals(peptide))
                    return i;

            return -1;
        }

        public Object getValueAt(int row, int col)
        {
            if (null == _peptides || row >= _peptides.length || col >= _colNames.length)
                return null;

            Feature nearestFeature = null;

            switch(col)
            {
                case 0:
                    return _peptides[row].toString();
                case 1:
                    return massFormat.format(_peptides[row].getMass());
                case 2:
                    return massFormat.format(_peptides[row].getMass() + _hMass);
                case 3:
                    return massFormat.format(_peptides[row].getMass()/2 + _hMass);
                case 4:
                    nearestFeature = getNearestFeature(row);
                    return null == nearestFeature ? null : massFormat.format(nearestFeature.mass - _peptides[row].getMass());
                case 5:
                    nearestFeature = getNearestFeature(row);
                    return null == nearestFeature ? null : massFormat.format(nearestFeature.mass);
                case 6:
                    nearestFeature = getNearestFeature(row);
                    return null == nearestFeature ? null : Integer.toString(nearestFeature.scan);
                default:
                    return null;
            }
        }

        public Feature getNearestFeature(int row)
        {
            if (null == _nearestNeighbors)
                return null;

            Feature nearestFeature = _nearestNeighbors[row];

            if (null != nearestFeature)
                return nearestFeature;

            double distance = Double.MAX_VALUE;
            java.util.List featureSets = (java.util.List) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
            if (null == featureSets)
                return null;

            Iterator featureSetIterator = featureSets.iterator();
            while (featureSetIterator.hasNext())
            {
                FeatureSet fs = (FeatureSet) featureSetIterator.next();
                if (!fs.isDisplayed())
                    continue;

                Feature[] features = fs.getFeatures();
                //Just find nearest by mass
                for (int i = 0; i < features.length; i++)
                {
                    Feature feature = features[i];
                    if (Math.abs(feature.mass - _peptides[row].getMass()) < distance)
                    {
                        nearestFeature = feature;
                        distance = Math.abs(feature.mass -  _peptides[row].getMass());
                    }
                }
            }

            _nearestNeighbors[row] = nearestFeature;

            return nearestFeature;
        }

        //TODO: Make this use a property listener approach
        //TODO: Do this on a background thread
        /*
        public Spectrum.Feature[] findHighlightedFeatures(boolean recalc)
        {
            if (recalc)
            {
                //Re-find everything...
                Arrays.fill(_nearestNeighbors, null);
                _highlightedFeatures = null;
            }
            else if (null != _highlightedFeatures)
                return _highlightedFeatures;

            if (null == _nearestNeighbors || null == _peptides)
                return null;

            ArrayList features = new ArrayList();
            double tolerance = getTolerance();

            for (int i = 0; i < getRowCount(); i++)
            {
                Spectrum.Feature feature = getNearestFeature(i);
                if (null != feature && Math.abs(feature.mass - _peptides[i].getMass()) <= tolerance)
                    features.add(feature);
            }

            _highlightedFeatures =(Spectrum.FeatureRange[]) features.toArray(new Spectrum.FeatureRange[features.size()]);
            return _highlightedFeatures;
        }
        */

        public String getColumnName(int col)
        {
            return _colNames[col];
        }

        public Feature[] findHighlightedFeatures()
        {
            if (null == _peptides)
                return null;
            List featureSets = (List) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
            if (null == featureSets)
                return null;

            int size = 0;
            for (int i = 0; i < featureSets.size(); i++)
            {
                FeatureSet fs = (FeatureSet) featureSets.get(i);
                if (!fs.isDisplayed())
                    continue;
                FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
                size += fs.getFeatures().length;
            }

            Feature[] features = new Feature[size];
            int offset = 0;
            for (int i = 0; i < featureSets.size(); i++)
            {
                FeatureSet fs = (FeatureSet) featureSets.get(i);
                if (!fs.isDisplayed())
                    continue;
                System.arraycopy(fs.getFeatures(), 0, features, offset, fs.getFeatures().length);
                offset += fs.getFeatures().length;
            }

            Comparator comp = new Feature.MassAscComparator();
            Arrays.sort(features, comp);
            ArrayList featureList = new ArrayList();
            double tolerance = getTolerance();
            for (int i = 0; i < _peptides.length; i++)
            {
                Peptide pep = _peptides[i];
                float minMass = (float) (pep.getMass() - tolerance);
                float maxMass = (float) (pep.getMass() + tolerance);
                Feature feat = new Feature();
                feat.mass = minMass;
                int indexStart = Arrays.binarySearch(features, feat, comp);
                if (indexStart < 0)
                    indexStart = -indexStart - 1;
                feat.mass = maxMass;
                int indexEnd = Arrays.binarySearch(features, feat, comp);
                if (indexEnd < 0)
                    indexEnd = -indexEnd - 1;
                if (indexEnd >= features.length)
                    indexEnd = features.length -1;

                for (int j = indexStart; j <= indexEnd; j++)
                {
                    //Still possible that nearest feature is not near enough
                    if (features[j].mass >= minMass && features[j].mass <= maxMass)
                    {
                        featureList.add(features[j]);
                    }
                }
            }


            _highlightedFeatures = (Feature[]) featureList.toArray(new Feature[featureList.size()]);

            return _highlightedFeatures;
        }

    }


    private class PeptideCollector implements PeptideGenerator.PeptideListener
    {
        HashMap hm = new HashMap();
        public void handlePeptide(Peptide peptide)
        {
            ArrayList protList = (ArrayList) hm.get(peptide);
            if (null == protList)
            {
                protList = new ArrayList();
                protList.add(peptide.getProtein());
                hm.put(peptide, protList);
            }
            else
                protList.add(peptide.getProtein());
        }

        public void handleDone()
        {
            final HashMap hm = this.hm;
            final double m = getMass();
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    MutableTreeNode rootNode = new DefaultMutableTreeNode("Peptides");
                    Collection peptides = hm.keySet();
                    Iterator iter = peptides.iterator();
                    while (iter.hasNext())
                    {
                        Peptide pep = (Peptide) iter.next();
                        ArrayList proteinList = (ArrayList) hm.get(pep);
                        MutableTreeNode node = new DefaultMutableTreeNode(new PepDisplay(pep, proteinList.size()));
                        for(int i = 0; i < proteinList.size(); i++)
                        {
                            Protein protein = (Protein) proteinList.get(i);
                            MutableTreeNode proteinNode = new DefaultMutableTreeNode(protein);
                            node.insert(proteinNode, i);
                        }
                        rootNode.insert(node, rootNode.getChildCount());
                    }
                    treePeptide.setModel(new DefaultTreeModel(rootNode));
                    treePeptide.addTreeSelectionListener(new TreeSelectionListener() {
                        public void valueChanged(TreeSelectionEvent e)
                        {
                            DefaultMutableTreeNode node = (DefaultMutableTreeNode) treePeptide.getLastSelectedPathComponent();
                            if (null == node)
                                return;
                            Object userObject = node.getUserObject();
                            if (userObject instanceof Protein)
                            {
                                DefaultMutableTreeNode pepNode = (DefaultMutableTreeNode) node.getParent();
                                Peptide pep = ((PepDisplay) pepNode.getUserObject()).peptide;
                                _protein= (Protein) userObject;
                                showProtein(_protein, pep);
                            }
                            else
                                labelSequence.setText("");
                        }
                    });
                    labelStatus.setText(String.valueOf(hm.size()) + " peptides found.");
                }
            });
        }
    }

    public class PepDisplay
    {
        protected Peptide peptide;
        private int nProteins;

        public PepDisplay(Peptide peptide, int nProteins)
        {
            this.nProteins = nProteins;
            this.peptide = peptide;
        }

        public String toString()
        {

           return peptide.toString() + " m=" + massFormat.format((double) peptide.getMass())
                   + " (" + massFormat.format(peptide.getMass() - getMass()) +  ") - " + nProteins + " proteins";
        }
    }


    public static class OpenFastaAction extends AbstractAction
    {
        JFileChooser chooser;


        public OpenFastaAction(JFileChooser chooser)
        {
            super(TextProvider.getText("OPEN_FASTA_DOTDOTDOT"));
            if (null == chooser)
                chooser = new WorkbenchFileChooser();
            this.chooser = chooser;
        }


        public void actionPerformed(ActionEvent evt)
        {
            JFrame frame = ApplicationContext.getFrame();
            int chooserStatus = chooser.showOpenDialog(frame);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            final File file = chooser.getSelectedFile();
            if (null == file)
                return;

            VerifyFile:
            {
                if (!file.exists())
                    break VerifyFile;

                OpenFastaDialog openFastaDialog =
                        new OpenFastaDialog();
                openFastaDialog.setFastaFile(file);
                openFastaDialog.setVisible(true);
                return;
            }

            //JOptionPane.showMessageDialog(frame, "Could not open file.", "Open File", JOptionPane.INFORMATION_MESSAGE, null);
            ApplicationContext.errorMessage("Could not open file: " + file.getPath(), null);
        }
    }
}
