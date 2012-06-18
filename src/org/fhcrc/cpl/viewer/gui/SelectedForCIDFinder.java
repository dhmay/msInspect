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

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.gui.FeatureSelectionFrame;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.table.*;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.awt.event.ActionEvent;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;


/**
 */
public class SelectedForCIDFinder
{
    protected static Logger _log = Logger.getLogger(SelectedForCIDFinder.class);

    protected double bestMassBucketSize = 0;
    protected double bestScanBucketSize = 0;
    protected int numMassBuckets = 4;
    protected int numScanBuckets = 4;

    public double getBestMassBucketSize()
    {
        return bestMassBucketSize;
    }

    public void setBestMassBucketSize(double bestMassBucketSize)
    {
        this.bestMassBucketSize = bestMassBucketSize;
    }

    public double getBestScanBucketSize()
    {
        return bestScanBucketSize;
    }

    public void setBestHydrophobicityBucketSize(double bestHydrophobicityBucketSize)
    {
        this.bestScanBucketSize = bestHydrophobicityBucketSize;
    }


    public SelectedForCIDFinder()
    {
    }

    /**
     * do the actual work
     */
    public static void identifyFeaturesSelectedForCID()
    {
        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        FeatureSet runMs2Features = run.getTandemFeatureSet(2);
        runMs2Features.setSourceFile(new File(run.getFile().getPath() + ".ms2.tsv"));
//        runMs2Features.setSourceFile(new File("generated_features_for_cid"));
        ClusteringFeatureSetMatcher cFSM =
                new ClusteringFeatureSetMatcher(.1f, FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE, 5);
        cFSM.setElutionMode(ClusteringFeatureSetMatcher.ELUTION_MODE_SCAN);
        cFSM.setUseMassInsteadOfMz(false);
        cFSM.setElutionBucketIncrement(1);
        cFSM.setMassBucketIncrement(.025);

        List<FeatureSet> featureSets =
                (List<FeatureSet>) ApplicationContext.getProperty(SharedProperties.FEATURE_SETS);
        FeatureSet ms1Features = null;

        boolean noFeatureSetToMine = false;
        if (featureSets == null)
            noFeatureSetToMine = true;
        if (!noFeatureSetToMine)
            for (FeatureSet fs : featureSets)
            {
                if (fs.isDisplayed())
                    ms1Features = fs;
            }
        if (ms1Features == null)
            noFeatureSetToMine = true;

        if (noFeatureSetToMine)
        {
            //make this translatable
            //not really necessary since we make this inaccessible if no features are displayed
            ApplicationContext.infoMessage("You must select and display a set of features in order to find the ones selected for CID");
            return;
        }


        FeatureSetMatcher.FeatureMatchingResult featureMatchingResult =
                cFSM.matchFeatures(ms1Features, runMs2Features);

        Feature[] ms1FeaturesSelectedForCID = new Feature[featureMatchingResult.size()];


        //first set CID scan of all ms1 features to not-found
        for (Feature feature : ms1Features.getFeatures())
        {
            feature.setProperty("cidscan",-1);
        }
        //set CID scan of all matched ms1 features to the MS2 feature's scan number
        int i=0;
        for (Feature ms1Feature : featureMatchingResult.getMasterSetFeatures())
        {
            ms1Feature.setProperty("cidscan",featureMatchingResult.getBestMatch(ms1Feature));
            ms1FeaturesSelectedForCID[i++] = ms1Feature;
        }

        FeatureSet ms1FeatureSetForCID = new FeatureSet(ms1FeaturesSelectedForCID,
                                                        FeatureSelectionFrame.FeatureSelectionDialog.nextColor());
//        ms1FeatureSetForCID.setSourceFile(new File("generated_features_for_cid"));

        //copy the feature set list because otherwise the main panel doesn't update
        ArrayList<FeatureSet> featureSetListCopy = new ArrayList<FeatureSet>(featureSets.size());
        featureSetListCopy.addAll(featureSets);
        featureSetListCopy.add(ms1FeatureSetForCID);
        featureSetListCopy.add(runMs2Features);

        ApplicationContext.setProperty(SharedProperties.FEATURE_SETS, featureSetListCopy);
        ApplicationContext.setProperty(SharedProperties.FEATURE_RANGES, featureSetListCopy);

        ResultsDialog resultsDialog =
                new ResultsDialog(ms1Features, ms1FeatureSetForCID, runMs2Features, featureMatchingResult);
        resultsDialog.setVisible(true);
    }

    //menu action
    public static class SelectedForCIDFinderAction extends AbstractAction
            implements PropertyChangeListener
    {
        public SelectedForCIDFinderAction()
        {
            super();
            ApplicationContext.addPropertyChangeListener(SharedProperties.FEATURE_RANGES, this);

            setEnabled(shouldEnable());
        }

        protected boolean shouldEnable()
        {
            List featureSets = (List) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);

            if (featureSets != null)
                for (int i = featureSets.size() - 1; i >= 0; i--)
                {
                    FeatureSet fs = (FeatureSet) featureSets.get(i);
                    if (fs.isDisplayed())
                        return true;
                }
            return false;
        }

        public void propertyChange(PropertyChangeEvent event)
        {
            setEnabled(shouldEnable());
        }

        public void actionPerformed(ActionEvent event)
        {
            identifyFeaturesSelectedForCID();
        }
    }

    protected static class ResultsDialog extends JDialog
    {
        protected ResultsDialogTableModel resultsTableModel = null;
        protected JTable tblResults;
        protected FeatureSet ms1Features;
        protected FeatureSet ms1FeatureSetForCID;
        protected FeatureSet runMs2Features;
        List<Pair<Feature, Feature>> featurePairs;

        public ResultsDialog(FeatureSet ms1Features, FeatureSet ms1FeatureSetForCID,
                             FeatureSet runMs2Features,
                             FeatureSetMatcher.FeatureMatchingResult featureMatchingResult)
        {
            setTitle(TextProvider.getText("FEATURE_SET_SUMMARY"));
            this.setSize(new Dimension(300,130));
            this.setLayout(new BorderLayout());

            resultsTableModel = new ResultsDialogTableModel();

            tblResults = new JTable(resultsTableModel);
            tblResults.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            tblResults.setModel(resultsTableModel);
            tblResults.setColumnModel(resultsTableModel.columnModel);
            tblResults.setTableHeader(new JTableHeader(resultsTableModel.columnModel));
//            tblResults.setSize(300,100);
            JScrollPane scrollpane = new JScrollPane(tblResults);
//            scrollpane.setSize(300,100);
            add(scrollpane, BorderLayout.CENTER);

            JPanel buttonPanel = new JPanel();
            JButton okButton = new JButton(TextProvider.getText("OK"));
//            JButton saveButton = new JButton(TextProvider.getText("SAVE_SELECTED_MS1"));
            buttonPanel.add(okButton);
//            buttonPanel.add(saveButton);
            ListenerHelper listenerHelper = new ListenerHelper(this);
            listenerHelper.addListener(okButton, "buttonOK_actionPerformed");
//            listenerHelper.addListener(saveButton, "buttonSave_actionPerformed");
            add(buttonPanel, BorderLayout.SOUTH);

            this.ms1Features = ms1Features;
            this.ms1FeatureSetForCID = ms1FeatureSetForCID;
            this.runMs2Features = runMs2Features;
            setLocation((int) this.getParent().getLocation().getX() + 300,
                        (int) this.getParent().getLocation().getY() + 300);
        }

        public void buttonOK_actionPerformed(ActionEvent event)
        {
            this.dispose();
        }
/*
        public void buttonSave_actionPerformed(ActionEvent event)
        {
            WorkbenchFileChooser wfc = new WorkbenchFileChooser();
            wfc.showSaveDialog(this);
            File saveFile = wfc.getSelectedFile();
            try
            {
                saveSelectedForCID(ms1FeatureSetForCID, featurePairs, saveFile);
            }
            catch (Exception e)
            {
                ApplicationContext.infoMessage("ERROR");
            }
        }
*/

        protected class ResultsDialogTableModel extends DefaultTableModel
        {
            protected TableColumnModel columnModel = new DefaultTableColumnModel();
            public TableColumn featureSetTypeColumn = new TableColumn();
            public TableColumn featureSetColorColumn = new TableColumn();
            public TableColumn featureNumberColumn = new TableColumn();

            public ResultsDialogTableModel()
            {
                featureSetTypeColumn.setHeaderValue(TextProvider.getText("FEATURES"));
                featureSetTypeColumn.setModelIndex(0);
//                featureSetTypeColumn.setMaxWidth(150);

                featureSetColorColumn.setHeaderValue(TextProvider.getText("COLOR"));
                featureSetColorColumn.setModelIndex(1);
                featureSetColorColumn.setCellRenderer(new FeatureSelectionFrame.FeatureSelectionDialog.ColorRenderer(false));
                featureSetColorColumn.setMaxWidth(50);

                featureNumberColumn.setHeaderValue(TextProvider.getText("NUMBER"));
                featureNumberColumn.setModelIndex(2);
                featureNumberColumn.setMaxWidth(50);

                columnModel.addColumn(featureSetTypeColumn);
                columnModel.addColumn(featureSetColorColumn);
                columnModel.addColumn(featureNumberColumn);
            }

            protected String getNameForColumn(int col)
            {
                switch(col)
                {
                    case 0:
                        return (String) featureSetTypeColumn.getHeaderValue();
                    case 1:
                        return (String) featureSetColorColumn.getHeaderValue();
                    case 2:
                        return (String) featureNumberColumn.getHeaderValue();
                    default:
                        return null;
                }
            }

            public int getRowCount()
            {
                return 3;
            }

            public int getColumnCount()
            {
                return 3;
            }

            public Class getColumnClass(int c)
            {
                switch(c)
                {
                    case 0:
                        return String.class;
                    case 1:
                        return Color.class;
                    case 2:
                        return Number.class;
                    default:
                        return Object.class;
                }
            }

            public Object getValueAt(int row, int col)
            {
                if (row >= 3)
                    return null;

                FeatureSet fs;
                switch (row)
                {
                    case 0:
                        fs = ms1Features;
                        break;
                    case 1:
                        fs = ms1FeatureSetForCID;
                        break;
                    case 2:
                        fs = runMs2Features;
                        break;
                    default:
                        fs = null;
                        break;
                }
                switch (col)
                {
                    case 0:
                        switch (row)
                        {
                            case 0:
                                return TextProvider.getText("MS1_FEATURES");
                            case 1:
                                return TextProvider.getText("MS1_FEATURES_FOR_CID");
                            case 2:
                                return TextProvider.getText("MS2_CIDS");
                            default:
                                return null;
                        }
                    case 1:
                        return fs.getColor();
                    case 2:
                        return fs.getFeatures().length;
                    default:
                        return null;
                }
            }
        }
    }

}
