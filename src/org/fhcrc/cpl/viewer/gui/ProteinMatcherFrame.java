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

import org.fhcrc.cpl.toolbox.datastructure.Pair;

import java.util.*;
import java.util.List;
import java.io.*;
import java.awt.event.ActionEvent;
import java.awt.*;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeaturePepXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.ClusteringFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.proteomics.ProteomicsRegressionUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.apache.log4j.Logger;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.chart.axis.NumberAxis;

import javax.swing.*;
import javax.swing.table.*;
import javax.swing.event.ListSelectionEvent;

/**
 * The idea is for this class to be as UI-oriented as possible.  Any heavy lifting should be
 * pushed down into ProteinMatcher or something
 */
public class ProteinMatcherFrame extends JDialog
{
    private static Logger _log = Logger.getLogger(ProteinMatcherFrame.class);

    public static final float DEFAULT_MIN_PEPTIDE_PROPHET_SCORE = .95f;

    //menu action
    public static class ProteinMatcherAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            ProteinMatcherFrame dialog = new ProteinMatcherFrame();
            dialog.setModal(false);
            dialog.setAlwaysOnTop(false);
        }
    }

    FeatureSet _ms1Features = null;
    FeatureSet _regressionMS2Features = null;
    boolean loadedRegressionMS2Features = false;
    FeatureSet _matchingMS2Features = null;
    File _matchingMS2FeatureFile = null;
    boolean loadedMatchingMS2Features = false;
    File _fastaFile = null;

    //variables related to hydrophobicity->retention time prediction
    Map<String,Double> _regressionLineHydroToTime = null;
    Map<String,Double> _regressionLineHydroToScan = null;
    Map<String,Double> _regressionLineTimeToHydro = null;
    Map<String,Double> _regressionLineScanToHydro = null;

    protected int _scanOrTimeMode = ProteomicsRegressionUtilities.REGRESSION_MODE_TIME;
    //TODO: actually control this somewhere. For now, always adaptive
    protected int _peptideTimeToleranceMode =
            PeptideMatcher.PEPTIDE_TIME_TOLERANCE_MODE_FIXED;
    protected int _featureDisplayMode = ProteinDisplay.DISPLAY_MATCHED_FEATURES_MODE;

    protected ArrayList<Protein> _matchedProteins = null;
    protected ArrayList<Protein> _displayedProteins = null;
    protected HashMap<Protein, Map<String,Feature>> _proteinMatchedMS1FeatureMap = null;
    protected Map<Protein, Double> _proteinPercentCoverageMap = null;
    protected HashMap<Protein, ArrayList<Feature>> _proteinSequenceMS2FeatureMap = null;
    protected ArrayList<Pair<Feature,Feature>> _ms2ms1MatchedFeaturePairs = null;
    protected HashMap<Protein,ArrayList<Feature>> _proteinMatchedMS1MS2FeatureMap = null;
    protected AmtDatabase _amtDatabase = null;
    protected ArrayList<Feature> _unmatchedMs2Features = null;
    protected ArrayList<Feature> _currentProteinMS1AndMS2MatchedFeatures = null;

    //initial filter is empty, show all proteins
    protected String _proteinFilterString = "";


    //menu
    public JMenuItem showMS2MatchedProteinsMenuItem;
//    public JMenuItem writePepXmlMs1Ms2MatchesMenuItem;
    public JMenuItem writePepXmlMs2RegressionMenuItem;
    public JMenuItem showPredictedHydroTimePlotMenuItem;
    public Action saveFastaAction = new SaveFastaAction();
//    public Action writePepXmlMs1Ms2Action = new WritePepXmlMatchedMs1Ms2Action();
    public Action writePepXmlMs2RegressionAction = new WritePepXmlMs2RegressionAction();
    public Action setParametersAction = new SetParametersAction();
    public Action loadMS2RegressionAction = new LoadMS2RegressionAction();
    public Action loadOffsetFeatureSetAction = new LoadOffsetFeatureSetAction();
    public Action matchProteinsFromFastaAction = new MatchProteinsFromFastaAction();
    public Action findProteinsMatchedInMS2Action = new FindProteinsMatchedInMS2Action();
    public Action showPredictedHydroTimePlotAction = new ShowPredictedHydroTimePlotAction();


    //protein display
    public JLabel labelNumProteins;
    public JLabel proteinLabel;
    public JTable tblProteins;
    public JTextField textProteinPrefix;
    public JButton buttonFilterProteins;
    protected ProteinTableModel _proteinTableModel = null;
    protected Protein _selectedProtein = null;

    //peptide display
    public JComboBox displayMatchedUnmatchedComboBox;
    protected ArrayList<Feature> _displayedFeatures = null;
    protected ArrayList<Feature> _currentProteinMS2MatchedFeatures  = null;

    public JTable tblFeatures;
    protected FeatureTableModel _featureTableModel = null;
    protected ArrayList<Feature> _selectedFeatures = null;

    //Transitory UI elements (referenced by CommandFileRunner)
    public ParametersDialog mParametersDialog = null;


    //matching parameters
    float deltaMass = FeatureSetMatcher.DEFAULT_DELTA_MASS_ABSOLUTE;
    int deltaMassType = FeatureSetMatcher.DEFAULT_DELTA_MASS_TYPE;
    int deltaScan = FeatureSetMatcher.DEFAULT_DELTA_SCAN;
    double deltaTime = FeatureSetMatcher.DEFAULT_DELTA_TIME;
    double deltaHydrophobicity = FeatureSetMatcher.DEFAULT_DELTA_HYDROPHOBICITY;
    float minPeptideProphet = DEFAULT_MIN_PEPTIDE_PROPHET_SCORE;
    int minMatchedFeatures = 1;
    int minPercentFeatureCoverage=0;
    FeatureSetMatcher featureSetMatcher = null;

    public ProteinMatcherFrame()
    {
        super(ApplicationContext.getFrame(), TextProvider.getText("AMT"));
        initialize();
    }


    /**
     * initial setup of UI and class variables
     */
    public void initialize()
    {
        _ms1Features = getMS1Features();
        if (_ms1Features == null)
        {
            ApplicationContext.infoMessage(TextProvider.getText(
                    "AMT_REQUIRES_DISPLAYED_FEATURES"));
            this.setVisible(false);
            this.dispose();
            return;
        }
        this.setVisible(true);

        //graphical stuff
        try
            {
            JMenuBar jmenu = (JMenuBar) Localizer.getSwingEngine(this).render(
                    "org/fhcrc/cpl/viewer/gui/ProteinMatcherMenu.xml");
            for (int i=0 ; i<jmenu.getMenuCount() ; i++)
                jmenu.getMenu(i).getPopupMenu().setLightWeightPopupEnabled(false);
            this.setJMenuBar(jmenu);
            }
        catch (Exception x)
            {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_LOADING_MENUS"), x);
            throw new RuntimeException(x);
            }

        Container contentPanel;
        try
        {
            contentPanel = Localizer.renderSwixml(
                    "org/fhcrc/cpl/viewer/gui/ProteinMatcherFrame.xml", this);
            setContentPane(contentPanel);
            pack();
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
            throw new RuntimeException(x);
        }

        Dimension d = contentPanel.getPreferredSize();
        setBounds(600, 100, (int)d.getWidth(), (int)d.getHeight());

        ListenerHelper helper = new ListenerHelper(this);
        helper.addListener(tblProteins.getSelectionModel(), "tblProteinsModel_valueChanged");
        helper.addListener(buttonFilterProteins,"buttonFilterProteins_actionPerformed");
        //hitting enter in the text field should act the same as hitting the filter button.  Hack, focus, whatever
        helper.addListener(textProteinPrefix,"buttonFilterProteins_actionPerformed");
        helper.addListener(tblFeatures.getSelectionModel(), "tblFeaturesModel_valueChanged");



        _proteinTableModel = new ProteinTableModel();
        tblProteins.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        tblProteins.setAutoCreateColumnsFromModel(false);
        tblProteins.setModel(_proteinTableModel);
        tblProteins.setColumnModel(_proteinTableModel.columnModel);


        if (null != displayMatchedUnmatchedComboBox)
        {
            //TODO: should really use TextProvider here, and use an internal value for determining state
            displayMatchedUnmatchedComboBox.addItem("matched");
            displayMatchedUnmatchedComboBox.addItem("unmatched peptides");
            displayMatchedUnmatchedComboBox.addItem("unmatched ms2");
            helper.addListener(displayMatchedUnmatchedComboBox,
                               "displayMatchedUnmatchedComboBox_actionPerformed");
        }

        _featureTableModel = new FeatureTableModel();
        tblFeatures.setModel(_featureTableModel);

        tblFeatures.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        tblFeatures.setAutoCreateColumnsFromModel(false);

        ProteinTablePopupMenu proteinPopup = new ProteinTablePopupMenu();

        tblProteins.setComponentPopupMenu(proteinPopup);
        proteinLabel.setComponentPopupMenu(proteinPopup);
    }


    /**
     * clean up and shut down this window
     */
    protected void quit()
    {
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES, null);
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES2, null);
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES3, null);

            this.setVisible(false);
            this.dispose();
    }



    /**
     * handler for the "display matched or unmatched features" combo box
     * @param e
     */
    public void displayMatchedUnmatchedComboBox_actionPerformed(ActionEvent e)
    {
        String comboBoxValue = (String) displayMatchedUnmatchedComboBox.getSelectedItem();
        //TODO: the actual strings "matched" and "unmatched" should be internal values, not displayed
        if ("matched".equals(comboBoxValue) &&
            _featureDisplayMode != ProteinDisplay.DISPLAY_MATCHED_FEATURES_MODE)
        {
            _featureDisplayMode = ProteinDisplay.DISPLAY_MATCHED_FEATURES_MODE;
            _displayedFeatures = new ArrayList<Feature>();
            _displayedFeatures.addAll(_proteinMatchedMS1FeatureMap.get(_selectedProtein).values());
            if (_displayedFeatures == null)
            {
                _displayedFeatures = new ArrayList<Feature>();
            }
            _featureTableModel.setDisplayedFeatures(_displayedFeatures);
            updateMainWindowHighlight();
        }
        else if ("unmatched peptides".equals(comboBoxValue) &&
                _featureDisplayMode != ProteinDisplay.DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE)
        {
//TODO: no longer doing this, should remove
            _featureDisplayMode = ProteinDisplay.DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE;
//            _displayedFeatures = _unmatchedProteinPeptideMap.get(_selectedProtein);
            if (_displayedFeatures == null)
            {
                _displayedFeatures = new ArrayList<Feature>();
            }
            _featureTableModel.setDisplayedFeatures(_displayedFeatures);
            updateMainWindowHighlight();
        }
        else if ("unmatched ms2".equals(comboBoxValue) &&
                _featureDisplayMode != ProteinDisplay.DISPLAY_UNMATCHED_MS2_FEATURES_MODE)
        {
            _featureDisplayMode = ProteinDisplay.DISPLAY_UNMATCHED_MS2_FEATURES_MODE;
            _displayedFeatures = _unmatchedMs2Features;
            _featureTableModel.setDisplayedFeatures(_displayedFeatures);
            updateMainWindowHighlight();
        }
    }


    /**
     * Handler for protein selection
     * @param e
     */
    public void tblProteinsModel_valueChanged(ListSelectionEvent e)
    {
        int selectedRow = tblProteins.getSelectedRow();
        if (selectedRow >= 0)
        {
            _selectedProtein =  _proteinTableModel.getSelectedProtein(tblProteins.getSelectedRow());
/*
System.err.println("**Name="+_selectedProtein.getTextCode());
if (_selectedProtein.getAliases() != null)
  for (int i=0; i<_selectedProtein.getAliases().length; i++) System.err.println("---alias="+_selectedProtein.getAliases()[i].getDescription());
System.err.println("**Header: "+ _selectedProtein.getHeader());
System.err.println("**OrigHeader: "+ _selectedProtein.getOrigHeader());
ProteinDisplay.openNCBIBrowserWindow(_selectedProtein);
*/
            _selectedFeatures = null;
            if (_featureDisplayMode == ProteinDisplay.DISPLAY_MATCHED_FEATURES_MODE)
            {
                _displayedFeatures = new ArrayList<Feature>();
                _displayedFeatures.addAll(_proteinMatchedMS1FeatureMap.get(_selectedProtein).values());

                if (_displayedFeatures == null)
                {
                    _displayedFeatures = new ArrayList<Feature>();
                }
                _featureTableModel.setDisplayedFeatures(_displayedFeatures);
            }
//TODO: no longer doing this, should remove
            else if (_featureDisplayMode == ProteinDisplay.DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE)
            {
//                _displayedFeatures = _unmatchedProteinPeptideMap.get(_selectedProtein);
                if (_displayedFeatures == null)
                {
                    _displayedFeatures = new ArrayList<Feature>();
                }
                _featureTableModel.setDisplayedFeatures(_displayedFeatures);
            }
            else if (_featureDisplayMode == ProteinDisplay.DISPLAY_UNMATCHED_MS2_FEATURES_MODE)
            {
                _displayedFeatures = _unmatchedMs2Features;
            }

            _currentProteinMS2MatchedFeatures = getMatchedMS2Features(_selectedProtein);

            //update protein sequence html
            updateProteinSequence();
            updateMainWindowHighlight();
        }
    }

    /**
     * update the main window feature highlighting
     */
    protected void updateMainWindowHighlight()
    {
        //update highlighted features in main window
        Feature[] displayedFeatureArray = null;
        if (_displayedFeatures != null)
            displayedFeatureArray = _displayedFeatures.toArray(new Feature[0]);

        if (_featureDisplayMode == ProteinDisplay.DISPLAY_MATCHED_FEATURES_MODE)
        {
            Feature[] currentProteinMS2MatchedFeatureArray = null;
            if (_currentProteinMS2MatchedFeatures != null)
                currentProteinMS2MatchedFeatureArray = _currentProteinMS2MatchedFeatures.toArray(new Feature[0]);

        ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES, displayedFeatureArray);
        ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES2, currentProteinMS2MatchedFeatureArray);
        ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES3, null);
        }
        else if (_featureDisplayMode == ProteinDisplay.DISPLAY_UNMATCHED_MS2_FEATURES_MODE ||
                 _featureDisplayMode == ProteinDisplay.DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE)
        {
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES,
                                           displayedFeatureArray);
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES2, null);
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES3, null);

        }
    }

    /**
     * Update the protein sequence that's displayed in the middle of the window
     */
    protected void updateProteinSequence()
    {
        _currentProteinMS1AndMS2MatchedFeatures =
                getMatchedMS1AndMS2Features(_selectedProtein);
        //update protein sequence html
        proteinLabel.setText(ProteinDisplay.getProteinSequenceHtml(_selectedProtein,_displayedFeatures,
                                                                   _currentProteinMS2MatchedFeatures,
                                                                   _currentProteinMS1AndMS2MatchedFeatures,
                                                                   _selectedFeatures, _featureDisplayMode));
    }

    /**
     * get the matched ms1 features for a given protein.  If none, return empty
     * arraylist
     * @param protein
     * @return
     */
    protected ArrayList<Feature> getMatchedMS1Features(Protein protein)
    {
        ArrayList<Feature> result = null;
        if (_proteinMatchedMS1FeatureMap != null)
        {
            result.addAll(_proteinMatchedMS1FeatureMap.get(protein).values());
        }
        if (result == null)
            result = new ArrayList<Feature>();
        return result;
    }

    /**
     * get the matched ms2 features for a given protein
     * @param protein
     * @return
     */
    protected ArrayList<Feature> getMatchedMS2Features(Protein protein)
    {
        ArrayList<Feature> result = null;
        if (_proteinSequenceMS2FeatureMap != null)
            result = _proteinSequenceMS2FeatureMap.get(protein);
        if (result == null)
            result = new ArrayList<Feature>();
        return result;
    }

    /**
     * get the matched features that match both ms1 and ms2 for a given protein.  Features are considered
     * identical if they have the same peptide.
     * NOTE: The actual Feature itself will be pulled from the MS2 feature list
     * @param protein
     * @return
     */
    protected ArrayList<Feature> getMatchedMS1AndMS2Features(Protein protein)
    {
/*
        ArrayList<Pair> result = new ArrayList<Pair>();
        ArrayList<Feature> ms1Features = getMatchedMS1Features(protein);
        if (ms1Features == null || ms1Features.size() == 0)
            return result;
        ArrayList<Feature> ms2Features = getMatchedMS2Features(protein);
                if (ms2Features == null || ms2Features.size() == 0)
            return result;


        //take no chances.  Don't use collection methods to compare features, go right to the peptide
        for (int i=0; i<ms1Features.size(); i++)
        {
            Feature ms1Feature = ms1Features.get(i);
            for (int j=0; j<ms2Features.size(); j++)
            {
                if (ms1Feature.getPeptide() != null &&
                    ms1Feature.getPeptide().equalsIgnoreCase(ms2Features.get(j).getPeptide()))
                {
                    result.add(ms2Features.get(j));
                    break;
                }
            }
        }
        return result;
*/
        if (_proteinMatchedMS1MS2FeatureMap == null)
            return new ArrayList<Feature>(0);
        return _proteinMatchedMS1MS2FeatureMap.get(protein);
    }

    /**
     * hide the choice of regression mode, time or scan
     * @return
     */
    protected Map<String,Double> getRegressionLineHydroToScanOrTime()
    {
        if (_scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME)
            return _regressionLineHydroToTime;
        else
            return _regressionLineHydroToScan;
    }

    /**
     * hide the choice of regression mode, time or scan
     * @return
     */
    protected Map<String,Double> getRegressionLineScanOrTimeToHydro()
    {
        if (_scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME)
            return _regressionLineTimeToHydro;
        else
            return _regressionLineScanToHydro;
    }

    /**
     * Query the user for an MS2 feature file to load for regression
     */
    public void loadMS2RegressionFeatures()
    {
        File regressionFeatureFile = WorkbenchFileChooser.chooseExistingFile(
                TextProvider.getText("MS2_FILE_FOR_REGRESSION"));
        loadMS2RegressionFeatures(regressionFeatureFile);
    }

    /**
     * load a specified MS2 feature file for regression.  Calculate the regression line for both scan and time
     * Note:  any feature offsets found in the feature file will NOT be used in regression line
     * calculation, as that line needs to be unmuddled by individual known feature offsets
     * @param regressionFeatureFile
     */
    public void loadMS2RegressionFeatures(File regressionFeatureFile)
    {
        if (regressionFeatureFile == null)
            return;
        try
        {
            _regressionMS2Features = new FeatureSet(regressionFeatureFile);
//Arrays.sort(_regressionMS2Features.getFeatures(),new Feature.MassAscComparator());
//for (int i=0; i<_regressionMS2Features.getFeatures().length; i++)
//  System.err.println(_regressionMS2Features.getFeatures()[i].getMass() + ", " +_regressionMS2Features.getFeatures()[i].getPProphet());

            FeatureSet.FeatureSelector peptideProphetFeatureSelector = new FeatureSet.FeatureSelector();
            peptideProphetFeatureSelector.setMinPProphet(minPeptideProphet);
            MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);

            _regressionMS2Features.populateTimesForMS2Features(run);
            _regressionMS2Features.filter(peptideProphetFeatureSelector);

            _regressionLineHydroToScan = AmtUtilities.calculateHydrophobicityScanOrTimeRelationship(
                    _regressionMS2Features.getFeatures(), ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN,
                     false);
            _regressionLineHydroToTime = AmtUtilities.calculateHydrophobicityScanOrTimeRelationship(
                    _regressionMS2Features.getFeatures(), ProteomicsRegressionUtilities.REGRESSION_MODE_TIME,
                     false);
            _regressionLineScanToHydro = AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                    _regressionMS2Features.getFeatures(), ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN,
                    false);
            _regressionLineTimeToHydro = AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                    _regressionMS2Features.getFeatures(), ProteomicsRegressionUtilities.REGRESSION_MODE_TIME,
                     false);


            if (_regressionLineHydroToScan == null ||
                    _regressionLineHydroToTime == null)
            {
                ApplicationContext.infoMessage(TextProvider.getText(
                        "ERROR_LOADING_FEATURE_FILE_FILENAME",
                        regressionFeatureFile.getAbsolutePath()));
                quit();
            }
            else
            {
                double[] coeffs = new double[2];
                coeffs[1] = AmtUtilities.getSlopeFromRegressionLine(getRegressionLineScanOrTimeToHydro());
                coeffs[0] = AmtUtilities.getInterceptFromRegressionLine(getRegressionLineScanOrTimeToHydro());
                //build a map of known offsets from the features
                AmtDatabaseBuilder.createAmtDatabaseForRun(_regressionMS2Features,
                        _scanOrTimeMode,
                        coeffs,
                        true, null, false
                );
                loadedRegressionMS2Features = true;
                writePepXmlMs2RegressionMenuItem.setEnabled(true);
                showPredictedHydroTimePlotMenuItem.setEnabled(true);
            }


//System.err.println("***Parameters: " + _ms1Features.getFeatures().length + ", " + _regressionMS2Features.getFeatures().length + ", " + deltaMass + ", " + deltaMassType + ", " + deltaTime + ", "  + Window2DFeatureSetMatcher.REGRESSION_MODE_TIME + ", " + minPeptideProphet);
            //enable matching menu item, since we now have a basis for prediction
//System.err.println("****placeholder*** Hydro-scan relationship:  Slope: " +
//                   _regressionLineHydroToScan.first + "; Intercept: " + _regressionLineHydroToScan.second +"\n" +
//                   "\nHydro-time relationship:  Slope: " +
//                   _regressionLineHydroToTime.first + "; Intercept: " + _regressionLineHydroToTime.second +"\n");

        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText(
                    "ERROR_LOADING_FEATURE_FILE_FILENAME",
                    regressionFeatureFile.getAbsolutePath()), e);
            quit();
        }
    }

    public void buttonFilterProteins_actionPerformed(ActionEvent e)
    {
        _proteinFilterString = textProteinPrefix.getText();

        filterDisplayedProteins();
    }

    /**
     * handler for peptide selection
     * @param e
     */
    public void tblFeaturesModel_valueChanged(ListSelectionEvent e)
    {
        int selectedRow = tblFeatures.getSelectedRow();
        if (selectedRow >= 0)
        {
            _selectedFeatures =  _featureTableModel.getSelectedFeatures(tblFeatures.getSelectedRows());
            updateProteinSequence();
            ApplicationContext.setProperty(SharedProperties.HIGHLIGHT_FEATURES3, _selectedFeatures.toArray(new Feature[0]));
            //select the feature -- scroll to it in the main window, and show its details
            if (_selectedFeatures.size() == 1)
                ((WorkbenchFrame) (Application.getInstance().getFrame())).imageComponent.selectFeature(_selectedFeatures.get(0));
        }
    }

    /**
     * Narrow down the displayed proteins, and remove protein and peptide selections
     */
    protected void filterDisplayedProteins()
    {
        if ("".equals(_proteinFilterString))
            _displayedProteins = _matchedProteins;
        else
        {
            String regexString = _proteinFilterString.toUpperCase();
//System.err.println("Before: " +regexString);
            //wildcards
            if (_proteinFilterString.contains("*"))
            {
                StringBuffer regexBuf = new StringBuffer();
                for (int i=0; i<_proteinFilterString.length(); i++)
                {
                    char currentChar = _proteinFilterString.charAt(i);
                    if (currentChar == '*')
                    {
                        regexBuf.append(".*");
                    }
                    else
                        regexBuf.append(Character.toUpperCase(currentChar));
                }
                regexString = regexBuf.toString();
            }
//System.err.println("after: " +regexString);

            _displayedProteins = new ArrayList<Protein>();

            for (Protein protein : _matchedProteins)
            {
                String currentHeader = protein.getHeader();
                if (currentHeader.toUpperCase().matches(regexString))
                {
                    _displayedProteins.add(protein);
                }
            }
            updateDisplayedProteins();

            //clear out dependent fields
            proteinLabel.setText("");
            _displayedFeatures = new ArrayList<Feature>(0);
            _featureTableModel.setDisplayedFeatures(_displayedFeatures);
        }
    }

    /**
     * Pick a featureset to represent the ms1 features.  By convention this is the most
     * recently chosen, visible featureset
     * @return
     */
    protected FeatureSet getMS1Features()
    {
        FeatureSet result = null;

        //TODO:  Does this do what it's supposed to?
        List<FeatureSet> featureSets = (List<FeatureSet>) ApplicationContext.getProperty(
                SharedProperties.FEATURE_SETS);

        if (null == featureSets)
            return null;
        else
            {
                for (FeatureSet fs : featureSets)
                {
                    if (fs.isDisplayed())
                        result = fs;
                }
            }
        return result;
    }

    /**
     * For CommandFileRunner
     * @param fastaFile
     */
    public void findMatchedMS1Proteins(File fastaFile)
    {
        _fastaFile = fastaFile;
        findMatchedMS1Proteins();
    }

    /**
     * Find all proteins in the loaded fasta file that match MS1 features
     */
    public void findMatchedMS1Proteins()
    {
        try
        {
            ArrayList<Protein> proteinsInFasta = ProteinUtilities.loadProteinsFromFasta(_fastaFile);
            //handling for no matches
            if (proteinsInFasta.size() == 0)
            {
                ApplicationContext.infoMessage(TextProvider.getText("NO_PROTEINS_FOUND"));
                return;
            }

            MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);

            double deltaScanOrTime = deltaScan;
            if (_scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME)
                deltaScanOrTime = deltaTime;
            Protein[] proteinArray = proteinsInFasta.toArray(new Protein[0]);
            //TODO: the last argument in this call indicates that we should be matching peptides not
            //TODO: found in the AMT database. this may not be the right thing to do in all cases...
            //TODO: might need to provide control over this in the UI
            _proteinMatchedMS1FeatureMap =
                ProteinMatcher.findMatchesForAllProteins(_ms1Features, proteinArray,
                                                         MS2ExtraInfoDef.getFeatureSetModifications(_regressionMS2Features),
                                                         _amtDatabase,
                                                         minMatchedFeatures);
//            _proteinPercentCoverageMap = ProteinMatcher.createProteinPercentCoverageMap(_proteinMatchedMS1FeatureMap);
            HashMap<Protein, Map<String,Feature>> filteredProteinMatchedMS1FeatureMap =
                    new HashMap<Protein, Map<String,Feature>>();
            for (Protein currentProtein : _proteinPercentCoverageMap.keySet())
            {
                if (_proteinPercentCoverageMap.get(currentProtein) >= minPercentFeatureCoverage)
                {
                    filteredProteinMatchedMS1FeatureMap.put(currentProtein,
                            _proteinMatchedMS1FeatureMap.get(currentProtein));
                }
            }
            _proteinMatchedMS1FeatureMap = filteredProteinMatchedMS1FeatureMap;

            //Create a map of features to the number of peptide masses that match them
//TODO: restore this            
//            createProteinMS1FeatureNumberMassMatchesMap(proteinArray);

//Iterator<ArrayList<Feature>> featureArrayIterator = _proteinMatchedMS1FeatureMap.values().iterator();
//while (featureArrayIterator.hasNext())
//{
//ArrayList<Feature> bunchOFeatures = featureArrayIterator.next();
//for (int i=0; i<bunchOFeatures.size(); i++)
//    System.err.println("Peptide: " + bunchOFeatures.get(i).getPeptide() +
//                        ", offset: " + bunchOFeatures.get(i).getProperty("PREDICTED_TIME_OFFSET"));
//}
            //update UI
            _displayedFeatures = new ArrayList<Feature>(0);
            _featureTableModel.setDisplayedFeatures(_displayedFeatures);

            Set<Protein> matchedProteinSet = _proteinMatchedMS1FeatureMap.keySet();
            Iterator<Protein> matchedProteinIterator = matchedProteinSet.iterator();
            _matchedProteins = new ArrayList<Protein>(matchedProteinSet.size());
            while (matchedProteinIterator.hasNext())
                _matchedProteins.add(matchedProteinIterator.next());

//            Collections.sort(_matchedProteins, new Protein.HeaderComparator());
            //sort by MS1 matches
            Collections.sort(_matchedProteins, new ProteinMS1PercentCoverageComparator());
            _displayedProteins = _matchedProteins;

            updateDisplayedProteins();

            _proteinTableModel.hideMS2Columns();

            tblProteins.getSelectionModel().setSelectionInterval(0, 0);
            tblProteinsModel_valueChanged(null);

            showMS2MatchedProteinsMenuItem.setEnabled(true);
            buttonFilterProteins.setEnabled(true);

        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_MATCHING_FEATURES"), e);

        }

    }


    /**
     * Write out a pepXml file containing the features in the regression MS2 file.
     * Each feature will contain the observed offset from predicted elution time.
     * NOTE: if a secondary offset file is loaded explicitly, those offsets will
     * also get written out along with the offsets from the regression file
     * @param outPepXmlFile
     */
    public void writePepXmlMs2RegressionOffsetFile(File outPepXmlFile)
    {
        if (!loadedRegressionMS2Features)
            return;

        try
        {
            FeaturePepXmlWriter pepXmlWriter =
                    new FeaturePepXmlWriter(_regressionMS2Features.getFeatures(),
                            MS2ExtraInfoDef.getFeatureSetModifications(_regressionMS2Features));
            pepXmlWriter.write(outPepXmlFile);
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
    }


    /**
     * Write out a pepXml file containing only the features that were matched in both MS1
     * and MS2.  Each feature will contain the observed offset from predicted elution time
     * @param outPepXmlFile
     */
/*
    public void writePepXmlMs1Ms2Matches(File outPepXmlFile)
    {
        if (!loadedMatchingMS2Features)
            return;

        ArrayList<Feature> allMatchedFeatures = new ArrayList<Feature>();

        Iterator<Protein> proteinIterator = _proteinMatchedMS1MS2FeatureMap.keySet().iterator();
        while (proteinIterator.hasNext())
        {
            ArrayList<Feature> currentProteinMatchedFeatures = getMatchedMS1AndMS2Features(proteinIterator.next());
            allMatchedFeatures.addAll(currentProteinMatchedFeatures);
        }

        try
        {

            FeaturePepXmlWriter pepXmlWriter = new FeaturePepXmlWriter(allMatchedFeatures.toArray(new Feature[0]),
                                                         _matchingMS2Features.getModifications());
            pepXmlWriter.write(outPepXmlFile);
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
    }
*/

    /**
     * Cover method passing in a file
     * @param matchingMS2FeatureFile
     */
    public void findMatchedMS2Proteins(File matchingMS2FeatureFile)
    {
        try
        {
            _matchingMS2Features = new FeatureSet(matchingMS2FeatureFile);
            if (_matchingMS2Features.getLoadStatus() != FeatureSet.FEATURESET_LOAD_SUCCESS)
            {
                ApplicationContext.infoMessage(TextProvider.getText("" +
                        "ERROR_LOADING_FEATURE_FILE_FILENAME",
                        matchingMS2FeatureFile.getAbsolutePath()) +
                        ": " + _matchingMS2Features.getLoadStatusMessage());
                return;
            }

            FeatureSet.FeatureSelector peptideProphetFeatureSelector = new FeatureSet.FeatureSelector();
            peptideProphetFeatureSelector.setMinPProphet(minPeptideProphet);
            _matchingMS2Features.filter(peptideProphetFeatureSelector);

            _matchingMS2FeatureFile = matchingMS2FeatureFile;

        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText(
                    "ERROR_LOADING_FEATURE_FILE_FILENAME",
                    matchingMS2FeatureFile.getAbsolutePath()), e);

        }
        findMatchedMS2Proteins();
    }

    /**
     * Identify all matchingMS2 features with peptides that occur in the matched proteins.
     * Match MS1 features with matchingMS2 Features
     */
    public void findMatchedMS2Proteins()
    {
        MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
        //using the ms2 feature set for _matching_, find the features with peptides that are
        //found in each protein.  This is used for display of matched peptides.
        //Only match if peptideprophet score above supplied cutoff
//TODO: restore this, or something
//        _proteinSequenceMS2FeatureMap =
//                ProteinMatcher.mapProteinsToFeatures(_matchedProteins,
//                                                     _matchingMS2Features.getFeatures(),
//                                                     minPeptideProphet);

        //this arraylist will hold all the unmatched ms2 features
        _unmatchedMs2Features = new ArrayList<Feature>();

        ClusteringFeatureSetMatcher featureSetMatcher =
                new ClusteringFeatureSetMatcher();
        featureSetMatcher.init(deltaMass, deltaMassType,
                        (float) deltaHydrophobicity);
        FeatureSetMatcher.FeatureMatchingResult featureMatchingResult =
                featureSetMatcher.matchFeatures(_ms1Features,
                        _matchingMS2Features);

        //associate ms1/ms2 matches with individual proteins based on peptide
        //(may be multiply associated)
        //TODO: should probably move this into ProteinMatcher for commandline use
        _proteinMatchedMS1MS2FeatureMap = new HashMap<Protein,ArrayList<Feature>>();
        Iterator<Protein> proteinIterator = _proteinMatchedMS1FeatureMap.keySet().iterator();
        while (proteinIterator.hasNext())
        {
            Protein currentProtein = proteinIterator.next();
            String currentProteinSequence = currentProtein.getSequenceAsString();
            ArrayList<Feature> currentProteinMS1MS2Matches = new ArrayList<Feature>();
            for (int i=0; i<_ms2ms1MatchedFeaturePairs.size(); i++)
            {
                Pair currentPair = _ms2ms1MatchedFeaturePairs.get(i);
                Feature ms2Feature = (Feature) currentPair.first;
                //null test shouldn't really be necessary.  Paranoia
                if (MS2ExtraInfoDef.getFirstPeptide(ms2Feature) != null &&
                    currentProteinSequence.contains(MS2ExtraInfoDef.getFirstPeptide(ms2Feature)))
                {
                    currentProteinMS1MS2Matches.add(ms2Feature);
                }
            }
            _proteinMatchedMS1MS2FeatureMap.put(currentProtein,currentProteinMS1MS2Matches);
        }

        //display changes
        if (_selectedProtein != null)
            _currentProteinMS2MatchedFeatures = getMatchedMS2Features(_selectedProtein);
        tblProteinsModel_valueChanged(null);
        loadedMatchingMS2Features = true;
        _proteinTableModel.showMS2Columns();
        _featureTableModel.fireTableDataChanged();
    }

    protected void updateDisplayedProteins()
    {
        _proteinTableModel.setDisplayedProteins(_displayedProteins);
        labelNumProteins.setText(""+_displayedProteins.size());
    }


    /**
     * Loads a FeatureSet from a file, solely for capturing a mapping between peptides and
     * their observed offsets.  Throws away the featureset
     * @param offsetFeatureFile
     */
    public void loadOffsetsFromFeatureFile(File offsetFeatureFile)
    {
/*
        FeatureSet features = null;
        int numPreviousOffsets = _peptideOffsetMap.size();
        try
        {
            features = new FeatureSet(offsetFeatureFile);
            if (features.getLoadStatus() != FeatureSet.FEATURESET_LOAD_SUCCESS)
            {
                ApplicationContext.infoMessage(TextProvider.getText("ERROR_LOADING_FEATURE_FILE_FILENAME", offsetFeatureFile.getAbsolutePath()) +
                        ": " + features.getLoadStatusMessage());
                return;
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_LOADING_FEATURE_FILE_FILENAME", offsetFeatureFile.getAbsolutePath()), e);
            return;
        }
        //build a map of known offsets from the features
        ProteinMatcher.augmentPeptideOffsetMap(_peptideOffsetMap, features.getFeatures());
        // _regressionMS2FeaturesPeptideOffsetMap = new HashMap<String,Double>();
        ApplicationContext.infoMessage(TextProvider.getText("LOADED_X_OFFSETS_FROM_FILE","" + (_peptideOffsetMap.size() - numPreviousOffsets)));
*/
    }

    protected void showPredictedHydroTimePlot()
    {
        //clone because we're sorting
        Feature[] regressionMS2FeatureArray = _regressionMS2Features.getFeatures().clone();
        Arrays.sort(regressionMS2FeatureArray, new Feature.ScanAscComparator());
        int numQualifyingFeatures = 0;
        for (int i=0; i<regressionMS2FeatureArray.length; i++)
        {
            if (MS2ExtraInfoDef.getPeptideProphet(regressionMS2FeatureArray[i]) >=
                    minPeptideProphet)
                numQualifyingFeatures++;
        }
        float[][] scatterPlotData = new float[2][numQualifyingFeatures];
        int[] sortedScanArray = new int[numQualifyingFeatures];
        int arrayIndex = 0;
        for (int i=0; i< regressionMS2FeatureArray.length; i++)
        {
            if (MS2ExtraInfoDef.getPeptideProphet(regressionMS2FeatureArray[i]) >=
                    minPeptideProphet)
            {
                Feature ms2Feature = regressionMS2FeatureArray[i];
                scatterPlotData[0][arrayIndex] = (float)
                    AmtUtilities.calculateNormalizedHydrophobicity(
                            MS2ExtraInfoDef.getFirstPeptide(ms2Feature));
//                scatterPlotData[1][arrayIndex] = ms2Feature.getScan();
                sortedScanArray[arrayIndex] = ms2Feature.getScan();
                arrayIndex++;
            }
        }
        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        double[] timesArrayTemp = AmtUtilities.getTimesForSortedScanArray(run,
                                                                            sortedScanArray);
        float[] timesArray = new float[timesArrayTemp.length];
        for (int i=0; i<timesArrayTemp.length; i++)
            timesArray[i] = (float) timesArrayTemp[i];
        scatterPlotData[1] = timesArray;

        FastScatterPlot scatterPlot = new FastScatterPlot(scatterPlotData,
              new NumberAxis(TextProvider.getText("CALCULATED_HYDROPHOBICITY")),
              new NumberAxis(TextProvider.getText("TIME")));
        scatterPlot.setPaint(Color.BLUE);
        scatterPlot.setOutlineStroke(new BasicStroke(10));
        ChartDialog chartDialog = new ChartDialog(scatterPlot);
        chartDialog.setLocation(getLocation());
        chartDialog.setVisible(true);
    }










    /**
     * tablemodel for the protein list table
     */
    private class ProteinTableModel extends DefaultTableModel
    {
        protected final String PROTEIN_COLUMN_NAME="PROTEIN";
        protected final String MS2_FEATURES_COLUMN_NAME="MS2_FEATURES";
        protected final String MS1_FEATURES_COLUMN_NAME="MS1_FEATURES";
        protected final String MS1_MS2_FEATURES_COLUMN_NAME="MS1_AND_MS2";
        protected final String MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME="PERCENT_COVERAGE";

        protected ArrayList<Protein> _displayedProteins = new ArrayList<Protein>();

        protected TableColumnModel columnModel = null;

        Hashtable<String,TableColumn> columnHash;
        Hashtable<TableColumn,String> columnReverseHash;

        public ProteinTableModel()
        {
            columnHash = new Hashtable<String,TableColumn>();
            columnReverseHash = new Hashtable<TableColumn,String>();

            TableColumn proteinColumn = new TableColumn();
            proteinColumn.setHeaderValue(TextProvider.getText(PROTEIN_COLUMN_NAME));
            proteinColumn.setModelIndex(0);
            columnHash.put(PROTEIN_COLUMN_NAME, proteinColumn);
            columnReverseHash.put(proteinColumn,PROTEIN_COLUMN_NAME);

            TableColumn ms1Column = new TableColumn();
            ms1Column.setHeaderValue(TextProvider.getText(MS1_FEATURES_COLUMN_NAME));
            ms1Column.setModelIndex(1);
            columnHash.put(MS1_FEATURES_COLUMN_NAME, ms1Column);
            columnReverseHash.put(ms1Column,MS1_FEATURES_COLUMN_NAME);

            TableColumn ms1FeatureCoverageColumn = new TableColumn();
            ms1FeatureCoverageColumn.setHeaderValue(
                    TextProvider.getText(MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME));
            ms1FeatureCoverageColumn.setModelIndex(2);
            columnHash.put(MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME,
                           ms1FeatureCoverageColumn);
            columnReverseHash.put(ms1FeatureCoverageColumn,
                                  MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME);

            TableColumn ms2Column = new TableColumn();
            ms2Column.setHeaderValue(TextProvider.getText(MS2_FEATURES_COLUMN_NAME));
            ms2Column.setModelIndex(3);
            columnHash.put(MS2_FEATURES_COLUMN_NAME, ms2Column);
            columnReverseHash.put(ms2Column,MS2_FEATURES_COLUMN_NAME);

            TableColumn ms1ms2Column = new TableColumn();
            ms1ms2Column.setHeaderValue(TextProvider.getText(MS1_MS2_FEATURES_COLUMN_NAME));
            ms1ms2Column.setModelIndex(4);
            columnHash.put(MS1_MS2_FEATURES_COLUMN_NAME, ms1ms2Column);
            columnReverseHash.put(ms1ms2Column,MS1_MS2_FEATURES_COLUMN_NAME);

            columnModel = new DefaultTableColumnModel();
            columnModel.addColumn(proteinColumn);
            columnModel.addColumn(ms1Column);
            columnModel.addColumn(ms1FeatureCoverageColumn);
        }

        public void showMS2Columns()
        {
            //just in case it's already there, don't show two copies
            columnModel.removeColumn(columnHash.get(MS2_FEATURES_COLUMN_NAME));
            columnModel.removeColumn(columnHash.get(MS1_MS2_FEATURES_COLUMN_NAME));
            columnModel.addColumn(columnHash.get(MS2_FEATURES_COLUMN_NAME));
            columnModel.addColumn(columnHash.get(MS1_MS2_FEATURES_COLUMN_NAME));
        }

        public void hideMS2Columns()
        {
            columnModel.removeColumn(columnHash.get(MS2_FEATURES_COLUMN_NAME));
            columnModel.removeColumn(columnHash.get(MS1_MS2_FEATURES_COLUMN_NAME));
        }

        public void setDisplayedProteins(ArrayList<Protein> displayedProteins)
        {
            _displayedProteins = displayedProteins;
            fireTableDataChanged();
        }

        public int getRowCount()
        {
            if (_displayedProteins == null)
                return 0;
            return _displayedProteins.size();
        }

        public Class getColumnClass(int c)
        {
            String columnName = getNameForColumn(c);
            if (PROTEIN_COLUMN_NAME.equals(columnName))
                return String.class;
            else if (MS2_FEATURES_COLUMN_NAME.equals(columnName))
                return Integer.class;
            else if (MS1_FEATURES_COLUMN_NAME.equals(columnName))
                return Integer.class;
            else if (MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME.equals(columnName))
                return String.class;
            else if (MS1_MS2_FEATURES_COLUMN_NAME.equals(columnName))
                return Integer.class;
            else
                return Object.class;
        }

        public boolean isCellEditable(int row, int col)
        {
            return false;
        }

        protected String getNameForColumn(int col)
        {
            return columnReverseHash.get(columnModel.getColumn(col));
        }

        public Object getValueAt(int row, int col)
        {

            if (row >= _displayedProteins.size())
                return null;

            Protein protein = _displayedProteins.get(row);
            String columnName = getNameForColumn(col);
            if (PROTEIN_COLUMN_NAME.equals(columnName))
                return protein.getLookup();
            else if (MS2_FEATURES_COLUMN_NAME.equals(columnName))
                return getMatchedMS2Features(protein).size();
            else if (MS1_FEATURES_COLUMN_NAME.equals(columnName))
                return getMatchedMS1Features(protein).size();
            else if (MS1_PERCENT_FEATURE_COVERAGE_COLUMN_NAME.equals(columnName))
                return round(ProteinDisplay.calculatePercentCovered(
                        getMatchedMS1Features(protein),
                        protein.getSequenceAsString()),2) + "%";
            else if (MS1_MS2_FEATURES_COLUMN_NAME.equals(columnName))
                return getMatchedMS1AndMS2Features(protein).size();
            else
                return "";
        }

        public Protein getSelectedProtein(int i)
        {
            return _displayedProteins.get(i);
        }

    }

    /**
     * tablemodel for the peptide display table
     */
    private class FeatureTableModel extends AbstractTableModel
    {
        protected ArrayList<Feature> _displayedFeatures = new ArrayList<Feature>();
        protected String[] _displayedFeatureHtmlStrings = null;

        public void setDisplayedFeatures(ArrayList<Feature> displayedFeatures)
        {
            _displayedFeatures = displayedFeatures;
            if (_displayedFeatures != null)
            {
                Collections.sort(_displayedFeatures,new Feature.MassAscComparator());
                _displayedFeatureHtmlStrings =
                        ProteinDisplay.getHtmlForFeatures(_displayedFeatures,
                                                  _currentProteinMS1AndMS2MatchedFeatures,
                                                  _selectedFeatures, _featureDisplayMode);
            }
            fireTableDataChanged();
        }

        public int getRowCount()
        {
            return _displayedFeatures.size();
        }

        public int getColumnCount()
        {
            return columnNames.length;
        }

        public Class getColumnClass(int c)
        {
            switch (c)
            {
                case 0:
                    return Object.class;
                case 1:
                    return Integer.class;
                case 2:
                    return Float.class;
                case 3:
                    return Float.class;
                case 4:
                    return Integer.class;
                default:
                    return Object.class;
            }
        }

        public boolean isCellEditable(int row, int col)
        {
            return false;
        }

        public Object getValueAt(int row, int col)
        {
            if (row >= _displayedFeatures.size())
                return null;

            Feature feature = _displayedFeatures.get(row);
            switch (col)
            {
                case 0:
                    return _displayedFeatureHtmlStrings[row];
//                    return feature.getPeptide();
                case 1:
                    return feature.getScan();
                case 2:
                    return feature.getMass();
                case 3:
                    return feature.getTotalIntensity();
                default:
                    return "";
            }
        }

        private String[] columnNames = new String[]{
                TextProvider.getText("PEPTIDE"),
                TextProvider.getText("SCAN"),
                TextProvider.getText("MASS"),
                TextProvider.getText("TOTAL_INTENSITY"),
                TextProvider.getText("MASS_MATCHES")
        };

        public String getColumnName(int i)
        {
            return i >= 0 && i < columnNames.length ? columnNames[i] : "";
        }

        public ArrayList<Feature> getSelectedFeatures(int[] indexes)
        {
            if (indexes == null || indexes.length == 0)
                return null;

            ArrayList<Feature> featureList = new ArrayList<Feature>(indexes.length);

            for (int i=0; i<indexes.length; i++)
            {
                featureList.add(_displayedFeatures.get(indexes[i]));
            }
            return featureList;
        }



    }

    /**
     *  menu action for saving a fasta file containing the currently displayed proteins
     */
    public class SaveFastaAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            if (_displayedProteins == null)
            {
                ApplicationContext.infoMessage(TextProvider.getText("NO_PROTEINS_TO_SAVE"));
                return;
            }

            WorkbenchFileChooser chooser = new WorkbenchFileChooser();
            int chooserStatus = chooser.showOpenDialog(ApplicationContext.getFrame());
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            File outFastaFile = chooser.getSelectedFile();

            try
            {
                Protein.saveProteinArrayToFasta(_displayedProteins.toArray(new Protein[0]),
                                                outFastaFile);
                ApplicationContext.infoMessage(TextProvider.getText("SAVED_FASTA_FILE"));
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage(
                        TextProvider.getText("ERROR_SAVING_FASTA"), e);
            }
        }
    }


    /**
     *  menu action for loading an ms2 feature file for regression
     */
    public class SetParametersAction  extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            mParametersDialog = new ParametersDialog();
            mParametersDialog.setVisible(true);
        }
    }
    //for accessing this action through CommandFileRunner
    public SetParametersAction getSetParametersAction()
    {
        return new SetParametersAction();
    }



    /**
     *  menu action for loading an ms2 feature file for regression
     */
    public class LoadMS2RegressionAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            loadMS2RegressionFeatures();
        }
    }

    /**
     * menu action for loading a fasta file for protein matching
     */
    public class MatchProteinsFromFastaAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            _fastaFile = WorkbenchFileChooser.chooseExistingFile("CHOOSE_FASTA_FILE");
            if (_fastaFile == null)
                return;
            findMatchedMS1Proteins();
        }
    }

    /**
     *  menu action for loading an ms2 feature file for regression
     */
    public class ShowPredictedHydroTimePlotAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            showPredictedHydroTimePlot();
        }
    }


    /**
     *  menu action for loading an ms2 feature file for regression
     */
    public class FindProteinsMatchedInMS2Action extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            File matchingFeatureFile =
                    WorkbenchFileChooser.chooseExistingFile("MATCHING_MS2_FILE");
            if (matchingFeatureFile == null)
                return;
            findMatchedMS2Proteins(matchingFeatureFile);

        }
    }

    /**
     *  menu action for
     */
    public class LoadOffsetFeatureSetAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            File offsetFeatureFile =
                    WorkbenchFileChooser.chooseExistingFile("OFFSET_FEATURE_FILE");
            if (offsetFeatureFile == null)
                return;
            loadOffsetsFromFeatureFile(offsetFeatureFile);
        }
    }


    /**
     * menu action for writing out a pepxml file containing all the features matched
     * between ms1 and ms2 feature sets
     */
/*
    public class WritePepXmlMatchedMs1Ms2Action extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            WorkbenchFileChooser chooser = new WorkbenchFileChooser();
            chooser.setDialogTitle(TextProvider.getText("OUTPUT_PEPXML_FILE"));
            chooser.showOpenDialog(ApplicationContext.getFrame());
            File outPepXmlFile = chooser.getSelectedFile();
            writePepXmlMs1Ms2Matches(outPepXmlFile);
        }
    }
*/

    public class WritePepXmlMs2RegressionAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent evt)
        {
            WorkbenchFileChooser chooser = new WorkbenchFileChooser();
            chooser.setDialogTitle(TextProvider.getText("OUTPUT_PEPXML_FILE"));
            int chooserStatus = chooser.showOpenDialog(ApplicationContext.getFrame());
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;                
            File outPepXmlFile = chooser.getSelectedFile();
            writePepXmlMs2RegressionOffsetFile(outPepXmlFile);
        }
    }


        /**
         * context menu for pulling up protein information
         */
        protected class ProteinTablePopupMenu extends JPopupMenu
        {
            JMenuItem menuItemNCBI;

            public ProteinTablePopupMenu()
            {
                super();

                menuItemNCBI =
                        new JMenuItem(TextProvider.getText("NCBI_WEB_LOOKUP_SELECTED_PROTEIN"));

                ListenerHelper helper = new ListenerHelper(this);
                helper.addListener(menuItemNCBI,"menuItemNCBI_actionPerformed");

                //layout
                add(menuItemNCBI);
            }

            public void menuItemNCBI_actionPerformed(ActionEvent event)
            {
                if (_selectedProtein != null)
                    ProteinDisplay.openNCBIBrowserWindow(_selectedProtein);
            }


        }

    /**
     * Compare proteins based on how many MS1 matches they have.  The one with fewer
     * matches is deemed the greater.  This is for display in the protein table
     */
    protected class ProteinMS1MatchesComparator implements Comparator
    {
         public int compare(Object o1, Object o2)
         {
            int o1Features = getMatchedMS1Features((Protein) o1).size();
            int o2Features = getMatchedMS1Features((Protein) o2).size();
            return o1Features < o2Features ? 1 : o2Features < o1Features ? -1 : 0;
         }
    }

    /*
     * Compare proteins based on percent feature coverage.  The one with fewer
     * matches is deemed the greater.  This is for display in the protein table
     */
    protected class ProteinMS1PercentCoverageComparator implements Comparator
    {
         public int compare(Object o1, Object o2)
         {
            double o1Percent = _proteinPercentCoverageMap.get(o1);
            double o2Percent = _proteinPercentCoverageMap.get(o2);

            return o1Percent < o2Percent ? 1 : o2Percent < o1Percent ? -1 : 0;
         }
    }



    /**
     * dialog box for specifying parameters.  This is also operated by CommandFileRunner
     */
    public class ParametersDialog extends JDialog
    {
        public JTextField textDeltaMass;
        public JComboBox comboBoxDaPpm;
        public JTextField textDeltaScanTime;
        public JComboBox comboBoxScanTime;
        public JTextField textMinMatchedFeatures;
        public JTextField textMinPercentFeatureCoverage;
        public JTextField textMinPeptideProphet;

        public JButton buttonCancel;
        public JButton buttonOK;

        protected ProteinMatcherFrame _proteinMatcherFrame = null;

        public ParametersDialog()
        {
            super(ApplicationContext.getFrame(), TextProvider.getText("SET_PARAMETERS"));

            //graphical stuff
            Container contentPanel = null;
            try
            {
                contentPanel =
                    Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/AMTParametersDialog.xml", this);
                setContentPane(contentPanel);
                pack();
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage(
                    TextProvider.getText("ERROR_CREATING_DIALOG"), x);
                throw new RuntimeException(x);
            }

            //TODO: should really use TextProvider here, and use an internal value for determining state
            comboBoxDaPpm.addItem("Daltons");
            comboBoxDaPpm.addItem("PPM");
            comboBoxScanTime.addItem("scans");
            comboBoxScanTime.addItem("seconds");

            textDeltaMass.setText(Double.toString(round((double)deltaMass,2)));

            if (deltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_PPM)
                comboBoxDaPpm.setSelectedItem(comboBoxDaPpm.getItemAt(1));
            else
                comboBoxDaPpm.setSelectedItem(comboBoxDaPpm.getItemAt(0));

            if (_scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME)
            {
                comboBoxScanTime.setSelectedItem(comboBoxScanTime.getItemAt(1));
                textDeltaScanTime.setText(Double.toString(round(deltaTime,2)));
            }
            else
            {
                comboBoxScanTime.setSelectedItem(comboBoxScanTime.getItemAt(0));
                textDeltaScanTime.setText(Integer.toString(deltaScan));
            }

            textMinMatchedFeatures.setText(Integer.toString(minMatchedFeatures));
            textMinPercentFeatureCoverage.setText(Integer.toString(minPercentFeatureCoverage));
            textMinPeptideProphet.setText(Double.toString(round((double) minPeptideProphet,2)));

            ListenerHelper helper = new ListenerHelper(this);
            helper.addListener(buttonCancel,"buttonCancel_actionPerformed");
            helper.addListener(buttonOK,"buttonOK_actionPerformed");

            getRootPane().setDefaultButton(buttonOK);
        }

        public void buttonOK_actionPerformed(ActionEvent event)
        {
            //one big try block to catch all parsing problems.  This could
            //be broken up.
            //Assign everything to temporary variables.  Then, if everything is parsed
            //ok, assign all the real variables
            try
            {
                float tempDeltaMass = (float) Double.parseDouble(textDeltaMass.getText());
                int tempDeltaMassType = FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE;
                if ("PPM".equals(comboBoxDaPpm.getSelectedItem()))
                    tempDeltaMassType = FeatureSetMatcher.DELTA_MASS_TYPE_PPM;
                int tempMinMatchedFeatures =
                        (int) Double.parseDouble(textMinMatchedFeatures.getText());
                int tempMinPercentFeatureCoverage =
                        (int) Double.parseDouble(textMinPercentFeatureCoverage.getText());
                float tempMinPeptideProphet =
                        (float) Double.parseDouble(textMinPeptideProphet.getText());

                String comboBoxValue = (String) comboBoxScanTime.getSelectedItem();

                //TODO: the actual strings "elution" and "scan" should be internal values, not displayed
                int tempScanOrTimeMode = ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN;
                double tempDeltaTime = deltaTime;
                int tempDeltaScan = deltaScan;
                if ("seconds".equals(comboBoxValue))
                {
                    tempScanOrTimeMode = ProteomicsRegressionUtilities.REGRESSION_MODE_TIME;
                    tempDeltaTime = (float) Double.parseDouble(textDeltaScanTime.getText());
                }
                else if ("scans".equals(comboBoxValue))
                {
                    tempScanOrTimeMode = ProteomicsRegressionUtilities.REGRESSION_MODE_SCAN;
                    tempDeltaScan = (int) Double.parseDouble(textDeltaScanTime.getText());
                }

                //if we got here, we're ok, assign the real variables
                deltaMass = tempDeltaMass;
                deltaMassType = tempDeltaMassType;
                minMatchedFeatures = tempMinMatchedFeatures;
                minPercentFeatureCoverage = tempMinPercentFeatureCoverage;
                minPeptideProphet = tempMinPeptideProphet;
                _scanOrTimeMode = tempScanOrTimeMode;
                if (_scanOrTimeMode == ProteomicsRegressionUtilities.REGRESSION_MODE_TIME)
                    deltaTime = tempDeltaTime;
                else
                    deltaScan = tempDeltaScan;

                if (_fastaFile != null)
                    findMatchedMS1Proteins();
                if (_matchingMS2Features != null)
                        findMatchedMS2Proteins();

                this.setVisible(false);
                this.dispose();
            }
            catch (NumberFormatException e)
            {
                ApplicationContext.errorMessage(TextProvider.getText("BAD_PARAMETER_ERROR_MESSAGE"), e);
            }
        }

        public void buttonCancel_actionPerformed(ActionEvent e)
        {
            this.setVisible(false);
            this.dispose();
        }
    }





/** can't find a good way to use this yet
    // This comparator is used to sort vectors of data
    public static class ColumnSorter implements Comparator {
        int colIndex;
        boolean ascending;
        ColumnSorter(int colIndex, boolean ascending) {
            this.colIndex = colIndex;
            this.ascending = ascending;
        }
        public int compare(Object a, Object b) {
            Vector v1 = (Vector)a;
            Vector v2 = (Vector)b;
            Object o1 = v1.get(colIndex);
            Object o2 = v2.get(colIndex);

            // Treat empty strains like nulls
            if (o1 instanceof String && ((String)o1).length() == 0) {
                o1 = null;
            }
            if (o2 instanceof String && ((String)o2).length() == 0) {
                o2 = null;
            }

            // Sort nulls so they appear last, regardless
            // of sort order
            if (o1 == null && o2 == null) {
                return 0;
            } else if (o1 == null) {
                return 1;
            } else if (o2 == null) {
                return -1;
            } else if (o1 instanceof Comparable) {
                if (ascending) {
                    return ((Comparable)o1).compareTo(o2);
                } else {
                    return ((Comparable)o2).compareTo(o1);
                }
            } else {
                if (ascending) {
                    return o1.toString().compareTo(o2.toString());
                } else {
                    return o2.toString().compareTo(o1.toString());
                }
            }
        }
    }

    // Regardless of sort order (ascending or descending), null values always appear last.
    // colIndex specifies a column in model.
    public static void sortAllRowsBy(DefaultTableModel model, int colIndex, boolean ascending) {
        Vector data = model.getDataVector();
        Collections.sort(data, new ColumnSorter(colIndex, ascending));
        model.fireTableStructureChanged();
    }
*/
    //TODO: this should really be a public static method in a common class somewhere
    protected static final double factors[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

    protected static double round(double n, int places)
    {
        if (places < 0 || places > 8) places = 8;

        return Math.round(n * factors[places]) / factors[places];
    }


}

