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

import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFindingBroker;
import org.fhcrc.cpl.viewer.feature.extraction.PeakExtractor;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.commandline.modules.FindPeptidesCommandLineModule;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.*;

/**
 * User: mbellew
 * Date: May 26, 2004
 * Time: 2:58:14 PM
 *
 * dhmay changing 2009/01/02: adding default filtering of minimum 2 peaks and maximum KL 3.0
 *
 */
public class ExtractFeatureRangesAction extends AbstractAction implements PropertyChangeListener
{
    public ExtractFeatureRangesAction()
    {
        super(TextProvider.getText("FIND_ALL_FEATURES_DOTDOTDOT"));
        ApplicationContext.addPropertyChangeListener(SharedProperties.MS_RUN, this);
    }


    public void actionPerformed(ActionEvent event)
    {
        MSRun run;
        run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
        if (null == run)
            return;

        FeatureScanner scanner = new FeatureScanner(run, null, true);

        ExtractFeatureRangesDialog dialog = new ExtractFeatureRangesDialog(scanner);
        dialog.setModal(true);
        dialog.setVisible(true);
        if (dialog.action == ExtractFeatureRangesDialog.ACTION_CANCEL || null == scanner._file)
            return;

        File file = scanner._file;
        if (dialog.action == ExtractFeatureRangesDialog.ACTION_LOAD)
        {
            if (file.exists())
            {
                FeatureSelectionFrame.FeatureSelectionDialog.getInstance().addFeatureSet(file);
                FeatureSelectionFrame.FeatureSelectionDialog.getInstance().setVisible(true);
            }
            else
                ApplicationContext.errorMessage(TextProvider.getText("FILE_FILE_DOES_NOT_EXIST","FILEPATH", file.getPath()), null);
        }
        else
        {
            //Don't have the information cached in a file, so do this on a background thread..
            Thread t = new Thread(scanner);
            t.setPriority(Thread.MIN_PRIORITY);
            t.start();
        }
    }


    public static File getFeatureRangeFile(MSRun run)
    {
        File runFile = run.getFile();
        return getFeatureRangeFile(runFile);
    }

    public static File getFeatureRangeFile(File runFile)
    {
        return new File(runFile.getAbsolutePath() + ".peptideRanges.tsv");
    }

    public static File getFeatureFile(MSRun run)
    {
        File featureRangeFile = run.getFile();
        String path = featureRangeFile.getAbsolutePath();
        if (path.toLowerCase().endsWith(".mzxml"))
            path = path.substring(0, path.length() - 6);
        String featureRangeFileName = path + ".peptides.tsv";

        return new File(featureRangeFileName);
    }

    public static File getFeatureFile(File runFile)
    {
        String featureRangeFileName = runFile.getAbsolutePath() + ".peptides.tsv";

        return new File(featureRangeFileName);
    }


    public static Feature[] ExtractScanFeatures(MSRun.MSScan scan, MSRun run)
    {
        try
        {
            FloatRange range = FeatureExtractor.getMzExtractionRange(scan);
            FeatureSet features = FeatureExtractor.getDefault(run, scan.getNum(), 1, 6, range, 2.0).analyze();
            return features.getFeatures();
        }
        catch (InterruptedException x)
        {
            return null;
        }
    }


    public void propertyChange(PropertyChangeEvent event)
    {
        setEnabled(null != ApplicationContext.getProperty(SharedProperties.MS_RUN));
    }


    public static class FeatureScanner implements Runnable
    {
        protected MSRun _run;
        protected boolean _updateUI = true;
        protected File _file = null;
        protected int _startScan = 1;
        protected int _endScan = -1;
        protected boolean _outputZeroCharge = false;
        protected boolean _deconvolute = false;
        protected boolean _quant = false;
        protected float _lightTagWeight = -1.f;
        protected float _heavyTagWeight = -1.f;
        protected char _labeledResidue = ' ';
        protected int _maxLabelCount = AnalyzeICAT.DEFAULT_MAX_LABEL_COUNT;

        public FeatureScanner(MSRun run, File file)
        {
            _run = run;
            _file = file;
            _endScan = run.getScanCount();
            _startScan = 1;
            // Default label params for ICAT
            _lightTagWeight = AnalyzeICAT.icatLabel.getLight();
            _heavyTagWeight = AnalyzeICAT.icatLabel.getHeavy() + _lightTagWeight; // getHeavy() gives delta
            _labeledResidue = AnalyzeICAT.icatLabel.getResidue();
            _maxLabelCount = AnalyzeICAT.icatLabel.getMaxLabelCount();
        }


        public FeatureScanner(MSRun run, File file, boolean updateUI)
        {
            this(run, file);
            _updateUI = updateUI;
        }


        public FeatureScanner(MSRun run, File file, boolean updateUI, int startScanNum, int endScanNum, boolean outputZeroCharge)
        {
            this(run, file, updateUI);
            setStartScan(startScanNum);
            setEndScan(endScanNum);
            _outputZeroCharge = outputZeroCharge;
        }


        public void run()
        {
            PrintWriter pw = null;

            FeatureSet fs = null;

            try
            {
                pw = new PrintWriter(new FileOutputStream(_file));
                pw.flush(); // test printwriter

                FloatRange range = FeatureExtractor.getMzExtractionRange(_run);

                int scanCount = _endScan - _startScan + 1;
                int maxCharge = PeakCombiner.DEFAULT_MAX_CHARGE;

                Class featureStrategyClass = FeatureExtractor.getDefaultClass();

                fs = FeatureFindingBroker.findPeptides(_run, _startScan, scanCount,
                        maxCharge, range, 0,
                        FeatureFinder.DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS,
                        featureStrategyClass, true,
                        PeakExtractor.DEFAULT_PEAK_RIDGE_WALK_SMOOTHED, false);
                //default filtering.  No way to turn this off in the GUI
                FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
                sel.setMinPeaks(FindPeptidesCommandLineModule.DEFAULT_MIN_PEAKS);
                sel.setMaxKL(FindPeptidesCommandLineModule.DEFAULT_MAX_KL);
                fs = fs.filter(sel);

                if (!_outputZeroCharge)
                {
                    Feature[] features = fs.getFeatures();
                    ArrayList collect = new ArrayList(features.length);
                    for (int i = 0; i < features.length; i++)
                    {
                        if (features[i].charge == 0)
                            continue;
                        collect.add(features[i]);
                    }
                    features = (Feature[])collect.toArray(new Feature[0]);
                    fs.setFeatures(features);
                }

                if (_deconvolute)
                {
                    // TODO: add to GUI as parameters
                    int scanWindow = 6;
                    double massWindow = .4;
                    boolean sumIntensities = true;
                    
                    fs = fs.deconvolute(scanWindow, massWindow, sumIntensities);
                }

//dhmay changing
                if (_quant)
                    fs = fs.quant(_lightTagWeight, _heavyTagWeight, _labeledResidue, _maxLabelCount, _run);

                fs.save(pw);

                pw.close();
                pw = null;
            }
            catch (Exception x1)
            {
                ApplicationContext.errorMessage(TextProvider.getText("ERROR_WRITING_FEATURE_FILE"), x1);
                return;
            }
            catch (java.lang.OutOfMemoryError ome)
            {
                ApplicationContext.errorMessage(TextProvider.getText("FEATURE_SCANNER_OME"), ome);
                return;
            }
            finally
            {
                if (null != pw)
                    pw.close();
            }

            if (_updateUI)
            {
                Arrays.sort(fs.getFeatures(), new Feature.MzScanAscComparator());
                fs.setColor(FeatureSelectionFrame.FeatureSelectionDialog.nextColor());
                final FeatureSet finalSet = fs;
                fs.setSourceFile(_file);
                EventQueue.invokeLater(new Runnable()
                {
                    public void run()
                    {
                        FeatureSelectionFrame.FeatureSelectionDialog.getInstance().addFeatureSet(finalSet);
                    }
                });
            }
            ApplicationContext.setMessage(TextProvider.getText("FINDING_FEATURES_COMPLETE_SEE_FILE","FILEPATH",_file.getAbsolutePath()));
        }

        public int getStartScan()
        {
            return _startScan;
        }

        public void setStartScan(int startScanNum)
        {
            _startScan = _run.getIndexForScanNum(startScanNum);
            if (_startScan < 0)
                _startScan = -(_startScan + 1);
        }

        public int getEndScan()
        {
            return _endScan;
        }

        public void setEndScan(int endScanNum)
        {
            _endScan = _run.getIndexForScanNum(endScanNum);
            if (_endScan < 0)
                _endScan = -(_endScan + 1);
            _endScan = Math.min(_run.getScanCount(),_endScan);
        }

        public void setOutputZeroCharge(boolean outputZeroCharge)
        {
            this._outputZeroCharge = outputZeroCharge;
        }

        public void setDeconvolute(boolean deconvolute)
        {
            this._deconvolute = deconvolute;
        }

        public void setLightTagWeight(float lightTagWeight)
        {
            this._lightTagWeight = lightTagWeight;
        }
        public void setHeavyTagWeight(float heavyTagWeight)
        {
            this._heavyTagWeight = heavyTagWeight;
        }
        public void setLabeledResidue(char labeledResidue)
        {
            this._labeledResidue = labeledResidue;
        }
        public void setMaxLabelCount(int maxLabelCount)
        {
            this._maxLabelCount = maxLabelCount;
        }
    }


    public static class ExtractFeatureRangesDialog extends JDialog
    {
        public static final int ACTION_CANCEL = 0;
        public static final int ACTION_LOAD = 1;
        public static final int ACTION_FIND = 2;

        public int action = ACTION_CANCEL;

        /*		public String fileName = null;
          public int startScan = 1;
          public int endScan = 2000;
          public boolean outputZeroCharge = false;
          public boolean icat = false;
          public boolean deconvolute = false; */
        FeatureScanner scanner = null;

        public JLabel labelAlgorithm;
        public JButton buttonCancel;
        public JButton buttonLoadFile;
        public JButton buttonFindFeatures;
        public JButton buttonBrowse;
        public JTextField textFileName;
        public JTextField textEndScan;
        public JTextField textStartScan;
        public JCheckBox checkboxZeroCharge;
        public JCheckBox checkboxDeconvolute;
        public JCheckBox checkboxQuant;
        public JComboBox comboboxTagSelect;
        public JTextField textLightTagWeight;
        public JTextField textHeavyTagWeight;
        public JTextField textLabeledResidue;
        public JTextField textMaxLabelCount;

        public ExtractFeatureRangesDialog(FeatureScanner scanner)
        {
            super(ApplicationContext.getFrame(), TextProvider.getText("EXTRACT_FEATURES"));
            this.scanner = scanner;

            try
            {
                Container content = Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/ExtractFeatureRanges.xml",this);
                setContentPane(content);
                pack();
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
                throw new RuntimeException(x);
            }

            final MSRun run = (MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null == run)
                throw new IllegalStateException();

            String start = "1";
            String end = "1";
            if (run.getScanCount() > 0)
            {
                start = String.valueOf(run.getScan(0).getNum());
                end = String.valueOf(run.getScan(run.getScanCount() - 1).getNum());
            }

            Class<?> c = (Class<?>)ApplicationContext.getProperty(FeatureExtractor.DEFAULT_EXTRACTOR_PROPERTYNAME);
            String name = c.getName();
            name = name.substring(name.lastIndexOf('.')+1);
            labelAlgorithm.setText(name);
            textStartScan.setText(start);
            textEndScan.setText(end);
            populateTagSelectComboBox(comboboxTagSelect);
            comboboxTagSelect.setEnabled(false);
            textLightTagWeight.setText(String.valueOf(scanner._lightTagWeight));
            textLightTagWeight.setEnabled(false);
            textHeavyTagWeight.setText(String.valueOf(scanner._heavyTagWeight));
            textHeavyTagWeight.setEnabled(false);
            textLabeledResidue.setText(String.valueOf(scanner._labeledResidue));
            textLabeledResidue.setEnabled(false);
            textMaxLabelCount.setText(String.valueOf(scanner._maxLabelCount));
            textMaxLabelCount.setEnabled(false);


            if (null == scanner._file)
                scanner._file = new File(ExtractFeatureRangesAction.getFeatureFile(run).getAbsolutePath());
            textFileName.setText(scanner._file.getPath());

            buttonLoadFile.setEnabled(scanner._file.exists());

            ListenerHelper helper = new ListenerHelper(this);
            helper.addListener(buttonCancel, "buttonCancel_actionPerformed");
            helper.addListener(buttonFindFeatures, "buttonFindFeatures_actionPerformed");
            helper.addListener(buttonLoadFile, "buttonLoadFile_actionPerformed");
            helper.addListener(buttonBrowse, "buttonBrowse_actionPerformed");
            helper.addListener(textFileName, "textFileName_keyTyped");
            helper.addListener(textFileName, "textFileName_actionPerformed");
            helper.addListener(textFileName.getDocument(), "textFileName_insertUpdate");
            helper.addListener(textFileName.getDocument(), "textFileName_removeUpdate");
            helper.addListener(textFileName.getDocument(), "textFileName_changedUpdate");
            helper.addListener(checkboxQuant, "checkboxQuant_actionPerformed");
            helper.addListener(comboboxTagSelect, "comboboxTagSelect_actionPerformed");
        }


        private void updateButtonState()
        {
            scanner._file = new File(textFileName.getText());
            buttonLoadFile.setEnabled(scanner._file.exists());
        }


        public void buttonCancel_actionPerformed(ActionEvent event)
        {
            action = ACTION_CANCEL;
            ExtractFeatureRangesDialog.this.dispose();
        }


        public void buttonFindFeatures_actionPerformed(ActionEvent event)
        {
            action = ACTION_FIND;
            try
            {
                scanner.setEndScan(Integer.parseInt(textEndScan.getText()));
                scanner.setStartScan(Integer.parseInt(textStartScan.getText()));
                scanner._outputZeroCharge = checkboxZeroCharge.isSelected();
                scanner._deconvolute = checkboxDeconvolute.isSelected();
                scanner._quant = checkboxQuant.isSelected();
                scanner.setLightTagWeight(Float.parseFloat(textLightTagWeight.getText()));
                scanner.setHeavyTagWeight(Float.parseFloat(textHeavyTagWeight.getText()));
                scanner.setLabeledResidue(textLabeledResidue.getText().charAt(0));
                try
                {
                    scanner.setMaxLabelCount(Integer.parseInt(textMaxLabelCount.getText()));
                }
                catch (Exception e)
                {
                    ApplicationContext.infoMessage(TextProvider.getText("BAD_VALUE_MAX_LABEL_USING_DEFAULT",
                            Integer.toString(AnalyzeICAT.DEFAULT_MAX_LABEL_COUNT)));
                    scanner.setMaxLabelCount(AnalyzeICAT.DEFAULT_MAX_LABEL_COUNT);
                }
            }
            catch (NumberFormatException x)
            {
            }
            ExtractFeatureRangesDialog.this.dispose();
        }


        public void buttonLoadFile_actionPerformed(ActionEvent event)
        {
            action = ACTION_LOAD;
            ExtractFeatureRangesDialog.this.dispose();
        }

        public void buttonBrowse_actionPerformed(ActionEvent event)
        {
            JFileChooser chooser = new WorkbenchFileChooser();
            int chooserStatus = chooser.showOpenDialog(ExtractFeatureRangesDialog.this);
            //if user didn't hit OK, ignore
            if (chooserStatus != JFileChooser.APPROVE_OPTION)
                return;
            File file = chooser.getSelectedFile();
            if (null != file)
            {
                textFileName.setText(file.getAbsolutePath());
            }
        }

        public void textFileName_keyTyped(KeyEvent event)
        {
            updateButtonState();
        }

        public void textFileName_actionPerformed(ActionEvent event)
        {
            updateButtonState();
        }

        public void textFileName_insertUpdate(DocumentEvent event)
        {
            updateButtonState();
        }

        public void textFileName_removeUpdate(DocumentEvent event)
        {
            updateButtonState();
        }

        public void textFileName_changedUpdate(DocumentEvent event)
        {
            updateButtonState();
        }

        public void checkboxQuant_actionPerformed(ActionEvent event)
        {
            boolean enabled = checkboxQuant.isSelected();
            textLightTagWeight.setEnabled(enabled);
            textHeavyTagWeight.setEnabled(enabled);
            comboboxTagSelect.setEnabled(enabled);
            textLabeledResidue.setEnabled(enabled);
            textMaxLabelCount.setEnabled(enabled);
        }

        /**
         * event handler for the prepopulated tag combobox.
         * Not very efficient at matching the selected label, but that list will never be big,
         * because it's displayed in the UI.
         * @param event
         */
        public void comboboxTagSelect_actionPerformed(ActionEvent event)
        {
            ArrayList tagArrayList = AnalyzeICAT.getPrePopulatedLabels();
            ListIterator tagArrayListIterator = tagArrayList.listIterator();
            String selectedTagName = (String) comboboxTagSelect.getSelectedItem();
            while (tagArrayListIterator.hasNext())
            {
                AnalyzeICAT.IsotopicLabel currentLabel =
                        (AnalyzeICAT.IsotopicLabel) tagArrayListIterator.next();
                if (currentLabel.getName().equals(selectedTagName))
                {
                    textLightTagWeight.setText(String.valueOf(currentLabel.getLight()));
                    textHeavyTagWeight.setText(String.valueOf(currentLabel.getHeavy() + currentLabel.getLight()));
                    textLabeledResidue.setText(String.valueOf(currentLabel.getResidue()));
                    textMaxLabelCount.setText(String.valueOf(currentLabel.getMaxLabelCount()));
                    break;
                }
            }
        }
    }

    /**
     * Populates the label name selection combobox from the prepopulated labels stored in
     * AnalyzeICAT
     * @param comboBox
     */
    protected static void populateTagSelectComboBox(JComboBox comboBox)
    {
        if (comboBox == null)
            return;
        ArrayList tagArrayList = AnalyzeICAT.getPrePopulatedLabels();
        ListIterator tagArrayListIterator = tagArrayList.listIterator();
        while (tagArrayListIterator.hasNext())
        {
            AnalyzeICAT.IsotopicLabel currentLabel =
                    (AnalyzeICAT.IsotopicLabel) tagArrayListIterator.next();
            comboBox.addItem(currentLabel.getName());
        }
    }
}
