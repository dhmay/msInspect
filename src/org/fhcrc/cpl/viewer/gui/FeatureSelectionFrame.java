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

import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.BucketedPeptideArray;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.proteomics.feature.*;

import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.chart.ChartPanel;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.table.*;
import javax.swing.event.ListSelectionEvent;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.text.SimpleDateFormat;
import java.text.NumberFormat;
import java.text.DecimalFormat;

public class FeatureSelectionFrame extends AbstractAction implements PropertyChangeListener
{
    protected static Logger _log = Logger.getLogger(FeatureSelectionFrame.class);

    public FeatureSelectionDialog dialog = null;


    public void propertyChange(PropertyChangeEvent event)
    {
//		setEnabled(null != ApplicationContext.getProperty(SharedProperties.MS_RUN));
    }

    public void actionPerformed(ActionEvent event)
    {
        dialog = new FeatureSelectionDialog();
        dialog.setVisible(true);
    }

    public static class FeatureSelectionDialog extends JDialog
    {
        public static final int ACTION_CANCEL = 0;
        public static final int ACTION_OK = 1;

        public int action = ACTION_CANCEL;

        public JTable tblFeatureSets;
        public JButton buttonAddFiles;
        public JButton buttonRemoveFiles;

        public JTabbedPane tabbedPane;

        public JTextField textMinMz;
        public JTextField textMaxMz;
        public JTextField textMinCharge;
        public JTextField textMaxCharge;
        public JTextField textMinScans;
        public JTextField textMaxKL;
        public JTextField textMinIntensity;
        public JTextField textMinPeaks;

        public JTextField textOutputFile;
        public JButton buttonBrowseOutputDir;
        public JTextField textPepArrayMzBucket;
        public JTextField textPepArrayScanBucket;
        public JCheckBox checkBoxNormalize;
        public JButton buttonOptimize;
        public JLabel lblNumBuckets;
        public JLabel lblExactMatches;
        public JButton buttonCreatePeptideArray;

        public JButton buttonCancel;
        public JButton buttonApply;
        public JButton buttonOK;

        public FileListTableModel model;

        //supporting components
        public static JFileChooser chooser = new WorkbenchFileChooser();
        public static JFileChooser saveChooser = new WorkbenchFileChooser();


        private static FeatureSelectionDialog _instance;

        private static Color[] featureSetColors = new Color[]{
                Color.BLUE,
                Color.GREEN,
                Color.CYAN,
                Color.MAGENTA,
                Color.ORANGE,
                Color.PINK,
                Color.RED,
                Color.YELLOW
        };

        synchronized public static FeatureSelectionDialog getInstance()
        {
            if (null == _instance)
                _instance = new FeatureSelectionDialog();
            return _instance;
        }

        static
        {
            chooser.setMultiSelectionEnabled(true);
        }

        public FeatureSelectionDialog()
        {

            super(ApplicationContext.getFrame(), TextProvider.getText("DISPLAY_FEATURES"));

            //graphical stuff

            Container contentPanel = null;
            try
            {
                contentPanel = Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/FeatureSelectionFrame.xml", this);
                setContentPane(contentPanel);
                pack();
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
                throw new RuntimeException(x);
            }

            Dimension d = contentPanel.getPreferredSize();
            setBounds(600, 100, 6 + (int) d.getWidth(), 50 + (int) d.getHeight());

            ListenerHelper helper = new ListenerHelper(this);

            buttonApply.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    FeatureSet.FeatureSelector selector = createFeatureSelector();
                    updateFeatureSelector(selector);
                }
            });
            buttonCancel.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    FeatureSelectionDialog.this.dispose();
                }
            });
            buttonOK.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    FeatureSet.FeatureSelector selector = createFeatureSelector();
                    updateFeatureSelector(selector);
                    FeatureSelectionDialog.this.dispose();
                }
            });
            helper.addListener(buttonAddFiles, "buttonAddFile_actionPerformed");
            helper.addListener(buttonRemoveFiles, "buttonRemoveFiles_actionPerformed");

            //end graphical stuff


            FeatureSet.FeatureSelector selector = (FeatureSet.FeatureSelector) ApplicationContext.getProperty("featureSelector");
            if (null == selector)
            {
                selector = new FeatureSet.FeatureSelector();
                updateFeatureSelector(selector);
            }

            textMaxKL.setText(String.valueOf(selector.getMaxKL()));
            textMinCharge.setText(String.valueOf(selector.getMinCharge()));
            textMaxCharge.setText(String.valueOf(selector.getMaxCharge()));
            textMaxMz.setText(String.valueOf(selector.getMaxMz()));
            textMinMz.setText(String.valueOf(selector.getMinMz()));
            textMinIntensity.setText(String.valueOf(selector.getMinIntensity()));
            textMinScans.setText(String.valueOf(selector.getMinScans()));

            File dir = null;
            MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
            if (null != run && null != run.getFile())
                dir = run.getFile().getParentFile();
            if (null == dir)
                dir = saveChooser.getCurrentDirectory();
            if (null != dir)
                textOutputFile.setText(dir.getAbsolutePath() + File.separatorChar + "pepArray.tsv");
            else
                textOutputFile.setText("pepArray.tsv");

            model = new FileListTableModel();
            tblFeatureSets.setModel(model);
            TableColumn displayColumn = tblFeatureSets.getColumnModel().getColumn(0);
            displayColumn.setPreferredWidth(50);
            displayColumn.setMaxWidth(50);
            TableColumn colorColumn = tblFeatureSets.getColumnModel().getColumn(2);
            colorColumn.setPreferredWidth(50);
            colorColumn.setMaxWidth(50);
            colorColumn.setCellRenderer(new ColorRenderer(false));
            colorColumn.setCellEditor(new ColorEditor());

            TableColumn browseColumn = tblFeatureSets.getColumnModel().getColumn(3);
            browseColumn.setPreferredWidth(30);
            browseColumn.setMaxWidth(30);
            final JButton button = new JButton("...");
            final TableCellRenderer browseButtonRenderer = new TableCellRenderer()
            {
                public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
                {
                    return button;
                }
            };
            browseColumn.setCellRenderer(browseButtonRenderer);
            browseColumn.setCellEditor(new BrowseEditor());

            List<FeatureSet> featureSets = (List<FeatureSet>) ApplicationContext.getProperty(SharedProperties.FEATURE_SETS);
            if (null != featureSets)
            {
                for (FeatureSet fs : featureSets)
                {
                    model.addFeatureSet(fs);
                }
            }

            buttonAddFiles.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    int chooserStatus = chooser.showOpenDialog(FeatureSelectionDialog.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;
                    File[] files = chooser.getSelectedFiles();
                    for (File file : files)
                    {
                        //String path = files[i].getAbsolutePath().toLowerCase();
                        FeatureSet featureSet = new FeatureSet(file, nextColor());
//int numWithoutComp = 0;
//for (Feature feature : featureSet.getFeatures()) if (feature.comprised == null) numWithoutComp++;
//System.err.println(numWithoutComp + " of " + featureSet.getFeatures().length + " features no comp");
                        _log.debug("Loaded " + featureSet.getFeatures().length + " features from file " +
                                   file.getAbsolutePath() + ".  Status: " + featureSet.getLoadStatus());
                        //check to see if features loaded correctly
                        if (featureSet.getLoadStatus() == FeatureSet.FEATURESET_LOAD_SUCCESS)
                        {
                            featureSet.setDisplayed(tabbedPane.getSelectedIndex() != 3);
                            addFeatureSet(featureSet);
                        }
                        else
                        {
                            ApplicationContext.infoMessage(featureSet.getLoadStatusMessage());
                        }
                    }

                }
            });

            buttonRemoveFiles.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    int[] rows = tblFeatureSets.getSelectedRows();
                    Arrays.sort(rows);
                    FileListTableModel model = (FileListTableModel) tblFeatureSets.getModel();
                    for (int i = rows.length - 1; i >= 0; i--)
                        model.removeFeatureSet(rows[i]);
                }
            });

            buttonBrowseOutputDir.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    File file = new File(textOutputFile.getText());
                    File dir = file.getParentFile();
                    if (dir.exists())
                    {
                        saveChooser.setCurrentDirectory(dir.getParentFile());
                        saveChooser.setSelectedFile(file);
                    }

                    int chooserStatus = saveChooser.showOpenDialog(FeatureSelectionDialog.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;                    
                    file = saveChooser.getSelectedFile();
                    if (null != file)
                        textOutputFile.setText(file.getAbsolutePath());
                }
            });

            buttonOptimize.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    double[] mzVal = new double[1];
                    mzVal[0] = Double.parseDouble(textPepArrayMzBucket.getText());
                    final Optimizer optimizer = new Optimizer(mzVal, new int[]{10, 25, 50, 75, 100, 150, 200, 300, 400}, getFeatureRanges());
                    optimizer.runLater = new OptimizeGrapher(optimizer);
                    Thread t = new Thread(optimizer);
                    t.start();
                }
            });

            buttonCreatePeptideArray.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    List featureSets = getFeatureSets();
                    if (null == featureSets || featureSets.size() < 2)
                    {
                        ApplicationContext.errorMessage(TextProvider.getText("TWO_FILES_TO_CREATE_PEPTIDE_ARRAY"), null);

                        return;
                    }

                    BucketedPeptideArray arr = new BucketedPeptideArray(featureSets, createFeatureSelector());
                    float mzBucketSize = Float.parseFloat(textPepArrayMzBucket.getText());
                    int scanBucketSize = Integer.parseInt(textPepArrayScanBucket.getText());
                    boolean normalize = checkBoxNormalize.isSelected();
                    arr.setElutionMode(FeatureClusterer.ELUTION_MODE_SCAN);
                    arr.setElutionBucket(scanBucketSize);
                    arr.setMassBucket(mzBucketSize);
                    arr.setNormalize(normalize);
                    File file = new File(textOutputFile.getText());
                    if (!file.getParentFile().exists())
                    {
                        ApplicationContext.errorMessage(TextProvider.getText("DIR_NOT_LEGAL_DIR_NAME",
                                file.getParentFile().getAbsolutePath()), null);
                        return;
                    }
                    arr.setOutFileName(file.getAbsolutePath());

                    Thread t = new Thread(arr);
                    t.setPriority(Thread.MIN_PRIORITY);
                    t.start();

//				//Create because user may have failed to hit apply button
//				FeatureSet.FeatureSelector sel = createFeatureSelector();
//				List ranges = generateMergedSets(featureSets, sel, false);
//				PeptideArray array = new PeptideArray(ranges);
//				array.setOutputDir(dir);
//				//TODO: Error checking
                }
            });
        }


        private static int nextColor;

        public static synchronized Color nextColor()
        {
            Color color = featureSetColors[nextColor];
            nextColor = ++nextColor % featureSetColors.length;
            return color;
        }

        public List getFeatureSets()
        {
            return (List) ApplicationContext.getProperty(SharedProperties.FEATURE_SETS);
        }

        public List getFeatureRanges()
        {
            return (List) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
        }

        public FeatureSet.FeatureSelector getFeatureSelector()
        {
            FeatureSet.FeatureSelector sel = (FeatureSet.FeatureSelector) ApplicationContext.getProperty("featureSelector");
            if (null == sel) //Should never happen...
            {
                sel = createFeatureSelector();
                ApplicationContext.setProperty("featureSelector", sel);
            }

            return sel;
        }

        private FeatureSet.FeatureSelector createFeatureSelector()
        {
            FeatureSet.FeatureSelector selector = new FeatureSet.FeatureSelector();
            //TODO: Error handling...
            selector.setMaxKL(Float.parseFloat(textMaxKL.getText()));
            selector.setMinCharge(Integer.parseInt(textMinCharge.getText()));
            selector.setMaxCharge(Integer.parseInt(textMaxCharge.getText()));
            selector.setMaxMz(Float.parseFloat(textMaxMz.getText()));
            selector.setMinMz(Float.parseFloat(textMinMz.getText()));
            selector.setMinIntensity(Float.parseFloat(textMinIntensity.getText()));
            selector.setMinPeaks(Integer.parseInt(textMinPeaks.getText()));
            selector.setMinScans(Integer.parseInt(textMinScans.getText()));

            return selector;
        }

        public void addFeatureSet(File file)
        {
            FeatureSet featureSet;
            try
            {
                featureSet = new FeatureSet(file);
            }
            catch (Exception x)
            {
                //JOptionPane.showMessageDialog(null, x.getMessage());
                ApplicationContext.errorMessage(null, x);
                return;
            }
            addFeatureSet(featureSet);
        }

        public void addFeatureSet(Feature[] features)
        {
            addFeatureSet(new FeatureSet(features));
        }

        public void addFeatureSet(FeatureSet featureSet)
        {
            featureSet.setColor(nextColor());
            _log.debug("Adding new featureset with " + featureSet.getFeatures().length +
                    " features to list of displayed featuresets");
            ((FileListTableModel) tblFeatureSets.getModel()).addFeatureSet(featureSet);
            ApplicationContext.setProperty(SharedProperties.SELECTED, featureSet);
        }

        /**
         * Updates the featureselector based on the currently loaded features
         *
         * @param sel
         */
        public void updateFeatureSelector(FeatureSet.FeatureSelector sel)
        {
            List featureSets = (List) ApplicationContext.getProperty(SharedProperties.FEATURE_SETS);
            if (null != featureSets)
            {
                List mergedFeatureSets = generateMergedSets(featureSets, sel);
                ApplicationContext.setProperty(SharedProperties.FEATURE_RANGES, mergedFeatureSets);
            }
            ApplicationContext.setProperty("featureSelector", sel);
        }

        /**
         * Generates a new merged feature list
         *
         * @param featureSets
         * @param sel
         * @return
         */
        protected List<FeatureSet> generateMergedSets(List<? extends FeatureSet> featureSets, FeatureSet.FeatureSelector sel)
        {
            List<FeatureSet> mergedFeatureSets = new ArrayList<FeatureSet>(featureSets.size());
            for (FeatureSet fs : featureSets)
            {
                FeatureSet fsMerged = fs.filter(sel);
                fsMerged.setDisplayed(fs.isDisplayed());
                fsMerged.setSourceFile(fs.getSourceFile());
                mergedFeatureSets.add(fsMerged);
            }
            return mergedFeatureSets;
        }

        public void updateFeatureSets(List<? extends FeatureSet> featureSets, boolean displayOnly)
        {
            FeatureSet.FeatureSelector sel = (FeatureSet.FeatureSelector) ApplicationContext.getProperty("featureSelector");
            if (null != sel)
            {
                List<FeatureSet> mergedFeatureSets = null;
                if (displayOnly)
                {
                    mergedFeatureSets = new ArrayList<FeatureSet>();
                    List<? extends FeatureSet> oldFeatureSets = (List<? extends FeatureSet>) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);
                    for (int i = 0; i < oldFeatureSets.size(); i++)
                    {
                        //Need to clone them or else property change support says nothing changed.
                        FeatureSet merged = (FeatureSet) oldFeatureSets.get(i).clone();
                        FeatureSet unmerged = featureSets.get(i);
                        merged.setColor(unmerged.getColor());
                        merged.setDisplayed(unmerged.isDisplayed());
                        mergedFeatureSets.add(merged);
                    }
                }
                else
                {
                    mergedFeatureSets = generateMergedSets(featureSets, sel);
                }
                ApplicationContext.setProperty(SharedProperties.FEATURE_RANGES, mergedFeatureSets);
            }
            ApplicationContext.setProperty(SharedProperties.FEATURE_SETS, featureSets);
        }

/*
*  UI event handlers
*/

        public void buttonAddFile_actionPerformed(ActionEvent event)
        {
        }


        public void buttonRemoveFiles_actionPerformed(ActionEvent event)
        {
        }

        /*
        * Inner classes
        */


        public static class Optimizer implements Runnable
        {
            //private static final int OPTIMIZE_STEPS = 15;
            //private static final int OPTIMIZE_INCREMENT = 20;

            double[] mzValues;
            int[] scanValues;
            int bestScanWindow = 50;
            double bestMzWindow = .3;
            int bestMatchMzIndex = 0;
            int bestMatchScanIndex = 0;
            //2D Arrays one subarray for each mzValue
            int[][] numBuckets = null;
            int[][] perfectMatches = null;
            private int topN = 0;
            XYSeries[] series;
            List<FeatureSet> featureSets = null;
            public boolean align = true;
            Runnable runLater = null;
            File outFile = null;
            PrintWriter out = new PrintWriter(System.out);


            public Optimizer(double[] mzValues, int[] scanValues, List<FeatureSet> featureSets)
            {
                this.mzValues = mzValues;
                this.scanValues = scanValues;
                this.featureSets = featureSets;
                series = new XYSeries[mzValues.length];
                for (int i = 0; i < series.length; i++)
                    series[i] = new XYSeries("-Points (" + mzValues[i] + ")");

                numBuckets = new int[mzValues.length][scanValues.length];
                perfectMatches = new int[mzValues.length][scanValues.length];
            }

            public Optimizer(double[] mzValues, int[] scanValues, List<FeatureSet> featureSets, File outFile)
            {
                this(mzValues, scanValues, featureSets);
                if (null == outFile)
                    return;

                this.outFile = outFile;
                try
                {
                    this.out = new PrintWriter(new FileOutputStream(outFile));
                }
                catch (Exception x)
                {
                    ApplicationContext.errorMessage(TextProvider.getText("COULD_NOT_OPEN_FILE_FILENAME", outFile.getPath()), x);
                }
            }

            public void run()
            {
                //use getFeatureRanges so that we have a filtered set.
                //Alignment was much too slow.
                if (null == featureSets || featureSets.size() < 2)
                {
                    //JOptionPane.showMessageDialog(FeatureSelectionFrame.this, "At least 2 files must be selected to create a peptide array");
                    ApplicationContext.errorMessage(TextProvider.getText("TWO_FILES_TO_CREATE_PEPTIDE_ARRAY"), null);
                    return;
                }

                FeatureSet[] alignedSets = null;
                if (align)
                {
                    List<FeatureSet> alignedList;
                    try
                    {
                        ApplicationContext.setMessage(TextProvider.getText("ALIGNING_DOTDOTDOT"));
                        Aligner aligner = new SplineAligner();
                        aligner.setBuildCharts(true);
                        Aligner.MzFeaturePairSelector pairSelector =
                                new Aligner.MzFeaturePairSelector();
                        pairSelector.setTopN(topN);
                        aligner.setFeaturePairSelector(Aligner.DEFAULT_FEATURE_PAIR_SELECTOR);
                        alignedList = aligner.alignFeatureSets(featureSets, true);
                    }
                    catch (Exception x)
                    {
                        ApplicationContext.errorMessage(TextProvider.getText("COULD_NOT_ALIGN_FEATURE_SETS"), x);
                        return;
                    }
                    alignedSets = alignedList.toArray(new FeatureSet[featureSets.size()]);
                }
                else
                    alignedSets = featureSets.toArray(new FeatureSet[featureSets.size()]);

                FeatureGrouper grouper = new FeatureGrouper();
                for (int i = 0; i < alignedSets.length; i++)
                    grouper.addSet(alignedSets[i]);

                out.println("mzBucketSize\tscanCount\tnumBuckets\tperfectMatches");
                for (int iMz = 0; iMz < mzValues.length; iMz++)
                {
                    for (int iScan = 0; iScan < scanValues.length; iScan++)
                    {
                        int scanCount = scanValues[iScan];
                        double mzBucketSize = mzValues[iMz];
                        ApplicationContext.setMessage(TextProvider.getText("CALCULATING_ARRAY_SCAN_MZ", "SCAN_COUNT", Integer.toString(scanCount), "MZ_WINDOW", Double.toString(mzBucketSize)));
                        grouper.split2D(mzBucketSize, scanCount);
                        numBuckets[iMz][iScan] = grouper.numBuckets();
                        perfectMatches[iMz][iScan] = grouper.rowsWithOneFromEach();
                        if (perfectMatches[iMz][iScan] > perfectMatches[bestMatchMzIndex][bestMatchScanIndex])
                        {
                            bestMatchMzIndex = iMz;
                            bestMatchScanIndex = iScan;
                            bestScanWindow = scanCount;
                            bestMzWindow = mzBucketSize;
                        }
                        series[iMz].add(scanCount, (double) perfectMatches[iMz][iScan]);
                        out.println(mzBucketSize + "\t" + scanCount + "\t" + numBuckets[iMz][iScan] + "\t" + perfectMatches[iMz][iScan]);
                    }
                }

                out.flush();
                if (null != this.outFile)
                    out.close();

                if (null != runLater)
                    EventQueue.invokeLater(runLater);
            }

            public int getTopN()
            {
                return topN;
            }

            public void setTopN(int topN)
            {
                this.topN = topN;
            }
            
        }


        private class FileListTableModel extends AbstractTableModel
        {
            ArrayList _featureSets = new ArrayList();

            synchronized public void addFeatureSet(FeatureSet fs)
            {
                _featureSets.add(fs);
                ArrayList changedFeatureSets = new ArrayList();
                changedFeatureSets.addAll(_featureSets);
                updateFeatureSets(changedFeatureSets, false);

                super.fireTableRowsInserted(_featureSets.size() - 1, _featureSets.size() - 1);
            }

            synchronized public void removeFeatureSet(int i)
            {
                _featureSets.remove(i);
                ArrayList changedFeatureSets = new ArrayList();
                changedFeatureSets.addAll(_featureSets);
                updateFeatureSets(changedFeatureSets, false);
                super.fireTableRowsDeleted(i, i);
            }

            public int getRowCount()
            {
                return _featureSets.size();
            }

            public int getColumnCount()
            {
                return 4;
            }

            public Class getColumnClass(int c)
            {
                switch (c)
                {
                    case 0:
                        return Boolean.class;
                    case 2:
                        return Color.class;
                    default:
                        return Object.class;
                }
            }

            public boolean isCellEditable(int row, int col)
            {
                return col == 0 || col == 2 || col == 3;
            }

            public Object getValueAt(int row, int col)
            {
                if (row >= _featureSets.size())
                    return null;

                FeatureSet fs = (FeatureSet) _featureSets.get(row);
                switch (col)
                {
                    case 0:
                        return new Boolean(fs.isDisplayed());
                    case 1:
                        return fs.getSourceFile() != null ? fs.getSourceFile().getAbsolutePath() : "(no file)";
                    case 2:
                        return fs.getColor();
                    default:
                        return fs;
                }
            }

            private String[] columnNames = new String[]{TextProvider.getText("DISPLAY"),
                    TextProvider.getText("FILE"),
                    TextProvider.getText("COLOR")};

            public String getColumnName(int i)
            {
                return i >= 0 && i < columnNames.length ? columnNames[i] : "";
            }

            public void setValueAt(Object o, int row, int col)
            {
                if (row >= _featureSets.size())
                    return;

                FeatureSet fs = (FeatureSet) _featureSets.get(row);
                switch (col)
                {
                    case 0:
                        fs.setDisplayed(((Boolean) o).booleanValue());
                        break;
                    case 1:
                        throw new IllegalArgumentException(TextProvider.getText("CANT_EDIT_FILE_NAME"));
                    case 2:
                        fs.setColor((Color) o);
                        break;
                    default:
                        return;
                }
                ArrayList changedFeatureSets = new ArrayList();
                changedFeatureSets.addAll(_featureSets);
                updateFeatureSets(changedFeatureSets, true);

                super.fireTableCellUpdated(row, col);
            }
        }


        private class OptimizeGrapher implements Runnable
        {
            FeatureSelectionFrame.FeatureSelectionDialog.Optimizer optimizer;

            private OptimizeGrapher(Optimizer optimizer)
            {
                this.optimizer = optimizer;
            }

            public void run()
            {
                textPepArrayScanBucket.setText(String.valueOf(optimizer.bestScanWindow));
                textPepArrayMzBucket.setText(String.valueOf(optimizer.bestMzWindow));
                lblNumBuckets.setText(String.valueOf(optimizer.numBuckets[optimizer.bestMatchMzIndex][optimizer.bestMatchScanIndex]));
                lblExactMatches.setText(String.valueOf(optimizer.perfectMatches[optimizer.bestMatchMzIndex][optimizer.bestMatchScanIndex]));
                XYSeriesCollection seriesColl = new XYSeriesCollection();
                for (int i = 0; i < optimizer.mzValues.length; i++)
                    seriesColl.addSeries(optimizer.series[i]);

                StringBuffer mzBucketSizeStr = new StringBuffer();
                String sep = "";
                for (int i = 0; i < optimizer.mzValues.length; i++)
                {
                    mzBucketSizeStr.append(sep);
                    mzBucketSizeStr.append(optimizer.mzValues[i]);
                    sep = ", ";
                }

                JFrame newFrame = new JFrame(TextProvider.getText("BUCKETS_FOR_MZ_OF_SIZE", new String(mzBucketSizeStr)));
                ChartPanel chart = SpectrumChartFactory.CreateChartPanel(seriesColl);
                chart.getChart().getXYPlot().getDomainAxis().setLabel(TextProvider.getText("SCAN_WINDOW"));
                chart.getChart().getXYPlot().getRangeAxis().setLabel(TextProvider.getText("QUOTE_PERFECT_BUCKETS"));
                newFrame.setContentPane(chart);
                newFrame.setSize(400, 500);
                newFrame.setVisible(true);
            }
        }

        public static class ColorRenderer extends JLabel
                implements TableCellRenderer
        {
            Border unselectedBorder = null;
            Border selectedBorder = null;
            boolean isBordered = true;

            public ColorRenderer(boolean isBordered)
            {
                this.isBordered = isBordered;
                setOpaque(true); //MUST do this for background to show up.
            }

            public Component getTableCellRendererComponent(JTable table, Object color,
                                                           boolean isSelected, boolean hasFocus,
                                                           int row, int column)
            {
                Color newColor = (Color) color;
                setBackground(newColor);
                if (isBordered)
                {
                    if (isSelected)
                    {
                        if (selectedBorder == null)
                        {
                            selectedBorder = BorderFactory.createMatteBorder(2, 5, 2, 5,
                                    table.getSelectionBackground());
                        }
                        setBorder(selectedBorder);
                    }
                    else
                    {
                        if (unselectedBorder == null)
                        {
                            unselectedBorder = BorderFactory.createMatteBorder(2, 5, 2, 5,
                                    table.getBackground());
                        }
                        setBorder(unselectedBorder);
                    }
                }

                setToolTipText(TextProvider.getText("RGB_VALUE_RED_GREEN_BLUE",
                        "RED", Integer.toString(newColor.getRed()),
                        "GREEN", Integer.toString(newColor.getGreen()),
                        "BLUE", Integer.toString(newColor.getBlue())));
                return this;
            }
        }

        public static class ColorEditor extends AbstractCellEditor
                implements TableCellEditor,
                ActionListener
        {
            Color currentColor;
            JButton button;
            JColorChooser colorChooser;
            JDialog dialog;
            protected static final String EDIT = "edit";

            public ColorEditor()
            {
                //Set up the editor (from the table's point of view),
                //which is a button.
                //This button brings up the color chooser dialog,
                //which is the editor from the user's point of view.
                button = new JButton();
                button.setActionCommand(EDIT);
                button.addActionListener(this);
                button.setBorderPainted(false);

                //Set up the dialog that the button brings up.
                colorChooser = new JColorChooser();
                dialog = JColorChooser.createDialog(button,
                        "Pick a Color",
                        true, //modal
                        colorChooser,
                        this, //OK button handler
                        null); //no CANCEL button handler
            }

            /**
             * Handles events from the editor button and from
             * the dialog's OK button.
             */
            public void actionPerformed(ActionEvent e)
            {
                if (EDIT.equals(e.getActionCommand()))
                {
                    //The user has clicked the cell, so
                    //bring up the dialog.
                    button.setBackground(currentColor);
                    colorChooser.setColor(currentColor);
                    dialog.setVisible(true);

                    //Make the renderer reappear.
                    fireEditingStopped();

                }
                else
                { //User pressed dialog's "OK" button.
                    currentColor = colorChooser.getColor();
                }
            }

            //Implement the one CellEditor method that AbstractCellEditor doesn't.
            public Object getCellEditorValue()
            {
                return currentColor;
            }

            //Implement the one method defined by TableCellEditor.
            public Component getTableCellEditorComponent(JTable table,
                                                         Object value,
                                                         boolean isSelected,
                                                         int row,
                                                         int column)
            {
                currentColor = (Color) value;
                return button;
            }
        }


        public static class BrowseEditor extends AbstractCellEditor
                implements TableCellEditor,
                ActionListener
        {
            FeatureSet fs;
            JButton button;

            public BrowseEditor()
            {
                button = new JButton("...");
                button.addActionListener(this);
            }

            public void actionPerformed(ActionEvent e)
            {
                if (null == fs)
                    return;

                FeatureFrame frame = new FeatureFrame(fs);
                frame.setVisible(true);
            }

            //Implement the one CellEditor method that AbstractCellEditor doesn't.
            public Object getCellEditorValue()
            {
                return fs;
            }

            //Implement the one method defined by TableCellEditor.
            public Component getTableCellEditorComponent(JTable table,
                                                         Object value,
                                                         boolean isSelected,
                                                         int row,
                                                         int column)
            {
                fs = (FeatureSet) value;
                return button;
            }
        }

        public static class FeatureFrame extends JFrame
        {
            FeatureSet _fs;
            FeatureTable _table;

            JDialog _popup = null;

            public FeatureFrame(FeatureSet fs)
            {
                super((fs.getSourceFile() != null ? fs.getSourceFile().getName() : "(no file)") + " -- " + fs.getFeatures().length);
                setSize(640, 400);
                setIconImage(ApplicationContext.getFrame().getIconImage());
                JPanel editPanel = new JPanel();
                JButton propertiesButton = new JButton(TextProvider.getText("PROPERTIES"));
                JButton deleteButton = new JButton(TextProvider.getText("DELETE"));
                JButton saveButton = new JButton(TextProvider.getText("SAVE"));
                JButton addButton = new JButton(TextProvider.getText("ADD"));
                JButton closeButton = new JButton(TextProvider.getText("CLOSE"));

                ListenerHelper helper = new ListenerHelper(this);
                helper.addListener(propertiesButton, "properties_actionPerformed");
                helper.addListener(deleteButton, "delete_actionPerformed");
                helper.addListener(addButton, "add_actionPerformed");
                helper.addListener(saveButton, "save_actionPerformed");
                helper.addListener(closeButton, "close_actionPerformed");
                editPanel.add(propertiesButton);
                editPanel.add(addButton);
                editPanel.add(deleteButton);
                editPanel.add(saveButton);
                editPanel.add(closeButton);

                _table = new FeatureTable(fs);
                JScrollPane pane = new JScrollPane(_table);
                getContentPane().add(pane);
                getContentPane().add(editPanel, BorderLayout.PAGE_END);
            }


            public void properties_actionPerformed(ActionEvent e)
            {
                if (true)
                {
                    if (ApplicationContext.getFrame() instanceof WorkbenchFrame)
                        ((WorkbenchFrame) ApplicationContext.getFrame()).showPropertiesPane();
                    ApplicationContext.setProperty(SharedProperties.SELECTED, _fs);
                }
                else
                {
                    PropertiesPane.popup(_fs.getProperties());
                }
/*            Map map = _fs.getProperties();
		    if (null == map)
			    map = new HashMap();

		    TableModel m = new DefaultTableModel(map.size(), 2);
		    int row = 0;
		    for (Iterator it = map.entrySet().iterator(); it.hasNext();)
			    {
			    Map.Entry entry = (Map.Entry)it.next();
				m.setValueAt(entry.getKey(), row, 0);
			    m.setValueAt(entry.getValue(), row, 1);
			    row++;
			    }

		    JTable table = new JTable(m);
		    table.getColumnModel().getColumn(0).setHeaderValue("property");
		    table.getColumnModel().getColumn(1).setHeaderValue("value");
		    String title = "" + _fs.getSourceFile() + " properties";
		    JDialog dialog = new JDialog(ApplicationContext.getFrame(), title, true);

		    ListenerHelper h = new ListenerHelper(this);
		    h.addListener(dialog.getContentPane(), "popup_keyPressed");
		    h.addListener(table, "popup_keyPressed");
		    dialog.getContentPane().add(new JScrollPane(table));
		    dialog.setSize(400, 300);
		    _popup = dialog;
		    dialog.setVisible(true);
				}
*/            }

/*	    public void popup_keyPressed(KeyEvent e)
            {
            if (e.getKeyCode() == KeyEvent.VK_ESCAPE && null != _popup)
                {
                _popup.dispose();
                _popup = null;
                e.consume();
                return;
                }
            } */


            public void delete_actionPerformed(ActionEvent e)
            {
                int row;
                int sortRow = _table.getSelectedRow();
                if (sortRow < 0 || sortRow >= _fs.getFeatures().length)
                    return;
                row = _table.convertRowIndexToModel(sortRow);

                Feature[] features = _fs.getFeatures();
                Feature[] newFeatures = new Feature[features.length - 1];
                for (int i = 0; i < row; i++)
                    newFeatures[i] = features[i];
                for (int i = row + 1; i < features.length; i++)
                    newFeatures[i - 1] = features[i];

                //TODO: Move fs into tablemodel so don't have to hack around like this...
                _fs.setFeatures(newFeatures);
                _table._model._features = newFeatures;
                _table._model.fireTableRowsDeleted(row, row);
                _table.getSelectionModel().setSelectionInterval(sortRow, sortRow);
                _table.requestFocus();
                FeatureSelectionDialog.getInstance().updateFeatureSets(FeatureSelectionDialog.getInstance().getFeatureSets(), false);

            }


            public void add_actionPerformed(ActionEvent e)
            {
                Spectrum.Peak p = (Spectrum.Peak) ApplicationContext.getProperty(SharedProperties.SELECTED_POINT);
                Feature f;
                if (null == p)
                    f = new Feature();
                else
                    f = new Feature(p);

                int row;
                int sortRow = _table.getSelectedRow();
                sortRow = Math.min(sortRow, _fs.getFeatures().length - 1);
                row = _table.convertRowIndexToModel(sortRow);

                Feature[] features = _fs.getFeatures();
                Feature[] newFeatures = new Feature[features.length + 1];
                for (int i = 0; i <= row; i++)
                    newFeatures[i] = features[i];
                newFeatures[row + 1] = f;
                for (int i = row + 1; i < features.length; i++)
                    newFeatures[i + 1] = features[i];

                String userName = System.getProperty("user.name");
                //not translating, since this will be stored at the feature level
                f.setDescription("Added=" + userName + ", " + new SimpleDateFormat().format(new Date()));
                //TODO: Move fs into tablemodel so don't have to hack around like this...
                _fs.setFeatures(newFeatures);
                _table._model._features = newFeatures;
                _table._model.fireTableRowsInserted(row + 1, row + 1);
                _table.getSelectionModel().setSelectionInterval(sortRow + 1, sortRow + 1);
                _table.requestFocus();
                FeatureSelectionDialog.getInstance().updateFeatureSets(FeatureSelectionDialog.getInstance().getFeatureSets(), false);
            }


            public void close_actionPerformed(ActionEvent e)
            {
                this.dispose();
            }


            public void save_actionPerformed(ActionEvent e)
            {
                try
                {
                    if (_fs.getSourceFile() == null)
                    {
                        WorkbenchFileChooser wfc = new WorkbenchFileChooser();
                        wfc.showSaveDialog(this);
                        File file = wfc.getSelectedFile();
                        _fs.setSourceFile(file);
                    }
                    _fs.save();
                }
                catch (Exception x)
                {
                    //JOptionPane.showMessageDialog(null, "Exception saving: " + x.getMessage());
                    ApplicationContext.errorMessage(null, x);
                }
            }


            private void selectRow()
            {
                if (null == _fs)
                    return;
                if (null == _table)
                    return;
                int row = _table.getSelectedRow();
                if (row < 0 || row >= _fs.getFeatures().length)
                    return;
                row = _table.convertRowIndexToModel(row);

                Feature f = _fs.getFeatures()[row];
                //System.err.println(f);
                float mz = f.getMz();
                if (mz <= 0)
                    return;
                // UNDONE: float time = f.getCPUTime();
                int scanNum = f.getScan();
                if (scanNum <= 0)
                    return;
                MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
                MSRun.MSScan selectedScan = null;
                if (run != null)
                {
                    if (scanNum >= run.getScanCount())
                        scanNum = run.getScanCount() - 1;
                    selectedScan = run.getScan(scanNum);
                }
                ApplicationContext.setProperty(SharedProperties.MS_SCAN, selectedScan);
                Feature c = (Feature) f.clone();
                ApplicationContext.setProperty(SharedProperties.SELECTED_POINT, c);
                ApplicationContext.setProperty(SharedProperties.SELECTED, c);
            }


            public void table_keyPressed(KeyEvent e)
            {
                ApplicationContext.setMessage("" + e.getKeyCode());
                if (e.getKeyCode() != 10)
                    return;
                e.consume();
                selectRow();
            }


            public void table_mouseClicked(MouseEvent e)
            {
                if (e.getClickCount() != 2)
                    return;
                e.consume();
                selectRow();
            }


            public class FeatureTable extends JTable
            {
                FeatureTableModel _model;
                JTextField _textEditor;

                public FeatureTable(FeatureSet fs)
                {
                    _fs = fs;
                    _model = new FeatureTableModel(_fs.getFeatures());
                    setModel(_model);

                    //dhmay replacing the old usage of our custom TableSorter class with the standard 1.6 syntax
                    TableRowSorter<TableModel> sorter
                            = new TableRowSorter<TableModel>(_model);
                    setRowSorter(sorter);
                    int count = _model.getColumnCount();

                    for (int i = 0; i < count - 1; i++)
                        getColumnModel().getColumn(i).setCellRenderer(new FeatureSelectionDialog.NumberRenderer());

                    getColumnModel().getColumn(count - 1).setPreferredWidth(200);
                    ListenerHelper helper = new ListenerHelper(FeatureFrame.this);
                    helper.addListener(this, "table_mouseClicked");
                    helper.addListener(this, "table_keyPressed");
                }

                public Component prepareEditor(TableCellEditor editor, int row, int column)
                {
                    Component c = super.prepareEditor(editor, row, column);

                    if (c instanceof JTextField)
                        ((JTextField) c).selectAll();

                    return c;
                }

                public void valueChanged(ListSelectionEvent e)
                {
                    super.valueChanged(e);
                    selectRow();
                }
            }
        }

        public static class NumberRenderer extends DefaultTableCellRenderer
        {
            NumberFormat format;
            int decimal = 2;

            public NumberRenderer()
            {
                super();
                setHorizontalAlignment(RIGHT);
                format = new DecimalFormat();
                format.setMaximumFractionDigits(decimal);
                format.setMinimumFractionDigits(decimal);
            }

            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
            {
                Object formatted = value;
                if (value instanceof Float || value instanceof Double)
                    formatted = format.format(value);
                return super.getTableCellRendererComponent(table, formatted, isSelected, hasFocus, row, column);
            }
        }

    }


}
