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
package org.fhcrc.cpl.viewer.align.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MoleculeRenderer2D;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAdd2HMod;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAddWaterMod;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.viewer.gui.FeatureViewerFrame;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.ms2.gui.MS2ScanViewer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.exception.CDKException;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.List;

/**
 */
public class ShowArrayRowMS2ScansCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ShowArrayRowMS2ScansCLM.class);

    protected PeptideArrayAnalyzer arrayAnalyzer;
    protected List<Integer> rowIds;

    protected int width = 400;
    protected int height = 400;

    protected float deltaPPM = 5;
    protected int deltaScans = 10;

    protected File mzXmlDir;

    protected Map<String, MSRun.MSScan[]> runNameMS2ScansMap = new HashMap<String, MSRun.MSScan[]>();


    JLabel scanInfoLabel = null;

    MS2ScanViewer.MultiMS2ScanViewer multiMS2ScanViewer;

    protected int numFragmentsPerScan = 6;

    protected float maxConsensusFragmentMzDiameter = 0.5f;

    protected boolean showCharts = true;

    protected File outFile;

    PrintWriter outPW = null;

    public ShowArrayRowMS2ScansCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "showarrayrowms2scans";
        mShortDescription = "showarrayrowms2scans";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Peptide array"),
                        new StringListArgumentDefinition("ids", true, "Row IDs to display, separated by commas"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "mzXML directory"),
                        new FileToReadArgumentDefinition("detailsfile", true, "Peptide array details file"),
                        new DecimalArgumentDefinition("deltappm", false, "Delta PPM", deltaPPM),
                        new DecimalArgumentDefinition("deltascans", false, "Delta seconds", deltaScans),
                        new IntegerArgumentDefinition("fragsperscan", false, "Fragments to consider per scan", numFragmentsPerScan),
                        new DecimalArgumentDefinition("maxfragbinsize", false, "Maximum fragment bin size", maxConsensusFragmentMzDiameter),
                        new BooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                        new FileToWriteArgumentDefinition("out",false,"output file with scan counts and peak lists"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        showCharts = getBooleanArgumentValue("showcharts");
        File pepArrayFile = this.getUnnamedFileArgumentValue();

        try
        {
            ApplicationContext.infoMessage("Loading array...");
            arrayAnalyzer = new PeptideArrayAnalyzer(pepArrayFile);
            if (hasArgumentValue("detailsfile"))
                arrayAnalyzer.setDetailsFile(getFileArgumentValue("detailsfile"));
            ApplicationContext.infoMessage("Loading details...");
            arrayAnalyzer.loadDetailsMap();
            ApplicationContext.infoMessage("Done loading array");
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to open array file",e);
        }

        mzXmlDir = getFileArgumentValue("mzxmldir");

        List<String> rowIdStrings = getStringListArgumentValue("ids");
        rowIds = new ArrayList<Integer>();
        for (String idString : rowIdStrings)
            rowIds.add(Integer.parseInt(idString));

        numFragmentsPerScan = getIntegerArgumentValue("fragsperscan");
        maxConsensusFragmentMzDiameter = getFloatArgumentValue("maxfragbinsize");
        outFile = getFileArgumentValue("out");

    }

    public void execute() throws CommandLineModuleExecutionException
    {
        if (outFile != null) {
            try {
            outPW = new PrintWriter(outFile);
                outPW.println("row\tms2scancount\tpeak1\tpeak2\tpeak3\tpeak4\tint1\tint2\tint3\tint4");
            } catch (Exception e) { throw new CommandLineModuleExecutionException("Can't write out file",e); }
        }
        for (int rowId : rowIds)  {
            ApplicationContext.infoMessage("\n\nProcessing row " + rowId);
            displayRow(rowId, arrayAnalyzer.getRowMapWithID(rowId));
        }
        if (outPW != null) {
            outPW.close();
        }
    }

    protected MSRun.MSScan[] displayMs2ScansNearFeature(Feature feature, String fileString, List<MSRun.MSScan> scansList)
             throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Examining feature " + feature + " ...");
        MSRun.MSScan[] ms2Scans = runNameMS2ScansMap.get(fileString);

        if (ms2Scans == null)
        {
            try
            {
                ApplicationContext.infoMessage("Finding scans for run " + fileString + "...");
//                File mzXmlFile = CommandLineModuleUtilities.findFileLikeFile(new File(fileString), mzXmlDir, ".mzXML");
//                ApplicationContext.infoMessage("Loading run from file " + mzXmlFile.getAbsolutePath());
//                ms2Scans = MSRun.load(mzXmlFile.getAbsolutePath()).getMS2Scans();
                List<MSRun.MSScan> ms2ScanList = new ArrayList<MSRun.MSScan>();
                boolean foundAFile = false;
                for (File potentialFile : mzXmlDir.listFiles())
                {
                    String potentialFileName = potentialFile.getName();
                    if (!potentialFileName.endsWith(".mzXML"))
                        continue;
                    if (potentialFileName.contains(fileString))
                    {
                        ApplicationContext.infoMessage("\tadding scans from file " + potentialFileName + "...");
                        ms2ScanList.addAll(Arrays.asList(MSRun.load(potentialFile.getAbsolutePath()).getMS2Scans()));
                        foundAFile = true;
                        break;
                    }
                }
                if (!foundAFile)
                    throw new CommandLineModuleExecutionException("No mzXML file found for run " + fileString);
                ms2Scans = ms2ScanList.toArray(new MSRun.MSScan[0]);
                runNameMS2ScansMap.put(fileString, ms2Scans);

            }
            catch (IOException e)
            {
                throw new  CommandLineModuleExecutionException("Failed to load run from file " + fileString,e);
            }
        }

        Arrays.sort(ms2Scans, new Comparator<MSRun.MSScan>() {

            int _compareAsc(float a, float b)
            {
                return a == b ? 0 : a < b ? -1 : 1;
            }
            public int compare(MSRun.MSScan o1, MSRun.MSScan o2)
            {
                return _compareAsc(o1.getPrecursorMz(), o2.getPrecursorMz());
            }

        });

        for (MSRun.MSScan ms2Scan : ms2Scans)
        {
            float precursorMz = ms2Scan.getPrecursorMz();
            float deltaMz =  precursorMz - feature.getMz();
            float deltaMass = deltaMz * feature.getCharge();
            float deltaMassPPM = (deltaMass) * 1000000 / feature.getMass();

            if (deltaMassPPM > deltaPPM)
                break;

            if (Math.abs(deltaMassPPM) <= deltaPPM)
            {
                if (ms2Scan.getNum() >= feature.getScanFirst()-deltaScans && ms2Scan.getNum() <= feature.getScanLast()+deltaScans) {
                    ApplicationContext.infoMessage("Adding scan " + ms2Scan.getNum() + " with precursor m/z " +
                            ms2Scan.getPrecursorMz()  + ", deltaMassPPM=" + deltaMassPPM);
                    scansList.add(ms2Scan);
                }
            }
        }

        return ms2Scans;
    }





    protected void displayRow(int rowId, Map<String, Object> rowMap)  throws CommandLineModuleExecutionException
    {
        List<MSRun.MSScan> scansList = new ArrayList<MSRun.MSScan>();

        Map<String, List<Feature>> fileFeaturesMapThisRow = arrayAnalyzer.createFileFeatureMapForRow(rowId);
        for (String fileString : fileFeaturesMapThisRow.keySet())
        {
            ApplicationContext.infoMessage("Displaying features from run " + fileString);
            for (Feature feature : fileFeaturesMapThisRow.get(fileString))
                displayMs2ScansNearFeature(feature, fileString, scansList);
        }
        if (scansList.isEmpty())
        {
            if (outPW != null)  {
                outPW.println(rowId + "\t0\t\t\t\t\t\t\t\t");
                outPW.flush();
            }
            ApplicationContext.infoMessage("No MS/MS scans.");
            return;
        }

        List<FragmentIonCluster> ionClusters = new ArrayList<FragmentIonCluster>();
        for (MSRun.MSScan scan : scansList) {
            float[] mzValues = scan.getSpectrum()[0];
            float[] intensityValues = scan.getSpectrum()[1];

            List<Pair<Float, Float>> peaksAsPairs = new ArrayList<Pair<Float, Float>>();
            for (int i=0; i<mzValues.length; i++)
                peaksAsPairs.add(new Pair<Float, Float>(mzValues[i], intensityValues[i]));
            Collections.sort(peaksAsPairs, new Comparator<Pair<Float, Float>>() {
                @Override
                public int compare(Pair<Float, Float> pair1, Pair<Float, Float> pair2) {
                    float diff = pair1.second - pair2.second;
                    if (diff > 0) return -1;
                    if (diff < 0) return 1;
                    return 0;
                }
            });

            for (int i=0; i<numFragmentsPerScan; i++) {
                Pair<Float, Float> pair = peaksAsPairs.get(i);
                ionClusters.add(new FragmentIonCluster(pair.first,pair.second));
            }
        }

        Collections.sort(ionClusters, new Comparator<FragmentIonCluster>() {
            int _compareAsc(float a, float b)
            {
                return a == b ? 0 : a < b ? -1 : 1;
            }
            public int compare(FragmentIonCluster o1, FragmentIonCluster o2)
            {
                return _compareAsc(o1.medianMz, o2.medianMz);
            }
        });

        List<Float> binDistances = new ArrayList<Float>();
        for (int i=0; i<ionClusters.size()-1; i++)
            binDistances.add(ionClusters.get(i+1).medianMz - ionClusters.get(i).medianMz);
        float minBinDistance = (float) BasicStatistics.min(binDistances);
        while (minBinDistance < maxConsensusFragmentMzDiameter) {
            for (int i=0; i<binDistances.size(); i++)
                if (binDistances.get(i) == minBinDistance)
                {
                    binDistances.remove(i);
                    ionClusters.get(i).agglomerateOtherCluster(ionClusters.get(i+1));
                    ionClusters.remove(i+1);
                    if (i<ionClusters.size()-1)
                        binDistances.set(i, ionClusters.get(i+1).medianMz - ionClusters.get(i).medianMz);
                    minBinDistance = (float) BasicStatistics.min(binDistances);
                    break;
                }
        }

        Collections.sort(ionClusters, new Comparator<FragmentIonCluster>() {
            int _compareDesc(float a, float b)
            {
                return a == b ? 0 : a > b ? -1 : 1;
            }
            public int compare(FragmentIonCluster o1, FragmentIonCluster o2)
            {
                return _compareDesc(o1.logIntensity, o2.logIntensity);
            }
        });

        String ion0String = "";
        String ion1String = "";
        String ion2String = "";
        String ion3String = "";

                String ion0intString = "";
        String ion1intString = "";
        String ion2intString = "";
        String ion3intString = "";

        ApplicationContext.infoMessage("Top fragment ions:");
        for (int i=0; i<Math.min(10,ionClusters.size()); i++) {
            ApplicationContext.infoMessage("\tmz=" + ionClusters.get(i).medianMz + ", logInt=" + ionClusters.get(i).logIntensity);
            if (i==0) ion0String = "" + ionClusters.get(i).medianMz;
            if (i==1) ion1String = "" + ionClusters.get(i).medianMz;
            if (i==2) ion2String = "" + ionClusters.get(i).medianMz;
            if (i==3) ion3String = "" + ionClusters.get(i).medianMz;

            if (i==0) ion0intString = "" + ionClusters.get(i).logIntensity;
            if (i==1) ion1intString = "" + ionClusters.get(i).logIntensity;
            if (i==2) ion2intString = "" + ionClusters.get(i).logIntensity;
            if (i==3) ion3intString = "" + ionClusters.get(i).logIntensity;

        }
        if (outPW != null)  {
            outPW.println(rowId + "\t" + scansList.size() + "\t" + ion0String  + "\t" + ion1String+ "\t" + ion2String+ "\t" + ion3String
                    + "\t" + ion0intString  + "\t" + ion1intString+ "\t" + ion2intString+ "\t" + ion3intString);
            outPW.flush();
        }


        List<String> scanNumbersList = new ArrayList<String>();
        for (MSRun.MSScan scan : scansList)
            scanNumbersList.add("" + scan.getNum());

        ApplicationContext.infoMessage("Row " + rowId + ": " + scansList.size() + "; bigions=" + ion0String + "," + ion1String + "," + ion2String + "," + ion3String  + "; scans=" +
                MS2ExtraInfoDef.convertStringListToString(scanNumbersList));

        if (showCharts) {
            scanInfoLabel = new JLabel("Scan , Precursor m/z: ");
            scanInfoLabel.setVisible(true);
            multiMS2ScanViewer =
                    new MS2ScanViewer.MultiMS2ScanViewer(scansList.toArray(new MSRun.MSScan[0]), 3);
            ChangeListener scanChangeListener = new MS2ScanChangedListener();
            multiMS2ScanViewer.addChangeListener(scanChangeListener);
            JDialog dialog = new JDialog();
            GridBagConstraints fullRowGBC = new GridBagConstraints();
            fullRowGBC.gridwidth = GridBagConstraints.REMAINDER;
            dialog.setLayout(new GridBagLayout());

            dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            dialog.setTitle("MS2 Scans");
            dialog.setSize(new Dimension(800, 600));
            dialog.add(scanInfoLabel, fullRowGBC);
            dialog.add(multiMS2ScanViewer, fullRowGBC);
            dialog.setVisible(true);
        }


    }

    protected class MS2ScanChangedListener implements ChangeListener
    {
        public void stateChanged(ChangeEvent e)
        {
            MSRun.MSScan scan = multiMS2ScanViewer.getMs2ScanViewer().getScanInViewer();
            if (showCharts) {
                scanInfoLabel.setText("Scan " + scan.getNum() + ", Precursor m/z: " + scan.getPrecursorMz());
                scanInfoLabel.updateUI();
            }
        }
    }

    protected class FragmentIonCluster
    {
        public List<Float> mzs;
        public float medianMz;
        public float logIntensity;

        public FragmentIonCluster(float mz, float intensity) {
            this.mzs = new ArrayList<Float>();
            mzs.add(mz);
            medianMz = mz;
            logIntensity = (float) Math.log(intensity);
        }

        public void agglomerateOtherCluster(FragmentIonCluster otherCluster) {
            mzs.addAll(otherCluster.mzs);
            medianMz = (float) BasicStatistics.median(mzs);
            logIntensity += otherCluster.logIntensity;
        }
    }
}


