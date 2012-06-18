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
package org.fhcrc.cpl.viewer.qa.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.modules.FilterFeaturesCommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.PlotMassCalibrationCLM;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.FeatureExtractor;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFindingBroker;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.viewer.gui.MSImageComponent;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 * test
 */
public class QAExperimentCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(QAExperimentCLM.class);

    protected File mzXmlDir;
    protected File allPepXmlFile;
    protected File allProtXmlFile;
    protected File protGeneFile;

    protected File qaDir;

    protected boolean force = false;

    protected float filteredMaxKL = 1.0f;
    protected int filteredMinPeaks = 3;
    protected float filteredMinMz = 400;
    protected float filteredMaxMz = 2500;

    protected float minPeptideProphet = .75f;
    protected float minProteinProphet = .9f;

    protected float labelMassDiff = 3f;


    protected List<Map<String,String>> runSummaryDataMaps = new ArrayList<Map<String,String>>();
    protected List<Set<String>> peptidesInEachRun = new ArrayList<Set<String>>();
    protected List<Set<String>> quantPeptidesInEachRun = new ArrayList<Set<String>>();

    protected List<Set<String>> proteinsInEachRun = new ArrayList<Set<String>>();
    protected List<Set<String>> genesInEachRun = new ArrayList<Set<String>>();
    protected List<Set<String>> quantProteinsInEachRun = new ArrayList<Set<String>>();
    protected List<Set<String>> quantGenesInEachRun = new ArrayList<Set<String>>();

    protected int numTotalPeptideIdentifications = 0;
    protected int numTotalQuantitatedPeptideIdentifications = 0;

    protected boolean shouldAnalyzeMS1 = true;

    public QAExperimentCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "qaexperiment";

        mHelpMessage ="Perform QA analysis on a single experiment";
        mShortDescription = "Perform QA analysis on a single experiment";

        CommandLineArgumentDefinition[] argDefs =
               {
                       new DirectoryToReadArgumentDefinition("mzxmldir", false, "mzXML Directory"),
                       new FileToReadArgumentDefinition("allpepxml", true, "all.pep.xml filepath"),
                       new FileToReadArgumentDefinition("allprotxml", true, "all.prot.xml filepath"),
                       new FileToReadArgumentDefinition("protgenefile", true,
                               "File associating gene symbols with protein accession numbers"),                       
                       new DirectoryToReadArgumentDefinition("qadir", true, "QA Output Root Directory"),
                       new DecimalArgumentDefinition("minpeptideprophet", false,
                               "Minimum PeptideProphet probability", minPeptideProphet),
                       new DecimalArgumentDefinition("minproteinprophet", false,
                               "Minimum ProteinProphet probability", minProteinProphet),
                       new DecimalArgumentDefinition("labelmassdiff", false,
                               "Isotopic label mass difference", labelMassDiff),
                       new BooleanArgumentDefinition("force", false,
                               "Force re-creation of output files if they exist?", force),
                       new BooleanArgumentDefinition("noms1", false,
                               "No MS1 analysis -- only pepXML and protXML", !shouldAnalyzeMS1)
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzXmlDir = getFileArgumentValue("mzxmldir");
        allPepXmlFile = getFileArgumentValue("allpepxml");
        allProtXmlFile = getFileArgumentValue("allprotxml");
        protGeneFile = getFileArgumentValue("protgenefile");

        minPeptideProphet =  getFloatArgumentValue("minpeptideprophet");
        minProteinProphet =  getFloatArgumentValue("minproteinprophet");

        labelMassDiff = getFloatArgumentValue("labelmassdiff");

        shouldAnalyzeMS1 = !getBooleanArgumentValue("noms1");

        if (shouldAnalyzeMS1)
            assertArgumentPresent("mzxmldir", "noms1");

        qaDir = getFileArgumentValue("qadir");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Starting MS1 QA Analysis...");
        List<String> runNames = null;

        if (shouldAnalyzeMS1)
        {
        try
        {
            ApplicationContext.infoMessage("Starting MS1 QA Analysis...");
            runNames = ms1QA();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        }

        ApplicationContext.infoMessage("Starting MS2 QA Analysis...");
        ms2QA();

    }

    protected List<String> ms1QA()
            throws Exception
    {
        File featuresDir = new File(qaDir,"features");
        if (!featuresDir.exists())
            featuresDir.mkdir();

        File filteredFeaturesDir = new File(featuresDir, "filtered");
        if (!filteredFeaturesDir.exists())
            filteredFeaturesDir.mkdir();        

        List<File> featureFiles = new ArrayList<File>();

        List<Double> featureCounts = new ArrayList<Double>();
        List<Double> filteredFeatureCounts = new ArrayList<Double>();
        List<String> runNames = new ArrayList<String>();

        File outScanCountsChartFile = new File(qaDir, "scancounts.png");
        boolean shouldWriteScanCounts = (!outScanCountsChartFile.exists() || force);
        List<Double> ms1ScanCounts = new ArrayList<Double>();
        List<Double> ms2ScanCounts = new ArrayList<Double>();

        for (File mzXmlFile : mzXmlDir.listFiles())
        {
            if (mzXmlFile.isDirectory())
                continue;

            String mzXmlFileName = mzXmlFile.getName();

            //dhmay removing a separate check for .xml extension, because that caused trouble
            if (!mzXmlFileName.toLowerCase().endsWith(".mzxml"))
                continue;

            int mzXmlFileNameLength = mzXmlFileName.length();
            ApplicationContext.setMessage("Processing file " + mzXmlFileName);

            String runName = mzXmlFileName;
            if (mzXmlFileName.toLowerCase().endsWith(".mzxml"))
            {
                runName =mzXmlFileName.substring(0, mzXmlFileNameLength -
                                ".mzxml".length());
            }
            else if (mzXmlFileName.endsWith(".xml"))
            {
                runName = mzXmlFileName.substring(0, mzXmlFileNameLength -
                                ".xml".length());
            }

            runNames.add(runName);
            String featuresFileName = runName + ".peptides.tsv";

            File outFeaturesFile = new File(featuresDir, featuresFileName);
            featureFiles.add(outFeaturesFile);

            MSRun run = null;
            FeatureSet featureSet = null;
            if (outFeaturesFile.exists() && !force)
            {
                ApplicationContext.setMessage("Feature file " + outFeaturesFile.getAbsolutePath() +
                        " exists, not overwriting.");
                featureSet = new FeatureSet(outFeaturesFile);
            }
            else
            {
                run = MSRun.load(mzXmlFile.getAbsolutePath());
                featureSet = FeatureFindingBroker.findPeptides(
                        run, 1, run.getScanCount(),
                        PeakCombiner.DEFAULT_MAX_CHARGE,
                        FeatureExtractor.getMzExtractionRange(run),
                        0, FeatureFinder.DEFAULT_ACCURATE_MASS_ADJUSTMENT_SCANS,
                        FeatureFinder.DEFAULT_FEATURE_FINDING_CLASS,
                        true, false, false);
                featureSet.save(outFeaturesFile);
                ApplicationContext.setMessage("Saved features file " + outFeaturesFile.getAbsolutePath());
            }

            featureCounts.add((double) featureSet.getFeatures().length);

            FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
            sel.setMinPeaks(filteredMinPeaks);
            sel.setMaxKL(filteredMaxKL);
            sel.setMinMz(filteredMinMz);
            sel.setMaxMz(filteredMaxMz);

            FeatureSet filteredFeatureSet = featureSet.filter(sel);
            filteredFeatureCounts.add((double) filteredFeatureSet.getFeatures().length);


            File outFilteredFeaturesFile = new File(filteredFeaturesDir,
                    FilterFeaturesCommandLineModule.createFilteredFeatureFileFilename(outFeaturesFile.getName(),
                            FilterFeaturesCommandLineModule.OUT_FORMAT_MSINSPECT));


            if (outFilteredFeaturesFile.exists())
            {
                ApplicationContext.setMessage("Filtered feature file " + outFilteredFeaturesFile.getAbsolutePath() +
                        " exists, not overwriting.");
            }
            else
            {
                filteredFeatureSet.save(outFilteredFeaturesFile);
                ApplicationContext.setMessage("Saved filtered features file " + outFilteredFeaturesFile.getAbsolutePath());
            }

            File outImageFile = new File(qaDir, mzXmlFileName + ".png");
            if (force || !outImageFile.exists())
            {
                if (run == null)
                    run = MSRun.load(mzXmlFile.getAbsolutePath());
                MSImageComponent comp =
                        new MSImageComponent(run.getImage(MSImageComponent.getPrefColorScheme()));
                comp.setRun(run);
            
                comp.saveImage(outImageFile, Integer.MAX_VALUE, Integer.MAX_VALUE, false);

                ApplicationContext.setMessage("Wrote run image file " + outImageFile.getAbsolutePath());
            }
            else
            {
                ApplicationContext.setMessage("Image file " + outImageFile.getAbsolutePath() +
                        " already exists, not overwriting");
            }

            if (shouldWriteScanCounts)
            {
                if (run == null)
                    run = MSRun.load(mzXmlFile.getAbsolutePath());
                ms1ScanCounts.add((double) run.getScanCount());
                ms2ScanCounts.add((double) run.getMS2Scans().length);
            }
        }


        List<Double> runNumbers = new ArrayList<Double>();
        for (int i=0; i<runNames.size(); i++)
            runNumbers.add((double)(i+1));

        if (shouldWriteScanCounts)
        {
            PanelWithLineChart lineChartScans = new PanelWithLineChart();
            lineChartScans.addData(runNumbers, ms1ScanCounts, "MS1 Scans");
            lineChartScans.addData(runNumbers, ms2ScanCounts, "MS2 Scans");
            lineChartScans.saveChartToImageFile(outScanCountsChartFile);
            ApplicationContext.setMessage("Saved chart file " + outScanCountsChartFile.getAbsolutePath());
        }

        PanelWithLineChart lineChartFeatureCounts = new PanelWithLineChart(runNumbers, featureCounts,
                                                                           "Feature Counts per Run");
        File featureCountsChartFile = new File(qaDir, "ms1_feature_counts.png");
        lineChartFeatureCounts.saveChartToImageFile(featureCountsChartFile);
        ApplicationContext.setMessage("Saved chart file " + featureCountsChartFile.getAbsolutePath());

        PanelWithLineChart lineChartFilteredFeatureCounts = new PanelWithLineChart(runNumbers, filteredFeatureCounts,
                                                                           "Feature Counts per Run");
        File filteredFeatureCountsChartFile = new File(qaDir, "ms1_filtered_feature_counts.png");
        lineChartFilteredFeatureCounts.saveChartToImageFile(filteredFeatureCountsChartFile);
        ApplicationContext.setMessage("Saved chart file " + filteredFeatureCountsChartFile.getAbsolutePath());
        
        PlotMassCalibrationCLM plotCalibModule = new PlotMassCalibrationCLM();
        File calibrationPlotFile = new File(qaDir,"mass_calibration.png");
        Map<String,String> argMap = new HashMap<String,String>();
        argMap.put("indir", filteredFeaturesDir.getAbsolutePath());
        argMap.put("outboxwhiskersplot",calibrationPlotFile.getAbsolutePath());
        argMap.put("showcharts","false");
        plotCalibModule.digestArguments(argMap);
        plotCalibModule.assignArgumentValues();
        plotCalibModule.execute();
        ApplicationContext.setMessage("Saved chart file " + calibrationPlotFile.getAbsolutePath());

        return runNames;
    }



    protected void ms2QA()
            throws CommandLineModuleExecutionException
    {
//        File outLegendFile = new File(qaDir,"legend.png");
//        try
//        {
//            String rPlotsCommand =
//                    RInterface.readResourceFile("/org/fhcrc/cpl/viewer/qa/QA_plots.R");
//            StringBuffer rCommandBuf = new StringBuffer("run_names<-c(");
//            for (int i=0; i<runNames.size(); i++)
//            {
//                rCommandBuf.append("'" + runNames.get(i) + "'");
//                if (i < runNames.size()-1)
//                    rCommandBuf.append(",");
//            }
//            rCommandBuf.append(");\n");
//            rCommandBuf.append(rPlotsCommand);
//
//            Map<String, Object> rScalarVarMap = new HashMap<String, Object>();
//            rScalarVarMap.put("out_legend_file", "'" + outLegendFile.getAbsolutePath() + "'");
//
//            RInterface.evaluateRExpression(rCommandBuf.toString(),
//                rScalarVarMap, null, null, null, 10000);
//
//        }
//        catch (IOException e)
//        {
//            throw new RuntimeException(e);
//        }

        try
        {
            File allProtTextFile = new File(allPepXmlFile.getParentFile(), "all.prot.txt");
            File allPepTextFile = new File(allPepXmlFile.getParentFile(), "all.pep.txt");
            File summaryTextFile = new File(allPepXmlFile.getParentFile(), "summary.txt");

            if (force || !allProtTextFile.exists() || !allPepTextFile.exists() || !summaryTextFile.exists() )
            {
                ApplicationContext.setMessage("Creating all.pep.txt...");
                createAllPepText(allPepTextFile);
                ApplicationContext.setMessage("Creating all.prot.txt...");
                createAllProtText(allProtTextFile);
                ApplicationContext.setMessage("Creating summary.txt...");
                createSummaryText(summaryTextFile);
            }
            else
            {
                ApplicationContext.setMessage("xtandem text files already exist.");
            }


            Map<String, File> rVarFileMap = new HashMap<String, File>();
            rVarFileMap.put("all_pep_text_file", allPepTextFile);
            rVarFileMap.put("all_prot_text_file", allProtTextFile);
            rVarFileMap.put("summary_text_file", summaryTextFile);

            rVarFileMap.put("out_legend_file", new File(qaDir,"legend.png"));
            rVarFileMap.put("delta_mass_plots_file", new File(qaDir, "delta_mass_plots.png"));
            rVarFileMap.put("delta_mass_boxplots_file", new File(qaDir, "delta_mass_boxplots.png"));
            rVarFileMap.put("quantified_protein_groups_file", new File(qaDir,"quantified_protein_groups.png"));
            rVarFileMap.put("peptides_file", new File(qaDir,"peptides.png"));
            rVarFileMap.put("gene_protein_groups_file", new File(qaDir,"gene_protein_groups.png"));
            rVarFileMap.put("quantified_peps_file", new File(qaDir,"quantified_peps.png"));
            rVarFileMap.put("quantified_peps_numlabels_file", new File(qaDir,"quantified_peps_numlabels.png"));
            rVarFileMap.put("quantified_mass_summary_file", new File(qaDir,"quantile_mass_summary.csv"));
            rVarFileMap.put("quantile_intensity_summary_file", new File(qaDir,"quantile_intensity_summary.csv"));
            rVarFileMap.put("quantile_cys_summary_file", new File(qaDir,"quantile_cys_summary.csv"));

            boolean chartMissing = false;

            if (!force)
            {
                for (File outChartFile : rVarFileMap.values())
                    if (!outChartFile.exists())
                    {
                        chartMissing = true;
                        break;
                    }
            }
            if (force || chartMissing)
            {
                ApplicationContext.setMessage("Creating charts in R...");
                String rPlotsCommand =
                        RInterface.readResourceFile("/org/fhcrc/cpl/viewer/qa/QA_plots.R");
                Map<String, Object> rScalarVarMap = new HashMap<String, Object>();
                for (String rVar : rVarFileMap.keySet())
                    rScalarVarMap.put(rVar, "'" + RInterface.generateRFriendlyPath(rVarFileMap.get(rVar)) + "'");
                rScalarVarMap.put("label_mass_diff", labelMassDiff);
                String qaDirForR = "'" + RInterface.generateRFriendlyPath(qaDir) +
                        RInterface.generateRFriendlyPath(File.separator) + "'";
                rScalarVarMap.put("qa_dir",qaDirForR);
//System.err.println("**" + qaDir.getAbsolutePath() + ", " + RInterface.generateRFriendlyPath(qaDir) + ", " +
 //                       RInterface.generateRFriendlyPath(File.separator));
//System.err.println("***" + qaDirForR);                

                RInterface.evaluateRExpression(rPlotsCommand,
                        rScalarVarMap, null, null, null, 1200000);
                ApplicationContext.setMessage("Done creating charts in R");
            }
            else
            {
                ApplicationContext.setMessage("R charts all exist already");

            }




        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Error processing MS/MS", e);
        }
    }

    protected File createSummaryText(File summaryTextFile)
            throws IOException, XMLStreamException
    {
        PrintWriter pw = new PrintWriter(summaryTextFile);
//        String headerLine = "*fractions\ttotal_peptide\tunique_peptide\tprotein_groups\tgene_symbols\ttotal_quantified_peptide\tunique_quantified_peptide\tquantified_protein\tquantified_genes\tlight_labeled\theavy_labeled\tpartial";
        String headerLine = "*fractions\ttotal_peptide\tunique_peptide\tprotein_groups\tgene_symbols\ttotal_quantified_peptide\tunique_quantified_peptide\tquantified_proteins\tquantified_genes";

        pw.println(headerLine);

        Set<String> allGenes = new HashSet<String>();
        Set<String> allProteins = new HashSet<String>();
        Set<String> allQuantGenes = new HashSet<String>();
        Set<String> allQuantProteins = new HashSet<String>();
        Set<String> allPeptides = new HashSet<String>();
        Set<String> allQuantPeptides = new HashSet<String>();


        for (int i=0; i< runSummaryDataMaps.size(); i++)
        {
            Map<String,String> runDataMap = runSummaryDataMaps.get(i);
            pw.println(runDataMap.get("fractions") + "\t" +
                    runDataMap.get("total_peptide") + "\t" +
                    runDataMap.get("unique_peptide") + "\t" +
                    runDataMap.get("protein_groups") + "\t" +
                    runDataMap.get("gene_symbols") + "\t" +
                    runDataMap.get("total_quantified_peptide") + "\t" +
                    runDataMap.get("unique_quantified_peptide") + "\t" +
                    runDataMap.get("quantified_proteins") + "\t" +
                    runDataMap.get("quantified_genes") + "\t" +
//                     runDataMap.get("light_labeled") + "\t" +
//                     runDataMap.get("heavy_labeled") + "\t" +
//                     runDataMap.get("partial")
                    "");
            pw.flush();

            allPeptides.addAll(peptidesInEachRun.get(i));
            allQuantPeptides.addAll(quantPeptidesInEachRun.get(i));


            allProteins.addAll(proteinsInEachRun.get(i));
            allGenes.addAll(genesInEachRun.get(i));
            allQuantProteins.addAll(quantProteinsInEachRun.get(i));
            allQuantGenes.addAll(quantGenesInEachRun.get(i));
        }

        pw.println("total\t" + numTotalPeptideIdentifications + "\t" + allPeptides.size() + "\t" +
                   allProteins.size() + "\t" + allGenes.size() + "\t" +  numTotalQuantitatedPeptideIdentifications
                   + "\t" + allQuantPeptides.size() + 
                   "\t" + allQuantProteins.size() + "\t" + allQuantGenes.size());


        pw.close();

        return summaryTextFile;
    }

    protected File createAllPepText(File allPepTextFile)
            throws IOException, XMLStreamException
    {

        PrintWriter pw = new PrintWriter(allPepTextFile);
        String headerLine = "*fraction\tZ\tPrecursor_Mass\tDelta_Mass\tPepProphet\tLight_Area\tHeavy_Area\tDecimal_Ratio\tHeavy_Mass\tLight_Mass\tPeptide\tProtein";

        pw.println(headerLine);

        PepXmlLoader loader = new PepXmlLoader(allPepXmlFile, null);
        PepXmlLoader.FractionIterator fi = loader.getFractionIterator();

        while (fi.hasNext())
        {
            PepXmlLoader.PepXmlFraction fraction = (PepXmlLoader.PepXmlFraction) fi.next();
            String runName = fraction.getDataBasename();

            PepXMLFeatureFileHandler fileHandler = new PepXMLFeatureFileHandler();

            FeatureSet blah = new FeatureSet();
            List<Feature> features = fileHandler.getFeaturesFromPepXmlFraction(fraction, loader, blah);

            Set<String> uniquePeptides = new HashSet<String>();
            Set<String> uniqueQuantifiedPeptides = new HashSet<String>();

            int numQuantifiedFeatures = 0;
            for (Feature feature : features)
            {
                float peptideProphet = (float) MS2ExtraInfoDef.getPeptideProphet(feature);
                if (peptideProphet < minPeptideProphet)
                    continue;
                numTotalPeptideIdentifications++;
                StringBuffer lineBuf = new StringBuffer(runName);

                lineBuf.append("\t" + feature.getCharge());
                float deltaMass = MS2ExtraInfoDef.getDeltaMass(feature);
                lineBuf.append("\t" + (feature.getMass() + deltaMass));
                lineBuf.append("\t" + deltaMass);
                lineBuf.append("\t" + peptideProphet);
                double lightIntensity = IsotopicLabelExtraInfoDef.getLightIntensity(feature);
                lineBuf.append("\t" + (lightIntensity > 0 ? lightIntensity : "-666"));
                double heavyIntensity = IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
                lineBuf.append("\t" + (heavyIntensity > 0 ? heavyIntensity : "-666"));
                double ratio = IsotopicLabelExtraInfoDef.getRatio(feature);
                lineBuf.append("\t" + (ratio > 0 ? ratio : "-666"));

                double heavyMass = IsotopicLabelExtraInfoDef.getHeavyMass(feature);
                double lightMass = IsotopicLabelExtraInfoDef.getLightMass(feature);

                if (heavyMass > 0 && lightMass > 0)
                {
                    lineBuf.append("\t" + heavyMass + "\t" + lightMass);
                }
                else
                {
                    lineBuf.append("\t-666\t-666");
                }



                String peptideSequence = MS2ExtraInfoDef.getFirstPeptide(feature);
                lineBuf.append("\t" + peptideSequence);
                uniquePeptides.add(peptideSequence);
                if (ratio>0)
                {
                    uniqueQuantifiedPeptides.add(peptideSequence);
                    numQuantifiedFeatures++;
                }

                //protein(s)
                lineBuf.append("\t");
                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                if (proteins != null)
                    for (int i=0; i<proteins.size(); i++)
                    {
                        if (i>0)
                            lineBuf.append(";");
                        lineBuf.append(proteins.get(i));
                    }
                pw.println(lineBuf);
                pw.flush();
            }

            Map<String,String> runDataMap = new HashMap<String,String>();
            runDataMap.put("fractions", runName);
            runDataMap.put("total_peptide", "" + features.size());
            runDataMap.put("unique_peptide", "" + uniquePeptides.size());
            runDataMap.put("total_quantified_peptide", "" + numQuantifiedFeatures);
            runDataMap.put("unique_quantified_peptide", "" + uniqueQuantifiedPeptides.size());
            runSummaryDataMaps.add(runDataMap);
            peptidesInEachRun.add(uniquePeptides);
            quantPeptidesInEachRun.add(uniqueQuantifiedPeptides);


            numTotalQuantitatedPeptideIdentifications += numQuantifiedFeatures;
        }

        pw.flush();
        pw.close();

        ApplicationContext.setMessage("Saved file " + allPepTextFile.getAbsolutePath());
        return allPepTextFile;
    }

    protected File createAllProtText(File allProtTextFile)
            throws Exception
    {

        ApplicationContext.setMessage("Loading gene mapping file...");
        Map<String,List<String>> ipiGeneListMap = QAUtilities.loadIpiGeneListMap(protGeneFile);

        ApplicationContext.setMessage("Gene map loaded. " + ipiGeneListMap.size() + " proteins with genes");



        String headerLine = "*Group\tGroup_Probability\tProtein_Probability\tL2H_Mean\tL2H_StdDev\tRatio_Peps\tH2L_Mean\tH2L_StdDev\tnum_Indistinguishable_Proteins\tIndistinguishable_Proteins\tnum_Genes\tGene_Symbol\tPeptides";

        PrintWriter pw = new PrintWriter(allProtTextFile);
        pw.println(headerLine);

        ProtXmlReader protXmlReader = new ProtXmlReader(allProtXmlFile);
        ProtXmlReader.ProteinGroupIterator groupIterator = protXmlReader.iterator();

        for (int i=0; i<peptidesInEachRun.size(); i++)
        {
            proteinsInEachRun.add(new HashSet<String>());
            genesInEachRun.add(new HashSet<String>());
            quantProteinsInEachRun.add(new HashSet<String>());
            quantGenesInEachRun.add(new HashSet<String>());
        }

        while (groupIterator.hasNext())
        {
            ProteinGroup proteinGroup = groupIterator.next();

            List<ProtXmlReader.Protein> proteins = proteinGroup.getProteins();

            for (int i=0; i<proteins.size(); i++)
            {
                //For proteins after the first, "group" number is actually group number with _# after it, where _# is
                //the count
                StringBuffer lineBuf = new StringBuffer("" + proteinGroup.getGroupNumber());
                if (i>0)
                    lineBuf.append("_" + i);

                ProtXmlReader.Protein protein = proteins.get(i);

                lineBuf.append("\t" + proteinGroup.getGroupProbability());
                lineBuf.append("\t" + protein.getProbability());

                ProtXmlReader.QuantitationRatio ratio = protein.getQuantitationRatio();
                if (ratio == null)
                    lineBuf.append("\t-666\t-666\t-666\t-666\t-666");
                else
                {
                    lineBuf.append("\t" +  ratio.getRatioMean());
                    lineBuf.append("\t" + ratio.getRatioStandardDev());
                    lineBuf.append("\t" + ratio.getRatioNumberPeptides());
                    lineBuf.append("\t" +  ratio.getHeavy2lightRatioMean());
                    lineBuf.append("\t" + ratio.getHeavy2lightRatioStandardDev());
                }

                StringBuffer indistProteinNamesBuf = new StringBuffer("");
                List<String> indistProteinNames = protein.getIndistinguishableProteinNames();
                if (indistProteinNames == null)
                    indistProteinNames = new ArrayList<String>(1);
                indistProteinNames.add(protein.getProteinName());
                for (int j=0; j<indistProteinNames.size(); j++)
                {
                    if (j>0)
                        indistProteinNamesBuf.append(";");
                    indistProteinNamesBuf.append(indistProteinNames.get(j));

                }
                lineBuf.append("\t" + indistProteinNames.size());                
                lineBuf.append("\t" + indistProteinNamesBuf);


                //num_genes, gene_symbol
                Set<String> genesAllIndistProteins = new HashSet<String>();
                for (String indistProteinName : indistProteinNames)
                {
                    List<String> geneList = ipiGeneListMap.get(indistProteinName);
                    if (geneList != null)
                        genesAllIndistProteins.addAll(geneList);
                }
                lineBuf.append("\t" + genesAllIndistProteins.size());
                StringBuffer geneSymbolBuf = new StringBuffer("");
                boolean firstGene = true;
                for (String gene : genesAllIndistProteins)
                {
                    if (!firstGene)
                        geneSymbolBuf.append(";");
                    geneSymbolBuf.append(gene);
                    firstGene = false;
                }
                if (geneSymbolBuf.length() == 0)
                    geneSymbolBuf.append("NA");
                lineBuf.append("\t" + geneSymbolBuf);


                StringBuffer peptidesBuf = new StringBuffer("");
                List<ProtXmlReader.Peptide> peptides = protein.getPeptides();
                int numWrittenPeptides = 0;
                Set<String> peptidesThisProtein = new HashSet<String>();
                for (ProtXmlReader.Peptide peptide : peptides)
                {
                    if (peptide.isContributingEvidence() && peptide.isNondegenerateEvidence())
                    {
                        if (numWrittenPeptides>0)
                            peptidesBuf.append(";");
                        peptidesBuf.append(peptide.getPeptideSequence());
                        peptidesThisProtein.add(peptide.getPeptideSequence());
                        numWrittenPeptides++;
                    }
                }

                for (int j=0; j<peptidesInEachRun.size(); j++)
                {
                    Set<String> peptidesInRun = peptidesInEachRun.get(j);
                    for (String peptide : peptidesThisProtein)
                    {
                        if (peptidesInRun.contains(peptide))
                        {
                            proteinsInEachRun.get(j).add(protein.getProteinName());
                            genesInEachRun.get(j).addAll(genesAllIndistProteins);
                            if (ratio != null)
                            {
                                quantProteinsInEachRun.get(j).add(protein.getProteinName());
                                quantGenesInEachRun.get(j).addAll(genesAllIndistProteins);
                            }
                            break;
                        }
                    }
                }

                if (peptidesBuf.length() == 0)
                    peptidesBuf.append("NA");
                lineBuf.append("\t" + peptidesBuf);



                pw.println(lineBuf);
                pw.flush();
            }
        }
        pw.flush();
        pw.close();

        for (int i=0; i<runSummaryDataMaps.size(); i++)
        {
            Map<String,String> runSummaryMap = runSummaryDataMaps.get(i);
            runSummaryMap.put("protein_groups", "" + proteinsInEachRun.get(i).size());
            runSummaryMap.put("gene_symbols", "" + genesInEachRun.get(i).size());
            runSummaryMap.put("quantified_proteins", "" + quantProteinsInEachRun.get(i).size());
            runSummaryMap.put("quantified_genes", "" + quantGenesInEachRun.get(i).size());            
        }

        ApplicationContext.setMessage("Saved file " + allProtTextFile.getAbsolutePath());
        return allProtTextFile;
    }

}
