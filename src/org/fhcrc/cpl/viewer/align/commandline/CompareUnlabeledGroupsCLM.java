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

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.viewer.align.BucketedPeptideArray;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * This commandline module bundles together several different actions (see help) involved with
 * comparing unlabeled data from unlabeled samples from multiple experimental groups in a peptide array.
 * Initially developed based on work for Teri Brentnall's group.
 *
 * dhmay 20100730 
 */
public class CompareUnlabeledGroupsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CompareUnlabeledGroupsCLM.class);

    protected File file;
    protected File outProteinQFile;
    protected File outGSEAFile;

    protected File outPepXMLFile;
    protected File outDir;
    
    protected boolean showCharts = false;
    protected boolean runRefreshParser = false;

    List<Protein> fastaProteins;

    Object[] arrayRows;

    protected File fastaFile;
    protected int maxMissedCleavages = 2;

    PeptideArrayAnalyzer arrayAnalyzer;

    double[] caseRunIndexes1Based;
    double[] controlRunIndexes1Based;

    protected Map<String, Collection<String>> peptideProteinMap = null;

    //Dummy isotopic label definition for creating features with ratios
    protected AnalyzeICAT.IsotopicLabel dummyLabel =
            new AnalyzeICAT.IsotopicLabel("Acrylamide (D0/D3)",71.03657F,3.0188F,'C',
                    AnalyzeICAT.DEFAULT_MAX_LABEL_COUNT);

    public CompareUnlabeledGroupsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "compareunlabeledgroups";

        mShortDescription = "Analyzes a peptide array, comparing two groups of unlabeled samples";
        mHelpMessage = "This single command pulls together several different steps of a complicated workflow " +
                "for analyzing a peptide array and comparing proteins in a 'case' and 'control' group based on " +
                "peptideintensity ratios.  The steps are: create a pepXML file describing ratios, create a " +
                "Geneset Enrichment file containing peptide-level t-scores, create a protein-level file containing " +
                "protein p-values and q-values.";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, null),
                        new FileToWriteArgumentDefinition("outgsea",false, "output GSEA file"),
                        new FileToWriteArgumentDefinition("outpepxml",false, "output pepXML file"),
                        new FileToWriteArgumentDefinition("outproteinq",false, "output protein q-value file"),


                        new FileToWriteArgumentDefinition("outdir",false, "output directory"),

                        new FileToReadArgumentDefinition("caserunlistfile",true,
                                "File containing the names of runs in the case group, one per line"),
                        new FileToReadArgumentDefinition("controlrunlistfile",true,
                                "File containing the names of runs in the control group, one per line"),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new IntegerArgumentDefinition("minrunspergroup", false,
                                "Minimum number of runs in each group in which a feature must be located to be counted",2),
                        new FileToReadArgumentDefinition("fasta", true,
                                "FASTA database"),
                        new FileToReadArgumentDefinition("protxml", false,
                                "protXML search results file"),
                        new IntegerArgumentDefinition("missedcleavages", false,"maximum missed cleavages for peptides (if fasta provided)", maxMissedCleavages),
                        new BooleanArgumentDefinition("refreshparser",false,"Run RefreshParser on pepXML file?", runRefreshParser),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        outGSEAFile = getFileArgumentValue("outgsea");
        outPepXMLFile = getFileArgumentValue("outpepxml");
        outProteinQFile = getFileArgumentValue("outproteinq");

        outDir = getFileArgumentValue("outdir");

        showCharts = getBooleanArgumentValue("showcharts");
        runRefreshParser = getBooleanArgumentValue("refreshparser");


        fastaFile = getFileArgumentValue("fasta");
        maxMissedCleavages = getIntegerArgumentValue("missedcleavages");
        File protXmlFile = getFileArgumentValue("protxml");

        try
        {
            ApplicationContext.infoMessage("Loading FASTA....");
            fastaProteins = ProteinUtilities.loadProteinsFromFasta(fastaFile);
            ApplicationContext.infoMessage("Mapping peptides to proteins...");
            if (protXmlFile == null)
            {
                //really stupid to have to do this, but for some reason map<string,set<string>>
                // doesn't cast to map<string,collection<string>> .  Same with list->collection
//                peptideProteinMap = new HashMap<String,Collection<String>>(ProteinUtilities.loadTrypticPeptideProteinMapFromFasta(fastaFile, maxMissedCleavages));
                peptideProteinMap = new HashMap<String,Collection<String>>();
                for (Protein protein : fastaProteins)
                {
                    for (String peptide : ProteinUtilities.getTrypticPeptidesForProtein(protein, maxMissedCleavages))
                    {
                        Collection<String> proteinsThisPeptide = peptideProteinMap.get(peptide);
                        if (proteinsThisPeptide == null)
                        {
                            proteinsThisPeptide = new ArrayList<String>();
                            peptideProteinMap.put(peptide, proteinsThisPeptide);
                        }
                        proteinsThisPeptide.add(protein.getLookup());
                    }
                }
                ApplicationContext.infoMessage("Loaded peptide-protein map from FASTA.  Warning: protein summary " +
                        "may include peptides that ProteinProphet would not assign");
            }
            else
            {
                ApplicationContext.infoMessage("Loading peptide-protein map from protxml....");                
                peptideProteinMap = new HashMap<String,Collection<String>>(ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFile, 0.1));
                ApplicationContext.infoMessage("Loaded peptide-protein map from protXML file");
            }
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("failed to load peptide-protein mapping",e);
        }

        try
        {
            arrayAnalyzer = new PeptideArrayAnalyzer(getUnnamedFileArgumentValue());
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("failed to load peptide array",e);
        }

        try
        {
            arrayAnalyzer.loadCaseControlRunListsFromFiles(getFileArgumentValue("caserunlistfile"),
                    getFileArgumentValue("controlrunlistfile"));
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load case or control run list file", e);
        }

        ApplicationContext.infoMessage("Case runs:");
        for (String caseRun : arrayAnalyzer.getCaseRunNames())
            ApplicationContext.setMessage("\t" + caseRun);
        ApplicationContext.infoMessage("Control runs:");
        for (String controlRun : arrayAnalyzer.getControlRunNames())
            ApplicationContext.setMessage("\t" + controlRun);
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            //first create feature file with case-control ratios
            List<Feature> features = new ArrayList<Feature>();
            for (Map<String, Object> rowMap : arrayAnalyzer.getRowMaps())
            {
                Feature feature = createFeatureWithPeptideAndRatio(rowMap);
                if (feature != null)
                    features.add(feature);
            }
            FeatureSet featureSet = new FeatureSet(features.toArray(new Feature[features.size()]));
            if (fastaFile != null)
            {
                MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSet, fastaFile.getAbsolutePath());
            }
            ProteinUtilities.guessProteinsForFeaturePeptides(new FeatureSet[] {featureSet}, fastaFile,
                    fastaProteins.toArray(new Protein[fastaProteins.size()]));
            List<Feature> featuresWithPeptides = new ArrayList<Feature>();
            for (Feature feature : featureSet.getFeatures())
                if (MS2ExtraInfoDef.getFirstProtein(feature) != null)
                    featuresWithPeptides.add(feature);
            if (featuresWithPeptides.size() < featureSet.getFeatures().length)
            {
                ApplicationContext.infoMessage("WARNING: " + (featureSet.getFeatures().length -
                        featuresWithPeptides.size()) + " peptides did not have proteins in fasta file, were removed");

                featureSet.setFeatures(featuresWithPeptides.toArray(new Feature[featuresWithPeptides.size()]));
            }
            try
            {
                featureSet.savePepXml(outPepXMLFile);
                ApplicationContext.infoMessage("Saved feature file " + outGSEAFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to save featureset",e);
            }

            //run RefreshParser
            //todo: parameterize this, since requires RefreshParser on PATH
            if (runRefreshParser)
            {
                ApplicationContext.setMessage("Running RefreshParser on all files...");

                String cmd = "RefreshParser " + outPepXMLFile.getAbsolutePath() + " " + fastaFile.getAbsolutePath();
                _log.debug("Running RefreshParser: " + cmd);
                Process p = Runtime.getRuntime().exec(cmd ,null);
                try
                {
                    int err = p.waitFor();
                    _log.debug("process returned, "+err);
                    if (err == 0)
                        ApplicationContext.setMessage("Successfully ran RefreshParser on " +
                                outPepXMLFile.getAbsolutePath());
                    else
                        throw new CommandLineModuleExecutionException("RefreshParser failed on " +
                                outPepXMLFile.getAbsolutePath() + " with error code " + err);
                }
                catch (InterruptedException e)
                {
                    ApplicationContext.infoMessage("Warning: InterruptedException while running RefreshParser.  Check that it completed.");
                }
            }



            float[][] resultMatrix = arrayAnalyzer.runTwoGroupTTestAllRows(2);

            Map<String, List<Float>> peptideTScoresMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideRatiosMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideCaseMeanIntensitiesMap = new HashMap<String, List<Float>>();
            Map<String, List<Float>> peptideControlMeanIntensitiesMap = new HashMap<String, List<Float>>();

            List<Float> noNaP = new ArrayList<Float>();
            List<Float> noNaQ = new ArrayList<Float>();
            List<Float> noNaT = new ArrayList<Float>();

            for (int i=0; i<resultMatrix.length; i++)
            {
                if (!Float.isNaN(resultMatrix[i][0]))
                {
                    noNaP.add(resultMatrix[i][0]);
                    noNaT.add(resultMatrix[i][1]);
                    noNaQ.add(resultMatrix[i][2]);
                }
            }

            if (showCharts)
            {
                new PanelWithHistogram(noNaP, "pep p-values").displayInTab();
                new PanelWithHistogram(noNaQ, "pep q-values").displayInTab();
                new PanelWithHistogram(noNaT, "pep t-scores").displayInTab();
            }

            Map<String, List<Integer>> peptideRowIndexesMap = arrayAnalyzer.loadPeptideRowIndexesMap();
            for (String peptide : peptideRowIndexesMap.keySet())
            {
                if (!peptideProteinMap.containsKey(peptide))
                {
                    _log.debug("Skipping unknown peptide " + peptide);
                    continue;
                }

                List<Float> tScoresThisPeptide = new ArrayList<Float>();
                peptideTScoresMap.put(peptide, tScoresThisPeptide);

                List<Float> ratiosThisPeptide = new ArrayList<Float>();
                peptideRatiosMap.put(peptide, ratiosThisPeptide);

                List<Float> caseMeanIntensitiesThisPeptide = new ArrayList<Float>();
                peptideCaseMeanIntensitiesMap.put(peptide, caseMeanIntensitiesThisPeptide);

                List<Float> controlMeanIntensitiesThisPeptide = new ArrayList<Float>();
                peptideControlMeanIntensitiesMap.put(peptide, controlMeanIntensitiesThisPeptide);

                for (int i : peptideRowIndexesMap.get(peptide))
                {
                    if (!Float.isNaN(resultMatrix[i][1]))
                    {
                        tScoresThisPeptide.add(resultMatrix[i][1]);
                        ratiosThisPeptide.add(
                                arrayAnalyzer.calcRowCaseControlRatioOfGeomMeans(
                                        arrayAnalyzer.getRowMap(i)).floatValue());
                        caseMeanIntensitiesThisPeptide.add((float)
                                BasicStatistics.geometricMean(arrayAnalyzer.getPresentCaseIntensities(
                                        arrayAnalyzer.getRowMap(i))));
                        controlMeanIntensitiesThisPeptide.add((float)
                                BasicStatistics.geometricMean(arrayAnalyzer.getPresentControlIntensities(
                                        arrayAnalyzer.getRowMap(i))));
                    }
                }
            }

            Set<String> observedProteinsSet = new HashSet<String>();
            for (String peptide : peptideTScoresMap.keySet())
            {
                if (!peptideTScoresMap.get(peptide).isEmpty())
                {
                    //This is a bit hacky.  Collapsing to a single peptide.  Rather than explicitly
                    //change the data model, I'm just holding a single value per peptide, the median of all the
                    //non-NaN values we had
                    List<Float> newList = new ArrayList<Float>();
                    newList.add((float) BasicStatistics.median(peptideTScoresMap.get(peptide)));
                    peptideTScoresMap.put(peptide, newList);

                    List<Float> newListR = new ArrayList<Float>();
                    newListR.add((float) BasicStatistics.geometricMean(peptideRatiosMap.get(peptide)));
                    peptideRatiosMap.put(peptide, newListR);

                    List<Float> newListCaseI = new ArrayList<Float>();
                    newListCaseI.add((float) BasicStatistics.geometricMean(peptideCaseMeanIntensitiesMap.get(peptide)));
                    peptideCaseMeanIntensitiesMap.put(peptide, newListCaseI);

                    List<Float> newListControlI = new ArrayList<Float>();
                    newListControlI.add((float) BasicStatistics.geometricMean(peptideControlMeanIntensitiesMap.get(peptide)));
                    peptideControlMeanIntensitiesMap.put(peptide, newListControlI);
                    observedProteinsSet.addAll(peptideProteinMap.get(peptide));
                }
            }
            List<String> observedProteins = new ArrayList<String>(observedProteinsSet);

            PrintWriter outPW = null;
            try
            {
                outPW = new PrintWriter(outGSEAFile);
            }
            catch (FileNotFoundException e)
            {
                throw new CommandLineModuleExecutionException("Failed to open output file " + outGSEAFile.getAbsolutePath());
            }
            StringBuffer headerLineBuf = new StringBuffer("peptide\tratio\ttscore\tintmeancase\tintmeancontrol");
            for (String protein : observedProteins)
            {
                headerLineBuf.append("\t" + protein);
            }
            outPW.println(headerLineBuf);
            outPW.flush();

            for (String peptide : peptideTScoresMap.keySet())
            {
                StringBuffer proteinsPartBuf = new StringBuffer();
                Collection<String> proteinsThisPeptide = peptideProteinMap.get(peptide);
                for (String protein : observedProteins)
                {
                    if (proteinsThisPeptide.contains(protein))
                        proteinsPartBuf.append("\t1");
                    else
                        proteinsPartBuf.append("\t0");
                }
                for (int i=0; i<peptideTScoresMap.get(peptide).size(); i++)
                {
                    float tscore = peptideTScoresMap.get(peptide).get(i);
                    float ratio = peptideRatiosMap.get(peptide).get(i);
                    float caseMean = peptideCaseMeanIntensitiesMap.get(peptide).get(i);
                    float controlMean = peptideControlMeanIntensitiesMap.get(peptide).get(i);

                    if (!Float.isNaN(tscore))
                        outPW.println(peptide + "\t" + ratio + "\t" + tscore + "\t" + caseMean + "\t" + controlMean + proteinsPartBuf.toString());
                }
            }
            outPW.close();
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Error running t test", e);
        }

        //Run protein overexpression analysis on peptide-level t-scores
        InputStream in = BucketedPeptideArray.class.getResourceAsStream("overrepresentation_groupcompare.R");
        StringBuffer scriptBuf = new StringBuffer();
        byte[] buf = new byte[1024];
        int len;
        try
        {
            while ((len = in.read(buf)) > 0)  scriptBuf.append(new String(buf, 0, len));
            in.close();
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(
                    "I/O Error running R script for protein-level Wilcox test", e);
        }
        Map<String, Object> scalarVariableValues = new HashMap<String, Object>();
        scalarVariableValues.put("file", "'" + RInterface.generateRFriendlyPath(outGSEAFile) + "'");
        scalarVariableValues.put("outfile", "'" + RInterface.generateRFriendlyPath(outProteinQFile) + "'");
        ApplicationContext.infoMessage("Running R script to calculate protein q-values...");
        RInterface.evaluateRExpression(scriptBuf.toString(),  scalarVariableValues, null, null,
                new String[] {"qvalue"}, 500000);
        ApplicationContext.infoMessage("R complete, wrote file " + outProteinQFile.getAbsolutePath()); 
    }

    /**
     *
     * @param rowMap
     * @return null if can't create feature
     */
    protected Feature createFeatureWithPeptideAndRatio(Map<String, Object> rowMap)
    {
        Feature rowFeature = new Feature();
        Set<String> rowPeptides = arrayAnalyzer.getAllRowPeptides(rowMap);

        //bomb out if no peptide or more than one
        if (rowPeptides.size() != 1)
            return null;

        String peptide = rowPeptides.iterator().next();
        rowFeature.setCharge(Integer.parseInt(rowMap.get("charge").toString()));
        if (rowMap.containsKey("minScan"))
            rowFeature.setScan((int) Double.parseDouble(rowMap.get("minScan").toString()));
        else
        {
            rowFeature.setScan((int) Double.parseDouble(rowMap.get("minTime").toString()));
        }
        MS2ExtraInfoDef.setPeptideList(rowFeature, peptide);
        MS2ExtraInfoDef.setPeptideProphet(rowFeature, 0.95);
        Protein dummyProtein = new Protein("asdf", peptide.getBytes());
        rowFeature.setMass((float) new Peptide(dummyProtein, 0, peptide.length()).getMonoisotopicMass());
        rowFeature.updateMz();

        List<Double> caseIntensities = arrayAnalyzer.getPresentCaseIntensities(rowMap);
        List<Double> controlIntensities = arrayAnalyzer.getPresentControlIntensities(rowMap);

//if (peptide.equals("FLGVAEQLHNEGFK")) System.err.println("FLGVAEQLHNEGFK, " + caseIntensities.size() +", " +  controlIntensities.size());
        if (caseIntensities.size() >= 1 && controlIntensities.size() >= 1)
        {
            double geoMeanCase = BasicStatistics.geometricMean(caseIntensities);
            double geoMeanControl = BasicStatistics.geometricMean(controlIntensities);
            //only add a ratio if at least one of the intensities being compared is nonzero
            if (Math.max(geoMeanCase, geoMeanControl) > 0)
            {
                double ratio = geoMeanCase / geoMeanControl;
                IsotopicLabelExtraInfoDef.setRatio(rowFeature, ratio);
                IsotopicLabelExtraInfoDef.setLabel(rowFeature, dummyLabel);
                IsotopicLabelExtraInfoDef.setHeavyFirstScan(rowFeature, rowFeature.getScan());
                IsotopicLabelExtraInfoDef.setHeavyLastScan(rowFeature, rowFeature.getScan());
                IsotopicLabelExtraInfoDef.setLightFirstScan(rowFeature, rowFeature.getScan());
                IsotopicLabelExtraInfoDef.setLightLastScan(rowFeature, rowFeature.getScan());
                IsotopicLabelExtraInfoDef.setLightIntensity(rowFeature, geoMeanCase);
                IsotopicLabelExtraInfoDef.setHeavyIntensity(rowFeature, geoMeanControl);
//if (peptide.equals("FLGVAEQLHNEGFK")) System.err.println("\tratio: " + ratio);
            }
        }
        return rowFeature;
    }

}

