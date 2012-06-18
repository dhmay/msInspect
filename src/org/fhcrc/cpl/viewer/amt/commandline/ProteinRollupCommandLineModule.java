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
package org.fhcrc.cpl.viewer.amt.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.viewer.ms2.GeneMappingUtilities;
import org.fhcrc.cpl.toolbox.filehandler.TabWriter;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule to help roll up peptide information to the protein level.
 * Includes some naive mechanisms for generating lists of possible proteins given
 * a list of peptides and a Fasta file.
 * More usefully, includes some mechanisms for generating pepXml files that can be
 * used by ProteinProphet, given some peptide IDs and a Fasta file
 */
public class ProteinRollupCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProteinRollupCommandLineModule.class);

    protected int mode=-1;                  

    protected final static String[] modeStrings =
            {"digestprotxml","tabfile","createpepxmlfromtabfile",
                    "collapsepeptides","rollupgenes","rollupgenesfrompeptides"};

    protected final static String[] modeExplanations =
            {
                    "Digest the output of ProteinProphet",
                    "Take as input a single tab-delimited ratio file in a special format",
                    "Create a pepXML file, using an input tab-delimited ratio file",
                    "collapse per-peptide row entries into one per protein",
                    "Roll up proteins to gene level",
                    "Rollup peptides to proteins and then to gene level in one step",
            };

    protected static final int MODE_DIGEST_PROTXML = 0;
    protected static final int MODE_TABFILE = 1;
    protected static final int MODE_CREATE_PEPXML_FOR_TABFILE = 2;
    protected static final int MODE_COLLAPSE_PEPTIDES = 3;
    protected static final int MODE_ROLLUP_GENES = 4;
    protected static final int MODE_ROLLUP_PEPTIDES_TO_GENES = 5;


    protected FeatureSet[] ms1FeatureSets;
    protected Protein[] proteinsInFasta;
    protected File fastaFile;
    protected File[] protXmlFiles;
    protected FeatureSet pepXmlFeatureSet;
    protected int minUniquePeptides=1;
    protected double minProtXmlProbability=0;

    protected File outFile;
    protected File outDir;
    protected File inDir;
    protected File inTabFile;
    protected File geneLookupFile;

    protected boolean expanded = false;


    protected int minPeptides = 1;


    boolean showCharts = false;


    protected AmtDatabase amtDB;

    //protected boolean countMode=false;

    public ProteinRollupCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "proteinrollup";

        mShortDescription = "Tools for rolling up peptide information to the protein level.";
        mHelpMessage = "Provides a naive mechanism for generating lists of possible proteins given\n" +
                " a list of peptides and a Fasta file.\n";


        CommandLineArgumentDefinition[] argDefs =
               {
                    createUnnamedSeriesFileArgumentDefinition(false, "Input feature file(s)"),
                    new EnumeratedValuesArgumentDefinition("mode",true,modeStrings, modeExplanations),
                    new FileToReadArgumentDefinition("fasta",false,"Input protein database"),
                    new FileToReadArgumentDefinition("peptidelistfile",false, "Input file containing one peptide sequence per line"),
                    new FileToWriteArgumentDefinition("out",false, "output filepath"),
                       new DirectoryToWriteArgumentDefinition("outdir",false,
                               "output directory"),
                    new IntegerArgumentDefinition("minuniquepeptides",false,"Minimum number of unique peptides for a protein to be considered to be identified in a protXml file",minUniquePeptides),
                    new DecimalArgumentDefinition("minProtXmlProbability",false,"Minimum probability assigned by protXml required for the protein to be considered to be identified",minProtXmlProbability),
                    new FileToReadArgumentDefinition("amtxml",false, "AMT database file (for spectral count info)"),
                       new DirectoryToReadArgumentDefinition("indir", false, "Directory of input files"),
//new BooleanArgumentDefinition("countmode",false, "count mode, for spectral counting"),
                       new BooleanArgumentDefinition("showcharts", false, "Show charts", showCharts),
                       new BooleanArgumentDefinition("expanded", false, "Peptide expanded format (for gene rollup)", showCharts),

                       new FeatureFileArgumentDefinition("pepxml",false,"pepXml file for use in digesting protXml file"),
                       new IntegerArgumentDefinition("minpeptides",false,"Minimum peptides for a protein to be collapsed",minPeptides),
                       new FileToReadArgumentDefinition("genelookupfile", false, "Gene lookup file for IPI numbers"),
                       new FileToReadListArgumentDefinition("protxml", false, "ProtXML file(s)"),

               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

//        countMode = getBooleanArgumentValue("countmode");

        minPeptides=getIntegerArgumentValue("minpeptides");

        File[] inFiles = getUnnamedSeriesFileArgumentValues();

        if (inFiles != null)
        {
            assertArgumentAbsent("indir");
        }
        else
        {
            assertArgumentPresent("indir");
            inDir = getFileArgumentValue("indir");
            List<File> inFileList = new ArrayList<File>();
            for (File file : inDir.listFiles())
            {
                if (!file.isDirectory())
                    inFileList.add(file);
            }
            inFiles = inFileList.toArray(new File[inFileList.size()]);
        }

        pepXmlFeatureSet = getFeatureSetArgumentValue("pepxml");
        if (mode == MODE_TABFILE)
            assertArgumentPresent("fasta");
        if (mode == MODE_ROLLUP_PEPTIDES_TO_GENES)
        {
            if (!hasArgumentValue("fasta") && !hasArgumentValue("protxml"))
            {
                throw new ArgumentValidationException("Either the 'fasta' or the 'protxml' argument is required for this mode");
            }
        }
        if (mode == MODE_TABFILE ||mode==MODE_COLLAPSE_PEPTIDES||
                mode == MODE_CREATE_PEPXML_FOR_TABFILE || mode == MODE_DIGEST_PROTXML || mode == MODE_ROLLUP_GENES ||
                mode == MODE_ROLLUP_PEPTIDES_TO_GENES)
        {
            if (inFiles.length > 1)
                throw new ArgumentValidationException("Only one input file allowed for this mode");
            if (mode == MODE_DIGEST_PROTXML)
            {
                assertArgumentPresent("protxml");
                assertArgumentPresent("pepxml");
            }
            else
            {
                inTabFile = inFiles[0];

            }



        }
        else
        {
            if (inFiles != null)
            {
                ms1FeatureSets = new FeatureSet[inFiles.length];
                for (int i=0; i<inFiles.length; i++)
                {
                    try
                    {
                        ms1FeatureSets[i] = new FeatureSet(inFiles[i]);
                    }
                    catch (Exception e)
                    {
                        throw new ArgumentValidationException(e);
                    }
                }
            }
        }

        protXmlFiles = getFileArrayArgumentValue("protxml");
        if (protXmlFiles != null)
        {
            ApplicationContext.infoMessage("protXML files:");
            for (File protXmlFile : protXmlFiles)
                ApplicationContext.infoMessage("\t" + protXmlFile.getAbsolutePath());
        }


        showCharts = getBooleanArgumentValue("showcharts");
        expanded = getBooleanArgumentValue("expanded");


        if (inFiles == null && mode != MODE_DIGEST_PROTXML)
        {
            assertArgumentPresent("peptidelistfile");
            File peptideListFile = getFileArgumentValue("peptidelistfile");
            List<Feature> peptideFeatureList = new ArrayList<Feature>();
            try
            {
                BufferedReader br = new BufferedReader(new FileReader(peptideListFile));
                int dummyScan = 1;
                while (true)
                {
                    String peptide = br.readLine();
                    if (null == peptide)
                        break;
                    if (peptide.length() == 0 || peptide.charAt(0) == '#')
                        continue;
                    Feature feature = new Feature();
                    MS2ExtraInfoDef.addPeptide(feature, peptide);
                    feature.setScan(dummyScan++);
                    feature.setCharge(1);
                    peptideFeatureList.add(feature);
                }
                ms1FeatureSets = new FeatureSet[1];
                ms1FeatureSets[0] = new FeatureSet(peptideFeatureList.toArray(new Feature[0]));
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }

        minUniquePeptides = getIntegerArgumentValue("minuniquepeptides");
        minProtXmlProbability = getDoubleArgumentValue("minprotxmlprobability");

        if (hasArgumentValue("fasta"))
        {
            fastaFile = getFileArgumentValue("fasta");
            try
            {
                proteinsInFasta = ProteinUtilities.loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Error loading fasta file " + fastaFile, true);
            }
            if (proteinsInFasta.length == 0)
            {
                throw new ArgumentValidationException("No proteins found in fasta file " + fastaFile, false);
            }
        }
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (inFiles.length > 2)
        {
            assertArgumentAbsent("out");
        }
        if (hasArgumentValue("out"))
            assertArgumentAbsent("outdir");
        else
            assertArgumentPresent("outdir");

        if (hasArgumentValue("amtxml"))
        {
            File amtDBFile = getFileArgumentValue("amtxml");
            try
            {
                amtDB = new AmtXmlReader(amtDBFile).getDatabase();
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }

        geneLookupFile = getFileArgumentValue("genelookupfile");
        if (mode == MODE_ROLLUP_GENES || mode == MODE_ROLLUP_PEPTIDES_TO_GENES)
            assertArgumentPresent("genelookupfile");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
List<Feature> ms1FeaturesWithPeptidesList = new ArrayList<Feature>();
Feature[] ms1FeaturesWithPeptides;
        switch(mode)
        {
            case MODE_TABFILE:
                rollupTabFile();
                break;
            case MODE_DIGEST_PROTXML:
//                digestProtXml();
                System.err.println("Mode deactivated");
                break;
            case MODE_CREATE_PEPXML_FOR_TABFILE:
                createPepXmlForTabFile();
                break;
            case MODE_COLLAPSE_PEPTIDES:
//                collapsePeptides();
                System.err.println("Mode deactivated");               
                break;
            case MODE_ROLLUP_GENES:
                System.err.println("****THIS MODE IS OFFLINE.  WANT IT?  PUT IT BACK TOGETHER*****");
//                proteinGeneRollup();
                break;
            case MODE_ROLLUP_PEPTIDES_TO_GENES:
                if (expanded)
                    peptideGeneRollupExpanded();
                else
                    peptideGeneRollup();
                break;
        }
    }

    protected void peptideGeneRollupExpanded() throws CommandLineModuleExecutionException
    {
        try
        {
            Map<String,List<String>> ipiGeneArrayMap =
                    GeneMappingUtilities.loadIPIGeneMap(geneLookupFile);

            TabLoader tl = new TabLoader(inTabFile);
            Map[] rowsAsMaps = (Map[]) tl.load();

            Set<String> allPeptides = new HashSet<String>();
            for (int i=0; i<rowsAsMaps.length; i++)
            {
                allPeptides.add((String) rowsAsMaps[i].get("peptide"));
            }

            Map<String,List<Protein>> peptideProteinMap =
                    ProteinUtilities.mapPeptidesToProteins(allPeptides, protXmlFiles, proteinsInFasta,0);
            Map<String,Object>[] outRows = new Map[rowsAsMaps.length];
            String[] outColumns = new String[]
                    {
                            "peptide",
                            "ratio",
                            "intensity1","intensity2",
                            "proteins",
                            "genes"
                    };
            int i=0;
            for (Map rowAsMap : rowsAsMaps)
            {
                String peptide = (String) rowAsMap.get("peptide");
                Double ratio = (Double) rowAsMap.get("ratio");
                Double intensity1 = (Double) rowAsMap.get("intensity1");
                Double intensity2 = (Double) rowAsMap.get("intensity2");

                List<Protein> proteins = peptideProteinMap.get(peptide);
                List<String> ipis = new ArrayList<String>();
                Set<String> genes = new HashSet<String>();
                if (proteins != null)
                {
                    for (Protein protein : proteins)
                    {
                        String proteinLookup = protein.getLookup();
                        ipis.add(proteinLookup);
                        if (ipiGeneArrayMap.containsKey(proteinLookup))
                        {
                            List<String> genesThisProtein = ipiGeneArrayMap.get(proteinLookup);
                            for (String gene : genesThisProtein)
                                genes.add(gene);
                        }
                    }
                }
                Map<String,Object> outRow = new HashMap<String,Object>(outColumns.length);
                String proteinsString = assembleStringFromCollection(ipis);
                String genesString = assembleStringFromCollection(genes);

                outRow.put("peptide",peptide);
                outRow.put("ratio",ratio);
                outRow.put("intensity1",intensity1);
                outRow.put("intensity2",intensity2);
                outRow.put("proteins",proteinsString);
                outRow.put("genes",genesString);

                outRows[i++] = outRow;
            }
            TabWriter tabWriter = new TabWriter(outColumns, outRows, outFile);
            tabWriter.write();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

    protected String assembleStringFromCollection(Collection<String> stringCollection)
    {
        if (stringCollection == null)
            return "";
        String[] stringCollectionArray =
                stringCollection.toArray(new String[stringCollection.size()]);
        StringBuffer resultBuf = new StringBuffer();
        for (int i=0; i<stringCollectionArray.length; i++)
        {
            if (i>0)
               resultBuf.append(",");
            resultBuf.append(stringCollectionArray[i]);
        }
        return resultBuf.toString();
    }



    protected void peptideGeneRollup() throws CommandLineModuleExecutionException
    {
        try
        {
            Map<String,List<String>> ipiGeneArrayMap =
                    GeneMappingUtilities.loadIPIGeneMap(geneLookupFile);

            TabLoader tl = new TabLoader(inTabFile);
            Map[] rowsAsMaps = (Map[]) tl.load();

            Set<String> allPeptides = new HashSet<String>();
            for (int i=0; i<rowsAsMaps.length; i++)
            {
                allPeptides.add((String) rowsAsMaps[i].get("peptide"));
            }

            Map<String,List<Protein>> peptideProteinMap =
                    ProteinUtilities.mapPeptidesToProteins(allPeptides, protXmlFiles, proteinsInFasta,0);

            Map<String,GeneMappingUtilities.InfoForGene> geneInfoMap =
                    new HashMap<String, GeneMappingUtilities.InfoForGene>();

            for (Map rowAsMap : rowsAsMaps)
            {
                String peptide = (String) rowAsMap.get("peptide");
                Double ratio = (Double) rowAsMap.get("ratio");
                Double intensity1 = (Double) rowAsMap.get("intensity1");
                Double intensity2 = (Double) rowAsMap.get("intensity2");


                if (peptideProteinMap.containsKey(peptide))
                {
                    for (Protein protein : peptideProteinMap.get(peptide))
                    {
                        String proteinLookup = protein.getLookup();
                        if (ipiGeneArrayMap.containsKey(proteinLookup))
                        {
                            List<String> genes = ipiGeneArrayMap.get(proteinLookup);

                            for (String gene : genes)
                            {
                                GeneMappingUtilities.InfoForGene infoForGene = geneInfoMap.get(gene);
                                if (infoForGene == null)
                                {
                                    infoForGene = new GeneMappingUtilities.InfoForGene(gene);
                                    geneInfoMap.put(gene, infoForGene);
                                }

                                //there should only be one ratio in this file per peptide, so if we've
                                //found it that means we've followed two different proteins for this gene,
                                //from the same peptide
                                if (infoForGene.getPeptides().contains(peptide))
                                    continue;

                                infoForGene.getPeptides().add(peptide);

                                infoForGene.getRatios().add(ratio);

                                infoForGene.getIntensities1().add(intensity1);
                                infoForGene.getIntensities2().add(intensity2);                                

                                infoForGene.getProteins().add(protein.getLookup());
                            }
                        }
                    }
                }
            }

            GeneMappingUtilities.InfoForGene.writeGeneRatioFile(geneInfoMap.values(), outFile);

            ApplicationContext.infoMessage(geneInfoMap.size() + " genes found in file");
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }





    protected void createPepXmlForTabFile()
            throws CommandLineModuleExecutionException
    {
        List<Feature> featuresForLines = new ArrayList<Feature>();

        BufferedReader br = null;

        try
        {

//            br = new BufferedReader(new FileReader(inTabFile));
//            String header = br.readLine();
                        TabLoader loader = new TabLoader(inTabFile, HashMap.class);

            HashMap[] rowsAsMaps = (HashMap[])loader.load();

            for (HashMap rowAsMap : rowsAsMaps)
            {
                Object ratioObject = rowAsMap.get("ratio");
                if (ratioObject == null)
                    ratioObject = rowAsMap.get("mean_ratio");

                double ratio = (Double) ratioObject;
                String peptide = (String) rowAsMap.get("peptide");

                Feature feature = new Feature();
                IsotopicLabelExtraInfoDef.setRatio(feature,ratio);
//                feature.setIntensity((float) ratio);
                IsotopicLabelExtraInfoDef.setHeavyIntensity(feature, ratio);
                IsotopicLabelExtraInfoDef.setLightIntensity(feature,1);
                MS2ExtraInfoDef.addPeptide(feature, peptide);
                feature.setCharge(2);

                featuresForLines.add(feature);
            }
            ProteinUtilities.createPepXml(featuresForLines.toArray(new Feature[featuresForLines.size()]), fastaFile, outFile);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (br != null)
                try {br.close();} catch(Exception e) {}
        }

    }

    /*
    protected void collapsePeptides()
            throws CommandLineModuleExecutionException

    {
        Map<String, List<Double>> proteinRatioListMap =
                new HashMap<String, List<Double>>();
        try
        {
            TabLoader loader = new TabLoader(inTabFile, HashMap.class);

            HashMap[] rowsAsMaps = (HashMap[])loader.load();

            ApplicationContext.infoMessage("Total lines: " + rowsAsMaps.length);
            for (HashMap rowMap :rowsAsMaps)
            {
                String protein = (String) rowMap.get("protein");
                List<Double> thisProteinRatioList = proteinRatioListMap.get(protein);
                if (thisProteinRatioList == null)
                {
                    thisProteinRatioList = new ArrayList<Double>();
                    proteinRatioListMap.put(protein, thisProteinRatioList);
                }

                thisProteinRatioList.add((Double) rowMap.get("ratio"));
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        PrintWriter outPW = null;

        int numProteinsPassing = 0;
        try
        {
            outPW = getPrintWriter(outFile);

            outPW.println("protein\tratio");
            for (String protein : proteinRatioListMap.keySet())
            {
                List<Double> ratios = proteinRatioListMap.get(protein);
                if (ratios.size() < minPeptides)
                    continue;
System.err.println("Protein: " + protein);                
                numProteinsPassing++;
                outPW.println(protein + "\t" + BasicStatistics.geometricMean(ratios));
                outPW.flush();
            }     
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (outPW != null)
                outPW.close();
        }
        ApplicationContext.infoMessage("Wrote out " + numProteinsPassing + " proteins");
    }
    */


    /**
     * roll up tab file to protein level naively (and natively)
     */
    protected void rollupTabFile()
            throws CommandLineModuleExecutionException
    {
        Map<Protein, Integer> proteinNumPeptidesMap = new HashMap<Protein,Integer>();
        PrintWriter outPW = null;
        try
        {

            outPW = getPrintWriter(outFile);
            TabLoader loader = new TabLoader(inTabFile, HashMap.class);

            StringBuffer firstLine = new StringBuffer();
            TabLoader.ColumnDescriptor[] columns = loader.getColumns();
            for (int i=0; i<columns.length; i++)
            {
                firstLine.append(columns[i].name);
                if (i<(columns.length-1))
                    firstLine.append("\t");
            }
            outPW.println("protein\t" + firstLine.toString());


            HashMap[] rowsAsMaps = (HashMap[])loader.load();

            ApplicationContext.infoMessage("Total lines: " + rowsAsMaps.length);

            Set<String> allPeptides = new HashSet<String>();
            for (int i=0; i<rowsAsMaps.length; i++)
            {
                allPeptides.add((String) rowsAsMaps[i].get("peptide"));
            }

            Map<String,List<Protein>> peptideProteinMap =
                    ProteinUtilities.mapPeptidesToProteins(allPeptides, protXmlFiles, proteinsInFasta,0);

            for (int i=0; i<rowsAsMaps.length; i++)
            {
                HashMap rowAsMap = rowsAsMaps[i];

                String peptide = (String) rowAsMap.get("peptide");

                List<Protein> proteinsThisPeptide = peptideProteinMap.get(peptide);
                if (proteinsThisPeptide == null)
                    continue;

                for (Protein protein : proteinsThisPeptide)
                {
                    outPW.print(protein.getLookup() + "\t");
                    for (int j=0; j<columns.length; j++)
                    {
                        outPW.print(rowAsMap.get(columns[j].name));
                        if (j<(columns.length-1))
                            outPW.print("\t");
                        else
                            outPW.print("\n");
                    }
                    outPW.flush();
                    if (proteinNumPeptidesMap.containsKey(protein))
                        proteinNumPeptidesMap.put(protein, proteinNumPeptidesMap.get(protein)+1);
                    else
                        proteinNumPeptidesMap.put(protein, 1);

                }
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (outPW != null)
                outPW.close();
        }

        double[] peptideCounts = new double[proteinNumPeptidesMap.size()];
        int i=0;
        for (Integer num : proteinNumPeptidesMap.values())
            peptideCounts[i++] = num;
        PanelWithHistogram pwh = new PanelWithHistogram(peptideCounts);
        pwh.setVisible(true);
    }



//
//            /**
//            * roll up to protein level naively (and natively)
//            */
//    protected void nativeRollup(Feature[] ms1FeaturesWithPeptides, File outputFile)
//            throws CommandLineModuleExecutionException
//    {
//        Set<String> alreadyUsedPeptides = new HashSet<String>();
//        int numAmbiguousPeptides = 0;
//        List<AmtMatchedProtein> matchedProteinList = new ArrayList<AmtMatchedProtein>();
//
//        int i=0;
//        for (Protein protein : proteinsInFasta)
//        {
//            if (i++%500 == 0)
//                ApplicationContext.setMessage(i + " proteins checked (out of " + proteinsInFasta.length + ")");
//            AmtMatchedProtein matchedProtein = new AmtMatchedProtein(protein);
//
//            int numMatchesThisProtein=0;
//            int numUniqueMatchesThisProtein = 0;
//            String proteinSequence = protein.getSequenceAsString();
//            for (Feature feature : ms1FeaturesWithPeptides)
//            {
//                String peptideSequence = MS2ExtraInfoDef.getFirstPeptide(feature);
//                if (proteinSequence.contains(peptideSequence))
//                {
//                    boolean unique = true;
//                    numMatchesThisProtein++;
//                    if (alreadyUsedPeptides.contains(peptideSequence))
//                    {
//                        numAmbiguousPeptides++;
//                        unique = false;
//                    }
//                    else
//                        numUniqueMatchesThisProtein++;
//                    alreadyUsedPeptides.add(peptideSequence);
//
//                    AmtMatchedProtein.AmtMatchedPeptide matchedPeptide =
//                            new AmtMatchedProtein.AmtMatchedPeptide(peptideSequence, feature.getIntensity(), unique);
//                    matchedProtein.addPeptideInfo(matchedPeptide);
//                }
//
//            }
//            if (numMatchesThisProtein > 0)
//            {
//                matchedProtein.setNumMatches(numMatchesThisProtein);
//                matchedProtein.setNumUniqueMatches(numUniqueMatchesThisProtein);
//                matchedProteinList.add(matchedProtein);
//            }
//        }
//
//
//
//        try
//        {
//            PrintWriter outPW = getPrintWriter(outputFile);
//            outPW.println(AmtMatchedProtein.getFileHeaderString());
//            for (AmtMatchedProtein matchedProtein : matchedProteinList)
//            {
//                outPW.print(matchedProtein.toString());
//            }
//            outPW.close();
//            if (outputFile != null)
//                ApplicationContext.setMessage("Wrote output file " +
//                        outputFile.getAbsolutePath());
//        }
//        catch (FileNotFoundException e)
//        {
//            throw new CommandLineModuleExecutionException(e);
//        }
//    }

}
