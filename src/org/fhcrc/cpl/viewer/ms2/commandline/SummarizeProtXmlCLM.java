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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBarChart;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class SummarizeProtXmlCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SummarizeProtXmlCLM.class);

    protected File protXmlFiles[];
    protected boolean showCharts = false;
    protected boolean listProteins = false;
    protected float minProteinProphet = 0.1f;
    protected File protGeneFile;
    protected Map<String,List<String>> protGenesMap = null;
	
	//map from proteins to groups, loaded from a groupified MasterGroup file
    protected Map<String, Integer> proteinGroupMap = null;

    protected File fastaFile;
    protected File pepXmlFile;

    protected File outProteinListFile = null;
    protected File outProteinPeptideFile = null;

    protected boolean shouldBarChartOrganism = false;

    int nMostAbundantToList = 20;

    public SummarizeProtXmlCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "summarizeprotxml";
        mShortDescription = "Summarize the contents of one or more protXML files";
        mHelpMessage = "Summarize the contents of one or more protXML files";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "protxml file(s)"),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new DecimalArgumentDefinition("minpprophet", false,
                                "Min proteinprophet for MA plot and protein list",
                                minProteinProphet),
                        new FileToReadArgumentDefinition("protgenefile", false,
                                "File associating gene symbols with protein accession numbers"),
                        new FileToReadArgumentDefinition("mastergroupfile", false,"'Groupified' MasterGroup file to use for protein grouping instead of protein groups in protxml file"),
                        new BooleanArgumentDefinition("listproteins", false, "List proteins to stderr?",
                                listProteins),
                        new BooleanArgumentDefinition("organism", false, "bar chart of organism? (only appropriate for " +
                                "SwissProt searches)",
                                shouldBarChartOrganism),
                        new FileToReadArgumentDefinition("fasta", false,
                                "FASTA file for examining protein sequences"),
                        new FileToReadArgumentDefinition("pepxml", false, "PepXML file"),
                        new FileToWriteArgumentDefinition("outproteinfile", false,
                                "Output file containing all proteins from this file, one per line"),
                        new FileToWriteArgumentDefinition("outproteinpeptidefile", false,
                                "Output file containing all protein and peptide information, one line per unique " +
                                        "peptide sequence"),
                        new IntegerArgumentDefinition("numabundanttolist", false,
                                "Number of the most-abundant proteins to list with counts", nMostAbundantToList),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();
        showCharts = getBooleanArgumentValue("showcharts");
        listProteins = getBooleanArgumentValue("listproteins");
        minProteinProphet = getFloatArgumentValue("minpprophet");
        protGeneFile = getFileArgumentValue("protgenefile");
        nMostAbundantToList = getIntegerArgumentValue("numabundanttolist");
        File masterGroupFile = getFileArgumentValue("mastergroupfile");
        if (masterGroupFile != null)
	{
 		proteinGroupMap = new HashMap<String, Integer>();
		try
		{
			Map[] rowMaps = TabLoader.loadMaps(masterGroupFile);
System.err.println("Loaded rows from mastergroup file");
			String groupColName = null;
			for (Object keyObj : rowMaps[0].keySet())
			{
				String key = keyObj.toString();
				if (key.startsWith("group"))
				{
					groupColName = key;
					break;
				}
			}

			for (Map rowMap : rowMaps)
			{
				if (rowMap.get(groupColName) == null)
					continue;
				int groupNumber = (Integer) rowMap.get(groupColName);
				String[] proteins = ((String)rowMap.get("protein")).split(";");
				for (String protein : proteins)
				{
					proteinGroupMap.put(protein, groupNumber);
				}
			}
		}
		catch (Exception e)
		{
			throw new ArgumentValidationException("Failed to load mastergroup file",e);
		}
 	}

        shouldBarChartOrganism = getBooleanArgumentValue("organism");
        fastaFile = getFileArgumentValue("fasta");
        pepXmlFile = getFileArgumentValue("pepxml");
        outProteinListFile = getFileArgumentValue("outproteinfile");
        outProteinPeptideFile = getFileArgumentValue("outproteinpeptidefile");
        if (outProteinPeptideFile != null)
            assertArgumentPresent("pepxml","outproteinpeptidefile");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        System.out.println("File\tGroups\tPoint1\tPoint5\tPoint75\tPoint9\tPoint95");
        if (protGeneFile != null)
        {
            try
            {
                ApplicationContext.infoMessage("Loading IPI-gene map...");
                protGenesMap = QAUtilities.loadIpiGeneListMap(protGeneFile);
                ApplicationContext.infoMessage("Done");

            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to load protein-gene map file",e);
            }
        }
        for (File protXmlFile : protXmlFiles)
            summarizeFile(protXmlFile);
    }

    protected void summarizeFile(File protXmlFile) throws CommandLineModuleExecutionException
    {
        try
        {
            if (showCharts)
            {
                PanelWithChart sensSpecChart =
                        ProteinUtilities.generateSensSpecChart(protXmlFile);
                sensSpecChart.setName("Sens/Spec");
                sensSpecChart.displayInTab();
            }

            ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

            List<Float> groupProbabilityList = new ArrayList<Float>();

            List<Float> proteinProbabilityList = new ArrayList<Float>();
            List<Float> proteinSpectralCountList = new ArrayList<Float>();
            List<Float> proteinRatioForMAList = new ArrayList<Float>();

            List<Float> proteinRatioCOVList = new ArrayList<Float>();
            List<Float> proteinWithRatioProbabilityList = new ArrayList<Float>();
            List<Float> proteinNumPeptidesList = new ArrayList<Float>();

            List<Float> proteinSpectralCountForMAList = new ArrayList<Float>();
            List<Float> proteinRatioNumPeptidesList = new ArrayList<Float>();
            List<Float> proteinRatioNumQuantPeptidesList = new ArrayList<Float>();


            //can only do these if we've got pepXML file
            List<Float> proteinRatioNumQuantEventsList = new ArrayList<Float>();
            List<Float> proteinRatioNumUniqueQuantPeptidesList = new ArrayList<Float>();

            Map<String, Integer> proteinSpectralCountMap = new HashMap<String, Integer>();

            List<Float> proteinCoveragePercentList = new ArrayList<Float>();

            Set<Integer> quantGroups = new HashSet<Integer>();


            Iterator<ProteinGroup> iterator = protXmlReader.iterator();

            Set<Integer> groupsWithProbabilityGreaterThanPoint1 = new HashSet<Integer>();
            Set<Integer> groupsWith2PlusPeptides = new HashSet<Integer>();
            Set<Integer> groupsWithZeroOrOneGenes = new HashSet<Integer>();

//            Map<String, Set<String>> peptideProteinMap = new HashMap<String, Set<String>>();
            Map<String, Set<String>> quantPeptideProteinMap = new HashMap<String, Set<String>>();

            //can only populate this if fastaFile != null
            List<Integer> proteinLengths = new ArrayList<Integer>();
            List<Double> proteinMasses = new ArrayList<Double>();
            Map<String, Protein> proteinNameProteinMap = null;
            if (fastaFile != null)
                proteinNameProteinMap = ProteinUtilities.loadProteinNameProteinMapFromFasta(fastaFile);

            Map<String, Integer> peptideQuantCountMap = null;
            Map<String, List<Feature>> peptideFeaturesMap = null;
            if (pepXmlFile != null)
            {
                List<FeatureSet> featureSets = PepXMLFeatureFileHandler.getSingletonInstance().loadAllFeatureSets(pepXmlFile);
                peptideQuantCountMap = createPeptideQuantEventCountMap(featureSets);

                peptideFeaturesMap = new HashMap<String, List<Feature>>();
                for (FeatureSet featureSet : featureSets)
                {
                    for (Feature feature : featureSet.getFeatures())
                    {
                        String peptide = MS2ExtraInfoDef.getPeptideList(feature).get(0);
                        List<Feature> featuresThisPeptide = peptideFeaturesMap.get(peptide);
                        if (featuresThisPeptide == null)
                        {
                            featuresThisPeptide = new ArrayList<Feature>();
                            peptideFeaturesMap.put(peptide, featuresThisPeptide);
                        }
                        featuresThisPeptide.add(feature);
                    }
                }
            }

            PrintWriter outProteinPeptidePW = null;
            if (outProteinPeptideFile != null)
            {
                try
                {
                    outProteinPeptidePW = new PrintWriter(outProteinPeptideFile);
                    outProteinPeptidePW.println(
                            "group\tprotein\tgene\tdesc\tprotein_probability\ttotal_peptides\tunique_peptides\tcoverage\tpeptide\tpeptide_probability\tmodifications\tmass_theoretical\tmass_difference\tcharge");
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Failed to write output file " +
                            outProteinPeptideFile.getAbsolutePath());
                }
            }

            Map<String, Float> organismCountMap = new HashMap<String, Float>();

            Set<String> allProteins = new HashSet<String>();

            while (iterator.hasNext())
            {
                ProteinGroup pg = iterator.next();
                Set<String> genesThisGroup = new HashSet<String>();

                Set<String> allPeptidesThisGroup = new HashSet<String>();

                groupProbabilityList.add(pg.getProbability());


                for (ProtXmlReader.Protein protein : pg.getProteins())
                {
                    String proteinName = protein.getProteinName();
                    String proteinPeptideFileProteinBits = null;
                    if (outProteinPeptideFile != null && protein.getProbability() >= minProteinProphet)
                    {
                        String description = "";
                        if (proteinNameProteinMap != null)
                        {
                            description = proteinNameProteinMap.get(protein.getProteinName()).getHeader();
                            int descriptionIndex = description.indexOf("Gene_Symbol=") + "Gene_Symbol=".length();
                            while (!Character.isSpaceChar(description.charAt(descriptionIndex)))
                                descriptionIndex++;
                            descriptionIndex++;
                            description = description.substring(descriptionIndex);
                        }
                        String gene = "";
                        if (protGenesMap != null && protGenesMap.containsKey(protein.getProteinName()))
                            gene = MS2ExtraInfoDef.convertStringListToString(protGenesMap.get(protein.getProteinName()));
                        int group = pg.getGroupNumber();
                        if (proteinGroupMap != null)
                        {
                            if (proteinGroupMap.containsKey(proteinName))
                            {
                                group = proteinGroupMap.get(proteinName);
                            }
                            else
                            {
                                boolean foundIt = false;
                                for (String name : protein.getIndistinguishableProteinNames())
                                    if (proteinGroupMap.containsKey(name))
                                    {
                                        group = proteinGroupMap.get(name);
                                        proteinName = name;
                                        foundIt = true;
                                    }
                                if (!foundIt)
                                    continue;
                            }


                        }
                        proteinPeptideFileProteinBits = group + ""
                                + "\t" + proteinName + "\t" + gene + "\t" + description + "\t" + protein.getProbability() + "\t"
                                + protein.getTotalNumberPeptides() + "\t" +
                                protein.getUniquePeptidesCount() + "\t" + (protein.getPercentCoverage()) + "\t";
                    }
                    if (protein.getProbability() >= minProteinProphet)
                    {
                        allProteins.add(protein.getProteinName());
                    }
                    if (proteinNameProteinMap != null)
                    {
                        try
                        {
//System.err.println("****" + proteinNameProteinMap.get(protein.getProteinName()));
                            Protein fastaProtein = proteinNameProteinMap.get(protein.getProteinName());
                            proteinLengths.add(fastaProtein.getSequenceAsString().length());
                            proteinMasses.add(fastaProtein.getMass());
                        }
                        catch(NullPointerException e)
                        {
                            ApplicationContext.infoMessage("Missing protein in FASTA: " + protein.getProteinName());
                        }
                    }
                    if (shouldBarChartOrganism)
                    {
                        String organism = protein.getProteinName().substring(protein.getProteinName().indexOf("_")+1);
                        if (!organismCountMap.containsKey(organism))
                            organismCountMap.put(organism, 0f);
                        organismCountMap.put(organism, organismCountMap.get(organism)+1);
                    }
                    proteinProbabilityList.add(protein.getProbability());
                    int specCount = 0;
                    proteinNumPeptidesList.add((float) protein.getUniquePeptidesCount());
                    if (protein.getPercentCoverage() != null)
                        proteinCoveragePercentList.add(protein.getPercentCoverage());
                    for (ProtXmlReader.Peptide peptide : protein.getPeptides())
                    {
                        allPeptidesThisGroup.add(peptide.getPeptideSequence());
                        specCount += peptide.getInstances();


                        if (outProteinPeptideFile != null && protein.getProbability() >= minProteinProphet && proteinPeptideFileProteinBits != null)
                        {
                            StringBuffer lineBuf = new StringBuffer(proteinPeptideFileProteinBits);
                            Feature bestProbFeature = null;
                            for (Feature feature : peptideFeaturesMap.get(peptide.getPeptideSequence()))
                            {
                                if (bestProbFeature == null || MS2ExtraInfoDef.getPeptideProphet(feature) >
                                        MS2ExtraInfoDef.getPeptideProphet(bestProbFeature))
                                    bestProbFeature = feature;
                            }



                            lineBuf.append(peptide.getPeptideSequence() + "\t" + MS2ExtraInfoDef.getPeptideProphet(bestProbFeature) + "\t" +
                                    MS2ExtraInfoDef.createModifiedSequenceString(peptide.getPeptideSequence(), MS2ExtraInfoDef.getModifiedAminoAcidsMap(bestProbFeature)) + "\t" +
                                    bestProbFeature.getMass() + "\t" + MS2ExtraInfoDef.getDeltaMass(bestProbFeature) + "\t" + bestProbFeature.getCharge());
                            outProteinPeptidePW.println(lineBuf.toString());
                        }
                    }
                    proteinSpectralCountList.add((float) specCount);
                    proteinSpectralCountMap.put(protein.getProteinName(), specCount);

                    if (protein.getProbability() >= minProteinProphet && protein.getQuantitationRatio() != null)
                    {
                        proteinRatioForMAList.add(protein.getQuantitationRatio().getRatioMean());
                        proteinSpectralCountForMAList.add((float) protein.getTotalNumberPeptides());
                        proteinRatioCOVList.add(protein.getQuantitationRatio().getRatioStandardDev() / protein.getQuantitationRatio().getRatioMean());
                        proteinRatioNumPeptidesList.add((float) protein.getQuantitationRatio().getPeptides().size());
                        proteinRatioNumQuantPeptidesList.add((float) protein.getQuantitationRatio().getRatioNumberPeptides());
                        int numQuantEvents = 0;
                        int numUniqueQuantPeptides = 0;
                        for (String peptide : protein.getQuantitationRatio().getPeptides())
                        {
                            if (peptideQuantCountMap != null && peptideQuantCountMap.containsKey(peptide))
                            {
                                numQuantEvents += peptideQuantCountMap.get(peptide);
                                numUniqueQuantPeptides++;
                            }
                        }
                        proteinRatioNumQuantEventsList.add((float) numQuantEvents);
                        proteinRatioNumUniqueQuantPeptidesList.add((float)numUniqueQuantPeptides);
                        proteinWithRatioProbabilityList.add(pg.getProbability());
                        quantGroups.add(pg.getGroupNumber());
                    }
                    List<String> genesThisProtein = null;
                    if (protGenesMap != null && protGenesMap.containsKey(protein.getProteinName()))
                    {
                        genesThisProtein = protGenesMap.get(protein.getProteinName());
                        genesThisGroup.addAll(genesThisProtein);
                    }
                    if (listProteins)
                    {
                        String genesString = "";
                        if (genesThisProtein != null)
                        {
                            StringBuffer genesStringBuf = new StringBuffer();
                            for (int i=0; i<genesThisProtein.size(); i++)
                            {
                                if (i>0) genesStringBuf.append(",");
                                genesStringBuf.append(genesThisProtein.get(i));
                            }
                            genesString = genesStringBuf.toString();
                        }
                        System.err.println(protein.getProteinName() + "\t" + genesString + "\t" +
                                protein.getProbability() + "\t" +
                                protein.getPeptides().size());
                    }
                }
                if (allPeptidesThisGroup.size() > 1)
                    groupsWith2PlusPeptides.add(pg.getGroupNumber());
                if (pg.getGroupProbability() >= minProteinProphet)
                    groupsWithProbabilityGreaterThanPoint1.add(pg.getGroupNumber());
                if (protGenesMap != null && genesThisGroup.size() < 2)
                    groupsWithZeroOrOneGenes.add(pg.getGroupNumber());
            }
            if (outProteinPeptidePW != null)
            {
                outProteinPeptidePW.close();
                ApplicationContext.infoMessage("Saved protein-peptide file to " + outProteinPeptideFile.getAbsolutePath());
            }

            ApplicationContext.infoMessage(nMostAbundantToList + " most abundant proteins, with counts: ");
            List<String> proteinsByCountDesc = new ArrayList<String>(proteinSpectralCountMap.keySet());
            Collections.sort(proteinsByCountDesc, new MapBasedComparatorDesc(proteinSpectralCountMap));
            for (int i=0; i<Math.min(nMostAbundantToList, proteinsByCountDesc.size()); i++)
            {
                String prot = proteinsByCountDesc.get(i);
                String genes = "?";
                if (protGenesMap != null && protGenesMap.containsKey(prot))
                {
                    StringBuffer genesBuf = new StringBuffer();
                    for (int j=0; j<protGenesMap.get(prot).size(); j++)
                    {
                        if (j>0) genesBuf.append(",");
                        genesBuf.append(protGenesMap.get(prot).get(j));
                    }
                    genes = genesBuf.toString();
                }
                ApplicationContext.infoMessage(prot + "\t" + genes + "\t" + proteinSpectralCountMap.get(prot));
            }


            ApplicationContext.infoMessage("Protein Groups with probability >=" + minProteinProphet + ": " +
                    groupsWithProbabilityGreaterThanPoint1.size());
            ApplicationContext.infoMessage("Single-peptide Protein Groups with probability >=" + minProteinProphet + ": " +
                    (groupsWithProbabilityGreaterThanPoint1.size() - setIntersection(groupsWith2PlusPeptides, groupsWithProbabilityGreaterThanPoint1).size()));
            ApplicationContext.infoMessage("Protein Groups with probability >=" + minProteinProphet + " and 2 or more peptides: " +
                    setIntersection(groupsWith2PlusPeptides, groupsWithProbabilityGreaterThanPoint1).size());
            if (protGenesMap != null)
            {
                ApplicationContext.infoMessage("Protein Groups with probability >=" + minProteinProphet + " and mapping to zero or one genes: " +
                        setIntersection(groupsWithZeroOrOneGenes, groupsWithProbabilityGreaterThanPoint1).size());
                ApplicationContext.infoMessage("Protein Groups with probability >=" + minProteinProphet + " and mapping to zero or one genes " +
                        " and 2 or more peptides: " +
                        setIntersection(groupsWith2PlusPeptides, setIntersection(groupsWithZeroOrOneGenes, groupsWithProbabilityGreaterThanPoint1)).size());                
            }

            ApplicationContext.infoMessage("Median spectral count per protein: " + BasicStatistics.median(proteinSpectralCountList) + ", mean: " + BasicStatistics.mean(proteinSpectralCountList));
            ApplicationContext.infoMessage("Median Unique Peptides per protein: " + BasicStatistics.median(proteinNumPeptidesList) + ", mean: " + BasicStatistics.mean(proteinNumPeptidesList));
            ApplicationContext.infoMessage("Median AA coverage: " + BasicStatistics.median(proteinCoveragePercentList));

            if (!proteinLengths.isEmpty())
                ApplicationContext.infoMessage("Median protein sequence length: " + BasicStatistics.median(proteinLengths) + ", mean: " + BasicStatistics.mean(proteinLengths));
            if (!proteinMasses.isEmpty())
                ApplicationContext.infoMessage("Median protein mass: " + BasicStatistics.median(proteinMasses) + ", mean: " + BasicStatistics.mean(proteinMasses));


            List<Float> logRatios = new ArrayList<Float>();
            List<Float> logSpectralCounts = new ArrayList<Float>();
            if (proteinRatioForMAList.size() > 0)
            {
                for (Float ratio : proteinRatioForMAList)
                    logRatios.add((float) Math.log(Math.max(ratio,0.000001)));
                for (Float count : proteinSpectralCountForMAList)
                    logSpectralCounts.add((float) Math.log(Math.max(1, count)));
                ApplicationContext.infoMessage("Median log ratio: " + BasicStatistics.median(logRatios));
                ApplicationContext.infoMessage("Median ratio: " + Math.exp(BasicStatistics.median(logRatios)));
                if (pepXmlFile != null)
                {
                ApplicationContext.infoMessage("Median unique quantitated peptides per protein: " +
                        BasicStatistics.median(proteinRatioNumUniqueQuantPeptidesList) + ", mean: " + BasicStatistics.mean(proteinRatioNumUniqueQuantPeptidesList));
                    ApplicationContext.infoMessage("Median quantification events per protein: " +
                            BasicStatistics.median(proteinRatioNumQuantEventsList) + ", mean: " + BasicStatistics.mean(proteinRatioNumQuantEventsList));

                }
            }


            if (showCharts)
            {
                if (shouldBarChartOrganism)
                {
                    Map<String, Float> abundantOrganismMap = new HashMap<String, Float>();
                    float totalProteinCount = 0;
                    for (String organism : organismCountMap.keySet())
                    {
                        float thisOrgCount = (float) organismCountMap.get(organism);
                        totalProteinCount += thisOrgCount;
                    }
                    for (String organism : organismCountMap.keySet())
                    {
                        float thisOrgCount = (float) organismCountMap.get(organism);
                        if (thisOrgCount >= totalProteinCount / 50f)
                        {
                            ApplicationContext.infoMessage("ORGANISM: " + organism + " with " + (thisOrgCount * 100 / totalProteinCount) + "%");
abundantOrganismMap.put(organism, thisOrgCount);
                        }
                    }
//                    ApplicationContext.infoMessage("ORGANISM: " + (humanCount * 100 / totalProteinCount) + "% Human");
//                    ApplicationContext.infoMessage("Highest Nonhuman: " + highestNonhuman + " with " + (highestNonhumanCount * 100 / totalProteinCount) + "%");

                    PanelWithBarChart pwbc = new PanelWithBarChart(abundantOrganismMap, "Abundant Organisms");
                    pwbc.displayInTab();
                }
                PanelWithHistogram pwh = new PanelWithHistogram(groupProbabilityList, "Probabilities of Proteins");
                pwh.displayInTab();

                if (!proteinCoveragePercentList.isEmpty())
                    new PanelWithHistogram(proteinCoveragePercentList, "% AA Coverage").displayInTab();

                if (!proteinMasses.isEmpty())
                    new PanelWithHistogram(proteinMasses, "Prot Masses").displayInTab();



                new PanelWithHistogram(proteinNumPeptidesList, "Unique Peptides", 200).displayInTab();


                if (proteinRatioForMAList.size() > 0)
                {

                    PanelWithScatterPlot pwsp = new PanelWithScatterPlot(logSpectralCounts,
                            logRatios, "MAPlot");
                    pwsp.displayInTab();



                    PanelWithHistogram pwhRatios = new PanelWithHistogram(logRatios, "Log Ratios");
                    pwhRatios.displayInTab();


                    PanelWithScatterPlot pwsp3 = new PanelWithScatterPlot(proteinWithRatioProbabilityList, logRatios, 
                            "Probability vs LogRatio");
                    pwsp3.setAxisLabels("Probability","Log Ratio");
                    pwsp3.displayInTab();

                    PanelWithScatterPlot pwsp4 = new PanelWithScatterPlot(proteinRatioNumQuantPeptidesList, proteinRatioCOVList,
                            "Quant Events vs COV");
                    pwsp4.setAxisLabels("Quant Events","COV");
                    pwsp4.displayInTab();

                    PanelWithScatterPlot pwsp5 = new PanelWithScatterPlot(proteinRatioNumPeptidesList, logRatios, 
                            "Quant Events vs LogRatio");
                    pwsp5.setAxisLabels("Quant Events","Log Ratio");
                    pwsp5.displayInTab();

                    if (pepXmlFile != null)
                    {
                    PanelWithHistogram pwhEvents = new PanelWithHistogram(proteinRatioNumQuantEventsList, "Quant Events");
                    pwhEvents.displayInTab();

                    new PanelWithHistogram(proteinRatioNumUniqueQuantPeptidesList, "Unique Quant Peptides").displayInTab();
                    }
                }
            }

            if (proteinRatioForMAList.size() > 0)
            {
                int numWithin25Percent=0;
                int numWithin50Percent=0;
                int numWithinTwofold=0;
                int numWithinThreefold=0;
                for (Float ratio : proteinRatioForMAList)
                {
                    Float ratioBigOverSmall = (float) Math.exp(Math.abs(Math.log(ratio)));
                    if (ratioBigOverSmall <= 3)
                    {
                        numWithinThreefold++;
                        if (ratioBigOverSmall <= 2)
                        {
                            numWithinTwofold++;
                            if (ratioBigOverSmall <= 1.5)
                            {
                                numWithin50Percent++;
                                if (ratioBigOverSmall <= 1.25)
                                    numWithin25Percent++;
                            }
                        }
                    }
                }
                int numRatios = proteinRatioForMAList.size();
                int percentWithin25Percent=100 * numWithin25Percent / numRatios;
                int percentWithin50Percent=100 * numWithin50Percent / numRatios;
                int percentWithinTwofold=100 * numWithinTwofold / numRatios;
                int percentWithinThreefold=100 * numWithinThreefold / numRatios;
                ApplicationContext.infoMessage("Quantitated proteins: " + numRatios +  "(" + 
                        (100f * numRatios / proteinProbabilityList.size()) + "%), groups: " + quantGroups.size() +
                        " (" + (100f * quantGroups.size() / groupsWithProbabilityGreaterThanPoint1.size()) + "%)");
                ApplicationContext.infoMessage("Ratios within 3fold: " + percentWithinThreefold + "%, 2fold: " +
                        numWithinTwofold + ", i.e., " + percentWithinTwofold + "%, 50%: " + percentWithin50Percent + "%, 25%: " + percentWithin25Percent + "%");

            }

            int numAbovePoint1 = 0;
            int numAbovePoint5 = 0;
            int numAbovePoint75 = 0;
            int numAbovePoint9 = 0;
            int numAbovePoint95 = 0;

            for (int i=0; i<groupProbabilityList.size(); i++)
            {
                float probability = groupProbabilityList.get(i);
                if (probability >= .5f)
                {
                    numAbovePoint1++;
                if (probability >= .5f)
                {
                    numAbovePoint5++;
                    if (probability >= .75f)
                    {
                        numAbovePoint75++;
                        if (probability >= .9f)
                        {
                            numAbovePoint9++;
                            if (probability >= .95f)
                            {
                                numAbovePoint95++;
                            }
                        }
                    }
                }
                }
            }

            List<Float> spectralCountsAbovePoint1 = new ArrayList<Float>();

            for (int i=0; i<proteinProbabilityList.size(); i++)
            {
                float probability = proteinProbabilityList.get(i);
                if (probability >= 0.1f)
                    spectralCountsAbovePoint1.add(proteinSpectralCountList.get(i));

            }

            System.out.println(protXmlFile.getName() + "\t" +  groupProbabilityList.size()+ "\t" + numAbovePoint1 + "\t" + numAbovePoint5 + "\t" +
                                        numAbovePoint75 + "\t" + numAbovePoint9 + "\t" + numAbovePoint95);

            if (showCharts)
            {
                List<Float> logSpectralCountsAbovePoint1 = new ArrayList<Float>(spectralCountsAbovePoint1.size());
                for (Float specCount : spectralCountsAbovePoint1)
                     logSpectralCountsAbovePoint1.add((float) Math.log(specCount));
                PanelWithHistogram pwhSpecCount = new PanelWithHistogram(logSpectralCountsAbovePoint1, "Log spec counts prob .1", 200);
                pwhSpecCount.displayInTab();
            }

            if (outProteinListFile != null)
            {
                PrintWriter outPW = null;
                try
                {
                    outPW = new PrintWriter(outProteinListFile);
                    for (String protein : allProteins)
                    {
                        outPW.println(protein);
                        outPW.flush();
                    }
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Failed to write output file " + outProteinListFile.getAbsolutePath());
                }
                finally
                {
                    if (outPW != null)
                        outPW.close();
                }
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

    public static <T> Set<T> setIntersection(Set<? extends T> a, Set<? extends T> b)
    {
        Set result = new HashSet<T>(a);
        result.retainAll(b);
        return result;
    }

    protected Map<String, Integer> createPeptideQuantEventCountMap(List<FeatureSet> featureSets)
    {
        Map<String, Integer> peptideEventCountMap = new HashMap<String, Integer>();

        for (int fraction=0; fraction<featureSets.size(); fraction++)
        {
            FeatureSet featureSet = featureSets.get(fraction);
            for (Feature feature : featureSet.getFeatures())
            {
                if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                    continue;
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (!peptideEventCountMap.containsKey(peptide))
                    peptideEventCountMap.put(peptide, 0);
                peptideEventCountMap.put(peptide, peptideEventCountMap.get(peptide)+1);
//
//                Map<Integer, Map<Integer, Map<String, Integer>>> fractionChargesMap =
//                        peptideFractionChargeModsQuantEventsMap.get(peptide);
//                if (fractionChargesMap == null)
//                {
//                    fractionChargesMap = new HashMap<Integer, Map<Integer, Map<String, Integer>>>();
//                    peptideFractionChargeModsQuantEventsMap.put(peptide, fractionChargesMap);
//                }
//
//
//                Map<Integer, Map<String, Integer>> chargeModificationsMap =
//                    fractionChargesMap.get(fraction);
//                if (chargeModificationsMap == null)
//                {
//                    chargeModificationsMap = new HashMap<Integer, Map<String, Integer>>();
//                    fractionChargesMap.put(fraction, chargeModificationsMap);
//                }
//
//                Map<String, Integer> modificationsEventsMap =
//                        chargeModificationsMap.get(feature.getCharge());
//                if (modificationsEventsMap == null)
//                {
//                    modificationsEventsMap = new HashMap<String, Integer>();
//                    chargeModificationsMap.put(feature.getCharge(), modificationsEventsMap);
//                }
//
//                String modState =  MS2ExtraInfoDef.convertModifiedAminoAcidsMapToString(MS2ExtraInfoDef.getModifiedAminoAcidsMap(feature));
//
//                if (!modificationsEventsMap.containsKey(modState))
//                    modificationsEventsMap.put(modState, 0);
//                modificationsEventsMap.put(modState, modificationsEventsMap.get(modState) + 1);
            }
        }
        return peptideEventCountMap;
    }

    /**
     * A comparator for comparing the keys of a String->Number map based on their associated values.
     * I.e., sort these strings by these associated numbers.  E.g., sort proteins by spectral count.
     *
     * Descending order
     */
    class MapBasedComparatorDesc implements Comparator
    {
        private Map<String, ? extends Number> map;

        public MapBasedComparatorDesc(Map<String, ? extends Number> map)
        {
            this.map = map;
        }

        public int compare(Object o1, Object o2)
        {
            double diff =  map.get(o2).doubleValue() - map.get(o1).doubleValue();
            return (diff > 0) ? 1 : ((diff == 0) ? 0 : -1);
        }
    }

}
