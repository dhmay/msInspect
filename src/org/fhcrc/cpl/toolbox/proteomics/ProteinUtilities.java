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
package org.fhcrc.cpl.toolbox.proteomics;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeaturePepXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLStreamReader;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.FastaLoader;
import org.apache.log4j.Logger;


import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;
import java.util.List;
import java.awt.*;


/**
 * This class doesn't do parsing of ProtXML files itself.  It uses ProtXmlReader for that.
 *
 * This is for utilities to work with the output of ProtXmlReader
 */
public class ProteinUtilities
{
    protected static Logger _log = Logger.getLogger(ProteinUtilities.class);


    public static Map<String, Set<String>> loadPeptideProteinMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        return loadPeptideProteinMapFromProtXML(protXmlFile, minProteinProphet, false);
    }

    /**
     * Create a mapping between all peptides noted in the protXml file, and all proteins
     * that they are associated with.  This means all indistinguishable proteins
     * @param protXmlFile
     * @return
     * @throws org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException
     */
    public static Map<String, Set<String>> loadPeptideProteinMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet,
                                                                            boolean quantProteinsOnly)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each peptide to all proteins it supports
        Map<String, Set<String>> peptideProteinMap = new HashMap<String, Set<String>>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        while (iterator.hasNext())
        {
            ProteinGroup group = iterator.next();
            if (minProteinProphet > 0 && group.getGroupProbability() < minProteinProphet)
                continue;
            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {
                if (quantProteinsOnly && protXmlReaderProtein.getQuantitationRatio() == null)
                    continue;
                Set<String> proteinNames = new HashSet<String>();
                //add first protein
                proteinNames.add(protXmlReaderProtein.getProteinName());
                //add indistinguishable proteins
                proteinNames.addAll(protXmlReaderProtein.getIndistinguishableProteinNames());

                for (ProtXmlReader.Peptide peptide : protXmlReaderProtein.getPeptides())
                {
                    String peptideSequence = peptide.getPeptideSequence();
                    Set<String> proteinsThisPeptide = peptideProteinMap.get(peptideSequence);
                    if (proteinsThisPeptide == null)
                    {
                        proteinsThisPeptide = new HashSet<String>();
                        peptideProteinMap.put(peptideSequence, proteinsThisPeptide);
                    }
                    proteinsThisPeptide.addAll(proteinNames);
                }
            }
        }
        return peptideProteinMap;
    }

    /**
     * Look inside a protXML file to find the source PepXML files
     */
    public static List<File> findSourcePepXMLFiles(File protXmlFile)
            throws FileNotFoundException, XMLStreamException
    {
        FileInputStream fIn = new FileInputStream(protXmlFile);
        SimpleXMLStreamReader parser = new SimpleXMLStreamReader(fIn);
        List<File> result = new ArrayList<File>();
        _log.debug("findSourcePepXMLFiles");
        while (parser.hasNext())
        {
            if (!parser.skipToStart("protein_summary_header"))
            {
                _log.debug("No protein_summary_header in this protXML file!");
                break;
            }
            _log.debug("findSourcePepXMLFiles: found protein_summary_header");
            String sourceFilesString = parser.getAttributeValue(null, "source_files");
            String[] sourceFilePathsArray = sourceFilesString.split(" ");
            _log.debug("findSourcePepXMLFiles: source_files='" + sourceFilesString + "'");
            for (String sourceFilePath : sourceFilePathsArray)
            {
                File pepXmlFile = new File(sourceFilePath);
                _log.debug("findSourcePepXMLFiles: source file " + sourceFilePath);
                if (!pepXmlFile.exists() || !pepXmlFile.canRead())
                    throw new FileNotFoundException("Can't read PepXML file " + pepXmlFile.getAbsolutePath());
                result.add(pepXmlFile);
            }
        }
        if (result.isEmpty())
            throw new FileNotFoundException("No PepXML file specified in ProtXML file " +
                    protXmlFile.getAbsolutePath());
        return result;
    }

    /**
     * Generate a sensitivity-and-specificity-curve chart
     * @param protXmlFile
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static PanelWithChart generateSensSpecChart(File protXmlFile)
            throws  FileNotFoundException, XMLStreamException
    {
        FileInputStream fIn = new FileInputStream(protXmlFile);
        SimpleXMLStreamReader parser = new SimpleXMLStreamReader(fIn);

        List<Double> minProbValues = new ArrayList<Double>();
        List<Double> sensValues = new ArrayList<Double>();
        List<Double> specValues = new ArrayList<Double>();

        while (parser.hasNext())
        {
            if (!parser.skipToStart("protein_summary_data_filter"))
            {
                break;
            }
            minProbValues.add((double) Float.parseFloat(parser.getAttributeValue(null, "min_probability")));
            sensValues.add((double) Float.parseFloat(parser.getAttributeValue(null, "sensitivity")));
            specValues.add((double) Float.parseFloat(parser.getAttributeValue(null, "false_positive_error_rate")));
        }

        if (minProbValues.size() == 0)
            throw new XMLStreamException("No sensitivity/specificity data found in file " +
                    protXmlFile.getAbsolutePath());
        PanelWithLineChart sensSpecChart = new PanelWithLineChart(minProbValues, sensValues, "Sensitivity");
        sensSpecChart.addData(minProbValues, specValues, "Specificity");
        sensSpecChart.setName("Sens/Spec");

        return sensSpecChart;
    }

    /**
     * Generate a sensitivity-and-specificity-curve chart for two files.
     * Second one will have dashed lines
     * @param protXmlFile1
     * @param protXmlFile2
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static PanelWithChart generateSensSpecChart(File protXmlFile1, File protXmlFile2)
            throws  FileNotFoundException, XMLStreamException
    {
        FileInputStream fIn1 = new FileInputStream(protXmlFile1);
        SimpleXMLStreamReader parser1 = new SimpleXMLStreamReader(fIn1);

        List<Double> minProbValues1 = new ArrayList<Double>();
        List<Double> sensValues1 = new ArrayList<Double>();
        List<Double> specValues1 = new ArrayList<Double>();

        while (parser1.hasNext())
        {
            if (!parser1.skipToStart("protein_summary_data_filter"))
            {
                break;
            }
            minProbValues1.add((double) Float.parseFloat(parser1.getAttributeValue(null, "min_probability")));
            sensValues1.add((double) Float.parseFloat(parser1.getAttributeValue(null, "sensitivity")));
            specValues1.add((double) Float.parseFloat(parser1.getAttributeValue(null, "false_positive_error_rate")));
        }

        if (minProbValues1.size() == 0)
            throw new XMLStreamException("No sensitivity/specificity data found in file " +
                    protXmlFile1.getAbsolutePath());

        FileInputStream fIn2 = new FileInputStream(protXmlFile2);
        SimpleXMLStreamReader parser2 = new SimpleXMLStreamReader(fIn2);

        List<Double> minProbValues2 = new ArrayList<Double>();
        List<Double> sensValues2 = new ArrayList<Double>();
        List<Double> specValues2 = new ArrayList<Double>();

        while (parser2.hasNext())
        {
            if (!parser2.skipToStart("protein_summary_data_filter"))
            {
                break;
            }
            minProbValues2.add((double) Float.parseFloat(parser2.getAttributeValue(null, "min_probability")));
            sensValues2.add((double) Float.parseFloat(parser2.getAttributeValue(null, "sensitivity")));
            specValues2.add((double) Float.parseFloat(parser2.getAttributeValue(null, "false_positive_error_rate")));
        }

        if (minProbValues2.size() == 0)
            throw new XMLStreamException("No sensitivity/specificity data found in file " +
                    protXmlFile2.getAbsolutePath());



        PanelWithLineChart sensSpecChart = new PanelWithLineChart(minProbValues1, sensValues1, "Sensitivity, MS2 Only");
        sensSpecChart.addData(minProbValues1, specValues1, "Specificity, MS2 Only");

        float dash[] = { 10.0f };
        Stroke dashedStroke = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
                BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);

        sensSpecChart.addData(minProbValues2, sensValues2, "Sensitivity, Augmented");
        sensSpecChart.addData(minProbValues2, specValues2, "Specificity, Augmented");

        sensSpecChart.setSeriesColor(0, Color.BLUE);
        sensSpecChart.setSeriesColor(1, Color.RED);



        sensSpecChart.setSeriesColor(2, Color.BLUE);
        sensSpecChart.getRenderer().setSeriesStroke(2, dashedStroke);
        sensSpecChart.setSeriesColor(3, Color.RED);
        sensSpecChart.getRenderer().setSeriesStroke(3, dashedStroke);


        sensSpecChart.setName("Sens/Spec");

        return sensSpecChart;
    }



    /**
     * Create a mapping between all peptides noted in the protXml file, and all protein groups
     * that they are associated with.
     * @param protXmlFile
     * @return
     * @throws org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException
     */
    public static Map<String, Set<Integer>> loadPeptideProteinGroupMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each peptide to all proteins it supports
        Map<String, Set<Integer>> peptideProteinGroupMap = new HashMap<String, Set<Integer>>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        while (iterator.hasNext())
        {
            ProteinGroup group = iterator.next();
            if (minProteinProphet > 0 && group.getGroupProbability() < minProteinProphet)
                continue;
            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {

                for (ProtXmlReader.Peptide peptide : protXmlReaderProtein.getPeptides())
                {
                    String peptideSequence = peptide.getPeptideSequence();
                    Set<Integer> proteinGroupsThisPeptide = peptideProteinGroupMap.get(peptideSequence);
                    if (proteinGroupsThisPeptide == null)
                    {
                        proteinGroupsThisPeptide = new HashSet<Integer>();
                        peptideProteinGroupMap.put(peptideSequence, proteinGroupsThisPeptide);
                    }
                    proteinGroupsThisPeptide.add(group.getGroupNumber());
                }
            }
        }
        return peptideProteinGroupMap;
    }

    /**
     * Create a mapping between all proteins in the protxml file, and all peptides
     * associated with them.  This means all indistinguishable proteins
     * @param protXmlFile
     * @return
     * @throws org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException
     */
    public static Map<String, List<ProtXmlReader.Peptide>> loadProteinPeptideMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each protein to its peptide support
        Map<String, List<ProtXmlReader.Peptide>> proteinPeptideMap =
                new HashMap<String, List<ProtXmlReader.Peptide>>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        int numGroups=0;
        while (iterator.hasNext())
        {
            numGroups++;
            ProteinGroup group = iterator.next();
            if (minProteinProphet > 0 && group.getGroupProbability() < minProteinProphet)
                continue;
            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {
                Set<String> proteinNames = new HashSet<String>();
                //add first protein
                proteinNames.add(protXmlReaderProtein.getProteinName());
                //add indistinguishable proteins
                proteinNames.addAll(protXmlReaderProtein.getIndistinguishableProteinNames());


                for (String proteinName : proteinNames)
                {
                    List<ProtXmlReader.Peptide> peptidesThisProtein = proteinPeptideMap.get(proteinName);
                    if (peptidesThisProtein == null)
                    {
                        peptidesThisProtein = new ArrayList<ProtXmlReader.Peptide>();
                        proteinPeptideMap.put(proteinName, peptidesThisProtein);
                    }
                    peptidesThisProtein.addAll(protXmlReaderProtein.getPeptides());
                }
            }
        }
        _log.debug("loaded " + numGroups + " protein groups from file " + protXmlFile.getAbsolutePath());        
        return proteinPeptideMap;
    }

    /**
     * Create a mapping between all proteins in the protxml file, and all peptides
     * associated with them.  This means all indistinguishable proteins
     * @param protXmlFile
     * @return
     * @throws org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException
     */
    public static Map<String, Set<String>> loadProteinPeptideSequenceMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each protein to its peptide support
        Map<String, Set<String>> proteinPeptideMap = new HashMap<String, Set<String>>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        Set<String> allPeptideSet = new HashSet<String>();
        int numGroups=0;
        while (iterator.hasNext())
        {
            numGroups++;
            ProteinGroup group = iterator.next();
            if (minProteinProphet > 0 && group.getGroupProbability() < minProteinProphet)
                continue;
            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {
                Set<String> proteinNames = new HashSet<String>();
                //add first protein
                proteinNames.add(protXmlReaderProtein.getProteinName());
                //add indistinguishable proteins
                proteinNames.addAll(protXmlReaderProtein.getIndistinguishableProteinNames());

                Set<String> peptideSet = new HashSet<String>();
                for (ProtXmlReader.Peptide peptide : protXmlReaderProtein.getPeptides())
                {
                    String peptideSequence = peptide.getPeptideSequence();
                    peptideSet.add(peptideSequence);
                }

                for (String proteinName : proteinNames)
                {
                    Set<String> peptidesThisProtein = proteinPeptideMap.get(proteinName);
                    if (peptidesThisProtein == null)
                    {
                        peptidesThisProtein = new HashSet<String>();
                        proteinPeptideMap.put(proteinName, peptidesThisProtein);
                    }
                    peptidesThisProtein.addAll(peptideSet);
                    allPeptideSet.addAll(peptidesThisProtein);
                }
            }
        }
        _log.debug("loaded " + numGroups + " protein groups, with " + allPeptideSet.size() +
                " unique peptides, from file " + protXmlFile.getAbsolutePath());        
        return proteinPeptideMap;
    }

    /**
     * Create a mapping between all proteins in the protxml file, and all peptides
     * associated with them.  This means all indistinguishable proteins
     * @param protXmlFile
     * @return
     * @throws org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException
     */
    public static Map<String, Float> loadProteinProbabilityMapFromProtXML(File protXmlFile)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        //build a map from each protein to its highest score
        Map<String, Float> proteinScoreMap = new HashMap<String, Float>();

        Iterator<ProteinGroup> iterator = protXmlReader.iterator();
        while (iterator.hasNext())
        {
            ProteinGroup group = iterator.next();
            float groupScore = group.getProbability();

            for (ProtXmlReader.Protein protXmlReaderProtein : group.getProteins())
            {
                Set<String> proteinNames = new HashSet<String>();
                //add first protein
                proteinNames.add(protXmlReaderProtein.getProteinName());
                //add indistinguishable proteins
                proteinNames.addAll(protXmlReaderProtein.getIndistinguishableProteinNames());


                for (String proteinName : proteinNames)
                {
                    Float scoreThisProtein = proteinScoreMap.get(proteinName);
                    if (scoreThisProtein == null)
                    {
                        scoreThisProtein = protXmlReaderProtein.getProbability();
                        proteinScoreMap.put(proteinName, scoreThisProtein);
                    }
                    if (groupScore > scoreThisProtein)
                        proteinScoreMap.put(proteinName, groupScore);
                }
            }
        }
        return proteinScoreMap;
    }

    /**
     * WARNING!  after running through the iterator, the proteins in each group disappear
     * @param protXmlFile
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static List<ProteinGroup> loadProteinGroupsFromProtXML(File protXmlFile)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);
        ProtXmlReader.ProteinGroupIterator pgIterator = protXmlReader.iterator();
        List<ProteinGroup> proteinGroupList = new ArrayList<ProteinGroup>();
        while (pgIterator.hasNext())
        {
            ProteinGroup pg = pgIterator.next();
            //if I don't reference the proteins list, it disappears
            //todo: maybe do this in some better way?
            pg.getProteins().size();
            proteinGroupList.add(pg);
        }
        return proteinGroupList;
    }

    public static List<ProtXmlReader.Protein> loadProtXmlProteinsFromProtXML(File protXmlFile)
            throws FileNotFoundException, XMLStreamException
    {
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);
        ProtXmlReader.ProteinGroupIterator pgIterator = protXmlReader.iterator();
        List<ProtXmlReader.Protein> result = new ArrayList<ProtXmlReader.Protein>();
        int numProteinGroups = 0;
        while (pgIterator.hasNext())
        {
            ProteinGroup proteinGroup = pgIterator.next();
            for (ProtXmlReader.Protein protein : proteinGroup.getProteins())
            {
                result.add(protein);
            }
            numProteinGroups++;
        }
        _log.debug("Loaded proteins from " + numProteinGroups + " protein groups");

        return result;
    }

    /**
     * returns null if protein not found with a minimum group probability
     * @param protXmlFile
     * @param proteinName
     * @param minProteinProphetGroupProbability
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static ProtXmlReader.Protein loadFirstProteinOccurrence(File protXmlFile, String proteinName,
                                                                   float minProteinProphetGroupProbability)
            throws FileNotFoundException, XMLStreamException
    {
        List<String> proteins = new ArrayList<String>();
        proteins.add(proteinName);
        return loadFirstProteinOccurrence(protXmlFile, proteins, minProteinProphetGroupProbability).get(proteinName);
    }

    /**
     * returns null if protein not found
     * @param protXmlFile
     * @param proteinName
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static ProtXmlReader.Protein loadFirstProteinOccurrence(File protXmlFile, String proteinName)
            throws FileNotFoundException, XMLStreamException
    {
        return loadFirstProteinOccurrence(protXmlFile, proteinName, 0f);
    }

    /**
     * Returns a map from protein names to first occurrences of proteins in the protXML file,
     * if they exist with minimum probability
     * @param protXmlFile
     * @param proteinNames
     * @param minProteinProphetGroupProbability
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static Map<String, ProtXmlReader.Protein> loadFirstProteinOccurrence(File protXmlFile,
                                                                                Collection<String> proteinNames,
                                                                   float minProteinProphetGroupProbability)
            throws FileNotFoundException, XMLStreamException
    {
        Map<String, ProtXmlReader.Protein> result = new HashMap<String, ProtXmlReader.Protein>();

        ProtXmlReader.ProteinGroupIterator proteinIterator = new ProtXmlReader(protXmlFile).iterator();
        while (proteinIterator.hasNext())
        {
            Set<String> remainingProteinNames = new HashSet<String>(proteinNames);
            ProteinGroup proteinGroup = proteinIterator.next();
            if (proteinGroup.getProbability() < minProteinProphetGroupProbability)
                continue;
            for (ProtXmlReader.Protein protein : proteinGroup.getProteins())
            {
                List<String> allProteinNamesThisProtein = new ArrayList<String>();
                allProteinNamesThisProtein.add(protein.getProteinName());
                if (protein.getIndistinguishableProteinNames() != null)
                {
                    allProteinNamesThisProtein.addAll(protein.getIndistinguishableProteinNames());
                }

                List<String> matchedNamesThisProtein = new ArrayList<String>();
                for (String proteinName : allProteinNamesThisProtein)
                {
                    if (proteinNames.contains(proteinName))
                        matchedNamesThisProtein.add(proteinName);
                }

                if (!matchedNamesThisProtein.isEmpty())
                {
                    for (String proteinName : matchedNamesThisProtein)
                    {
                        result.put(proteinName, protein);
                        remainingProteinNames.remove(proteinName);
                    }
                }
            }
        }
        return result;
    }


    /**
     * Create a pepxml file that can be used with proteinprophet
     */
    public static void createPepXml(Feature[] featuresWithPeptides, File fastaFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        for (Feature feature : featuresWithPeptides)
        {
            //assign all features a peptide ID probability of 1
            MS2ExtraInfoDef.setPeptideProphet(feature, 1f);
            feature.setProperty("protein",null);

            //proteinprophet can't handle charges >3
            if (feature.getCharge() > 3)
                feature.setCharge(3);
        }


        assignContainingProteinsToFeatures(featuresWithPeptides, fastaFile);

        FeaturePepXmlWriter pepXmlWriter =
                new FeaturePepXmlWriter(featuresWithPeptides,
                                        null);
        pepXmlWriter.setSearchDatabase(fastaFile.getAbsolutePath());
        try
        {
            pepXmlWriter.write(outputFile);
            ApplicationContext.infoMessage("Wrote output file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            _log.error("Failed to save pepXML",e);
        }
    }

    public static Map<String, Protein> loadProteinNameProteinMapFromFasta(File fastaFile)
    {
        Map<String, Protein> result = new HashMap<String, Protein>();
        FastaLoader fastaLoader = new FastaLoader(fastaFile);
        FastaLoader.ProteinIterator iterator = fastaLoader.iterator();

        while (iterator.hasNext())
        {
            Protein protein = iterator.next();
            result.put(protein.getLookup(), protein);
        }
        return result;
    }

    /**
     * Load all proteins from a fasta file
     * @param fastaFile
     * @return
     */
    public static ArrayList<Protein> loadProteinsFromFasta(File fastaFile)
    {
        ArrayList<Protein> proteinArray = new ArrayList<Protein>();
        FastaLoader fastaLoader = new FastaLoader(fastaFile);
        FastaLoader.ProteinIterator iterator = fastaLoader.iterator();

        while (iterator.hasNext())
        {
            Protein protein = iterator.next();
            proteinArray.add(protein);
        }
        return proteinArray;
    }

    /**
     * Load a map from protein ID to protein sequence
     * @param fastaFile
     * @return
     */
    public static Map<String, String> loadProteinSequenceMapFromFasta(File fastaFile)
    {
        List<Protein> proteins = loadProteinsFromFasta(fastaFile);
        Map<String, String> result = new HashMap<String, String>();
        for (Protein protein : proteins)
        {
            result.put(protein.getLookup(), protein.getSequenceAsString());
        }
        return result;
    }

    public static Set<String> loadTrypticPeptidesFromFasta(File fastaFile)
    {
        FastaLoader fastaLoader = new FastaLoader(fastaFile);
        FastaLoader.ProteinIterator iterator = fastaLoader.iterator();

        Set<String> peptides = new HashSet<String>();
        PeptideGenerator pg = new PeptideGenerator();

        while (iterator.hasNext())
        {
            Protein protein = iterator.next();
            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
                peptides.add(new String(peptideThisProtein.getChars()));
        }
        return peptides;
    }

    public static Map<String, List<String>> loadTrypticPeptideProteinMapFromFasta(File fastaFile)
    {
        return loadTrypticPeptideProteinMapFromFasta(fastaFile, 0);
    }

    public static Map<String, List<String>> loadTrypticPeptideProteinMapFromFasta(File fastaFile, int maxMissedCleavages)
    {
        FastaLoader fastaLoader = new FastaLoader(fastaFile);
        FastaLoader.ProteinIterator iterator = fastaLoader.iterator();

        PeptideGenerator pg = new PeptideGenerator();
        pg.setMaxMissedCleavages(maxMissedCleavages);
        Map<String, List<String>> result = new HashMap<String, List<String>>();

        while (iterator.hasNext())
        {
            Protein protein = iterator.next();
            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
            {
                String peptide = new String(peptideThisProtein.getChars());
                List<String> proteinsThisPeptide = result.get(peptide);
                if (proteinsThisPeptide == null)
                {
                    proteinsThisPeptide = new ArrayList<String>();
                    result.put(peptide, proteinsThisPeptide);
                }
                proteinsThisPeptide.add(protein.getLookup());
            }
        }
        return result;
    }

    public static List<String> getTrypticPeptidesForProtein(Protein protein, int maxMissedCleavages)
    {
        PeptideGenerator pg = new PeptideGenerator();
        pg.setMaxMissedCleavages(maxMissedCleavages);

        List<String> result = new ArrayList<String>();
        Peptide[] peptidesThisProtein = pg.digestProtein(protein);
        for (Peptide peptideThisProtein : peptidesThisProtein)
        {
            result.add(new String(peptideThisProtein.getChars()));
        }
        return result;
    }


    
    public static Map<String, Set<String>> findFastaProteinsForPeptides(Collection<String> peptideList, File fastaFile)
    {
        FastaLoader fastaLoader = new FastaLoader(fastaFile);
        FastaLoader.ProteinIterator iterator = fastaLoader.iterator();

        Map<String, Set<String>> result = new HashMap<String, Set<String>>();

        PeptideGenerator pg = new PeptideGenerator();
        while (iterator.hasNext())
        {
            Protein protein = iterator.next();

            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
            {
                String peptideSequence = new String(peptideThisProtein.getChars());
                if (peptideList.contains(peptideSequence))
                {
                    Set<String> proteinsThisPeptide = result.get(peptideSequence);
                    if (proteinsThisPeptide == null)
                    {
                        proteinsThisPeptide = new HashSet<String>();
                        result.put(peptideSequence, proteinsThisPeptide);
                    }
                    proteinsThisPeptide.add(protein.getName());
                }
            }
        }
        return result;
    }

    /**
     *
     * @param ms1FeaturesWithPeptides
     * @param fastaFile
     */
    public static void assignContainingProteinsToFeatures(Feature[] ms1FeaturesWithPeptides, File fastaFile)
    {
        Protein[] proteinsInFasta = loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);

        for (int i=0; i<ms1FeaturesWithPeptides.length; i++)
        {
            Feature feature = ms1FeaturesWithPeptides[i];


            if (i % (ms1FeaturesWithPeptides.length/100) == 0)
                ApplicationContext.setMessage(Rounder.round((((float) i) * 100f/ms1FeaturesWithPeptides.length),0) + "% complete");

            String peptideSequence = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptideSequence == null)
                continue;
            for (Protein protein : proteinsInFasta)
            {
                String proteinSequence = protein.getSequenceAsString();
                if (proteinSequence.contains(peptideSequence))
                {
                    MS2ExtraInfoDef.addProtein(feature, protein.getLookup());
                }
            }
        }
    }

    /**
     * Map peptides to proteins using multiple protxml files
     * @param peptides
     * @param protXmlFiles
     * @param proteinsInFasta
     * @return
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public static Map<String,List<Protein>> mapPeptidesToProteins(Set<String> peptides,
                                                              File[] protXmlFiles,
                                                              Protein[] proteinsInFasta,
                                                              double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        Map<String,List<Protein>> result = new HashMap<String,List<Protein>>();
        for (File protXmlFile : protXmlFiles)
        {
            Map<String,List<Protein>> mapThisFile =
                    mapPeptidesToProteins(peptides, protXmlFile, proteinsInFasta, minProteinProphet);
            for (String peptide : mapThisFile.keySet())
            {
                if (result.containsKey(peptide))
                {
                    List<Protein> existingProteins = mapThisFile.get(peptide);
                    for (Protein newProtein : mapThisFile.get(peptide))
                    {
                        boolean foundIt = false;
                        for (Protein existingProtein : existingProteins)
                        {
                             if (existingProtein.getLookup().equals(newProtein.getLookup()))
                             {
                                 foundIt = true;
                                 break;
                             }
                        }
                        if (!foundIt)
                            existingProteins.add(newProtein);
                    }
                }
                else
                {
                    result.put(peptide, mapThisFile.get(peptide));
                }
            }
        }
        return result;
    }

    /**
     * Map peptides to proteins.  If the user has supplied a fasta file, do this based on
     * sequence.  If the user has supplied a protxml file, do this based on that.
     *
     * protxml file takes precedence
     * @param peptides
     * @return
     */
    public static Map<String,List<Protein>> mapPeptidesToProteins(Set<String> peptides,
                                                              File protXmlFile,
                                                              Protein[] proteinsInFasta,
                                                              double minProteinProphet)
            throws FileNotFoundException, XMLStreamException
    {
        Map<String,List<Protein>> result = new HashMap<String,List<Protein>>();

        int i=0;

        if (protXmlFile != null)
        {
                Map<String, Set<String>> peptideIPIMap =
                        loadPeptideProteinMapFromProtXML(protXmlFile,minProteinProphet);
                for (String peptide : peptideIPIMap.keySet())
                {
                    List<Protein> protList = new ArrayList<Protein>();
                    for(String ipi : peptideIPIMap.get(peptide))
                    {
                        Protein protein = new Protein(ipi, "".getBytes());
                        protList.add(protein);
                    }
                    result.put(peptide, protList);
                }
        }
        else
        {
            for (String peptide : peptides)
            {

                if (i++ % 100 == 0)
                    ApplicationContext.setMessage("Protein mapping peptide " + i + " / " + peptides.size());
                for (Protein protein : proteinsInFasta)
                {
                    String proteinSequence = protein.getSequenceAsString();

                    if (proteinSequence.contains(peptide))
                    {
                        List<Protein> proteinList = result.get(peptide);
                        if (proteinList == null)
                        {
                            proteinList = new ArrayList<Protein>();
                            result.put(peptide, proteinList);
                        }
                        proteinList.add(protein);
                    }
                }
            }
        }

        return result;
    }


    /**
     * helper method for one featureset
     * @param featureSet
     * @param fastaFile
     */
    public static void guessProteinsForFeaturePeptides(FeatureSet featureSet, File fastaFile)
    {
         guessProteinsForFeaturePeptides(new FeatureSet[] { featureSet }, fastaFile);
    }

    /**
     * helper method for one featureset
     * @param featureSet
     * @param fastaProteins
     */
    public static void guessProteinsForFeaturePeptides(FeatureSet featureSet, Protein[] fastaProteins)
    {
         guessProteinsForFeaturePeptides(new FeatureSet[] { featureSet }, null, fastaProteins);
    }



    /**
     * cover method
     * @param featureSets
     * @param fastaFile
     */
    public static void guessProteinsForFeaturePeptides(FeatureSet[] featureSets, File fastaFile)
    {
        List<Protein> fastaProteins = ProteinUtilities.loadProteinsFromFasta(fastaFile);
        guessProteinsForFeaturePeptides(featureSets, fastaFile, fastaProteins.toArray(new Protein[fastaProteins.size()]));
    }


    /**
     * For every feature with a peptide in every featureset passed in, find some protein in the fasta
     * file that contains that peptide, and assign it
     * @param featureSets
     * @param fastaProteins
     */
    public static void guessProteinsForFeaturePeptides(FeatureSet[] featureSets, File fastaFile, Protein[] fastaProteins)
    {

        Set<String> peptidesRemainingInAllFiles = new HashSet<String>();

        for (FeatureSet featureSet : featureSets)
        {
            String path = "unknown";
            if (featureSet.getSourceFile() != null)
                path = featureSet.getSourceFile().getAbsolutePath();
            ApplicationContext.setMessage("Getting peptides from file " + path);
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null)
                    peptidesRemainingInAllFiles.add(peptide);
            }
        }
        ApplicationContext.setMessage("Total peptides: " + peptidesRemainingInAllFiles.size());
        ApplicationContext.setMessage("Looking for tryptic peptides...");

        Map<String, Protein> peptideProteinMap =   new HashMap<String, Protein>();

        PeptideGenerator pg = new PeptideGenerator();
        pg.setMaxMissedCleavages(2);

        for (Protein protein : fastaProteins)
        {
            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
            {
                String peptideSequence = new String(peptideThisProtein.getChars());
                if (peptidesRemainingInAllFiles.contains(peptideSequence))
                {
                    peptideProteinMap.put(peptideSequence, protein);
                    peptidesRemainingInAllFiles.remove(peptideSequence);
                }
            }
        }

        ApplicationContext.setMessage("After tryptic digest, " +
                peptidesRemainingInAllFiles.size() + " peptides remain.");

        int numProteins = fastaProteins.length;


        int currentProteinIndex = 0;


        //the PeptideGenerator method won't always work -- semitryptic searches, etc.
        if (peptidesRemainingInAllFiles.size() > 0)
        {
            ApplicationContext.setMessage("Doing nontryptic on " + peptidesRemainingInAllFiles.size() +
                    " remaining peptides...");

            for (Protein protein : fastaProteins)
            {
                if (peptidesRemainingInAllFiles.size() == 0)
                    break;
                if (currentProteinIndex % (Math.max(numProteins / 100,1)) == 0)
                    ApplicationContext.setMessage("\t" + (currentProteinIndex * 100f / numProteins) + "% complete... remaining peptides: " + peptidesRemainingInAllFiles.size());
                currentProteinIndex++;
                Set<String> peptidesFound = new HashSet<String>();

                String proteinSequence = protein.getSequenceAsString();
                for (String peptide : peptidesRemainingInAllFiles)
                {
                    if (proteinSequence.contains(peptide))
                    {
                        peptidesFound.add(peptide);
                        peptideProteinMap.put(peptide, protein);
                    }
                }
                peptidesRemainingInAllFiles.removeAll(peptidesFound);
            }
        }
        ApplicationContext.setMessage("All peptides assigned proteins");

        for (FeatureSet featureSet : featureSets)
        {
            if (fastaFile != null)
                MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSet, fastaFile.getAbsolutePath());
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null)
                {
                    Protein protein = peptideProteinMap.get(peptide);
                    if (protein == null)
                    {
                        ApplicationContext.infoMessage("\tWARNING!!! no protein for peptide " + peptide);
                    }
                    else
                    {
                        List<String> proteinList = new ArrayList<String>();
                        proteinList.add(protein.getLookup());
                        MS2ExtraInfoDef.setProteinList(feature, proteinList);
                        Pair<Character, Character> prevNextAAs = getPrevNextAAs(peptide, protein);
                        char prevAA = prevNextAAs.first;
                        char nextAA = prevNextAAs.second;
                        MS2ExtraInfoDef.setPrevAminoAcid(feature, prevAA);
                        MS2ExtraInfoDef.setNextAminoAcid(feature, nextAA);

                        //WARNING WARNING WARNING!!!
                        //This behavior is trypsin-specific.  If another enzyme is used, number of enzymatic
                        //ends will be set incorrectly.
                        //check for start of protein sequence or trypsin digestion at start of peptide  (remember proline)
                        int numTrypticEnds = 0;
                        if (prevAA == '-' ||
                                (prevAA == 'K' || prevAA == 'R') && !peptide.startsWith("P"))
                            numTrypticEnds++;
                        //check for end of protein sequence or trypsin digestion at end of peptide
                        if (nextAA == ('-') ||
                                (nextAA != 'P' && (peptide.endsWith("K") ||
                                        peptide.endsWith("R"))))
                            numTrypticEnds++;
//if (numTrypticEnds != 2) System.err.println("*** " + numTrypticEnds + ", " + peptide + ", " + prevAA + ", " + nextAA);
                        MS2ExtraInfoDef.setNumEnzymaticEnds(feature, numTrypticEnds);
                    }
                }
            }
        }
    }



    /**
     * For every feature with a peptide in every featureset passed in, find ALL proteins in the fasta
     * file that contains that peptide, and assign them all.  This is much more computationally intensive than
     * finding just one
     * @param featureSets
     * @param fastaProteins
     */
    public static void guessAllProteinsForFeaturePeptides(FeatureSet[] featureSets, File fastaFile, Protein[] fastaProteins)
    {

        Set<String> peptidesWithNoProteins = new HashSet<String>();

        for (FeatureSet featureSet : featureSets)
        {
            ApplicationContext.setMessage("Getting peptides from file " +
                    featureSet.getSourceFile().getAbsolutePath());
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null)
                    peptidesWithNoProteins.add(peptide);
            }
        }
        ApplicationContext.setMessage("Total peptides: " + peptidesWithNoProteins.size());
        ApplicationContext.setMessage("Looking for tryptic peptides...");

        Map<String, List<Protein>> peptideProteinMap =   new HashMap<String, List<Protein>>();

        PeptideGenerator pg = new PeptideGenerator();
        pg.setMaxMissedCleavages(2);

        for (Protein protein : fastaProteins)
        {
            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
            {
                String peptideSequence = new String(peptideThisProtein.getChars());
                List<Protein> proteinsThisPeptide = peptideProteinMap.get(peptideSequence);
                if (proteinsThisPeptide == null)
                {
                    proteinsThisPeptide = new ArrayList<Protein>();
                    peptideProteinMap.put(peptideSequence, proteinsThisPeptide);
                }
                proteinsThisPeptide.add(protein);
                if (peptidesWithNoProteins.contains(peptideSequence))
                {
                    peptidesWithNoProteins.remove(peptideSequence);
                }
            }
        }

        ApplicationContext.setMessage("After tryptic digest, " +
                peptidesWithNoProteins.size() + " peptides remain.");

        int numProteins = fastaProteins.length;


        int currentProteinIndex = 0;


        //the PeptideGenerator method won't always work -- semitryptic searches, etc.
        if (peptidesWithNoProteins.size() > 0)
        {
            ApplicationContext.setMessage("Doing nontryptic on " + peptidesWithNoProteins.size() +
                    " remaining peptides...");

            for (Protein protein : fastaProteins)
            {
                if (peptidesWithNoProteins.size() == 0)
                    break;
                if (currentProteinIndex % (numProteins / 100) == 0)
                    ApplicationContext.setMessage("\t" + (currentProteinIndex * 100f / numProteins) + "% complete... remaining peptides: " + peptidesWithNoProteins.size());
                currentProteinIndex++;
                Set<String> peptidesFound = new HashSet<String>();

                String proteinSequence = protein.getSequenceAsString();
                for (String peptide : peptidesWithNoProteins)
                {
                    if (proteinSequence.contains(peptide))
                    {
                        peptidesFound.add(peptide);
                        List<Protein> proteinsThisPeptide = peptideProteinMap.get(peptide);
                        if (proteinsThisPeptide == null)
                        {
                            proteinsThisPeptide = new ArrayList<Protein>();
                            peptideProteinMap.put(peptide, proteinsThisPeptide);
                        }
                        proteinsThisPeptide.add(protein);
                    }
                }
            }
        }
        ApplicationContext.setMessage("All peptides assigned proteins");

        for (FeatureSet featureSet : featureSets)
        {
            if (fastaFile != null)
                MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSet, fastaFile.getAbsolutePath());
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide != null)
                {
                    List<Protein> proteins = peptideProteinMap.get(peptide);
                    if (proteins == null)
                    {
                        ApplicationContext.infoMessage("\tWARNING!!! no protein for peptide " + peptide);
                    }
                    else
                    {
                        List<String> proteinList = new ArrayList<String>();
                        for (Protein protein : proteins)
                            proteinList.add(protein.getLookup());
                        MS2ExtraInfoDef.setProteinList(feature, proteinList);
                        Pair<Character, Character> prevNextAAs = getPrevNextAAs(peptide, proteins.get(0));
                        char prevAA = prevNextAAs.first;
                        char nextAA = prevNextAAs.second;
                        MS2ExtraInfoDef.setPrevAminoAcid(feature, prevAA);
                        MS2ExtraInfoDef.setNextAminoAcid(feature, nextAA);

                        //WARNING WARNING WARNING!!!
                        //This behavior is trypsin-specific.  If another enzyme is used, number of enzymatic
                        //ends will be set incorrectly.
                        //check for start of protein sequence or trypsin digestion at start of peptide  (remember proline)
                        int numTrypticEnds = 0;
                        if (prevAA == '-' ||
                                (prevAA == 'K' || prevAA == 'R') && !peptide.startsWith("P"))
                            numTrypticEnds++;
                        //check for end of protein sequence or trypsin digestion at end of peptide
                        if (nextAA == ('-') ||
                                (nextAA != 'P' && (peptide.endsWith("K") ||
                                        peptide.endsWith("R"))))
                            numTrypticEnds++;
//if (numTrypticEnds != 2) System.err.println("*** " + numTrypticEnds + ", " + peptide + ", " + prevAA + ", " + nextAA);
                        MS2ExtraInfoDef.setNumEnzymaticEnds(feature, numTrypticEnds);
                    }
                }
            }
        }
    }



    protected static Pair<Character, Character> getPrevNextAAs(String peptide, Protein protein)
    {
        Character prevAA = '-';
        Character nextAA = '-';

        String proteinSequence = protein.getSequenceAsString();
        int prevAAIndex = proteinSequence.indexOf(peptide) - 1;
        if (prevAAIndex > 0)
            prevAA = proteinSequence.charAt(prevAAIndex);
        int nextAAIndex = proteinSequence.indexOf(peptide) + peptide.length();
        if (nextAAIndex >= peptide.length() && nextAAIndex < proteinSequence.length())
            nextAA = proteinSequence.charAt(nextAAIndex);

        return new Pair<Character, Character>(prevAA, nextAA);
    }
}
