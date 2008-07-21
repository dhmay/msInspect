/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.ms2;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeaturePepXmlWriter;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithLineChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.SimpleXMLStreamReader;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
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

    /**
     * Create a mapping between all peptides noted in the protXml file, and all proteins
     * that they are associated with.  This means all indistinguishable proteins
     * @param protXmlFile
     * @return
     * @throws CommandLineModuleExecutionException
     */
    public static Map<String, Set<String>> loadPeptideProteinMapFromProtXML(File protXmlFile,
                                                                            double minProteinProphet)
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
     * @throws CommandLineModuleExecutionException
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
     * @throws CommandLineModuleExecutionException
     */
    public static Map<String, Set<String>> loadProteinPeptideMapFromProtXML(File protXmlFile,
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
     * @throws CommandLineModuleExecutionException
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
            ApplicationContext.setMessage("Getting peptides from file " +
                    featureSet.getSourceFile().getAbsolutePath());
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

        ApplicationContext.setMessage("Doing nontryptic...");

        //the PeptideGenerator method won't always work -- semitryptic searches, etc.
        if (peptidesRemainingInAllFiles.size() > 0)
        {
            for (Protein protein : fastaProteins)
            {
                if (peptidesRemainingInAllFiles.size() == 0)
                    break;
                if (currentProteinIndex % (numProteins / 100) == 0)
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
                    }
                }
            }
        }

    }
}