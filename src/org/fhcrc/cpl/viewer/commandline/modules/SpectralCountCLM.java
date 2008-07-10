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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.viewer.ms2.GeneMappingUtilities;
import org.fhcrc.cpl.toolbox.TabWriter;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;


/**
 *
 */
public class SpectralCountCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SpectralCountCLM.class);

    protected File protXmlFile = null;
    protected FeatureSet[] ms2FeatureSets = null;
    protected File geneLookupFile = null;
    protected File outFile = null;

    protected boolean showCharts = false;

    protected int mode;

    protected final static String[] modeStrings =
            {
                    "peptide",
                    "gene",
                    "proteingroup",
            };

    protected final static String[] modeExplanations =
            {
                    "Peptide-level counts",
                    "Gene-level counts",
                    "protein group-level counts",
            };

    protected static final int MODE_PEPTIDE = 0;
    protected static final int MODE_GENE = 1;
    protected static final int MODE_PROTEIN_GROUP = 2;



    public SpectralCountCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "spectralcount";
        mShortDescription = "Create a spreadsheet with spectral count information";
        mHelpMessage = "Create a spreadsheet with spectral count information";
        CommandLineArgumentDefinition[] argDefs =
                {
                    createEnumeratedArgumentDefinition("mode",true,modeStrings, modeExplanations),
                    createFileToReadArgumentDefinition("protxml", false, "ProtXML file"),
                    createUnnamedSeriesFileArgumentDefinition(true, "MS2 Feature file(s)"),
                    createFileToReadArgumentDefinition("genelookupfile", false,
                            "Gene lookup file for IPI numbers"),
                    createFileToWriteArgumentDefinition("out", false, "output file"),
                    createBooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFile = getFileArgumentValue("protxml");
        File[] ms2FeatureFiles = this.getUnnamedSeriesFileArgumentValues();

        ms2FeatureSets = new FeatureSet[ms2FeatureFiles.length];
        try
        {
            for (int i=0; i<ms2FeatureFiles.length; i++)
                ms2FeatureSets[i] = new FeatureSet(ms2FeatureFiles[i]);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }


        geneLookupFile = getFileArgumentValue("genelookupfile");

        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        switch(mode)
        {
            case MODE_PEPTIDE:
                assertArgumentAbsent("protxml","mode");
                assertArgumentAbsent("genelookupfile","mode");
                break;
            case MODE_GENE:
                assertArgumentPresent("protxml","mode");
                assertArgumentPresent("genelookupfile","mode");
                break;
            case MODE_PROTEIN_GROUP:
                assertArgumentPresent("protxml","mode");
                break;
        }

        outFile = getFileArgumentValue("out");
        showCharts = getBooleanArgumentValue("showcharts");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Map<String,Integer>[] peptideSpectralCountMaps = new Map[ms2FeatureSets.length];


        for (int i=0; i<ms2FeatureSets.length; i++)
        {
            Map<String,Integer> peptideSpectralCountMap = new HashMap<String,Integer>();
            for (Feature feature : ms2FeatureSets[i].getFeatures())
            {
                String featurePeptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (featurePeptide == null)
                    continue;
                Integer currentCount = peptideSpectralCountMap.get(featurePeptide);
                if (currentCount == null)
                {
                    currentCount = 0;
                }
                peptideSpectralCountMap.put(featurePeptide, currentCount+1);
            }
            peptideSpectralCountMaps[i] = peptideSpectralCountMap;
        }

        switch(mode)
        {
            case MODE_PEPTIDE:
                doPeptides(peptideSpectralCountMaps);
                break;
            case MODE_GENE:
                doGenes(peptideSpectralCountMaps);
                break;
            case MODE_PROTEIN_GROUP:
                doProteinGroups(peptideSpectralCountMaps);
                break;
        }

    }



    protected void doPeptides(Map<String,Integer>[] peptideSpectralCountMaps)
            throws CommandLineModuleExecutionException
    {


        if (outFile != null)
        {
            Map<String,Integer> peptideSpectralCountMap = peptideSpectralCountMaps[0];
            String[] columns = new String[] {"peptide","spectra"};

            List<Map<String,Object>> rowsList =
                    new ArrayList<Map<String,Object>>(peptideSpectralCountMap.size());
            for (String peptide : peptideSpectralCountMap.keySet())
            {
                Map<String,Object> row = new HashMap<String,Object>();
                row.put("peptide",peptide);
                row.put("spectra",peptideSpectralCountMap.get(peptide));
                rowsList.add(row);
            }
            TabWriter tabWriter = new TabWriter(columns, rowsList, outFile);
            try
            {
                tabWriter.write();
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        if (showCharts)
        {

            PanelWithHistogram pwh = new PanelWithHistogram("peptide spectral counts");
            pwh.setOffsetSeries(true);
            for (int i=0; i<peptideSpectralCountMaps.length; i++)
            {
                List<Double> specCounts = new ArrayList<Double>();
                for (Integer specCount : peptideSpectralCountMaps[i].values())
                    specCounts.add((double) specCount);
                pwh.addData(specCounts, ms2FeatureSets[i].getSourceFile().getName());
            }

            ChartDialog cd = new ChartDialog(pwh);
            cd.setVisible(true);
        }

    }

    protected void doProteinGroups(Map<String,Integer>[] peptideSpectralCountMaps)
            throws CommandLineModuleExecutionException
    {
        Map<String,Integer> peptideSpectralCountMap = peptideSpectralCountMaps[0];

        Map<String, Set<Integer>> peptideGroupMap;
        try
        {
            peptideGroupMap = ProteinUtilities.loadPeptideProteinGroupMapFromProtXML(protXmlFile,0);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Error parsing protxml file",e);
        }


        Map<Integer, Integer> proteinGroupSpectralCountMap = new HashMap<Integer,Integer>();
        for (String peptide : peptideGroupMap.keySet())
        {
            if (!peptideSpectralCountMap.containsKey(peptide))
                continue;
            Set<Integer> proteinGroupsThisPeptide = peptideGroupMap.get(peptide);
            for (Integer proteinGroup : proteinGroupsThisPeptide)
            {
                Integer proteinGroupSpectralCount = proteinGroupSpectralCountMap.get(proteinGroup);
                if (proteinGroupSpectralCount == null)
                {
                    proteinGroupSpectralCount = 0;
                }
                proteinGroupSpectralCountMap.put(proteinGroup,
                        proteinGroupSpectralCount + peptideSpectralCountMap.get(peptide));                
            }
        }

        String[] columns = new String[] {"proteingroup","spectra"};

        List<Map<String,Object>> rowsList =
                new ArrayList<Map<String,Object>>(proteinGroupSpectralCountMap.size());
        for (Integer proteinGroup : proteinGroupSpectralCountMap.keySet())
        {         
            Map<String,Object> row = new HashMap<String,Object>();
            row.put("proteingroup",proteinGroup);
            row.put("spectra",proteinGroupSpectralCountMap.get(proteinGroup));
            rowsList.add(row);
        }

        if (outFile != null)
        {
            TabWriter tabWriter = new TabWriter(columns, rowsList, outFile);
            try
            {
                tabWriter.write();
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
    }

    protected void doGenes(Map<String,Integer>[] peptideSpectralCountMaps)
            throws CommandLineModuleExecutionException
    {
        Map<String,Integer> peptideSpectralCountMap = peptideSpectralCountMaps[0];

        Map<String, Set<String>> peptideIPIMap;
        try
        {
            peptideIPIMap = ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFile,0);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Error parsing protxml file",e);
        }

        Map<String,String[]> ipiGeneArrayMap =
                GeneMappingUtilities.loadIPIGeneMap(geneLookupFile);
        Set<String> geneSet = new HashSet<String>();
        for (String[] genearray : ipiGeneArrayMap.values())
            for (String gene : genearray)
                geneSet.add(gene);

        Map<String, Integer> geneSpectralCountMap = new HashMap<String,Integer>();
        for (String peptide : peptideSpectralCountMap.keySet())
        {
            Set<String> genesAlreadyUpdated = new HashSet<String>();
            Set<String> proteinsThisPeptide = peptideIPIMap.get(peptide);
            if (proteinsThisPeptide == null)
                continue;
            for (String ipi : proteinsThisPeptide)
            {
                String[] genesThisIPI = ipiGeneArrayMap.get(ipi);
                if (genesThisIPI == null)
                    continue;
                for (String gene : genesThisIPI)
                {
                    if (genesAlreadyUpdated.contains(gene))
                        continue;
                    Integer geneSpectralCount = geneSpectralCountMap.get(gene);
                    if (geneSpectralCount == null)
                        geneSpectralCount = 0;
                    geneSpectralCountMap.put(gene, geneSpectralCount + peptideSpectralCountMap.get(peptide));
                    genesAlreadyUpdated.add(gene);
                }

            }
        }

        String[] columns = new String[] {"gene","spectra"};

        List<Map<String,Object>> rowsList =
                new ArrayList<Map<String,Object>>(geneSpectralCountMap.size());
        for (String peptide : peptideSpectralCountMap.keySet())
        {
            Map<String,Object> row = new HashMap<String,Object>();
            row.put("gene",peptide);
            row.put("spectra",peptideSpectralCountMap.get(peptide));
            rowsList.add(row);
        }

        if (outFile != null)
        {
            TabWriter tabWriter = new TabWriter(columns, rowsList, outFile);
            try
            {
                tabWriter.write();
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
    }

}
