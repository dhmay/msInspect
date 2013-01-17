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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.quant.gui.ProteinQuantSummaryFrame;
import org.fhcrc.cpl.viewer.quant.gui.ProteinSummarySelectorFrame;
import org.fhcrc.cpl.viewer.quant.gui.QuantitationReviewer;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


/**
 * test
 */
public class ProteinQuantChartsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProteinQuantChartsCLM.class);
    public File protXmlFile;
    public File pepXmlFile;
    public File outDir;
    public File mzXmlDir;
    public Boolean appendOutput = true;
    public float minProteinProphet = 0.9f;
    public float minHighRatio = 0f;
    public float maxLowRatio = 999f;

    public Map<String, List<String>> proteinGeneListMap;
    public List<ProtXmlReader.Protein> proteins;
    public File outFile;

    protected ProteinQuantSummaryFrame quantSummaryFrame;



    protected boolean hasRun = false;

    public ProteinQuantChartsCLM()
    {
        init(true);
    }

    public ProteinQuantChartsCLM(boolean shouldCaptureOutputArgs)
    {
        init(shouldCaptureOutputArgs);
    }

    protected void init(boolean shouldCaptureOutputArgs)
    {
        mCommandName = "proteinquantcharts";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("protxml", true, "ProtXML file with protein " +
                                "identifications"),
                        new FileToReadArgumentDefinition("pepxml", false,
                                "PepXML file containing peptide identifications.  If absent, will look in " +
                                        "ProtXML file for location"),

                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "Directory with mzXML " +
                                "files from the runs that generated the database search results"),

                        new DecimalArgumentDefinition("minproteinprophet", false,
                                "Minimum ProteinProphet group probability for proteins",
                                minProteinProphet),
                        new FileToReadArgumentDefinition("protgenefile", false,
                                "Tab-delimited file associating gene symbols with protein accession numbers"),
                        new StringListArgumentDefinition("proteins", false,
                                "Protein(s) whose events you wish to survey.  If specifying multiple proteins, " +
                                        "separate names with ','.  Leave blank for a table of all proteins.)"),
                        new FileToReadArgumentDefinition("proteinfile", false,
                                                        "File containing protein accessions, one per line"),
                        new DecimalArgumentDefinition("minhighratio", false,
                                "Ratios must be higher than this, or lower than maxratio, or both",
                                minHighRatio),
                        new DecimalArgumentDefinition("maxlowratio", false,
                                "Ratios must be lower than this, or higher than minratio, or both",
                                maxLowRatio),
                        new StringListArgumentDefinition("genes", false,
                                     "Gene(s) whose events you wish to survey (requires 'protgenefile', " +
                                     "can't be used with 'proteins').  If specifying multiple proteins, " +
                                     "separate names with ','.  Leave blank for a table of all proteins.)"),
                };
        addArgumentDefinitions(argDefs);
        if (shouldCaptureOutputArgs)
        {
            addArgumentDefinition(new DirectoryToWriteArgumentDefinition("outdir", false,
                                "Base output directory for charts (protein-specific charts will be created in " +
                                        "protein-specific subdirectories)"));
            addArgumentDefinition(new FileToWriteArgumentDefinition("out", false,
                                "Output .tsv file location (if blank, output will be written to a temporary file)"));
            addArgumentDefinition(new BooleanArgumentDefinition("appendoutput", false,
                                "Append output to a file, if that file already exists? (otherwise, remove " +
                                        "existing file)", appendOutput), true);
        }

    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzXmlDir = getFileArgumentValue("mzxmldir");
        protXmlFile = getFileArgumentValue("protxml");
        pepXmlFile = getFileArgumentValue("pepxml");
        if (pepXmlFile == null)
        {
            try
            {
                ApplicationContext.infoMessage("Finding source PepXML file in ProtXML file " +
                        protXmlFile.getAbsolutePath() + "...");
                List<File> pepXmlFiles = ProteinUtilities.findSourcePepXMLFiles(protXmlFile);
                if (pepXmlFiles.size() > 1)
                    throw new ArgumentValidationException("Multiple PepXML files specified in ProtXML file " +
                            protXmlFile.getAbsolutePath() +
                            ".  Multiple PepXML files per ProtXML file are not currently supported by Qurate.");
                pepXmlFile = pepXmlFiles.get(0);
                ApplicationContext.infoMessage("Located PepXML file " + pepXmlFile.getAbsolutePath());
            }
            catch (FileNotFoundException e)
            {
                throw new ArgumentValidationException("Can't open PepXML file specified in ProtXML file " +
                        protXmlFile.getAbsolutePath());
            }
            catch (XMLStreamException e)
            {
                throw new ArgumentValidationException("Can't open PepXML file specified in ProtXML file " +
                        protXmlFile.getAbsolutePath());
            }
        }

        if (hasArgumentValue("outdir"))
            outDir = getFileArgumentValue("outdir");
        if (hasArgumentValue("appendoutput"))
            appendOutput = getBooleanArgumentValue("appendoutput");

        if (hasArgumentValue("genes"))
        {
            assertArgumentPresent("protgenefile", "genes");
            assertArgumentAbsent("proteins", "genes");
            assertArgumentAbsent("proteinfile", "genes");

        }

        List<String> proteinNames = (List<String>) getArgumentValue("proteins");


        if ((proteinNames == null || proteinNames.isEmpty()) && hasArgumentValue("proteinfile"))
        {
            proteinNames = new ArrayList<String>();
            File proteinFile = getFileArgumentValue("proteinfile");
            try {
                BufferedReader br = new BufferedReader(new FileReader(proteinFile));
                while (true) {
                    String line = br.readLine();
                    if (line == null)
                        break;
                    proteinNames.add(line);
                }
            } catch(Exception e) {
                throw new ArgumentValidationException("Failed to read proteinfile.",e);
            }
        }

        if (proteinNames == null || proteinNames.isEmpty())
        {
            ApplicationContext.infoMessage("No protein names provided, searching all proteins");
        }
        else
        {
        }

        minProteinProphet = getFloatArgumentValue("minproteinprophet");
        if (hasArgumentValue("out"))
            outFile = getFileArgumentValue("out");

        File protGeneFile = getFileArgumentValue("protgenefile");
        if (protGeneFile != null)
        {
            try
            {
                proteinGeneListMap = QAUtilities.loadIpiGeneListMap(protGeneFile);
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Failed to load protein-gene map file",e);
            }
        }

        List<String> geneNames = (List<String>) getArgumentValue("genes");
        if (geneNames != null && !geneNames.isEmpty())
        {
            //todo: this is inefficient -- reopens the file, rather than using the earlier map to construct this one
            try
            {
                ApplicationContext.infoMessage("Loading gene-protein map....");
                Map<String,List<String>> geneIpiListMap = QAUtilities.loadGeneIpiListMap(protGeneFile);
                Set<String> proteinNamesTheseGenes = new HashSet<String>();
                for (String geneName : geneNames)
                {
                    List<String> proteinsThisGene = null;
                    if (geneIpiListMap.containsKey(geneName))
                    {
                        proteinsThisGene = geneIpiListMap.get(geneName);
                    }
                    else
                    {
                        //check different case
                        for (String geneInMap : geneIpiListMap.keySet())
                        {
                            if (geneInMap.equalsIgnoreCase(geneName))
                            {
                                ApplicationContext.infoMessage("NOTE: for gene " + geneName + ", found similar gene " +
                                        geneInMap + " with different case.  Using " + geneInMap);
                                proteinsThisGene = geneIpiListMap.get(geneInMap);
                            }
                        }
                    }
                    if (proteinsThisGene == null)
                    {
                        ApplicationContext.infoMessage("WARNING: No proteins found for gene " + geneName);
                    }
                    else
                    {
                        proteinNamesTheseGenes.addAll(proteinsThisGene);
                        for (String protein : proteinsThisGene)
                        {
                            ApplicationContext.infoMessage("Using protein " + protein + " for gene " + geneName);
                        }
                    }
                }
                proteinNames = new ArrayList<String>(proteinNamesTheseGenes);

            }
            catch (IOException e)
            {
                throw new ArgumentValidationException("Failed to load protein-gene map file",e);
            }
        }

        if (proteinNames != null && !proteinNames.isEmpty())
            loadProteins(proteinNames);        

        minHighRatio = getFloatArgumentValue("minhighratio");
        maxLowRatio = getFloatArgumentValue("maxlowratio");


    }

    protected void loadProteins(List<String> proteinNames) throws ArgumentValidationException
    {
        proteins = new ArrayList<ProtXmlReader.Protein>();
        try
        {
            //todo: this will only check the first protein occurrence for quantitation.  We could be
            //throwing away quantitated proteins that are quantitated in another occurrence
            Map<String, ProtXmlReader.Protein> proteinNameProteinMap =
                    ProteinUtilities.loadFirstProteinOccurrence(protXmlFile, proteinNames, minProteinProphet);
            List<String> notFoundProteinNames = new ArrayList<String>();
            List<String> unquantitatedProteinNames = new ArrayList<String>();
            List<String> okProteinNames = new ArrayList<String>();
            for (String proteinName : proteinNames)
            {
                if (!proteinNameProteinMap.containsKey(proteinName))
                    notFoundProteinNames.add(proteinName);
                else
                {
                    if (proteinNameProteinMap.get(proteinName).getQuantitationRatio() == null)
                        unquantitatedProteinNames.add(proteinName);
                    else okProteinNames.add(proteinName);
                }
            }
            if (okProteinNames.isEmpty())
                throw new ArgumentValidationException("None of the specified proteins were found, and quantitated, in the protXml file " +
                        protXmlFile.getAbsolutePath() + " with probability >= " + minProteinProphet);
            if (!notFoundProteinNames.isEmpty())
            {
                StringBuffer message = new StringBuffer("WARNING! Some specified proteins were not found" +
                        " in file " +
                        protXmlFile.getAbsolutePath()  + " with probability >= " + minProteinProphet + ".  Missing protein(s): ");
                for (String proteinName : notFoundProteinNames)
                    message.append(" " + proteinName);
                QuantitationReviewer.infoMessage(message.toString());
            }
            if (!unquantitatedProteinNames.isEmpty())
            {
                StringBuffer message = new StringBuffer("WARNING!!! Some specified proteins were not quantitated in file " +
                        protXmlFile.getAbsolutePath() + " with probability >= " + minProteinProphet + ".  Unquantitated protein(s): ");
                for (String proteinName : unquantitatedProteinNames)
                    message.append(" " + proteinName);
                QuantitationReviewer.infoMessage(message.toString());
            }


            for (String proteinName : okProteinNames)
            {
                ProtXmlReader.Protein thisProtein = proteinNameProteinMap.get(proteinName);
                if (proteins.contains(thisProtein))
                {
                    ApplicationContext.infoMessage("NOTE: Indistinguishable protein " +
                            thisProtein.getProteinName() + " used for protein " + proteinName);
                }
                else
                {
                    String origProteinName = thisProtein.getProteinName();
                    if (!origProteinName.equals(proteinName))
                    {
                        thisProtein.setProteinName(proteinName);
                        ApplicationContext.infoMessage("NOTE: Using specified protein name " + proteinName +
                                " instead of original name " + origProteinName + " for indistinguishable protein");
                    }
                    proteins.add(thisProtein);
                }

            }
            System.err.println("Done figuring out proteins.");
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Error loading proteins from file " +
                    protXmlFile.getAbsolutePath(), e);
        }
        System.err.println("Done parsing args.");
        if (proteins != null)
            System.err.println("Proteins to manage: " + proteins.size());

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        System.err.println("execute 1");
        final QuantitationReviewer quantReviewer = new QuantitationReviewer(true, false);
        System.err.println("execute 2");
        quantReviewer.settingsCLM = this;
        if (proteins != null)
        {
            System.err.println("Showing quant summary frame for " + proteins.size() + "  proteins...");

            quantReviewer.showProteinQuantSummaryFrame(proteins, proteinGeneListMap);
        }
        else
        {
            try
            {
                final ProteinSummarySelectorFrame proteinSummarySelector = new ProteinSummarySelectorFrame();
                proteinSummarySelector.setMinProteinProphet(minProteinProphet);
                proteinSummarySelector.setMinHighRatio(minHighRatio);
                proteinSummarySelector.setMaxLowRatio(maxLowRatio);


                proteinSummarySelector.setProteinGeneMap(proteinGeneListMap);
                proteinSummarySelector.addSelectionListener(
                        new ActionListener()
                        {
                            public void actionPerformed(ActionEvent e)
                            {
                                if (proteinSummarySelector.getSelectedProteins() != null &&
                                        !proteinSummarySelector.getSelectedProteins().isEmpty())
                                {
                                    quantReviewer.showProteinQuantSummaryFrame(proteinSummarySelector.getSelectedProteins());
                                }
                            }
                        }
                );

                proteinSummarySelector.displayProteins(protXmlFile);
                proteinSummarySelector.setVisible(true);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Error opening ProtXML file " +
                        protXmlFile.getAbsolutePath(),e);
            }
        }
    }

    public Map<String, List<String>> getProteinGeneListMap()
    {
        return proteinGeneListMap;
    }

    public void setProteinGeneListMap(Map<String, List<String>> proteinGeneListMap)
    {
        this.proteinGeneListMap = proteinGeneListMap;
    }
}
