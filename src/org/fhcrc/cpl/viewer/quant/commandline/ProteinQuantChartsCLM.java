/* 
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
import org.fhcrc.cpl.viewer.quant.QuantEvent;
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.FileNotFoundException;
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
    public Map<String, List<String>> proteinGeneListMap;
    public List<ProtXmlReader.Protein> proteins;
    public File outFile;

    protected ProteinQuantSummaryFrame quantSummaryFrame;

    protected boolean hasRun = false;

    public ProteinQuantChartsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "proteinquantcharts";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("protxml", true, "ProtXML file with protein " +
                                "identifications"),
                        new FileToReadArgumentDefinition("pepxml", false,
                                "PepXML file containing peptide identifications.  If absent, will look in " +
                                        "ProtXML file for location"),
                        new DirectoryToReadArgumentDefinition("outdir", true,
                                "Base output directory for charts (protein-specific charts will be created in " +
                                        "protein-specific subdirectories)"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "Directory with mzXML " +
                                "files from the runs that generated the database search results"),
                        new BooleanArgumentDefinition("appendoutput", false,
                                "Append output to a file, if that file already exists? (otherwise, remove " +
                                        "existing file)", appendOutput),
                        new DecimalArgumentDefinition("minproteinprophet", false,
                                "Minimum ProteinProphet group probability for proteins",
                                minProteinProphet),
                        new FileToReadArgumentDefinition("protgenefile", false,
                                "Tab-delimited file associating gene symbols with protein accession numbers"),
                        new StringArgumentDefinition("proteins", false,
                                "Protein(s) whose events you wish to survey.  If specifying multiple proteins, " +
                                        "separate names with ','.  Leave blank for a table of all proteins.)"),
                        new FileToWriteArgumentDefinition("out", false,
                                "Output .tsv file location (if blank, output will be written to a temporary file)"),
                };

        addArgumentDefinitions(argDefs);
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

        outDir = getFileArgumentValue("outdir");
        appendOutput = getBooleanArgumentValue("appendoutput");

        List<String> proteinNames = new ArrayList<String>();
        if (hasArgumentValue("proteins"))
        {
            String[] rawProteins = getStringArgumentValue("proteins").split(",");
            for (String protein : rawProteins)
            {
                String trimmedProtein = protein.trim();
                if (trimmedProtein.length() > 0)
                {

                    proteinNames.add(trimmedProtein);
                    ApplicationContext.infoMessage("\tUsing protein " + trimmedProtein);
                }
            }
            if (proteinNames.isEmpty())
            {
                ApplicationContext.infoMessage("No protein names provided, searching all proteins");
            }
        }

        if (!proteinNames.isEmpty())
        {
            proteins = new ArrayList<ProtXmlReader.Protein>();
            try
            {
                //todo: this will only check the first protein occurrence for quantitation.  We could be
                //throwing away quantitated proteins that are quantitated in another occurrence
                Map<String, ProtXmlReader.Protein> proteinNameProteinMap =
                        ProteinUtilities.loadFirstProteinOccurrence(protXmlFile, proteinNames, minProteinProphet);
                if (proteinNameProteinMap.size() < proteinNames.size())
                {
                    StringBuffer message = new StringBuffer("Not all specified proteins were found in file " +
                            protXmlFile.getAbsolutePath() + ".  Missing protein(s): ");
                    for (String proteinName : proteinNames)
                        if (!proteinNameProteinMap.containsKey(proteinName))
                            message.append(" " + proteinName);
                    throw new ArgumentValidationException(message.toString());
                }
                List<String> unquantitatedProteinNames = new ArrayList<String>();
                for (ProtXmlReader.Protein protein : proteinNameProteinMap.values())
                    if (protein.getQuantitationRatio() == null)
                        unquantitatedProteinNames.add(protein.getProteinName());
                if (!unquantitatedProteinNames.isEmpty())
                {
                    StringBuffer message = new StringBuffer("Some proteins were not quantitated in file " +
                            protXmlFile.getAbsolutePath() + ".  Unquantitated protein(s): ");
                    for (String proteinName : unquantitatedProteinNames)
                        message.append(" " + proteinName);
                    throw new ArgumentValidationException(message.toString());
                }

                for (String proteinName : proteinNames)
                    proteins.add(proteinNameProteinMap.get(proteinName));
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Error loading proteins from file " +
                        protXmlFile.getAbsolutePath(), e);
            }
        }

        minProteinProphet = getFloatArgumentValue("minproteinprophet");
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

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        final QuantitationReviewer quantReviewer = new QuantitationReviewer(false);
        quantReviewer.settingsCLM = this;
        if (proteins != null)
        {
            quantReviewer.showProteinQuantSummaryFrame(proteins);
        }
        else
        {
            try
            {
                final ProteinSummarySelectorFrame proteinSummarySelector = new ProteinSummarySelectorFrame();
                proteinSummarySelector.setMinProteinProphet(minProteinProphet);
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
}
