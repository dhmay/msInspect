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
import org.fhcrc.cpl.viewer.qa.QAUtilities;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.File;


/**
 * test
 */
public class ProteinQuantChartsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProteinQuantChartsCLM.class);

    protected File protXmlFile;
    protected File pepXmlFile;
    protected String proteinName;
    protected File outDir;
    protected File outFile;
    protected File mzXmlDir;
    protected Boolean appendOutput = true;
    protected float minProteinProphet = 0.9f;
    protected File protGeneFile;

    protected ProteinQuantSummaryFrame quantSummaryFrame;

    public ProteinQuantChartsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "proteinquantcharts";

        mHelpMessage ="proteinquantcharts";
        mShortDescription = "proteinquantcharts";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("protxml", true, "ProtXML file"),
                        new FileToReadArgumentDefinition("pepxml", true, "PepXML file"),
                        new StringArgumentDefinition("protein", false, "Protein name"),
                        new DirectoryToReadArgumentDefinition("outdir", true, "Output Directory"),
                        new FileToWriteArgumentDefinition("out", false, "Output File"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "Directory with mzXML files"),
                        new BooleanArgumentDefinition("appendoutput", false,
                                "Append output to file, if already exists?", appendOutput),
                        new DecimalArgumentDefinition("minproteinprophet", false,
                                "Minimum ProteinProphet score for proteins (if protein not specified)",
                                minProteinProphet),
                       new FileToReadArgumentDefinition("protgenefile", false,
                               "File associating gene symbols with protein accession numbers"),
                };

        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzXmlDir = getFileArgumentValue("mzxmldir");
        protXmlFile = getFileArgumentValue("protxml");
        pepXmlFile = getFileArgumentValue("pepxml");
        proteinName = getStringArgumentValue("protein");

        outDir = getFileArgumentValue("outdir");
        outFile = getFileArgumentValue("out");
        appendOutput = getBooleanArgumentValue("appendoutput");

        minProteinProphet = getFloatArgumentValue("minproteinprophet");

        protGeneFile = getFileArgumentValue("protgenefile");

    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (proteinName == null)
        {
            try
            {
                final ProteinSummarySelectorFrame proteinSummarySelector = new ProteinSummarySelectorFrame(false);
                proteinSummarySelector.setMinProteinProphet(minProteinProphet);
                if (protGeneFile != null)
                    proteinSummarySelector.setProteinGeneMap(QAUtilities.loadIpiGeneListMap(protGeneFile));
                proteinSummarySelector.displayProteins(protXmlFile);
                proteinSummarySelector.setVisible(true);
                //dialog is modal, so next lines happen after close
                if (proteinSummarySelector.getSelectedProteins() == null || 
                        proteinSummarySelector.getSelectedProteins().isEmpty())
                    return;
                proteinName = proteinSummarySelector.getSelectedProteins().get(0).getProteinName();
                proteinSummarySelector.dispose();
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Error opening ProtXML file " + protXmlFile.getAbsolutePath(),e);
            }
        }
        if (proteinName == null)
        {
            ApplicationContext.infoMessage("No protein selected");
            return;
        }

        quantSummaryFrame = new ProteinQuantSummaryFrame(outDir, mzXmlDir, outFile, appendOutput);

        while(true)
        {
            if (!quantSummaryFrame.isVisible())
                return;
            try
            {
                Thread.sleep(500);
            }
            catch (InterruptedException e)
            {
                return;
            }
        }

    }

    public ProteinQuantSummaryFrame getQuantSummaryFrame()
    {
        return quantSummaryFrame;
    }

}
