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

import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithProteinRatiosChart;
import org.fhcrc.cpl.viewer.quant.gui.ProteinQuantSummaryFrame;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;


/**
 * test
 */
public class ProteinQuantChartsCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProteinQuantChartsCLM.class);

    protected File protXmlFile;
    protected File pepXmlFile;
    protected String proteinName;
    protected File outDir;
    protected File mzXmlDir;

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
                        this.createFileToReadArgumentDefinition("protxml", true, "ProtXML file"),
                        this.createFileToReadArgumentDefinition("pepxml", true, "PepXML file"),
                        this.createStringArgumentDefinition("protein", true, "Protein name"),
                        this.createDirectoryToReadArgumentDefinition("outdir", true, "Output Directory"),
                        this.createDirectoryToReadArgumentDefinition("mzxmldir", true, "Directory with mzXML files"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFile = getFileArgumentValue("protxml");
        pepXmlFile = getFileArgumentValue("pepxml");
        proteinName = getStringArgumentValue("protein");
        outDir = getFileArgumentValue("outdir");
        mzXmlDir = getFileArgumentValue("mzxmldir");         
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
         quantSummaryFrame =
                new ProteinQuantSummaryFrame(protXmlFile, pepXmlFile, proteinName, outDir, mzXmlDir);
//        Thread doneCheckThread = new Thread(new Runnable()
//        {
//            public void run()
//            {
//                while(true)
//                {
//  System.err.println("@@@Checking done, " + isDone());
//
//
//            }
//        });
        while(true)
        {
                    if (isDone())
                    {
                        return;
                    }
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

    public boolean isDone()
    {
        if (quantSummaryFrame == null)
            return false;
        return quantSummaryFrame.isDone();
    }
}
