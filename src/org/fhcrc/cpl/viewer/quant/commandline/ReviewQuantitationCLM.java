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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.toolbox.TabLoader;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.quant.QuantitationVisualizer;
import org.fhcrc.cpl.viewer.quant.QuantitationReviewer;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.gui.util.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;


/**
 * test
 */
public class ReviewQuantitationCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ReviewQuantitationCLM.class);

    protected File quantSummaryFile;

    public ReviewQuantitationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "qurate";

        mHelpMessage ="Review quantitation charts";
        mShortDescription = "Review quantitation charts";

        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Quantitation summary file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        quantSummaryFile = this.getUnnamedFileArgumentValue();
    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        QuantitationReviewer quantReviewer = null;
        try
        {
            quantReviewer = new QuantitationReviewer(quantSummaryFile);

        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to open quantitation summary file",e);
        }

        quantReviewer.setTitle("Quantitation");
        quantReviewer.setVisible(true);


    }


}
