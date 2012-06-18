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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.quant.gui.QuantitationReviewer;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.awt.*;


/**
 * test
 */
public class ReviewQuantitationCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ReviewQuantitationCLM.class);

    protected File quantSummaryFile;

    protected boolean done = false;

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
                        this.createUnnamedFileArgumentDefinition(false, "Quantitation summary file"),
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
            if (quantSummaryFile != null)
                quantReviewer = new QuantitationReviewer(quantSummaryFile);
            else
            {
                quantReviewer = new QuantitationReviewer();
            }
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to open quantitation summary file", e);
        }

        quantReviewer.setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        quantReviewer.setVisible(true);
    }


}
