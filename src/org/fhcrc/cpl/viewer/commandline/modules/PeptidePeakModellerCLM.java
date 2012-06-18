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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.DecimalArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBarChart;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;


/**
 * Show a plot that lets the user determine the mass accuracy of the instrumentation
 */
public class PeptidePeakModellerCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptidePeakModellerCLM.class);

    protected float monoIsotopicMass = 0;



    public PeptidePeakModellerCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "modelpeptide";
        mShortDescription = "";
        mUsageMessage = CommandLineModule.MODULE_USAGE_AUTOMATIC;
        mHelpMessage =
                "";

        CommandLineArgumentDefinition[] argDefs =
            {
                    new DecimalArgumentDefinition("daltons", true, "peptide mass"),
            };
        addArgumentDefinitions(argDefs);
    }



    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        monoIsotopicMass = getFloatArgumentValue("daltons");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        float[] intensities = Spectrum.Poisson(monoIsotopicMass);
        String[] masses = new String[intensities.length];
        for (int i=0; i<masses.length; i++)
            masses[i] = "" + (monoIsotopicMass + i);
        PanelWithBarChart chartPanel = new PanelWithBarChart(masses, intensities, "Peak Intensities");

        ChartDialog chartDialog = new ChartDialog(chartPanel);
        chartDialog.setVisible(true);

        ApplicationContext.infoMessage("Mass\tIntensity");
        for (int i=0; i<masses.length; i++)
            ApplicationContext.infoMessage(masses[i] + "\t" + intensities[i]);
    }

}
