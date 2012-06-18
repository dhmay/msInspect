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
package org.fhcrc.cpl.viewer.metabologna.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.chem.ChemCalcs;
import org.fhcrc.cpl.toolbox.chem.ChemicalFormula;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteDatabaseMatcher;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 */
public class GuessFormulasForMassCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(GuessFormulasForMassCLM.class);

    protected float massTolerancePPM = 2f;
    protected float mass;


    protected List<ChemicalCompound> databaseCompoundsByMass = null;

    protected TabLoader loader;

    protected boolean shouldCollapseByMax = false;

    public GuessFormulasForMassCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "guessformulasformass";
        mShortDescription = "guessformulasformass";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        new DecimalArgumentDefinition("mass", true, "mass"),
                        new DecimalArgumentDefinition("deltappm", false, "Delta mass (ppm)", massTolerancePPM),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        massTolerancePPM = getFloatArgumentValue("deltappm");
        mass = getFloatArgumentValue("mass");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            List<ChemicalFormula> formulas =
                    ChemCalcs.CDKFormulaSet2ChemFormList(ChemCalcs.calcMass2Formulas(mass, massTolerancePPM));
            ApplicationContext.infoMessage("Found " + formulas.size() + " formulas");
            for (ChemicalFormula formula : formulas)
            {
                System.err.println(formula + ", mass=" + formula.getCommonestIsotopeMass() + ", deltaPPM="
                    + MassUtilities.convertDaToPPM((float) formula.getCommonestIsotopeMass() - mass, mass));
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }



    }
}


