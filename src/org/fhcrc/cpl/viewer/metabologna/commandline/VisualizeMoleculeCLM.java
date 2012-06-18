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
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MoleculeRenderer2D;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAdd2HMod;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAddWaterMod;
import org.apache.log4j.Logger;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 */
public class VisualizeMoleculeCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(VisualizeMoleculeCLM.class);

    protected String smilesString = null;
    protected int width = 400;
    protected int height = 400;

    protected boolean showHydrogens = false;
    protected boolean showCarbons = false;

    public VisualizeMoleculeCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "visualizemolecule";
        mShortDescription = "Visualize a molecule";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "File with SMILES string"),
                        new BooleanArgumentDefinition("showhydrogens", false, "Show hydrogens?", showHydrogens),
                        new BooleanArgumentDefinition("showcarbons", false, "Show carbons?", showCarbons),

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        showHydrogens = getBooleanArgumentValue("showhydrogens");
        showCarbons = getBooleanArgumentValue("showcarbons");

        File smilesFile = this.getUnnamedFileArgumentValue();

        try
        {
        smilesString = new BufferedReader(new FileReader(smilesFile)).readLine();
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException(e);
        }
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            ChemicalCompound compound = ChemicalCompound.createFromSmiles("dummyname", smilesString);
            if (showHydrogens)
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(compound.getCdkMolecule());

            ApplicationContext.infoMessage("Formula: " + compound.getFormula());
            ApplicationContext.infoMessage("Mass: " + compound.getCommonestIsotopeMass());

            System.err.println("\t" + compound.getFormula());

            MoleculeRenderer2D renderer = new MoleculeRenderer2D();
            renderer.setWidth(width);
            renderer.setHeight(height);
            renderer.setShouldShowHydrogens(showHydrogens);
            renderer.setShouldShowCarbons(showCarbons);


            Image image = renderer.renderMolecule(compound.getCdkMolecule());
            PanelWithBlindImageChart chart = new PanelWithBlindImageChart((BufferedImage) image, smilesString);
            chart.displayInTab();

Adduct adduct = new Adduct(compound);
ReduceDoubleBondAddWaterMod mod = new ReduceDoubleBondAddWaterMod();
if (mod.canPerform(adduct))
{
    mod.perform(adduct);
            if (showHydrogens)
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(adduct.getMolecule());    
            image = renderer.renderMolecule(adduct.getMolecule());
            chart = new PanelWithBlindImageChart((BufferedImage) image, "adduct");
            chart.displayInTab();
    System.err.println("Adduct formula: " + adduct.getFormula());
    System.err.println("Adduct mass: " + adduct.getCommonestIsotopeMass());
}
else
{
    System.err.println("CAN'T PERFORM MODIFICATION!!!");
}

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }
}


