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
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteUtilities;
import org.fhcrc.cpl.viewer.metabologna.MoleculeRenderer2D;
import org.apache.log4j.Logger;
import org.openscience.cdk.reaction.type.*;
import org.openscience.cdk.reaction.type.parameters.IParameterReact;
import org.openscience.cdk.reaction.type.parameters.SetReactionCenter;
import org.openscience.cdk.reaction.ReactionEngine;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.dict.DictionaryDatabase;

import java.io.*;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;

/**
 */
public class AddReactionProductsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AddReactionProductsCLM.class);

    protected File file;
    protected File outFile;


    public AddReactionProductsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "addreactionproducts";
        mShortDescription = "Add reaction products to a database";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Input file"),
                        new FileToWriteArgumentDefinition("out", true, "Output file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        file = this.getUnnamedFileArgumentValue();
        outFile = getFileArgumentValue("out");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        List<ChemicalCompound> compounds = null;
        try
        {
            compounds = ChemicalCompound.loadCompoundsFromFile(file, 2);
//    System.err.println(formula + ", " + numOfElement);
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to load file " + file.getAbsolutePath(),e);
        }
    }


    public IMolecule createCOOH()
    {
        IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
        IMolecule molecule = builder.newMolecule();
        molecule.addAtom(builder.newAtom("C"));
        molecule.addAtom(builder.newAtom("O"));
        molecule.addBond(0, 1, IBond.Order.DOUBLE);
        molecule.addAtom(builder.newAtom("O"));
        molecule.addBond(0, 2, IBond.Order.SINGLE);

        try {

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
            LonePairElectronChecker lpcheck = new LonePairElectronChecker();
            lpcheck.saturate(molecule);

        }
        catch (CDKException e)
        {
            e.printStackTrace();
        }
        return molecule;
    }
}


