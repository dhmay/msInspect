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
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.chem.ChemCalcs;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteUtilities;
import org.apache.log4j.Logger;
import org.openscience.cdk.reaction.type.*;
import org.openscience.cdk.reaction.type.parameters.IParameterReact;
import org.openscience.cdk.reaction.type.parameters.SetReactionCenter;
import org.openscience.cdk.reaction.ReactionEngine;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.dict.DictionaryDatabase;

import java.io.*;
import java.util.*;

/**
 */
public class AnalyzeFormulasCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AnalyzeFormulasCLM.class);

    protected File file;

    protected boolean shouldCollapseByMax = false;

    protected String compoundClass = null;

    public AnalyzeFormulasCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "analyzeformulas";
        mShortDescription = "Analyze chemical formulas";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Input file"),
                        new StringArgumentDefinition("class", false, "Class of compounds to consider")
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        file = this.getUnnamedFileArgumentValue();
        compoundClass = getStringArgumentValue("class");
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

        if (compoundClass != null)
        {
            List<ChemicalCompound> newCompounds = new ArrayList<ChemicalCompound>();
            for (ChemicalCompound compound: compounds)
            {
                if (compound.getCompoundClass().equals(compoundClass))
                     newCompounds.add(compound);
            }
            compounds = newCompounds;
        }

        ApplicationContext.setMessage("Loaded " + compounds.size() + " compounds");
        List<Float> log12PeakRatios = new ArrayList<Float>();
        List<Float> masses = new ArrayList<Float>();
        List<Float> unitMassDefects = new ArrayList<Float>();
        List<Float> totalMassDefects = new ArrayList<Float>();
        List<Float> totalMassOffsets = new ArrayList<Float>();
        List<Float> unitMassOffsets = new ArrayList<Float>();

        List<IParameterReact> paramList1 = new ArrayList<IParameterReact>();
        //Reaction is not localized
        IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList1.add(param);
      

        IReactionProcess reduceDoubleBondReaction = new AdductionProtonPBReaction();
        IReactionProcess addHydrogenReaction = new AdductionProtonLPReaction();
        List<IReactionProcess> allReactions = new ArrayList<IReactionProcess>();
        allReactions.add(reduceDoubleBondReaction);
//        allReactions.add(addHydrogenReaction);


        IMolecule oneHydrogenMol = new Molecule();
        oneHydrogenMol.addAtom(new Atom("H"));

        List<IMolecule> otherReactants = new ArrayList<IMolecule>();
//        otherReactants.add(oneHydrogenMol);

        try
        {
            for (IReactionProcess reaction : allReactions)
                reaction.setParameterList(paramList1);
        }
        catch (CDKException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        for (ChemicalCompound compound : compounds)
        {
            try
            {
                IMolecularFormula cdkFormula = MolecularFormulaManipulator.getMolecularFormula(compound.getCdkMolecule());
ApplicationContext.infoMessage(compound.toString() + ", " + compound.getCommonestIsotopeMass() + ", " +
                        MolecularFormulaManipulator.getString(cdkFormula) + ", " +
                        MolecularFormulaManipulator.getMajorIsotopeMass(cdkFormula));
                List<ChemicalCompound> products = compound.applyReactions(allReactions, otherReactants);
System.err.println("\tProducts: " + products.size());
                for (ChemicalCompound productCompound : products)
                {
                    ApplicationContext.infoMessage("\t\t" + productCompound.getFormula());
                }
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }


            if (compound.getPeakMasses().length > 1)
            {
                float mass = (float)compound.getCommonestIsotopeMass();

                if (compound.getPeakFrequencies()[0] == 0 || compound.getPeakFrequencies()[1] == 0)
                    continue;                   

                float peakRatio = (float) (compound.getPeakFrequencies()[0] / compound.getPeakFrequencies()[1]);
                float peakLogRatio = (float) Math.log(peakRatio);

                if (peakLogRatio < 0) continue;
                if (mass < 100 && peakLogRatio < 1.5)
                            continue;
                if (mass < 50 || mass > 1000)
                            continue;
//                System.err.println(peakLogRatio);
                if (Float.isInfinite(peakLogRatio)) System.err.println(compound.getPeakFrequencies()[0] + "," + compound.getPeakFrequencies()[1]);
                log12PeakRatios.add(peakLogRatio);
                masses.add(mass);
                float totalMassDefect = (float) compound.getFormula().calcTotalMassDefect();
                float unitMassDefect = (float) compound.getFormula().calcUnitMassDefect();
                unitMassDefects.add(unitMassDefect);
                totalMassDefects.add(totalMassDefect);
                float totalMassOffset = (float) MetaboliteUtilities.calcTotalMassDefectOffset(mass);
                totalMassOffsets.add(totalMassOffset);
                unitMassOffsets.add((float) MetaboliteUtilities.calcUnitMassDefectOffset(mass));
            }
        }
        new PanelWithHistogram(log12PeakRatios, "log peak ratios").displayInTab();
        new PanelWithHistogram(unitMassDefects, "unit defects").displayInTab();
        System.err.println("Mean unit mass defect: " + BasicStatistics.mean(unitMassDefects));

        new PanelWithScatterPlot(masses, unitMassDefects, "mass vs unit mass defect").displayInTab();
        new PanelWithScatterPlot(masses, totalMassDefects, "mass vs total mass defect").displayInTab();
        new PanelWithScatterPlot(masses, totalMassOffsets, "mass vs total mass offset").displayInTab();
        new PanelWithHistogram(totalMassOffsets, "total mass offsets").displayInTab();
        new PanelWithHistogram(unitMassOffsets, "unit mass offsets").displayInTab();
        new PanelWithScatterPlot(masses, unitMassOffsets, "mass vs unit mass offset").displayInTab();





        PanelWithScatterPlot pwsp = new PanelWithScatterPlot(masses, log12PeakRatios, "mass v logratio");

        try
        {
            double[] coeffs = RegressionUtilities.modalRegression(masses,
                    log12PeakRatios, 4);

            pwsp.addLineOrCurve(coeffs, BasicStatistics.min(masses), BasicStatistics.max(masses));
            pwsp.displayInTab();
            System.err.print("Coefficients: ");
            for (int i=0; i<coeffs.length; i++)
                System.err.print(", " + coeffs[i]);
            System.err.println();

            List<Double> logRatioResiduals = new ArrayList<Double>();
            for (int i=0; i<masses.size(); i++)
            {
                 logRatioResiduals.add(log12PeakRatios.get(i) - RegressionUtilities.mapValueUsingCoefficients(coeffs, masses.get(i)));
            }
            double sd = BasicStatistics.standardDeviation(logRatioResiduals);
            System.err.println("SD residual: " + sd);
            coeffs[0] += sd;
            pwsp.addLineOrCurve(coeffs, BasicStatistics.min(masses), BasicStatistics.max(masses));
            coeffs[0] -= 2*sd;
            pwsp.addLineOrCurve(coeffs, BasicStatistics.min(masses), BasicStatistics.max(masses));

            new PanelWithHistogram(logRatioResiduals, "SD errors").displayInTab();
        }
        catch (Exception e) {
            e.printStackTrace(System.err);  
            throw new CommandLineModuleExecutionException(e);
        }

//        PanelWithScatterPlot pwsp = new PanelWithScatterPlot();
//        pwsp.setName("mass vs " + element + "s");
//        pwsp.addData(massesInRange, numsInRange, "in range");
//        pwsp.addData(masses, numsOfElements, "mass vs " + element + "s");
//        pwsp.displayInTab();
//        new PanelWithHistogram(numsOfElements, "num " + element).displayInTab();
    }
}


