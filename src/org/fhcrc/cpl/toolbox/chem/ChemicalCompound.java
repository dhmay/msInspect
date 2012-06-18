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

package org.fhcrc.cpl.toolbox.chem;

import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.IOException;

/**
 * Represents a named chemical compound.  Contains a chemical formula and a name.  Lots of passing through
 * of methods into the chemical formula
* User: dhmay
* Date: Apr 6, 2010
* Time: 4:06:07 PM
* To change this template use File | Settings | File Templates.
*/
public class ChemicalCompound
{
    protected static Logger _log = Logger.getLogger(ChemicalCompound.class);


    protected String name;
    protected ChemicalFormula formula;


    protected IMolecule cdkMolecule;

    //this will need to get more sophisticated once we figure out what's useful
    protected String compoundClass;

    /**
     * Creates a ChemicalCompound, does not populate peak masses and intensities
     * @param name
     * @param formulaString
     * @throws IllegalArgumentException
     */
    public ChemicalCompound(String name, String formulaString) throws IllegalArgumentException
    {
        this(name, formulaString, 0);
    }

    /**
     * Creates a ChemicalCompound, populates peak masses and intensities as specified
     * @param name
     * @param formulaString
     * @param numPeaksToPopulate
     * @throws IllegalArgumentException
     */
    public ChemicalCompound(String name, String formulaString, int numPeaksToPopulate)
            throws IllegalArgumentException
    {
        this(name, new ChemicalFormula(formulaString, numPeaksToPopulate));
    }

    public static ChemicalCompound createFromSmiles(String name, String smilesString)
            throws CDKException
    {
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule cdkMolecule = smilesParser.parseSmiles(smilesString);
        return new ChemicalCompound(name, cdkMolecule);
    }

    /**
     *
     * @param name
     * @param formula
     */
    public ChemicalCompound(String name, ChemicalFormula formula)
    {
        this.name = name;
        this.formula = formula;
    }

    /**
     * SIDE EFFECT: converts implicit to explicit hydrogens in molecule, then converts them /back/
     * @param name
     * @param molecule
     */
    public ChemicalCompound(String name, IMolecule molecule)
            throws CDKException
    {
        this.name = name;
        this.cdkMolecule = molecule;

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cdkMolecule);


        //When SmilesParser gives us a molecule, the hydrogens are implicit.  If you,
        //e.g., create a MolecularFormula object from the Molecule and ask it for its mass
        //and its formula string, they will be wrong, because they have no H's.
        //So I add the hydrogens, create the formula with the hydrogens, and then remove them
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(cdkMolecule);

        formula = ChemCalcs.CDKMolForm2ChemForm(MolecularFormulaManipulator.getMolecularFormula(cdkMolecule));
        //This doesn't actually remove hydrogens in the argument, just returns a new IMolecule
        cdkMolecule = (IMolecule) AtomContainerManipulator.removeHydrogens(cdkMolecule);
        try
        {
            new LonePairElectronChecker().saturate(cdkMolecule);
        }
        catch (Exception e) //failed to saturate.  Not sure why this happens sometimes on valid SMILES strings
        {}
    }

    /**
     * Create a new ChemicalCompound identical to this one with the additional elements added and the specified name.
     * Do not populate peaks
     * @param additionFormula
     * @return
     */    
    public ChemicalCompound createCompoundWithAddition(ChemicalFormula additionFormula,
                                                      String newCompoundName)
    {
        return new ChemicalCompound(newCompoundName, formula.createFormulaWithAddition(additionFormula));
    }

    /**
     * Create a new ChemicalCompound identical to this one with the specified elements removed and the specified name.  
     * Do not populate peaks  Throw IllegalArgumentException if the formula doesn't have the specified elements
     * @param subtractionFormula
     * @return
     */
    public ChemicalCompound createCompoundWithSubtraction(ChemicalFormula subtractionFormula,
                                                      String newCompoundName)
            throws IllegalArgumentException

    {
        return new ChemicalCompound(newCompoundName, formula.createFormulaWithSubtraction(subtractionFormula));
    }

    //comparators

    public static class ComparatorMassAsc implements Comparator<ChemicalCompound>
    {
        public int compare(ChemicalCompound o1, ChemicalCompound o2)
        {
            if (o1.getCommonestIsotopeMass() == o2.getCommonestIsotopeMass())
                return 0;
            return o1.getCommonestIsotopeMass() == o2.getCommonestIsotopeMass() ? 0 :
                    o1.getCommonestIsotopeMass() < o2.getCommonestIsotopeMass() ? -1 : 1;
        }
    }

    public static class ComparatorNameAsc implements Comparator<ChemicalCompound>
    {
        public int compare(ChemicalCompound o1, ChemicalCompound o2)
        {
            return o1.name.compareTo(o2.getName());
        }
    }

    public String toString()
    {
        StringBuffer buf = new StringBuffer("Compound: " + name + ", formula=" + formula);
        return buf.toString();
    }


    //getters and setters

    public String getCompoundClass()
    {
        return compoundClass;
    }

    public void setCompoundClass(String compoundClass)
    {
        this.compoundClass = compoundClass;
    }

    public double getCommonestIsotopeMass()
    {
        return formula.getCommonestIsotopeMass();
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }

    public ChemicalFormula getFormula()
    {
        return formula;
    }

    public void setFormula(ChemicalFormula formula)
    {
        this.formula = formula;
    }

    public double[] getPeakFrequencies()
    {
        return formula.getPeakFrequencies();
    }

    public double[] getPeakMasses()
    {
        return formula.getPeakMasses();
    }

    public double[] getPeakFrequencies(int numPeaksToCalculate)
    {
        return formula.getPeakFrequencies(numPeaksToCalculate);
    }

    public double[] getPeakMasses(int numPeaksToCalculate)
    {
        return formula.getPeakMasses(numPeaksToCalculate);
    }

    public Map<String, Integer> getElementCountMap()
    {
        return formula.getElementCountMap();
    }


    /**
     * Load chemical compounds from a tsv file containing columns called name and formula (at least)
     * @param file
     * @param numPeaksToPopulate
     * @return
     * @throws IOException
     */
    public static final List<ChemicalCompound> loadCompoundsFromFile(File file, int numPeaksToPopulate)
            throws IOException
    {
        return loadCompoundsFromFile(file, numPeaksToPopulate, "name","formula","smiles");
    }


    /**
     * Load chemical compounds from a tsv file containing columns for name and formula (at least)
     *
     * @param file
     * @param numPeaksToPopulate
     * @param nameColName
     * @param formulaColName
     * @return
     * @throws IOException
     */
    public static List<ChemicalCompound> loadCompoundsFromFile(File file, int numPeaksToPopulate,
                                                                     String nameColName, String formulaColName,
                                                                     String smilesColumnName)
            throws IOException
    {
        TabLoader loader = null;


        Map[] rowsAsMaps = null;

        try {
            loader = new TabLoader(file);
            rowsAsMaps = (Map[])loader.load();
        }
        catch (Exception e) { throw new IOException("Error loading tab-delimited compound file: " + e.getMessage()); }

        boolean foundNameCol = false;
        boolean foundFormulaCol = false;
        boolean foundSmilesCol = false;
        for (TabLoader.ColumnDescriptor col : loader.getColumns()) {
            if (col.name.equals(nameColName))
                foundNameCol = true;
            else if (col.name.equals(formulaColName))
                foundFormulaCol = true;
            else if (col.name.equals(smilesColumnName))
                foundSmilesCol = true;
        }
        if (!foundNameCol || !foundFormulaCol || !foundSmilesCol) {
            throw new IOException("Compounds file is missing at least one column of " + nameColName +
                    ", " + formulaColName + ", " + smilesColumnName);
        }

        List<ChemicalCompound> result = new ArrayList<ChemicalCompound>();
        _log.debug("Loading " + rowsAsMaps.length + " compounds...");
        for (int i=0; i<rowsAsMaps.length; i++)
        {
            Map rowMap = rowsAsMaps[i];
            if (i % 500 == 0)
                _log.debug("\tLoaded " + i + " of " + rowsAsMaps.length + " compounds");
            try
            {
                String name = (String) rowMap.get(nameColName);
                String formulaString = rowMap.get(formulaColName).toString();
                ChemicalCompound compound = null;



                //Load IMolecule using the SMILES formula string
                if (rowMap.containsKey(smilesColumnName) && rowMap.get(smilesColumnName) != null)
                {

                    String smilesString = rowMap.get(smilesColumnName).toString();
                    try
                    {
                        compound = createFromSmiles(name, smilesString);
                    }
                    catch (Exception e)
                    {
                        ApplicationContext.infoMessage("WARNING: Failed to load SMILES formula " + smilesString + " for compound " + name + ": " + e.getMessage());
                        compound = new ChemicalCompound(name,
                                formulaString, numPeaksToPopulate);
//                        throw new IOException(e);
                    }
                }
                else
                {
                    compound = new ChemicalCompound(name,
                            formulaString, numPeaksToPopulate);
                }

                //experimental
                if (rowMap.containsKey("class") && rowMap.get("class") != null)
                {
//                    System.err.println("Class: " + rowMap.get("class").toString());
                    compound.setCompoundClass(rowMap.get("class").toString());
                }

                result.add(compound);
            }
            catch (IllegalArgumentException e)
            {
                ApplicationContext.setMessage("Skipping bad compound: " + e.getMessage());
            }

        }
        _log.debug("Loaded " + rowsAsMaps.length + " compounds.");

        return result;
    }

    /**
     * This is kind of a weird one.  Apply multiple reactions.  After each reaction, apply the following
     * reaction to /all/ the products of the first reaction
     * @param reactions
     * @param otherReactants
     * @return
     * @throws CDKException
     */
    public List<ChemicalCompound> applyReactions(List<IReactionProcess> reactions, List<IMolecule> otherReactants)
            throws CDKException
    {
        List<ChemicalCompound> reactants = new ArrayList<ChemicalCompound>();
        reactants.add(this);

        for (IReactionProcess reaction : reactions)
        {
            List<ChemicalCompound> newReactants = new ArrayList<ChemicalCompound>();
            for (ChemicalCompound reactant : reactants)
            {
                 newReactants.addAll(reactant.applyReaction(reaction, otherReactants));
            }
            reactants = newReactants;
        }
        return reactants;
    }

    /**
     * Apply a chemical reaction to the cdkMolecule of this compound, possibly with other reactants involved.
     *
     * Returns ChemicalCompounds for all possible products
     *
     * todo: do I need to group the products by IReaction? 
     * @param reaction
     * @param otherReactants
     * @return
     * @throws CDKException
     */
    public List<ChemicalCompound> applyReaction(IReactionProcess reaction, List<IMolecule> otherReactants)
            throws CDKException
    {
        IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
        setOfReactants.addMolecule(cdkMolecule);
        if (otherReactants != null)
        {
            for (IMolecule otherReactant : otherReactants)
                setOfReactants.addMolecule(otherReactant);
        }
        IReactionSet setOfReactions = reaction.initiate(setOfReactants, null);

        List<ChemicalCompound> result = new ArrayList<ChemicalCompound>();

        for (IReaction outputReaction : setOfReactions.reactions())
        {
            IMoleculeSet products = outputReaction.getProducts();
            for (int i=0; i<products.getMoleculeCount(); i++)
            {
                result.add(new ChemicalCompound(name, products.getMolecule(i)));
            }
        }
        return result;
    }

    public IMolecule getCdkMolecule()
    {
        return cdkMolecule;
    }

    public void setCdkMolecule(IMolecule cdkMolecule)
    {
        this.cdkMolecule = cdkMolecule;
    }
}
