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


import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

/**
 * Represents a chemical compound with optional ChemicalModifications.
* User: dhmay
* Date: Apr 6, 2010
* Time: 4:06:07 PM
* To change this template use File | Settings | File Templates.
*/
public class Adduct
{
    protected ChemicalFormula formula;

    //None of these variables should ever be set without using the setters.  setting any of these with the
    //setters will update formula

    protected ChemicalCompound compound;

    //Molecule structure that's correct for this adduct.  Hopefully we can always populate this
    protected IMolecule molecule;

    //These are guaranteed never to be null unless explicitly set that way by pathological code
    protected List<ChemicalModification> modifications = new ArrayList<ChemicalModification>();

    /**
     * Create a fully-defined Adduct based on the supplied compound, with modifications supplied to be applied in order
     * @param compound
     * @param modifications
     *
     * @throws IllegalArgumentException if the modifications can't be performed
     */
    public Adduct(ChemicalCompound compound,
                  List<ChemicalModification> modifications)
            throws IllegalArgumentException
    {
        setCompound(compound);
        if (compound.getCdkMolecule() != null)
        {
            try
            {
                setMolecule((IMolecule) compound.getCdkMolecule().clone());
            }
            catch (CloneNotSupportedException e)
            {
                //cloning /is/ supported
            }
        }

        if (modifications != null)
            for (ChemicalModification mod : modifications)
            {
                if (!mod.canPerform(this))
                    throw new IllegalArgumentException("Can't perform modification " + mod.getSymbol());
                mod.perform(this);
            }
    }

    /**
     * Takes a slightly cheaper route of simply copying the compound, formula and modifications from the other
     * adduct, rather than applying all the modifications one by one
     * @param adduct
     */
    public Adduct(Adduct adduct)
    {
        setCompound(adduct.getCompound());
        this.formula = adduct.getFormula();
        this.modifications = new ArrayList<ChemicalModification>(modifications);
        try
        {            
            setMolecule((IMolecule) compound.getCdkMolecule().clone());
        }
        catch (Exception e)
        {
            //cloning /is/ supported
        }
    }

    /**
     * Create an adduct equivalent to the compound with no additions or subtractions
     * @param compound
     */
    public Adduct(ChemicalCompound compound)
    {
        this(compound, null);
    }

    public ChemicalCompound getCompound()
    {
        return compound;
    }

    public void setCompound(ChemicalCompound compound)
    {
        this.compound = compound;
        formula = new ChemicalFormula(compound.getFormula());
    }

    //convenience methods for getting at a few things in ChemicalFormula

    public double getCommonestIsotopeMass()
    {
        return formula.getCommonestIsotopeMass();
    }

    public double[] getPeakFrequencies()
    {
        return formula.getPeakFrequencies();
    }

    public double[] getPeakMasses()
    {
        return formula.getPeakMasses();
    }

    public ChemicalFormula getFormula() {
        return formula;
    }

    /**
     * Use with EXTREME CAUTION, because it can make formula disagree with molecule.
     * This is used by, e.g., ChemicalModification
     * @param formula
     */
    public void setFormula(ChemicalFormula formula) {
        this.formula = formula;
    }

    public String getCompoundNameAndIonTypeString()
    {
        return compound.getName() + ":" + getIonTypeString();
    }

    /**
     * Returns a string that describes the modifications done to the base compound to produce this adduct
     * @return
     */
    public String getIonTypeString()
    {
        if (modifications.isEmpty())
            return "[M]";
        StringBuffer resultBuf = new StringBuffer("[");
        boolean first = true;
        for (ChemicalModification mod : modifications)
        {
            if (!first)
                resultBuf.append(" ");
            first=false;
            resultBuf.append(mod.getSymbol());
        }
        resultBuf.append("]");
        return resultBuf.toString();
    }

    public String toString()
    {
        StringBuffer resultBuf = new StringBuffer(compound.getName() + "\t" + formula.toString() + "\t" +
            getIonTypeString());
        return resultBuf.toString();
    }

    public static class ComparatorMassAsc implements Comparator<Adduct>
    {
        public int compare(Adduct o1, Adduct o2)
        {
            if (o1.getCommonestIsotopeMass() == o2.getCommonestIsotopeMass())
                return 0;
            return o1.getCommonestIsotopeMass() == o2.getCommonestIsotopeMass() ? 0 :
                    o1.getCommonestIsotopeMass() < o2.getCommonestIsotopeMass() ? -1 : 1;
        }
    }

    public static class ComparatorNameAndIonTypeAsc implements Comparator<Adduct>
    {
        public int compare(Adduct o1, Adduct o2)
        {
            return o1.getCompoundNameAndIonTypeString().compareTo(o2.getCompoundNameAndIonTypeString());
        }
    }

    public List<ChemicalModification> getModifications() {
        return modifications;
    }

    public IMolecule getMolecule() {
        return molecule;
    }

    public void setMolecule(IMolecule molecule) {
        this.molecule = molecule;
        updateFormula();
    }

    /**
     * Update formula, mass, etc., because molecule has changed
     */
    public void updateFormula()
    {
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);

        formula = ChemCalcs.CDKMolForm2ChemForm(MolecularFormulaManipulator.getMolecularFormula(molecule));
        //This doesn't actually remove hydrogens in the argument, just returns a new IMolecule
        this.molecule = (IMolecule) AtomContainerManipulator.removeHydrogens(molecule);
    }
}
