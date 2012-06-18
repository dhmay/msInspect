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


/**
 * Represents a simple chemical subtraction that can happen to an Adduct.  Just subtracts the specified ChemicalFormula
 * atoms from the Adduct.  canPerform returns false if the atoms aren't there
*/
public class SimpleChemicalSubtraction implements ChemicalModification
{
    protected ChemicalFormula formulaToSubtract;

    public SimpleChemicalSubtraction(ChemicalFormula formulaToSubtract)
    {
        this.formulaToSubtract = formulaToSubtract;
    }

    /**
     * Creates a new Adduct with the modification performed
     * @param adduct
     * @return
     */
    public void perform(Adduct adduct)
    {
        adduct.setFormula(adduct.getFormula().createFormulaWithSubtraction(formulaToSubtract));
        adduct.modifications.add(this);
    }

    /**
     * We can only subtract if the elements are present
     * @param adduct
     * @return
     */
    public boolean canPerform(Adduct adduct)
    {
        try
        {
            adduct.getFormula().createFormulaWithSubtraction(formulaToSubtract);
        }
        catch (IllegalArgumentException e)
        {
            return false;
        }
        return true;
    }

    public String getSymbol()
    {
        return "-" + formulaToSubtract;
    }

    public String getName()
    {
        return "Subtraction of " + formulaToSubtract;
    }
}
