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
 * Represents a simple chemical addition that can happen to an Adduct.  Just adds the specified ChemicalFormula
 * atoms to the Adduct.  Does no checking to see if this is reasonable
*/
public class SimpleChemicalAddition implements ChemicalModification
{
    protected ChemicalFormula formulaToAdd;

    public SimpleChemicalAddition()
    {
    }

    public SimpleChemicalAddition(ChemicalFormula formulaToAdd)
    {
        this.formulaToAdd = formulaToAdd;
    }

    /**
     * Creates a new Adduct with the modification performed
     * @param adduct
     * @return
     */
    public void perform(Adduct adduct)
    {
        adduct.setFormula(adduct.getFormula().createFormulaWithAddition(formulaToAdd));
        adduct.modifications.add(this);
    }

    /**
     * Assume we can always add stuff
     * @param adduct
     * @return
     */
    public boolean canPerform(Adduct adduct)
    {
         return true;
    }

    public String getSymbol()
    {
        return "+" + formulaToAdd;
    }

    public String getName()
    {
        return "Addition of " + formulaToAdd;
    }
}
