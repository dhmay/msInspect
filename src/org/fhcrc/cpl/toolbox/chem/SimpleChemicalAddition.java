package org.fhcrc.cpl.toolbox.chem;


/**
 * Represents a simple chemical addition that can happen to an Adduct.  Just adds the specified ChemicalFormula
 * atoms to the Adduct.  Does no checking to see if this is reasonable
*/
public class SimpleChemicalAddition implements ChemicalModification
{
    protected ChemicalFormula formulaToAdd;

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
