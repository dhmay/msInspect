package org.fhcrc.cpl.toolbox.chem;

/**
 * Represents a chemical modification that can happen to an Adduct.
 *
 * It happens to an Adduct, rather than to a Compound, because an Adduct is a more general case -- an
 * Adduct with no modifications is equivalent to its Compound.
*/
public interface ChemicalModification
{
    /**
     * Perform the modification on the adduct.  To be absolutely clear, adduct will be altered
     * @param adduct
     * @return
     */
    public void perform(Adduct adduct);

    /**
     * Determines whether the modification can be performed on a particular Adduct
     * @param adduct
     * @return
     */
    public boolean canPerform(Adduct adduct);

    /**
     * Return a symbol to represent this modification, e.g., +2H
     * @return
     */
    public String getSymbol();

    /**
     * A human-readable name
     * @return
     */
    public String getName();
}
