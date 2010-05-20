package org.fhcrc.cpl.toolbox.chem;

/**
 * Represents a chemical modification that can create a new Adduct based on an existing Adduct, if the
 * input Adduct meets certain criteria.
 *
 * This is /not/ meant to model a chemical reaction, exactly.  A chemical modification doesn't imply that
 * you can biologically start with the input Adduct and end up with the output.  Rather, this represents
 * a way to model compounds that we think are likely to exist, based on the existence of another compound.
 *
 * The modification 'happens' to an Adduct, rather than to a Compound, because an Adduct is a more general case --
 * an Adduct with no modifications is equivalent to its Compound.
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
