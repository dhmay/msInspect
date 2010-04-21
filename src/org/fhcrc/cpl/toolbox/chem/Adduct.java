package org.fhcrc.cpl.toolbox.chem;


import java.util.*;

/**
 * Represents a chemical compound with optional added AND/OR subtracted elements.
 * todo: revisit terminology?
* User: dhmay
* Date: Apr 6, 2010
* Time: 4:06:07 PM
* To change this template use File | Settings | File Templates.
*/
public class Adduct
{
    //None of these variables should ever be set without using the setters.  setting any of these with the
    //setters will update formula

    protected ChemicalCompound compound;
    protected ChemicalFormula formula;

    //These are guaranteed never to be null unless explicitly set that way by pathological code
    protected List<ChemicalModification> modifications = new ArrayList<ChemicalModification>();

    /**
     * Create a fully-defined Adduct.  If nulls passed for added or subtracted formulas, instantiate empty lists.
     * Throws IllegalArgumentException if all subtractedFormulas can't be subtracted
     * @param compound
     * @param modifications
     */
    public Adduct(ChemicalCompound compound,
                  List<ChemicalModification> modifications)
            throws IllegalArgumentException
    {
        setCompound(compound);

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
     * Use with EXTREME CAUTION.  This is used by, e.g., ChemicalModification
     * @param formula
     */
    public void setFormula(ChemicalFormula formula) {
        this.formula = formula;
    }

    public String getCompoundNameAndIonTypeString()
    {
        return compound.getName() + ":" + getIonTypeString();
    }

    public String getIonTypeString()
    {
        StringBuffer resultBuf = new StringBuffer("[M");
        for (ChemicalModification mod : modifications)
        {
            resultBuf.append(" ");
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


}
