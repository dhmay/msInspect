package org.fhcrc.cpl.toolbox.chem;

import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;

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
    protected String name;
    protected ChemicalFormula formula;    

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
        return loadCompoundsFromFile(file, numPeaksToPopulate, "name","formula");
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
                                                                     String nameColName, String formulaColName)
            throws IOException
    {
        TabLoader loader = new TabLoader(file);

        Map[] rowsAsMaps = (Map[])loader.load();

        List<ChemicalCompound> result = new ArrayList<ChemicalCompound>();


        for (Map rowMap : rowsAsMaps)
        {
            try
            {
                ChemicalCompound compound = new ChemicalCompound((String) rowMap.get(nameColName),
                        rowMap.get(formulaColName).toString(), numPeaksToPopulate);
                //experimental
                if (rowMap.containsKey("class"))
                {
                    compound.setCompoundClass(rowMap.get("class").toString());
                }
                result.add(compound);
            }
            catch (IllegalArgumentException e)
            {
                ApplicationContext.setMessage("Skipping bad compound: " + e.getMessage());
            }

        }
        return result;
    }


}
