package org.fhcrc.cpl.toolbox.chem;

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import java.util.Map;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.IOException;

/**
 * Represents a named chemical compound.  Chemical formula need not be unique
* User: dhmay
* Date: Apr 6, 2010
* Time: 4:06:07 PM
* To change this template use File | Settings | File Templates.
*/
public class ChemicalCompound
{
    protected double mass;
    protected String name;
    protected String formula;
    protected double[] peakFrequencies;
    protected double[] peakMasses;
    protected Map<String, Integer> elementCountMap;

    public ChemicalCompound(String name, String formula) throws IllegalArgumentException
    {
        this.name = name;
        this.formula = formula;
        elementCountMap = ChemCalcs.chemicalFormula2AtomCount(formula);
        Pair<double[], double[]> massesAndProbabilities = ChemCalcs.calcPeakMassesAndProbabilities(elementCountMap, 2);
        peakMasses = massesAndProbabilities.first;
        peakFrequencies = massesAndProbabilities.second;
        mass = peakMasses[0];
    }

    public double getMass()
    {
        return mass;
    }

    public void setMass(double mass)
    {
        this.mass = mass;
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }

    public String getFormula() {
        return formula;
    }

    public void setFormula(String formula) {
        this.formula = formula;
    }

    public static class ComparatorMassAsc implements Comparator<ChemicalCompound>
    {
        public int compare(ChemicalCompound o1, ChemicalCompound o2)
        {
            if (o1.mass == o2.mass)
                return 0;
            return o1.mass == o2.mass ? 0 : o1.mass < o2.mass ? -1 : 1;
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
        StringBuffer buf = new StringBuffer("Compound: " + name + ", formula=" + formula +
                ", monomass=" + mass);
        if (peakMasses != null)
        {
            buf.append(" all masses (frequencies): ");
            for (int i=0; i< peakMasses.length; i++)
            {
                if (i>0)
                    buf.append(", ");
                buf.append(peakMasses[i] + " (" + peakFrequencies[i] + ")");
            }
        }

        return buf.toString();
    }

    public double[] getPeakFrequencies() {
        return peakFrequencies;
    }

    public double[] getPeakMasses() {
        return peakMasses;
    }

    public Map<String, Integer> getElementCountMap() {
        return elementCountMap;
    }


    /**
     * Load chemical compounds from a tsv file containing columns name and formula (at least)
     * @param file
     * @return
     * @throws IOException
     */
    public static final List<ChemicalCompound> loadCompoundsFromFile(File file) throws IOException
    {
        TabLoader loader = new TabLoader(file);

        Map[] rowsAsMaps = (Map[])loader.load();

        List<ChemicalCompound> result = new ArrayList<ChemicalCompound>();

        for (Map rowMap : rowsAsMaps)
        {
            try
            {
                ChemicalCompound compound = new ChemicalCompound((String) rowMap.get("name"),
                        rowMap.get("formula").toString());
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
