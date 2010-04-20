package org.fhcrc.cpl.toolbox.chem;

import org.fhcrc.cpl.toolbox.datastructure.Pair;

import java.util.*;

/**
* Represents a chemical formula, with useful helper methods.
* User: dhmay
* Date: Apr 6, 2010
* Time: 4:06:07 PM
* To change this template use File | Settings | File Templates.
*/
public class ChemicalFormula
{
    //mass of the "monoisotope", or technically "commonest isotope"
    protected double commonestIsotopeMass;

    //formula as a String.  Preserve this in case the caller wants the order of elements preserved
    protected String formulaString;

    //this is the most important thing.  How much of what elements
    protected Map<String, Integer> elementCountMap;

    //peak frequencies and masses are not necessarily populated, because that's somewhat costly
    protected double[] peakFrequencies;
    protected double[] peakMasses;

    //default number of peaks whose freqs and masses we should calculate, /if/ we calculate
    public static final int DEFAULT_NUMPEAKS_TOCALC = 3;
    //controls the number of peaks calculated if number not supplied
    public static int numPeaksToCalculate = DEFAULT_NUMPEAKS_TOCALC;

    /**
     * Creates a ChemicalFormula and does not populate peak masses and frequencies
     * @param formulaString
     * @throws IllegalArgumentException
     */
    public ChemicalFormula(String formulaString) throws IllegalArgumentException
    {
        this(formulaString, false);
    }

    /**
     * Creates a ChemicalFormula, only populating peak masses and frequencies if specified
     * @param elementCountMap
     * @param numPeaksToPopulate
     * @throws IllegalArgumentException
     */
    public ChemicalFormula(Map<String, Integer> elementCountMap, int numPeaksToPopulate)
            throws IllegalArgumentException
    {
        this.elementCountMap = elementCountMap;
        commonestIsotopeMass = ChemCalcs.calcCommonestIsotopeMass(elementCountMap);
    }

    /**
     * Creates a ChemicalFormula, only populating peak masses and frequencies if specified
     * @param formulaString
     * @param shouldPopulatePeaks
     * @throws IllegalArgumentException
     */
    public ChemicalFormula(String formulaString, boolean shouldPopulatePeaks) throws IllegalArgumentException
    {
        this(ChemCalcs.chemicalFormula2AtomCount(formulaString),
                shouldPopulatePeaks ? numPeaksToCalculate : 0);

    }

    /**
     * Creates a ChemicalFormula and populates peak masses and frequencies out to the specified # of peaks
     * @param formulaString
     * @param numPeaksToPopulate
     * @throws IllegalArgumentException
     */
    public ChemicalFormula(String formulaString, int numPeaksToPopulate) throws IllegalArgumentException
    {
        this.formulaString = formulaString;
        elementCountMap = ChemCalcs.chemicalFormula2AtomCount(formulaString);
        if (numPeaksToPopulate > 0)
            populatePeakMassesAndFrequencies(numPeaksToPopulate);
    }

    /**
     * Create a new ChemicalFormula identical to this one with the additional elements added.  Do not
     * populate peaks
     * @param additionElementCountMap
     * @return
     */
    public ChemicalFormula createFormulaWithAddition(Map<String, Integer> additionElementCountMap)
    {
        HashMap<String, Integer> newElementCountMap = new HashMap<String, Integer>(elementCountMap);
        for (String atom : additionElementCountMap.keySet())
            newElementCountMap.put(atom, additionElementCountMap.get(atom) +
                    (newElementCountMap.containsKey(atom) ? newElementCountMap.get(atom) : 0));
        return new ChemicalFormula(newElementCountMap, 0);
    }

    /**
     * Create a new ChemicalFormula identical to this one with the specified elements removed.  Do not
     * populate peaks.  Throw IllegalArgumentException if the formula doesn't have the specified elements
     * @param subtractionElementCountMap
     * @return
     */
    public ChemicalFormula createFormulaWithSubtraction(Map<String, Integer> subtractionElementCountMap)
            throws IllegalArgumentException
    {
        HashMap<String, Integer> newElementCountMap = new HashMap<String, Integer>(elementCountMap);
        for (String atom : subtractionElementCountMap.keySet())
        {
            if (!elementCountMap.containsKey(atom))
                throw new IllegalArgumentException("Can't remove nonpresent element " + atom);
            int numPresent = elementCountMap.get(atom);
            int numForNewMap = numPresent - subtractionElementCountMap.get(atom);
            if (numForNewMap < 0)
                throw new IllegalArgumentException("Can't remove " + subtractionElementCountMap.get(atom) +
                        " of element " + atom + ", only " + numPresent + " present");
            if (numForNewMap == 0)
                newElementCountMap.remove(atom);
            else
                newElementCountMap.put(atom, numForNewMap);                     
        }
        return new ChemicalFormula(newElementCountMap, 0);
    }


    public double getCommonestIsotopeMass()
    {
        return commonestIsotopeMass;
    }

    public void setCommonestIsotopeMass(double mass)
    {
        this.commonestIsotopeMass = mass;
    }

    public String getFormula() {
        return formulaString;
    }

    public void setFormula(String formula) {
        this.formulaString = formula;
    }

    public static class ComparatorMassAsc implements Comparator<ChemicalFormula>
    {
        public int compare(ChemicalFormula o1, ChemicalFormula o2)
        {
            if (o1.commonestIsotopeMass == o2.commonestIsotopeMass)
                return 0;
            return o1.commonestIsotopeMass == o2.commonestIsotopeMass ? 0 :
                    o1.commonestIsotopeMass < o2.commonestIsotopeMass ? -1 : 1;
        }
    }

    public String toString()
    {
        StringBuffer buf = new StringBuffer("formulaString: " + formulaString +
                ", monomass=" + commonestIsotopeMass);
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

    public double[] getPeakFrequencies()
    {
         return getPeakFrequencies(numPeaksToCalculate);
    }

    public double[] getPeakMasses()
    {
         return getPeakMasses(numPeaksToCalculate);
    }

    /**
     * Can be slightly costly if masses or frequencies not already calculated
     * @return
     */
    public double[] getPeakFrequencies(int numPeaksToCalculate)
    {
        if (peakFrequencies == null || peakFrequencies.length < numPeaksToCalculate)
            populatePeakMassesAndFrequencies(numPeaksToCalculate);
        return peakFrequencies;
    }

    /**
     * Can be slightly costly if masses or frequencies not already calculated
     * @return
     */
    public double[] getPeakMasses(int numPeaksToCalculate)
    {
        if (peakMasses == null || peakMasses.length < numPeaksToCalculate)
            populatePeakMassesAndFrequencies(numPeaksToCalculate);
        return peakMasses;
    }

    public void populatePeakMassesAndFrequencies()
    {
        populatePeakMassesAndFrequencies(numPeaksToCalculate);
    }

    public void populatePeakMassesAndFrequencies(int numPeaksToCalc)
    {
        Pair<double[], double[]> massesAndFrequencies =
                ChemCalcs.calcPeakMassesAndProbabilities(elementCountMap, numPeaksToCalc);
        peakMasses = massesAndFrequencies.first;
        peakFrequencies = massesAndFrequencies.second;
    }

    public Map<String, Integer> getElementCountMap()
    {
        return elementCountMap;
    }

    public static int getNumPeaksToCalculate()
    {
        return numPeaksToCalculate;
    }

    public static void setNumPeaksToCalculate(int numPeaksToCalculate)
    {
        ChemicalFormula.numPeaksToCalculate = numPeaksToCalculate;
    }
}
