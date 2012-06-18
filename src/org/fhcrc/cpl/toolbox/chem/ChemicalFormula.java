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
    //mass of the "monoisotope", or technically "commonest isotope". This must be updated when any change happens
    protected double commonestIsotopeMass;


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
     * Create a ChemicalFormula that's the same as another ChemicalFormula, but without peaks populated
     * @param otherFormula
     */
    public ChemicalFormula(ChemicalFormula otherFormula)
    {
        this(otherFormula.getElementCountMap(), 0);
    }

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
        this.elementCountMap = new HashMap<String, Integer>(elementCountMap);
        commonestIsotopeMass = ChemCalcs.calcCommonestIsotopeMass(elementCountMap);
        if (numPeaksToPopulate > 0)
            populatePeakMassesAndFrequencies(numPeaksToPopulate);
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
        this(ChemCalcs.chemicalFormula2AtomCount(formulaString), numPeaksToPopulate);
    }

    /**
     * Add additional elements to this formula.  If peak masses, freqs already populated, update them
     * @param additionFormula
     */
    public void addFormula(ChemicalFormula additionFormula)
    {
        for (String atom : additionFormula.getElementCountMap().keySet())
            elementCountMap.put(atom, additionFormula.getElementCountMap().get(atom) +
                    (elementCountMap.containsKey(atom) ? elementCountMap.get(atom) : 0));
        commonestIsotopeMass = ChemCalcs.calcCommonestIsotopeMass(elementCountMap);
        if (peakFrequencies != null)
            populatePeakMassesAndFrequencies(peakFrequencies.length);
    }

    /**
     * Create a new ChemicalFormula identical to this one with the additional elements added.  Do not
     * populate peaks
     * @param additionFormula
     * @return
     */
    public ChemicalFormula createFormulaWithAddition(ChemicalFormula additionFormula)
    {
        ChemicalFormula newFormula = new ChemicalFormula(this);
        newFormula.addFormula(additionFormula);
        return newFormula;
    }

    /**
     * Remove specified elements.  Populate peaks if they were already populated.  Throw IllegalArgumentException
     * if the formula doesn't have the specified elements
     * @param subtractionFormula
     * @return
     */
    public ChemicalFormula createFormulaWithSubtraction(ChemicalFormula subtractionFormula)
            throws IllegalArgumentException
    {
        ChemicalFormula newFormula = new ChemicalFormula(this);
        newFormula.subtractFormula(subtractionFormula);
        return newFormula;
    }

    /**
     * Create a formula identical to this one with elements removed.  Populate peaks if they were already populated.
     * Throw IllegalArgumentException if the formula doesn't have the specified elements
     * @param subtractionFormula
     * @return
     */
    public void subtractFormula(ChemicalFormula subtractionFormula)
            throws IllegalArgumentException
    {
        for (String atom : subtractionFormula.getElementCountMap().keySet())
        {
            if (!elementCountMap.containsKey(atom))
                throw new IllegalArgumentException("Can't remove nonpresent element " + atom);
            int numPresent = elementCountMap.get(atom);
            int numForNewMap = numPresent - subtractionFormula.getElementCountMap().get(atom);
            if (numForNewMap < 0)
                throw new IllegalArgumentException("Can't remove " + subtractionFormula.getElementCountMap().get(atom) +
                        " of element " + atom + ", only " + numPresent + " present");
            if (numForNewMap == 0)
                elementCountMap.remove(atom);
            else
                elementCountMap.put(atom, numForNewMap);
        }
        commonestIsotopeMass = ChemCalcs.calcCommonestIsotopeMass(elementCountMap);
        if (peakFrequencies != null)
            populatePeakMassesAndFrequencies(peakFrequencies.length);
    }


    public double getCommonestIsotopeMass()
    {
        return commonestIsotopeMass;
    }

    public void setCommonestIsotopeMass(double mass)
    {
        this.commonestIsotopeMass = mass;
    }

    public String toString()
    {
        return ChemCalcs.atomCount2FormulaString(getElementCountMap());
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

    public String getPeakMassesFrequenciesString()
    {
        StringBuffer buf = new StringBuffer("");
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

    /**
     * Equality check.  This is for storing things in hashtables by ChemicalFormula.
     * NOTE! ChemicalFormulas with identical formulas will hash to the same thing
     * @param otherFormula
     * @return
     */
    public boolean equals(Object otherFormula)
    {
         return otherFormula instanceof ChemicalFormula && toString().equals(otherFormula.toString());
    }

    /**
     * Hash code generation.  WARNING! This will generate the same code as a String with the formula contents.
     * So don't put ChemicalFormulas and Strings in the same HashMap.
     * @return
     */
    public int hashCode()
    {
        return toString().hashCode();
    }

    /**
     * Calculate the nominal mass, defined as the number of nucleons
     * todo: is this the standard way to calculate it?
     * @return
     */
    public int getNominalMass()
    {
        int nominalMass = 0;
        for (String symbol : elementCountMap.keySet())
        {
            nominalMass += elementCountMap.get(symbol) * Elements.get(symbol).getNominalMass();
//if (toString().equals("C5H7N3O1")) System.err.println(symbol + ", " + elementCountMap.get(symbol) + ", " + Elements.get(symbol).getNominalMass());            
        }
        return nominalMass;
    }

    /**
     * Calculate the mass defect of the entire molecule.  This is NOT the unit mass defect
     * @return
     */
    public double calcTotalMassDefect()
    {
//System.err.println(toString() + ", " + getNominalMass() + ", " + commonestIsotopeMass);        
        //todo: is this standard?
        return commonestIsotopeMass - getNominalMass();
    }

    /**
     * Calculate the unit mass defect.
     * todo: is this standard?
     * @return
     */
    public double calcUnitMassDefect()
    {
        return calcTotalMassDefect() / getNominalMass();
    }
}
