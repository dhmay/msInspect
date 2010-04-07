package org.fhcrc.cpl.toolbox.chem;

import org.fhcrc.cpl.toolbox.datastructure.Pair;

import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 25, 2010
 * Time: 12:39:19 PM
 * To change this template use File | Settings | File Templates.
 */
public final class ChemCalcs {
    protected static final Pattern atomCounter = Pattern.compile("[A-Z][a-z]?[0-9]*");

    /**
     * Calculate the "commonest isotope" (i.e., "monoisotope" mass for a chemical formula
     * @param emp
     * @return
     */
    public static double calcCommonestIsotopeMass(String emp) {
        Double retval = 0D;
        Matcher formulaSplit = atomCounter.matcher(emp);
        while(formulaSplit.find()) {
           String curAtomCount = formulaSplit.group();
           if(!Character.isDigit(curAtomCount.charAt(curAtomCount.length()-1))) curAtomCount += "1";
           String atom = curAtomCount.split("[0-9]")[0];
           Double curMass = Elements.get(atom).getCommonestIsotopeMass();
           int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
           retval += (curMass*curRep);
        }

        return retval;
    }

    /**
     * Take a formula and return a map from element symbols to counts of atoms.
     * AtomCounts will have CommonestIsotopeMass
     * @param emp
     * @return
     */
    public static HashMap<String, Integer> chemicalFormula2AtomCount(String emp)
            throws IllegalArgumentException
    {
        HashMap<String, Integer> retVal = new HashMap<String,Integer>();
        try
        {
            Matcher formulaSplit = atomCounter.matcher(emp);
            while(formulaSplit.find())
            {
                String curAtomCount = formulaSplit.group();
                if(!Character.isDigit(curAtomCount.charAt(curAtomCount.length()-1))) curAtomCount += "1";
                String atom = curAtomCount.split("[0-9]")[0];
                Element curEl = Elements.get(atom);
                if (curEl == null)
                    throw new IllegalArgumentException("Bad formula " + emp + " contains unknown element " + curEl);
                int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
                retVal.put(atom, curRep);
            }
        }
        catch (IllegalArgumentException iae)
        {
            throw iae;
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException ("Bad formula " + emp);
        }
        return retVal;
    }

    public static void  main(String argv[]) {
        System.out.print(argv[0]);
        Double d = calcCommonestIsotopeMass(argv[0]);
        System.out.println("\t"+d);
        HashMap<String, Integer> formula = chemicalFormula2AtomCount(argv[0]);
        System.err.println(formula);
    }

    public static double[] calcPeakProbabilities(String formula, int maxPeaksToCalc)
            throws IllegalArgumentException
    {
        return calcPeakMassesAndProbabilities(formula, maxPeaksToCalc).second;
    }

    public static double[] calcPeakMasses(String formula, int maxPeaksToCalc)
            throws IllegalArgumentException
    {
        return calcPeakMassesAndProbabilities(formula, maxPeaksToCalc).first;
    }

    /**
     * This method calculates the relative abundance of isotopic peaks within the specified
     * chemical formula.  It uses the algorithm described here:
     *
     * Efficient Calculation of Accurate Masses of Isotopic Peaks
     * Alan L. Rockwood, Perttu Haimi
     * Journal of the American Society for Mass Spectrometry
     * Volume 17, Issue 3, March 2006, Pages 415-419
     *
     * Specifically, I begin with one single atom from the chemical formula and then
     * add each remaining atom, one at a time.  Adding each atom involves creating a
     * 'super-atom' from the cumulative peak frequencies and the new atom's peak frequencies.
     *
     * Throws an IllegalArgumentException if the formula isn't valid
     * @param formula
     * @param maxPeaksToCalc
     * @return
     */
    public static Pair<double[],double[]> calcPeakMassesAndProbabilities(String formula, int maxPeaksToCalc)
            throws IllegalArgumentException
    {
        HashMap<String, Integer> elementAtomCountMap = chemicalFormula2AtomCount(formula);
        return calcPeakMassesAndProbabilities(elementAtomCountMap, maxPeaksToCalc);
    }

    public static Pair<double[],double[]> calcPeakMassesAndProbabilities(Map<String, Integer> elementAtomCountMap,
                                                                         int maxPeaksToCalc)
    {
        //parse the formula into a map of atoms to counts


        //This array holds our peak probabilities.  It changes with each atom added
        double[] peakMassesCum = new double[maxPeaksToCalc];
        double[] peakProbsCum = new double[maxPeaksToCalc];

        //special handling for first atom considered
        boolean firstAtom = true;
        for (String element : elementAtomCountMap.keySet())
        {
            double[] elementFrequencies = Elements.get(element).getIsotopicPeakFrequenciesWithMissing();
            double[] elementMasses = Elements.get(element).getIsotopicPeakMassesWithMissing();


            //how many of this element do we have?
            int elementCount = elementAtomCountMap.get(element);

            //one by one, add the atoms of this element
            for (int elemCount=0; elemCount < elementCount; elemCount++)
            {
                //If first atom, start off with the frequencies for that atom
                if (firstAtom)
                {
                    //simply set the first elements of peakProbsCum to that atom's freqs
                    System.arraycopy(elementFrequencies, 0, peakProbsCum, 0,
                            Math.min(elementFrequencies.length, maxPeaksToCalc));
                    System.arraycopy(elementMasses, 0, peakMassesCum, 0,
                            Math.min(elementMasses.length, maxPeaksToCalc));
                    firstAtom = false;
                }
                else
                {
                    //new peak probabilities array values will replace old ones
                    double[] newPeakProbsCum = new double[maxPeaksToCalc];
                    double[] newPeakMassesCum = new double[maxPeaksToCalc];

                    //for every atom after the first, "combine" it with the accumulated 'super-atom'.
                    //Define gp(i) as the new atom's probability for peak i.
                    //Define fp(i) as the super-atom's probability for peak i.
                    //The formula for each peak's probability hp(k) is:
                    //hp(k) = sum over all i (gp(i) * fp(k-i))
                    //
                    //masses:
                    //Define gm(i) as the new atom's mass for peak i.
                    //Define fm(i) as the super-atom's mass for peak i.
                    //formula for each peak's mass hm(k) is:
                    // (sum over all i (gp(i) * fp(k-i) * (gm(i) + fm(k-i)) ) ) / hp(k)
                    for (int k=0; k<maxPeaksToCalc; k++)
                    {
                        //begin with 0, sum up components
                        double hpk = 0;
                        double massFormulaNumerator = 0;
                        for (int i=0; i<Math.min(maxPeaksToCalc, elementFrequencies.length); i++)
                        {
                            int kminusi = k-i;
                            //check indexes both valid
                            if (kminusi < 0 || kminusi >= maxPeaksToCalc)
                                continue;
                            double gpi = elementFrequencies[i];
                            double fpkminusi = peakProbsCum[kminusi];

                            hpk += gpi * fpkminusi;

                            massFormulaNumerator +=
                                    gpi * fpkminusi * (elementMasses[i] + peakMassesCum[kminusi]);
                        }
                        newPeakProbsCum[k] = hpk;
                        newPeakMassesCum[k] = massFormulaNumerator / hpk;
                    }
                    peakProbsCum = newPeakProbsCum;
                    peakMassesCum = newPeakMassesCum;
                }
            }
        }

        return new Pair<double[],double[]>(peakMassesCum, peakProbsCum);
    }









}
