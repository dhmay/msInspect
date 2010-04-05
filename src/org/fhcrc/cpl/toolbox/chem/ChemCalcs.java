package org.fhcrc.cpl.toolbox.chem;

import java.util.HashMap;
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
     * Take a formula and return a map from element symbols to AtomCounts.
     * AtomCounts will have CommonestIsotopeMass
     * @param emp
     * @return
     */
    public static HashMap<String,AtomCount> emp2atomCount(String emp)
            throws IllegalArgumentException
    {
        HashMap<String,AtomCount> retVal = new HashMap<String,AtomCount>();
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
                Double curMass = curEl.getCommonestIsotopeMass();
                int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
                retVal.put(atom,new AtomCount(curMass,curRep,curEl));
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
        HashMap<String,AtomCount> formula = emp2atomCount(argv[0]);
        System.err.println(formula);
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
    public static double[] calcPeakProbabilities(String formula, int maxPeaksToCalc)
            throws IllegalArgumentException
    {
        //parse the formula into a map of atoms to counts
        HashMap<String, AtomCount> elementAtomCountMap = emp2atomCount(formula);

        //This array holds our peak probabilities.  It changes with each atom added
        double[] peakProbsCum = new double[maxPeaksToCalc];

        //special handling for first atom considered
        boolean firstAtom = true;
        for (String element : elementAtomCountMap.keySet())
        {
            double[] elementFrequencies = Elements.get(element).getIsotopicPeakFrequenciesWithMissing();

            //how many of this element do we have?
            int elementCount = elementAtomCountMap.get(element).getCount();

            //one by one, add the atoms of this element
            for (int elemCount=0; elemCount < elementCount; elemCount++)
            {
                //If first atom, start off with the frequencies for that atom
                if (firstAtom)
                {
                    //simply set the first elements of peakProbsCum to that atom's freqs
                    System.arraycopy(elementFrequencies, 0, peakProbsCum, 0,
                            Math.min(elementFrequencies.length, maxPeaksToCalc));
                    firstAtom = false;
                }
                else
                {
                    //new peak probabilities array values will replace old ones
                    double[] newPeakProbsCum = new double[maxPeaksToCalc];

                    //for every atom after the first, "combine" it with the accumulated 'super-atom'.
                    //Define gp(i) as the new atom's probability for peak i.
                    //Define fp(i) as the super-atom's probability for peak i.
                    //The formula for each peak's probability hp(k) is:
                    //hp(k) = sum over all i (gp(i) * fp(k-i))
                    for (int k=0; k<maxPeaksToCalc; k++)
                    {
                        //begin with 0, sum up components
                        double newProbThisPeak = 0;
                        for (int i=0; i<Math.min(maxPeaksToCalc, elementFrequencies.length); i++)
                        {
                            int kminusi = k-i;
                            //check indexes both valid
                            if (kminusi < 0 || kminusi >= maxPeaksToCalc)
                                continue;
                            double gpi = elementFrequencies[i];
                            double fpkminusi = peakProbsCum[kminusi];

                            newProbThisPeak += gpi * fpkminusi;
                        }
                        newPeakProbsCum[k] = newProbThisPeak;
                    }
                    peakProbsCum = newPeakProbsCum;
                }
            }
        }

        return peakProbsCum;
    }









}
