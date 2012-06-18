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
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.formula.MassToFormulaTool;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecularFormulaSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.formula.MolecularFormulaRange;
import org.openscience.cdk.formula.rules.*;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;


import java.util.*;
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
    protected static final IChemObjectBuilder CDKObjBuilder = NoNotificationChemObjectBuilder.getInstance();
    protected static IsotopeFactory ifac= null; 
    protected static final MassToFormulaTool mfTool = new MassToFormulaTool(NoNotificationChemObjectBuilder.getInstance());

    //Default constraints on numbers of each atom in formulas.
    //todo: is this redundant with Seven Golden Rules?
    protected static IRule _defaultMFRangeRule = null;
    //this will contain default IRules that are appropriate for the following ranges of masses, in order:
    //0-500
    //500-1000
    //1000-2000
    //2000-3000
    //3000+
    protected static List<IRule>[] _defaultRulesForMasses = null;

    //Initialize CDK isotope factory and default IRules
    static
    {
        try
        {
            ifac = IsotopeFactory.getInstance(CDKObjBuilder);
        }
        catch (Exception e)
        {
            System.err.println("CDK initialization failure: "+e);
        }

        MolecularFormulaRange mfRange = new MolecularFormulaRange();
        mfRange.addIsotope( ifac.getMajorIsotope("C"), 1, 50);
        mfRange.addIsotope( ifac.getMajorIsotope("H"), 1, 100);
        mfRange.addIsotope( ifac.getMajorIsotope("O"), 0, 50);
        mfRange.addIsotope( ifac.getMajorIsotope("N"), 0, 50);
        mfRange.addIsotope( ifac.getMajorIsotope("S"), 0, 5);
        mfRange.addIsotope( ifac.getMajorIsotope("P"), 0, 5);
        mfRange.addIsotope( ifac.getMajorIsotope("Cl"), 0, 2);
        _defaultMFRangeRule  = new ElementRule();
        Object[] params = new Object[1];
        params[0] = mfRange;
        try
        {
            _defaultMFRangeRule.setParameters(params);
        }
        catch (Exception e)
        {
            throw new RuntimeException("Failed to create default MolecularFormulaRange rule",e);
        }

        Object[] mmParams = new Object[]
                {
                        MMElementRule.RangeMass.Minus500,
                        MMElementRule.RangeMass.Minus1000,
                        MMElementRule.RangeMass.Minus2000,
                        MMElementRule.RangeMass.Minus3000
                };
        _defaultRulesForMasses = (List<IRule>[]) new ArrayList[mmParams.length + 1];

        for (int i=0; i<mmParams.length; i++)
        {
            MMElementRule sevenGolden = new MMElementRule();
            Object[] paramsMM = new Object[2];
            paramsMM[0] = MMElementRule.Database.WILEY;
            paramsMM[1] = mmParams[i];
            try
            {
                sevenGolden.setParameters(paramsMM);
            }
            catch (Exception e)
            {
                throw new RuntimeException("Failed to create default SevenGolden rule",e);
            }

            List<IRule> rules = new ArrayList<IRule>();
            rules.add(_defaultMFRangeRule);
            rules.add(sevenGolden);
            _defaultRulesForMasses[i] = rules;
        }
        _defaultRulesForMasses[_defaultRulesForMasses.length-1] = new ArrayList<IRule>();
        _defaultRulesForMasses[_defaultRulesForMasses.length-1].add(_defaultMFRangeRule);
    }

    /**
     * Calculate the "commonest isotope" (i.e., "monoisotope" mass for a chemical formula
     * @param atomCountMap
     * @return
     * @throws IllegalArgumentException
     */
    public static double calcCommonestIsotopeMass(Map<String, Integer> atomCountMap)
            throws IllegalArgumentException
    {
        Double retval = 0D;

        for (String atom : atomCountMap.keySet())
        {
            Element curEl = Elements.get(atom);
            if (curEl == null)
                throw new IllegalArgumentException("Bad formula " + atomCount2FormulaString(atomCountMap) +
                        " contains unknown element " + curEl);
            retval += (curEl.getCommonestIsotopeMass() * atomCountMap.get(atom));
        }

        return retval;
    }

    public static String atomCount2FormulaString(Map<String, Integer> atomCountMap)
    {
        List<String> atoms = new ArrayList<String>(atomCountMap.keySet());
        Collections.sort(atoms);
        StringBuffer resultBuf = new StringBuffer();
        for (String atom : atoms)
            resultBuf.append(atom + atomCountMap.get(atom));
        return resultBuf.toString();
    }

    /**
     * Take a formula and return a map from element symbols to counts of atoms.
     * AtomCounts will have CommonestIsotopeMass
     * @param emp
     * @return
     * @throws IllegalArgumentException
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

    public static ChemicalFormula CDKMolForm2ChemForm(IMolecularFormula f, boolean shouldPopulatePeaks) {
       return new ChemicalFormula(MolecularFormulaManipulator.getString(f),shouldPopulatePeaks);
    }

    public static ChemicalFormula CDKMolForm2ChemForm(IMolecularFormula f) {
        return new ChemicalFormula(MolecularFormulaManipulator.getString(f),false);
    }

    public static IMolecularFormula ChemForm2CDKMolForm(ChemicalFormula f){
        return MolecularFormulaManipulator.getMajorIsotopeMolecularFormula(f.toString(),CDKObjBuilder);
    }



    /**
     * Create an IRule that constrains formula results to a PPM tolerance around a mass
     * @param mass
     * @param ppmerr
     * @return
     */
    public static IRule createIRuleForMassAndTolerance(double mass, double ppmerr)
    {
        ToleranceRangeRule ruleToleran = new ToleranceRangeRule();
        Object[] paramsT = new Object[2];
        double ppm = ppmerr*mass/1.0E6;
        paramsT[0] = mass;
        paramsT[1] = ppm;
        try
        {
            ruleToleran.setParameters(paramsT);
        }
        catch (Exception e)
        {
            throw new RuntimeException(e);
        }
        return ruleToleran;
    }

    /**
     * Get the default IRules for a given mass.  Note: two of these are defaults that are used over and over
     * again, and the third (Seven Golden Rules) is generated uniquely for this call.  So don't do anything
     * to change the first two rules of this result.
     * @param mass
     * @param ppmerr
     * @return
     * @throws Exception
     */
    public static List<IRule> defaultRules(double mass, double ppmerr) throws Exception
    {                        
        List<IRule> basicRules = null;
        int masscode = (int) Math.floor(mass/500);
        switch (masscode) {
            case 0: basicRules = _defaultRulesForMasses[0];
                break;
            case 1: basicRules = _defaultRulesForMasses[1];
                break;
            case 2:
            case 3: basicRules = _defaultRulesForMasses[2];
                break;
            case 4:
            case 5: basicRules = _defaultRulesForMasses[3];
                break;
            default: basicRules = _defaultRulesForMasses[4];
        }

        List<IRule> result = new ArrayList<IRule>(basicRules);

        result.add(createIRuleForMassAndTolerance(mass, ppmerr));

        return result;
    }

    /**
     * Generate a set of formulas for a given mass and set of constraints.  Note: time-consuming
     * @param mass
     * @param constraints
     * @return
     * @throws Exception
     */
    public static IMolecularFormulaSet calcMass2Formulas(double mass, ArrayList<IRule> constraints) throws Exception {
       mfTool.setRestrictions(constraints);
       return mfTool.generate(mass);
    }

    /**
     * Generate a default set of constraints for the specified mass and tolerance, and return all formulas
     * that satisfy those constraints.  Note: time-consuming
     * @param mass
     * @param ppm
     * @return
     * @throws Exception
     */
    public static IMolecularFormulaSet calcMass2Formulas(double mass, double ppm) throws Exception {
       mfTool.setRestrictions(defaultRules(mass,ppm));
       return mfTool.generate(mass);
    }

    public static List<ChemicalFormula> CDKFormulaSet2ChemFormList(IMolecularFormulaSet imfs) {
       ArrayList<ChemicalFormula> retVal = new ArrayList<ChemicalFormula>();
       if(imfs == null) return retVal;
       for(IMolecularFormula imf: imfs.molecularFormulas()) {
           retVal.add(CDKMolForm2ChemForm(imf));
       }
       return retVal; 
    }

    public static void  main(String argv[]) {
        System.out.print(argv[0]);
        Double d = calcCommonestIsotopeMass(chemicalFormula2AtomCount(argv[0]));
        System.out.println("\t"+d);
        HashMap<String, Integer> formula = chemicalFormula2AtomCount(argv[0]);
        System.err.println(formula);
        ChemicalFormula cf = new ChemicalFormula(argv[0]);
        System.out.println("A new chemical formula: "+cf);
        IMolecularFormula imf = ChemForm2CDKMolForm(cf);
        System.out.println("A new ImolForm: "+imf);
        try {
           IMolecularFormulaSet imfs = calcMass2Formulas(MolecularFormulaManipulator.getMajorIsotopeMass(imf),2.0);
           System.out.println("Formulas calculated for natural mass of "+argv[0]+" ("+MolecularFormulaManipulator.getMajorIsotopeMass(imf)+")");
           for(IMolecularFormula i: imfs.molecularFormulas()) {
             System.out.println(MolecularFormulaManipulator.getString(i));
           }
        } catch (Exception e) {
            System.err.println("Formula calculation from mass failed: "+e);
        }
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
                        for (int i=0; i<Math.min(maxPeaksToCalc,elementFrequencies.length); i++)
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
                        newPeakMassesCum[k] = (hpk == 0) ? peakMassesCum[k] : massFormulaNumerator / hpk;
                    }
                    peakProbsCum = newPeakProbsCum;
                    peakMassesCum = newPeakMassesCum;
                }
            }
        }

        int lastNonzeroMassPeak = 0;
        for (int i=0; i<maxPeaksToCalc; i++)
        {
            if (peakMassesCum[i] != 0)
               lastNonzeroMassPeak = i;
        }
        if (lastNonzeroMassPeak < maxPeaksToCalc-1)
        {
            double[] peakProbsCumShorter = new double[lastNonzeroMassPeak+1];
            double[] peakMassesCumShorter = new double[lastNonzeroMassPeak+1];
            System.arraycopy(peakMassesCum, 0, peakMassesCumShorter, 0, lastNonzeroMassPeak+1);
            System.arraycopy(peakProbsCum, 0, peakProbsCumShorter, 0, lastNonzeroMassPeak+1);
            peakMassesCum = peakMassesCumShorter;
            peakProbsCum = peakProbsCumShorter;            
        }
//if (peakMassesCum.length > 1 && (peakMassesCum[1] - peakMassesCum[0] > 1.5))
//{
//    for (String elem : elementAtomCountMap.keySet()) System.err.print(elem);
//    System.err.println("******" + peakMassesCum[0] + ", " + peakMassesCum[1]);
//}
        return new Pair<double[],double[]>(peakMassesCum, peakProbsCum);
    }


    public static String createSMILESString(IMolecule molecule)
    {
        SmilesGenerator sg = new SmilesGenerator();
        return sg.createSMILES(molecule);
    }






}
