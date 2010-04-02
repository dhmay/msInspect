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
    public static HashMap<String,AtomCount> emp2atomCount(String emp) {
       HashMap<String,AtomCount> retVal = new HashMap<String,AtomCount>();
       Matcher formulaSplit = atomCounter.matcher(emp);
       while(formulaSplit.find()) {
           String curAtomCount = formulaSplit.group();
           if(!Character.isDigit(curAtomCount.charAt(curAtomCount.length()-1))) curAtomCount += "1";
           String atom = curAtomCount.split("[0-9]")[0];
           Element curEl = Elements.get(atom);
           Double curMass = curEl.getCommonestIsotopeMass();
           int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
           retVal.put(atom,new AtomCount(curMass,curRep,curEl));
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
}
