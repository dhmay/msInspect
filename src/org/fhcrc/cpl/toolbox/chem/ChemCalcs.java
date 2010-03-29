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

    public static Double emp2mass(String emp) {
        return emp2mass(emp, Elements.CommonestIsotopes);
    }

    public static Double emp2mass(String emp,HashMap<String, Element> massMap) {
        Double retval = 0D;
        Matcher formulaSplit = atomCounter.matcher(emp);
        while(formulaSplit.find()) {
           String curAtomCount = formulaSplit.group();
           if(!Character.isDigit(curAtomCount.charAt(curAtomCount.length()-1))) curAtomCount += "1";
           String atom = curAtomCount.split("[0-9]")[0];
           Double curMass = massMap.get(atom).getCanonical_mass();
           int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
           retval += (curMass*curRep);
        }

        return retval;
    }

    public static HashMap<String,atomCount> emp2atomCount(String emp) {
        return emp2atomCount(emp,Elements.CommonestIsotopes);
    }

    public static HashMap<String,atomCount> emp2atomCount(String emp, HashMap<String,Element>massMap) {
       HashMap<String,atomCount> retVal = new HashMap<String,atomCount>();
       Matcher formulaSplit = atomCounter.matcher(emp);
       while(formulaSplit.find()) {
           String curAtomCount = formulaSplit.group();
           if(!Character.isDigit(curAtomCount.charAt(curAtomCount.length()-1))) curAtomCount += "1";
           String atom = curAtomCount.split("[0-9]")[0];
           Element curEl = massMap.get(atom);
           Double curMass = curEl.getCanonical_mass();
           int curRep = Integer.parseInt(curAtomCount.split("[A-Za-z]+")[1]);
           retVal.put(atom,new atomCount(curMass,curRep,curEl));
       }    
       return retVal;
    }

    public static void  main(String argv[]) {
        System.out.print(argv[0]);
        Double d = calcs.emp2mass(argv[0]);
        System.out.println("\t"+d);
        HashMap<String,atomCount> formula = emp2atomCount(argv[0]);
        System.err.println(formula);
    }
}
