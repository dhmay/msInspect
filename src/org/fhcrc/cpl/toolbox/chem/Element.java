package org.fhcrc.cpl.toolbox.chem;
/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 25, 2010
 * Time: 11:56:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class Element {
    String symbol;
    Double canonical_mass;
    int atomic_number;

    public String getSymbol() {
        return symbol;
    }

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    public Double getCanonical_mass() {
        return canonical_mass;
    }

    public void setCanonical_mass(Double canonical_mass) {
        this.canonical_mass = canonical_mass;
    }

    public int getAtomic_number() {
        return atomic_number;
    }

    public void setAtomic_number(int atomic_number) {
        this.atomic_number = atomic_number;
    }


    public Element(String symbol, Double canonical_mass, int atomic_number) {
        this.symbol = symbol;
        this.canonical_mass = canonical_mass;
        this.atomic_number = atomic_number;
    }
}
