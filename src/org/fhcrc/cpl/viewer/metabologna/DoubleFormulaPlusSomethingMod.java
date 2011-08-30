package org.fhcrc.cpl.viewer.metabologna;

import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.chem.ChemicalFormula;
import org.fhcrc.cpl.toolbox.chem.ChemicalModification;

/**
 * Double the original formula and then add addFormula
 */
public class DoubleFormulaPlusSomethingMod implements ChemicalModification {
    protected String symbol;
    protected ChemicalFormula addFormula;
    protected String name;

    public DoubleFormulaPlusSomethingMod(ChemicalFormula addFormula, String symbol, String name) {
        this.symbol = symbol;
        this.addFormula = addFormula;
        this.name = name;
    }

    public void perform(Adduct adduct) {
        ChemicalFormula newFormula = new ChemicalFormula(adduct.getFormula());
        newFormula.addFormula(adduct.getFormula());
        newFormula.addFormula(addFormula);
        adduct.setFormula(newFormula);
        adduct.getModifications().add(this);
    }

    public boolean canPerform(Adduct adduct) {
        return true;
    }

    public String getSymbol() {
        return symbol;
    }

    public String getName() {
        return name;
    }
}
