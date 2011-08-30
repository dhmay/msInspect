package org.fhcrc.cpl.viewer.metabologna;

import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.chem.ChemicalFormula;
import org.fhcrc.cpl.toolbox.chem.ChemicalModification;

/**
 * Created by IntelliJ IDEA.
 * User: dhmay
 * Date: 8/30/11
 * Time: 12:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class UnknownFormulaAdditionMod implements ChemicalModification {
    protected String symbol;
    protected ChemicalFormula formula;
    protected String name;

    public UnknownFormulaAdditionMod(ChemicalFormula formula, String symbol, String name) {
        this.symbol = symbol;
        this.formula = formula;
        this.name = name;
    }

    public void perform(Adduct adduct) {
        ChemicalFormula newFormula = new ChemicalFormula(adduct.getFormula());
        newFormula.addFormula(formula);
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