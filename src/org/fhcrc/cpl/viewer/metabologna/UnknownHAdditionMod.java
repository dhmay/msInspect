package org.fhcrc.cpl.viewer.metabologna;

import org.fhcrc.cpl.toolbox.chem.ChemicalModification;
import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.chem.ChemicalFormula;

/**
 * Created by IntelliJ IDEA.
 * User: dhmay
 * Date: May 27, 2011
 * Time: 10:15:41 AM
 * Cheat, and add an H that gets reflected in the formula but not in the chemical structure.  Pretend we can
 * always do this.
 */
public class UnknownHAdditionMod implements ChemicalModification {
    public void perform(Adduct adduct) {
        ChemicalFormula newFormula = new ChemicalFormula(adduct.getFormula());
        newFormula.addFormula(new ChemicalFormula("H1"));
        adduct.setFormula(newFormula);
        adduct.getModifications().add(this);
    }

    public boolean canPerform(Adduct adduct) {
        return true;
    }

    public String getSymbol() {
        return "+ H";
    }

    public String getName() {
        return "UnknownHAdditionMod";
    }
}
