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