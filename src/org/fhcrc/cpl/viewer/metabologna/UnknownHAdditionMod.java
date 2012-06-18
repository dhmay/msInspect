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
public class UnknownHAdditionMod extends UnknownFormulaAdditionMod {
    public UnknownHAdditionMod() {
        super(new ChemicalFormula("H1"),"M+H","UnknownHAdditionMod");
    }

}
