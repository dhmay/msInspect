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

/**
 * Represents a chemical modification that can create a new Adduct based on an existing Adduct, if the
 * input Adduct meets certain criteria.
 *
 * This is /not/ meant to model a chemical reaction, exactly.  A chemical modification doesn't imply that
 * you can biologically start with the input Adduct and end up with the output.  Rather, this represents
 * a way to model compounds that we think are likely to exist, based on the existence of another compound.
 *
 * The modification 'happens' to an Adduct, rather than to a Compound, because an Adduct is a more general case --
 * an Adduct with no modifications is equivalent to its Compound.
*/
public interface ChemicalModification
{
    /**
     * Perform the modification on the adduct.  To be absolutely clear, adduct will be altered
     * @param adduct
     * @return
     */
    public void perform(Adduct adduct);

    /**
     * Determines whether the modification can be performed on a particular Adduct
     * @param adduct
     * @return
     */
    public boolean canPerform(Adduct adduct);

    /**
     * Return a symbol to represent this modification, e.g., +2H
     * @return
     */
    public String getSymbol();

    /**
     * A human-readable name
     * @return
     */
    public String getName();

}
