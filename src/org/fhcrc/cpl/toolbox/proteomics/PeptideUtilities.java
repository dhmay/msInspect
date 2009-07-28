/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.toolbox.proteomics;

import org.apache.log4j.Logger;

import java.util.List;

/**
 * Utilities for working with peptide sequences and modifications
 */
public class PeptideUtilities
{
    protected static Logger _log = Logger.getLogger(PeptideUtilities.class);

    public static final float MOD_MASS_SLOP = 0.1f;


    /**
     * Given a peptide sequence and a list of modifications for every position in the sequence,
     * and a specified modification residue and mass, checks to see whether the peptide contains the
     * modification.
     *
     * Returns true if the peptide does NOT contain the modification on any residue, false if it
     * does.  Special case: returns false if the peptide does not contain the residue.
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    public static boolean checkForModNoResidues(String peptide, List<ModifiedAminoAcid>[] mods, 
                                            char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return true;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    continue;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                        return false;
            }
        }
        return true;
    }

    /**
     * Given a peptide sequence and a list of modifications for every position in the sequence,
     * and a specified modification residue and mass, checks to see whether the peptide contains the
     * modification.
     *
     * Returns true if the peptide contains the modification on every single occurrence of the residue,
     * false if it does not.  Special case: returns false if the peptide does not contain the residue.
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    public static boolean checkForModAllResidues(String peptide, List<ModifiedAminoAcid>[] mods,
                                                 char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return false;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                boolean foundIt = false;
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    return false;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                    {
                        foundIt = true;
                        break;
                    }
                if (!foundIt)
                    return false;
            }
        }
        return true;
    }
}
