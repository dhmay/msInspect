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
package org.fhcrc.cpl.toolbox.proteomics;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.datastructure.Pair;

import java.util.*;

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

    /**
     * Calculate peptide masses for all possible modification states for this peptide, given the list of
     * possible static and variable modifications.  All variable modifications are either applied to
     * all of the appropriate residue, or none.  So, e.g., ELVISMMM either has all oxidized M's or none.
     *
     * This is done recursively.  This would be somewhat more
     * efficient if it reordered the list to put the static modifications first.
     * @param peptideSequence
     * @param mods
     * @return
     */
    public static List<Pair<List<ModifiedAminoAcid>[], Float>> calculatePeptideMassesForMods(
            String peptideSequence, List<MS2Modification> mods)
    {
        List<MS2Modification> modsThisPeptideHas = new ArrayList<MS2Modification>();
        for (MS2Modification mod : mods)
            if (peptideSequence.contains(mod.getAminoAcid()))
            {
                modsThisPeptideHas.add(mod);
            }

        if (modsThisPeptideHas.isEmpty())
        {
            float mass = (float) new PeptideGenerator().createPeptideForFullyTrypticPeptideSequence(
                    peptideSequence).getMonoisotopicMass();
            Pair<List<ModifiedAminoAcid>[], Float> noModsPair =
                    new Pair<List<ModifiedAminoAcid>[], Float>(new List[peptideSequence.length()], mass);
            List<Pair<List<ModifiedAminoAcid>[], Float>> result = new ArrayList<Pair<List<ModifiedAminoAcid>[], Float>>(1);
            result.add(noModsPair);
            return result;
        }
        else
            return recursivelyAddMassesForMods(peptideSequence, new ArrayList<MS2Modification>(),
                                                          modsThisPeptideHas);
    }

    /**
     * Recursively generate modified sequences and masses for all possible masses, given the list of
     * variable modifications known to exist in this peptide and the list of modifications
     * already applied
     * @param peptideSequence
     * @param appliedMods
     * @param modsRemaining
     * @return
     */
    protected static List<Pair<List<ModifiedAminoAcid>[], Float>> recursivelyAddMassesForMods(
                                                          String peptideSequence,
                                                          List<MS2Modification> appliedMods,
                                                          List<MS2Modification> modsRemaining)
    {
        List<Pair<List<ModifiedAminoAcid>[], Float>> result = new ArrayList<Pair<List<ModifiedAminoAcid>[], Float>>();
        if (modsRemaining.size() == 0)
        {
            Pair<List<ModifiedAminoAcid>[], Float> modsAndMass =
                    calcMassForPeptideWithMods(peptideSequence, appliedMods);
            result.add(modsAndMass);
            return result;
        }

        //This is hacky and wasteful.
        //TODO: find a better way to deal with keeping the integrity of these lists through recursion
        List<MS2Modification> modsRemainingCopy =
                new ArrayList<MS2Modification>(modsRemaining.size());
        modsRemainingCopy.addAll(modsRemaining);
        modsRemaining = modsRemainingCopy;

        List<MS2Modification> appliedModsCopy =
                new ArrayList<MS2Modification>(appliedMods.size());
        appliedModsCopy.addAll(appliedMods);
        appliedMods = appliedModsCopy;


        MS2Modification mod = modsRemaining.get(0);
        modsRemainingCopy.remove(mod);

        if (mod.getVariable())
            result.addAll(recursivelyAddMassesForMods(peptideSequence,
                                                    appliedMods, modsRemaining));
        appliedMods.add(mod);
        result.addAll(recursivelyAddMassesForMods(peptideSequence,
                                                    appliedMods,
                                                    modsRemaining));
        return result;
    }

    /**
     * Modifications are applied all-or-nothing.  If we modify one residue, we modify all
     * of that residue
     * @param peptideSequence
     * @param modifications
     * @return
     */
    protected static Pair<List<ModifiedAminoAcid>[], Float> calcMassForPeptideWithMods(String peptideSequence,
                                                    List<MS2Modification> modifications)
    {
        double mass = PeptideGenerator.computeMass(peptideSequence.getBytes(),
                                                   0, peptideSequence.length(),
                                                   PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
        List<ModifiedAminoAcid>[] modifiedAminoAcids = new ArrayList[peptideSequence.length()];
        if (modifications.size() > 0)
        {
            Map<String, List<MS2Modification>> acidModListMap =
                    new HashMap<String, List<MS2Modification>>();

            for (MS2Modification mod : modifications)
            {
                String acidString = mod.getAminoAcid();
                List<MS2Modification> modListThisAcid = acidModListMap.get(acidString);
                if (modListThisAcid == null)
                {
                    modListThisAcid = new ArrayList<MS2Modification>();
                    acidModListMap.put(acidString, modListThisAcid);
                }
                modListThisAcid.add(mod);
            }


            for (int i=0; i<peptideSequence.length(); i++)
            {
                char thisResidue = peptideSequence.charAt(i);
                double massThisResidue =
                        PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[thisResidue];

                List<MS2Modification> modListThisAcid = acidModListMap.get("" + thisResidue);
                if (modListThisAcid != null)
                {
                    for (MS2Modification mod : modListThisAcid)
                    {
                        mass += mod.getMassDiff();
                        List<ModifiedAminoAcid> modsThisResidue = modifiedAminoAcids[i];
                        if (modsThisResidue == null)
                        {
                            modsThisResidue = new ArrayList<ModifiedAminoAcid>();
                            modifiedAminoAcids[i] = modsThisResidue;
                        }
                        massThisResidue += mod.getMassDiff();
                        ModifiedAminoAcid moddedAcid = new ModifiedAminoAcid(thisResidue,
                                massThisResidue);
                        modsThisResidue.add(moddedAcid);
                    }
                }
            }
        }

        return new Pair<List<ModifiedAminoAcid>[], Float>(modifiedAminoAcids, (float) mass);
    }

    /*
        protected static Pair<List<ModifiedAminoAcid>[], Float> calcMassForPeptideWithMods(String peptideSequence,
                                                    List<MS2Modification> modifications)
    {
        double mass = PeptideGenerator.computeMass(peptideSequence.getBytes(),
                                                   0, peptideSequence.length(),
                                                   PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
        List<ModifiedAminoAcid>[] modifiedAminoAcids = new ArrayList[peptideSequence.length()];
        if (modifications.size() > 0)
        {
            Map<String, List<MS2Modification>> acidStaticModListMap =
                    new HashMap<String, List<MS2Modification>>();
            Map<String, List<MS2Modification>> acidVarModListMap =
                    new HashMap<String, List<MS2Modification>>();

            for (MS2Modification mod : modifications)
            {
                Map<String, List<MS2Modification>> appropriateMap = acidStaticModListMap;
                if (mod.getVariable())
                    appropriateMap = acidVarModListMap;
                String acidString = mod.getAminoAcid();
                List<MS2Modification> modListThisAcid = appropriateMap.get(acidString);
                if (modListThisAcid == null)
                {
                    modListThisAcid = new ArrayList<MS2Modification>();
                    appropriateMap.put(acidString, modListThisAcid);
                }
                modListThisAcid.add(mod);
            }


            for (int i=0; i<peptideSequence.length(); i++)
            {
                char thisResidue = peptideSequence.charAt(i);
                double staticMassThisResidue =
                        PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[thisResidue];

                List<MS2Modification> staticModListThisAcid = acidStaticModListMap.get("" + thisResidue);
                if (staticModListThisAcid != null)
                {
                    for (MS2Modification staticMod : staticModListThisAcid)
                    {
                        mass += staticMod.getMassDiff();
                        List<ModifiedAminoAcid> modsThisResidue = modifiedAminoAcids[i];
                        if (modsThisResidue == null)
                        {
                            modsThisResidue = new ArrayList<ModifiedAminoAcid>();
                            modifiedAminoAcids[i] = modsThisResidue;
                        }
                        staticMassThisResidue += staticMod.getMassDiff();
                        ModifiedAminoAcid moddedAcid = new ModifiedAminoAcid(thisResidue,
                                staticMassThisResidue);
                        modsThisResidue.add(moddedAcid);
                    }
                }

                List<MS2Modification> varModListThisAcid = acidVarModListMap.get("" + thisResidue);
                if (varModListThisAcid != null)
                {
                    for (MS2Modification varMod : varModListThisAcid)
                    {
                        mass += varMod.getMassDiff();
                        List<ModifiedAminoAcid> modsThisResidue = modifiedAminoAcids[i];
                        if (modsThisResidue == null)
                        {
                            modsThisResidue = new ArrayList<ModifiedAminoAcid>();
                            modifiedAminoAcids[i] = modsThisResidue;
                        }
                        ModifiedAminoAcid moddedAcid = new ModifiedAminoAcid(thisResidue,
                                staticMassThisResidue + varMod.getMassDiff());
                        modsThisResidue.add(moddedAcid);
                    }
                }
            }
        }

        return new Pair<List<ModifiedAminoAcid>[], Float>(modifiedAminoAcids, (float) mass);
    }
     */
}
