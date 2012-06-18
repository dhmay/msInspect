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
package org.fhcrc.cpl.viewer.amt;

import java.util.*;
import java.util.List;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;


/**
 * This class Generates FeatureSets based on AMT Database entries.  These featuresets
 * are what's actually used in matching.
 */
public class AmtDatabaseFeatureSetGenerator
{
    static Logger _log = Logger.getLogger(AmtDatabaseFeatureSetGenerator.class);

    protected AmtDatabase amtDatabase;

    public AmtDatabaseFeatureSetGenerator(AmtDatabase amtDatabase)
    {
        this.amtDatabase = amtDatabase;
    }

    /**
     * Create a feature set based on this AMT database, accounting for expected
     * modifications.
     * That is, bump up the mass by the mass of any static modifications, and
     * create duplicate features for any variable modifictions
     * TODO: maybe interrogate modifications[] to see if any of these match up with
     * TODO: the mods we've already got in the db, and get the hydrophobicity accordingly
     *
     * @param modifications
     * @return
     */
    public Feature[] createFeaturesForModifications(MS2Modification[] modifications)
    {
        List<MS2Modification> varModList = new ArrayList<MS2Modification>();
        List<MS2Modification> staticModList = new ArrayList<MS2Modification>();
        if (modifications != null)
            for (MS2Modification modification : modifications)
            {
                if (modification.getVariable())
                    varModList.add(modification);
                else
                    staticModList.add(modification);
            }
//for (MS2Modification mod : varModList) System.err.println(mod);


        List<Feature> resultList = new ArrayList<Feature>();


        for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
        {
            //This call is a bit confusing.  We start off the recursive call declaring that all the static
            //mods should be applied, and supplying the variable mod list as the list of mods to potentially
            //apply.  Later recursive calls will or will not add variable mods to this list, one by one.
            List<Feature> featuresForThisPeptide =
                    generateModFeaturesForPeptide(peptideEntry,
                                                  peptideEntry.getMedianObservedHydrophobicity(),
                                                  staticModList, varModList);

            resultList.addAll(featuresForThisPeptide);
            
        }

        return resultList.toArray(new Feature[resultList.size()]);
    }



    /**
     * Note: Assumes that a variable modification is either fully applied, to all
     * residues in the peptide, or not applied at all. So, e.g., ELVISMMM either has all
     * oxidized M's or none.
     * @param peptideEntry
     * @param observedHydrophobicity
     * @param staticMods
     * @param varMods
     * @return
     */
    public static List<Feature>generateModFeaturesForPeptide(AmtPeptideEntry peptideEntry,
                                                      double observedHydrophobicity,
                                                      List<MS2Modification> staticMods,
                                                      List<MS2Modification> varMods)
    {
        String peptideSequence = peptideEntry.getPeptideSequence();

        List<MS2Modification> staticModsThisFeature =
                new ArrayList<MS2Modification>();
        for (MS2Modification mod : staticMods)
        {
            if (peptideSequence.contains(mod.getAminoAcid()))
                staticModsThisFeature.add(mod);
        }

        List<MS2Modification> varModsThisFeature =
                new ArrayList<MS2Modification>();
        for (MS2Modification mod : varMods)
        {
            if (peptideSequence.contains(mod.getAminoAcid()))
                varModsThisFeature.add(mod);
        }

//System.err.println(peptideSequence + ", static: " + staticModsThisFeature + ", var: " + varModsThisFeature);        
        return recursivelyAddFeaturesForMods(peptideEntry,
                observedHydrophobicity, staticModsThisFeature, varModsThisFeature);
    }

    /**
     * Recursively generate features for all possible masses, given the list of
     * /variable/ modifications known to exist in this peptide and the list of modifications
     * already applied
     * @param peptideEntry
     * @param observedHydrophobicity
     * @param appliedMods already-applied modifications.  This will include all static mods and a growing
     * list of variable mods
     * @param varModsRemaining  variable mods remaining.  None of these mods should ever be static
     * @return
     */
    protected static List<Feature> recursivelyAddFeaturesForMods(AmtPeptideEntry peptideEntry,
                                                          double observedHydrophobicity,
                                                          List<MS2Modification> appliedMods,
                                                          List<MS2Modification> varModsRemaining)
    {
        List<Feature> result = new ArrayList<Feature>();
        if (varModsRemaining.size() == 0)
        {
            Feature feature = createFeatureForPeptideWithMods(peptideEntry,
                    observedHydrophobicity, appliedMods);
            result.add(feature);
            return result;
        }

        //This is hacky and wasteful.
        //TODO: find a better way to deal with keeping the integrity of these lists through recursion
        List<MS2Modification> modsRemainingCopy =
                new ArrayList<MS2Modification>(varModsRemaining.size());
        modsRemainingCopy.addAll(varModsRemaining);
        varModsRemaining = modsRemainingCopy;

        List<MS2Modification> appliedModsCopy =
                new ArrayList<MS2Modification>(appliedMods.size());
        appliedModsCopy.addAll(appliedMods);
        appliedMods = appliedModsCopy;


        MS2Modification mod = varModsRemaining.get(0);
        modsRemainingCopy.remove(mod);

        result.addAll(recursivelyAddFeaturesForMods(peptideEntry, observedHydrophobicity,
                                                    appliedMods, varModsRemaining));
        appliedMods.add(mod);
        result.addAll(recursivelyAddFeaturesForMods(peptideEntry, observedHydrophobicity,
                                                    appliedMods,
                                                    varModsRemaining));
        return result;
    }

    /**
     * Modifications are applied all-or-nothing.  If we modify one residue, we modify all
     * of that residue
     * @param peptideEntry
     * @param observedHydrophobicity
     * @param modifications
     * @return
     */
    protected static Feature createFeatureForPeptideWithMods(AmtPeptideEntry peptideEntry,
                                                    double observedHydrophobicity,
                                                    List<MS2Modification> modifications)
    {
        String peptideSequence = peptideEntry.getPeptideSequence();


        Feature feature = new Feature();
        MS2ExtraInfoDef.addPeptide(feature, peptideSequence);
        MS2ExtraInfoDef.setPeptideProphet(feature, peptideEntry.calculateIDProbability());
        feature.setPeaks(1);
        feature.setScanCount(1);

//Hack! Hack! Hack!  I've been using this to be able to match two of these sets
// to each other.  Really there's no good value for scan, though, so, hey, whatever.         
        feature.setScan((int) ((1000 * (observedHydrophobicity)) + 2000));
        feature.setIntensity(1000);
        AmtExtraInfoDef.setObservedHydrophobicity(feature, observedHydrophobicity);

        double mass = PeptideGenerator.computeMass(peptideSequence.getBytes(),
                                                   0, peptideSequence.length(),
                                                   PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
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

            List<ModifiedAminoAcid>[] modifiedAminoAcids = new ArrayList[peptideSequence.length()];
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
            MS2ExtraInfoDef.setModifiedAminoAcids(feature, modifiedAminoAcids);
        }
        feature.setMass((float) mass);
        feature.updateMz();

        return feature;
    }

    /**
     * create a featureset that represents all the features from a run within an
     * amt database, with their original times.
     *
     * Note: scan will be incorrectB.
     * @param runEntry
     * @param modifications
     * @return
     */
    public FeatureSet createFeatureSetForRun(AmtRunEntry runEntry,
                                             MS2Modification[] modifications)
    {

        List<MS2Modification> varModList = new ArrayList<MS2Modification>();
        List<MS2Modification> staticModList = new ArrayList<MS2Modification>();

        if (modifications != null)
            for (MS2Modification modification : modifications)
            {
                if (modification.getVariable())
                    varModList.add(modification);
                else
                    staticModList.add(modification);
            }


        List<Feature> resultList = new ArrayList<Feature>();
        for (AmtPeptideEntry peptideEntry : amtDatabase.getPeptideEntriesForRun(runEntry))
        {
            AmtPeptideEntry.AmtPeptideObservation obs =
                    peptideEntry.getObservationForRun(runEntry);
            List<Feature> featuresForThisObservation =
                    generateModFeaturesForPeptide(peptideEntry,
                                                  obs.getObservedHydrophobicity(),
                                                  staticModList, varModList);
            for (Feature feature : featuresForThisObservation)
            {
                feature.setTime((float) obs.getTimeInRun());
                resultList.add(feature);
            }
        }

        return new FeatureSet(resultList.toArray(new Feature[resultList.size()]));
    }



}
