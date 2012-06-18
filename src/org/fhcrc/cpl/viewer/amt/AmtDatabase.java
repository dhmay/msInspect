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

import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;

/**
 * A representation of a full AMT database.  Stores everything we'd want
 * to read from or write to amtXml files, and allows us to manipulate the
 * database.
 *
 * An AMT database contains a hierarchical structure, exactly equivalent
 * to the structure defined in the amtXml schema.  Briefly:
 * -An AMT database contains a number of runs.
 * -An AMT database also contains a number of known Aminoacid Modifications,
 * which are referenced in the runs and in the modification states (see below)
 * -An AMT database also contains a number of peptide entries.
 * -Each peptide entry contains multiple Modification states.
 * -Each Modification State contains multiple observations.
 * -Each observation knows its hydrophobicity, quality score, and the run
 * it came from.
 *
 * Retention times for observations can come from MS/MS observations, or from matched
 * MS1 feature times (preferred).  This is controlled by AmtDatabaseBuilder and isn't
 * tracked in the database, so you'd best keep track yourself.
 *
 * There are obvious connections between AMT databases and peptides, and
 * Features and FeatureSets, but all references to Features are kept
 * out of this class.
 */
public class AmtDatabase
    implements Cloneable
{
    static Logger _log = Logger.getLogger(AmtDatabase.class);

    public static int DEFAULT_PRECISION=10;

    //Map of peptide entries.  Key = peptide sequence
    protected HashMap<String, AmtPeptideEntry> mAmtPeptideEntryMap = null;

    //List of runs.  Note:  The sequence of each run is its index + 1
    protected List<AmtRunEntry> mAmtRunEntries = null;

    //tolerance used when determining whether two modifications are the same.
    //Units are Da
    //TODO: should be able to make this MUCH smaller.  Can't.  What the heck?
    public static final float MODIFICATION_EQUALITY_MASS_TOLERANCE = .1f;

    //list of modifications.  Note:  The sequence of each mod is its index + 1
    protected List<MS2Modification> mAminoacidModifications = null;

    //Map from AMT run entries to their sequence number within this database.
    //For conversion to XML, which has no other way of linking observations
    //to runs
    protected HashMap<AmtRunEntry,Integer> mAmtRunSequenceMap = null;

    //Map from modifications to their sequence number within this database.
    //For conversion to XML, which has no other way of linking observations
    //to runs
    protected HashMap<MS2Modification,Integer> mAminoacidModificationSequenceMap = null;

    //File that the database was loaded from
    protected File mAmtDBSourceFile = null;

    //modification states' modified residue masses will be rounded to the nearest
    //multiple of this value when stored.  For bucketing slightly different masses
    public static final double DEFAULT_MODIFICATION_MASS_ROUNDING_FACTOR = 1;

    //defaults for hydrophobicity algorithm name and version
    public static final String DEFAULT_HYDROPHOBICITY_ALGORITHM_NAME = "krokhin";
    public static final double DEFAULT_HYDROPHOBICITY_ALGORITHM_VERSION = 3;

    protected String mHydrophobicityAlgorithmName =
            DEFAULT_HYDROPHOBICITY_ALGORITHM_NAME;
    protected double mHydrophobicityAlgorithmVersion =
            DEFAULT_HYDROPHOBICITY_ALGORITHM_VERSION;




    public AmtDatabase()
    {
        init();
    }

    /**
     * Initialize hashtables, etc
     */
    protected void init()
    {
        mAmtPeptideEntryMap = new HashMap<String, AmtPeptideEntry>();
        mAmtRunEntries = new ArrayList<AmtRunEntry>();
        mAmtRunSequenceMap = new HashMap<AmtRunEntry,Integer>();
        mAminoacidModificationSequenceMap = new HashMap<MS2Modification,Integer>();
        mAminoacidModifications = new ArrayList<MS2Modification>();
    }

    /**
     * VERY basic summary information
     * @return
     */
    public String toString()
    {
        return "(AMT Database with " + numRuns() + " runs, " + numAminoacidModifications() +
                " distinct modifications, and " + numEntries() + " entries)";
    }

    /**
     * Structure of database will be copied.  Individual entries, however, will be the same,
     * so don't mess with them.  Same with run entries.
     *
     * Hence, not a deep copy, not completely shallow.  Waist-deep.
     *
     * @return
     */
    public Object waistDeepCopy()
    {
        AmtDatabase cloneDb = new AmtDatabase();
        cloneDb.addOrOverrideEntriesWithAnotherDatabase(this);
        return cloneDb;
    }

    public MS2Modification findExistingEquivalentModification(MS2Modification newMS2Mod)
    {
//System.err.println("findequiv 1, newms2mod=" + newMS2Mod);
        if (mAminoacidModifications == null)
            return null;
        for (MS2Modification ms2Mod : mAminoacidModifications)
        {
//System.err.println("   findequiv: " + ms2Mod);
            if (ms2Mod.getAminoAcid().equalsIgnoreCase(newMS2Mod.getAminoAcid()) &&
               (Math.abs(ms2Mod.getMassDiff() - newMS2Mod.getMassDiff()) < MODIFICATION_EQUALITY_MASS_TOLERANCE) &&
               (ms2Mod.getVariable() == newMS2Mod.getVariable()))
            {
//System.err.println("  found equivalent");
                return ms2Mod;
            }
        }
        return null;
    }


    /**
     * Given an instance of a modification, figure out which MS2Modification within
     * this run it represents.  If none, return null
     * @return
     */
    public MS2Modification resolveMS2VariableModification(String residue, float massDiff,
                                                          AmtRunEntry runEntry)
    {
        _log.debug("resolveMS2VarMod 1, " + residue + ", " + massDiff);

        for (MS2Modification ms2Mod : runEntry.getVariableModifications())
        {
            _log.debug("  checking against " + ms2Mod);
            if (ms2Mod.getAminoAcid().equalsIgnoreCase(residue) &&
               (Math.abs(ms2Mod.getMassDiff() - massDiff) <
                       MODIFICATION_EQUALITY_MASS_TOLERANCE))
            {
                _log.debug("resolveMS2VarMods, found match: " + ms2Mod);
                return ms2Mod;
            }
        }
        return null;
    }




    /**
     * Add an observation.  Requires all the things that an observation needs
     * to know about, including the AmtRunEntry that it's associated with --
     * and that run entry must be a valid entry in THIS database
     * @param peptideSequence
     * @param modifiedAminoAcids
     * @param qualityScore
     * @param hydrophobicity
     * @param runEntry
     */
    public void resolveModsAndAddObservation(String peptideSequence,
                                             List<ModifiedAminoAcid>[] modifiedAminoAcids,
                                             double qualityScore,
                                             double hydrophobicity,
                                             AmtRunEntry runEntry,
                                             Map<String, Integer> spectralCountsMap,
                                             double timeInRun)
    {
        if (peptideSequence.contains("X"))
        {
            _log.debug("Skipping peptide sequence with aminoacid 'X'.  Peptide: " +
                    peptideSequence);
            return;
        }
        List<MS2Modification>[] ms2Modifications =
                resolveMods(peptideSequence, modifiedAminoAcids,
                            runEntry);
        addObservation(peptideSequence, ms2Modifications, qualityScore,
                       hydrophobicity, runEntry, spectralCountsMap, timeInRun);
    }

    /**
     * add an observation, having already resolved the modifications
     * @param peptideSequence
     * @param ms2Modifications
     * @param qualityScore
     * @param hydrophobicity
     * @param runEntry
     * @param spectralCountsMap
     * @param timeInRun
     */
    public void addObservation(String peptideSequence,
                                             List<MS2Modification>[] ms2Modifications,
                                             double qualityScore,
                                             double hydrophobicity,
                                             AmtRunEntry runEntry,
                                             Map<String, Integer> spectralCountsMap,
                                             double timeInRun)
    {
        if (spectralCountsMap == null)
            addObservation(peptideSequence, ms2Modifications, qualityScore,
                           hydrophobicity, runEntry,
                           AmtPeptideEntry.AmtPeptideObservation.SPECTRAL_COUNT_UNKNOWN,
                           timeInRun);
        else
            addObservation(peptideSequence, ms2Modifications, qualityScore,
                           hydrophobicity, runEntry, spectralCountsMap.get(peptideSequence),
                           timeInRun);
    }

    /**
     * Add an observation.  Requires all the things that an observation needs
     * to know about, including the AmtRunEntry that it's associated with --
     * and that run entry must be a valid entry in THIS database
     * @param peptideSequence
     * @param modifiedAminoAcids
     * @param runEntry
     */
    public List<MS2Modification>[] resolveMods(String peptideSequence,
                                             List<ModifiedAminoAcid>[] modifiedAminoAcids,
                                             AmtRunEntry runEntry)
    {
        _log.debug("resolveMods 1");

        List<MS2Modification>[] ms2Modifications =
                (List<MS2Modification>[]) new List[peptideSequence.length()];
        for (int i=0; i<ms2Modifications.length; i++)
            ms2Modifications[i] = new ArrayList<MS2Modification>();

        Map<String, Float> residueBaseMassMap = new HashMap<String, Float>();
        Map<String, List<MS2Modification>> residueStaticModMap =
                new HashMap<String, List<MS2Modification>>();

        //TODO: this should be moved up higher, not done for each peptide
        _log.debug("resolveMods, peptide " + peptideSequence + ", calculating static mods for each residue");
        for (int i=0; i<peptideSequence.length(); i++)
        {
            String residueString = "" + peptideSequence.charAt(i);
            if (!residueBaseMassMap.containsKey(residueString))
            {
                float residueBaseMass =
                        (float) PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[peptideSequence.charAt(i)];
                _log.debug("  Residue " + residueString + ", monoisotopic mass " + residueBaseMass);

                List<MS2Modification> staticModsThisResidue =
                        new ArrayList<MS2Modification>();
                residueStaticModMap.put(residueString, staticModsThisResidue);
                //resolve static modifications
                MS2Modification[] staticMods = runEntry.getStaticModifications();
                if (staticMods != null)
                {
                    for (MS2Modification staticMod : staticMods)
                    {
                        if (residueString.equalsIgnoreCase(staticMod.getAminoAcid()))
                        {
                            if (ms2Modifications[i] == null)
                                ms2Modifications[i] = new ArrayList<MS2Modification>();
                            staticModsThisResidue.add(staticMod);
                            residueBaseMass += staticMod.getMassDiff();
                            _log.debug("      Adding static mod of " + staticMod.getMassDiff());
                        }
                    }
                }
                _log.debug("    Added total " + staticModsThisResidue.size() + " modifications, residue base mass=" + residueBaseMass);
                residueBaseMassMap.put(residueString, residueBaseMass);
            }

            //add static mods calculated above
            ms2Modifications[i].addAll(residueStaticModMap.get(residueString));
        }

        //variable mods
        if (modifiedAminoAcids != null && modifiedAminoAcids.length > 0)
        {
            //variable modifications
            for (int i = 0; i < modifiedAminoAcids.length; i++)
            {
                if (modifiedAminoAcids[i] == null)
                    continue;

                for (ModifiedAminoAcid modifiedAminoAcid : modifiedAminoAcids[i])
                {
                    _log.debug("  Trying to add modification " + modifiedAminoAcids[i]);
                    String aminoAcidAsString = modifiedAminoAcid.getAminoAcidAsString();
                    float effectiveMassDiff =
                            (float) (modifiedAminoAcid.getMass() -
                                    residueBaseMassMap.get(aminoAcidAsString));
                    //if the effective mass diff is small, no need to account for it
                    if (Math.abs(effectiveMassDiff) > MODIFICATION_EQUALITY_MASS_TOLERANCE)
                    {
                        MS2Modification runMod =
                                resolveMS2VariableModification(aminoAcidAsString,
                                        effectiveMassDiff,
                                        runEntry);
                        if (runMod != null)
                            ms2Modifications[i].add(runMod);
                        else
                            throw new IllegalArgumentException("Unable to resolve variable modification on peptide " +
                                    peptideSequence + " : " + modifiedAminoAcid + ", base residue mass: " +
                                    residueBaseMassMap.get(aminoAcidAsString) + ", effective mass diff: " + effectiveMassDiff);
                    }

                }
            }
        }

        //if no mods at all, null out result
        //This is hokey.
        boolean atLeastOnePopulatedMod = false;
        for (List<MS2Modification> modList : ms2Modifications)
        {
            if (modList.size() > 0)
                atLeastOnePopulatedMod = true;
        }

        if (_log.isDebugEnabled())
        {
            StringBuffer annotatedPeptideSequence = new StringBuffer();
            for (int i=0; i<peptideSequence.length(); i++)
            {
                annotatedPeptideSequence.append(peptideSequence.charAt(i));
                for (MS2Modification mod : ms2Modifications[i])
                {
                    annotatedPeptideSequence.append("[+" + mod.getMassDiff());
                    if (mod.getVariable())
                        annotatedPeptideSequence.append("(V)");
                    annotatedPeptideSequence.append("]");
                }
            }
            _log.debug("*** " + annotatedPeptideSequence);
        }

        if (!atLeastOnePopulatedMod)
            ms2Modifications = null;


        return ms2Modifications;
    }


    /**
     * Add an observation.  Requires all the things that an observation needs
     * to know about, including the AmtRunEntry that it's associated with --
     * and that run entry must be a valid entry in THIS database
     * @param peptideSequence
     * @param qualityScore
     * @param hydrophobicity
     * @param runEntry
     */
    public void addObservation(String peptideSequence,
                               List<MS2Modification>[] ms2Modifications,
                               double qualityScore,
                               double hydrophobicity,
                               AmtRunEntry runEntry,
                               int spectralCount,
                               double timeInRun)
    {

        AmtPeptideEntry entryForPeptide = mAmtPeptideEntryMap.get(peptideSequence);

        AmtPeptideEntry.AmtPeptideObservation observation =
                AmtPeptideEntry.AmtPeptideObservation.createObservation(
                        hydrophobicity,
                        qualityScore,
                        runEntry,
                        timeInRun);

        observation.setSpectralCount(spectralCount);

        if (entryForPeptide == null)
        {
            entryForPeptide =
                    AmtPeptideEntry.createEntryFromObservation(
                            peptideSequence, ms2Modifications,
                            observation);
            mAmtPeptideEntryMap.put(peptideSequence,entryForPeptide);
        }
        else
        {
            entryForPeptide.addObservation(peptideSequence,
                                           ms2Modifications,
                                           observation);
        }
    }


    public void addObservationsFromEntry(AmtPeptideEntry newEntry)
    {
        addObservationsFromEntry(newEntry, null);
    }

    /**
     * Given an entry, add all the observations from all modification states in that entry
     * to the database.
     * If there's no existing entry for this peptide, create one
     * @param newEntry
     */
    public void addObservationsFromEntry(AmtPeptideEntry newEntry,
                                         Map<MS2Modification, MS2Modification> oldNewModMap)
    {
        String peptideSequence = newEntry.getPeptideSequence();
        AmtPeptideEntry entryForPeptide = mAmtPeptideEntryMap.get(peptideSequence);
        for (AmtPeptideEntry.AmtPeptideModificationStateEntry modEntry :
                newEntry.getModificationStateEntries())
        {
            if (modEntry.getModifications() != null)
            {
                for (List<MS2Modification> mods : modEntry.getModifications())
                {
                    if (mods == null)
                        continue;
                    if (oldNewModMap != null)
                        for (int i=0; i<mods.size(); i++)
                        {
                            if (oldNewModMap.containsKey(mods.get(i)))
                            {
                                mods.set(i, oldNewModMap.get(mods.get(i)));
                            }
                        }

                }
            }
        }
        if (entryForPeptide == null)
            mAmtPeptideEntryMap.put(peptideSequence, newEntry);
        else
        {
            for (AmtPeptideEntry.AmtPeptideModificationStateEntry modEntry :
                    newEntry.getModificationStateEntries())
            {
                entryForPeptide.addModificationStateEntry(modEntry);
            }
        }
    }

    /**
     * For each entry in another AmtDatabase, add all that entry's observations
     * to this database.  If this involves creating new entries here, so be it.  If entries
     * already exist, augment them with the new data
     * @param otherDatabase
     */
    public void addObservationsFromAnotherDatabase(AmtDatabase otherDatabase)
    {
        Map<MS2Modification, MS2Modification> oldNewModMap =
                new HashMap<MS2Modification, MS2Modification>();
        for (AmtRunEntry runEntry : otherDatabase.getRuns())
        {
            Map<MS2Modification, MS2Modification> oldNewModMapThisRun = 
                    addRunEntry(runEntry);
            oldNewModMap.putAll(oldNewModMapThisRun);
        }

        if (_log.isDebugEnabled())
        {
            for (MS2Modification oldMod : oldNewModMap.keySet())
            {
                _log.debug("Adding new run, modification map:  ");
                _log.debug("  " + oldMod + " -> " + oldNewModMap.get(oldMod) +
                           "  (ID " + getSequenceForAminoacidModification(oldNewModMap.get(oldMod)) + ")");
            }
        }

        for (AmtPeptideEntry otherDatabasePeptideEntry : otherDatabase.getEntries())
        {
            addObservationsFromEntry(otherDatabasePeptideEntry, oldNewModMap);
        }
    }

    /**
     * For each entry in another AmtDatabase, add it to this database,
     * overriding an existing entry if one exists.  All data from existing
     * overridden entries will be lost.
     *
     * No checking is done to determine if we have any orphan runs.
     * @param otherDatabase
     */
    public void addOrOverrideEntriesWithAnotherDatabase(AmtDatabase otherDatabase)
    {
        for (AmtRunEntry runEntry : otherDatabase.getRuns())
            addRunEntry(runEntry);
        for (AmtPeptideEntry otherDatabasePeptideEntry : otherDatabase.getEntries())
        {
            addOrOverrideEntry(otherDatabasePeptideEntry);
        }
    }



    /**
     * Save to a tsv file in some hokey format
     * @param tsvFile
     * @throws FileNotFoundException
     */
    public void saveToTsvSpreadsheet(File tsvFile) throws FileNotFoundException
    {
        PrintWriter pw = new PrintWriter(tsvFile);

        AmtRunEntry[] runs = getRuns();
        int numRuns = runs.length;

        Map<AmtRunEntry,Integer> runToIndexMap = new HashMap<AmtRunEntry,Integer>(numRuns);
        for (int i=0; i<numRuns; i++)
            runToIndexMap.put(runs[i],i);

        //write header line
        pw.write("sequence\tmass\tcalch\thaverage");

        for (int i=0; i<numRuns; i++)
            pw.write("\th" + i + "\tt" + i);
        pw.write("\n");

        double sentinelHValue = -999999;
        double[] hydrophobicitiesInRuns = new double[numRuns];
        double[] timesInRuns = new double[numRuns];

        for (AmtPeptideEntry peptideEntry : getEntries())
        {
            pw.write(peptideEntry.getPeptideSequence() + "\t" +
                    peptideEntry.getMass(PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES) + "\t" +
                    peptideEntry.getPredictedHydrophobicity() +
                     //"\t" + peptideEntry.getMass() +
                     "\t" + peptideEntry.getMedianObservedHydrophobicity());
            Arrays.fill(hydrophobicitiesInRuns,sentinelHValue);
            for (AmtPeptideEntry.AmtPeptideObservation peptideObservation : peptideEntry.getObservations())
            {
                hydrophobicitiesInRuns[runToIndexMap.get(peptideObservation.getRunEntry())] =
                        peptideObservation.getObservedHydrophobicity();
                timesInRuns[runToIndexMap.get(peptideObservation.getRunEntry())] =
                        peptideObservation.getTimeInRun();
            }
            for (int i=0; i<hydrophobicitiesInRuns.length; i++)
            {
                pw.write("\t");
                if (hydrophobicitiesInRuns[i] != sentinelHValue)
                    pw.write("" + hydrophobicitiesInRuns[i]);
                else
                    pw.write("NA");
                pw.write("\t");
                if (hydrophobicitiesInRuns[i] != sentinelHValue)
                {
                    pw.write("" + timesInRuns[i]);
                }
                else
                    pw.write("NA");
            }
            pw.write("\n");
        }
        pw.close();
    }

    //stuff related to peptide entries

    /**
     * This isn't free, because internally we store these as a HashMap
     * @return
     */
    public AmtPeptideEntry[] getEntries()
    {
        return mAmtPeptideEntryMap.values().toArray(new AmtPeptideEntry[mAmtPeptideEntryMap.size()]);
    }

    public AmtPeptideEntry[] getPeptideEntriesForRun(AmtRunEntry runEntry)
    {
        List<AmtPeptideEntry> resultList = new ArrayList<AmtPeptideEntry>();
        for (AmtPeptideEntry peptideEntry : getEntries())
        {
            if (peptideEntry.getObservationForRun(runEntry) != null)
               resultList.add(peptideEntry);
        }
        return resultList.toArray(new AmtPeptideEntry[resultList.size()]);
    }

    /**
     * this REALLY isn't free
     * @param runEntry
     * @return
     */
    public AmtPeptideEntry.AmtPeptideObservation[] getObservationsForRun(AmtRunEntry runEntry)
    {
        List<AmtPeptideEntry.AmtPeptideObservation> observations =
                new ArrayList<AmtPeptideEntry.AmtPeptideObservation>();
        for (AmtPeptideEntry peptideEntry : getEntries())
        {
            AmtPeptideEntry.AmtPeptideObservation obs =
                    peptideEntry.getObservationForRun(runEntry);
            if (obs != null)
                observations.add(obs);
        }
        return observations.toArray(new AmtPeptideEntry.AmtPeptideObservation[observations.size()]);
    }

    /**
     * Return the minimum time-in-run value for any observation in this run
     * @param runEntry
     * @return
     */
    public double getMinTimeInRun(AmtRunEntry runEntry)
    {
        double result = Double.MAX_VALUE;
        for (AmtPeptideEntry.AmtPeptideObservation obs : getObservationsForRun(runEntry))
        {
            if (obs.getTimeInRun() < result)
                result = obs.getTimeInRun();
        }
        return result;
    }

    /**
     * Return the maximum time-in-run value for any observation in this run
     * @param runEntry
     * @return
     */
    public double getMaxTimeInRun(AmtRunEntry runEntry)
    {
        double result = Double.MIN_VALUE;
        for (AmtPeptideEntry.AmtPeptideObservation obs : getObservationsForRun(runEntry))
        {
            if (obs.getTimeInRun() > result)
                result = obs.getTimeInRun();
        }
        return result;
    }

    public String[] getPeptides()
    {
        return mAmtPeptideEntryMap.keySet().toArray(new String[0]);
    }

    public AmtPeptideEntry getEntry(String peptideSequence)
    {
        return mAmtPeptideEntryMap.get(peptideSequence);
    }

    public boolean contains(String peptideSequence)
    {
        return getEntry(peptideSequence) != null;
    }

    /**
     * Add a peptide entry, blowing away the existing entry if
     * it was there.  Make no attempt to reconcile runs.
     * @param overridingEntry
     */
    protected void addOrOverrideEntry(AmtPeptideEntry overridingEntry)
    {
        mAmtPeptideEntryMap.put(overridingEntry.getPeptideSequence(),
                                overridingEntry);
    }

    /**
     * Remove an entry for a given sequence
     * @param peptideSequence
     */
    public void removeEntry(String peptideSequence)
    {
        mAmtPeptideEntryMap.remove(peptideSequence);
    }

    /**
     * count 'em
     * @return
     */
    public int numEntries()
    {
        return mAmtPeptideEntryMap.keySet().size();
    }


    //stuff related to aminoacid modifications
    /**
     * count 'em
     * @return
     */
    public int numAminoacidModifications()
    {
        if (mAminoacidModifications == null)
            return 0;
        return mAminoacidModifications.size();
    }

    /**
     * Add a run.  The sequence of this mod will be the new size of the ArrayList after addition
     */
    public void addAminoacidModification(MS2Modification newMod)
    {
//System.err.println("Adding mod: " + newMod);
        mAminoacidModifications.add(newMod);
        mAminoacidModificationSequenceMap.put(newMod,numAminoacidModifications());
    }

    /**
     * Get the mod with the specified sequence.
     * Note: sequence is one-based and ArrayLists are zero-based, so we subtract one
     * when referencing the ArrayList;
     * @param sequence (one-based)
     * @return
     */
    public MS2Modification getAminoacidModificationBySequence(int sequence)
    {
        if (sequence > numAminoacidModifications())
            return null;
        return mAminoacidModifications.get(sequence-1);
    }

    /**
     * Return an array containing all runs
     * @return
     */
    public MS2Modification[] getAminoacidModifications()
    {
        if (mAminoacidModifications == null)
            return null;
        return mAminoacidModifications.toArray(new MS2Modification[0]);
    }

    /**
     *
     * @return the sequence of this run entry, or -1 if it's not in the database
     */
    public int getSequenceForAminoacidModification(MS2Modification mod)
    {
        Integer modSequence = mAminoacidModificationSequenceMap.get(mod);
        if (modSequence == null)
            return -1;
        return modSequence;
    }


    

    //stuff related to runs

    /**
     * count 'em
     * @return
     */
    public int numRuns()
    {
        return mAmtRunEntries.size();
    }

    /**
     * Add a run.  The sequence of this run will be the new size of the ArrayList after addition
     * @param newRunEntry
     */
    public Map<MS2Modification, MS2Modification> addRunEntry(AmtRunEntry newRunEntry)
    {
        _log.debug("Adding run entry");
        _log.debug("Before add, num mods: " + this.numAminoacidModifications());

        //override all modifications on the run entry that are already in the
        // database with their equivalents in the database.
        //So that everything agrees BEFORE
        //we start adding observations
        Map<MS2Modification, MS2Modification> oldNewModMap =
                newRunEntry.overrideDuplicateModifications(this);



        //check for any modifications that aren't already in the database,
        //add them all
        //For ones that were there, update the object reference in the run
        if (newRunEntry.getModifications() != null)
        {
            MS2Modification[] newRunModArray = newRunEntry.getModifications();
            _log.debug("Run entry has " + newRunEntry.getModifications().length + " modifications");
            _log.debug("Existing database has " + this.getAminoacidModifications().length + " modifications");

            for (int i=0; i<newRunModArray.length; i++)
            {
                MS2Modification newMs2Mod = newRunModArray[i];
                _log.debug("Evaluating modification for inclusion: " + newMs2Mod);

                boolean alreadyHasMod = false;
                for (MS2Modification existingMod : getAminoacidModifications())
                {
                    _log.debug("  Comparing against: " + existingMod);
                    if (existingMod == newMs2Mod)
                    {
//                        oldNewModMap.put(newRunModArray[i], existingMod);
                        newRunModArray[i] = existingMod;
                        alreadyHasMod = true;
                        _log.debug(" Found it already there");
                        break;
                    }
                }
                if (!alreadyHasMod)
                {
                    addAminoacidModification(newMs2Mod);
                    _log.debug("Inserting it, now we've got " + this.numAminoacidModifications());
                }
            }
            newRunEntry.setModifications(newRunModArray);
        }
        mAmtRunEntries.add(newRunEntry);
        mAmtRunSequenceMap.put(newRunEntry,numRuns());

        _log.debug("After run add, num mods: " + this.numAminoacidModifications());


        return oldNewModMap;
    }

    /**
     * Get the run with the specified sequence.
     * Note: sequence is one-based and ArrayLists are zero-based, so we subtract one
     * when referencing the ArrayList;
     * @param sequence (one-based)
     * @return
     */
    public AmtRunEntry getRunBySequence(int sequence)
    {
        if (sequence > numRuns())
            return null;
        return mAmtRunEntries.get(sequence-1);
    }

    /**
     * Return an array containing all runs
     * @return
     */
    public AmtRunEntry[] getRuns()
    {
        return mAmtRunEntries.toArray(new AmtRunEntry[0]);
    }

    /**
     *
     * @param runEntry
     * @return the sequence of this run entry, or -1 if it's not in the database
     */
    public int getSequenceForRun(AmtRunEntry runEntry)
    {
        Integer runSequence = mAmtRunSequenceMap.get(runEntry);
        if (runSequence == null)
            return -1;
        return runSequence;
    }

    //Consider moving this stuff to a new diagnostic class

    /**
     * Calculate the mean difference of all median peptide H observations from prediction
     * @return
     */
    public double calculateMeanDifferenceFromPredictedHydro()
    {
        double[] medianHydDifferencesFromPredicted = new double[numEntries()];
        int i = 0;
        for (AmtPeptideEntry entry : getEntries())
        {
            medianHydDifferencesFromPredicted[i++] =
                    Math.abs(entry.getMedianObservedHydrophobicity() - entry.getPredictedHydrophobicity());
        }
        return BasicStatistics.mean(medianHydDifferencesFromPredicted);
    }

    /**
     * Calculate the standard deviation of all deviations of median peptide H observations from prediction
     * @return
     */
    public double calculateStandardDeviationDifferenceFromPredictedHydro()
    {
        double[] medianHydDeviationsFromPredicted = new double[numEntries()];
        int i = 0;
        for (AmtPeptideEntry entry : getEntries())
        {
            medianHydDeviationsFromPredicted[i++] =
                    Math.abs(entry.getMedianObservedHydrophobicity() -
                             entry.getPredictedHydrophobicity());
        }
        return BasicStatistics.standardDeviation(medianHydDeviationsFromPredicted);
    }





    public String getHydrophobicityAlgorithmName()
    {
        return mHydrophobicityAlgorithmName;
    }

    public void setHydrophobicityAlgorithmName(String hydrophobicityAlgorithmName)
    {
        this.mHydrophobicityAlgorithmName = hydrophobicityAlgorithmName;
    }

    public double getHydrophobicityAlgorithmVersion()
    {
        return mHydrophobicityAlgorithmVersion;
    }

    public void setHydrophobicityAlgorithmVersion(double hydrophobicityAlgorithmVersion)
    {
        this.mHydrophobicityAlgorithmVersion = hydrophobicityAlgorithmVersion;
    }


    public File getAmtDBSourceFile()
    {
        return mAmtDBSourceFile;
    }

    public void setAmtDBSourceFile(File mAmtDBSourceFile)
    {
        this.mAmtDBSourceFile = mAmtDBSourceFile;
    }
}
