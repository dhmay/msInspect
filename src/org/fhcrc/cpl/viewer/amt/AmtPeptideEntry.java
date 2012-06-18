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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
//import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
//import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;

/**
 * Encapsulates everything we need to retrieve from / store into an AMT database.
 * Awfully redundant with the XmlBeans amtXml representation.  It's worth it to have a
 * separate class so that we have complete control over the code and can do things like
 * calculate predicted hydrophobicities.
 */
public class AmtPeptideEntry
{
    static Logger _log = Logger.getLogger(AmtPeptideEntry.class);

    protected String mPeptideSequence = null;
    protected Peptide mPeptide = null;
    protected double mPredictedHydrophobicity = 0.0;
    protected double mMedianObservedHydrophobicity = 0.0;
    protected double mMedianPeptideProphet = 0.0;
    protected double mHydrophobicityStandardDeviation = 0.0;
    protected double mMassUsingMonoisotopicMassTable = -1;
    protected double mMeanObservedHydrophobicity = MEAN_HYDRO_UNCALCULATED_SENTINEL;

    //value that indicates that the mean observed hydro hasn't been calculated yet
    protected static final double MEAN_HYDRO_UNCALCULATED_SENTINEL=-99999999;

    //a map of modified peptides to entries for those modified peptides
    protected HashMap<String,AmtPeptideModificationStateEntry> mModificationStateEntryMap = null;

    public AmtPeptideEntry()
    {
        mModificationStateEntryMap = new HashMap<String,AmtPeptideModificationStateEntry>();
    }


    public String toString()
    {
        return "AmtPeptideEntry: sequence="+ mPeptideSequence +
               ", medianObservedHydro="+ getMedianObservedHydrophobicity() +
               ", observations=" + getNumObservations();
    }

    public static AmtPeptideEntry createEntryFromObservation(String peptideSequence,
                                                             List<MS2Modification>[] modifications,
                                                             AmtPeptideObservation observation)

    {
        return createEntryFromObservation(peptideSequence, modifications,
                                          AmtUtilities.calculateNormalizedHydrophobicity(peptideSequence),
                                          observation);
    }

    /**
     * Create a new entry with one observation based on a calculated hydrophobicity and
     * that observation
     * @param predictedHydrophobicity
     * @param observation
     * @return
     */
    public static AmtPeptideEntry createEntryFromObservation(String peptideSequence,
                                                             List<MS2Modification>[] modifications,
                                                             double predictedHydrophobicity,
                                                             AmtPeptideObservation observation)

    {
        AmtPeptideEntry newEntry = new AmtPeptideEntry();
        newEntry.setPeptideSequence(peptideSequence);
        newEntry.setPredictedHydrophobicity(predictedHydrophobicity);

        AmtPeptideModificationStateEntry modEntry =
                AmtPeptideModificationStateEntry.createEntryFromObservation(peptideSequence,
                                                                            modifications,
                                                                            observation,
                                                                            newEntry);
        newEntry.addModificationStateEntry(modEntry);

        return newEntry;
    }


    /**
     * Round to the nearest multiple of modificationMassRoundingFactor.
     * Do it by dividing by the factor, rounding to nearest 1, then multiplying
     * by the factor
     * @param inputMass
     * @return
     */
    protected static double roundModifiedMass(double inputMass)
    {
        double result = inputMass / AmtDatabase.DEFAULT_MODIFICATION_MASS_ROUNDING_FACTOR;
        result = Math.round(result);
        result *= AmtDatabase.DEFAULT_MODIFICATION_MASS_ROUNDING_FACTOR;
        return result;
    }

        /**
         * Add a bunch of observations for this peptide, without recalculating stats,
         * then recalculate once.  Because recalculation is expensive
         *
         * @param observations
         */
        public void addObservations(String peptide,
                                    List<MS2Modification>[] modifications,
                                    AmtPeptideObservation[] observations,
                                    double modificationMassRoundingFactor)
        {
            for (AmtPeptideObservation observation : observations)
                addObservationNoRecalc(peptide, modifications, observation);
            recalculateStats(AmtDatabase.DEFAULT_PRECISION);
        }

        /**
         * Add an observation for this peptide and recalculate the stats
         *
         * @param observation
         */
        public void addObservation(String peptide,
                                   List<MS2Modification>[] modifications,
                                   AmtPeptideObservation observation)
        {
            addObservationNoRecalc(peptide, modifications, observation);
            recalculateStats();
        }


        /**
         * Add an observation for this peptide.  Don't recalculate stats
         *
         * @param observation
         */
        protected void addObservationNoRecalc(String peptide,
                                              List<MS2Modification>[] modifications,
                                              AmtPeptideObservation observation)
        {
        String modifiedPeptideSequence =
                AmtPeptideModificationStateEntry.calculateModifiedPeptideSequence(peptide,
                        modifications);
            AmtPeptideModificationStateEntry modEntry =
                    getModificationStateEntry(modifiedPeptideSequence);
            if (modEntry == null)
            {

                modEntry =
                        AmtPeptideModificationStateEntry.createModificationStateEntryFromObservation(
                                modifiedPeptideSequence,
                                modifications,
                                observation,
                                this);
                mModificationStateEntryMap.put(modifiedPeptideSequence, modEntry);
            }
            else
                modEntry.addObservationNoRecalc(observation);
        }

    /**
     * Build this on demand only
     * @return
     */
    public Peptide getPeptide()
    {
        if (mPeptide == null)
            mPeptide = getPeptideForPeptideSequence(this.getPeptideSequence());
        return mPeptide;
    }

    /**
     * Calculate this on demand only
     * @return
     */
    public double getMass()
    {
        if (mMassUsingMonoisotopicMassTable < 0)
            mMassUsingMonoisotopicMassTable = getPeptide().getMass(PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
        return mMassUsingMonoisotopicMassTable;
    }

    /**
     * Assume this is a table other than monoisotopic mass table, actually calculate
     * @param massTable
     * @return
     */
    public double getMass(double[] massTable)
    {
        return getPeptide().getMass(massTable);
    }

    public static Peptide getPeptideForPeptideSequence(String peptideSequence)
    {
        //this is lame.  Really peptide should have a constructor that doesn't require a protein
        Protein fakeProtein = new Protein("",peptideSequence.getBytes());
        return new Peptide(fakeProtein,0,fakeProtein.getBytes().length);
    }



    /**
     * predict hydrophobicity based on a peptide sequence.  This isn't free,
     * so we want to avoid doing it too often, which explains the structure of some of this
     * code
     * @param peptideSequence
     * @return
     */
    public static double predictHydrophobicityForPeptideSequence(String peptideSequence)
    {
        return getPeptideForPeptideSequence(peptideSequence).getHydrophobicity();
    }


    /**
     * Recalculate the mean observed hydrophobicity and other stats, munging together all
     * observations, giving them equal weight.
     * Possibly, in the future, if we store probability on observations, weight them
     */
    public void recalculateStats(int precision)
    {
        AmtPeptideObservation[] observations = getObservations();

        double[] observedPeptideProphetScores = new double[observations.length];
        double[] observedHydrophobicities = new double[observations.length];
        for (int i = 0; i < observations.length; i++)
        {
            observedHydrophobicities[i] = observations[i].getObservedHydrophobicity();
            observedPeptideProphetScores[i] = observations[i].getPeptideProphet();
        }
        mMedianObservedHydrophobicity =
                Rounder.round(BasicStatistics.median(observedHydrophobicities), precision);
        mMedianPeptideProphet =
                Rounder.round(BasicStatistics.median(observedPeptideProphetScores), precision);
        mHydrophobicityStandardDeviation =
                Rounder.round(BasicStatistics.standardDeviation(observedHydrophobicities), precision);
    }

    public void recalculateStats()
    {
        recalculateStats(AmtDatabase.DEFAULT_PRECISION);
    }

    /**
     * Calculate the probability that the ID for this peptide entry is correct.
     * This is calculated simply as:
     * 1 - the product of (1 - each peptide observation's probability)
     *
     * @return
     */
    public float calculateIDProbability()
    {
        float productOfBadIDProbabilities = 1;

        for (AmtPeptideObservation obs : getObservations())
            productOfBadIDProbabilities *= (1 - obs.getPeptideProphet());
        return 1 - productOfBadIDProbabilities;
    }

    /**
     * Quick and dirty method to find which observation, if any, came from this run entry.
     * If nothing, return null.
     * TODO: If this is something we do often, I should build a datastructure for this.
     * TODO: This is totally broken with respect  to multiple modification states
     * @param runEntry
     * @return
     */
    public AmtPeptideObservation getObservationForRun(AmtRunEntry runEntry)
    {
        for (AmtPeptideObservation observation : getObservations())
        {
            if (observation.getRunEntry() == runEntry)
                return observation;
        }
        return null;
    }

    /**
     * how many observations for this peptide entry?
     *
     * @return
     */
    public int getNumObservations()
    {
        int numObservations = 0;
        for (AmtPeptideModificationStateEntry modEntry : getModificationStateEntries())
        {
            numObservations += modEntry.getNumObservations();
        }
        return numObservations;
    }

    /**
     * how many observations for this peptide entry?
     * @return
     */
    public int getNumModificationStates()
    {
        return mModificationStateEntryMap.size();
    }

    //getters and setters


    public double getHydrophobicityStandardDeviation()
    {
        return mHydrophobicityStandardDeviation;
    }

    public void setHydrophobicityStandardDeviation(double newValue)
    {
        mHydrophobicityStandardDeviation = newValue;
    }

    /**
     * This getter is computationally intensive, since observations are stored one level down
     * @return
     */
    public AmtPeptideObservation[] getObservations()
    {
        ArrayList<AmtPeptideObservation> observationList = new ArrayList<AmtPeptideObservation>();
        for (AmtPeptideModificationStateEntry modStateEntry : mModificationStateEntryMap.values())
        {
            for (AmtPeptideObservation observation : modStateEntry.getObservations())
            observationList.add(observation);
        }
        return observationList.toArray(new AmtPeptideObservation[0]);
    }

    /**
     * Removes an observation
     * @return true if the observation was in this entry already
     */
    public boolean removeObservation(AmtPeptideObservation observationToRemove)
    {
        for (String modStateString : mModificationStateEntryMap.keySet())
        {
            AmtPeptideModificationStateEntry modStateEntry = mModificationStateEntryMap.get(modStateString);
            if (modStateEntry.removeObservation(observationToRemove))
            {
                AmtPeptideObservation[] observations = modStateEntry.getObservations();
                if (observations == null || observations.length == 0)
                    mModificationStateEntryMap.remove(modStateString);
                recalculateStats();
                return true;
            }
        }
        return false;
    }


    public AmtPeptideModificationStateEntry getModificationStateEntry(String modifiedSequence)
    {
        return mModificationStateEntryMap.get(modifiedSequence);
    }

    public AmtPeptideModificationStateEntry[] getModificationStateEntries()
    {
        return mModificationStateEntryMap.values().toArray(new AmtPeptideModificationStateEntry[0]);
    }

    public AmtPeptideModificationStateEntry
            addModificationStateEntry(String modifiedPeptideSequence,
                                      double modifiedMass,
                                      List<MS2Modification>[] modifications)
    {
        AmtPeptideModificationStateEntry newModificationStateEntry =
                AmtPeptideModificationStateEntry.createModificationStateEntry(
                        modifiedPeptideSequence, modifiedMass, modifications, this);
        addModificationStateEntry(newModificationStateEntry);
        return newModificationStateEntry;
    }

    /**
     * Add a modification state entry.  If we've already got one for this modified peptide
     * (after re-rounding), add all of the new entry's observations to the old entry
     * @param newModificationStateEntry
     */
    public void addModificationStateEntry(AmtPeptideModificationStateEntry newModificationStateEntry)
    {
        newModificationStateEntry.recalculateModifiedSequenceAndMass(getPeptideSequence());
        newModificationStateEntry.setParentPeptideEntry(this);

        AmtPeptideModificationStateEntry modEntry =
                getModificationStateEntry(newModificationStateEntry.getModifiedSequence());
        if (modEntry == null)
        {
            mModificationStateEntryMap.put(newModificationStateEntry.getModifiedSequence(),
                                           newModificationStateEntry);
        }
        else
        {
            modEntry.addObservations(newModificationStateEntry.getObservations());
        }

        recalculateStats();
    }

    public String getPeptideSequence()
    {
        return mPeptideSequence;
    }

    public void setPeptideSequence(String peptideSequence)
    {
        this.mPeptideSequence = peptideSequence;
    }

    public double getPredictedHydrophobicity()
    {
        return mPredictedHydrophobicity;
    }

    public void setPredictedHydrophobicity(double predictedHydrophobicity)
    {
        this.mPredictedHydrophobicity = predictedHydrophobicity;
    }

    public double getMedianObservedHydrophobicity()
    {
        return mMedianObservedHydrophobicity;
    }

    public double getMedianPeptideProphet()
    {
        return mMedianPeptideProphet;
    }

    /**
     * This is not a cost-free getter.  We store median, not mean, at this level, so if we need
     * the mean we have to calculate it.  Once calculated, store it.
     * @return
     */
    public double getMeanObservedHydrophobicity()
    {
        if (mMeanObservedHydrophobicity == MEAN_HYDRO_UNCALCULATED_SENTINEL)
        {
            double[] observedHydrophobicities = new double[getNumObservations()];
            for (int i=0; i<getNumObservations(); i++)
            {
                observedHydrophobicities[i] = getObservations()[i].getObservedHydrophobicity();
            }
            mMeanObservedHydrophobicity = BasicStatistics.mean(observedHydrophobicities);
        }
        return mMeanObservedHydrophobicity;
    }

    public int getSpectralCount()
    {
        int count = 0;
        for (AmtPeptideObservation observation : getObservations())
        {
            if (!observation.hasSpectralCount())
                return AmtPeptideObservation.SPECTRAL_COUNT_UNKNOWN;
            count += observation.getSpectralCount();
        }
        return count;
    }



    /**
     * Represents an entry for a specific modification state
     */
    public static class AmtPeptideModificationStateEntry
    {
        ArrayList<AmtPeptideObservation> mObservations = new ArrayList<AmtPeptideObservation>();

        List<MS2Modification>[] mModifications = null;
        double mModifiedMass = 0.0;
        double mMedianPeptideProphet = 0.0;
        double mHydrophobicityStandardDeviation = 0.0;
        double mMedianObservedHydrophobicity = 0.0;
        String mModifiedSequence = null;
        AmtPeptideEntry mParentPeptideEntry = null;

        public String toString()
        {
            return "AmtPeptideModificationStateEntry: modifiedSequence=" + mModifiedSequence +
                    ", numObservations = " + getNumObservations();
        }

//        public boolean hasSameModifications(List<MS2Modification>[] otherModArray)
//        {
//            if (mModifications == null && otherModArray == null)
//                return true;
//            //both are not null, so if one is, not same
//            if ((mModifications == null || otherModArray == null))
//                return false;
//            if (otherModArray.length != mModifications.length)
//            for (int i=0; i<mModifications.length; i++)
//            {
//                if (mModifications[i] == null && otherModArray[i] == null)
//                    continue;
//                //both are not null, so if one is, not same
//                if ((mModifications[i] == null || otherModArray[i] == null))
//                    return false;
//                if (mModifications[i].size() != otherModArray[i].size())
//                    return false;
//                for (int j=0; j<mModifications[i].size(); i++)
//                    if (!(mModifications[i].equals(otherModArray[i])))
//                        return false;
//            }
//            //if we got here, every index compared fine
//            return true;
//        }

        public AmtPeptideObservation getObservationForRun(AmtRunEntry runEntry)
        {
            for (AmtPeptideObservation observation : getObservations())
            {
                if (observation.getRunEntry() == runEntry)
                    return observation;
            }
            return null;
        }


        /**
         * Create a new entry with one observation based on a calculated
         * hydrophobicity and that observation
         * @param modifications
         * @param observation
         * @return
         */
        public static AmtPeptideModificationStateEntry
                createModificationStateEntryFromObservation(String modifiedPeptideSequence,
                                                            List<MS2Modification>[] modifications,
                                                            AmtPeptideObservation observation,
                                                            AmtPeptideEntry parentPeptideEntry)

        {
            AmtPeptideModificationStateEntry newModificationStateEntry =
                    createModificationStateEntry(modifiedPeptideSequence,
                                             -1,
                                             modifications,
                                             parentPeptideEntry);
            newModificationStateEntry.addObservation(observation);
            newModificationStateEntry.recalculateMass();

            return newModificationStateEntry;
        }

        /**
         * Create a new entry with one observation based on a calculated hydrophobicity and
         * that observation
         * @param modifications
         * @param observation
         * @return
         */
        public static AmtPeptideModificationStateEntry
                createEntryFromObservation(String peptideSequence,
                                           List<MS2Modification>[] modifications,
                                           AmtPeptideObservation observation,
                                           AmtPeptideEntry parentPeptideEntry)

        {
            String modifiedPeptideSequence =
                    calculateModifiedPeptideSequence(peptideSequence,
                    modifications);
            return createModificationStateEntryFromObservation(modifiedPeptideSequence, modifications,
                                              observation, parentPeptideEntry);
        }

        /**
         * Create a new entry with no observations
         * @param modifications
         * @return
         */
        public static AmtPeptideModificationStateEntry
                createModificationStateEntry(String modifiedPeptideSequence,
                                             double modifiedMass,
                                             List<MS2Modification>[] modifications,
                                             AmtPeptideEntry parentPeptideEntry)

        {
            AmtPeptideModificationStateEntry newModificationStateEntry = new AmtPeptideModificationStateEntry();
            newModificationStateEntry.setModifiedSequence(modifiedPeptideSequence);
            newModificationStateEntry.setModifications(modifications);
            newModificationStateEntry.setParentPeptideEntry(parentPeptideEntry);
            if (modifiedMass > 0)
                newModificationStateEntry.setModifiedMass(modifiedMass);
            else
                newModificationStateEntry.recalculateMass();
            return newModificationStateEntry;
        }

        public int getSpectralCount()
        {
            int count = 0;
            for (AmtPeptideObservation observation : getObservations())
            {
                if (!observation.hasSpectralCount())
                    return AmtPeptideObservation.SPECTRAL_COUNT_UNKNOWN;
                else count++;
            }
            return count;
        }

        /**
         * Look through all the modifications on this peptide and account
         * for all of them
         * @param peptideSequence
         * @param modifications
         * @return
         */
        protected static String calculateModifiedPeptideSequence(String peptideSequence,
                                                                 List<MS2Modification>[] modifications)
        {
            if (modifications == null)
                return peptideSequence;
            StringBuffer modifiedPeptideSequenceBuf = new StringBuffer();

            for (int i = 0; i < peptideSequence.length(); i++)
            {
                char residue = peptideSequence.charAt(i);
                modifiedPeptideSequenceBuf.append(residue);
                if (modifications[i] != null && modifications[i].size() > 0)
                {
                    double newMass = PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[residue];
                    for (MS2Modification mod : modifications[i])
                    {
                        newMass += mod.getMassDiff();
                    }
                    double roundedMass = roundModifiedMass(newMass);
                    //TODO: will this ever produce trailing zeroes?  Shouldn't
                    modifiedPeptideSequenceBuf.append("[" + roundedMass + "]");
                }
            }
            return modifiedPeptideSequenceBuf.toString();
        }

        /**
         * Recalculate the modified sequence and mass.  Used when wanting to change the rounding
         * factor
         *
         * @param peptideSequence
         */
        public void recalculateModifiedSequenceAndMass(String peptideSequence)
        {
            mModifiedSequence =
                    calculateModifiedPeptideSequence(peptideSequence,
                    mModifications);
            recalculateMass();
        }

        /**
         * Recalculate the mass of this modification state in the following complicated way:
         * -For all modified aminoacids, total up the masses
         * -For all unmodified aminoacids, total up the masses by creating a fake Peptide and
         * asking it for its mass
         * -add the two totals together
         */
        public void recalculateMass()
        {
            StringBuffer sequenceCopyBuf = new StringBuffer();
            double massAdditionFromModifiedMasses = 0;
            String parentPeptideSequence = getParentPeptideEntry().getPeptideSequence();

            if (mModifications == null)
            {
                mModifiedMass = getParentPeptideEntry().getMass();
                return;
            }

            for (int i=0; i<mModifications.length; i++)
            {
                if (mModifications[i] != null && mModifications[i].size() > 0)
                {
                    for (MS2Modification mod : mModifications[i])
                        massAdditionFromModifiedMasses += mod.getMassDiff();
                }
                else
                    sequenceCopyBuf.append(parentPeptideSequence.charAt(i));
            }
            Peptide fakePeptideForUnmodifiedMasses =
                    AmtPeptideEntry.getPeptideForPeptideSequence(parentPeptideSequence);
            mModifiedMass =
                    fakePeptideForUnmodifiedMasses.getMass(PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES) +
                                                           massAdditionFromModifiedMasses;
        }


        /**
         * Add a bunch of observations for this peptide, without recalculating stats,
         * then recalculate once.  Because recalculation is expensive
         *
         * @param observations
         */
        public void addObservations(AmtPeptideObservation[] observations)
        {
            for (AmtPeptideObservation observation : observations)
                addObservationNoRecalc(observation);
            recalculateStats();
        }

        /**
         * Add an observation for this peptide and recalculate the stats
         *
         * @param observation
         */
        public void addObservation(AmtPeptideObservation observation)
        {
            addObservationNoRecalc(observation);
            recalculateStats();
        }


        /**
         * Add an observation for this peptide.  Don't recalculate stats
         *
         * @param observation
         */
        protected void addObservationNoRecalc(AmtPeptideObservation observation)
        {
            if (mObservations == null)
                mObservations = new ArrayList<AmtPeptideObservation>();
            mObservations.add(observation);
        }

        /**
         * Remove an observation.  To do this you have to have the actual AmtPeptideObservation
         * object from this entry.  Might want to add a way to do this based on observed
         * value, too, but if so, need to figure out how close it needs to be.
         *
         * @param observation
         */
        public boolean removeObservation(AmtPeptideObservation observation)
        {
            if (mObservations == null)
                return false;

            boolean result = mObservations.remove(observation);
            recalculateStats();
            return result;
        }

        public void recalculateStats()
        {
            recalculateStats(AmtDatabase.DEFAULT_PRECISION);
        }


        /**
         * Recalculate the mean observed hydrophobicity, munging together all
         * observations, giving them equal weight.
         * Possibly, in the future, do other things.
         * Possibly, in the future, if we store probability on observations, weight them
         */
        public void recalculateStats(int precision)
        {
//        double sumObservedHydrophobicity = 0.0;
            double sumPeptideProphet = 0.0;
            for (AmtPeptideObservation obs : mObservations)
            {
//            sumObservedHydrophobicity += obs.getObservedHydrophobicity();
                sumPeptideProphet += obs.getPeptideProphet();
            }
            double[] observedHydrophobicities = new double[mObservations.size()];
            for (int i = 0; i < mObservations.size(); i++)
                observedHydrophobicities[i] = mObservations.get(i).getObservedHydrophobicity();
            double medianHydrophobicity = BasicStatistics.median(observedHydrophobicities);
            setMedianObservedHydrophobicity(Rounder.round(medianHydrophobicity, precision));
            setMedianPeptideProphet(Rounder.round(sumPeptideProphet / (double) mObservations.size(), precision));

            mHydrophobicityStandardDeviation =
                    Rounder.round(BasicStatistics.standardDeviation(observedHydrophobicities), precision);
        }



        /**
         * how many observations for this peptide entry?
         * @return
         */
        public int getNumObservations()
        {
            return mObservations.size();
        }

        public void setModifications(List<MS2Modification>[] modifications)
        {
            this.mModifications = modifications;
        }


        //simple getters and setters
        public String getModifiedSequence()
        {
            return mModifiedSequence;
        }

        public void setModifiedSequence(String mModifiedSequence)
        {
            this.mModifiedSequence = mModifiedSequence;
        }

        public AmtPeptideEntry getParentPeptideEntry()
        {
            return mParentPeptideEntry;
        }

        public void setParentPeptideEntry(AmtPeptideEntry mParentPeptideEntry)
        {
            this.mParentPeptideEntry = mParentPeptideEntry;
        }

        public List<MS2Modification>[] getModifications()
        {
            return mModifications;
        }

        public double getHydrophobicityStandardDeviation()
        {
            return mHydrophobicityStandardDeviation;
        }

        public void setHydrophobicityStandardDeviation(double newValue)
        {
            mHydrophobicityStandardDeviation = newValue;
        }

        public AmtPeptideObservation[] getObservations()
        {
            return mObservations.toArray(new AmtPeptideObservation[0]);
        }



        public double getMedianObservedHydrophobicity()
        {
            return mMedianObservedHydrophobicity;
        }

        public void setMedianObservedHydrophobicity(double meanObservedHydrophobicity)
        {
            this.mMedianObservedHydrophobicity = meanObservedHydrophobicity;
        }

        public double getMedianPeptideProphet()
        {
            return mMedianPeptideProphet;
        }

        public void setMedianPeptideProphet(double meanPeptideProphet)
        {
            this.mMedianPeptideProphet = meanPeptideProphet;
        }

        public double getModifiedMass()
        {
            return mModifiedMass;
        }

        public void setModifiedMass(double modifiedMass)
        {
            this.mModifiedMass = modifiedMass;
        }
    }

    /**
     * Represents a single peptide observation.  Carries a reference to the AmtRunEntry that
     * the observation came from
     */
    public static class AmtPeptideObservation
    {
        //a single hydrophobicity observation
        protected double mObservedHydrophobicity = 0.0;

        protected double mPeptideProphet = 0.0;

        protected AmtRunEntry runEntry = null;

        protected double mTimeInRun = 0;


        //Number of times the peptide was seen in this run.
        protected int mSpectralCount = SPECTRAL_COUNT_UNKNOWN;

        public static final int SPECTRAL_COUNT_UNKNOWN = -1;

        public AmtPeptideObservation()
        {
        }

        public String toString()
        {
            return "AmtPeptideObservation: observedHydrophobicity=" + getObservedHydrophobicity();
        }



        /**
         * Given an observed hydrophobicity, create an observation
         *
         * @param observedHydrophobicity
         * @return
         */
        public static AmtPeptideObservation createObservation(double observedHydrophobicity,
                                                              double qualityScore,
                                                              AmtRunEntry runEntry,
                                                              double timeInRun)
        {
            AmtPeptideObservation newObservation = new AmtPeptideObservation();
            newObservation.setObservedHydrophobicity(Rounder.round(observedHydrophobicity, AmtDatabase.DEFAULT_PRECISION));
            newObservation.setPeptideProphet(Rounder.round(qualityScore, AmtDatabase.DEFAULT_PRECISION));
            newObservation.setRunEntry(runEntry);
            newObservation.setTimeInRun(timeInRun);
            return newObservation;
        }

        //getters and setters

        public double getObservedHydrophobicity()
        {
            return mObservedHydrophobicity;
        }

        public void setObservedHydrophobicity(double observedHydrophobicity)
        {
            this.mObservedHydrophobicity = observedHydrophobicity;
        }


        public double getPeptideProphet()
        {
            return mPeptideProphet;
        }

        public void setPeptideProphet(double peptideProphet)
        {
            this.mPeptideProphet = peptideProphet;
        }

        public AmtRunEntry getRunEntry()
        {
            return runEntry;
        }

        public void setRunEntry(AmtRunEntry runEntry)
        {
            this.runEntry = runEntry;
        }

        public int getSpectralCount()
        {
            return mSpectralCount;
        }

        public void setSpectralCount(int spectralCount)
        {
            this.mSpectralCount = spectralCount;
        }

        public boolean hasSpectralCount()
        {
            return getSpectralCount() != SPECTRAL_COUNT_UNKNOWN;
        }


        public double getTimeInRun()
        {
            return mTimeInRun;
        }

        public void setTimeInRun(double timeInRun)
        {
            mTimeInRun = timeInRun;
        }
        
    }

}
