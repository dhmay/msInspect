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
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;


/**
 * Utilities related to peptide mass
 */
public class MassUtilities
{
    private static Logger _log = Logger.getLogger(MassUtilities.class);


    public static float calcMassForMzAndCharge(float mz, int charge)
    {
        return (mz - Spectrum.HYDROGEN_ION_MASS) * charge;   
    }

    public static float calcMzForMassAndCharge(float mass, int charge)
    {
        return (mass / charge) + Spectrum.HYDROGEN_ION_MASS;
    }

        /**
         * Utility method to calculate the absolute mass tolerance, given a mass tolerance
         * parameter that may be absolute or relative
         * @param centerMass
         * @param deltaMass
         * @param deltaMassType
         * @return
         */
    public static float calculateAbsoluteDeltaMass(float centerMass,
                                                   float deltaMass,
                                                   int deltaMassType)
    {
        if (deltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE)
            return deltaMass;
        //deltamass must be in ppm
        return convertPPMToDa(deltaMass, centerMass);
    }

    /**
     *
     * @param deltaMass
     * @param centerMass
     * @return
     */
    public static float convertPPMToDa(float deltaMass, float centerMass)
    {
        return (deltaMass * centerMass) / 1000000;
    }

    /**
     * Utility method to calculate the PPM mass tolerance, given a mass tolerance
     * parameter that may be absolute or relative
     * @param centerMass
     * @param deltaMass
     * @param deltaMassType
     * @return
     */
    public static float calculatePPMDeltaMass(float centerMass,
                                              float deltaMass,
                                              int deltaMassType)
    {
        if (deltaMassType == FeatureSetMatcher.DELTA_MASS_TYPE_PPM)
            return deltaMass;
        //deltamass must be in ppm
        return convertDaToPPM(deltaMass, centerMass);
    }

    /**
     * Convert Daltons to PPM
     * @param centerMass
     * @param deltaMass
     * @return
     */
    public static float convertDaToPPM(float deltaMass , float centerMass)
    {
        return deltaMass * 1000000 / centerMass;
    }

    /**
     * Calculate the neutral mass of a given peptide sequence with no modifications
     * @param peptideSequence
     * @return
     */
    public static float calcPeptideNeutralMass(String peptideSequence)
    {
        //this is lame.  Really peptide should have a constructor that doesn't require a protein
        Protein dummyProtein = new Protein("",peptideSequence.getBytes());
        Peptide dummyPeptide = new Peptide(dummyProtein,0,dummyProtein.getBytes().length);

        return (float) dummyPeptide.getMass(PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES);
    }

    /**
     * Calculate the neutral mass of a given peptide sequence with the specified STATIC modifications.
     *
     * @param peptideSequence
     * @param modifications
     */
    public static float calcModifiedPeptideNeutralMass(String peptideSequence, MS2Modification[] modifications)
    {
        float unmodifiedMass = calcPeptideNeutralMass(peptideSequence);
        if (modifications == null || modifications.length == 0)
        {
            return unmodifiedMass;
        }

        float massAdditionFromModifiedMasses = 0;

        for (int i=0; i<peptideSequence.length(); i++)
        {
            for (MS2Modification mod : modifications)
            {
                if (!mod.getVariable() && mod.getAminoAcid().charAt(0) == peptideSequence.charAt(i))
                    massAdditionFromModifiedMasses += mod.getMassDiff();
            }
        }
        return unmodifiedMass + massAdditionFromModifiedMasses;
    }
}
