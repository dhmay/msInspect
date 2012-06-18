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

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;


/**
 * This class knows how to adjust normalized hydrophobicities based on the effect of acrylamide,
 * vs. iodoacetamide, on peptide hydrophobicity.
 *
 * These techniques were developed based on a fractionated run of CSF by Wendy Law, containing
 * a pool of unlabeled CSF with iodoacetamide protecting the Cysteines, and the same CSF with
 * heavy-labeled acrylamide.
 *
 * Executive summary:  acrylamide increases hydrophobicity slightly over iodoacetamide.  More
 * Cysteines -> bigger effect
 *
 * WARNING:  there is actually a rather more complex relationship between peptide composition
 * and the acrylamide hydrophobicity contribution than we capture here.  It looks like there's
 * a nonlinear relationship between peptide hydrophobicity and the amount of contribution, and of
 * course that relationship is different for different numbers of Cysteines.  What you see here
 * is just a good-enough shorthand.
 *
 */
public class AcrylamideHydrophobicityAdjuster
{
    private static Logger _log = Logger.getLogger(AcrylamideHydrophobicityAdjuster.class);

    //these values are the median values that we observed in our experiment
    public static final double[] ACRYLAMIDE_HYDRO_DIFF_COEFFS_0_RESIDUES =
            new double[] { 0 , 0 };
    public static final double[] ACRYLAMIDE_HYDRO_DIFF_COEFFS_1_RESIDUE =
            new double[] { 0.00825, -0.0060 };
    public static final double[] ACRYLAMIDE_HYDRO_DIFF_COEFFS_2_RESIDUES =
            new double[] { 0.0155, -0.0120 };
    public static final double[] ACRYLAMIDE_HYDRO_DIFF_COEFFS_3PLUS_RESIDUES =
            new double[] { 0.0268, -0.0055 };

    public static double calculateAcrylamideHDifference(String peptideSequence, double oldH)
    {
        int numCysteines = 0;
        for (int i=0; i<peptideSequence.length(); i++)
        {
            if (peptideSequence.charAt(i) == 'C')
                numCysteines++;
        }

        double[] coefficientsThisPeptide;
        switch (numCysteines)
        {
            case 0:
                coefficientsThisPeptide = ACRYLAMIDE_HYDRO_DIFF_COEFFS_0_RESIDUES;
                break;
            case 1:
                coefficientsThisPeptide = ACRYLAMIDE_HYDRO_DIFF_COEFFS_1_RESIDUE;
                break;
            case 2:
                coefficientsThisPeptide = ACRYLAMIDE_HYDRO_DIFF_COEFFS_2_RESIDUES;
                break;                
            default:
                coefficientsThisPeptide = ACRYLAMIDE_HYDRO_DIFF_COEFFS_3PLUS_RESIDUES;
        }

        double deltaH =
                RegressionUtilities.predictYFromX(coefficientsThisPeptide[1], coefficientsThisPeptide[0], oldH);

        //Here's where I don't really trust our coefficients.  I don't think acrylamide ever leads to
        //a negative H contribution
        if (deltaH < 0)
            deltaH = 0;
        return deltaH;
    }

    public static double adjustHForAcrylamide(String peptideSequence, double oldH)
    {
        return oldH + calculateAcrylamideHDifference(peptideSequence, oldH);
    }

    public static double removeAcrylamideEffect(String peptideSequence, double oldH)
    {
        return oldH - calculateAcrylamideHDifference(peptideSequence, oldH);
    }
}
