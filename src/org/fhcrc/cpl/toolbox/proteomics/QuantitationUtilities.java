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
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;


/**
 * Utilities for working with quantitation
 */
public class QuantitationUtilities
{
    protected static Logger _log = Logger.getLogger(QuantitationUtilities.class);

    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;
    public static final String ALGORITHM_Q3 = "Q3";
    public static final String ALGORITHM_XPRESS = "XPRESS";
    public static final int LABEL_ACRYLAMIDE = 0;
    public static final int LABEL_LYCINE = 1;
    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
    public static final float SILAC_LABEL_MASS = 134.115092f;
    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;


    /**
     * Calculate the difference between heavy and light masses of a completely heavy and completely light
     * labeled version of this peptide
     * @param peptide
     * @param labelType
     * @return
     */
    public static float calcHeavyLightMassDiff(String peptide, int labelType)
    {
        String residue = "";
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                residue = "C";
                break;
            case LABEL_LYCINE:
                residue = "K";
                break;
        }
        int numLabels = 0;
        for (int i=0; i<peptide.length(); i++)
        {
            if (residue.equals("" + peptide.charAt(i)))
                numLabels++;
        }

        return numLabels * (labelType == LABEL_ACRYLAMIDE ? ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE :
                SILAC_LABEL_MASSDIFF_PERRESIDUE);
    }
}
