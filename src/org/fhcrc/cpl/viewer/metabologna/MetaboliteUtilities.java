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

public class MetaboliteUtilities
{
    public static final double DEFAULT_MEAN_UNIT_MASS_DEFECT = 0.00045;
    

    public static double calcTotalMassDefectOffset(double mass)
    {
        return calcTotalMassDefectOffset(mass, DEFAULT_MEAN_UNIT_MASS_DEFECT);
    }

    public static double calcTotalMassDefectOffset(double mass, double meanUnitMassDefect)
    {
        double totalMassOffset = mass % (1d + meanUnitMassDefect);
        if (totalMassOffset > 0.5)
            totalMassOffset -= 1d;
        return totalMassOffset;
    }

    public static double calcUnitMassDefectOffset(double mass)
    {
        return calcUnitMassDefectOffset(mass, DEFAULT_MEAN_UNIT_MASS_DEFECT);
    }

    public static double calcUnitMassDefectOffset(double mass, double meanUnitMassDefect)
    {
        return calcTotalMassDefectOffset(mass, meanUnitMassDefect) / mass;
    }
}
