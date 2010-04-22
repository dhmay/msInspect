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
