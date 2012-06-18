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

package org.fhcrc.cpl.toolbox.chem;

import java.util.List;
import java.util.ArrayList;

/**
 *
 * Describes a chemical element
 *
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 25, 2010
 * Time: 11:56:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class Element {
    protected String symbol;
    protected double commonestIsotopeMass;
    protected double averageMass;
    protected int atomicNumber;
    protected int nominalMass;

    protected double[] isotopicPeakMassesWithoutMissing;
    protected double[] isotopicPeakFrequenciesWithoutMissing;

    protected double[] isotopicPeakMassesWithMissing;
    protected double[] isotopicPeakFrequenciesWithMissing;


    public String getSymbol() {
        return symbol;
    }

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    public int getAtomicNumber() {
        return atomicNumber;
    }

    public void setAtomicNumber(int atomicNumber) {
        this.atomicNumber = atomicNumber;
    }

    public Element(String symbol, double commonestIsotopeMass, double averageMass, int atomicNumber)
    {
        this.symbol = symbol;
        this.commonestIsotopeMass = commonestIsotopeMass;
        this.averageMass = averageMass;
        this.atomicNumber = atomicNumber;
        nominalMass = (int) Math.round(commonestIsotopeMass);
    }

    /**
     * Constructor that takes isotope masses and frequencies and calculates commonest isotope mass
     * and average mass.
     * isotope masses should /not/ include 0-frequency (missing) peaks.
     * @param symbol
     * @param atomicNumber
     * @param isotopeMasses
     * @param isotopeFrequencies
     */
    public Element(String symbol, int atomicNumber, double[] isotopeMasses, double[] isotopeFrequencies)
    {
        this.symbol = symbol;
        this.atomicNumber = atomicNumber;
        setIsotopicPeakMassesAndFrequencies(isotopeMasses, isotopeFrequencies);

        double massesSum = 0;
        double weightsSum = 0;
        double highestFrequency = 0;
        for (int i=0; i<isotopeMasses.length; i++)
        {
            massesSum += isotopeMasses[i] * isotopeFrequencies[i];
            weightsSum += isotopeFrequencies[i];

            if (isotopeFrequencies[i] > highestFrequency)
            {
                commonestIsotopeMass = isotopeMasses[i];
                highestFrequency = isotopeFrequencies[i];
            }
        }

        this.averageMass = massesSum / weightsSum;
        nominalMass = (int) Math.round(commonestIsotopeMass);
    }

    public double getCommonestIsotopeMass() {
        return commonestIsotopeMass;
    }

    public void setCommonestIsotopeMass(double commonestIsotopeMass) {
        this.commonestIsotopeMass = commonestIsotopeMass;
    }

    public double getAverageMass() {
        return averageMass;
    }

    public void setAverageMass(double averageMass) {
        this.averageMass = averageMass;
    }



    /**
     * Created arrays of masses and frequencies including non-occurring masses with frequency 0
     * @param isotopicPeakMasses
     * @param isotopicPeakFrequencies
     */
    public void setIsotopicPeakMassesAndFrequencies(double[] isotopicPeakMasses, double[] isotopicPeakFrequencies)
    {
        this.isotopicPeakMassesWithoutMissing = isotopicPeakMasses;
        this.isotopicPeakFrequenciesWithoutMissing = isotopicPeakFrequencies;

        List<Double> isotopicPeakMassesPadded = new ArrayList<Double>();
        List<Double> isotopicPeakFrequenciesPadded = new ArrayList<Double>();

        double previousMass = 0;
        for (int i=0; i<isotopicPeakMasses.length; i++)
        {
            double consecutivePeakMassDiff = isotopicPeakMasses[i] - previousMass;
            if (i>0 && consecutivePeakMassDiff > 1.3)
            {
                int numMissingPeaks = (int) Math.round(consecutivePeakMassDiff / 1.0078) - 1;
                for (int j=1; j <= numMissingPeaks; j++)
                {
                    isotopicPeakMassesPadded.add(previousMass + (j * (consecutivePeakMassDiff / (double) (numMissingPeaks+1))));
                    isotopicPeakFrequenciesPadded.add(0.0);
                }
            }
            isotopicPeakMassesPadded.add(isotopicPeakMasses[i]);
            isotopicPeakFrequenciesPadded.add(isotopicPeakFrequencies[i]);
            previousMass = isotopicPeakMasses[i];
        }
        double[] freqsPadded = new double[isotopicPeakFrequenciesPadded.size()];
        double[] massesPadded = new double[isotopicPeakFrequenciesPadded.size()];
        for (int i=0; i<freqsPadded.length; i++)
        {
            freqsPadded[i] = isotopicPeakFrequenciesPadded.get(i);
            massesPadded[i] = isotopicPeakMassesPadded.get(i);
        }

        setIsotopicPeakMassesWithMissing(massesPadded);
        setIsotopicPeakFrequenciesWithMissing(freqsPadded);
    }

    /**
     * IMPORTANT!  These isotopic peak masses do NOT include 0-frequency peaks!  If you
     * want those, use getIsotopicPeakMassesPadded()
     * @return
     */
    public double[] getIsotopicPeakMassesWithoutMissing() {
        return isotopicPeakMassesWithoutMissing;
    }


    /**
     * IMPORTANT!  These isotopic peak frequencies do NOT include 0-frequency peaks!  If you
     * want those, use getIsotopicPeakFrequenciesPadded()
     * @return
     */
    public double[] getIsotopicPeakFrequenciesWithoutMissing() {
        return isotopicPeakFrequenciesWithoutMissing;
    }

    public void setIsotopicPeakFrequenciesWithoutMissing(double[] isotopicPeakFrequencies)
    {
        this.isotopicPeakFrequenciesWithoutMissing = isotopicPeakFrequencies;
    }

    public double[] getIsotopicPeakMassesWithMissing() {
        return isotopicPeakMassesWithMissing;
    }


    /**
     * Includes 0-frequency masses
     * @param isotopicPeakMassesWithMissing
     */
    public void setIsotopicPeakMassesWithMissing(double[] isotopicPeakMassesWithMissing) {
        this.isotopicPeakMassesWithMissing = isotopicPeakMassesWithMissing;
    }

    /**
     * Includes 0-frequency peaks
     * @return
     */
    public double[] getIsotopicPeakFrequenciesWithMissing() {
        return isotopicPeakFrequenciesWithMissing;
    }

    public void setIsotopicPeakFrequenciesWithMissing(double[] isotopicPeakFrequenciesWithMissing) {
        this.isotopicPeakFrequenciesWithMissing = isotopicPeakFrequenciesWithMissing;
    }

    public String toString()
    {
        StringBuffer buf = new StringBuffer("Element: " + symbol + ", anumber=" + atomicNumber +
                ", monomass=" + getCommonestIsotopeMass() + ", avemass=" + getAverageMass());
        if (isotopicPeakFrequenciesWithMissing != null)
        {
            buf.append(" all masses (frequencies): ");
            for (int i=0; i< isotopicPeakFrequenciesWithMissing.length; i++)
            {
                if (i>0)
                    buf.append(", ");
               buf.append(isotopicPeakMassesWithMissing[i] + " (" + isotopicPeakFrequenciesWithMissing[i] + ")");
            }
        }

        return buf.toString();
    }

    public int getNominalMass() {
        return nominalMass;
    }
}
