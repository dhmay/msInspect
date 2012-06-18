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
import org.fhcrc.cpl.toolbox.datastructure.Pair;

import java.util.Map;
import java.util.HashMap;

/**
 * Converts hydrophobicity values from various algorithms to normalized
 * values.  Stores information about all known hydrophobicity algorithms.
 *
 * Normalized values have a mean of 0 and a standard deviation of 1 when
 * calculated for the entire IPI human database.
 *
 * Normalization of an algorithm is performed by calculating the mean and
 * standard deviation of that algorithm on the all peptides in theentire
 * IPI human database.  Minimum peptide length 6, one missed cleavage allowed.
 * 
 * Individual values are normalized by scaling and shifting the value
 * appropriately.
 *
 * The IPI human database used for these calculations is the human IPI
 * database released on July 13, 2006, filename ipi.HUMAN.fasta.20060713
 *
 */
public class HydrophobicityNormalizer
{
    private static Logger _log = Logger.getLogger(HydrophobicityNormalizer.class);

    //constants for easy reference

    //Krokhin's algorithm
    public static final double KROKHIN_1_MEAN=30.763700485229492;
    public static final double KROKHIN_1_STDDEV=21.886646343829856;

    public static final double KROKHIN_3_MEAN=29.575716018676758;
    public static final double KROKHIN_3_STDDEV=17.596638381015442;

//redone with no missed cleavages
//public static final double KROKHIN_3_MEAN=27.470834732055664;
//public static final double KROKHIN_3_STDDEV=17.091617040482113;


    //predictNET
    public static final double PREDICTNET_1_MEAN=0.337;
    public static final double PREDICTNET_1_STDDEV=0.210;


    //Key is based on algorithm name and version, as built by constructKey().
    //Value is a pair containing mean and standard deviation, in that order.
    protected static Map<String, Pair<Double,Double>> _algStatsMap = null;



    static
    {
        _algStatsMap = new HashMap<String, Pair<Double,Double>>();

        addAlgorithmStatistics("krokhin", 1.0,
                               KROKHIN_1_MEAN,KROKHIN_1_STDDEV);
        addAlgorithmStatistics("krokhin", 3.0,
                               KROKHIN_3_MEAN,KROKHIN_3_STDDEV);
        addAlgorithmStatistics("predictnet", 1.0,
                               PREDICTNET_1_MEAN,PREDICTNET_1_STDDEV);
    }

    /**
     * For adding stats on new algorithms.  Public, so it can be added by any
     * standard or custom piece of code
     * @param algorithmName
     * @param algorithmVersion
     * @param mean
     * @param stddev
     */
    public static void addAlgorithmStatistics(String algorithmName,
                                              double algorithmVersion,
                                              double mean, double stddev)
    {
        _algStatsMap.put(constructKey(algorithmName, algorithmVersion),
                         new Pair<Double,Double>(mean, stddev));
    }

    protected static String constructKey(String algorithmName,
                                         double algorithmVersion)
    {
        return algorithmName.toLowerCase() + "_" +
              algorithmVersion;
    }

    public static boolean algorithmVersionKnown(String algorithmName,
                                                double algorithmVersion)
    {
        return _algStatsMap.containsKey(constructKey(algorithmName, algorithmVersion));
    }

    /**
     * Normalize a hydrophobicity value calculated using the specified
     * algorithm name and version: subtract the mean, divide by the
     * standard deviation.
     *
     * Will throw NullPointerException if algorithm info isn't there
     * @param inputHydrophobicity
     * @param algorithmName
     * @param algorithmVersion
     * @return
     */
    public static double normalize(double inputHydrophobicity,
                                   String algorithmName,
                                   double algorithmVersion)
    {
        double result;
        try
        {
            Pair<Double,Double> algInfo =
                    _algStatsMap.get(constructKey(algorithmName, algorithmVersion));
            //shift by the mean, to center on 0
            result = inputHydrophobicity - algInfo.first;
            //divide by the standard deviation, so that the standard deviation of
            //normalized values will be 1
            result = result / algInfo.second;
        }
        catch (NullPointerException npe)
        {
            _log.error("normalize: Failed to find algorithm information for algorithm " +
                       algorithmName + ", version " + algorithmVersion);
            throw npe;
        }
        return result;
    }

}
