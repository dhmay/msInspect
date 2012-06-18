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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.fhcrc.cpl.toolbox.proteomics.Clusterer2D;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * 2-dimensional clusterer for clustering Features
 */
public class FeatureClusterer extends Clusterer2D
{
    private static Logger _log = Logger.getLogger(FeatureClusterer.class);

    private List<FeatureSet> _featureSets = new ArrayList<FeatureSet>();

    //different modes for clustering by mass or mz
    public static final int MASS_MZ_MODE_MASS=0;
    public static final int MASS_MZ_MODE_MZ=1;
    public static final int DEFAULT_MASS_MZ_MODE=MASS_MZ_MODE_MZ;

    //constants used to specify how mass tolerance should be calculated
    public static final int DELTA_MASS_TYPE_ABSOLUTE = 0;
    public static final int DELTA_MASS_TYPE_PPM = 1;
    public static final int DEFAULT_DELTA_MASS_TYPE = DELTA_MASS_TYPE_ABSOLUTE;

    //different modes for clustering by elution time
    public static final int ELUTION_MODE_TIME=0;
    public static final int ELUTION_MODE_SCAN=1;
    public static final int DEFAULT_ELUTION_MODE=ELUTION_MODE_TIME;

    protected int _massMzMode = DEFAULT_MASS_MZ_MODE;
    protected int _elutionMode = DEFAULT_ELUTION_MODE;

    protected int _massType = DEFAULT_DELTA_MASS_TYPE;




    /**
     * Initialize the list of featuresets
     */
    public FeatureClusterer()
    {
        _featureSets = new ArrayList<FeatureSet>();
    }

    public FeatureClusterer(int massMzMode, int elutionMode)
    {
        this();
        _massMzMode = massMzMode;
        _elutionMode = elutionMode;
    }

    public FeatureClusterer(int massMzMode, int elutionMode, FeatureSet featureSet)
    {
        this();
        _massMzMode = massMzMode;
        _elutionMode = elutionMode;
        addSet(featureSet);
    }

    /**
     * Add a featureset to be clustered
     * @param featureSet
     */
    public void addSet(FeatureSet featureSet)
    {
        _featureSets.add(featureSet);

            FeatureClusterable[] featureClusterables =
                    createClusterablesForFeatures(featureSet.getFeatures());
            super.addSet(featureClusterables);
    }

    public FeatureSet getSet(int i)
    {
        return _featureSets.get(i);
    }

    /**
     * Given an array of features, create an array of Clusterables.  The Clusterables should
     * be of the appropriate type for the mass/mz and elution mode
     * @param features
     * @return
     */
    public FeatureClusterable[] createClusterablesForFeatures(Feature[] features)
    {
        FeatureClusterable[] result = new FeatureClusterable[features.length];
        for (int i = 0; i < features.length; i++)
        {
            switch (_massMzMode)
            {
                case MASS_MZ_MODE_MASS:
                    switch (_elutionMode)
                    {
                        case ELUTION_MODE_TIME:
                            result[i] = new FeatureMassTimeClusterable(features[i]);
                            break;
                        case ELUTION_MODE_SCAN:
                            result[i] = new FeatureMassScanClusterable(features[i]);
                            break;
                    }
                    break;
                case MASS_MZ_MODE_MZ:
                    switch (_elutionMode)
                    {
                        case ELUTION_MODE_TIME:
                            result[i] = new FeatureMzTimeClusterable(features[i]);
                            break;
                        case ELUTION_MODE_SCAN:
                            result[i] = new FeatureMzScanClusterable(features[i]);
                            break;
                    }
                    break;
            }
        }
        return result;
    }


    //getters and setters

    public int getElutionMode()
    {
        return _elutionMode;
    }

    public void setElutionMode(int _elutionMode)
    {
        this._elutionMode = _elutionMode;
        setDimension2IsInt(true);
    }

    public int getMassMzMode()
    {
        return _massMzMode;
    }

    public void setMassMzMode(int _massMzMode)
    {
        this._massMzMode = _massMzMode;
    }


    //Clusterable classes specific to clustering Features




    /**
     * Pragmatic abstract class to help you cluster sets of features.
     */
    public abstract static class FeatureClusterable implements Clusterer2D.Clusterable
    {
        public Feature parentFeature;

        public FeatureClusterable(Feature feature)
        {
            setParentFeature(feature);
        }

        public void setParentFeature(Feature feature)
        {
            parentFeature = feature;
        }

        public Feature getParentFeature()
        {
            return parentFeature;
        }
    }

    /**
     * Clusterable for clustering based on feature mass and scan.
     * Allows access to the original Feature
     */
    public static class FeatureMassScanClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMassScanClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mass;
        }
        public double getDimension2Value()
        {
            return parentFeature.getScan();
        }
    }

    /**
     * Clusterable for clustering based on feature mass and time.
     */
    public static class FeatureMassTimeClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMassTimeClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mass;
        }
        public double getDimension2Value()
        {
            return parentFeature.getTime();
        }
    }

    /**
     * Clusterable for clustering based on feature mz and scan.
     */
    public static class FeatureMzScanClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMzScanClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mz;
        }
        public double getDimension2Value()
        {
            return parentFeature.getScan();
        }
    }

    /**
     * Clusterable for clustering based on feature mz and time.
     */
    public static class FeatureMzTimeClusterable
        extends FeatureClusterer.FeatureClusterable
    {
        public FeatureMzTimeClusterable(Feature feature)
        {
            super(feature);
        }

        public double getDimension1Value()
        {
            return parentFeature.mz;
        }
        public double getDimension2Value()
        {
            return parentFeature.getTime();
        }
    }


    public int getMassType()
    {
        return _massType;
    }

    public void setMassType(int _massType)
    {
        if (_massMzMode == MASS_MZ_MODE_MZ)
            throw new RuntimeException("FeatureClusterer: Tried to set mass type when in m/z mode");

        this._massType = _massType;


        if (_massType == DELTA_MASS_TYPE_PPM)
        {
            setDimensionSplitCalculator(new Clusterer2D.ClusterDimensionSplitCalculator()
            {
                public double calculateDimension1ForSplit(double referenceValue, double dimensionValue)
                {
                    double dim1 = calculateAbsoluteDeltaMass((float) referenceValue,
                            (float) dimensionValue,
                            DELTA_MASS_TYPE_PPM);
                    return dim1;
                }

                public double calculateDimension2ForSplit(double referenceValue, double dimensionValue)
                {
                    return dimensionValue;
                }

                public String toString()
                {
                    return "THE RIGHT ONE";
                }
            });
        }
        else
        {
            setDimensionSplitCalculator(new DefaultClusterDimensionSplitCalculator());
        }           
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
        if (deltaMassType == DELTA_MASS_TYPE_ABSOLUTE)
            return deltaMass;
        //deltamass must be in ppm
        return (deltaMass * centerMass) / 1000000;
    }


}
