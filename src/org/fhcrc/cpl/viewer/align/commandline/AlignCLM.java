/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.align.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.arguments.*;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.viewer.align.QuantileRegressionAligner;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.apache.log4j.Logger;


import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class AlignCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AlignCLM.class);


    protected FeatureSet[] featureSets;

    protected float massMatchDeltaMass = 10;
    protected int massMatchDeltaMassType = AmtFeatureSetMatcher.DELTA_MASS_TYPE_PPM;

    protected int nonlinearMappingPolynomialDegree =
            AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;


    protected double maxStudRes = AmtDatabaseMatcher.DEFAULT_MAX_STUDENTIZED_RESIDUAL;

    protected double maxLeverageNumerator = AmtDatabaseMatcher.DEFAULT_LEVERAGE_NUMERATOR;

    protected int topN = -1;

    protected int mode = -1;

    protected static final int MODE_SPLINE = 0;
    protected static final int MODE_QUANTILE = 1;

    protected final static String[] modeStrings = {
            "spline",
            "quantile"
    };

    protected final static String[] modeExplanations = {
            "spline",
            "quantile"
    };

    public AlignCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "align";
        mShortDescription = "Align feature sets";
        mHelpMessage = "Align feature sets";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createEnumeratedArgumentDefinition("mode",true,modeStrings,
                                modeExplanations),
                        createUnnamedSeriesArgumentDefinition(ViewerArgumentDefinitionFactory.FEATURE_FILE,true, null),
                        createIntegerArgumentDefinition("mappingpolynomialdegree", false,
                                "The degree of the polynomial to fit when mapping time to hydrophobicity nonlinearly",
                                nonlinearMappingPolynomialDegree),
                        createDecimalArgumentDefinition("maxstudres", false,
                                "Maximum studentized residual for regression", maxStudRes),
                        createDecimalArgumentDefinition("maxleverage", false,
                                "Maximum NUMERATOR of the leverage of features used for regression", maxLeverageNumerator),
                        createIntegerArgumentDefinition("topn", false,
                                "topN argument for spline-based mapping",
                                topN),
                        createDeltaMassArgumentDefinition("deltamass", false,
                                "delta-mass for matching features for alignment"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureSets = new FeatureSet[this.getUnnamedSeriesArgumentValues().length];

        if (featureSets.length < 2)
            throw new ArgumentValidationException("More files, please");
        for (int i=0; i<this.getUnnamedSeriesArgumentValues().length; i++)
            featureSets[i] = (FeatureSet) this.getUnnamedSeriesArgumentValues()[i];
        nonlinearMappingPolynomialDegree = getIntegerArgumentValue("mappingpolynomialdegree");
        maxStudRes = getDoubleArgumentValue("maxstudres");
        maxLeverageNumerator = getDoubleArgumentValue("maxleverage");

        if (hasArgumentValue("deltamass"))
        {
            DeltaMassArgumentDefinition.DeltaMassWithType deltaMassWithType =
                    getDeltaMassArgumentValue("deltamass");
            massMatchDeltaMass = deltaMassWithType.getDeltaMass();
            massMatchDeltaMassType = deltaMassWithType.getDeltaMassType();
        }

        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        topN = getIntegerArgumentValue("topn");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        FeatureSet baseFeatureSet = featureSets[0];

        for (int i=1; i<featureSets.length; i++)
        {
            FeatureSet featureSetToAlign = featureSets[i];

            Aligner aligner = null;
            switch (mode)
            {
                case MODE_SPLINE:
                    SplineAligner splineAligner = new SplineAligner();
                    aligner = splineAligner;
                    break;
                case MODE_QUANTILE:
                    QuantileRegressionAligner qrAligner = new QuantileRegressionAligner();
                    aligner = qrAligner;
                    break;
                default:
                    throw new CommandLineModuleExecutionException("Unknown mode");
                    
            }

            Aligner.MassOrMzFeaturePairSelector pairSelector = new Aligner.MassFeaturePairSelector();
            aligner.setFeaturePairSelector(pairSelector);
            pairSelector.setTopN(topN);
            pairSelector.setDeltaMassOrMz(massMatchDeltaMass);
            pairSelector.setDeltaMassType(massMatchDeltaMassType);
            aligner.setMaxLeverageNumerator(maxLeverageNumerator);
            aligner.setMaxStudRes(maxStudRes);

            aligner.setBuildCharts(true);

            List<FeatureSet> featureSetList = new ArrayList<FeatureSet>(2);
            featureSetList.add(baseFeatureSet);
            featureSetList.add(featureSetToAlign);
            aligner.alignFeatureSets(featureSetList, true);
        }
       
    }


}
