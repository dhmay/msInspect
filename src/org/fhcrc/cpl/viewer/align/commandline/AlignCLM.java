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
package org.fhcrc.cpl.viewer.align.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.viewer.align.QuantileRegressionAligner;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;


import java.util.*;
import java.io.File;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class AlignCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AlignCLM.class);


    protected float massMatchDeltaMass = 10;
    protected int massMatchDeltaMassType = FeatureSetMatcher.DELTA_MASS_TYPE_PPM;

    protected int nonlinearMappingPolynomialDegree =
            AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;

    protected File[] featureFiles;

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
                        new EnumeratedValuesArgumentDefinition("mode",true,modeStrings,
                                modeExplanations),
                        createUnnamedSeriesFileArgumentDefinition(true, null),
                        new IntegerArgumentDefinition("mappingpolynomialdegree", false,
                                "The degree of the polynomial to fit when mapping time to hydrophobicity nonlinearly",
                                nonlinearMappingPolynomialDegree),
                        new DecimalArgumentDefinition("maxstudres", false,
                                "Maximum studentized residual for regression", maxStudRes),
                        new DecimalArgumentDefinition("maxleverage", false,
                                "Maximum NUMERATOR of the leverage of features used for regression", maxLeverageNumerator),
                        new IntegerArgumentDefinition("topn", false,
                                "topN argument for spline-based mapping",
                                topN),
                        new DeltaMassArgumentDefinition("deltamass", false,
                                "delta-mass for matching features for alignment"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = getUnnamedSeriesFileArgumentValues();

        if (featureFiles.length < 2)
            throw new ArgumentValidationException("More files, please");
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
        FeatureSet baseFeatureSet = null;

        for (int i=0; i<featureFiles.length; i++)
        {
            File featureFile = featureFiles[i];
            FeatureSet featureSetToAlign = null;
            try
            {
                featureSetToAlign = new FeatureSet(featureFile);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException("Failed to load feature file " +
                        featureFile.getAbsolutePath());
            }
            if (i==0)
            {
                baseFeatureSet = featureSetToAlign;
                continue;
            }

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
