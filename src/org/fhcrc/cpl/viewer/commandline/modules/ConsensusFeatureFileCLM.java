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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.viewer.align.*;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule for feature finding
 */
public class ConsensusFeatureFileCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ConsensusFeatureFileCLM.class);

    protected File[] featureFiles;
    protected int minRunsPerFeature = 2;
    protected File outFile;
    protected File pepArrayFile;

    protected boolean shouldRequireSamePeptide = false;

    protected boolean showCharts = false;
    protected float maxStudRes = 1;
    protected float maxLeverageNumerator=8;

    protected static final int ALIGNMENT_MODE_SPLINE = 0;
    protected static final int ALIGNMENT_MODE_QUANTILE = 1;

    protected int quantilePolynomialDegree = AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;
    protected int degreesOfFreedom = SplineAligner.DEFAULT_DEGREES_OF_FREEDOM;

    protected int intensityMode = PeptideArrayAnalyzer.CONSENSUS_INTENSITY_MODE_MEDIAN;

    protected boolean shouldAllowNegScanAndTime = false;


    protected final static String[] alignmentModeStrings = {
            "spline",
            "quantile"
    };

    protected final static String[] alignmentModeExplanations = {
            "Use the spline regression algorithm for alignment.  This is the original algorithm implemented. " +
            "It runs a bit faster than the quantile regression and behaves less oddly at the very extremes of the data.",
            "Use quantile regression algorithm for alignment.  This algorithm performs better with very noisy data."
    };

    protected int alignmentMode = ALIGNMENT_MODE_SPLINE;

    public int outFormat = FilterFeaturesCommandLineModule.OUT_FORMAT_MSINSPECT;


       public ConsensusFeatureFileCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "consensusfeaturefile";

        mShortDescription = "Create a 'consensus' feature file representing all of the features in multiple MS runs.";
        mHelpMessage = "Create a 'consensus' feature file representing all of the features in multiple MS runs. This is " +
                "done by creating a peptide array (or specifying an existing peptide array) and then creating one feature" +
                "for each cell in the array.";

        //todo: stop ignoring outformat.  This requires redoing the way we do consensus --
        //I'm currently using the details file to build features.
        CommandLineArgumentDefinition[] argDefs =
               {
                       new FileToWriteArgumentDefinition("out",false, "output file"),
                       createUnnamedSeriesFileArgumentDefinition(false, "input feature files"),
                       new FileToReadArgumentDefinition("peparray", false, "input peptide array"),
                       new IntegerArgumentDefinition("minfeatureruns", false,
                               "minimun number of runs a feature must appear in",
                               minRunsPerFeature),
                       new BooleanArgumentDefinition("samepeptide", false,
                               "Require all features in row to have non-conflicting peptide assignment??", shouldRequireSamePeptide),
                       new IntegerArgumentDefinition("topn", false, "number of features to use for alignment"),
                       new BooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                       new DecimalArgumentDefinition("maxstudres", false, "max studentized residual", maxStudRes),
                       new DecimalArgumentDefinition("maxleverage", false, "max leverage numerator", maxLeverageNumerator),
                       new EnumeratedValuesArgumentDefinition("alignmentmode",false,alignmentModeStrings,
                               alignmentModeExplanations, "spline"),
                       new EnumeratedValuesArgumentDefinition("intensitymode",false,PeptideArrayAnalyzer.INTENSITY_MODE_STRINGS,
                               PeptideArrayAnalyzer.INTENSITY_MODE_EXPLANATIONS, "median"),
                       new EnumeratedValuesArgumentDefinition("outformat",false,FilterFeaturesCommandLineModule.outFormatStrings,
                               FilterFeaturesCommandLineModule.outFormatExplanations),
                       new BooleanArgumentDefinition("allownegscanandtime", false, "Allow scan and time values to be " +
                               "negative? (otherwise, peg to 1 and 0, respectively",
                               shouldAllowNegScanAndTime),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = getUnnamedSeriesFileArgumentValues();
        pepArrayFile = getFileArgumentValue("peparray");
        if ((featureFiles == null) == (pepArrayFile == null))
            throw new ArgumentValidationException("Please specify feature files, or peptide array file, but not both");

        if (hasArgumentValue("outformat"))
        {
            outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(getStringArgumentValue("outformat"));
        }


        outFile = getFileArgumentValue("out");
        minRunsPerFeature = getIntegerArgumentValue("minfeatureruns");
        shouldRequireSamePeptide = getBooleanArgumentValue("samepeptide");
        showCharts = getBooleanArgumentValue("showcharts");
        maxStudRes = getFloatArgumentValue("maxstudres");
        maxLeverageNumerator = getFloatArgumentValue("maxleverage");
        alignmentMode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("alignmentmode")).getIndexForArgumentValue(
                getStringArgumentValue("alignmentmode"));
        intensityMode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("intensitymode")).getIndexForArgumentValue(
                getStringArgumentValue("intensitymode"));
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (pepArrayFile == null)
        {
            pepArrayFile = new File(outFile.getParent(), "temp.array.tsv");
//        File pepArrayDetailsFile = new File(outFile.getParent(), "temp.array.details.tsv");


            List<File> featureFileList = new ArrayList<File>(featureFiles.length);
            for (File featureFile : featureFiles)
                featureFileList.add(featureFile);
            BucketedPeptideArray arr = new BucketedPeptideArray(featureFileList, new FeatureSet.FeatureSelector());


            Aligner aligner = null;
            switch (alignmentMode)
            {
                case ALIGNMENT_MODE_SPLINE:
                    SplineAligner splineAligner = new SplineAligner();
                    splineAligner.setDegreesOfFreedom(degreesOfFreedom);
                    aligner = splineAligner;
                    break;
                case ALIGNMENT_MODE_QUANTILE:
                    QuantileRegressionAligner qrAligner = new QuantileRegressionAligner();
                    qrAligner.setNonlinearMappingPolynomialDegree(quantilePolynomialDegree);
                    aligner = qrAligner;
                    break;
                default:
                    throw new CommandLineModuleExecutionException("Unknown mode");

            }
            arr.setAligner(aligner);

            if (hasArgumentValue("topn"))
            {
                Aligner.FeaturePairSelector featurePairSelector =
                        Aligner.DEFAULT_FEATURE_PAIR_SELECTOR;
                arr.setFeaturePairSelector(featurePairSelector);
                ((Aligner.MassOrMzFeaturePairSelector) featurePairSelector).setTopN(getIntegerArgumentValue("topn"));
            }



            arr.setOutFileName(pepArrayFile.getAbsolutePath());
            arr.setAlign(true);
            arr.getAligner().setBuildCharts(showCharts);
            arr.getAligner().setMaxStudRes(maxStudRes);
            arr.getAligner().setMaxLeverageNumerator(maxLeverageNumerator);
            arr.run(true, FeatureGrouper.BUCKET_EVALUATION_MODE_ONE_FROM_EACH, showCharts);
        }
//        arr.getAligner().getWarpingMaps()

        try
        {
            PeptideArrayAnalyzer peptideArrayAnalyzer =
                    new PeptideArrayAnalyzer(pepArrayFile);
            File detailsFile = new File(BucketedPeptideArray.calcDetailsFilepath(pepArrayFile.getAbsolutePath()));
            ApplicationContext.infoMessage("Using details file " + detailsFile.getAbsolutePath());
            FeatureSet consensusFeatureSet =
                    peptideArrayAnalyzer.createConsensusFeatureSet(detailsFile,
                            minRunsPerFeature, intensityMode,
                            shouldRequireSamePeptide, shouldAllowNegScanAndTime);
            ApplicationContext.infoMessage("Created consensus feature set with " +
                    consensusFeatureSet.getFeatures().length + " features");
            consensusFeatureSet.save(outFile);
            ApplicationContext.infoMessage("Saved consensus feature set to file " +
                    outFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }



    }



}
