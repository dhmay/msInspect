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

import org.fhcrc.cpl.viewer.commandline.modules.FeatureSelectionParamsCommandLineModule;
import org.fhcrc.cpl.viewer.align.BucketedPeptideArray;
import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.viewer.align.QuantileRegressionAligner;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureClusterer;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseMatcher;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;


/**
 * Commandline module for creating peptide arrays.
 * dhmay changing 20100104.  Big changes: separating features by charge before clustering.  Alignment is still done based
 * on a deconvoluted featureset, and normalization is done on a deconvoluted featureset in which all feature charges
 * are assigned 1.  But clustering is done separately for each charge state.  That way, feature intensity comparison
 * actually means something.
 *
 * In support of these changes, adding '--deconvolute' argument, so that the user can get behavior that's more
 * like the earlier functionality
 *
 * dhmay 20100310: changing default to alignment by retention time, and adding option for alignment by scan
 */
public class PeptideArrayCommandLineModule extends FeatureSelectionParamsCommandLineModule
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptideArrayCommandLineModule.class);

    protected File outFile = null;
    protected File tagFile = null;
    protected double elutionBucket=50;
    protected double massBucket=.2;

    protected int massType = FeatureClusterer.DELTA_MASS_TYPE_PPM;

    protected double[] elutionBuckets = null;
    protected double[] massBuckets = null;

    protected int topN = 0;

    protected boolean normalize = false;
    protected boolean align = true;
    protected boolean optimize = false;
    protected int optimizationMode = BucketedPeptideArray.DEFAULT_OPTIMIZATION_MODE;
    protected boolean optimizeOnPeptideIds = false;
    protected File[] featureFiles = null;

    //for spline-based alignment
    protected int degreesOfFreedom = SplineAligner.DEFAULT_DEGREES_OF_FREEDOM;

    protected static final String[] alignByTagsStrings = { "none", "loose", "strict", };
    protected String alignByTags = alignByTagsStrings[0];

    protected int peptideMatchScore = FeatureGrouper.DEFAULT_PEPTIDE_AGREEMENT_MATCH_SCORE;
    protected int peptideMismatchPenalty = FeatureGrouper.DEFAULT_PEPTIDE_AGREEMENT_MISMATCH_PENALTY;

    protected EnumeratedValuesArgumentDefinition featurePairSelectorArgDef =
                new EnumeratedValuesArgumentDefinition("featurepairselector", false,
                                "How features should be paired when performing alignment.  Default is MZ",
                                featurePairSelectorStrings);

    protected Aligner.FeaturePairSelector featurePairSelector =
            Aligner.DEFAULT_FEATURE_PAIR_SELECTOR;

    protected int alignmentOrderMode = Aligner.ALIGNMENT_ORDER_MODE_DEFAULT;

    protected float alignmentDeltaMz = Aligner.MzFeaturePairSelector.DEFAULT_DELTA_MZ;
    protected float alignmentMinIntensity = Aligner.MzFeaturePairSelector.DEFAULT_MIN_INTENSITY;

    protected double alignmentMaxStudRes = 100;
    protected double alignmentMaxLeverageNumerator = 100;

    protected int quantilePolynomialDegree = AmtDatabaseMatcher.DEFAULT_NONLINEAR_MAPPING_DEGREE;

    //Align by scan? If not, align by RT
    protected boolean alignByScan = false;

    protected static final int ALIGNMENT_MODE_SPLINE = 0;
    protected static final int ALIGNMENT_MODE_QUANTILE = 1;

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

    protected boolean showCharts = false;

    protected boolean shouldDeconvolute = false;

    protected int elutionMode = FeatureClusterer.DEFAULT_ELUTION_MODE;

    protected String[] elutionModeNames = new String[] { "time", "scan" };
    EnumeratedValuesArgumentDefinition elutionModeArg =
            new EnumeratedValuesArgumentDefinition("elutionmode", false, "Cluster features by retention time or by scan",
                                elutionModeNames, "time");

    protected String[] optimizationModeNames = new String[] { "em", "perfectbuckets" };
    EnumeratedValuesArgumentDefinition optimizationModeArg =
            new EnumeratedValuesArgumentDefinition("optimizationmode", false, "Optimize cluster diameters using the " +
                    "'perfect buckets' strategy (original msInspect strategy) or the newer strategy using EM algorithm",
                                optimizationModeNames, "em");

    protected String[] massModeNames = new String[] { "da", "ppm" };
    EnumeratedValuesArgumentDefinition massModeArg =
            new EnumeratedValuesArgumentDefinition("masstype", false, "For mass clustering, use Daltons or PPM?",
                                massModeNames, "ppm");

    protected float maxOptimizeMatchFDR = BucketedPeptideArray.DEFAULT_MAX_MATCH_FDR_EM_OPT;

    EnumeratedValuesArgumentDefinition alignOrderModeArg = new EnumeratedValuesArgumentDefinition("alignOrderMode", false,
                    "Method of selecting run order for feature alignment.  \"alltofirst\" will align all runs " +
                            "to the first run.  \"daisychain\" will align each run to the run before it (e.g., for " +
                            "fractionated samples).  \"cumulative\" will grow the 'run' being aligned to incrementally",
                            Aligner.alignmentOrderModeDescs,
                    Aligner.alignmentOrderModeDescs[alignmentOrderMode]);
    public CommandLineArgumentDefinition[] arrayParameterArgDefs =
    {
            new FileToReadArgumentDefinition("tags",false, "optional file of tags for each run"),
            new EnumeratedValuesArgumentDefinition("alignByTags", false,
                                               "Request a pre-alignment within each group of tags. \"None\" indicates " +
                                                       "that no pre-alignment be done. \"Loose\" requests pre-alignment, " +
                                                       "preserving those features present in at least 3/4 of the " +
                                                       "runs in each group. \"Strict\" preserves those features " +
                                                       "present in all runs in each group. Unless \"None\", a tags " +
                                                       "file *must* be specified.",
                                               alignByTagsStrings, alignByTags),
            alignOrderModeArg,
            new DecimalArgumentDefinition("elutionwindow", false,
                    "number of seconds or scans to use as a window when aligning features", elutionBucket),
            new DecimalArgumentDefinition("masswindow", false,
                    "number of Daltons to use as a window when aligning features", massBucket),
            new BooleanArgumentDefinition("normalize", false,
                    "should the intensities be normalized after alignment?", normalize),
            new IntegerArgumentDefinition("topN", false,
                    "if > 0, sets the number of highest-intensity features from each run to use when constructing " +
                            "the non-linear mapping", topN),
            new BooleanArgumentDefinition("align", false, "should alignment be performed?",align),
            new BooleanArgumentDefinition("optimize", false,
                    "Should we optimize the size of the mass and scan buckets based on the number of 'perfect buckets'?",
                    optimize),
            new EnumeratedValuesArgumentDefinition("intensitytype", false,
                    "What type of intensity should we include in the array when there are conflicts?  " +
                    "A sum of all matching features in the bucket, or the intensity of the feature with lowest kl," +
                            "or the maximum intensity?",
                    intensityTypeStrings, intensityTypeStrings[FeatureGrouper.DEFAULT_CONFLICT_RESOLVER]),
            new StringArgumentDefinition("massbuckets", false,
                    "comma-separated list of decimal values for the maximum mass bucket size"),
            new StringArgumentDefinition("elutionbuckets", false,
                    "comma-separated list of values for the maximum time or scan bucket size"),
            new BooleanArgumentDefinition("optimizeonpeptideids", false,
                    "If optimizing, optimize based on the number of rows with agreeing peptide IDs, rather than the " +
                            "number of rows with one peptide (of any ID) from each row",
                    optimizeOnPeptideIds),
            new IntegerArgumentDefinition("peptidematchscore", false,
                    "If optimizing based on rows with agreeing peptides, give this score for each agreeing row",
                    peptideMatchScore),
            new IntegerArgumentDefinition("peptidemismatchpenalty", false,
                    "If optimizing based on rows with agreeing peptides, give this penalty for each mismatched row",
                    peptideMismatchPenalty),
            featurePairSelectorArgDef,
            new BooleanArgumentDefinition("showcharts", false,
                    "Optionally plot the warping functions for each aligned run",
                    showCharts),
            new DecimalArgumentDefinition("alignmztolerance", false,
                    "M/Z tolerance used in picking pairs of features for alignment (if applicable)",
                    alignmentDeltaMz),
            new DecimalArgumentDefinition("alignminintensity", false,
                    "Minimum intensity used in picking pairs of features for alignment (if applicable)",
                    alignmentMinIntensity),
            new IntegerArgumentDefinition("df", false,
                    "Degrees of Freedom for alignment (spline mode)",
                    degreesOfFreedom),
            new EnumeratedValuesArgumentDefinition("alignmentmode",false,alignmentModeStrings,
                    alignmentModeExplanations, "spline"),
            new DecimalArgumentDefinition("maxstudres", false,
                    "Maximum studentized residual for alignment.  Default value is effectively no filtering.  " +
                            "Use values < 3 for tighter filtering", alignmentMaxStudRes),
            new DecimalArgumentDefinition("maxleverage", false,
                    "Maximum NUMERATOR of the leverage of features used for alignment (denominator is N). Default " +
                            "value is effectively no filtering.  Use values < 8 for tighter filtering",
                    alignmentMaxLeverageNumerator),
            new IntegerArgumentDefinition("polynomialdegree", false,
                    "The degree of the polynomial to fit (for quantile mode)",
                    quantilePolynomialDegree),
            new BooleanArgumentDefinition("deconvolute", false,
                    "Should features be deconvoluted (collapsed to lowest charge state) before array creation?",
                    shouldDeconvolute),
            new BooleanArgumentDefinition("alignbyscan", false,
                    "Align by scan? (If not, align by retention time)", alignByScan),
            new DecimalArgumentDefinition("maxoptimizematchfdr", false, "For EM optimization.  Max FDR " +
                    "that every match within the bounding box of tolerances is a good match",
                    BucketedPeptideArray.DEFAULT_MAX_MATCH_FDR_EM_OPT),
            elutionModeArg,
            optimizationModeArg,
            massModeArg,
    };


    protected static final String[] intensityTypeStrings =
        {
            "sum",
            "best",
                "max",
        };

    protected static final String[] featurePairSelectorStrings =
        {
            "mz",
            "peptide",
            "hybrid"
        };

    protected int conflictResolver = FeatureGrouper.DEFAULT_CONFLICT_RESOLVER;

    public PeptideArrayCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        super.init();
        mCommandName = "peptidearray";

        mHelpMessage =
                "The peptidearray command creates a peptide array from two or more feature files.";

        mShortDescription = "Create a peptide array from two or more feature files";

        List<CommandLineArgumentDefinition> argDefList =
                new ArrayList<CommandLineArgumentDefinition>();
        for (CommandLineArgumentDefinition def : arrayParameterArgDefs)
            argDefList.add(def);
        argDefList.add(new FileToWriteArgumentDefinition("out",true, "output file"));
        argDefList.add(createUnnamedSeriesFileArgumentDefinition(true,
                                "A series of feature files to align"));


        addArgumentDefinitions(argDefList.toArray(new CommandLineArgumentDefinition[0]));
    }

    protected int translateIntensityTypeString(String intensityType)
    {
        int result = -1;
        for (int i = 0; i < intensityTypeStrings.length; i++)
        {
            String modeString = intensityTypeStrings[i];
            if (modeString.equalsIgnoreCase(intensityType))
            {
                result = i;
                break;
            }
        }
        return result;
    }

    protected double[] parseBucketValues(String paramString)
            throws ArgumentValidationException
    {
        if (paramString == null)
            return null;
        double[] result;

        try
        {
            List<Double> resultList = new ArrayList<Double>();
            String[] chunks = paramString.split(",");
            for (String chunk : chunks)
            {
                if (chunk != null && chunk.length() > 0)
                {
                    Double bucketsize = new Double(chunk);
                    if (bucketsize <= 0)
                        throw new ArgumentValidationException("All bucket sizes must be positive");
                    resultList.add(bucketsize);
                }
            }
            result = new double[resultList.size()];
            for (int i=0; i<resultList.size(); i++)
            {
                result[i] = resultList.get(i);
            }
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Error while parsing bucket value parameter string " + paramString,e);
        }

        return result;
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        super.assignArgumentValues();

        outFile = getFileArgumentValue("out");
        tagFile = getFileArgumentValue("tags");

        featureFiles = getUnnamedSeriesFileArgumentValues();

        massBucket = getDoubleArgumentValue("masswindow");
        massType = massModeArg.getIndexForArgumentValue(getStringArgumentValue("masstype"));
//System.err.println("massType = " + massType);

        elutionBucket = getDoubleArgumentValue("elutionwindow");
        topN = getIntegerArgumentValue("topN");
        normalize = getBooleanArgumentValue("normalize");
        align = getBooleanArgumentValue("align");
        optimize = getBooleanArgumentValue("optimize");
        optimizeOnPeptideIds = getBooleanArgumentValue("optimizeonpeptideids");
        maxOptimizeMatchFDR = getFloatArgumentValue("maxoptimizematchfdr");
        peptideMatchScore = getIntegerArgumentValue("peptidematchscore");
        peptideMismatchPenalty = getIntegerArgumentValue("peptidemismatchpenalty");
        degreesOfFreedom = getIntegerArgumentValue("df");

        shouldDeconvolute = getBooleanArgumentValue("deconvolute");

        alignmentMode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("alignmentmode")).getIndexForArgumentValue(
                getStringArgumentValue("alignmentmode"));


        alignmentMaxStudRes = getDoubleArgumentValue("maxstudres");
        alignmentMaxLeverageNumerator = getDoubleArgumentValue("maxleverage");

        alignByScan = getBooleanArgumentValue("alignbyscan");

        alignByTags = getStringArgumentValue("alignByTags");
        if (!"none".equals(alignByTags) && null == tagFile)
            throw new ArgumentValidationException("When aligning by tags, a tag file must be specified with --tag");

        alignmentOrderMode = alignOrderModeArg.getIndexForArgumentValue(getStringArgumentValue("alignOrderMode")); 

        elutionMode = elutionModeArg.getIndexForArgumentValue(getStringArgumentValue("elutionmode"));

        conflictResolver =
                translateIntensityTypeString(getStringArgumentValue("intensitytype"));

        //if we're not optimizing, then the -buckets arguments are meaningless
        if (!optimize)
        {
            assertArgumentAbsent("massbuckets", "optimize");
            assertArgumentAbsent("scanbuckets", "optimize");
            assertArgumentAbsent("optimizationmode", "optimize");
        }
        optimizationMode = optimizationModeArg.getIndexForArgumentValue(getStringArgumentValue("optimizationmode"));
        if (optimizationMode == BucketedPeptideArray.OPTIMIZE_MODE_ERRORDIST)
        {
            assertArgumentAbsent("massbuckets", "optimizationmode");
            assertArgumentAbsent("scanbuckets", "optimizationmode");
        }

        massBuckets = parseBucketValues(getStringArgumentValue("massbuckets"));
        elutionBuckets = parseBucketValues(getStringArgumentValue("elutionbuckets"));

        quantilePolynomialDegree = getIntegerArgumentValue("polynomialdegree");


        alignmentDeltaMz = getFloatArgumentValue("alignmztolerance");
//System.err.println("Alignment delta mz: " + alignmentDeltaMz);        
        alignmentMinIntensity = getFloatArgumentValue("alignminintensity");


//System.err.println("adsf");
        if (hasArgumentValue(featurePairSelectorArgDef.getArgumentName()))
        {
//System.err.println("has pair arg, which is " +
            int argIndex = featurePairSelectorArgDef.getIndexForArgumentValue(getStringArgumentValue(featurePairSelectorArgDef.getArgumentName()));
            if (argIndex != 0 && topN > 0)
                throw new ArgumentValidationException("topN may only be speicified if aligning based on mz");
            switch (argIndex)
            {
                case 0:
                    //same as default
                    break;
                case 1:
                    featurePairSelector = new Aligner.PeptideFeaturePairSelector();
                    break;
                case 2:
                    featurePairSelector = new Aligner.HybridFeaturePairSelector();
                    break;
            }
        }
        
        // Set alignment deltaMz and minIntensity explicitly in case the user has changed them on the command-line
        if (featurePairSelector instanceof Aligner.MzFeaturePairSelector)
        {
            ((Aligner.MzFeaturePairSelector) featurePairSelector).setDeltaMz(alignmentDeltaMz);
            ((Aligner.MzFeaturePairSelector) featurePairSelector).setMinIntensity(alignmentMinIntensity);
        }

        if (topN > 0)
        {
            if (!(featurePairSelector instanceof Aligner.MassOrMzFeaturePairSelector))
                    throw new ArgumentValidationException("To specify topN, featureselector must be mass or mz");
            //if topN > 0, then we must have an MzFeaturePairSelector
            ((Aligner.MassOrMzFeaturePairSelector) featurePairSelector).setTopN(topN);
        }

        showCharts = getBooleanArgumentValue("showcharts");
    }

    /**
     * Read the tags file and associate with each input filename. Tags should contain one row
     * per FeatureSet, each row containing a filename prefix and a tag separated by a tab.
     */
    private Map<String, String> assignTags()
        throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(tagFile));
        HashMap<String,String> tags = new HashMap<String,String>();
        String line = null;
        while(null != (line = br.readLine()))
        {
            if (line.startsWith("#"))
                continue;
            String[] fields;
            if (line.indexOf("\t") >= 0)
                fields = line.split("\t");
            else
                fields = line.split(" ");
            if (fields.length != 2)
            {
                _log.warn("Skipping unrecognized line in tag file: '" + line + "'");
                continue;
            }
            for (File f : featureFiles)
                if (f.getName().startsWith(fields[0]))
                    tags.put(f.getName(), fields[1]);
        }
        return tags;
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            List<File> featureFileList = new ArrayList<File>(featureFiles.length);

            for (File featureFile : featureFiles)
            {
                featureFileList.add(featureFile);
            }
            BucketedPeptideArray arr = new BucketedPeptideArray(featureFileList, featureSelector);
            arr.setElutionMode(elutionMode);            
            arr.setElutionBucket(elutionBucket);
            arr.setMassType(massType);
            arr.setMassBucket(massBucket);
            arr.setOutFileName(outFile.getAbsolutePath());
            arr.setAlign(align);
            arr.setNormalize(normalize);
            arr.setConflictResolver(conflictResolver);
            arr.setFeaturePairSelector(featurePairSelector);
            arr.setShouldDeconvolute(shouldDeconvolute);
            arr.setOptimizationMode(optimizationMode);
            arr.setMaxMatchFDRForToleranceBoxCalc(maxOptimizeMatchFDR);


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

            if (alignByScan)
                aligner.setAlignmentDatasource(Aligner.ALIGNMENT_DATASOURCE_SCAN);
            else
                aligner.setAlignmentDatasource(Aligner.ALIGNMENT_DATASOURCE_TIME);

            ((Aligner.MassOrMzFeaturePairSelector) aligner.getFeaturePairSelector()).setTopN(topN);
            aligner.setMaxLeverageNumerator(alignmentMaxLeverageNumerator);
            aligner.setMaxStudRes(alignmentMaxStudRes);
            aligner.setBuildCharts(showCharts);

            aligner.setAlignmentOrderMode(alignmentOrderMode);

            arr.setAligner(aligner);

            if (optimize)
            {
                if (optimizeOnPeptideIds)
                {
                    arr.getFeatureGrouper().setBucketEvaluationPeptideAgreementMatchScore(peptideMatchScore);
                    arr.getFeatureGrouper().setBucketEvaluationPeptideAgreementMismatchPenalty(peptideMismatchPenalty);
                }
                if (massBuckets != null)
                    arr.setMassBuckets(massBuckets);
                if (elutionBuckets != null)
                    arr.setScanBuckets(elutionBuckets);
            }

            if ("none".equals(alignByTags))
            {
                arr.run(optimize, optimizeOnPeptideIds ?
                            FeatureGrouper.BUCKET_EVALUATION_MODE_PEPTIDE_AGREEMENT :
                            FeatureGrouper.BUCKET_EVALUATION_MODE_ONE_FROM_EACH, showCharts);
            }
            else
            {
                boolean strict = "strict".equals(alignByTags);
                Map<String, String> tags = assignTags();
                arr.alignByGroup(tags, strict);
            }
//System.err.println("Done with array, plotwarping: " + plotWarping + ", align: " + align);

        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Error while generating peptide array: " + e.getMessage());
            throw new CommandLineModuleExecutionException(e);
        }
    }

    public File getOutFile()
    {
        return outFile;
    }

    public void setOutFile(File outFile)
    {
        this.outFile = outFile;
    }

    public void setInputFiles(File[] inputFiles)
    {
        featureFiles = inputFiles;
    }
}
