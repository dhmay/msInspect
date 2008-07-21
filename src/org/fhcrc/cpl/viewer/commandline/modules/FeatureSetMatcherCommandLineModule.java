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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.DeltaMassArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.util.MsInspectRegressionUtilities;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;
import org.fhcrc.cpl.viewer.align.Aligner;
import org.fhcrc.cpl.viewer.align.SplineAligner;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.apache.log4j.Logger;
import org.jfree.data.xy.XYSeriesCollection;

import java.util.*;
import java.io.File;
import java.io.IOException;


/**
 */
public class FeatureSetMatcherCommandLineModule extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FeatureSetMatcherCommandLineModule.class);

    protected FeatureSet ms1Features = null;
    protected FeatureSet ms2Features = null;
    protected File outFile = null;
    protected File ms2directory = null;
    protected MSRun run1 = null;
    protected double minPeptideProphet = 0;
    float deltaMass = AmtFeatureSetMatcher.DEFAULT_DELTA_MASS_PPM;
    int deltaMassType = AmtFeatureSetMatcher.DELTA_MASS_TYPE_PPM;
    int deltaScan = AmtFeatureSetMatcher.DEFAULT_DELTA_SCAN;
    double deltaTime = AmtFeatureSetMatcher.DEFAULT_DELTA_TIME;
    double deltaHydrophobicity = AmtFeatureSetMatcher.DEFAULT_DELTA_HYDROPHOBICITY;
    protected boolean useTime = false;
    protected boolean matchOnHydro = false;
    protected double deltaHydro = 0;
    protected boolean writeUnmatched = true;
    protected boolean stripMultipleMS2 = true;

    protected FeatureSet[] ms2FeatureSets;
    boolean alignMS2=false;

    ClusteringFeatureSetMatcher cFSM = null;

    protected File outUnmatchedMS2File = null;
    protected File outAllMS2MarkedFile = null;

    protected boolean showCharts = false;

    protected Protein[] proteinsInFasta = null;


    protected static final int MODE_MS1_MS2 = 0;
    protected static final int MODE_MS1_MS2_DIR = 1;

    protected int mode=MODE_MS1_MS2;

    protected final static String[] modeStrings = {"ms1ms2","ms1ms2dir"};

    protected List<MS2Modification> ms2ModificationsForMatching;
    protected AmtDatabase amtDatabaseForPeptideExclusion = null;


    public FeatureSetMatcherCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matchfeatures";
        mShortDescription = "perform simple feature-matching";

        mHelpMessage =
                "perform feature-matching between different sets of fe" +
                        "atures";

        CommandLineArgumentDefinition[] argDefs =
            {
                    createEnumeratedArgumentDefinition("mode",false,
                            "Mode of operation, default ms1ms2",modeStrings, "ms1ms2"),
                    createFeatureFileArgumentDefinition("ms1features",true,
                            "MS1 feature file"),
                    createFeatureFileArgumentDefinition("ms2features",false,
                            "MS2 feature file (usually pepXml)"),
                    createDirectoryToReadArgumentDefinition("ms2dir",false,
                            "MS2 feature file directory"),
                    createBooleanArgumentDefinition("align", false,
                            "align all MS2 runs to the MS1 run being matched", alignMS2),
                    createFileToReadArgumentDefinition("mzxml",false,
                            "mzXML file"),
                    createFileToWriteArgumentDefinition("out", false,
                            "Output File"),
                    createDecimalArgumentDefinition("minpprophet", false,
                            "Minimum PeptideProphet score", minPeptideProphet),
                    createDecimalArgumentDefinition("deltatime", false,
                            "Maximum time between matched features", deltaTime),
                    createIntegerArgumentDefinition("deltascan", false,
                            "Maximum number of scans between matched features", deltaScan),
                    createDeltaMassArgumentDefinition("deltamass", false,
                              "Maximum mass difference between matched features (in units of da [Daltons] or ppm [parts per million]",
                            new DeltaMassArgumentDefinition.DeltaMassWithType(deltaMass, deltaMassType)),
                    createBooleanArgumentDefinition("matchonhydro", false,
                            "under the hood, perform matching based on hydrophobicity", matchOnHydro),
                    createBooleanArgumentDefinition("writeunmatched", false,
                            "Write out unmatched features", writeUnmatched),
                    createBooleanArgumentDefinition("stripmultiplems2", false,
                            "Strip subsequent MS2 identifications for the same peptide out of the file when matching",
                            stripMultipleMS2),
                    createFileToWriteArgumentDefinition("outunmatchedms2", false,
                            "Output File for unmatched MS2"),
                    createFileToWriteArgumentDefinition("outallms2marked", false,
                            "Output File for all MS2, with unmatched having 'unmatched' in description"),
                    createBooleanArgumentDefinition("showcharts", false,
                            "show useful charts created when matching", showCharts),
                    createFastaFileArgumentDefinition("fasta", false,
                            "Fasta database for matching"),
                    createModificationListArgumentDefinition("modifications", false,
                            "a list of modifications to use when creating features to represent peptide sequences"),
                    createFileToReadArgumentDefinition("amtdbtoexclude", false,
                            "an AMT database whose peptides should be excluded from protein matching")
            };
        addArgumentDefinitions(argDefs);
    }

//    public void digestArguments(Map<String,String> argumentMap)
//            throws ArgumentValidationException
//    {
//        ApplicationContext.infoMessage("Damon's little matching tool v12");
//        //remember, month is (ridiculously) zero-based. Day is not.  Year is not.  Just month.
//        Date deathDate = new GregorianCalendar(2007, 5-1, 1).getTime();
//        if (new Date().getTime() > deathDate.getTime())
//        {
//            throw new RuntimeException("Sorry, this feature expired on " + deathDate +
//                    ".  Please go back to Damon and get a new version");
//        }
//        else
//            ApplicationContext.infoMessage("Note: this tool expires on " + deathDate + ", please check with Damon if you need to use it after that date");
//
//        super.digestArguments(argumentMap);
//    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        ms1Features = getFeatureSetArgumentValue("ms1features");

        //under the hood, match based on hydrophobicity
        matchOnHydro = getBooleanArgumentValue("matchonhydro");

        if (hasArgumentValue("mode"))
            mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        if (hasArgumentValue("minpprophet"))
            minPeptideProphet = (float) getDoubleArgumentValue("minpprophet");
        FeatureSet.FeatureSelector peptideProphetFeatureSelector = new FeatureSet.FeatureSelector();
        peptideProphetFeatureSelector.setMinPProphet((float) minPeptideProphet);

        stripMultipleMS2 = getBooleanArgumentValue("stripmultiplems2");



        if (hasArgumentValue("mzxml"))
        {
            try
            {
                run1 = MSRun.load((getFileArgumentValue("mzxml")).getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }        

        switch(mode)
        {
            case MODE_MS1_MS2:
                assertArgumentPresent("ms2features");
                assertArgumentAbsent("ms2dir");
                ms2Features = getFeatureSetArgumentValue("ms2features");
                ms2Features = ms2Features.filter(peptideProphetFeatureSelector);
                if (ms2Features.getFeatures().length == 0)
                    throw new ArgumentValidationException("There are no MS2 features with PeptideProphet score >= " + minPeptideProphet + ", quitting.");
                if (run1 != null)
                    ms2Features.populateTimesForMS2Features(run1);
                if (stripMultipleMS2)
                    MS2ExtraInfoDef.removeAllButFirstFeatureForEachPeptide(ms2Features);
                break;
            case MODE_MS1_MS2_DIR:
                assertArgumentPresent("ms2dir");
                assertArgumentAbsent("ms2features");

                matchOnHydro = true;
                try
                {
                    ms2directory = getFileArgumentValue("ms2dir");
                    File[] ms2files = ms2directory.listFiles();
                    ms2FeatureSets = new FeatureSet[ms2files.length];
                    for (int i = 0; i < ms2files.length; i++)
                    {
                        ms2FeatureSets[i] = new FeatureSet(ms2files[i]);
                        ms2FeatureSets[i] = ms2FeatureSets[i].filter(peptideProphetFeatureSelector);
                        if (stripMultipleMS2)
                            MS2ExtraInfoDef.removeAllButFirstFeatureForEachPeptide(ms2FeatureSets[i]);
                    }
                }
                catch (Exception e)
                {
                    throw new ArgumentValidationException(e);
                }
                break;
        }


        outUnmatchedMS2File = getFileArgumentValue("outunmatchedms2");
        outAllMS2MarkedFile = getFileArgumentValue("outallms2marked");

        alignMS2 = getBooleanArgumentValue("align");







        ApplicationContext.infoMessage("MS1 features: " + ms1Features.getFeatures().length);
//        ApplicationContext.infoMessage("MS2 features: " + ms2Features.getFeatures().length);

        outFile = getFileArgumentValue("out");

        if (hasArgumentValue("deltatime"))
        {
            if (hasArgumentValue("deltascan"))
                throw new ArgumentValidationException("deltatime and deltascan are mutually exclusive");
            deltaTime = getDoubleArgumentValue("deltatime");
            if (hasArgumentValue("mzxml"))
            {
              ms2Features.populateTimesForMS2Features(run1);
              ms1Features.populateTimesForMS1Features(run1);
            }
            useTime = true;
        }
        if (hasArgumentValue("deltascan"))
        {
            if (hasArgumentValue("deltatime"))
                throw new ArgumentValidationException("deltatime and deltascan are mutually exclusive");
            deltaScan = getIntegerArgumentValue("deltascan");
        }


        if (matchOnHydro)
        {
            List<Feature> featuresWithPeptides = new ArrayList<Feature>();
            for (Feature feature : ms2Features.getFeatures())
                if (MS2ExtraInfoDef.getFirstPeptide(feature) != null)
                    featuresWithPeptides.add(feature);
            Feature[] featuresWithPeptidesArray = featuresWithPeptides.toArray(new Feature[featuresWithPeptides.size()]);
            double[] studentizedResiduals =
                    AmtDatabaseBuilder.calculateStudentizedResiduals(featuresWithPeptidesArray,
                            MsInspectRegressionUtilities.REGRESSION_MODE_TIME);

            Feature[] featuresForRegression =
                    AmtDatabaseBuilder.chooseFeaturesWithMaxStudentizedResidual(featuresWithPeptidesArray,
                            studentizedResiduals,
                            AmtDatabaseBuilder.DEFAULT_MAX_STUDENTIZED_RESIDUAL_FOR_REGRESSION);
            Map<String, Double> ms2HydroLineMap =
                    AmtUtilities.calculateScanOrTimeHydrophobicityRelationship(
                            featuresForRegression,
                            MsInspectRegressionUtilities.REGRESSION_MODE_TIME,
                            false);
            AmtUtilities.recordHydrophobicities(ms2Features,
                    ms2HydroLineMap,
                    MsInspectRegressionUtilities.REGRESSION_MODE_TIME);
            double[] hydros = new double[ms2Features.getFeatures().length];
            for (int i=0; i<hydros.length; i++) hydros[i] = AmtExtraInfoDef.getObservedHydrophobicity(ms2Features.getFeatures()[i]);

            _log.debug("recording hydros for ms1 features");
            AmtUtilities.recordHydrophobicities(ms1Features,
                    ms2HydroLineMap,
                    useTime? MsInspectRegressionUtilities.REGRESSION_MODE_TIME :
                            MsInspectRegressionUtilities.REGRESSION_MODE_SCAN);

            deltaHydro =
                    AmtUtilities.convertDeltaScanOrTimeToHydro(ms2HydroLineMap,
                            useTime? deltaTime : deltaScan) -
                    AmtUtilities.convertDeltaScanOrTimeToHydro(ms2HydroLineMap,
                                                               0);
            _log.debug("deltaHydro: " + deltaHydro);
        }


        writeUnmatched = getBooleanArgumentValue("writeunmatched");



//        deltaElution = AmtCommandLineModule.convertDeltaScanOrTimeToHydro(
//                AmtUtilities.getSlopeFromRegressionLine(regressionLineMap),
//                AmtUtilities.getInterceptFromRegressionLine(regressionLineMap),
//                deltaTime);
//System.err.println("converted " + deltaTime + " seconds to " + deltaElution + " hydrophobicity units");

        //this could also be broken out into a new ArgumentDefinition class
        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassAndType =
                getDeltaMassArgumentValue("deltamass");
        deltaMass = deltaMassAndType.getDeltaMass();
        deltaMassType = deltaMassAndType.getDeltaMassType();


        showCharts = getBooleanArgumentValue("showCharts");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        int elutionMode = matchOnHydro? ClusteringFeatureSetMatcher.ELUTION_MODE_HYDROPHOBICITY :
                (useTime? ClusteringFeatureSetMatcher.ELUTION_MODE_TIME :
                          ClusteringFeatureSetMatcher.ELUTION_MODE_SCAN);
        float deltaElution = -1;
        double elutionBucketIncrement = -1;
        switch (elutionMode)
        {
            case ClusteringFeatureSetMatcher.ELUTION_MODE_HYDROPHOBICITY:
                deltaElution = (float) deltaHydro;
                elutionBucketIncrement =
                        ClusteringFeatureSetMatcher.DEFAULT_HYDRO_ELUTION_BUCKET_INCREMENT;
                break;
            case ClusteringFeatureSetMatcher.ELUTION_MODE_TIME:
                deltaElution = (float) deltaTime;
                elutionBucketIncrement = 100;
                break;
            case ClusteringFeatureSetMatcher.ELUTION_MODE_SCAN:
                deltaElution = (float) deltaScan;
                elutionBucketIncrement = 1;
                break;
        }

        cFSM = new ClusteringFeatureSetMatcher(deltaMass, deltaMassType,
                        deltaElution);

        cFSM.setElutionMode(elutionMode);
        cFSM.setElutionBucketIncrement(elutionBucketIncrement);



        switch(mode)
        {
            case MODE_MS1_MS2:
                matchOneToOne(ms1Features, ms2Features);
                break;
            case MODE_MS1_MS2_DIR:
                matchOneToMany();
                if (outFile != null)
                {
                    try
                    {
                        ms1Features.save(outFile);
                        ApplicationContext.infoMessage("Saved " +  ms1Features.getFeatures().length +
                                " features, with annotations for matches, to file " + outFile.getAbsolutePath());
                    }
                    catch (Exception e)
                    {
                        throw new CommandLineModuleExecutionException(e);
                    }
                }
                break;
        }



    }

    protected void matchOneToMany()
            throws CommandLineModuleExecutionException
    {
        Set<String> peptidesMatched = new HashSet<String>();

        for (FeatureSet ms2FeatureSet : ms2FeatureSets)
        {
            matchOneToOne(ms1Features, ms2FeatureSet);
        }

        for (Feature ms1Feature : ms1Features.getFeatures())
        {
            List<String> ms1PeptideList = MS2ExtraInfoDef.getPeptideList(ms1Feature);
            if (ms1PeptideList != null)
                for (String peptide : ms1PeptideList)
                    peptidesMatched.add(peptide);
        }

        System.err.println("Distinct peptides matched: " + peptidesMatched.size());
    }

    protected void matchOneToOne(FeatureSet ms1Features, FeatureSet ms2Features)
             throws CommandLineModuleExecutionException
    {
        ms1Features.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());

        if (alignMS2)
        {
            Aligner aligner = new SplineAligner();
            List<FeatureSet> featureSets = new ArrayList<FeatureSet>();
            featureSets.add(ms1Features);
            featureSets.add(ms2Features);
            try
            {
                aligner.alignFeatureSets(featureSets, showCharts);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }

            if (showCharts)
            {
                aligner.plotWarpings();
            }
        }

        _log.debug("ms1Features: " + ms1Features.getFeatures().length + ", ms2Features: " + ms2Features.getFeatures().length);
        Set<String> allMS2Peptides = new HashSet<String>();
        for (Feature ms2Feature : ms2Features.getFeatures())
        {
            allMS2Peptides.add(MS2ExtraInfoDef.getFirstPeptide(ms2Feature));
        }


        XYSeriesCollection dataset = new XYSeriesCollection();


        AmtFeatureSetMatcher.FeatureMatchingResult featureMatchingResult =
                cFSM.matchFeatures(ms1Features, ms2Features);

        Set<Feature> matchedMs1FeatureHashSet = new HashSet<Feature>();
        Set<String> matchedPeptides = new HashSet<String>();
        Set<Feature> matchedMs2FeatureHashSet = new HashSet<Feature>();

        int numUnmatchedMatched=0;
        int numMatchedMatched=0;
        int numUnmatched=0;
        int numMatched=0;

        for (Feature ms2Feature : ms2Features.getFeatures())
        {
    if (ms2Feature.getDescription() != null && ms2Feature.getDescription().contains("unmatched"))
        numUnmatched++;
    else
        numMatched++;
        }

//        Set<Feature> ms2MatchedSet = new HashSet<Feature>();

        for (Feature ms1Feature : featureMatchingResult.getMasterSetFeatures())
        {
            matchedMs1FeatureHashSet.add(ms1Feature);

            for (Feature ms2Feature : featureMatchingResult.get(ms1Feature))
            {
                matchedMs2FeatureHashSet.add(ms2Feature);

                String ms2Peptide = MS2ExtraInfoDef.getFirstPeptide(ms2Feature);
                List<String> ms1PeptideList = MS2ExtraInfoDef.getPeptideList(ms1Feature);
//if (!ms2MatchedSet.contains(ms2Feature))
//{
//    ms2MatchedSet.add(ms2Feature);
//    if (ms2Feature.getDescription() != null && ms2Feature.getDescription().contains("unmatched"))
//        numUnmatchedMatched++;
//    else
//        numMatchedMatched++;
//}
                if (ms1PeptideList != null && ms1PeptideList.contains(ms2Peptide))
                    continue;

                MS2ExtraInfoDef.addPeptideWithProtein(ms1Feature,
                        ms2Peptide,
                        MS2ExtraInfoDef.getFirstProtein(ms2Feature));
                //TODO: what the heck do I do about this?
                MS2ExtraInfoDef.setPeptideProphet(ms1Feature,
                        MS2ExtraInfoDef.getPeptideProphet(ms2Feature));
                matchedPeptides.add(MS2ExtraInfoDef.getFirstPeptide(ms2Feature));

            }
//if (featureMatchingResult.get(ms1Feature).size() > 1) System.err.println("Matched MS1 with " + featureMatchingResult.get(ms1Feature).size() + " MS2 , list=" + MS2ExtraInfoDef.convertStringListToString(MS2ExtraInfoDef.getPeptideList(ms1Feature)));

        }

        ApplicationContext.infoMessage("Matched " + matchedMs1FeatureHashSet.size() + " out of " + ms1Features.getFeatures().length + " MS1 features.");
        ApplicationContext.infoMessage("\t(" + matchedMs2FeatureHashSet.size() + " out of " + ms2Features.getFeatures().length + " MS2 features, (" + (matchedMs2FeatureHashSet.size()* 100 / ms2Features.getFeatures().length) +"%))");

        ApplicationContext.infoMessage("\t(" + matchedPeptides.size() +
                " distinct new peptides, out of " +
                allMS2Peptides.size() + ", " +
                (100 * matchedPeptides.size() / allMS2Peptides.size()) + "%)");
//System.err.println("Self-matched matched: " + numMatchedMatched + " / " + numMatched + " ( " + 100 * (double)((double)numMatchedMatched / (double)numMatched) + "%)");
//System.err.println("Non-self-matched matched: " + numUnmatchedMatched + " / " + numUnmatched + " ( " + 100 * (double)((double)numUnmatchedMatched / (double)numUnmatched) + "%)");


        double[] matchedIntensities = new double[matchedMs1FeatureHashSet.size()];
        int i=0;
        for (Feature matchedMs1Feature : matchedMs1FeatureHashSet)
            matchedIntensities[i++] = matchedMs1Feature.getIntensity();
        double meanMatchedIntensity = BasicStatistics.mean(matchedIntensities);

        System.err.println("Mean intensity of matched features: " + meanMatchedIntensity);

        Set<Feature> unmatchedMs1Features = new HashSet<Feature>();
        for (Feature ms1Feature : ms1Features.getFeatures())
            if (!matchedMs1FeatureHashSet.contains(ms1Feature))
                unmatchedMs1Features.add(ms1Feature);
        double[] unmatchedIntensities = new double[unmatchedMs1Features.size()];
        i=0;
        for (Feature unmatchedMs1Feature : unmatchedMs1Features)
            unmatchedIntensities[i++] = unmatchedMs1Feature.getIntensity();
        double meanUnmatchedIntensity = BasicStatistics.mean(unmatchedIntensities);

        System.err.println("Mean intensity of unmatched features: " + meanUnmatchedIntensity);



        if (writeUnmatched)
        {
            for (Feature feature : ms1Features.getFeatures())
            {
                matchedMs1FeatureHashSet.add(feature);
            }
        }

        Feature[] annotatedFeatureArray = new Feature[matchedMs1FeatureHashSet.size()];
        int afaIndex = 0;
        for (Feature feature : matchedMs1FeatureHashSet)
            annotatedFeatureArray[afaIndex++] = feature;

        FeatureSet annotatedFeatureSet = new FeatureSet(annotatedFeatureArray);
        annotatedFeatureSet.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());

        if (outFile != null)
        {
            try
            {
                annotatedFeatureSet.save(outFile);
                ApplicationContext.infoMessage("Saved " +  annotatedFeatureSet.getFeatures().length +
                        " features, with annotations for matches, to file " + outFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        if (outUnmatchedMS2File != null)
        {
            _log.debug("Writing unmatched MS2 features");
            List<Feature> unmatchedMs2FeatureList = new ArrayList<Feature>();
            for (Feature ms2Feature : ms2Features.getFeatures())
            {
                if (!matchedPeptides.contains(MS2ExtraInfoDef.getFirstPeptide(ms2Feature)))
                    unmatchedMs2FeatureList.add(ms2Feature);
            }
            try
            {
                new FeatureSet(unmatchedMs2FeatureList.toArray(new Feature[0])).save(outUnmatchedMS2File);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        if (outAllMS2MarkedFile != null)
        {
            _log.debug("Writing all MS2 features, marking unmatched with description 'unmatched'");
            FeatureSet ms2FeaturesClone = (FeatureSet) ms2Features.clone();
            for (Feature ms2Feature : ms2FeaturesClone.getFeatures())
            {
                if (!matchedPeptides.contains(MS2ExtraInfoDef.getFirstPeptide(ms2Feature)))
                    ms2Feature.setDescription("unmatched");
            }
            try
            {
                ms2FeaturesClone.save(outAllMS2MarkedFile);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }

        if (showCharts)
        {
            List<Pair<Float, Float>> matchedFeatureTimes =
                    new ArrayList<Pair<Float,Float>>();
            for (Feature ms1Feature : featureMatchingResult.getMasterSetFeatures())
            {
                matchedMs1FeatureHashSet.add(ms1Feature);

                for (Feature ms2Feature : featureMatchingResult.get(ms1Feature))
                {
                    matchedFeatureTimes.add(new Pair<Float,Float>(ms1Feature.getTime(), ms2Feature.getTime()));
                }
            }
            float[][] scatterPlotData = new float[2][matchedFeatureTimes.size()];
            double[] histData = new double[matchedFeatureTimes.size()];

            for (int j=0; j<matchedFeatureTimes.size(); j++)
            {
                Pair<Float,Float> pair = matchedFeatureTimes.get(j);
                scatterPlotData[0][j] = pair.first;
                scatterPlotData[1][j] = pair.second;

                histData[j] = Math.abs(pair.first-pair.second);
            }

            ScatterPlotDialog spd =
                    new ScatterPlotDialog(scatterPlotData[0], scatterPlotData[1], "");
            spd.setAxisLabels("MS1 time","MS2 Time");
            spd.setVisible(true);
        }


        //if there are peptides in the ms1 features, explore the agreement and disagreement
        //with ms2 peptides
        Set<String> ms1Peptides = new HashSet<String>();
        for (Feature ms1Feature : ms1Features.getFeatures())
            if (MS2ExtraInfoDef.getFirstPeptide(ms1Feature) != null)
                ms1Peptides.add(MS2ExtraInfoDef.getFirstPeptide(ms1Feature));
        if (ms1Peptides.size() > 0)
        {
            int numAgreement=0;
            int numConflict=0;
            int numAmbiguous=0;
            for (Feature ms1Feature : featureMatchingResult.getMasterSetFeatures())
            {
                boolean ambiguous = false;
                List<Feature> matchedMS2Features = featureMatchingResult.get(ms1Feature);
                Set<String> ms2PeptidesSet = new HashSet<String>();
                for (Feature feature : matchedMS2Features)
                {
                    ms2PeptidesSet.add(MS2ExtraInfoDef.getFirstPeptide(feature));
                }
                if (ms2PeptidesSet.size() > 1)
                {
                    ambiguous=true;
                }
                List<String> ms2PeptideList = new ArrayList<String>();
                ms2PeptideList.addAll(ms2PeptidesSet);

                List<String> ms1PeptideList = MS2ExtraInfoDef.getPeptideList(ms1Feature);

                if (ms1PeptideList == null)
                    continue;

                if (ms1PeptideList.size() > 1)
                    ambiguous=true;

                Set<String> commonPeptides = new HashSet<String>();
                boolean agreement = false;

                for (String ms1Peptide : ms1PeptideList)
                {
                    if (ms2PeptidesSet.contains(ms1Peptide))
                        commonPeptides.add(ms1Peptide);
                }
                if (commonPeptides.size() > 0)
                    agreement=true;

                if (ambiguous)
                {
                    numAmbiguous++;
                    ApplicationContext.infoMessage("Ambiguous:  MS1: " +
                            MS2ExtraInfoDef.convertStringListToString(ms1PeptideList) + "\n\tamtmsw4 MS2: " +
                            MS2ExtraInfoDef.convertStringListToString(ms2PeptideList));
                }
                if (agreement)
                    numAgreement++;
                else
                    numConflict++;
            }
            ApplicationContext.infoMessage("Peptide comparison, MS1 to MS2:");
            ApplicationContext.infoMessage("Agreement: " + numAgreement + ", conflict: " +
                    numConflict + ", ambiguous: " + numAmbiguous);
        }


    }

/*
    protected void matchManyToMany()
            throws CommandLineModuleExecutionException
    {
        for (int i=0; i<featureSets.length; i++)
        {
            Map<Feature, String> newPeptideAssignments =
                    new HashMap<Feature,String>();
            Set<Feature> conflictedFeatures = new HashSet<Feature>();

            for (int j=0; j<featureSets.length; j++)
            {
                if (i == j)
                    continue;

                AmtFeatureSetMatcher.FeatureMatchingResult featureMatchingResult =
                    cfsm.matchFeatures(featureSets[i], featureSets[j]);

                for (Feature set1Feature : featureMatchingResult.getMasterSetFeatures())
                {
                    if (MS2ExtraInfoDef.getFirstPeptide(set1Feature) != null)
                        continue;
                    if (conflictedFeatures.contains(set1Feature))
                        continue;

                    for (Feature set2Feature : featureMatchingResult.get(set1Feature))
                    {
                        String set2Peptide =
                                MS2ExtraInfoDef.getFirstPeptide(set2Feature);
                        if (set2Peptide != null)
                        {
                            if (newPeptideAssignments.containsKey(set1Feature))
                            {
                                if (!newPeptideAssignments.get(set1Feature).equals(set2Peptide))
                                    conflictedFeatures.add(set1Feature);
                            }
                            else
                                newPeptideAssignments.put(set1Feature, set2Peptide);
                        }
                    }
                }
            }

            for (Feature featureWithNewPeptide : newPeptideAssignments.keySet())
            {
                MS2ExtraInfoDef.addPeptide(featureWithNewPeptide,
                        newPeptideAssignments.get(featureWithNewPeptide));
            }

            try
            {
                File out = new File(outDir.getAbsolutePath() + File.separatorChar + featureSets[i].getSourceFile().getName());
                System.err.println("Saving file " + out.getName() + " with " + newPeptideAssignments.size() + " peptide additions");
                featureSets[i].save(out);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }

        }
    }

    */


}