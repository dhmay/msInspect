/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandModuleUtilities;
import org.fhcrc.cpl.viewer.feature.extraction.WaveletPeakExtractor;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategy;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.Window2DFeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

import java.io.File;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class FlagQuantEventsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FlagQuantEventsCLM.class);

    protected File[] featureFiles;
    protected File outFile;
    protected File outDir;
    protected File mzXmlDir;
//    protected File ms1FeaturesDir;

    public static final int FLAG_REASON_NONE = 0;
    public static final int FLAG_REASON_COELUTING_ABOVEBELOW = 1;
    public static final int FLAG_REASON_DISSIMILAR_KL = 2;
    public static final int FLAG_REASON_DISSIMILAR_MS1_RATIO = 3;
    public static final int FLAG_REASON_DIFF_HIGHESTPEAK = 4;
    public static final int FLAG_REASON_DIFF_PEAKDIST = 5;




    protected int numScansSlopForCoeluting = 2;
    protected float massSlopPPM = 25f;
    protected float massSlopDa = 0.05f;


    public static final float MOD_MASS_SLOP = 0.1f;

    //Minimum proportion of the quantitated feature's intensity (if that's determinable) " +
    //for a coeluting feature to have to cause an event to be flagged.
    protected float minCoeluteIntensityProportion = 0.1f;


    protected int labelType = LABEL_LYCINE;

    protected boolean showCharts = false;

    public static final int LABEL_ACRYLAMIDE = 0;
    public static final int LABEL_LYCINE = 1;

    public static final float SILAC_LABEL_MASS = 134.115092f;
    public static final float SILAC_LABEL_MASSDIFF_PERRESIDUE = 6.020129f;


    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;
    public static final float ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE =
            ACRYLAMIDE_LABEL_HEAVYMASS - ACRYLAMIDE_LABEL_LIGHTMASS;


    protected String[] labelStrings = new String[]
            {
                    "acrylamide",
                    "silac"
            };

    protected String[] labelExplanations = new String[]
            {
                    "Acrylamide (3.0106Da on C)",
                    "SILAC Lycine labeling (134.115092 on K, i.e. 6Da SILAC)"
            };

    public static final String[] flagReasonDescs = new String[]
            {
                    "None",
                    "Coeluting peptide",
                    "Dissimilar KL",
                    "MS1 Ratio different from MS2",
                    "Light/Heavy have different high peak",
                    "Light/Heavy have different peak distributions",
            };

    protected int maxDaltonsUpCoelute = 3;
    protected int maxDaltonsDownCoelute = 3;





//scaffolding
    private List<Float> daltonsOff = new ArrayList<Float>();


    public FlagQuantEventsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "flagquant";
        mShortDescription = "Flag questionable quantitation events";
        mHelpMessage = "Flag questionable quantitation events";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "MS2 Feature files"),
//                        new DirectoryToReadArgumentDefinition("ms1featuresdir", true, "MS1 Feature directory"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "mzXML directory"),
                        new FileToWriteArgumentDefinition("out", false, "Output file"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
                       new EnumeratedValuesArgumentDefinition("label", false, labelStrings, labelExplanations,
                               "silac"),
                        new IntegerArgumentDefinition("maxdaltonsupcoelute", false,
                                "Number of Daltons up from the event that we care about coeluting features",
                                maxDaltonsUpCoelute),
                        new IntegerArgumentDefinition("maxdaltonsdowncoelute", false,
                                "Number of Daltons up from the event that we care about coeluting features",
                                maxDaltonsDownCoelute),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new DecimalArgumentDefinition("coeluteminintproportion", false,
                                "Minimum proportion of the quantitated feature's intensity (if that's determinable) " +
                                        "for a coeluting feature to have to cause an event to be flagged.",
                                minCoeluteIntensityProportion),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (featureFiles.length > 1)
        {
            assertArgumentPresent("outdir");
            assertArgumentAbsent("out");
        }

        if (outFile != null)
            assertArgumentAbsent("outdir");
        if (outDir != null)
            assertArgumentAbsent("out");

        mzXmlDir = getFileArgumentValue("mzxmldir");
//        ms1FeaturesDir = getFileArgumentValue("ms1featuresdir");

        labelType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label"));

        maxDaltonsUpCoelute = getIntegerArgumentValue("maxdaltonsupcoelute");
        maxDaltonsDownCoelute = getIntegerArgumentValue("maxdaltonsdowncoelute");

        minCoeluteIntensityProportion = getFloatArgumentValue("coeluteminintproportion");

        showCharts = getBooleanArgumentValue("showcharts");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File featureFile : featureFiles)
        {
            File outputFile = outFile;
            if (outputFile == null)
                outputFile = new File(outDir, featureFile.getName());

            Iterator<FeatureSet> featureSetIterator = null;
            try
            {
                featureSetIterator =
                        new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failure loading pepXML file " + featureFile.getAbsolutePath(),e);
            }

            int flaggedFeaturesThisFile = 0;
            int totalFeaturesThisFile = 0;

            List<File> tempFeatureFiles = new ArrayList<File>();
            int numSetsProcessed = 0;
            while (featureSetIterator.hasNext())
            {
                FeatureSet featureSet = featureSetIterator.next();
                totalFeaturesThisFile += featureSet.getFeatures().length;

//                try
//                {
//                    featureSet = new FeatureSet(featureFile);
//                }
//                catch (Exception e)
//                {
//                    throw new CommandLineModuleExecutionException("Failure processing file " +
//                            featureFile.getAbsolutePath(),e);
//                }

                MSRun run = null;
                try
                {
                    File featureSetFile = featureSet.getSourceFile();
                    if (MS2ExtraInfoDef.getFeatureSetBaseName(featureSet) != null)
                        featureSetFile = new File(MS2ExtraInfoDef.getFeatureSetBaseName(featureSet) + ".pep.xml");
                    File mzXmlFile = ViewerCommandModuleUtilities.findCorrespondingMzXmlFile(
                            featureSetFile, mzXmlDir);
                    ApplicationContext.infoMessage("Loading mzXml file " + mzXmlFile.getAbsolutePath());
                    run = MSRun.load(mzXmlFile.getAbsolutePath());
                    ApplicationContext.infoMessage("Loaded.");
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Can't find or open mzXml file for file " +
                            featureFile.getAbsolutePath(), e);
                }


                ApplicationContext.infoMessage("\tProcessing fraction " + (numSetsProcessed+1) + "...");

                String baseName = MS2ExtraInfoDef.getFeatureSetBaseName(featureSet);
                if (baseName == null)
                {
                    baseName = featureFile.getName();
                    if (numSetsProcessed > 0 || featureSetIterator.hasNext())
                        baseName = baseName + "_" + numSetsProcessed;
                }
                //if PeptideProphet was run from a directory below the directory containing the
                //mzXML files, we may have ../ in the baseName, which causes trouble in saving
                //the temporary files
                while (baseName.contains(".." + File.separator))
                    baseName.replaceFirst(".." + File.separator, "");

                File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);
//
//                FeatureSet ms1FeatureSet = null;
//                try
//                {
//                    ms1FeatureSet = ViewerCommandModuleUtilities.findCorrespondingFeatureSet(thisFractionFile, ms1FeaturesDir);
//                }
//                catch (IOException e)
//                {
//                    throw new CommandLineModuleExecutionException("Can't find or open MS1 feature file for file " +
//                            featureFile.getAbsolutePath(), e);
//                }

                processFeatureSet(featureSet, run);

                flaggedFeaturesThisFile += featureSet.getFeatures().length;

                try
                {
                    featureSet.savePepXml(thisFractionFile);
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failure writing file " +
                            thisFractionFile.getAbsolutePath(),e);
                }

                _log.debug("Saved fraction file as " + thisFractionFile.getAbsolutePath());
                tempFeatureFiles.add(thisFractionFile);

                numSetsProcessed++;
            }

            ApplicationContext.infoMessage("Flagged " + flaggedFeaturesThisFile + " out of " + totalFeaturesThisFile +
                       " features (" + ((float) flaggedFeaturesThisFile *100f / (float) totalFeaturesThisFile) +
                       "%) in this file");

            ApplicationContext.infoMessage("Saving output file " + outputFile.getAbsolutePath() + "...");

            try
            {
                if (numSetsProcessed == 1)
                {
                    FileReader in = new FileReader(tempFeatureFiles.get(0));
                    FileWriter out = new FileWriter(outputFile);
                    int c;

                    while ((c = in.read()) != -1)
                        out.write(c);

                    in.close();
                    out.close();
                }
                else
                {
                    ApplicationContext.infoMessage("\tCombining individual fraction files... " +
                            outputFile.getAbsolutePath() + "...");
                    new PepXMLFeatureFileHandler().combinePepXmlFiles(tempFeatureFiles, outputFile);
                }
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to save output file " + outputFile.getAbsolutePath(), e);
            }
            ApplicationContext.infoMessage("Done.");
        }
    }

    /**
     * Side effect: ms2FeatureSet retains only flagged features
     * @param ms2FeatureSet
     * @param run
     * @throws CommandLineModuleExecutionException
     */
    protected void processFeatureSet(FeatureSet ms2FeatureSet, MSRun run)
            throws CommandLineModuleExecutionException
    {
        List<Feature> flaggedFeatures = new ArrayList<Feature>();

        Window2DFeatureSetMatcher fsm = new Window2DFeatureSetMatcher();
        fsm.setMatchWithinChargeOnly(true);
        fsm.setMassDiffType(FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);
        fsm.setMinMassDiff(-maxDaltonsUpCoelute - massSlopDa);
        fsm.setMaxMassDiff(maxDaltonsDownCoelute + massSlopDa);

        fsm.setElutionMode(Window2DFeatureSetMatcher.ELUTION_MODE_SCAN);
        fsm.setElutionRangeMode(Window2DFeatureSetMatcher.ELUTION_RANGE_MODE_RANGE);
        fsm.setMinElutionDiff(-numScansSlopForCoeluting);
        fsm.setMaxElutionDiff(numScansSlopForCoeluting);

        Map<Feature, Feature> origOppositeFeatureMap = new HashMap<Feature, Feature>();
        List<Feature> origAndOppositeFeatures = new ArrayList<Feature>();
        Map<Feature, Boolean> featureIsLightMap = new HashMap<Feature, Boolean>();
        for (Feature feature : ms2FeatureSet.getFeatures())
        {
            if (!IsotopicLabelExtraInfoDef.hasRatio(feature))
                continue;
            Feature oppFeature = (Feature) feature.clone();
            boolean isLight = false;
            if (isHeavyLabeled(feature))
                oppFeature.setMass(feature.getMass() - calcHeavyLightMassDiff(feature));
            else if (isLightLabeled(feature))
            {
                isLight = true;
                oppFeature.setMass(feature.getMass() + calcHeavyLightMassDiff(feature));
            }
            else throw new RuntimeException("Neither light nor heavy labeled?!\n" + feature.toString());
            oppFeature.updateMz();

            featureIsLightMap.put(feature, isLight);
            origOppositeFeatureMap.put(feature, oppFeature);
            origAndOppositeFeatures.add(feature);
            origAndOppositeFeatures.add(oppFeature);
        }

        //WARNING! messing with ms2FeatureSet features here
        ms2FeatureSet.setFeatures(origAndOppositeFeatures.toArray(new Feature[origAndOppositeFeatures.size()]));
        List<Float> flaggedReasons = new ArrayList<Float>();


//        FeatureSetMatcher.FeatureMatchingResult matchingResult = fsm.matchFeatures(ms2FeatureSet, ms1FeatureSet);

//        Arrays.sort(ms2FeatureSet.getFeatures(), new Feature.ScanAscComparator());
//
//        List<Float> numsMatchedFeatures = new ArrayList<Float>();
//        List<Float> numsOppMatchedFeatures = new ArrayList<Float>();
//        List<Float> numFeaturesWithMatch = new ArrayList<Float>();
        for (Feature feature : origOppositeFeatureMap.keySet())
        {
            Feature oppFeature = origOppositeFeatureMap.get(feature);
//            List<Feature> matchedMs1Features = matchingResult.getSlaveSetFeatures(feature);
//            List<Feature> oppMatchedMs1Features = matchingResult.getSlaveSetFeatures(oppFeature);
//
//            numsMatchedFeatures.add(matchedMs1Features == null ? 0 : (float) matchedMs1Features.size());
//            numsOppMatchedFeatures.add(oppMatchedMs1Features == null ? 0 : (float) oppMatchedMs1Features.size());
//            numFeaturesWithMatch.add(0f + (matchedMs1Features == null ? 0 : 1) + (oppMatchedMs1Features == null ? 0 : 1));
//
            int flagReason = processFeature(feature, oppFeature,
                    featureIsLightMap.get(feature), run);
            flaggedReasons.add((float) flagReason);
            if (flagReason != FLAG_REASON_NONE)
            {
                feature.setDescription(flagReasonDescs[flagReason]);
                flaggedFeatures.add(feature);
            }
        }

        if (showCharts)
        {
//            new PanelWithHistogram(numsMatchedFeatures, "Orig matches").displayInTab();
//            new PanelWithHistogram(numsOppMatchedFeatures, "Opp matches").displayInTab();
//            new PanelWithHistogram(numFeaturesWithMatch, "num with match").displayInTab();
//            new PanelWithHistogram(daltonsOff, "Daltons off").displayInTab();

            new PanelWithHistogram(flaggedReasons, "flag reason").displayInTab();
        }
        ms2FeatureSet.setFeatures(flaggedFeatures.toArray(new Feature[flaggedFeatures.size()]));
    }

    protected int processFeature(Feature feature,// List<Feature> ms1MatchedFeatures,
                                 Feature oppFeature,// List<Feature> oppMs1MatchedFeatures,
                                 boolean baseIsLight,
                                 MSRun run)
    {



//List<Integer> badScans = Arrays.asList(new Integer[] {3424, 3516,8436,2304,3638,7493,8978,5545,2829,4879});
//Level origLevel = _log.getLevel();
//if (badScans.contains(feature.getScan()))
//{
//_log.setLevel(Level.DEBUG);
//_log.debug("***** scan " + feature.getScan());
//}
if (baseIsLight)
    _log.debug("Feature is light");
else
    _log.debug("Feature is heavy");
        int reason = FLAG_REASON_NONE;
//        Feature featureHighIntensityMs1Match = findHighestIntensityMs1Match(feature, ms1MatchedFeatures);
//        Feature oppFeatureHighIntensityMs1Match = findHighestIntensityMs1Match(oppFeature, oppMs1MatchedFeatures);
//        _log.debug("Has Feature? " + (featureHighIntensityMs1Match != null) +
//                ", oppFeature? " + (oppFeatureHighIntensityMs1Match != null));
//
////        if (checkCoelutingAboveBelow(feature, featureHighIntensityMs1Match, ms1MatchedFeatures) ||
////            checkCoelutingAboveBelow(oppFeature, oppFeatureHighIntensityMs1Match, oppMs1MatchedFeatures))
////            reason = FLAG_REASON_COELUTING_ABOVEBELOW;
//        if (reason == FLAG_REASON_NONE &&
//                checkDissimilarKL(featureHighIntensityMs1Match, oppFeatureHighIntensityMs1Match))
//            reason =  FLAG_REASON_DISSIMILAR_KL;
//        if (reason == FLAG_REASON_NONE &&
//                checkMs1RatioDifferent((float)IsotopicLabelExtraInfoDef.getRatio(feature),
//                featureHighIntensityMs1Match, oppFeatureHighIntensityMs1Match,
//                baseIsLight))
//            reason =  FLAG_REASON_DISSIMILAR_MS1_RATIO;
        if (reason == FLAG_REASON_NONE)
        {
//            _log.debug("Finding peaks...");
//            Feature[] peaksNearFeature = findPeaksNearFeature(feature, run);
//            Feature[] peaksNearOppFeature = findPeaksNearFeature(oppFeature, run);
//            _log.debug("\tnear: " + peaksNearFeature.length + ", opp: " + peaksNearOppFeature.length);

            List<Float> peakIntensities = calcPeakIntensities(feature, run);
            List<Float> oppPeakIntensities = calcPeakIntensities(oppFeature, run);

            float intensityBelow = peakIntensities.remove(0);
            float intensityBelowOpp = oppPeakIntensities.remove(0);


_log.debug("**" + peakIntensities.get(0) + ", " + peakIntensities.get(1) + ", " + peakIntensities.get(2) + ", " + peakIntensities.get(3));
_log.debug("**" + oppPeakIntensities.get(0) + ", " + oppPeakIntensities.get(1) + ", " + oppPeakIntensities.get(2) + ", " + oppPeakIntensities.get(3));            

_log.debug("************************" + calcKL(feature.getMass(), peakIntensities) + ", " + calcKL(oppFeature.getMass(), oppPeakIntensities));            


            float maxPeakInt = 0;
            int highPeakIndex = 0;
            float secondHighPeakInt = 0f;
            float maxPeakIntOpp = 0;
            int highPeakIndexOpp = 0;
            float secondHighPeakIntOpp = 0f;
            for (int i=0; i<peakIntensities.size(); i++)
            {
                if (peakIntensities.get(i) > maxPeakInt)
                {
                    if (maxPeakInt > secondHighPeakInt)
                        secondHighPeakInt = maxPeakInt;
                    maxPeakInt = peakIntensities.get(i);
                    highPeakIndex = i;
                }
                else if (peakIntensities.get(i) > secondHighPeakInt)
                    secondHighPeakInt = peakIntensities.get(i);
                if (oppPeakIntensities.get(i) > maxPeakIntOpp)
                {
                    if (maxPeakIntOpp > secondHighPeakIntOpp)
                        secondHighPeakIntOpp = maxPeakIntOpp;
                    maxPeakIntOpp = oppPeakIntensities.get(i);
                    highPeakIndexOpp = i;
                }
                else if (oppPeakIntensities.get(i) > secondHighPeakIntOpp)
                    secondHighPeakIntOpp = oppPeakIntensities.get(i);
            }


            float kl = calcKL(feature.getMass(), peakIntensities);
            float klOpp = calcKL(oppFeature.getMass(), oppPeakIntensities);
            float klDiff = Math.abs(kl - klOpp);
            float klRatio = kl / klOpp;
            if (kl > klOpp)
                klRatio = 1/klRatio;
_log.debug("KL: " + kl + ", klOPP: " + klOpp + ": diff=" + klDiff + ", ratio=" + klRatio);
            if (klDiff > 0.25 && klRatio < 0.7)
            {
                reason = FLAG_REASON_DIFF_PEAKDIST;
_log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
                return reason;
            }

            float maxPeakIntensity = BasicStatistics.max(peakIntensities);
            float maxPeakIntensityOpp = BasicStatistics.max(oppPeakIntensities);

            float belowIntensityRatio = intensityBelow / maxPeakIntensity;
            float belowIntensityRatioOpp = intensityBelowOpp / maxPeakIntensityOpp;

            _log.debug("BELOW_RATIO: " + belowIntensityRatio + ", OPP: " + belowIntensityRatioOpp);
            if (Math.max(belowIntensityRatio, belowIntensityRatioOpp) > 0.5)
            {
                reason = FLAG_REASON_COELUTING_ABOVEBELOW;
_log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
                return reason;
            }




//            float peakRatio = secondHighPeakInt / maxPeakInt;
//            float peakRatioOpp = secondHighPeakIntOpp / maxPeakIntOpp;



//            if (highPeakIndex != highPeakIndexOpp)
//            {
//                if ((peakRatio < 0.7) ||
//                        (peakRatioOpp < 0.7))
//                {
//_log.debug("Feature: high: " + maxPeakInt + ", second: " + secondHighPeakInt);
//_log.debug("Opp: high: " + maxPeakIntOpp + ", second: " + secondHighPeakIntOpp);
//
//_log.debug("DIFF HIGH: " + highPeakIndex + ", " + highPeakIndexOpp + ", ratios: " + (secondHighPeakInt / maxPeakInt) + ", " + (secondHighPeakIntOpp / maxPeakIntOpp));
//                    reason = FLAG_REASON_DIFF_HIGHESTPEAK;
//                }
//            }
//            else
//            {
//_log.debug("PEAKDIST: ratio1: " + peakRatio + ", Ratio2: " + peakRatioOpp);
//
//                if (Math.abs(Math.log(peakRatio) - Math.log(peakRatioOpp)) >=
//                    Math.log(0.9) - Math.log(0.6))
//                {
//_log.debug("DIFF PEAKDIST");
//                    reason = FLAG_REASON_DIFF_PEAKDIST;
//                }
//            }
        }
_log.debug("REASON: " + flagReasonDescs[reason]);
//_log.setLevel(origLevel);
        return reason;
    }

    protected float calcKL(float mass, List<Float> peakIntensities)
    {
        int signalLength = peakIntensities.size();
        float[] signal = new float[signalLength];
        int countPeaks = 0;
        float sum = 0;
        for (int i = 0; i < peakIntensities.size(); i++)
        {
            if (i < signal.length)
            {
                signal[i] = Math.max(0.1F, peakIntensities.get(i));
                sum += signal[i];
            }
        }
        for (int i = 0; i < signal.length; i++)
            signal[i] /= sum;

        return Spectrum.KLPoissonDistance(mass, signal);
    }


    /**
     * Return the list of peak intensities, starting with 1Da below mass
     * @param feature
     * @param run
     * @return
     */
    protected List<Float> calcPeakIntensities(Feature feature, MSRun run)
    {

        FloatRange mzRange = new FloatRange(
                feature.getMz() - ((1 * 1.000495f + 0.2f) / feature.getCharge()),
                feature.getMz() + ((6 * 1.000495f + 0.2f) / feature.getCharge()));
        //Resampling adds some slop on each end
        float minResampledMz = (int) mzRange.min;
        if (minResampledMz > mzRange.min) minResampledMz--;
        float maxResampledMz = (int) mzRange.max;
        if (maxResampledMz < mzRange.max) maxResampledMz++;
        FloatRange mzWindowResampled = new FloatRange(minResampledMz, maxResampledMz);

//        int numScansAround = 20;
//        int precursorIndex = run.getIndexForScanNum(feature.getScan(), true);
//System.err.println("Scans: " + IsotopicLabelExtraInfoDef.getHeavyFirstScan(feature) + ", " + IsotopicLabelExtraInfoDef.getHeavyLastScan(feature));
        int firstScanIndex = Math.max(Math.abs(run.getIndexForScanNum(IsotopicLabelExtraInfoDef.getHeavyFirstScan(feature))) - 2, 0);
        int lastScanIndex = Math.min(Math.abs(run.getIndexForScanNum(IsotopicLabelExtraInfoDef.getHeavyLastScan(feature))) + 2, run.getScanCount()-1);

//        int startScanIndex = Math.max(0, precursorIndex - numScansAround);
//        int endScanIndex = Math.max(precursorIndex + numScansAround, run.getScans().length-1);

        Scan[] scans = FeatureFinder.getScans(run, firstScanIndex,
                lastScanIndex - firstScanIndex + 1);
        _log.debug("Scans: " + scans.length);
//for (Scan scan : scans) _log.debug(scan.getNum());        
//System.err.println("SCANFIRST: " + firstScanIndex + ", SCANLAST: " + lastScanIndex + ", SCAN COUNT: " + (lastScanIndex - firstScanIndex + 1));
        //need fake mzRange?
        SpectrumResampler spectrumResampler = new SpectrumResampler(mzRange);
        int resolution = PanelWithSpectrumChart.DEFAULT_RESOLUTION;
        spectrumResampler.setResamplingFrequency(resolution);
        spectrumResampler.setUseMedianSmooth(false);
        float[][] resampledSpectra = null;
//System.err.println("Range resampled: " + mzRange.min + "-" + mzRange.max + ", buckets should be: " + ((mzRange.max -mzRange.min) * resolution));

        try
        {
            resampledSpectra =
                    spectrumResampler.resampleSpectra(scans);
//System.err.println("Actual buckets: " + resampledSpectra[0].length + ", max mz should be " + (mzRange.min + (resampledSpectra[0].length / resolution)));            
        }
        catch (InterruptedException e)
        {}



        List<Float> peakIntensities = new ArrayList<Float>();
        for (int i=-1; i<4; i++)
        {
            float intensitySum = 0;
            float mz = feature.getMz() + ((float) i * 1.000495f / (float) feature.getCharge());
            float mzTol = (MassUtilities.calculateAbsoluteDeltaMass(
                    feature.getMass(), 50, FeatureSetMatcher.DELTA_MASS_TYPE_PPM)) / feature.getCharge();
            float minPeakMz = mz - mzTol;
            float maxPeakMz = mz + mzTol;
            float lowBucket = (minPeakMz - mzRange.min) * resolution;
            int lowBucketIndex = (int)Math.floor(lowBucket);
            float highBucket = (maxPeakMz - mzRange.min) * resolution;
            int highBucketIndex = (int)Math.floor(highBucket);
//System.err.println(i + ", " + mz + ", " + minPeakMz + ", " + maxPeakMz + ", lowb=" + lowBucket + ", highb=" + highBucket + ", buckets=" + resampledSpectra[0].length +  ", highmzrecalc=" + (mzRange.min + (highBucketIndex / (float) resolution)));

            for (int scanIndex = 0; scanIndex < resampledSpectra.length; scanIndex++)
            {
                for (int bucket = lowBucketIndex; bucket <= highBucketIndex; bucket++)
                    intensitySum += resampledSpectra[scanIndex][bucket];
            }
//_log.debug(i + ", " + mz + ", " + minPeakMz + ", " + maxPeakMz + ", lowb=" + lowBucket + ", highb=" + highBucket  + ", int=" + intensitySum);

            peakIntensities.add(intensitySum);
        }
        return peakIntensities;
    }

    protected Feature[] findPeaksNearFeature(Feature feature, MSRun run)
    {
            WaveletPeakExtractor peakExtractor = new WaveletPeakExtractor();
            FloatRange mzRange = new FloatRange(
                    feature.getMz() - ((maxDaltonsDownCoelute + massSlopDa) / feature.getCharge()),
                    feature.getMz() + ((maxDaltonsUpCoelute + massSlopDa) / feature.getCharge()));

        int numScansAround = 20;
        int precursorIndex = run.getIndexForScanNum(feature.getScan(), true);
        int startScanIndex = Math.max(0, precursorIndex - numScansAround);
        int scanCountWithWindow = Math.max(precursorIndex + numScansAround, run.getScans().length-1);

        Scan[] scans = FeatureFinder.getScans(run, startScanIndex,
                scanCountWithWindow);
        //need fake mzRange?
        SpectrumResampler spectrumResampler = new SpectrumResampler(mzRange);
        spectrumResampler.setResamplingFrequency(FeatureStrategy.DEFAULT_RESAMPLING_FREQUENCY);
        spectrumResampler.setUseMedianSmooth(false);
        Feature[] peaks = null;
        try
        {
            float[][] resampledSpectra =
                    spectrumResampler.resampleSpectra(scans);
            peaks = peakExtractor.extractPeakFeatures(scans, resampledSpectra, mzRange);

        }
        catch (InterruptedException e)
        {}
//        FeatureSet.FeatureSelector intensitySel = new FeatureSet.FeatureSelector();
//        intensitySel.setMinIntensity(feature.getIntensity() / 20f);
//        peaks = new FeatureSet(peaks).filter(intensitySel).getFeatures();
        return peaks;
    }

    protected boolean checkDissimilarKL(Feature ms1BestMatchFeature, Feature oppMs1BestMatchFeature)
    {
        if (ms1BestMatchFeature == null || oppMs1BestMatchFeature == null)
            return false;
        float kl1 = ms1BestMatchFeature.kl;
        float kl2 = oppMs1BestMatchFeature.kl;
        _log.debug("\tKL1: " + kl1 + ", KL2: " + kl2 + ", ratio: " + (kl1/kl2));
        //can pick up more by lowering the abs difference below 0.2, but lotsa false pos
        if (Math.abs(kl1 - kl2) > 1 &&
            ((kl1 / kl2) > 1.5 || (kl1 / kl2) < (.67)))
        {
            _log.debug("\t***Bad KL: " + kl1 + ", " + kl2);
            return true;
        }
        return false;
    }

    protected boolean checkMs1RatioDifferent(float ms2Ratio, Feature ms1BestMatchFeature, Feature oppMs1BestMatchFeature,
                                             boolean baseFeatureIsLight)
    {
        if (ms1BestMatchFeature == null || oppMs1BestMatchFeature == null)
            return false;
        float ms1Ratio = ms1BestMatchFeature.intensity / oppMs1BestMatchFeature.intensity;
        if (!baseFeatureIsLight)
            ms1Ratio = 1.0f / ms1Ratio;
        float ratioDiff = (float) Math.abs(Math.log(ms1Ratio) - Math.log(ms2Ratio));
        _log.debug("\tRatios: MS2: " + ms2Ratio + ", MS1: " + ms1Ratio);

        //we care more, and are willing to accept less of a difference, if one ratio goes one way and the other goes the other
        if (ratioDiff > Math.log(2) - Math.log(1.0) ||
            (((ms1Ratio > 1 && ms2Ratio < 1) || (ms1Ratio < 1 && ms2Ratio > 1)) && ratioDiff > Math.log(1.3) - Math.log(1.0)))
        {
            _log.debug("\t***Bad ratio.  MS2: " + ms2Ratio + ", MS1: " + ms1Ratio);
            return true;
        }
        return false;
    }

    /**
     * Null if nothing matches
     * @param ms2Feature
     * @param ms1MatchedFeatures
     * @return
     */
    protected Feature findHighestIntensityMs1Match(Feature ms2Feature, List<Feature> ms1MatchedFeatures)
    {
        if (ms1MatchedFeatures == null)
            return null;
        Feature highestIntensityMs1TightMatch = null;
        for (Feature matchedFeature : ms1MatchedFeatures)
        {
            float massDiffPPM =
                    MassUtilities.calculatePPMDeltaMass(ms2Feature.getMass(), matchedFeature.getMass() - ms2Feature.getMass(),
                    FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);
            if (massDiffPPM < massSlopPPM)
            {
                if (highestIntensityMs1TightMatch == null ||
                        matchedFeature.getIntensity() > highestIntensityMs1TightMatch.getIntensity())
                    highestIntensityMs1TightMatch = matchedFeature;
            }
        }
        return highestIntensityMs1TightMatch;
    }

    /**
     * Returns true if there's a feature above or below, with intensity greater than
     * minCoeluteIntensityProportion times the intensity of the event feature itself (if there is one) 
     *
     * @param feature
     * @param highIntensityMs1Match the highest-intensity feature that matches feature within massSlopPPM
     * @param ms1MatchedFeatures
     * @return
     */
    protected boolean checkCoelutingAboveBelow(Feature feature, Feature highIntensityMs1Match,
                                               List<Feature> ms1MatchedFeatures)
    {
        if (ms1MatchedFeatures == null)
            return false;
        float intensityToBeat = 0f;
        float thisFeatureIntensity = 0f;

        if (highIntensityMs1Match != null)
            thisFeatureIntensity = highIntensityMs1Match.getIntensity();
        intensityToBeat = minCoeluteIntensityProportion * thisFeatureIntensity;
        _log.debug("Intensity to beat: " + intensityToBeat +", feature int: " + thisFeatureIntensity);
        for (Feature matchedFeature : ms1MatchedFeatures)
        {
            float massDiff = matchedFeature.getMass() - feature.getMass();
            float fractionalMassDiff = massDiff - (int) massDiff;//(massDiff + 0.5f) % 1.0f - 0.5f;
            float fractionaMassDiffPPM = MassUtilities.calculatePPMDeltaMass(matchedFeature.getMass(), fractionalMassDiff, FeatureSetMatcher.DELTA_MASS_TYPE_ABSOLUTE);

            if (Math.abs(massDiff) < 0.1f)
            {
                if (Math.abs(massDiff) < massSlopPPM)
                {
                    _log.debug("Identity feature, moving along....");
                    continue;

                }
                else if (matchedFeature.getIntensity() > intensityToBeat)
                {
                    _log.debug("Coeluting non-identity feature!");
                    return true;
                }
                else
                {
                    _log.debug("Too low intensity, moving along....");
                    continue;
                }
            }
            _log.debug("\tcoelut check, massdiff=" + massDiff + ", match int: " + matchedFeature.getIntensity());

            if (Math.abs(massDiff) > 0.5f && Math.abs(fractionaMassDiffPPM) < massSlopPPM  &&
                    matchedFeature.getIntensity() >= intensityToBeat)
            {
                daltonsOff.add((float) ((int) massDiff));
                _log.debug("\tCoeluting, massdiff=" + massDiff);

                return true;
            }
        }

        return false;
    }


    /**
     * Returns false if the peptide does not contain the residue
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    protected boolean checkForModAllResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return false;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                boolean foundIt = false;
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    return false;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                    {
                        foundIt = true;
                        break;
                    }
                if (!foundIt)
                    return false;
            }
        }
        return true;
    }

    protected float calcHeavyLightMassDiff(Feature feature)
    {
        String residue = "";
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                residue = "C";
                break;
            case LABEL_LYCINE:
                residue = "K";
                break;
        }
        int numLabels = 0;
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
        for (int i=0; i<peptide.length(); i++)
        {
            if (residue.equals("" + peptide.charAt(i)))
                numLabels++;
        }

        return numLabels * (labelType == LABEL_ACRYLAMIDE ? ACRYLAMIDE_LABEL_MASSDIFF_PERRESIDUE :
                SILAC_LABEL_MASSDIFF_PERRESIDUE);
    }


    /**
     * Returns false if the peptide does not contain the residue
     * @param peptide
     * @param mods
     * @param residue
     * @param modMass
     * @return
     */
    protected boolean checkForModNoResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {
        if (!peptide.contains("" + residue))
            return false;
        if (mods == null)
            return true;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == residue)
            {
                List<ModifiedAminoAcid> modsThisResidue = mods[i];
                if (modsThisResidue == null)
                    continue;
                for (ModifiedAminoAcid mod : modsThisResidue)
                    if (Math.abs(mod.getMass() - modMass) < MOD_MASS_SLOP)
                        return false;
            }
        }
        return true;
    }


    protected boolean isLightLabeled(Feature feature)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean lightLabeled = false;
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                lightLabeled =  checkForModAllResidues(peptide, mods, 'C', ACRYLAMIDE_LABEL_LIGHTMASS);
                break;
            case LABEL_LYCINE:
                lightLabeled =  checkForModNoResidues(peptide, mods, 'K', SILAC_LABEL_MASS);
                break;
        }
        return lightLabeled;
    }

    protected boolean isHeavyLabeled(Feature feature)
    {
        List<ModifiedAminoAcid>[] mods = MS2ExtraInfoDef.getModifiedAminoAcids(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

        boolean heavyLabeled = false;
        switch(labelType)
        {
            case LABEL_ACRYLAMIDE:
                heavyLabeled = checkForModAllResidues(peptide, mods, 'C', ACRYLAMIDE_LABEL_HEAVYMASS);
                break;
            case LABEL_LYCINE:
                heavyLabeled = checkForModAllResidues(peptide, mods, 'K', SILAC_LABEL_MASS);
                break;
        }
        return heavyLabeled;
    }
}
