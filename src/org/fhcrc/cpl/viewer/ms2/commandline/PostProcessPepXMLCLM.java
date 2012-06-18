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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;
import org.apache.commons.lang.StringUtils;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.QuantitationUtilities;

import java.io.*;
import java.util.*;


/**
 * post-process a pepxml file
 */
public class PostProcessPepXMLCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PostProcessPepXMLCLM.class);

    protected File[] pepXmlFiles;

    protected boolean medianCenter = false;
    protected boolean medianCenterByNumCysteines = false;
    protected boolean medianCenterAllRunsTogether = false;
    protected float medianRatioForCentering = -1;

    protected boolean stripQuantMissingLightOrHeavyWithinRun = false;
    protected boolean stripQuantMissingLightOrHeavyAcrossAll = false;
    protected boolean stripQuantNotInHeavyAcrossAll = false;
    protected boolean stripLightIDs = false;                                            
    protected boolean adjustQuantZeroAreas = false;
    protected boolean stripQuantZeroAreas = false;
    protected boolean stripQuantSingleScans = false;
    protected boolean stripQuantMatchingRegexp = false;

    protected String regexp;

    protected Map<File, Float> fileMedianLogRatioMap;
    protected Map<File, Map<Integer, Float>> fileNumCysteinesMedianLogRatioMap;

    protected boolean filterByProteinPrefix = false;

    protected int percentileForQuantZeroAreaAdjustment = 1;

    protected String badProteinPrefix;
    protected String goodProteinPrefix;
    protected boolean excludeProteinPrefixQuantOnly;

    protected File outFile;
    protected File outDir;

    protected Set<String> peptidesToStrip = null;
    protected Set<String> proteinsToStrip = null;
    protected Set<String> proteinsToStripQuant = null;


    protected Set<String> proteinsToKeep = null;


    protected boolean showCharts = false;

    protected int minRatiosForMedianCenter = 10;

    protected float minPeptideProphetForMedian = .75f;

    protected float minPeptideProphet = 0f;
    protected float minQuantPeptideProphet = 0f;

    protected float maxFDR = Float.MAX_VALUE;


    protected float maxExpect = Float.MAX_VALUE;
    protected float maxQuantExpect = Float.MAX_VALUE;

    protected float minRatio = 0f;
    protected float maxRatio = Float.MAX_VALUE;


    protected boolean requirePepXmlExtension = false;

    protected DeltaMassArgumentDefinition.DeltaMassWithType maxFracDeltaMass = null;




    protected int labelType = QuantitationUtilities.LABEL_ACRYLAMIDE;


    //for filtering out peptides not seen in both light and heavy in some run
    protected Set<String> lightPeptidesAllRuns = new HashSet<String>();
    protected Set<String> heavyPeptidesAllRuns = new HashSet<String>();




    protected String[] labelStrings = new String[]
            {
                    "acrylamide",
                    "silac",
                    "silacrk"
            };

    protected String[] labelExplanations = new String[]
            {
                    "Acrylamide (3.0106Da on C)",
                    "SILAC Lycine labeling (134.115092 on K, i.e. 6Da SILAC)",
                    "SILAC Lycine and Arginine labeling (6Da SILAC)"
            };


    public PostProcessPepXMLCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "postprocesspepxml";

        mHelpMessage ="Post-process PepXML.  This provides tools for stripping out peptides and median-centering log ratios.";
        mShortDescription = "Post-process PepXML.  This provides tools for stripping out peptides and median-centering log ratios.";

        CommandLineArgumentDefinition[] argDefs =
               {
                       this.createUnnamedSeriesFileArgumentDefinition(true, "PepXML files to process"),
                       new BooleanArgumentDefinition("mediancenter", false, "Median-center ratios?", medianCenter),
                       new BooleanArgumentDefinition("bynumcysteines", false,
                               "Median-center ratios separately by number of Cysteines?", medianCenter),
                       new BooleanArgumentDefinition("mediancenterallrunstogether", false,
                               "Median-center ratios all runs together, across all files?  If false, median-centers " +
                               "separately for each file (all fractions together)", medianCenterAllRunsTogether),
                       new FileToReadArgumentDefinition("strippeptidefile", false,
                               "File containing a list of peptides to strip from results, one per line, all caps"),
                       new FileToReadArgumentDefinition("stripproteinfile", false,
                               "File containing a list of protein identifiers to strip from results, one per line"),
                       new FileToReadArgumentDefinition("stripproteinquantfile", false,
                               "File containing a list of protein identifiers to strip /quantitation/ from results, one per line"),
                       new FileToReadArgumentDefinition("keepproteinfile", false,
                               "File containing a list of protein identifiers to keep in results (strip all others), one per line"),
                       new FileToWriteArgumentDefinition("out", false, "Output file"),
                       new DirectoryToWriteArgumentDefinition("outdir", false, "Output directory"),
                       new BooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                       new IntegerArgumentDefinition("minmediancenterratios", false,
                               "Minimum number of ratios necessary in order to perform median-centering",
                               minRatiosForMedianCenter),
                       new DecimalArgumentDefinition("minmedianpprophet", false,
                               "Minimum PeptideProphet score to be counted in median calculation",
                               minPeptideProphetForMedian),
                       new BooleanArgumentDefinition("stripquantmissingheavy", false,
                               "Strip quantitation events in which the heavy " +
                               "isotope was never identified, in any run",
                               stripQuantNotInHeavyAcrossAll),
                       new BooleanArgumentDefinition("stripquantmissinglightorheavywithinrun", false,
                               "Strip peptides that we haven't seen in both light and heavy states, " +
                               "within a single run",
                               stripQuantMissingLightOrHeavyWithinRun),
                       new BooleanArgumentDefinition("stripquantmissinglightorheavyacrossruns", false,
                               "Strip peptides that we haven't seen in both light and heavy states, " +
                               "across all runs.  This ONLY makes sense for multiple metabolically-labeled " +
                               "experiments with a label flip",
                               stripQuantMissingLightOrHeavyAcrossAll),
                       new EnumeratedValuesArgumentDefinition("label", false, labelStrings, labelExplanations,
                               "acrylamide"),
                       new BooleanArgumentDefinition("filterbyproteinprefix", false,
                               "Filter peptides based on prefixes of the protein names that they're associated with?",
                               filterByProteinPrefix),
                       new StringArgumentDefinition("badproteinprefix",false,
                               "Exclude any peptides with any associated proteins with this prefix to their names"),
                       new StringArgumentDefinition("goodproteinprefix",false,
                               "Include any peptides with any associated proteins with this prefix to their names"),
                       new BooleanArgumentDefinition("protprefixexcludequantonly", false,
                               "When excluding peptides based on protein prefix, exclude only quantitation?  " +
                                       "If false, excludes entire ID",
                               excludeProteinPrefixQuantOnly),
                       new BooleanArgumentDefinition("requirepepxmlextension", false,
                               "When looking for files in a pepxmldir, require that they end with .pep.xml?",
                               requirePepXmlExtension),
                       new DecimalArgumentDefinition("minpprophet", false,
                               "Minimum PeptideProphet score to keep", minPeptideProphet),
                       new DecimalArgumentDefinition("maxfdr", false,
                               "Maximum FDR (determined using PeptideProphet probabilities, " +
                                       "separately for each file) to keep", maxFDR),
                       new DecimalArgumentDefinition("minquantpprophet", false,
                               "Minimum PeptideProphet score for quantitation", minQuantPeptideProphet),
                       new DecimalArgumentDefinition("maxexpect", false,
                               "Maximum expect score to keep", maxExpect),
                       new DecimalArgumentDefinition("maxquantexpect", false,
                               "Maximum expect score for quantitation", maxQuantExpect),
                       new DeltaMassArgumentDefinition("maxfracdeltamass", false,
                               "Maximum fractional delta mass"),
                       new BooleanArgumentDefinition("adjustquantzeroareas", false,
                               "Adjust zero values for light or heavy areas in quantitation (and ratios) to the " +
                               percentileForQuantZeroAreaAdjustment + " percentile of all the (nonzero) values",
                               adjustQuantZeroAreas),
                       new BooleanArgumentDefinition("stripquantzeroareas", false,
                               "Strip quantitation with zero values for light or heavy areas",
                               stripQuantZeroAreas),
                       new BooleanArgumentDefinition("stripquantsinglescan", false,
                               "Strip quantitation with light OR heavy scan extents of just one scan",
                               stripQuantSingleScans),
                       new BooleanArgumentDefinition("striplightids", false,
                               "Strip all light-labeled IDs (for aminoacid labels, light on all residues), regardless " +
                                       "of whether there is a corresponding heavy ID",
                               stripLightIDs),
                       new BooleanArgumentDefinition("stripquantmatchingregexp", false,
                               "Strip quantitation from all peptides whose sequences match regexp", stripQuantMatchingRegexp),
                       new StringArgumentDefinition("regexp",false, "Regular expression"),
                       new DecimalArgumentDefinition("minratio", false,
                               "Minimum quantitative ratio to keep", minRatio),
                       new DecimalArgumentDefinition("maxratio", false,
                               "Maximum quantitative ratio to keep", maxRatio),
                       new DecimalArgumentDefinition("medianratioforcenter", false,
                               "Median ratio to use explicitly when median-centering " +
                                       "(i.e., median will not be calculated from data)", medianRatioForCentering),


               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        pepXmlFiles = getUnnamedSeriesFileArgumentValues();

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (pepXmlFiles.length > 1)
        {
            assertArgumentAbsent("out", "(unnamed)");
            assertArgumentPresent("outdir", "(unnamed)");
        }

        if (hasArgumentValue("maxfdr"))
            assertArgumentAbsent("minpprophet","maxfdr");

        medianCenter = getBooleanArgumentValue("mediancenter");
        medianCenterByNumCysteines = getBooleanArgumentValue("bynumcysteines");
        medianCenterAllRunsTogether = getBooleanArgumentValue("mediancenterallrunstogether");
        medianRatioForCentering = getFloatArgumentValue("medianratioforcenter");
        if (hasArgumentValue("medianratioforcenter") && (medianCenterByNumCysteines || medianCenterAllRunsTogether ||
            !medianCenter)) {
            throw new ArgumentValidationException("medianratioforcenter was provided along with incompatible args");
        }

        requirePepXmlExtension = getBooleanArgumentValue("requirepepxmlextension");

        stripQuantMissingLightOrHeavyWithinRun = getBooleanArgumentValue("stripquantmissinglightorheavywithinrun");
        stripQuantMissingLightOrHeavyAcrossAll = getBooleanArgumentValue("stripquantmissinglightorheavyacrossruns");
        stripQuantNotInHeavyAcrossAll = getBooleanArgumentValue("stripquantmissingheavy");
        stripLightIDs = getBooleanArgumentValue("striplightids");
        if (stripLightIDs)
            assertArgumentPresent("label", "striplightids");

        stripQuantMatchingRegexp = getBooleanArgumentValue("stripquantmatchingregexp");
        if (stripQuantMatchingRegexp)
            assertArgumentPresent("regexp", "stripquantmatchingregexp");
        regexp = getStringArgumentValue("regexp");

        int numIncompatibleFlags = 0;
        if (stripQuantMissingLightOrHeavyWithinRun) numIncompatibleFlags++;
        if (stripQuantMissingLightOrHeavyAcrossAll) numIncompatibleFlags++;
        if (stripQuantNotInHeavyAcrossAll) numIncompatibleFlags++;

        if (numIncompatibleFlags > 1)
        {
            throw new ArgumentValidationException("Only one of these flags may be specified: " +
                    "'stripquantmissinglightorheavywithinrun', 'stripquantmissinglightorheavyacrossruns', " +
                    "'stripquantmissingheavy'");
        }

        adjustQuantZeroAreas = getBooleanArgumentValue("adjustquantzeroareas");
        stripQuantZeroAreas = getBooleanArgumentValue("stripquantzeroareas");
        if (adjustQuantZeroAreas && stripQuantZeroAreas)
            throw new ArgumentValidationException("Can't both adjust /and/ strip zero areas!");

        stripQuantSingleScans = getBooleanArgumentValue("stripquantsinglescan");

        filterByProteinPrefix = getBooleanArgumentValue("filterbyproteinprefix");

        if (filterByProteinPrefix && !hasArgumentValue("badproteinprefix") && !hasArgumentValue("goodproteinprefix"))
        {                                                                        
            throw new ArgumentValidationException("If filtering by protein prefix, must specify either a bad or good protein prefix");
        }

        badProteinPrefix = getStringArgumentValue("badproteinprefix");
        goodProteinPrefix = getStringArgumentValue("goodproteinprefix");
        excludeProteinPrefixQuantOnly = getBooleanArgumentValue("protprefixexcludequantonly");

        if (badProteinPrefix != null && goodProteinPrefix != null)
            throw new ArgumentValidationException("Can't have it both ways (good and bad protein prefixes), sorry");

        if (!medianCenter) {
            assertArgumentAbsent("bynumcysteines","mediancenter");
            assertArgumentAbsent("medianratioforcenter", "mediancenter");
        }

        labelType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label"));

        if (hasArgumentValue("stripproteinfile"))
        {
            assertArgumentAbsent("keepproteinfile", "stripproteinfile");
            proteinsToStrip = new HashSet<String>(readOneStringPerLine(getFileArgumentValue("stripproteinfile")));
        }

        if (hasArgumentValue("stripproteinquantfile"))
        {
            proteinsToStripQuant = new HashSet<String>(readOneStringPerLine(getFileArgumentValue("stripproteinquantfile")));
        }

        if (hasArgumentValue("keepproteinfile"))
        {
            assertArgumentAbsent("stripproteinfile", "keepproteinfile");
            proteinsToKeep = new HashSet<String>(readOneStringPerLine(getFileArgumentValue("keepproteinfile")));
        }

        if (hasArgumentValue("strippeptidefile"))
        {
            Set<String> rawPeptides = new HashSet<String>(readOneStringPerLine(getFileArgumentValue("strippeptidefile")));
            peptidesToStrip = new HashSet<String>();
            //special processing for peptides.  Make sure they're actually peptides
            for (String peptide : rawPeptides)
            {
                boolean ok = true;
                if (peptide.length() == 0)
                    ok = false;
                for (int i = 0; i < peptide.length(); i++)
                {
                    char c = peptide.charAt(i);
                    if ((c >= 'A') && (c <= 'Z')) continue; // uppercase
                    ok = false;
                }
                if (!ok)
                    throw new RuntimeException();
                peptidesToStrip.add(peptide);
            }
        }



        minRatiosForMedianCenter = getIntegerArgumentValue("minmediancenterratios");

        showCharts = getBooleanArgumentValue("showcharts");

        minPeptideProphetForMedian = getFloatArgumentValue("minmedianpprophet");
        minPeptideProphet = getFloatArgumentValue("minpprophet");
        minQuantPeptideProphet = getFloatArgumentValue("minquantpprophet");

        maxFDR = getFloatArgumentValue("maxfdr");

        minRatio = getFloatArgumentValue("minratio");
        maxRatio = getFloatArgumentValue("maxratio");


        maxExpect = getFloatArgumentValue("maxexpect");
        maxQuantExpect = getFloatArgumentValue("maxquantexpect");

        if (hasArgumentValue("maxfracdeltamass"))
            maxFracDeltaMass = getDeltaMassArgumentValue("maxfracdeltamass");

        if (peptidesToStrip == null && proteinsToStrip == null && proteinsToStripQuant == null && proteinsToKeep == null &&
                !medianCenter && !stripQuantMissingLightOrHeavyWithinRun && !stripQuantMissingLightOrHeavyAcrossAll &&
                !filterByProteinPrefix &&
                !stripQuantNotInHeavyAcrossAll && !adjustQuantZeroAreas && !stripQuantZeroAreas &&
                !stripQuantSingleScans &&
                !hasArgumentValue("maxexpect") && minPeptideProphet == 0 && minQuantPeptideProphet == 0 &&
                !stripLightIDs && !stripQuantMatchingRegexp && minRatio == 0 && maxRatio == Float.MAX_VALUE &&
                maxFDR == Float.MAX_VALUE)
        {
            throw new ArgumentValidationException("Nothing to do!  Quitting");
        }

    }

    /**
     * Read each line as a String, stopping at first whitespace
     * @param file
     * @return
     * @throws ArgumentValidationException
     */
    protected List<String> readOneStringPerLine(File file)
            throws ArgumentValidationException
    {
        List<String> result =  new ArrayList<String>();
        try
        {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = null;
            while ((line = br.readLine()) != null)
            {
                if (line.startsWith("#"))
                    continue;
                String protein = StringUtils.strip(line);
                //if there's more than one column in the file, take first
                protein = protein.replaceFirst("\\s.*", "");

                result.add(protein);
            }

        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Failed to retrieve list from file " +
                    file + ".  Please make sure file contains only a list of identifiers, one per line");
        }
        return result;
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {

        if (stripQuantNotInHeavyAcrossAll || stripQuantMissingLightOrHeavyAcrossAll)
            loadLightHeavyPeptidesAcrossAll();
        if (medianCenter && medianRatioForCentering == -1)
            calcLogMedianRatiosAllFiles();

        for (File file : pepXmlFiles)
        {
            try
            {
                File outputFile = outFile;
                if (outFile == null)
                    outputFile = new File(outDir, calcOutputFilename(file.getName()));
                ApplicationContext.infoMessage("Processing file " + file.getAbsolutePath() + ", output file " +
                        outputFile.getAbsolutePath());
                handleFeatureFile(file, outputFile);
            }
            catch (Exception e)
            {
                ApplicationContext.setMessage("WARNING: Failed to process file " + file.getAbsolutePath() + " as a PepXML file");
            }
        }
        ApplicationContext.setMessage("Done.");
    }

    /**
     * TODO: fold this in with loadLightHeavyPeptides
     *
     * This is some really weird stuff, right here.  This method handles calculation of median log ratios, regardless
     * of whether you're doing it separately by file or all together, by number of cysteines or not.
     * It creates a map data structure that gets used down below, with a key for each file.  If we're doing global
     * centering, the value gets repeated for each key.  This is a bit silly, but it makes things simpler below.
     *
     * @return
     * @throws CommandLineModuleExecutionException
     */
    protected void calcLogMedianRatiosAllFiles()
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Calculating median log ratio(s) for all files, " +
                "this may take a while...");

        if (medianCenterByNumCysteines)
        {
            fileNumCysteinesMedianLogRatioMap = new HashMap<File, Map<Integer, Float>>();
            Map<Integer, List<Float>> numCysLogRatiosMapAllFiles = new HashMap<Integer, List<Float>>();
            for (File featureFile : pepXmlFiles)
            {
                Map<Integer, List<Float>> numCysLogRatiosMapThisFile = new HashMap<Integer, List<Float>>();
                ApplicationContext.infoMessage("\tProcessing file " + featureFile.getAbsolutePath() + "...");
                try
                {
                    Iterator<FeatureSet> featureSetIterator =
                            new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);
                    while (featureSetIterator.hasNext())
                    {
                        FeatureSet featureSet = featureSetIterator.next();
                        filterOnQualityScores(featureSet);

                        for (Feature feature : featureSet.getFeatures())
                        {
                            if (IsotopicLabelExtraInfoDef.hasRatio(feature) &&
                                    MS2ExtraInfoDef.getPeptideProphet(feature) >= minPeptideProphetForMedian)
                            {
                                float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                                if (Float.isInfinite(ratio) || (ratio == 0) || Float.isNaN(ratio))
                                    continue;

                                //don't consider to-be-stripped proteins in median ratio calculation
                                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                                boolean hasBadProtein = false;
                                if (proteins != null)
                                {
                                    for (String protein : proteins)
                                    {
                                        if ((proteinsToStripQuant != null && proteinsToStripQuant.contains(protein)) ||
                                            (proteinsToStrip != null &&  proteinsToStrip.contains(protein)))
                                        {
                                            hasBadProtein = true;
                                            break;
                                        }
                                    }
                                }
                                if (hasBadProtein)
                                    continue;

                                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

                                if (peptide != null)
                                {
                                if (peptidesToStrip != null && peptidesToStrip.contains(peptide))
                                    continue;
                                    int numCysteines = 0;
                                    for (int i=0; i<peptide.length(); i++)
                                    {
                                        if (peptide.charAt(i) == 'C')
                                            numCysteines++;
                                    }
                                    if (numCysteines > 0)
                                    {
                                        List<Float> logRatiosList = numCysLogRatiosMapThisFile.get(numCysteines);
                                        if (logRatiosList == null)
                                        {
                                            logRatiosList = new ArrayList<Float>();
                                            numCysLogRatiosMapThisFile.put(numCysteines, logRatiosList);
                                        }
                                        logRatiosList.add((float) Math.log(ratio));
                                    }
                                }
                            }
                        }
                    }

                    if (medianCenterAllRunsTogether)
                    {
                        for (int numCys : numCysLogRatiosMapThisFile.keySet())
                        {
                            List<Float> logRatiosThisNumCys = numCysLogRatiosMapAllFiles.get(numCys);
                            if (logRatiosThisNumCys == null)
                            {
                                logRatiosThisNumCys = new ArrayList<Float>();
                                numCysLogRatiosMapAllFiles.put(numCys, logRatiosThisNumCys);
                            }
                            logRatiosThisNumCys.addAll(numCysLogRatiosMapThisFile.get(numCys));
                        }
                    }
                    else
                    {
                        HashMap<Integer, Float> numCysteinesMedianThisFile =
                                new HashMap<Integer, Float>();
                        for (int numCys : numCysLogRatiosMapThisFile.keySet())
                        {
                            numCysteinesMedianThisFile.put(numCys,
                                    (float)BasicStatistics.median(numCysLogRatiosMapThisFile.get(numCys)));
                        }
                        fileNumCysteinesMedianLogRatioMap.put(featureFile, numCysteinesMedianThisFile);
                    }

                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failed to load feature file " + featureFile);
                }
            }

            if (medianCenterAllRunsTogether)
            {
                Map<Integer, Float> numCysteinesMedianMap = new HashMap<Integer, Float>();
                for (int numCys = 0; numCys<10; numCys++)
                {
                    if (!numCysLogRatiosMapAllFiles.containsKey(numCys))
                        continue;
                    numCysteinesMedianMap.put(numCys, (float)BasicStatistics.median(numCysLogRatiosMapAllFiles.get(numCys)));
                    if (showCharts)
                    {
                        new PanelWithHistogram(numCysLogRatiosMapAllFiles.get(numCys),
                                "LogRatiosCys" + numCys, 200).displayInTab();
                    }
                }
                for (File file : pepXmlFiles)
                {
                    fileNumCysteinesMedianLogRatioMap.put(file, numCysteinesMedianMap);
                }
                ApplicationContext.infoMessage("Median log ratios by num Cysteines:");
                for (int i=0; i<20; i++)
                {
                    if (numCysteinesMedianMap.containsKey(i))
                    {
                        ApplicationContext.infoMessage(i + ": " + numCysteinesMedianMap.get(i) + " (" + numCysLogRatiosMapAllFiles.get(i).size() + " events)");
                    }
                }
            }
            else ApplicationContext.infoMessage("Separate median log ratio (by #Cysteines) per file, not displaying.");
        } //end if (medianCenterByNumCysteines)
        else
        {
            List<Float> logRatiosForMedianCalc = new ArrayList<Float>();
            fileMedianLogRatioMap = new HashMap<File, Float>();
            int fileIndex = 1;
            for (File featureFile : pepXmlFiles)
            {
                List<Float> logRatiosForMedianCalcThisFile = new ArrayList<Float>();

                ApplicationContext.infoMessage("\tProcessing file " + featureFile.getAbsolutePath() + "...");
                try
                {
                    Iterator<FeatureSet> featureSetIterator =
                            new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);

                    while (featureSetIterator.hasNext())
                    {
                        FeatureSet featureSet = featureSetIterator.next();
                        filterOnQualityScores(featureSet);
                        for (Feature feature : featureSet.getFeatures())
                        {
                            if (IsotopicLabelExtraInfoDef.hasRatio(feature) &&
                                    MS2ExtraInfoDef.getPeptideProphet(feature) >= minPeptideProphetForMedian)
                            {
                                float ratio =  (float) IsotopicLabelExtraInfoDef.getRatio(feature);

                                if (!Float.isInfinite(ratio) && !Float.isNaN(ratio) && ratio != 0)
                                {
                                    logRatiosForMedianCalcThisFile.add((float) Math.log(ratio));
                                }
                            }
                        }
                    }
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failed to load feature file " + featureFile);
                }
                if (!medianCenterAllRunsTogether)
                {
                    if (logRatiosForMedianCalcThisFile.size() < minRatiosForMedianCenter)
                        throw new CommandLineModuleExecutionException(
                                "Not enough ratios to calculate median for file "
                                + featureFile.getAbsolutePath() + " (only " + logRatiosForMedianCalcThisFile.size() +
                                        " ratios, needed " + minRatiosForMedianCenter + ")");
                    else
                    {
                        float medianLogRatioThisFile = (float)BasicStatistics.median(logRatiosForMedianCalcThisFile);
                        ApplicationContext.infoMessage("Median log ratio for file " +
                                featureFile.getName() + ": " +
                                medianLogRatioThisFile);
                        fileMedianLogRatioMap.put(featureFile, medianLogRatioThisFile);
                        if (showCharts)
                        {
                            PanelWithHistogram pwh =
                                    new PanelWithHistogram(logRatiosForMedianCalcThisFile,
                                            "RAW Log Ratios " + fileIndex++, 200);
                            pwh.displayInTab();
                        }
                    }
                }
                else logRatiosForMedianCalc.addAll(logRatiosForMedianCalcThisFile);
            }

            if (medianCenterAllRunsTogether)
            {
                if (logRatiosForMedianCalc.size() < minRatiosForMedianCenter)
                        throw new CommandLineModuleExecutionException(
                                "Not enough ratios to calculate median (only " +
                                        logRatiosForMedianCalc.size() +
                                        " ratios, needed " + minRatiosForMedianCenter + ")");
                //assign the same median to each file.  This keeps code same down below
                float medianLogRatioAcrossAll = (float)BasicStatistics.median(logRatiosForMedianCalc);
                ApplicationContext.infoMessage("Median log ratio across all runs: " + medianLogRatioAcrossAll);
                for (File featureFile : pepXmlFiles)
                    fileMedianLogRatioMap.put(featureFile, medianLogRatioAcrossAll);
                if (showCharts)
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(logRatiosForMedianCalc, "RAW Log Ratios", 200);
                    pwh.displayInTab();
                }
            }
        } //end not-by-cysteines behavior
        ApplicationContext.infoMessage("Done calculating median log ratio across all files.");
    }



    protected void loadLightHeavyPeptidesAcrossAll()
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Loading light and heavy peptide occurrences across all files, " +
                "this may take a while...");

        for (File featureFile : pepXmlFiles)
        {
            ApplicationContext.infoMessage("\tProcessing file " + featureFile.getAbsolutePath() + "...");
            try
            {
                Iterator<FeatureSet> featureSetIterator =
                        new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);

                while (featureSetIterator.hasNext())
                {
                    FeatureSet featureSet = featureSetIterator.next();
                    //todo: adding this here for Lynn.  Is this always appropriate?  Sometimes, might want to know about low-quality stuff here
                    filterOnQualityScores(featureSet);
                    addLightHeavyPeptides(featureSet);
                }
            }
            catch (IOException e)
            {
                throw new CommandLineModuleExecutionException("Failed to load feature file " + featureFile);
            }
        }
        ApplicationContext.infoMessage("Done loading light (" + lightPeptidesAllRuns.size() + ") and heavy (" +
                heavyPeptidesAllRuns.size() + ") peptide occurrences across all files.");
    }

    protected String calcOutputFilename(String inputFilename)
    {
        if (inputFilename.toLowerCase().endsWith(".pep.xml"))
            return inputFilename.substring(0, inputFilename.length()-".pep.xml".length()) + ".mod.pep.xml";
        else if (inputFilename.toLowerCase().endsWith(".xml"))
            return inputFilename.substring(0, inputFilename.length()-".xml".length()) + ".mod.xml";
        else return inputFilename + ".mod.pep.xml";
    }

    /**
     *
     * @param featureFile
     * @param outputFile
     * @throws CommandLineModuleExecutionException
     */
    protected void handleFeatureFile(File featureFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Handling file " + featureFile.getAbsolutePath() + "...");

        if (maxFDR < Float.MAX_VALUE)
        {
            minPeptideProphet = (float) calcMinPeptideProphetForMaxFDR(featureFile);
            ApplicationContext.infoMessage("Minimum PeptideProphet for FDR " + maxFDR + ": " + minPeptideProphet);
        }

        try
        {
            Iterator<FeatureSet> featureSetIterator =
                    new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);

            List<File> tempFeatureFiles = new ArrayList<File>();
            int numSetsProcessed = 0;
            while (featureSetIterator.hasNext())
            {
                FeatureSet featureSet = featureSetIterator.next();

                ApplicationContext.infoMessage("\tProcessing fraction " + (numSetsProcessed+1) + "...");

                processFeatureSet(featureSet);
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
//                while (baseName.contains(".." + File.separator))
//                    baseName.replaceFirst(".." + File.separator, "");
                if (baseName.contains(File.separator))
                    baseName = baseName.substring(baseName.lastIndexOf(File.separator) + 1);

                File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);

                featureSet.savePepXml(thisFractionFile);

                _log.debug("Saved fraction file as " + thisFractionFile.getAbsolutePath());
                tempFeatureFiles.add(thisFractionFile);

                numSetsProcessed++;
            }

            ApplicationContext.infoMessage("Saving output file " + outputFile.getAbsolutePath() + " ...");

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
            ApplicationContext.infoMessage("Done.");
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to process features from file " +
                    featureFile.getAbsolutePath(),e);
        }
        finally
        {
            TempFileManager.deleteTempFiles(this);            
        }
    }

    protected void filterOnQualityScores(FeatureSet featureSet)
    {
        //filter on feature attributes
        if (minPeptideProphet > 0f)
        {
            List<Feature> featuresToKeep = new ArrayList<Feature>();

            int numFeaturesStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (MS2ExtraInfoDef.getPeptideProphet(feature) >= minPeptideProphet)
                    featuresToKeep.add(feature);
                else
                    numFeaturesStripped++;
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));
            ApplicationContext.infoMessage("\tStripped " + numFeaturesStripped +
                    " features with PeptideProphet < " + minPeptideProphet);
        }



        if (maxExpect < Float.MAX_VALUE)
        {
            List<Feature> featuresToKeep = new ArrayList<Feature>();

            int numFeaturesStripped = 0;

            for (Feature feature : featureSet.getFeatures())
            {
                try
                {
                    Float expectScore =
                            Float.parseFloat(MS2ExtraInfoDef.getSearchScore(feature,"expect"));
                    if (expectScore <= maxExpect)
                        featuresToKeep.add(feature);
                    else
                        numFeaturesStripped++;
                }
                catch (Exception e)
                {}
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));
            ApplicationContext.infoMessage("\tStripped " + numFeaturesStripped +
                                " features with expect > " + maxExpect);
        }
    }

    protected int countCysteines(String peptide)
    {
        int numCysteines = 0;
        for (int i=0; i<peptide.length(); i++)
        {
            if (peptide.charAt(i) == 'C')
                numCysteines++;
        }
        return numCysteines;
    }

    protected void processFeatureSet(FeatureSet featureSet)
    {

        if (filterByProteinPrefix)
        {
            filterByProteinPrefix(featureSet);
        }



        filterOnQualityScores(featureSet);

        if (stripLightIDs)
        {
            int numFeaturesBefore = featureSet.getFeatures().length;
            List<Feature> featuresToKeep = new ArrayList<Feature>();

            for (Feature feature : featureSet.getFeatures())
            {
                if (!IsotopicLabelExtraInfoDef.isLightLabeled(feature, labelType))
                    featuresToKeep.add(feature);
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));

            ApplicationContext.setMessage("\tStripped " + (numFeaturesBefore - featureSet.getFeatures().length) +
                    " light-identified peptides. Kept " +
                    featureSet.getFeatures().length +
                    " out of " + numFeaturesBefore + " identifications.");
        }

        if (peptidesToStrip != null)
        {
            int numFeaturesBefore = featureSet.getFeatures().length;
            List<Feature> featuresToKeep = new ArrayList<Feature>();
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide == null)
                    continue;
                if (!peptidesToStrip.contains(peptide))
                    featuresToKeep.add(feature);
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));

            ApplicationContext.setMessage("\tStripped indicated peptides. Kept " + featureSet.getFeatures().length +
                    " out of " + numFeaturesBefore + " identifications.");
        }
        if (proteinsToStrip != null)
        {
            int numFeaturesBefore = featureSet.getFeatures().length;
            Set<String> proteinsNotOnStripListWithStrippedPeptides = new HashSet<String>();
            List<Feature> featuresToKeep = new ArrayList<Feature>();
            for (Feature feature : featureSet.getFeatures())
            {
                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                if (proteins == null)
                    continue;
                Set<String> proteinsNotOnStripListThisFeature = new HashSet<String>();
                boolean foundBadProtein = false;
                for (String protein : proteins)
                {
                    if (proteinsToStrip.contains(protein))
                    {
                        foundBadProtein = true;
                    }
                    else
                        proteinsNotOnStripListThisFeature.add(protein);

                }
                if (foundBadProtein)
                {
                    proteinsNotOnStripListWithStrippedPeptides.addAll(proteinsNotOnStripListThisFeature);
                }
                else
                {
                    featuresToKeep.add(feature);
                }
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));
//            StringBuffer proteinsNotOnListWithStrippedBuf =
//                    new StringBuffer("Proteins not on list that have stripped peptides");
//            for (String protein : proteinsNotOnStripListWithStrippedPeptides)
//                proteinsNotOnListWithStrippedBuf.append("," + protein);
//            ApplicationContext.setMessage(proteinsNotOnListWithStrippedBuf.toString());
            ApplicationContext.setMessage("\tStripped indicated proteins. Kept " + featureSet.getFeatures().length +
                    " out of " + numFeaturesBefore + " identifications.");
        }

        if (stripQuantMatchingRegexp)
        {
            int numStrippedQuant = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide == null)
                    continue;
                if (peptide.matches(regexp) && IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    numStrippedQuant++;
                }
            }
            ApplicationContext.setMessage("\tStripped quantitation from " + numStrippedQuant +
                    " features with peptides matching regexp " + regexp);
        }


        if (proteinsToStripQuant != null)
        {
            Set<String> proteinsNotOnStripQuantListWithStrippedPeptides = new HashSet<String>();
            for (Feature feature : featureSet.getFeatures())
            {
                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                if (proteins == null)
                    continue;
                Set<String> proteinsNotOnStripQuantListThisFeature = new HashSet<String>();
                boolean foundBadProtein = false;
                for (String protein : proteins)
                {
                    if (proteinsToStripQuant.contains(protein))
                    {
                        foundBadProtein = true;
                    }
                    else
                        proteinsNotOnStripQuantListThisFeature.add(protein);

                }
                if (foundBadProtein)
                {
                    proteinsNotOnStripQuantListWithStrippedPeptides.addAll(proteinsNotOnStripQuantListThisFeature);
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                }
            }

            StringBuffer proteinsNotOnListWithStrippedBuf =
                    new StringBuffer("Proteins not on list that have quantitation-stripped peptides");
            for (String protein : proteinsNotOnStripQuantListWithStrippedPeptides)
                proteinsNotOnListWithStrippedBuf.append("," + protein);
            ApplicationContext.setMessage(proteinsNotOnListWithStrippedBuf.toString());
            ApplicationContext.setMessage("\tStripped quantitative events from indicated proteins.");
        }

        if (proteinsToKeep != null)
        {
            int numFeaturesBefore = featureSet.getFeatures().length;
            List<Feature> featuresToKeep = new ArrayList<Feature>();
            for (Feature feature : featureSet.getFeatures())
            {
                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                if (proteins == null)
                    continue;
                boolean shouldKeep = false;
                for (String protein : proteins)
                {
                    if (proteinsToKeep.contains(protein))
                    {
                        shouldKeep = true;
                        break;
                    }
                }
                if (shouldKeep)
                {
                    featuresToKeep.add(feature);
                }
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));
            ApplicationContext.setMessage("\tKept only indicated proteins. Kept " + featureSet.getFeatures().length +
                    " out of " + numFeaturesBefore + " identifications.");
        }

        if (stripQuantMissingLightOrHeavyWithinRun)
        {
            stripQuantWithoutLightOrHeavyIDWithinSet(featureSet);
        }

        if (minQuantPeptideProphet > 0f)
        {
            int numFeaturesQuantStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (MS2ExtraInfoDef.getPeptideProphet(feature) < minQuantPeptideProphet
                        && IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    numFeaturesQuantStripped++;
                }
            }
            ApplicationContext.infoMessage("\tStripped quantitation from " + numFeaturesQuantStripped +
                    " features with PeptideProphet < " + minQuantPeptideProphet);
        }

        if (minRatio > 0f)
        {
            int numFeaturesQuantStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature) &&
                    IsotopicLabelExtraInfoDef.getRatio(feature) < minRatio)
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    numFeaturesQuantStripped++;
                }
            }
            ApplicationContext.infoMessage("\tStripped quantitation from " + numFeaturesQuantStripped +
                    " features with ratio < " + minRatio);
        }

        if (maxRatio < Float.MAX_VALUE)
        {
            int numFeaturesQuantStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature) &&
                    IsotopicLabelExtraInfoDef.getRatio(feature) > maxRatio)
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    numFeaturesQuantStripped++;
                }
            }
            ApplicationContext.infoMessage("\tStripped quantitation from " + numFeaturesQuantStripped +
                    " features with ratio > " + maxRatio);
        }

        if (maxQuantExpect < Float.MAX_VALUE)
        {
            int numFeaturesQuantStripped = 0;
            int numQuantFeaturesWithNoExpect = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                String expectString = MS2ExtraInfoDef.getSearchScore(feature,"expect");
                if (expectString == null)
                {
                    if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                        numQuantFeaturesWithNoExpect++;
                    continue;
                }
                Float expectScore =
                        Float.parseFloat(expectString);

                if (expectScore > maxQuantExpect && IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    numFeaturesQuantStripped++;
                }
            }
            ApplicationContext.infoMessage("\tStripped quantitation from " + numFeaturesQuantStripped +
                    " features with expect > " + maxQuantExpect);
            if (numQuantFeaturesWithNoExpect > 0)
                 ApplicationContext.infoMessage("\t\tNote: " + numQuantFeaturesWithNoExpect +
                         " quantified features had no expect score, were left alone");
        }


        if (stripQuantNotInHeavyAcrossAll || stripQuantMissingLightOrHeavyAcrossAll)
        {
            Set<String> strippedPeptidesThisRun = new HashSet<String>();
            Set<String> allPeptides = new HashSet<String>();
            Set<String> allQuantifiedPeptides = new HashSet<String>();

            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide == null)
                    continue;
                allPeptides.add(peptide);
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    allQuantifiedPeptides.add(peptide);
                else continue;

                boolean shouldRemove = false;

                //only one of the two variables will be true, so take the appropriate action
                if (stripQuantMissingLightOrHeavyAcrossAll)
                {
                    if (!heavyPeptidesAllRuns.contains(peptide) || !lightPeptidesAllRuns.contains(peptide))
                        shouldRemove = true;
                }
                //not MissingLightOrHeavy, so must be NotInHeavy
                else
                {
                    if (!heavyPeptidesAllRuns.contains(peptide))
                        shouldRemove = true;
                }
                if (shouldRemove)
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    strippedPeptidesThisRun.add(peptide);
                }
            }
            ApplicationContext.infoMessage("\tStripped quant from " + strippedPeptidesThisRun.size() +
                    " peptides (out of " + allPeptides.size() + " total, " + allQuantifiedPeptides.size() +
                    " quantified) not found in appropriate states in any run");
        }



        //filter on fractional deltamass
        if (maxFracDeltaMass != null)
        {
            List<Feature> featuresToKeep = new ArrayList<Feature>();

            int numFeaturesStripped = 0;

            for (Feature feature : featureSet.getFeatures())
            {
                float deltaMass = MS2ExtraInfoDef.getDeltaMass(feature);

                float fracDeltaMass = (float) ((deltaMass + .5) %1) - .5f;
                if (fracDeltaMass < -.5f)
                    fracDeltaMass += 1f;

                if (maxFracDeltaMass.getDeltaMassType() == AnalyzeICAT.DELTA_MASS_PPM)
                    fracDeltaMass = fracDeltaMass * 1000000f / feature.getMass();

                if (fracDeltaMass <= maxFracDeltaMass.getDeltaMass())
                    featuresToKeep.add(feature);
                    else
                        numFeaturesStripped++;
            }
            featureSet.setFeatures(featuresToKeep.toArray(new Feature[0]));
            ApplicationContext.infoMessage("\tStripped " + numFeaturesStripped +
                                " features with fractional deltamass > " + maxFracDeltaMass);
        }        

        if (adjustQuantZeroAreas)
        {
            adjustQuantZeroAreas(featureSet);
        }

        if (stripQuantZeroAreas)
        {
            int numQuantStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature);
                    float heavyIntensity = (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);

                    if (lightIntensity == 0 || heavyIntensity == 0)
                    {
                       IsotopicLabelExtraInfoDef.removeRatio(feature);
                        numQuantStripped++;
                    }
                }
            }
            ApplicationContext.infoMessage("Stripped quantitation from " + numQuantStripped +
                    " features with zero light and/or heavy area");
        }

        if (stripQuantSingleScans)
        {
            int numQuantStripped = 0;
            for (Feature feature : featureSet.getFeatures())
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    if ((IsotopicLabelExtraInfoDef.getLightFirstScan(feature) ==
                         IsotopicLabelExtraInfoDef.getLightLastScan(feature)) ||
                        (IsotopicLabelExtraInfoDef.getHeavyFirstScan(feature) ==
                         IsotopicLabelExtraInfoDef.getHeavyLastScan(feature)))
                    {
                        IsotopicLabelExtraInfoDef.removeRatio(feature);
                        numQuantStripped++;
                    }
                }
            }
            ApplicationContext.infoMessage("Stripped quantitation from " + numQuantStripped +
                    " features with a single scan for light or heavy species");
        }


        if (medianCenter)
        {
            if (medianCenterByNumCysteines)
            {
                Map<Integer, Float> numCysteineMedianLogRatioMapThisFile =
                        fileNumCysteinesMedianLogRatioMap.get(featureSet.getSourceFile());
                for (Feature feature : featureSet.getFeatures())
                {
                    if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        float ratio =  (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                        if (peptide != null)
                        {
                            int numCysteines = countCysteines(peptide);
//if (peptide.equals("QSSGENCDVVVNTLGK")) throw new RuntimeException("heavy_area is NaN! Light: " +
//IsotopicLabelExtraInfoDef.getLightIntensity(feature) + ", old heavy: " +
//        IsotopicLabelExtraInfoDef.getHeavyIntensity(feature) + ", old ratio: " + ratio + ", numcys: " + numCysteines +
//     ", containsKey? " + numCysteineMedianLogRatioMapThisFile.containsKey(numCysteines));
                            if (numCysteineMedianLogRatioMapThisFile.containsKey(numCysteines))
                            {
                                float newRatio = (float) Math.exp(Math.log(ratio) -
                                        numCysteineMedianLogRatioMapThisFile.get(numCysteines));
                                IsotopicLabelExtraInfoDef.setRatio(feature,newRatio);
                                //dhmay 20091215, fixing divide-by-zero area that left NaN in output pepXML heavy_areas
                                double newHeavyIntensity = IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
                                if (newRatio > 0)
                                    newHeavyIntensity =
                                            IsotopicLabelExtraInfoDef.getLightIntensity(feature) / newRatio;
//if (peptide.equals("QSSGENCDVVVNTLGK")) throw new RuntimeException("heavy_area is NaN! Light: " +
//IsotopicLabelExtraInfoDef.getLightIntensity(feature) + ", old heavy: " +
//        IsotopicLabelExtraInfoDef.getHeavyIntensity(feature) + ", new heavy: " + newHeavyIntensity + ", new ratio: " + newRatio);
                                IsotopicLabelExtraInfoDef.setHeavyIntensity(feature, newHeavyIntensity);
                            }
                        }
                    }
                }
            }
            else
            {
                float medianLogRatioThisFile = 0;
                if (medianRatioForCentering != -1)
                    medianLogRatioThisFile = (float) Math.log(medianRatioForCentering);
                else medianLogRatioThisFile = fileMedianLogRatioMap.get(featureSet.getSourceFile());
                ApplicationContext.infoMessage("Median centering file " + featureSet.getSourceFile().getName() +
                        ", median: " + medianLogRatioThisFile);
//List<Float> newLogRatios = new ArrayList<Float>();
//List<Float> oldLogRatios = new ArrayList<Float>();

                for (Feature feature : featureSet.getFeatures())
                {
                    if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        float ratio =  (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                        float newRatio = (float) Math.exp(Math.log(ratio) - medianLogRatioThisFile);
//newLogRatios.add((float)Math.log(newRatio));
//oldLogRatios.add((float)Math.log(ratio));

                        IsotopicLabelExtraInfoDef.setRatio(feature, newRatio);
                        //dhmay 20091215, fixing divide-by-zero area that left NaN in output pepXML heavy_areas
                        double newHeavyIntensity = IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
                        if (newRatio > 0) newHeavyIntensity =
                                IsotopicLabelExtraInfoDef.getLightIntensity(feature) / newRatio;
                        IsotopicLabelExtraInfoDef.setHeavyIntensity(feature, newHeavyIntensity);
                    }
                }
//System.err.println("Old median log ratio: " + BasicStatistics.median(oldLogRatios));
//System.err.println("New median log ratio: " + BasicStatistics.median(newLogRatios));

            }
        }

    }

    protected void adjustQuantZeroAreas(FeatureSet featureSet)
    {
        List<Float> allNonzeroAreas = new ArrayList<Float>();
        for (Feature feature : featureSet.getFeatures())
        {
            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
                if (IsotopicLabelExtraInfoDef.getLightIntensity(feature) > 0)
                    allNonzeroAreas.add((float) IsotopicLabelExtraInfoDef.getLightIntensity(feature));
                if (IsotopicLabelExtraInfoDef.getHeavyIntensity(feature) > 0)
                   allNonzeroAreas.add((float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature));
            }
        }

        float adjustedPercentileValue = (float)BasicStatistics.percentile(allNonzeroAreas,
                percentileForQuantZeroAreaAdjustment);

        int numFeaturesQuantAdjusted = 0;
        
        for (Feature feature : featureSet.getFeatures())
        {
            boolean adjusted = false;
            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
                float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature);
                float heavyIntensity = (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);

                if (lightIntensity == 0)
                {
                    IsotopicLabelExtraInfoDef.setLightIntensity(feature,
                            Math.min(heavyIntensity, adjustedPercentileValue));
                    adjusted = true;
                }
                if (heavyIntensity == 0)
                {
                    IsotopicLabelExtraInfoDef.setHeavyIntensity(feature,
                            Math.min(lightIntensity, adjustedPercentileValue));
                    adjusted = true;
                }
            }
            if (adjusted)
            {
                float newRatio = (float) (IsotopicLabelExtraInfoDef.getLightIntensity(feature) /
                        IsotopicLabelExtraInfoDef.getHeavyIntensity(feature));

                IsotopicLabelExtraInfoDef.setRatio(feature, newRatio);
                numFeaturesQuantAdjusted++;
            }

        }

        if (numFeaturesQuantAdjusted == 0)
            ApplicationContext.infoMessage("\tNo quantitated features with zero light or heavy areas to adjust");
        else
            ApplicationContext.infoMessage("\tAdjusted zero values for light or heavy quantitation areas to " +
                adjustedPercentileValue + " (percentile " + percentileForQuantZeroAreaAdjustment  +
                ") for " +numFeaturesQuantAdjusted + " features");
    }

    /**
     * Does filtering by protein prefix. if excludeProteinPrefixQuantOnly is true, then this
     * "filtering" just means removing quantitation
     * @param featureSet
     */
    protected void filterByProteinPrefix(FeatureSet featureSet)
    {
        List<Feature> outFeatureList = new ArrayList<Feature>();
        for (Feature feature : featureSet.getFeatures())
        {
            boolean exclude = false;

            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);

            if (peptide == null)
                exclude=true;
            else
            {
                List<String> proteins = MS2ExtraInfoDef.getProteinList(feature);
                if (proteins != null)
                {
                    if (badProteinPrefix != null)
                    {
                        for (String protein : proteins)
                        {
                            if (protein.startsWith(badProteinPrefix))
                            {
                                exclude=true;
                                break;
                            }
                        }
                    }
                    else
                    {
                        exclude=true;
                        for (String protein : proteins)
                            if (protein.startsWith(goodProteinPrefix))
                            {
                                exclude=false;
                                break;
                            }
                    }
                }
            }
            if (excludeProteinPrefixQuantOnly)
            {
                outFeatureList.add(feature);
                if (exclude)
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
            }
            else
            {
                if (!exclude)
                    outFeatureList.add(feature);
            }
        }
        featureSet.setFeatures(outFeatureList.toArray(new Feature[outFeatureList.size()]));
    }

    protected Map<String, Map<Integer, List<Feature>>> createPeptideChargeFeatureListMap(FeatureSet featureSet)
    {
        Map<String, Map<Integer, List<Feature>>> peptideChargeFeatureListMap =
                new HashMap<String, Map<Integer, List<Feature>>>();
        for (Feature feature : featureSet.getFeatures())
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            Map<Integer, List<Feature>> thisPeptideChargeFeatureListMap =
                peptideChargeFeatureListMap.get(peptide);
            if (thisPeptideChargeFeatureListMap == null)
            {
                thisPeptideChargeFeatureListMap = new HashMap<Integer, List<Feature>>();
                peptideChargeFeatureListMap.put(peptide, thisPeptideChargeFeatureListMap);
            }

            int charge = feature.getCharge();
            List<Feature> featureList = thisPeptideChargeFeatureListMap.get(charge);
            if (featureList == null)
            {
                featureList = new ArrayList<Feature>();
                thisPeptideChargeFeatureListMap.put(charge, featureList);
            }
            featureList.add(feature);
        }
        return peptideChargeFeatureListMap;
    }

    protected void addLightHeavyPeptides(FeatureSet featureSet)
    {
//        Map<Integer, Pair<Set<String>,Set<String>>> result = new HashMap<Integer, Pair<Set<String>,Set<String>>>();

        for (Feature feature : featureSet.getFeatures())
        {
            String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
            if (peptide == null)
                continue;
            if (IsotopicLabelExtraInfoDef.isLightLabeled(feature, labelType))
            {
                lightPeptidesAllRuns.add(peptide);
            }
            else if (IsotopicLabelExtraInfoDef.isHeavyLabeled(feature, labelType))
            {
                heavyPeptidesAllRuns.add(peptide);
            }
        }
    }

    protected void stripQuantWithoutLightOrHeavyIDWithinSet(FeatureSet featureSet)
    {
        Map<String, Map<Integer, List<Feature>>> peptideChargeFeatureListMap =
                createPeptideChargeFeatureListMap(featureSet);


        int numFeaturesQuantStripped=0;
        for (String peptide : peptideChargeFeatureListMap.keySet())
        {
            for (List<Feature> featuresThisPeptideThisCharge : peptideChargeFeatureListMap.get(peptide).values())
            {
                boolean foundLight = false;
                boolean foundHeavy = false;

                for (Feature feature : featuresThisPeptideThisCharge)
                {
                    if (IsotopicLabelExtraInfoDef.isLightLabeled(feature, labelType))
                    {
                        foundLight=true;
                        continue;
                    }
                    if (IsotopicLabelExtraInfoDef.isHeavyLabeled(feature, labelType))
                    {
                        foundHeavy=true;
                        continue;
                    }
                }

                if (!foundLight || !foundHeavy)
                {
//                    ApplicationContext.setMessage("Stripping ratios from peptide " + peptide + ", charge " +
//                            featuresThisPeptideThisCharge.get(0).getCharge());
                    for (Feature feature : featuresThisPeptideThisCharge)
                    {
                        if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                        {
                            IsotopicLabelExtraInfoDef.removeRatio(feature);
                        }
                        numFeaturesQuantStripped++;
                    }
                }
            }
        }
        ApplicationContext.setMessage("\tStripped ratios from " + numFeaturesQuantStripped +
                " features for peptides not sequenced in both light and heavy");
    }

    /**
     * median-centers the natural log on 0, not in place
     * @param allInputList
     * @param inputListForMedianCalc this gets munged
     * @return
     */
    protected List<Float> logMedianCenterOn0(List<Float> allInputList, List<Float> inputListForMedianCalc,
                                             String chartTitleSuffix)
    {
        for (int i=0; i<inputListForMedianCalc.size(); i++)
            inputListForMedianCalc.set(i, (float) Math.log(inputListForMedianCalc.get(i)));
        if (showCharts)
        {
            PanelWithHistogram pwhBefore =
                    new PanelWithHistogram(inputListForMedianCalc, "Mediancalc logratios" + chartTitleSuffix);
            pwhBefore.displayInTab();
        }

        float medianLog = (float)BasicStatistics.median(inputListForMedianCalc);

        ApplicationContext.setMessage("\t\tSubtracting " + medianLog + " (old median) from all LOG ratios");

        List<Float> outputList = new ArrayList<Float>();
        for (float input : allInputList)
            outputList.add((float) Math.exp(Math.log(input) - medianLog));
        return outputList;
    }

    protected double calcMinPeptideProphetForMaxFDR(File pepXmlFile)
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Calculating min PeptideProphet probability for file " + pepXmlFile.getAbsolutePath());
        try
        {
            Iterator<FeatureSet> featureSetIterator =
                    new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
            List<Double> allProbabilities = new ArrayList<Double>();
            while (featureSetIterator.hasNext())
            {
                for (Feature feature : featureSetIterator.next().getFeatures())
                    allProbabilities.add(MS2ExtraInfoDef.getPeptideProphet(feature));
            }
            Collections.sort(allProbabilities);

            double cumulativeFalseMatches = 0;
            int count = 0;

            for (int i=allProbabilities.size()-1; i>=0; i--)
            {
                cumulativeFalseMatches += (1 - allProbabilities.get(i));
                count++;

                if (cumulativeFalseMatches / count > maxFDR)
                {
                    return allProbabilities.get(i);
                }
            }
            return 0;
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to process pepXML file " +
                    pepXmlFile.getAbsolutePath(),e);
        }
    }

}
