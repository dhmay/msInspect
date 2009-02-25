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
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;

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
    protected boolean stripQuantMissingLightOrHeavyWithinRun = false;
    protected boolean stripQuantMissingLightOrHeavyAcrossAll = false;
    protected boolean stripQuantNotInHeavyAcrossAll = false;
    protected boolean adjustQuantZeroAreas = false;
    protected boolean stripQuantZeroAreas = false;

    protected boolean filterByProteinPrefix = false;

    protected int percentileForQuantZeroAreaAdjustment = 1;

    protected String badProteinPrefix;
    protected String goodProteinPrefix;
    protected boolean excludeProteinPrefixQuantOnly;

    protected File outFile;
    protected File outDir;

    protected Set<String> peptidesToStrip = null;

    protected boolean showCharts = false;

    protected int minRatiosForMedianCenter = 10;

    protected float minPeptideProphetForMedian = .75f;

    protected float minPeptideProphet = 0f;
    protected float minQuantPeptideProphet = 0f;

    protected float maxExpect = Float.MAX_VALUE;
    protected float maxQuantExpect = Float.MAX_VALUE;

    protected boolean requirePepXmlExtension = false;

    protected DeltaMassArgumentDefinition.DeltaMassWithType maxFracDeltaMass = null;


    public static final float MOD_MASS_SLOP = 0.1f;


    protected int labelType = LABEL_ACRYLAMIDE;

    public static final int LABEL_ACRYLAMIDE = 0;
    public static final int LABEL_LYCINE = 1;

    public static final float ACRYLAMIDE_LABEL_LIGHTMASS = 174.0458f;
    public static final float ACRYLAMIDE_LABEL_HEAVYMASS = 177.05591f;

    public static final float SILAC_LABEL_MASS = 134.115092f;
    


    //for filtering out peptides not seen in both light and heavy in some run
    protected Set<String> lightPeptidesAllRuns = new HashSet<String>();
    protected Set<String> heavyPeptidesAllRuns = new HashSet<String>();




    protected String[] labelStrings = new String[]
            {
                    "acrylamide",
                    "silac"
            };

    protected String[] labelExplanations = new String[]
            {
                    "Acrylamide (3.0106Da on C)",
                    "SILAC Lycine labeling (134.115092 on K)"
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
                       new FileToReadArgumentDefinition("strippeptidefile", false,
                               "File containing a list of peptides to strip from results, one per line, all caps"),
                       new FileToWriteArgumentDefinition("out", false, "Output file"),
                       new DirectoryToReadArgumentDefinition("outdir", false, "Output directory"),
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

        medianCenter = getBooleanArgumentValue("mediancenter");
        medianCenterByNumCysteines = getBooleanArgumentValue("bynumcysteines");

        requirePepXmlExtension = getBooleanArgumentValue("requirepepxmlextension");


        stripQuantMissingLightOrHeavyWithinRun = getBooleanArgumentValue("stripquantmissinglightorheavywithinrun");
        stripQuantMissingLightOrHeavyAcrossAll = getBooleanArgumentValue("stripquantmissinglightorheavyacrossruns");
        stripQuantNotInHeavyAcrossAll = getBooleanArgumentValue("stripquantmissingheavy");

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

        if (!medianCenter)
            assertArgumentAbsent("bynumcysteines","mediancenter");

        labelType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("label")).getIndexForArgumentValue(getStringArgumentValue("label"));


        if (hasArgumentValue("strippeptidefile"))
        {
            File peptidesToStripFile = getFileArgumentValue("strippeptidefile");
            try
            {
                FileReader fr = new FileReader(peptidesToStripFile);
                BufferedReader br = new BufferedReader(fr);
                peptidesToStrip = new HashSet<String>();
                String line = null;
                while ((line = br.readLine()) != null)
                {
                    String peptide = StringUtils.strip(line);
                    boolean ok = true;
                    if (peptide.length() == 0)
                        ok = false;
                    for (int i = 0; i < peptide.length(); i++) {
                        char c = peptide.charAt(i);
                        if ((c >= 'A') && (c <= 'Z')) continue; // uppercase
                        ok = false;
                    }
                    if (!ok)
                        throw new RuntimeException();

                    peptidesToStrip.add(peptide);
                }

            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Failed to retrieve list of peptides to strip from file " +
                        peptidesToStripFile + ".  Please make sure file contains only a list of peptides, one per line, all caps");
            }
        }



        minRatiosForMedianCenter = getIntegerArgumentValue("minmediancenterratios");

        showCharts = getBooleanArgumentValue("showcharts");

        minPeptideProphetForMedian = (float) getDoubleArgumentValue("minmedianpprophet");
        minPeptideProphet = (float) getDoubleArgumentValue("minpprophet");
        minQuantPeptideProphet = (float) getDoubleArgumentValue("minquantpprophet");

        maxExpect = (float) getDoubleArgumentValue("maxexpect");
        maxQuantExpect = (float) getDoubleArgumentValue("maxquantexpect");

        if (hasArgumentValue("maxfracdeltamass"))
            maxFracDeltaMass = getDeltaMassArgumentValue("maxfracdeltamass");

        if (peptidesToStrip == null &&
                !medianCenter && !stripQuantMissingLightOrHeavyWithinRun && !stripQuantMissingLightOrHeavyAcrossAll &&
                !filterByProteinPrefix &&
                !stripQuantNotInHeavyAcrossAll && !adjustQuantZeroAreas && !stripQuantZeroAreas &&
                !hasArgumentValue("maxexpect") && minPeptideProphet == 0 && minQuantPeptideProphet == 0)
        {
            throw new ArgumentValidationException("Nothing to do!  Quitting");
        }

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {

        if (stripQuantNotInHeavyAcrossAll || stripQuantMissingLightOrHeavyAcrossAll)
            loadLightHeavyPeptidesAcrossAll();

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
        ApplicationContext.infoMessage("Done loading light and heavy peptide occurrences across all files.");
    }

    protected String calcOutputFilename(String inputFilename)
    {
        if (inputFilename.toLowerCase().endsWith(".pep.xml"))
            return inputFilename.substring(0, inputFilename.length()-".pep.xml".length()) + ".mod.pep.xml";
        else if (inputFilename.toLowerCase().endsWith(".xml"))
            return inputFilename.substring(0, inputFilename.length()-".xml".length()) + ".mod.xml";
        else return inputFilename + ".mod.pep.xml";
    }


    protected void handleFeatureFile(File featureFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        try
        {
            ApplicationContext.infoMessage("Loading file " + featureFile.getAbsolutePath() + "...");
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
                while (baseName.contains(".." + File.separator))
                    baseName.replaceFirst(".." + File.separator, "");

                File thisFractionFile = TempFileManager.createTempFile(baseName + ".pep.xml", this);

                featureSet.savePepXml(thisFractionFile);

                tempFeatureFiles.add(thisFractionFile);

                numSetsProcessed++;
            }

            ApplicationContext.infoMessage("Saving output file " + outputFile.getAbsolutePath() + "...");

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

    protected void processFeatureSet(FeatureSet featureSet)
    {
        if (filterByProteinPrefix)
        {
            filterByProteinPrefix(featureSet);
        }

        filterOnQualityScores(featureSet);
                 

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
                    " out of " + numFeaturesBefore + " peptide hits.");
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
            ApplicationContext.infoMessage("\tStripped " + strippedPeptidesThisRun.size() +
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


        if (medianCenter)
        {
            if (medianCenterByNumCysteines)
            {
                Map<Integer, List<Feature>> numCysteinesFeatureListMap =
                        new HashMap<Integer, List<Feature>>();
                for (Feature feature : featureSet.getFeatures())
                {
                    String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                    if (peptide != null)
                    {
                        int numCysteines = 0;
                        for (int i=0; i<peptide.length(); i++)
                        {
                            if (peptide.charAt(i) == 'C')
                                numCysteines++;
                        }
                        if (numCysteines > 0)
                        {
                            List<Feature> featureList = numCysteinesFeatureListMap.get(numCysteines);
                            if (featureList == null)
                            {
                                featureList = new ArrayList<Feature>();
                                numCysteinesFeatureListMap.put(numCysteines, featureList);
                            }

                            featureList.add(feature);
                        }
                    }
                }

                for (int i=0; i<20; i++)
                {
                    if (numCysteinesFeatureListMap.containsKey(i))
                    {
                        ApplicationContext.setMessage("\tCentering ratios for " + i + " Cysteines...");
                        logMedianCenterPeptideRatios(numCysteinesFeatureListMap.get(i).toArray(new Feature[0]),
                                ", " + i + " C's");
                    }
                }
            }
            else
            {
                ApplicationContext.setMessage("\tCentering ratios all peptides...");
                logMedianCenterPeptideRatios(featureSet.getFeatures(), "");
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

        float adjustedPercentileValue = BasicStatistics.percentile(allNonzeroAreas,
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
                            if (protein.startsWith(badProteinPrefix))
                            {
                                exclude=true;
                                break;
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
            if (isLightLabeled(feature))
            {
                lightPeptidesAllRuns.add(peptide);
            }
            else if (isHeavyLabeled(feature))
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
                    if (isLightLabeled(feature))
                    {
                        foundLight=true;
                        continue;
                    }
                    if (isHeavyLabeled(feature))
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

    protected boolean checkForModAllResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {
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

    protected boolean checkForModNoResidues(String peptide, List<ModifiedAminoAcid>[] mods, char residue, float modMass)
    {

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

    protected boolean logMedianCenterPeptideRatios(Feature[] features, String chartTitleSuffix)
    {        
        List<Float> allRatios = new ArrayList<Float>();
        List<Float> ratiosForMedianCalc = new ArrayList<Float>();

        List<Feature> featuresWithRatios = new ArrayList<Feature>();

        for (Feature feature : features)
        {
            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
            {
                float ratio =  (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                allRatios.add(ratio);
                featuresWithRatios.add(feature);

                if (MS2ExtraInfoDef.getPeptideProphet(feature) >= minPeptideProphetForMedian)
                {
                     ratiosForMedianCalc.add(ratio);
                }
            }
        }

        if (ratiosForMedianCalc.size() >= minRatiosForMedianCenter)
        {
            List<Float> newRatios = logMedianCenterOn0(allRatios, ratiosForMedianCalc, chartTitleSuffix);
            for (int i=0; i<allRatios.size(); i++)
            {
                Feature feature = featuresWithRatios.get(i);
                IsotopicLabelExtraInfoDef.setRatio(feature, newRatios.get(i));
                float newHeavyIntensity =
                        (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature) / newRatios.get(i);
                IsotopicLabelExtraInfoDef.setHeavyIntensity(feature, newHeavyIntensity);
            }

            if (showCharts)
            {
                List<Float> logOldRatios = new ArrayList<Float>();
                for (Float oldRatio : allRatios)
                    logOldRatios.add((float) Math.log(oldRatio));

                List<Float> logNewRatios = new ArrayList<Float>();
                for (Float newRatio : newRatios)
                    logNewRatios.add((float) Math.log(newRatio));



                PanelWithHistogram pwhAfter = new PanelWithHistogram(logNewRatios, "Log ratios after" + chartTitleSuffix);
                pwhAfter.displayInTab();
            }

            return true;
        }
        else
        {
            ApplicationContext.setMessage("Not enough ratios (" + ratiosForMedianCalc.size() + ").  Not recentering.");

            return false;
        }
    }

    /**
     * median-centers the log on 0, not in place
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

        float medianLog = BasicStatistics.median(inputListForMedianCalc);

        List<Float> outputList = new ArrayList<Float>();
        for (float input : allInputList)
            outputList.add((float) Math.exp(Math.log(input) - medianLog));
        return outputList;
    }

}
