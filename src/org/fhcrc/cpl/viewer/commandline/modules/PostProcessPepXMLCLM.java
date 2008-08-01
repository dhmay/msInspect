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
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.DeltaMassArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.AnalyzeICAT;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;
import org.apache.commons.lang.StringUtils;
import org.fhcrc.cpl.toolbox.BasicStatistics;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;

import java.io.*;
import java.util.*;


/**
 * post-process a pepxml file
 */
public class PostProcessPepXMLCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PostProcessPepXMLCLM.class);

    protected File pepXmlFile;
    protected File pepXmlDir;

    protected boolean medianCenter = false;
    protected boolean medianCenterByNumCysteines = false;
    protected boolean stripQuantMissingLightOrHeavy = false;
    protected boolean stripQuantNotInHeavyAcrossAll = false;
    protected boolean adjustQuantZeroAreas = false;
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
                       createFileToReadArgumentDefinition("pepxml", false, "PepXML file to process"),
                       createDirectoryToReadArgumentDefinition("pepxmldir", false, "Directory of PepXML files to process"),
                       createBooleanArgumentDefinition("mediancenter", false, "Median-center ratios?", medianCenter),
                       createBooleanArgumentDefinition("bynumcysteines", false,
                               "Median-center ratios separately by number of Cysteines?", medianCenter),
                       createFileToReadArgumentDefinition("strippeptidefile", false,
                               "File containing a list of peptides to strip from results, one per line, all caps"),
                       createFileToWriteArgumentDefinition("out", false, "Output file"),
                       createDirectoryToReadArgumentDefinition("outdir", false, "Output directory"),
                       createBooleanArgumentDefinition("showcharts", false, "Show charts?", showCharts),
                       createIntegerArgumentDefinition("minmediancenterratios", false,
                               "Minimum ratios necessary in order to perform median-centering",
                               minRatiosForMedianCenter),
                       createDecimalArgumentDefinition("minmedianpprophet", false,
                               "Minimum PeptideProphet score to be counted in median calculation",
                               minPeptideProphetForMedian),
                       createBooleanArgumentDefinition("stripquantmissingheavy", false,
                               "Strip quantitation events in which either the heavy " +
                               "isotope was never identified",
                               stripQuantNotInHeavyAcrossAll),
                       createBooleanArgumentDefinition("stripquantmissinglightorheavy", false,
                               "Strip peptides that we haven't seen in both light and heavy states, " +
                               "across all runs.  This ONLY makes sense for multiple metabolically-labeled " +
                               "experiments with a label flip",
                               stripQuantMissingLightOrHeavy),
                       createEnumeratedArgumentDefinition("label", false, labelStrings, labelExplanations,
                               "acrylamide"),
                       createBooleanArgumentDefinition("filterbyproteinprefix", false,
                               "Filter peptides based on prefixes of the protein names that they're associated with?",
                               filterByProteinPrefix),
                       createStringArgumentDefinition("badproteinprefix",false,
                               "Exclude any peptides with any associated proteins with this prefix to their names"),
                       createStringArgumentDefinition("goodproteinprefix",false,
                               "Include any peptides with any associated proteins with this prefix to their names"),
                       createBooleanArgumentDefinition("protprefixexcludequantonly", false,
                               "When excluding peptides based on protein prefix, exclude only quantitation?  " +
                                       "If false, excludes entire ID",
                               excludeProteinPrefixQuantOnly),
                       createDecimalArgumentDefinition("minpprophet", false,
                               "Minimum PeptideProphet score to keep", minPeptideProphet),
                       createDecimalArgumentDefinition("minquantpprophet", false,
                               "Minimum PeptideProphet score for quantitation", minQuantPeptideProphet),
                       createDecimalArgumentDefinition("maxexpect", false,
                               "Maximum expect score to keep", maxExpect),
                       createDecimalArgumentDefinition("maxquantexpect", false,
                               "Maximum expect score for quantitation", maxQuantExpect),
                       createDeltaMassArgumentDefinition("maxfracdeltamass", false,
                               "Maximum fractional delta mass"),
                       createBooleanArgumentDefinition("adjustquantzeroareas", false,
                               "Adjust zero values for light or heavy areas in quantitation (and ratios) to the " +
                               percentileForQuantZeroAreaAdjustment + " percentile of all the (nonzero) values",
                               adjustQuantZeroAreas)
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        pepXmlFile = getFileArgumentValue("pepxml");
        pepXmlDir = getFileArgumentValue("pepxmldir");

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");


        if (pepXmlFile == null)
        {
            if (!hasArgumentValue("pepxmldir"))
                throw new ArgumentValidationException("Either pepxml or pepxmldir is required");
        }
        else
        {
            assertArgumentAbsent("pepxmldir","pepxml");
            assertArgumentPresent("out","pepxml");
        }

        if (pepXmlDir != null)
        {
            assertArgumentAbsent("pepxml","pepxmldir");
            assertArgumentPresent("outdir","pepxmldir");
        }

        medianCenter = getBooleanArgumentValue("mediancenter");
        medianCenterByNumCysteines = getBooleanArgumentValue("bynumcysteines");

        stripQuantMissingLightOrHeavy = getBooleanArgumentValue("stripquantmissinglightorheavy");
        stripQuantNotInHeavyAcrossAll = getBooleanArgumentValue("stripquantmissingheavy");
        adjustQuantZeroAreas = getBooleanArgumentValue("adjustquantzeroareas");


        filterByProteinPrefix = getBooleanArgumentValue("filterbyproteinprefix");

        if (filterByProteinPrefix && !hasArgumentValue("badproteinprefix") && !hasArgumentValue("goodproteinprefix"))
        {
            throw new ArgumentValidationException("If filtering by protein prefix, must specify either a bad or good protein prefix");
        }

        badProteinPrefix = getStringArgumentValue("badproteinprefix");
        goodProteinPrefix = getStringArgumentValue("goodproteinprefix");
        excludeProteinPrefixQuantOnly = getBooleanArgumentValue("protprefixexcludequantonly");

        if (badProteinPrefix != null && goodProteinPrefix != null)
            throw new ArgumentValidationException("Can't have it both ways, sorry");

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
                !medianCenter && !stripQuantMissingLightOrHeavy && !filterByProteinPrefix &&
                !stripQuantNotInHeavyAcrossAll && !adjustQuantZeroAreas &&
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
        if (pepXmlFile != null)
            handleFeatureFile(pepXmlFile, outFile);
        else
        {
            for (File file : pepXmlDir.listFiles())
            {
                if (file.isDirectory() || !file.canRead())
                    continue;
                try
                {
                    File outputFile = new File(outDir, calcOutputFilename(file.getName()));

                    handleFeatureFile(file, outputFile);
                }
                catch (Exception e)
                {
                    ApplicationContext.setMessage("WARNING: Failed to process file " + file.getAbsolutePath() + " as a PepXML file");
                }
            }
        }
        ApplicationContext.setMessage("Done.");
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

            if (stripQuantNotInHeavyAcrossAll)
            {
                while (featureSetIterator.hasNext())
                    addLightHeavyPeptides(featureSetIterator.next());
                featureSetIterator = new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(featureFile);
            }

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

    protected void processFeatureSet(FeatureSet featureSet)
    {
        if (filterByProteinPrefix)
        {
            filterByProteinPrefix(featureSet);
        }

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

        if (stripQuantMissingLightOrHeavy)
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


        if (stripQuantNotInHeavyAcrossAll)
        {
            Set<String> strippedPeptidesThisRun = new HashSet<String>();
            Set<String> allPeptides = new HashSet<String>();

            for (Feature feature : featureSet.getFeatures())
            {
                String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                if (peptide == null)
                    continue;
                allPeptides.add(peptide);

                if (!heavyPeptidesAllRuns.contains(peptide))
                {
                    IsotopicLabelExtraInfoDef.removeRatio(feature);
                    strippedPeptidesThisRun.add(peptide);
                }
            }
            ApplicationContext.infoMessage("\tStripped " + strippedPeptidesThisRun.size() +
                    " peptides (out of " + allPeptides.size() + ") not found in heavy in any run");
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
