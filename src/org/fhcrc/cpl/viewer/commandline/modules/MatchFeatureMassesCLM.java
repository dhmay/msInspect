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
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.PeptideUtilities;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.amt.AmtPeptideEntry;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseFeatureSetGenerator;
import org.apache.log4j.Logger;



import java.io.*;
import java.util.*;


/**
 * Calculates mz values for peptides for a targeted MS/MS experiment
 */
public class MatchFeatureMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MatchFeatureMassesCLM.class);

    protected File[] inputFiles = null;

    protected File outFile = null;
    protected File outDir = null;


    public static final float DEFAULT_DELTA_PPM = 3;
    public static final float DEFAULT_CALIBRATE_DELTA_PPM = 10;


    protected float deltaPPM = DEFAULT_DELTA_PPM;
    protected float calibrateDeltaPPM = DEFAULT_CALIBRATE_DELTA_PPM;

    List<MapWithMass> mapsWithMassesAsc = new ArrayList<MapWithMass>();

    List<String> otherMassFileColNames = new ArrayList<String>();

    protected boolean showCharts = false;
    protected boolean shouldCalibrate = false;

    public MatchFeatureMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matchfeaturemasses";
        mShortDescription = "Match a feature set(s) by mass to a specified list of masses";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "Feature file(s)"),
                        new FileToReadArgumentDefinition("massesfile", true,
                                "File of masses to match, in column 'mass'. " +
                                "Any other columns in the file will be carried forward into output"),
                        new FileToWriteArgumentDefinition("out", false, null),
                        new DirectoryToWriteArgumentDefinition("outdir", false, null),
                        new DecimalArgumentDefinition("deltappm", false, "Mass tolerance (ppm)",
                                DEFAULT_DELTA_PPM),
                        new DecimalArgumentDefinition("calibratedeltappm", false, "Mass tolerance (ppm) for calibration",
                                DEFAULT_CALIBRATE_DELTA_PPM),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", showCharts),
                        new BooleanArgumentDefinition("calibrate", false, 
                            "calibrate masses based on matches? Beware, this will increase bad matches if " +
                            "there isn't a strong signal of good matches. Calibration is very simple, " +
                            "a ppm linear transformation based on median of uncalibrated match delta mass ppm.", false)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inputFiles = getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");
        if (!hasArgumentValue("out") && !hasArgumentValue("outdir"))
            throw new ArgumentValidationException("No output location specified");
        if (hasArgumentValue("out") && hasArgumentValue("outdir"))
            throw new ArgumentValidationException("out and outdir both specified");
        if (inputFiles.length > 1) {
            assertArgumentPresent("outdir");
        }

        deltaPPM = getFloatArgumentValue("deltappm");
        calibrateDeltaPPM = getFloatArgumentValue("calibratedeltappm");


        File massesFile = getFileArgumentValue("massesfile");
        try {
            TabLoader loader = new TabLoader(massesFile);

            boolean foundMass = false;
            for (TabLoader.ColumnDescriptor column : loader.getColumns()) {
                if ("mass".equals(column.name))
                    foundMass = true;
                else otherMassFileColNames.add(column.name);
            }
            if (!foundMass) throw new ArgumentValidationException("masses file doesn't have 'mass' column.");

            for (Object mapObj : loader.load())
                mapsWithMassesAsc.add(new MapWithMass( (Map<String, Object>) mapObj));
        } catch (IOException e) {
            throw new ArgumentValidationException("Failed to load file " + massesFile.getAbsolutePath(),e);
        }
        Collections.sort(mapsWithMassesAsc, new MapWithMassComparatorAsc());

        _log.info("Loaded " + mapsWithMassesAsc.size() + ". Min=" + mapsWithMassesAsc.get(0).mass +
                ", Max=" + mapsWithMassesAsc.get(mapsWithMassesAsc.size()-1).mass);

        showCharts = getBooleanArgumentValue("showcharts");
        shouldCalibrate = getBooleanArgumentValue("calibrate");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File file : inputFiles)
            processFile(file);
    }

    public void processFile(File inputFile) throws CommandLineModuleExecutionException {

        _log.info("Processing file " + inputFile.getAbsolutePath());

        FeatureSet featureSet = null;
        try {
            featureSet = new FeatureSet(inputFile);
        }     catch (Exception e) {
            throw new CommandLineModuleExecutionException(e);
        }

        List<Feature> featuresMassAsc = Arrays.asList(featureSet.getFeatures());
        Collections.sort(featuresMassAsc, new Feature.MassAscComparator());


        if (shouldCalibrate) {
            _log.info("Calibrating...");
            Map<MapWithMass, List<Feature>> matchingResult = massMatchFull(featuresMassAsc, mapsWithMassesAsc, calibrateDeltaPPM);
            List<Float> matchDeltaMassesPPM = new ArrayList<Float>();
            for (MapWithMass map : matchingResult.keySet()) {
                for (Feature feature : matchingResult.get(map)) {
                    float deltaMass = feature.mass - map.mass;
                    float deltaMassPPM = MassUtilities.convertDaToPPM(deltaMass, map.mass);
                    matchDeltaMassesPPM.add(deltaMassPPM);
                }
            }
            float meanDeltaMassPPM = (float) BasicStatistics.median(matchDeltaMassesPPM);

            for (Feature feature : featuresMassAsc) {
                float adjustmentThisFeature = MassUtilities.convertPPMToDa(meanDeltaMassPPM, feature.getMass()); 
                feature.setMass(feature.getMass() - adjustmentThisFeature);
                feature.updateMz();
            }
            _log.info("Done calibrating. Changed all feature masses by " + -meanDeltaMassPPM + "ppm.");


            if (showCharts && !matchDeltaMassesPPM.isEmpty()) {
                new PanelWithHistogram(matchDeltaMassesPPM, "deltaMassPPM preCalibrate").displayInTab();
            }
        }

        Map<MapWithMass, List<Feature>> matchingResult = massMatchFull(featuresMassAsc,
               mapsWithMassesAsc, deltaPPM);

        _log.info("Masses matched: " + matchingResult.size());

        File outputFile;
        if (outFile != null)
            outputFile = outFile;
        else {
            outputFile = new File(outDir, inputFile.getName() + ".matched.tsv");
        }

        List<Float> matchDeltaMassesPPM = new ArrayList<Float>();
        List<Float> matchMasses = new ArrayList<Float>();
        List<Float> matchLogInts = new ArrayList<Float>();
        List<Integer> massMatchedFeatureCounts = new ArrayList<Integer>();
        List<Integer> massMatchedChargeStateCounts = new ArrayList<Integer>();
        Map<Feature, Integer> featureMatchcountMap = new HashMap<Feature, Integer>();
        try {
            PrintWriter outPW = new PrintWriter(outputFile);
            StringBuffer headerLineBuf = new StringBuffer("targetmass");
            for (String colName : otherMassFileColNames)
                headerLineBuf.append("\t" + colName);
            headerLineBuf.append("\tdeltamass\tdeltappm\t");
            headerLineBuf.append(Feature.getFeatureHeader(featureSet.getExtraInformationTypesArray()));
            outPW.println(headerLineBuf.toString());
            outPW.flush();

            for (MapWithMass map : matchingResult.keySet()) {
                massMatchedFeatureCounts.add(matchingResult.get(map).size());
                Set<Integer> chargeStatesMatched = new HashSet<Integer>();
                for (Feature feature : matchingResult.get(map)) {
                    float deltaMass = feature.mass - map.mass;
                    float deltaMassPPM = MassUtilities.convertDaToPPM(deltaMass, map.mass);
                    matchDeltaMassesPPM.add(deltaMassPPM);
                    matchMasses.add(map.mass);
                    matchLogInts.add((float) Math.log(feature.intensity));
                    chargeStatesMatched.add(feature.charge);

                    int matchCountThisFeature = 0;
                    if (featureMatchcountMap.containsKey(feature))
                        matchCountThisFeature = featureMatchcountMap.get(feature);
                    featureMatchcountMap.put(feature, matchCountThisFeature + 1);
                    List<String> chunks = new ArrayList<String>();
                    chunks.add("" + map.mass);
                    for (String otherCol : otherMassFileColNames)
                        chunks.add("" + map.map.get(otherCol));
                    chunks.add("" + deltaMass);
                    chunks.add("" + deltaMassPPM);
                    String lineStart = MS2ExtraInfoDef.convertStringListToString(chunks, "\t");
                    String line = lineStart + "\t" + feature.toString(featureSet.getExtraInformationTypesArray());

                    outPW.println(line);
                    outPW.flush();
                }
                massMatchedChargeStateCounts.add(chargeStatesMatched.size());
            }
            outPW.close();
            _log.info("Wrote output file " + outputFile.getAbsolutePath());
        }
        catch (Exception e) {
            throw new CommandLineModuleExecutionException("Failed to write output file " + outFile.getAbsolutePath());
        }

        if (showCharts && !matchDeltaMassesPPM.isEmpty()) {
            new PanelWithHistogram(matchDeltaMassesPPM, "deltaMassPPM").displayInTab();
            new PanelWithScatterPlot(matchMasses, matchDeltaMassesPPM, "mass_v_deltappm").displayInTab();
            new PanelWithHistogram(new ArrayList<Integer>(featureMatchcountMap.values()), "MatchesPerFeature").displayInTab();

            new PanelWithHistogram(new ArrayList<Integer>(massMatchedFeatureCounts), "MatchesPerMass").displayInTab();
            new PanelWithHistogram(new ArrayList<Integer>(massMatchedChargeStateCounts), "ChargeStatesPerMass").displayInTab();

            new PanelWithScatterPlot(matchDeltaMassesPPM, matchLogInts, "deltappm_v_logint").displayInTab();

        }
    }

    protected static class MapWithMass {
        public float mass;
        public Map<String, Object> map;

        public MapWithMass(float mass, Map<String, Object> map) {
            this.mass = mass;
            this.map = map;
        }

        public MapWithMass(Map<String, Object> map) {
            float mass = ((Double) map.get("mass")).floatValue();
            this.mass = mass;
            this.map = map;
        }
    }

    public static class MapWithMassComparatorAsc implements Comparator<MapWithMass>
    {
        public int compare(MapWithMass o1, MapWithMass o2)
        {
            return o1.mass == o2.mass ? 0 : o1.mass < o2.mass ? -1 : 1;

        }
    }


    /**
     * This version takes a map from formulas to adducts
     * @param featuresSortedMassAsc
     * * @return
     */
    public Map<MapWithMass, List<Feature>> massMatchFull(List<Feature> featuresSortedMassAsc,
                                                              List<MapWithMass> mapsSortedMassAsc,
                                                              float deltaPPMThisMatch)
    {

        int minPossibleMassIndex = 0;

        Map<MapWithMass, List<Feature>> result =
                new HashMap<MapWithMass, List<Feature>>();

        for (int featureIndex = 0; featureIndex<featuresSortedMassAsc.size(); featureIndex++)
        {
            Feature feature = featuresSortedMassAsc.get(featureIndex);
            //not using calculated feature mass, because that presumes M+H and subtracts the H
            float featureMass = feature.getMass();
            float massToleranceDa = MassUtilities.calculateAbsoluteDeltaMass(featureMass, deltaPPMThisMatch,
                    FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
            float minMatchMass = featureMass - massToleranceDa;
            float maxMatchMass = featureMass + massToleranceDa;

            while (mapsSortedMassAsc.get(minPossibleMassIndex).mass < minMatchMass &&
                   minPossibleMassIndex < mapsSortedMassAsc.size()-1)
            {
                minPossibleMassIndex++;
            }


            for (int i=minPossibleMassIndex ; i< mapsSortedMassAsc.size(); i++)
            {
                MapWithMass mapWithMass = mapsSortedMassAsc.get(i);
                float mass = mapsSortedMassAsc.get(i).mass;

                if (mass < minMatchMass)
                    continue;
                else if (maxMatchMass < mass)
                    break;
                else
                {
                    List<Feature> featuresMatchedThisMap = result.get(mapWithMass);
                    if (featuresMatchedThisMap == null)
                    {
                        featuresMatchedThisMap = new ArrayList<Feature>();
                        result.put(mapWithMass, featuresMatchedThisMap);
                    }
                    featuresMatchedThisMap.add(feature);
                }
            }

        }
        return result;
    }

}
