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
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.apache.log4j.Logger;


import java.util.*;
import java.io.File;
import java.io.IOException;


public class CaseControlRatioFeatureSetCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CaseControlRatioFeatureSetCLM.class);

      protected File file;
    protected File outFile;
    protected int minRunsPerGroup = 1;
    Object[] arrayRows;
    PeptideArrayAnalyzer peptideArrayAnalyzer = null;

    protected float minGreaterIntensityForRatio = 0f;

    protected File fastaFile;

    protected boolean shouldAddRatiosForMissing = false;
    protected int minPresentToAddMissingRatio = 2;

    protected boolean winsorize = false;

    //90% Winsorization means locking values to 5th percentile and 95th percentile
    protected float winsorizePercent = 96;



    public CaseControlRatioFeatureSetCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "casecontrolratiofeatureset";
        mShortDescription = "Create a featureset with one feature per peptide array row with a peptide assignment " +
                "and ratios for those rows with both case and control features";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true, "peptide array"),
                        new FileToReadArgumentDefinition("caserunlistfile",false,
                                "File containing the names of runs in the case group, one per line"),
                        new FileToReadArgumentDefinition("controlrunlistfile",false,
                                "File containing the names of runs in the control group, one per line"),


                        new IntegerArgumentDefinition("minrunspergroup",false,
                                "Minimum number of runs required for a feature to be included in the feature set",
                                minRunsPerGroup),
                        new FileToWriteArgumentDefinition("out", true, "Output pepXML file"),
                        new DecimalArgumentDefinition("mingreaterintensity", false,
                                "Minimum intensity for the /greater/ of the two mean intensity values, case or control",
                                minGreaterIntensityForRatio),
                        new FileToReadArgumentDefinition("fasta", false, "FASTA database to pretend the file came from"),
                        new BooleanArgumentDefinition("addmissing", false, "Add extreme ratios for peptides with at " +
                                "least mintoaddmissing values in one group and none in the other?", shouldAddRatiosForMissing),
                        new IntegerArgumentDefinition("mintoaddmissing", false, "if addmissing specified, minimum " +
                                "number of intensity values in the present group in order to add an extreme ratio",
                                minPresentToAddMissingRatio),
                        new BooleanArgumentDefinition("winsorize", false, "Winsorize missing ratios at 95% (default: " +
                                "use 0 and 999)", winsorize)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        file = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        fastaFile = getFileArgumentValue("fasta");
        outFile = getFileArgumentValue("out");                                                       

        minRunsPerGroup = getIntegerArgumentValue("minrunspergroup");

        minGreaterIntensityForRatio = getFloatArgumentValue("mingreaterintensity");

        minPresentToAddMissingRatio = getIntegerArgumentValue("mintoaddmissing");
        shouldAddRatiosForMissing = getBooleanArgumentValue("addmissing");

        winsorize = getBooleanArgumentValue("winsorize");
        if (!shouldAddRatiosForMissing) assertArgumentAbsent("winsorize", "addmissing");

        try
        {
            peptideArrayAnalyzer = new PeptideArrayAnalyzer(file);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }

        File caseRunListFile = getFileArgumentValue("caserunlistfile");
        File controlRunListFile = getFileArgumentValue("controlrunlistfile");

        if (peptideArrayAnalyzer.getRunNames().size() > 2)
        {
            assertArgumentPresent("caserunlistfile");
            assertArgumentPresent("controlrunlistfile");
        }
        else
        {
            if (caseRunListFile == null)
            {
                ApplicationContext.infoMessage("Only two runs this file.  First will be 'case', second 'control'");
                peptideArrayAnalyzer.setCaseRunNames(new String[] { peptideArrayAnalyzer.getRunNames().get(0) });
                peptideArrayAnalyzer.setControlRunNames(new String[] { peptideArrayAnalyzer.getRunNames().get(1) });
            }
        }

        if (caseRunListFile != null)
        {
            try
            {
                peptideArrayAnalyzer.loadCaseControlRunListsFromFiles(caseRunListFile, controlRunListFile);
            }
            catch (IOException e)
            {
                throw new ArgumentValidationException("Failed to load case or control run list file", e);
            }

            ApplicationContext.infoMessage("Case runs:");
            for (String caseRun : peptideArrayAnalyzer.getCaseRunNames())
                ApplicationContext.setMessage("\t" + caseRun);
            ApplicationContext.infoMessage("Control runs:");
            for (String controlRun : peptideArrayAnalyzer.getControlRunNames())
                ApplicationContext.setMessage("\t" + controlRun);
        }
        else
        {
            peptideArrayAnalyzer.setCaseRunNames(new String[] { peptideArrayAnalyzer.getRunNames().get(0)});
            peptideArrayAnalyzer.setControlRunNames(new String[] { peptideArrayAnalyzer.getRunNames().get(1)});
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<Feature> features = new ArrayList<Feature>();
        //todo: un-hardcode
        AnalyzeICAT.IsotopicLabel label = new AnalyzeICAT.IsotopicLabel("Acrylamide (D0/D3)",71.03657F,3.0188F,'C',
                    AnalyzeICAT.DEFAULT_MAX_LABEL_COUNT);
        int numFeaturesWithRatios = 0;
        int numFeaturesRatiosInsufficientIntensity = 0;
        int numFeaturesWithPeptides = 0;
        Set<String> allPeptides = new HashSet<String>();

        List<Feature> extremeRatioFeatures = new ArrayList<Feature>();
        List<Float> nonMissingRatios = new ArrayList<Float>();
        for (Map<String, Object> rowMap : peptideArrayAnalyzer.getRowMaps())
        {
             Set<String> rowPeptides = peptideArrayAnalyzer.getAllRowPeptides(rowMap);
            if (rowPeptides.size() != 1)
                continue;

            List<Double> caseIntensities = new ArrayList<Double>();
            for (String caseRun : peptideArrayAnalyzer.getCaseRunNames())
            {
                Double intensity = peptideArrayAnalyzer.getRunIntensity(rowMap, caseRun);
                if (intensity != null)
                    caseIntensities.add(intensity);
            }

            List<Double> controlIntensities = new ArrayList<Double>();
            for (String run : peptideArrayAnalyzer.getControlRunNames())
            {
                Double intensity = peptideArrayAnalyzer.getRunIntensity(rowMap, run);
                if (intensity != null)
                    controlIntensities.add(intensity);
            }
            if (caseIntensities.size() == 0 && controlIntensities.size() == 0)
                continue;

            numFeaturesWithPeptides++;
            Feature rowFeature = new Feature();
            String peptide = rowPeptides.iterator().next();
            allPeptides.add(peptide);
            rowFeature.setCharge(Integer.parseInt(rowMap.get("charge").toString()));
            if (rowMap.containsKey("minScan"))
                rowFeature.setScan((int) Double.parseDouble(rowMap.get("minScan").toString()));
            else
            {
                rowFeature.setScan((int) Double.parseDouble(rowMap.get("minTime").toString()));                
            }
            MS2ExtraInfoDef.setPeptideList(rowFeature, peptide);
            MS2ExtraInfoDef.setPeptideProphet(rowFeature, 0.95);
            Protein dummyProtein = new Protein("asdf", peptide.getBytes());
            rowFeature.setMass((float) new Peptide(dummyProtein, 0, peptide.length()).getMonoisotopicMass());
            rowFeature.updateMz();            

//if (peptide.equals("FLGVAEQLHNEGFK")) System.err.println("FLGVAEQLHNEGFK, " + caseIntensities.size() +", " +  controlIntensities.size());
            if (caseIntensities.size() >= minRunsPerGroup && controlIntensities.size() >= minRunsPerGroup)
            {
                double geoMeanCase = BasicStatistics.geometricMean(caseIntensities);
                double geoMeanControl = BasicStatistics.geometricMean(controlIntensities);
                //only add a ratio if at least one of the intensities being compared is above threshold
                if (Math.max(geoMeanCase, geoMeanControl) > minGreaterIntensityForRatio)
                {
                    double ratio = geoMeanCase / geoMeanControl;
                    IsotopicLabelExtraInfoDef.setRatio(rowFeature, ratio);
                    IsotopicLabelExtraInfoDef.setLabel(rowFeature, label);
                    IsotopicLabelExtraInfoDef.setHeavyFirstScan(rowFeature, rowFeature.getScan());
                    IsotopicLabelExtraInfoDef.setHeavyLastScan(rowFeature, rowFeature.getScan());
                    IsotopicLabelExtraInfoDef.setLightFirstScan(rowFeature, rowFeature.getScan());
                    IsotopicLabelExtraInfoDef.setLightLastScan(rowFeature, rowFeature.getScan());
                    IsotopicLabelExtraInfoDef.setLightIntensity(rowFeature, geoMeanCase);
                    IsotopicLabelExtraInfoDef.setHeavyIntensity(rowFeature, geoMeanControl);
                    nonMissingRatios.add((float) ratio);
                    numFeaturesWithRatios++;
//if (peptide.equals("FLGVAEQLHNEGFK")) System.err.println("\tratio: " + ratio);                    
                }
                else
                    numFeaturesRatiosInsufficientIntensity++;
            }
            else
            {
                //check if we should add a fake ratio for missing data in one group or the other.
                //if case missing, use 0.  If control missing, use 999
                if (shouldAddRatiosForMissing)
                {
                    if ((caseIntensities.size() >= minPresentToAddMissingRatio && controlIntensities.size() == 0) ||
                        (controlIntensities.size() >= minPresentToAddMissingRatio && caseIntensities.size() == 0))
                    {

                        IsotopicLabelExtraInfoDef.setLabel(rowFeature, label);
                        IsotopicLabelExtraInfoDef.setHeavyFirstScan(rowFeature, rowFeature.getScan());
                        IsotopicLabelExtraInfoDef.setHeavyLastScan(rowFeature, rowFeature.getScan());
                        IsotopicLabelExtraInfoDef.setLightFirstScan(rowFeature, rowFeature.getScan());
                        IsotopicLabelExtraInfoDef.setLightLastScan(rowFeature, rowFeature.getScan());
                        numFeaturesWithRatios++;

                        if (caseIntensities.size() == 0)
                        {
                            IsotopicLabelExtraInfoDef.setRatio(rowFeature, 0);
                            IsotopicLabelExtraInfoDef.setLightIntensity(rowFeature, 0);
                            IsotopicLabelExtraInfoDef.setHeavyIntensity(rowFeature, BasicStatistics.geometricMean(controlIntensities));
                        }
                        else
                        {
                            IsotopicLabelExtraInfoDef.setRatio(rowFeature, 999);
                            IsotopicLabelExtraInfoDef.setLightIntensity(rowFeature, BasicStatistics.geometricMean(caseIntensities));
                            IsotopicLabelExtraInfoDef.setHeavyIntensity(rowFeature, 0);
                        }

                    }
                    extremeRatioFeatures.add(rowFeature);
                }
            }
            features.add(rowFeature);

        }

        if (winsorize)
        {
            float winsorMinRatio = (float) BasicStatistics.percentile(nonMissingRatios, (100-winsorizePercent)/2);
            float winsorMaxRatio = (float) BasicStatistics.percentile(nonMissingRatios, winsorizePercent + (100-winsorizePercent)/2);
            ApplicationContext.infoMessage("Winsorizing: min=" + winsorMinRatio + ", max=" + winsorMaxRatio);
            List<Float> logRatiosBeforeWinsor = new ArrayList<Float>();
            for (float ratio : nonMissingRatios)
                logRatiosBeforeWinsor.add((float) Math.log(ratio));
            new PanelWithHistogram(logRatiosBeforeWinsor, "LogRatio preWinsor", 100).displayInTab();

            List<Float> logRatiosAfterWinsor = new ArrayList<Float>();
            for (Feature feature : features)
            {
                if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                {
                    float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
                    float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature);
                    float heavyIntensity = (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
// if (lightIntensity == 0 && heavyIntensity == 0)
//             throw new RuntimeException ("light & heavy 0");
                    if (ratio < winsorMinRatio)
                    {
                        ratio = winsorMinRatio;
                        lightIntensity = ratio * heavyIntensity;

                    }
                    if (ratio > winsorMaxRatio)
                    {
                        ratio = winsorMaxRatio;
                        heavyIntensity = lightIntensity / ratio;

                    }
                    logRatiosAfterWinsor.add((float) Math.log(ratio));
                    IsotopicLabelExtraInfoDef.setRatio(feature, ratio);
                    IsotopicLabelExtraInfoDef.setLightIntensity(feature, lightIntensity);
                    IsotopicLabelExtraInfoDef.setHeavyIntensity(feature, heavyIntensity);

                }
            }
            new PanelWithHistogram(logRatiosAfterWinsor, "LogRatio postWinsor", 100).displayInTab();

        }
        ApplicationContext.infoMessage("Rows with one peptide: " + numFeaturesWithPeptides);
        ApplicationContext.infoMessage("Unique peptides: " + allPeptides.size());
        ApplicationContext.infoMessage("Ratios created for " + numFeaturesWithRatios + " features.  Without: " + (numFeaturesWithPeptides - numFeaturesWithRatios));
        ApplicationContext.infoMessage("Features with intensities too low for ratio: " + numFeaturesRatiosInsufficientIntensity);
        FeatureSet featureSet = new FeatureSet(features.toArray(new Feature[features.size()]));
        if (fastaFile != null)
            MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSet, fastaFile.getAbsolutePath());
        try
        {
            featureSet.savePepXml(outFile);
            ApplicationContext.infoMessage("Saved feature file " + outFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to save featureset",e);
        }
    }


}
