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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToWriteArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.DecimalArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBoxAndWhiskerChart;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;


/**
 * test
 */
public class PeptideRatioVariationCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptideRatioVariationCLM.class);

    protected File[] featureFiles;
    protected File outFile;
    protected float minPeptideProphet = .9f;


    public PeptideRatioVariationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "peptideratiovariation";

        mHelpMessage ="peptideratiovariation";
        mShortDescription = "peptideratiovariation";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedSeriesFileArgumentDefinition(true, "input files"),
                       new FileToWriteArgumentDefinition("out", true, "output file"),
                       new DecimalArgumentDefinition("minpprophet", false, "min peptideprophet", minPeptideProphet)
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        outFile = this.getFileArgumentValue("out");

        minPeptideProphet = getFloatArgumentValue("minpprophet");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Map<String, List<List<Double>>> peptidePerFileLogRatiosMap = new HashMap<String, List<List<Double>>>();

        try
        {
            for (File featureFile : featureFiles)
            {
                FeatureSet featureSet = new FeatureSet(featureFile);

                Map<String, List<Double>> peptideLogRatiosMapThisFile = new HashMap<String, List<Double>>();
                for (Feature feature : featureSet.getFeatures())
                {
                    if (MS2ExtraInfoDef.hasPeptideProphet(feature) &&
                            MS2ExtraInfoDef.getPeptideProphet(feature) >= minPeptideProphet &&
                            IsotopicLabelExtraInfoDef.hasRatio(feature))
                    {
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                        List<Double> thisPeptideLogRatios = peptideLogRatiosMapThisFile.get(peptide);
                        if (thisPeptideLogRatios == null)
                        {
                            thisPeptideLogRatios = new ArrayList<Double>();
                            peptideLogRatiosMapThisFile.put(peptide, thisPeptideLogRatios);
                        }
                        thisPeptideLogRatios.add( Math.log(IsotopicLabelExtraInfoDef.getRatio(feature)));
                    }
                }

                for (String peptide : peptideLogRatiosMapThisFile.keySet())
                {
                    List<List<Double>> peptideRatioLists = peptidePerFileLogRatiosMap.get(peptide);
                    if (peptideRatioLists == null)
                    {
                        peptideRatioLists = new ArrayList<List<Double>>();
                        peptidePerFileLogRatiosMap.put(peptide, peptideRatioLists);
                    }
                    peptideRatioLists.add(peptideLogRatiosMapThisFile.get(peptide));
                }
            }

            Map<String, List<Double>> peptideMeanLogRatios = new HashMap<String, List<Double>>();
            Map<String, Integer> peptideIdCountMap = new HashMap<String, Integer>();
            for (String peptide : peptidePerFileLogRatiosMap.keySet())
            {
                List<Double> thisPeptideLogRatios = new ArrayList<Double>();
                int idCount = 0;
                for (List<Double> perFileLogRatios : peptidePerFileLogRatiosMap.get(peptide))
                {
                    thisPeptideLogRatios.add(BasicStatistics.mean(perFileLogRatios));
                    idCount += perFileLogRatios.size();
                }
                peptideMeanLogRatios.put(peptide, thisPeptideLogRatios);
                peptideIdCountMap.put(peptide, idCount);
            }

            Map<Integer, List<Float>> perFractionCountStandardDevs = new HashMap<Integer, List<Float>>();
            List<Float> allPeptideStandardDevs = new ArrayList<Float>();
            Map<String, Float> peptideStdDevMap = new HashMap<String, Float>();
            Map<String, Float> peptideUberMeanMap = new HashMap<String, Float>();


            for (String peptide : peptideMeanLogRatios.keySet())
            {
                List<Double> logRatioMeans = peptideMeanLogRatios.get(peptide);
                int fractionCount = logRatioMeans.size();

                double stdDevLogRatios = BasicStatistics.standardDeviation(logRatioMeans);
                allPeptideStandardDevs.add((float) stdDevLogRatios);
                peptideStdDevMap.put(peptide, (float) stdDevLogRatios);
                List<Float> thisNumFractionsDevList = perFractionCountStandardDevs.get(fractionCount);
                if (thisNumFractionsDevList == null)
                {
                    thisNumFractionsDevList = new ArrayList<Float>();
                     perFractionCountStandardDevs.put(fractionCount,thisNumFractionsDevList );
                }
                thisNumFractionsDevList.add((float) stdDevLogRatios);

                peptideUberMeanMap.put(peptide, (float)BasicStatistics.mean(logRatioMeans));
            }

            PanelWithHistogram pwh = new PanelWithHistogram(allPeptideStandardDevs, "Fraction Std Devs", 200);
            pwh.displayInTab();

            PanelWithBoxAndWhiskerChart boxWhiskers = new PanelWithBoxAndWhiskerChart("Std Devs Per Frac Count");
            for (int i=0; i<featureFiles.length; i++)
            {
                List<Float> thisNumFractionsStdDevs = perFractionCountStandardDevs.get(i);
                if (thisNumFractionsStdDevs != null)
                {
                    double[] asDouble = new double[thisNumFractionsStdDevs.size()];
                    for (int j=0; j<thisNumFractionsStdDevs.size(); j++)
                        asDouble[j] = thisNumFractionsStdDevs.get(j);
                    boxWhiskers.addData(asDouble, "" + (i));
                }
                boxWhiskers.displayInTab();
            }


            for (int i=0; i<featureFiles.length; i++)
            {
                List<Float> thisNumFractionsStdDevs = perFractionCountStandardDevs.get(i);
                if (thisNumFractionsStdDevs != null)
                {
                    PanelWithHistogram pwhi = new PanelWithHistogram(thisNumFractionsStdDevs,
                            "Std Devs " + i + " fracs", 200);
                    pwhi.displayInTab();
                }
            }


            PrintWriter pw = new PrintWriter(outFile);
            pw.println("peptide\tmean_frac_log_ratio\tstddev_frac_log_ratio\tnum_fractions\tnum_times_idd");
            for (String peptide : peptidePerFileLogRatiosMap.keySet())
            {
                pw.println(peptide + "\t" + peptideUberMeanMap.get(peptide) + "\t" +
                        peptideStdDevMap.get(peptide) + "\t" + peptideIdCountMap.get(peptide));
                pw.flush();
            }
            pw.close();

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
