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
package org.fhcrc.cpl.viewer.amt.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.HistogramDataset;

import javax.imageio.ImageIO;
import java.util.*;
import java.io.File;
import java.io.PrintWriter;


/**
 * Command line module for getting information and charts on an AMT database.
 * Has an extra mode for saving all charts to a directory
 */
public class AmtDatabaseDiagnosticCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AmtDatabaseDiagnosticCommandLineModule.class);

    protected File outFile;
    protected File outDir;

    protected PrintWriter outPW;

    protected int mode=-1;

    protected AmtDatabase amtDB;
    protected boolean showCharts;

    protected FeatureSet ms1FeatureSet;
    
    protected String peptide = null;

    protected final static String[] modeStrings = {"compareobshcalch",
             "comparecalchtime","plotthmaps","histhydrostddev",
             "plotrunstddevs",
             "saveallcharts","plotmassvsobshydro","basicinfo",
             "histmasses","plotrunhvsmedian","peptidedetails","plotmassh",
             "histidprobs"};

    protected static final String[] modeExplanations =
            {
                    "Compare Observed vs. Calculated H, per run",
                    "Compare Calculated H vs. Time, per run",
                    "Plot the T->H mappings for each run",
                    "Histogram of the hydrophobicity standard deviations across all observations for each peptide",
                    "Bar chart of the mean difference between predicted and observed H, per run",
                    "Save all charts to image files in a specified directory",
                    "Plot peptide mass vs. observed hydrophobicity, per run",
                    "Return basic AMT database information",
                    "Histogram all peptide masses",
                    "Plot observed hydrophobicity against median peptide hydrophobicity, by run",
                    "Show details about a single peptide",
                    "Plot mass and H for every entry in the DB",
                    "Histogram the probabilities of all peptide entry IDs' correctness",
            };

    protected static final int MODE_COMPARE_OBS_H_CALC_H = 0;
    protected static final int MODE_COMPARE_CALC_H_TIME = 1;
    protected static final int MODE_PLOT_T_H_MAPS = 2;
    protected static final int MODE_HIST_HYDRO_STD_DEV = 3;
    protected static final int MODE_PLOT_RUN_HYDRO_DIFFS = 4;
    protected static final int MODE_SAVE_ALL_CHARTS=5;
    protected static final int MODE_PLOT_MASS_VS_OBS_HYDRO=6;
    protected static final int MODE_BASIC_INFO=7;
    protected static final int MODE_HIST_MASSES=8;
    protected static final int MODE_PLOT_RUN_H_VS_MEDIAN=9;
    protected static final int MODE_PEPTIDE_DETAILS=10;
    protected static final int MODE_PLOT_MASS_H=11;
    protected static final int MODE_HIST_ID_PROBS=12;




    public AmtDatabaseDiagnosticCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "amtdiagnostic";
        mShortDescription = "Get various kinds of information about an AMT database.";
        mHelpMessage = "Get various kinds of information (mostly in chart form) about an AMT database.  The 'saveallcharts' mode will generate all possible charts and save them to a directory";

        CommandLineArgumentDefinition[] argDefs =
               {
                    new EnumeratedValuesArgumentDefinition("mode",true,modeStrings, modeExplanations),
                    createUnnamedFileArgumentDefinition(true, "AMT database file"),
                    new FileToWriteArgumentDefinition("out",false, "output filepath (for individual charts)"),
                    new DirectoryToWriteArgumentDefinition("outdir",false, "output directory (for saving all charts)"),
                    new BooleanArgumentDefinition("showcharts",false, "show charts?", false),
                    new StringArgumentDefinition("peptide",false,"Peptide to get details about (for mode peptidedetails only)"),
                    new FeatureFileArgumentDefinition("ms1features", false, "MS1 feature file to show along with database entries (plotmassandh mode only)"),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        File dbFile = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        showCharts = getBooleanArgumentValue("showcharts");
        

        switch(mode)
        {
            case MODE_SAVE_ALL_CHARTS:
                assertArgumentPresent("outdir");
                showCharts = false;
                break;
            case MODE_PEPTIDE_DETAILS:
                assertArgumentPresent("peptide");
                peptide = getStringArgumentValue("peptide");
                break;
            default:
                showCharts = true;
                break;
        }


        try
        {
            amtDB = new AmtXmlReader(dbFile).getDatabase();
        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Failure!! " + e.getMessage() + ", type: " + e.getClass().getName());
            e.printStackTrace(System.err);
            throw new ArgumentValidationException(e);
        }

        ms1FeatureSet = getFeatureSetArgumentValue("ms1features");
    }


    /**
     * a great big switch statement
     */
    public void execute() throws CommandLineModuleExecutionException
    {

        try
        {
            ApplicationContext.infoMessage(amtDB.toString());
              JFreeChart chart = null;
            switch(mode)
            {
                case  MODE_COMPARE_OBS_H_CALC_H:
                    chart = compareObservedHAgainstCalculatedH();
                    break;
                case MODE_COMPARE_CALC_H_TIME:
                    chart = compareCalculatedHAgainstTime();
                    break;
                case MODE_PLOT_T_H_MAPS:
                    chart = showTimeHydMaps();
                    break;
                case MODE_HIST_HYDRO_STD_DEV:
                    chart = histogramHydroStdDev();
                    break;
                case MODE_PLOT_RUN_HYDRO_DIFFS:
                    chart = showMeanHydroDiffByRun();
                    break;
                case MODE_SAVE_ALL_CHARTS:

                    Map<String,JFreeChart> filenameChartMap =
                            new HashMap<String,JFreeChart>();
                    filenameChartMap.put("obsH_calcH.png",compareObservedHAgainstCalculatedH());

                    filenameChartMap.put("calcH_time.png",compareCalculatedHAgainstTime());
                    filenameChartMap.put("hist_peptide_hydro_stddev.png",histogramHydroStdDev());
                    filenameChartMap.put("mean_diff_by_run.png",showMeanHydroDiffByRun());
                    filenameChartMap.put("mass_obs_hydro.png",plotMassVersusObsHydroByRun());
                    filenameChartMap.put("run_h_vs_median.png", plotRunHvsMedian());
                    filenameChartMap.put("mass_h_all_entries.png", plotMassAndH());
                    for (String filename : filenameChartMap.keySet())
                    {
                        chart = filenameChartMap.get(filename);
                        File outChartFile = new File(outDir.getAbsolutePath() + File.separatorChar + filename);
                        ImageIO.write(chart.createBufferedImage(500,400),"png",outChartFile);
                        ApplicationContext.infoMessage("Wrote image file " + filename);
                    }
                    ApplicationContext.infoMessage("Done.");
                    break;
                case MODE_PLOT_MASS_VS_OBS_HYDRO:
                    chart = plotMassVersusObsHydroByRun();
                    break;
                case MODE_BASIC_INFO:
                    basicInfo();
                    break;
                case MODE_HIST_MASSES:
                    chart = histogramMasses();
                    break;
                case MODE_PLOT_RUN_H_VS_MEDIAN:
                    chart = plotRunHvsMedian();
                    break;
                case MODE_PEPTIDE_DETAILS:
                    peptideDetails();
                    break;
                case MODE_PLOT_MASS_H:
                    chart = plotMassAndH();
                    break;
                case MODE_HIST_ID_PROBS:
                    chart = histogramIDProbs();
                    break;
            }
            if (chart != null && showCharts)
            {
                ChartDialog cd = new ChartDialog(chart);
                cd.setVisible(true);
            }


        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

    }

    protected void basicInfo()
    {
        ApplicationContext.infoMessage("Run\tPeptides");
        int sumPeptidesPerRun = 0;
        int minPeptidesPerRun = Integer.MAX_VALUE;
        int maxPeptidesPerRun = 0;

        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            int numPeptidesThisRun = amtDB.getPeptideEntriesForRun(runEntry).length;
            String runFilename = runEntry.getPepXmlFilename();
            ApplicationContext.infoMessage(runFilename + "\t" + numPeptidesThisRun);
            sumPeptidesPerRun += numPeptidesThisRun;
            if (numPeptidesThisRun < minPeptidesPerRun)
                minPeptidesPerRun = numPeptidesThisRun;
            if (numPeptidesThisRun > maxPeptidesPerRun)
                maxPeptidesPerRun = numPeptidesThisRun;
        }
        ApplicationContext.infoMessage(amtDB.toString());

        ApplicationContext.infoMessage("\nTotal\t" + amtDB.getEntries().length);
        ApplicationContext.infoMessage("Mean Peptides Per Run: " +
                                       (double) sumPeptidesPerRun / amtDB.numRuns() +
                                       ", Min: " + minPeptidesPerRun + ", Max: " + maxPeptidesPerRun);
        int sumModsPerPeptide = 0;
        int numTotalObservations = 0;

        for (AmtPeptideEntry peptideEntry : amtDB.getEntries())
        {
            sumModsPerPeptide += peptideEntry.getNumModificationStates();
            numTotalObservations += peptideEntry.getNumObservations();
        }
        ApplicationContext.infoMessage("Mean modification states per peptide: " +
                (((double) sumModsPerPeptide) / (double) amtDB.numEntries()));
        ApplicationContext.infoMessage("Total number of peptide observations: " + numTotalObservations);
    }

    protected void peptideDetails()
    {
        ApplicationContext.infoMessage("Peptide: " + peptide);
        AmtPeptideEntry entry = amtDB.getEntry(peptide);
        if (entry == null)
        {
            ApplicationContext.infoMessage("Not in database.");
            return;
        }

        ApplicationContext.infoMessage("Predicted H: " + entry.getPredictedHydrophobicity());        
        ApplicationContext.infoMessage("Median observed H: " + entry.getMedianObservedHydrophobicity());
        ApplicationContext.infoMessage("Total observations: " + entry.getNumObservations());
        double[] obsHs = new double[entry.getNumObservations()];
        ApplicationContext.infoMessage("H standard deviation: " +
                entry.getHydrophobicityStandardDeviation());

        AmtPeptideEntry.AmtPeptideModificationStateEntry[] modStates =
                entry.getModificationStateEntries();
        ApplicationContext.infoMessage("Modification States: " + modStates.length + "\n");

        for (AmtPeptideEntry.AmtPeptideModificationStateEntry modState : modStates)
        {
            ApplicationContext.infoMessage("Mod State " + modState.getModifiedSequence() + ", median H: " + modState.getMedianObservedHydrophobicity());
        }

    }

    /**
     * Plot observed H values in each individual run vs. the median H value
     * @return
     */
    protected JFreeChart plotRunHvsMedian()
    {
        PanelWithScatterPlot spd =  new PanelWithScatterPlot();
        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            List<Double> medianValues = new ArrayList<Double>();
            List<Double> runValues = new ArrayList<Double>();

            for (AmtPeptideEntry entry : amtDB.getEntries())
            {
                AmtPeptideEntry.AmtPeptideObservation obs =
                        entry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    medianValues.add(entry.getMedianObservedHydrophobicity());
                    runValues.add(obs.getObservedHydrophobicity());
                }
            }

            double[] medianValuesArray = new double[medianValues.size()];
            double[] runValuesArray = new double[medianValues.size()];
            for (int i = 0; i < medianValuesArray.length; i++)
            {
                medianValuesArray[i] = medianValues.get(i);
                runValuesArray[i] = runValues.get(i);
            }
            spd.addData(medianValuesArray, runValuesArray,
                    "" + amtDB.getSequenceForRun(runEntry));
        }
        if (outPW != null)
            outPW.close();
        spd.setAxisLabels("Median Observed Hydrophobicity", "Run Observed Hydrophobicity");
        if (showCharts)
            spd.displayInTab();
        return spd.getChart();
    }

    /**
     * Show mean hydrophobicity deviation from prediction, on a handy bar chart.
     * Useful for finding runs that are really wacked out.
     * @return
     */
    protected JFreeChart showMeanHydroDiffByRun()
    {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        int i=1;
        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            double sumDifference = 0;
            int numObsThisRun = 0;
            for (AmtPeptideEntry entry : amtDB.getEntries())
            {
                AmtPeptideEntry.AmtPeptideObservation obs =
                        entry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    sumDifference += (entry.getMedianObservedHydrophobicity() - obs.getObservedHydrophobicity());
                    numObsThisRun++;
                }
            }
            sumDifference = sumDifference /numObsThisRun;
            dataset.setValue(sumDifference, "Mean diff", "" +i);
            i++;
        }

        JFreeChart chart = ChartFactory.createBarChart(null, null, null, dataset,
                PlotOrientation.VERTICAL, true, false, false);

        PanelWithChart pwc = new PanelWithChart(chart);
        pwc.setName("mean hydro diff");
        pwc.displayInTab();

        return chart;
    }

    /**
     * Plot the T->H mappings of all runs.  Useful for showing up odd runs
     * @return
     */
    protected JFreeChart showTimeHydMaps()
    {

        double maxTimeAnyRun = 0;
        if (amtDB.numRuns() > 300)
        {
            for (AmtPeptideEntry.AmtPeptideObservation obs : amtDB.getObservationsForRun(amtDB.getRuns()[1]))
                if (obs.getTimeInRun() > maxTimeAnyRun)
                    maxTimeAnyRun = obs.getTimeInRun();
            maxTimeAnyRun += 50;
            ApplicationContext.setMessage("Punting on figuring out max time, too many runs.  Using " + maxTimeAnyRun);
        }
        else
        for (AmtRunEntry run : amtDB.getRuns())
        {
            for (AmtPeptideEntry.AmtPeptideObservation obs : amtDB.getObservationsForRun(run))
                if (obs.getTimeInRun() > maxTimeAnyRun)
                    maxTimeAnyRun = obs.getTimeInRun();
        }

        PanelWithScatterPlot spd = new PanelWithScatterPlot();
        if (amtDB.numRuns() > 50)
            spd.setShowLegend(false);
        spd.setName("Realignment of DB runs");
        int i=0;
        for (AmtRunEntry run : amtDB.getRuns())
        {
            int numDotsOnChart = ((int) maxTimeAnyRun+1) /5;
            double[] mappingXVals = new double[numDotsOnChart];
            double[] mappingYVals = new double[numDotsOnChart];
            for (int j=0; j<numDotsOnChart; j++)
            {
                mappingXVals[j] = 5 * j;
                mappingYVals[j] = RegressionUtilities.mapValueUsingCoefficients(run.getTimeHydMapCoefficients(), j);
            }
            spd.addData(mappingXVals, mappingYVals, "");
            if (amtDB.numRuns() > 100) ApplicationContext.setMessage("Added run " + i++);
        }

//            spd.addData(reverseRegressionLineXVals, reverseRegressionLineYVals, "reverse regression line");
//            spd.addData(regressionLineXVals, regressionLineYVals, "regression line");
//            spd.addData(lowStudResMs1Times, lowStudResMs2H, "low studentized residual");

        spd.setAxisLabels("MS1 time","AMT Hydrophobicity");
        spd.displayInTab();
        return spd.getChart();

    }

    /**
     * Plot peptide mass vs. observed hydrophobicity
     * @return
     */
    protected JFreeChart plotMassVersusObsHydroByRun()
    {
        PanelWithScatterPlot spd =
                new PanelWithScatterPlot();
        spd.setName("mass vs obs h");

        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            List<Double> xValues = new ArrayList<Double>();
            List<Double> yValues = new ArrayList<Double>();

            for (AmtPeptideEntry entry : amtDB.getEntries())
            {
                AmtPeptideEntry.AmtPeptideObservation obs =
                        entry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    xValues.add(entry.getMass());
                    yValues.add(obs.getObservedHydrophobicity());
                }
            }
            double[] xValuesArray = new double[xValues.size()];
            double[] yValuesArray = new double[xValues.size()];
            for (int i = 0; i < xValuesArray.length; i++)
            {
                xValuesArray[i] = xValues.get(i);
                yValuesArray[i] = yValues.get(i);
            }
            spd.addData(xValuesArray, yValuesArray,
                    "" + amtDB.getSequenceForRun(runEntry));
            if (outPW != null)
                for (int i = 0; i < xValuesArray.length; i++)
                {
                    outPW.println(xValuesArray[i] + "\t" + yValuesArray[i] + "\t" + amtDB.getSequenceForRun(runEntry));
                }
        }
        if (outPW != null)
            outPW.close();
        spd.setAxisLabels("Mass", "observed hydro");
        if (showCharts)
            spd.setVisible(true);
        return spd.getChart();

    }

    /**
     * Scatterplot observed hydrophobicity against calculated H, by run
     * @return
     */
    protected JFreeChart compareObservedHAgainstCalculatedH()
    {
        ApplicationContext.infoMessage("Mean of deviations from predicted: " + amtDB.calculateMeanDifferenceFromPredictedHydro());
        ApplicationContext.infoMessage("Standard deviation of deviations from predicted: " + amtDB.calculateStandardDeviationDifferenceFromPredictedHydro());

        PanelWithScatterPlot spd =
                new PanelWithScatterPlot();
        spd.setName("obs vs calc h");
        if (outFile != null)
        {
            try
            {
                outPW = new PrintWriter(outFile);
                outPW.println("predictedH\tobservedH\tfraction");
            }
            catch (Exception e)
            {
                throw new RuntimeException(e);
            }
        }
        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            List<Double> xValues = new ArrayList<Double>();
            List<Double> yValues = new ArrayList<Double>();

            for (AmtPeptideEntry entry : amtDB.getEntries())
            {

                AmtPeptideEntry.AmtPeptideObservation obs =
                        entry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    xValues.add(entry.getPredictedHydrophobicity());
                    yValues.add(obs.getObservedHydrophobicity());
                }
            }
            double[] xValuesArray = new double[xValues.size()];
            double[] yValuesArray = new double[xValues.size()];
            for (int i = 0; i < xValuesArray.length; i++)
            {
                xValuesArray[i] = xValues.get(i);
                yValuesArray[i] = yValues.get(i);
            }
            spd.addData(xValuesArray, yValuesArray,
                    "" + amtDB.getSequenceForRun(runEntry));
            if (outPW != null)
                for (int i = 0; i < xValuesArray.length; i++)
                {
                    outPW.println(xValuesArray[i] + "\t" + yValuesArray[i] + "\t" + amtDB.getSequenceForRun(runEntry));
                }
        }
        if (outPW != null)
            outPW.close();
        spd.setAxisLabels("predicted hydro", "observed hydro");
        if (showCharts)
            spd.setVisible(true);
        return spd.getChart();

    }

    /**
     * Plot time vs. calculated H, per run
     * @return
     */
    protected JFreeChart compareCalculatedHAgainstTime()
    {
        PanelWithScatterPlot spd =
                new PanelWithScatterPlot();
        spd.setName("calc h vs time");

        if (outFile != null)
        {
            try
            {
                outPW = new PrintWriter(outFile);
                outPW.println("timeInRun\tpredictedH\tfraction");
            }
            catch (Exception e)
            {
                throw new RuntimeException(e);
            }
        }
        for (AmtRunEntry runEntry : amtDB.getRuns())
        {
            List<Double> xValues = new ArrayList<Double>();
            List<Double> yValues = new ArrayList<Double>();

            for (AmtPeptideEntry entry : amtDB.getEntries())
            {

                AmtPeptideEntry.AmtPeptideObservation obs =
                        entry.getObservationForRun(runEntry);
                if (obs != null)
                {
                    yValues.add(entry.getPredictedHydrophobicity());
                    xValues.add(obs.getTimeInRun());
                }
            }
            double[] xValuesArray = new double[xValues.size()];
            double[] yValuesArray = new double[xValues.size()];
            for (int i = 0; i < xValuesArray.length; i++)
            {
                xValuesArray[i] = xValues.get(i);
                yValuesArray[i] = yValues.get(i);
            }
            spd.addData(xValuesArray, yValuesArray,
                    "" + amtDB.getSequenceForRun(runEntry));
            if (outPW != null)
                for (int i = 0; i < xValuesArray.length; i++)
                {
                    outPW.println(xValuesArray[i] + "\t" + yValuesArray[i] + "\t" + amtDB.getSequenceForRun(runEntry));
                }
        }
        if (outPW != null)
            outPW.close();
        spd.setAxisLabels("Time", "Predicted Hydrophobicity");
        if (showCharts)
            spd.setVisible(true);
        return spd.getChart();

    }


    /**
     *
     * @return
     */
    protected JFreeChart histogramIDProbs()
    {
        HistogramDataset dataset = new HistogramDataset();

        List<Double> histDataList = new ArrayList<Double>();


        AmtPeptideEntry[] entries = amtDB.getEntries();
        for (AmtPeptideEntry entry : entries)
        {
            if (entry.getNumObservations() > 1)
                histDataList.add((double) entry.calculateIDProbability());
        }
        double[] histData = new double[histDataList.size()];

        for (int i=0; i<histData.length; i++)
            histData[i] = histDataList.get(i);

        ApplicationContext.setMessage("Mean ID probability: " +
                BasicStatistics.mean(histData));

        dataset.addSeries("ID probablity",histData,200);


        JFreeChart chart = ChartFactory.createHistogram(
                "id_probability","id_probability","id_probability", dataset, PlotOrientation.VERTICAL,
                false,false,false);
        PanelWithChart pwc = new PanelWithChart(chart);
        pwc.setName("id probs");

        if (showCharts)
            pwc.displayInTab();
        return pwc.getChart();
    }



    /**
     * Show a histogram of the standard deviations of hydrophobicity observations per peptide
     * @return
     */
    protected JFreeChart histogramHydroStdDev()
    {
        HistogramDataset dataset = new HistogramDataset();

        List<Double> histDataList = new ArrayList<Double>();


            AmtPeptideEntry[] entries = amtDB.getEntries();
        for (int i=0; i<amtDB.numEntries(); i++)
        {
            if (entries[i].getNumObservations() > 1)
                histDataList.add(entries[i].getHydrophobicityStandardDeviation());
        }
        double[] histData = new double[histDataList.size()];

        for (int i=0; i<histData.length; i++)
            histData[i] = histDataList.get(i);

        ApplicationContext.setMessage("Mean H standard deviation: " +
                BasicStatistics.mean(histData));        

        dataset.addSeries("hyd_std_dev",histData,100);


        JFreeChart chart = ChartFactory.createHistogram(
                "hydstddev","hydstddev","hydstddev", dataset, PlotOrientation.VERTICAL,
                false,false,false);
        PanelWithChart pwc = new PanelWithChart(chart);
        pwc.setName("mean h std devs");

        if (showCharts)
            pwc.displayInTab();
        return pwc.getChart();
    }

    /**
     * Histogram all masses in the database
     * @return
     */
    protected JFreeChart histogramMasses()
    {
        HistogramDataset dataset = new HistogramDataset();

        double[] histData = new double[amtDB.numEntries()];


            AmtPeptideEntry[] entries = amtDB.getEntries();
        for (int i=0; i<amtDB.numEntries(); i++)
        {
           histData[i] = entries[i].getMass();


        }
        dataset.addSeries("masses",histData,500);

        JFreeChart chart = ChartFactory.createHistogram(
                "masses","masses","masses", dataset, PlotOrientation.VERTICAL,
                false,false,false);
        ChartDialog cd2 = new ChartDialog(chart);
        cd2.setSize(800,600);
        PanelWithChart pwc = new PanelWithChart(chart);

        pwc.setName("masses");

        if (showCharts)
            pwc.displayInTab();
        return pwc.getChart();
    }

    protected JFreeChart plotMassAndH()
    {
        AmtPeptideEntry[] entries = amtDB.getEntries();
        double[] entryMasses = new double[entries.length];
        double[] entryHs = new double[entries.length];

        for (int i=0; i<entries.length; i++)
        {
            entryMasses[i] = entries[i].getMass();
            entryHs[i] = entries[i].getMedianObservedHydrophobicity();
        }

        PanelWithScatterPlot spd = new PanelWithScatterPlot(entryHs, entryMasses, "AMT Database Entries");
        spd.setAxisLabels("H","Mass");
        spd.setPointSize(3);

        if (ms1FeatureSet != null)
        {
            AmtDatabaseMatcher matcher = new AmtDatabaseMatcher();
            AmtDatabaseFeatureSetGenerator featureGen =
                    new AmtDatabaseFeatureSetGenerator(amtDB);
            Feature[] dbFeatures =
                    featureGen.createFeaturesForModifications(
                            amtDB.getAminoacidModifications());
            Feature[] ms1Features = ms1FeatureSet.getFeatures();
            matcher.calculateFeatureHydrophobicities(ms1Features, null,
                    new FeatureSet(dbFeatures), 1);
            double[] ms1FeatureMasses = new double[ms1Features.length];
            double[] ms1FeatureHs = new double[ms1Features.length];
            for (int i=0; i<ms1Features.length; i++)
            {
                ms1FeatureMasses[i] = ms1Features[i].getMass();
                ms1FeatureHs[i] = AmtExtraInfoDef.getObservedHydrophobicity(ms1Features[i]);
            }
            spd.addData(ms1FeatureHs, ms1FeatureMasses, "MS1 Features");
        }

        spd.setVisible(true);

        return spd.getChart();
    }
}
