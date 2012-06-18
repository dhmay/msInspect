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

import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.ChartDialog;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.util.List;
import java.util.ArrayList;


/**
 * Command linemodule for deconvolution.  The "quant" and "icat" commands subclass
 * this module... all the work is done here.
 */
public class DeconvoluteCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(DeconvoluteCommandLineModule.class);

    protected File[] files;
    protected File outFile;
    protected File outDir;


    protected boolean deconvolute = true;
    protected boolean quant = false;
    protected int scanWindow = 6;
    protected double massWindow = .4;
    protected float lightTagWeight = -1.f;
    protected float heavyTagWeight = -1.f;
    protected char labeledResidue = ' ';
    protected int maxLabelCount = 3;
    protected MSRun run = null;
    protected int intensityType = FeatureSet.DEFAULT_INTENSITY_TYPE;
    protected float massTolerance = AnalyzeICAT.DEFAULT_DELTA_MASS;
    protected int massToleranceType = AnalyzeICAT.DEFAULT_DELTA_MASS_TYPE;
    protected float timeTolerance = AnalyzeICAT.DEFAULT_DELTA_TIME;
    protected boolean sumIntensities = true;

    protected boolean showCharts = false;

    AnalyzeICAT.IsotopicLabel label = AnalyzeICAT.icatLabel; // Default to ICAT parameters for quantitation



    protected static final String[] intensityTypeStrings =
            {
                    "total",
                    "max",
                    "recalculated"
            };

    public DeconvoluteCommandLineModule()
    {
        init();
    }

    protected CommandLineArgumentDefinition[] createCommonArgDefs()
    {
        CommandLineArgumentDefinition[] argDefs =
                {
                    createUnnamedSeriesFileArgumentDefinition(true, "Input File(s)"),
                    new FileToWriteArgumentDefinition("out",false, "Output File"),
                    new DirectoryToWriteArgumentDefinition("outdir",false, "Output Directory (for multiple files)"),
                    new DecimalArgumentDefinition("masswindow",false,"Mass Window", massWindow),
                    new IntegerArgumentDefinition("scanwindow",false,"Scan Window", scanWindow),
                    new DecimalArgumentDefinition("lighttagweight",false,"Light tag weight",
                            lightTagWeight),
                    new DecimalArgumentDefinition("heavytagweight",false,"Heavy tag weight",
                            heavyTagWeight),
                    new IntegerArgumentDefinition("maxlabelcount",false,"Maximum Label Count",
                            maxLabelCount),
                    new StringArgumentDefinition("labeledresidue",false,"Labeled Residue"),
                    new FileToReadArgumentDefinition("msfile",false,"mzXml File"),
                    new EnumeratedValuesArgumentDefinition("intensitytype", false, "Intensity type",intensityTypeStrings),
                    new DecimalArgumentDefinition("deltatime",false,"Time Tolerance", timeTolerance),
                    new DeltaMassArgumentDefinition("deltamass",false,"Mass Tolerance",
                            new DeltaMassArgumentDefinition.DeltaMassWithType(massTolerance, massToleranceType)),
                    new BooleanArgumentDefinition("showcharts", false,"Show Charts", showCharts),
                    new BooleanArgumentDefinition("sumintensities", false,
                            "If true, deconvoluted feature intensities reflect the sum of all component feature intensities.  If false, intensity of most-intense feature is kept.",
                            sumIntensities),

               };
        return argDefs;
    }

    protected void init()
    {
        mCommandName = "deconvolute";
        mShortDescription = "Deconvolute";
        mHelpMessage = "Deconvolute";

        addArgumentDefinitions(createCommonArgDefs());
        addArgumentDefinition(new BooleanArgumentDefinition("quant",false, "Quantitate", false));
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        files = this.getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        if (outFile == null)
        {
            assertArgumentPresent("outdir");
            outDir = getFileArgumentValue("outdir");
        }
        else
            assertArgumentAbsent("outdir");

        if (files.length > 1)
            assertArgumentPresent("outdir");

        if (hasArgumentValue("quant"))
            quant = getBooleanArgumentValue("quant");
        massWindow = getDoubleArgumentValue("masswindow");
        scanWindow = getIntegerArgumentValue("scanwindow");

        lightTagWeight = getFloatArgumentValue("lighttagweight");
        heavyTagWeight = getFloatArgumentValue("heavytagweight");

        maxLabelCount = getIntegerArgumentValue("maxlabelcount");

        String labeledResidueString = getStringArgumentValue("labeledresidue");
        if (hasArgumentValue("labeledresidue"))
        {
            if (labeledResidueString.length() > 1 ||
                !Character.isLetter(labeledResidueString.charAt(0)))
                throw new ArgumentValidationException("Bad residue value " + labeledResidueString);
            labeledResidue = labeledResidueString.charAt(0);
        }

        File runFile = getFileArgumentValue("msfile");
        if (runFile != null)
        {
            try
            {
                run = MSRun.load(runFile.getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException("Failed to load file " +
                        runFile.getAbsolutePath(),e);
            }
        }

        if (hasArgumentValue("intensitytype"))
        {
            int rawIntensityType = ((EnumeratedValuesArgumentDefinition)
                getArgumentDefinition("intensitytype")).getIndexForArgumentValue(getStringArgumentValue("intensitytype"));
            switch(rawIntensityType)
            {
                case 0:
                    intensityType = FeatureSet.TOTAL_INTENSITY;
                    break;
                case 1:
                    intensityType = FeatureSet.MAX_INTENSITY;
                    break;
                case 2:
                    intensityType = FeatureSet.RECALCULATED_INTENSITY;
                    if (run == null)
                        throw new ArgumentValidationException("When specifying intensityType=recalc, you must specify an msFile as well");
                    break;
            }
        }

        timeTolerance = getFloatArgumentValue("deltatime");

        DeltaMassArgumentDefinition.DeltaMassWithType deltaMassResult =
                getDeltaMassArgumentValue("deltamass");
        massTolerance = deltaMassResult.getDeltaMass();
        massToleranceType = deltaMassResult.getDeltaMassType();

        // ???? TODO: Expand these checks; deltaMass and deltaTime and maxLabelCount have no meaning if ! --quant
        if (!quant)
        {
            if ( lightTagWeight != -1.f || heavyTagWeight != -1.f )
                throw new ArgumentValidationException("Tag weights can only be specified along with --quant");
        }
        else if ((lightTagWeight > -1.f && heavyTagWeight == -1.f ) ||
                (heavyTagWeight > -1.f && lightTagWeight == -1.f ))
            throw new ArgumentValidationException("Must specify both light and heavy tag weights explicitly, or leave both to default");
        else if (lightTagWeight > -1.f && heavyTagWeight > -1.f) {
            // Build an Isotopic label from the given params. NB: we pass in the mass *delta*, not the heavy weight
            label = new AnalyzeICAT.IsotopicLabel(lightTagWeight, heavyTagWeight - lightTagWeight, labeledResidue, maxLabelCount);
        }

        showCharts = getBooleanArgumentValue("showcharts");
        sumIntensities = getBooleanArgumentValue("sumintensities");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        for (File file : files)
        {
            ApplicationContext.setMessage("Processing file " + file.getName());
            File outputFile = outFile;
            if (outFile == null)
            {
                outputFile =
                        new File(outDir, file.getName().substring(0, file.getName().indexOf(".")) + ".quant.tsv");
            }
            processFile(file, outputFile);
        }                          
    }

    /**
     * do the actual work
     */
    public void processFile(File inputFile, File outputFile)
            throws CommandLineModuleExecutionException
    {
        try
        {

            PrintWriter pw = null;
            if (null != outputFile)
            {
                try
                {
                    pw = new PrintWriter(new FileOutputStream(outputFile));
                }
                catch (java.io.FileNotFoundException e)
                {
                    throw new CommandLineModuleExecutionException("Error creating PrintWriter from file " +
                            outputFile.getAbsolutePath() + ", file not found");
                }
            }
            else
                pw = new PrintWriter(System.out);


            try
            {
                FeatureSet fs = new FeatureSet(inputFile);
                if (deconvolute)
                    fs = fs.deconvolute(scanWindow, massWindow, sumIntensities);
                if (quant)
                    fs = fs.quant(label, intensityType, run, massTolerance, massToleranceType,
                            timeTolerance);
                fs.save(pw);
//PrintWriter deconPW = new PrintWriter(new File("/home/dhmay/temp/decon.tsv"));
//FeatureSet fsDecon = (FeatureSet)fs.clone();
//List<Feature> deconOnlyFeatures = new ArrayList<Feature>();
//for (Feature feature : fs.getFeatures())
//  if (feature.comprised != null && feature.comprised.length > 1) {for (Spectrum.Peak comprisedFeature : feature.comprised) deconOnlyFeatures.add((Feature) comprisedFeature);}
//System.err.println("Decon only: " + deconOnlyFeatures.size());
//fsDecon.setFeatures(deconOnlyFeatures.toArray(new Feature[deconOnlyFeatures.size()]));
//fsDecon.save(deconPW);
                if (quant && showCharts)
                {
                    List<Float> ratioList = new ArrayList<Float>();
                    for (Feature feature : fs.getFeatures())
                    {
                        double ratio = IsotopicLabelExtraInfoDef.getRatio(feature);
                        if (ratio != IsotopicLabelExtraInfoDef.NO_RATIO_FOR_FEATURE)
                            ratioList.add((float) ratio);
                    }

                    PanelWithHistogram pwh =
                            new PanelWithHistogram(ratioList, "Quantitation Ratios");
                    ChartDialog histDialog = new ChartDialog(pwh);
                    histDialog.setVisible(true);
                }
            }
            finally
            {
                if (null != outputFile && null != pw)
                    pw.close();
            }
        }
        catch (Exception x)
        {
            throw new CommandLineModuleExecutionException(x);
        }
    }

}
