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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.quant.gui.QuantitationVisualizer;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;


/**
 * test
 */
public class PeptideQuantVisualizationCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PeptideQuantVisualizationCLM.class);


    protected File pepXmlFile;

    //The amount by which two features can differ and still be combined into the same single representation.
    //TODO: convert this to log space, parameterize
    protected float maxCombineFeatureRatioDiff = 0.2f;


    protected boolean writeHTMLAndText = true;
    protected Set<String> peptidesFound = new HashSet<String>();
    protected int sidebarWidth = 180;

    protected static final String DUMMY_PROTEIN_NAME = "DUMMY_PROTEIN";

    protected QuantitationVisualizer quantVisualizer;

    protected File outTsvFile;


    public PeptideQuantVisualizationCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "quantcharts";

        mHelpMessage ="Create charts showing the MS1 spectra near quantitative events for specified peptides or " +
                "proteins.  Includes an HTML index to all of the charts, broken down by protein, peptide, " +
                "fraction, charge, etc.";
        mShortDescription = "Create charts showing the MS1 spectra near quantitative events for specified peptides " +
                "or proteins";

        CommandLineArgumentDefinition[] argDefs =
                {
                        new DirectoryToWriteArgumentDefinition("outdir", true, "Output directory"),

                        new FileToReadArgumentDefinition("pepxml", true, "pepXML file"),
                        new FileToReadArgumentDefinition("mzxml", false, "mzXML file"),

                        new DirectoryToReadArgumentDefinition("mzxmldir", false, "mzXML directory"),
                        new StringArgumentDefinition("peptides", false, "comma-separated list of peptides to examine"),
                        new StringArgumentDefinition("proteins", false, "comma-separated list of proteins to examine"),
                        new IntegerArgumentDefinition("scan", false,
                                "Scan number of desired quantitation event (leave at 0 for all scans)",0),
                        new StringArgumentDefinition("fractions", false,
                                "Fraction containing desired quantitation event"),
                        new DecimalArgumentDefinition("minpprophet", false, "minimum PeptideProphet value",
                                0),
                        new BooleanArgumentDefinition("show3dplots", false,
                                "Show 3D plot? (takes more time)", true),

                        new FileToWriteArgumentDefinition("outtsv", false,
                                "Output TSV file"),
                        new FileToWriteArgumentDefinition("outhtml", false,
                                "Output HTML file"),
                };
        CommandLineArgumentDefinition[] advancedArgDefs =
                {
                        new IntegerArgumentDefinition("paddingscans", false,
                                "number of scans before and after quant envelope to display", 5),
                        new DecimalArgumentDefinition("mzpadding", false,
                                "amount of m/z space to display around quant", 1.5),
                        new IntegerArgumentDefinition("numpeaksaboveheavy", false,
                                "number of peaks above the heavy-ion monoisotope to display", 4),
                        new IntegerArgumentDefinition("maxscansimageheight", false,
                                "Maximum overall height for the all-scans line plot image " +
                                        "(overrides scansfileimageheight)",
                                QuantitationVisualizer.DEFAULT_MAX_SINGLE_SCANS_TOTAL_IMAGE_HEIGHT),
                        new IntegerArgumentDefinition("spectrumimageheight", false,
                                "Image height  (used for spectrum, scans, and sum scan intensities charts)",
                                QuantitationVisualizer.DEFAULT_SPECTRUM_IMAGE_HEIGHT),
                        new IntegerArgumentDefinition("resolution", false, "resolution (number of breaks per Thompson",
                                PanelWithSpectrumChart.DEFAULT_RESOLUTION),
                        new IntegerArgumentDefinition("width", false,
                                "Image width (used for spectrum, scans, and sum scan intensities charts)",
                                QuantitationVisualizer.DEFAULT_SPECTRUM_IMAGE_WIDTH),
                        new IntegerArgumentDefinition("scanimageheight", false,
                                "Height of EACH per-scan image, in the output file",
                                QuantitationVisualizer.DEFAULT_SINGLE_SCAN_IMAGE_HEIGHT),
                        new IntegerArgumentDefinition("3drotation", false,
                                "Rotation angle for 3D plot", PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_ROTATION),
                        new IntegerArgumentDefinition("3dtilt", false,
                                "Tilt angle for 3D plot", PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_TILT),
                        new IntegerArgumentDefinition("3dwidth", false,
                                "Image width for 3D plot", QuantitationVisualizer.DEFAULT_IMAGE_WIDTH_3D),
                        new IntegerArgumentDefinition("3dheight", false,
                                "Image height for 3D plot", QuantitationVisualizer.DEFAULT_IMAGE_HEIGHT_3D),
                        new BooleanArgumentDefinition("3dshowaxes", false,
                                "Include axes on 3D plot?", true),
                        new BooleanArgumentDefinition("infooncharts", false,
                                "Write quantitation information directly on the charts?", false),
                        new DecimalArgumentDefinition("peakdistance", false,
                                "Distance, in Daltons, between peaks.  This is configurable in Q3, " +
                                "so it has to be configurable here.  Used in generating the intensity sum chart",
                                PanelWithSpectrumChart.DEFAULT_PEAK_SEPARATION_MASS),
                        new DecimalArgumentDefinition("peakmasstoleranceppm", false,
                                "Mass tolerance, in PPM, around each theoretical peak to consider part of " +
                                        "the peptide being quantitated.  Used in generating the intensity sum chart",
                                PanelWithSpectrumChart.DEFAULT_PEAK_TOLERANCE_PPM),
                        new BooleanArgumentDefinition("nocharts", false,
                                "Don't create charts, just write out HTML and/or TSV files", false),
                        new BooleanArgumentDefinition("markallbad", false,
                                "Mark all events as Bad, rather than Unknown", false),
                };
        addArgumentDefinitions(argDefs);
        this.addArgumentDefinitions(advancedArgDefs, true);

    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        pepXmlFile = getFileArgumentValue("pepxml");

        if (!hasArgumentValue("mzxml") && !hasArgumentValue("mzxmldir"))
            throw new ArgumentValidationException("Must supply mzxml or mzxmldir argument");


        quantVisualizer = new QuantitationVisualizer();

        quantVisualizer.setResolution(getIntegerArgumentValue("resolution"));
        quantVisualizer.setSpectrumImageHeight(getIntegerArgumentValue("spectrumimageheight"));
        quantVisualizer.setImageWidth(getIntegerArgumentValue("width"));
        quantVisualizer.setScanImageHeight(getIntegerArgumentValue("scanimageheight"));
        quantVisualizer.setMaxScansImageHeight(getIntegerArgumentValue("maxscansimageheight"));
        quantVisualizer.setMzXmlFile(getFileArgumentValue("mzxml"));
        quantVisualizer.setMzXmlDir(getFileArgumentValue("mzxmldir"));
        quantVisualizer.setOutDir(getFileArgumentValue("outdir"));
        if (hasArgumentValue("outtsv"))
            quantVisualizer.setOutTsvFile(getFileArgumentValue("outtsv"));
        if (hasArgumentValue("outhtml"))
            quantVisualizer.setOutHtmlFile(getFileArgumentValue("outhtml"));
        if (hasArgumentValue("outtsv") || hasArgumentValue("outhtml"))
        {
            quantVisualizer.setWriteHTMLAndText(true);
            if (!hasArgumentValue("outtsv"))
                quantVisualizer.setOutTsvFile(TempFileManager.createTempFile("qurate_dummy_out.tsv", this));

            if (!hasArgumentValue("outhtml"))
                quantVisualizer.setOutHtmlFile(TempFileManager.createTempFile("qurate_dummy_out.html", this));

        }

        quantVisualizer.setMinPeptideProphet(getFloatArgumentValue("minpprophet"));
        quantVisualizer.setPeakSeparationMass(getFloatArgumentValue("peakdistance"));
        quantVisualizer.setPeakTolerancePPM(getFloatArgumentValue("peakmasstoleranceppm"));

        quantVisualizer.setShouldCreateCharts(!getBooleanArgumentValue("nocharts"));
        quantVisualizer.setMarkAllEventsBad(getBooleanArgumentValue("markallbad"));

        if (quantVisualizer.isMarkAllEventsBad())
            ApplicationContext.infoMessage("Note: all events will be marked as 'Bad'");
        if (!quantVisualizer.isShouldCreateCharts())
            ApplicationContext.infoMessage("Note: no charts will be created");


        int numModeArgs = 0;
        if (hasArgumentValue("peptides")) numModeArgs++;
        if (hasArgumentValue("proteins")) numModeArgs++;
        if (hasArgumentValue("scan")) numModeArgs++;

        if (numModeArgs > 1)
            throw new ArgumentValidationException(
                    "Please supply either peptides, or proteins, or scan (but no more than one");
        if (numModeArgs == 0)
            ApplicationContext.infoMessage("WARNING! Generating charts for ALL events in input file(s).");

        if (hasArgumentValue("peptides"))
        {
            String peptidesString = getStringArgumentValue("peptides");
            String[] peptidesArray = peptidesString.split(",");
            Set<String> peptidesToExamine = new HashSet<String>();
            peptidesToExamine.addAll(Arrays.asList(peptidesArray));
            quantVisualizer.setPeptidesToExamine(peptidesToExamine);
        }

        if (hasArgumentValue("proteins"))
        {
            String proteinsString = getStringArgumentValue("proteins");
            String[] proteinsArray = proteinsString.split(",");
            Set<String> proteinsToExamine = new HashSet<String>();
            proteinsToExamine.addAll(Arrays.asList(proteinsArray));
            quantVisualizer.setProteinsToExamine(proteinsToExamine);
        }

        quantVisualizer.setScan(getIntegerArgumentValue("scan"));

        if (hasArgumentValue("fractions"))
        {
            String fractionsString = getStringArgumentValue("fractions");
            String[] fractionsArray = fractionsString.split(",");
            Set<String> fractionsToExamine = new HashSet<String>();
            fractionsToExamine.addAll(Arrays.asList(fractionsArray));
            quantVisualizer.setFractionsToExamine(fractionsToExamine);
//for (String fraction : fractionsToExamine) System.err.println(fraction);            
        }

        quantVisualizer.setNumPaddingScans(getIntegerArgumentValue("paddingscans"));
        quantVisualizer.setMzPadding(getFloatArgumentValue("mzpadding"));
        quantVisualizer.setNumHeavyPeaksToPlot(getIntegerArgumentValue("numpeaksaboveheavy"));
        quantVisualizer.setShow3DPlots(getBooleanArgumentValue("show3dplots"));
        quantVisualizer.setShow3DAxes(getBooleanArgumentValue("3dshowaxes"));
        quantVisualizer.setRotationAngle3D(getIntegerArgumentValue("3drotation"));
        quantVisualizer.setTiltAngle3D(getIntegerArgumentValue("3dtilt"));
        quantVisualizer.setImageWidth3D(getIntegerArgumentValue("3dwidth"));
        quantVisualizer.setImageHeight3D(getIntegerArgumentValue("3dheight"));
        quantVisualizer.setWriteInfoOnCharts(getBooleanArgumentValue("infooncharts"));

    }



    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Iterator<FeatureSet> fsi;
        try
        {
            fsi = new PepXMLFeatureFileHandler.PepXMLFeatureSetIterator(pepXmlFile);
        }
        catch (IOException e)
        {
            try
            {
                if (quantVisualizer.getFractionsToExamine() != null)
                    throw new CommandLineModuleExecutionException("'fractions' argument provided on a non-pepxml file.  Quitting");
                FeatureSet featureSet = new FeatureSet(pepXmlFile);
                List<FeatureSet> dummyList = new ArrayList<FeatureSet>();
                dummyList.add(featureSet);
                fsi = dummyList.iterator();
            }
            catch (Exception e2)
            {
                throw new CommandLineModuleExecutionException("Failed to open file ",e2);
            }
        }
        quantVisualizer.setFeatureSetIterator(fsi);
        try
        {
            quantVisualizer.visualizeQuantEvents();
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        outTsvFile = quantVisualizer.getOutTsvFile();
    }

    public File getOutTsvFile()
    {
        return outTsvFile;
    }
}
