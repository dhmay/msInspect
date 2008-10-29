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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.gui.chart.*;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.quant.QuantitationVisualizer;
import org.fhcrc.cpl.viewer.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.viewer.gui.util.PanelWithSpectrumChart;
import org.apache.log4j.Logger;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;


/**
 * test
 */
public class PeptideQuantVisualizationCLM extends BaseCommandLineModuleImpl
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

    protected Map<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>> proteinPeptideFractionChargeFilesMap =
            new HashMap<String, Map<String, Map<String, Map<Integer, List<Pair<File, File>>>>>>();

    protected static final String DUMMY_PROTEIN_NAME = "DUMMY_PROTEIN";

    protected QuantitationVisualizer quantVisualizer;



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
                        createIntegerArgumentDefinition("spectrumimageheight", false, "Spectrum image height",
                                800),
                        createIntegerArgumentDefinition("resolution", false, "resolution (number of breaks per Thompson",
                                PanelWithSpectrumChart.DEFAULT_RESOLUTION),
                        createIntegerArgumentDefinition("width", false, "Image width (used for both)", 1200),
                        createIntegerArgumentDefinition("scanimageheight", false,
                                "Height of EACH per-scan image, in the output file",
                                100),
                        createDirectoryToReadArgumentDefinition("outdir", true, "Output directory"),
                        createIntegerArgumentDefinition("maxscansimageheight", false,
                                "Maximum overall height for the all-scans line plot image (overrides scansfileimageheight)",
                                4000),
                        createFileToReadArgumentDefinition("pepxml", true, "pepXML file"),
                        createFileToReadArgumentDefinition("mzxml", false, "mzXML file"),

                        createDirectoryToReadArgumentDefinition("mzxmldir", false, "mzXML directory"),
                        createStringArgumentDefinition("peptides", false, "comma-separated list of peptides to examine"),
                        createStringArgumentDefinition("proteins", false, "comma-separated list of proteins to examine"),
                        createIntegerArgumentDefinition("scan", false, "Scan number of desired quantitation event",0),
                        createStringArgumentDefinition("fractions", false, "Fraction containing desired quantitation event"),

                        createDecimalArgumentDefinition("minpprophet", false, "minimum PeptideProphet value",
                                0),
                        createIntegerArgumentDefinition("paddingscans", false,
                                "number of scans before and after quant envelope to display", 5),
                        createDecimalArgumentDefinition("mzpadding", false,
                                "amount of m/z space to display around quant", 1.5),
                        createIntegerArgumentDefinition("numpeaksaboveheavy", false,
                                "number of peaks above the heavy-ion monoisotope to display", 4),
                        createBooleanArgumentDefinition("show3dplots", false,
                                "Show 3D plot? (takes more time)", true),
                        createIntegerArgumentDefinition("3drotation", false,
                                "Rotation angle for 3D plot", PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_ROTATION),
                        createIntegerArgumentDefinition("3dtilt", false,
                                "Tilt angle for 3D plot", PanelWithSpectrumChart.DEFAULT_CONTOUR_PLOT_TILT),
                        createIntegerArgumentDefinition("3dwidth", false,
                                "Image width for 3D plot", 1000),
                        createIntegerArgumentDefinition("3dheight", false,
                                "Image height for 3D plot", 1000),
                        createBooleanArgumentDefinition("3dshowaxes", false,
                                "Include axes on 3D plot?", true),
                        createBooleanArgumentDefinition("infooncharts", false,
                                "Write quantitation information directly on the charts?", false),
                        createFileToWriteArgumentDefinition("outtsv", false,
                                "Output TSV file"),
                        createFileToWriteArgumentDefinition("outhtml", false,
                                "Output HTML file"),
                };
        addArgumentDefinitions(argDefs);
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

        quantVisualizer.setMinPeptideProphet(getFloatArgumentValue("minpprophet"));




        int numModeArgs = 0;
        if (hasArgumentValue("peptides")) numModeArgs++;
        if (hasArgumentValue("proteins")) numModeArgs++;
        if (hasArgumentValue("scan")) numModeArgs++;

        if (numModeArgs != 1)
            throw new ArgumentValidationException(
                    "Please supply either peptides, or proteins, or scan (but no more than one");

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
        quantVisualizer.setWriteHTMLAndText(!hasArgumentValue("scan"));

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
        quantVisualizer.visualizeQuantEvents();
    }


}
