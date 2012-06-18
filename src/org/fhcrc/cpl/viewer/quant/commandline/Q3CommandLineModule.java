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
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.viewer.quant.Q3;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import javax.xml.stream.XMLStreamException;

/**
 * Command linemodule for Q3 quantitation
 */
public class Q3CommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(Q3CommandLineModule.class);

    protected boolean forceOutput = false;
    protected boolean mimicXpress = false;
    protected boolean noSentinels = false;
    protected boolean compatMode = false;
    protected boolean debugMode = false;
    protected boolean stripExistingQ3 = false; // If there are existing analysis_results, remove them


    String alternateMzXmlDir = null;

    protected char labeledResidue = '\0';
    protected float massDiff = 0f;

    protected float minPeptideProphet = -1f;

    protected float maxFracDeltaMass = -1f;
    protected boolean maxFracDeltaMassIsPPM = true;

    protected FeatureSet inFeatureSet = null;
    protected File[] inFiles = null;
    protected File outFile = null;
    protected File outDir = null;


    public Q3CommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "q3new";

        mHelpMessage = "Run Q3 quantitation on a pepXml file. Labeled residue and mass difference must be specified.";

        mShortDescription = "Run Q3 quantitation on a pepXml file";

        CommandLineArgumentDefinition[] argDefs =
                {
                        //mandatory
                        this.createUnnamedSeriesFileArgumentDefinition(true, "Input file(s)"),

                        //optional
                        new StringArgumentDefinition("labeledresidue", false, "Labeled residue"),
                        new DecimalArgumentDefinition("massdiff", false, "Mass difference between heavy and light forms of the labeled residue"),
                        new FileToWriteArgumentDefinition("out", false, "Output file"),
                        new BooleanArgumentDefinition("mimicxpress", false, "Mimic XPress? (default false)"),
                        new BooleanArgumentDefinition("forceoutput", false, "Force output (default false)"),
                        new BooleanArgumentDefinition("nosentinels", false, "No sentinels (default false)"),
                        new BooleanArgumentDefinition("compat", false, "Match behavior of the original R code when the center scan has fewer than three matching isotopes (default false)"),
                        new BooleanArgumentDefinition("debug", false, "Output extra debugging information in the pepXML"),
                        new BooleanArgumentDefinition("stripoldq3", false,
                                "Strip existing analysis_results and Q3 analysis_summary elements from the file", stripExistingQ3),

                        new DecimalArgumentDefinition("minpeptideprophet", false, "Minimum PeptideProphet score"),
                        new DeltaMassArgumentDefinition("maxfracdeltamass", false, "Maximum fractional delta mass"),
                        new StringArgumentDefinition("alternatemzxmldir", false, "Alternate mzXML directory"),
                        // new DeltaMassArgumentDefinition("masstolerance", false, "mass tolerance"),

                        // optional short-form arguments
                        new StringArgumentDefinition("n", false, "Label definition (e.g. -nC,3.0100645"),
                        new StringArgumentDefinition("d", false, "Alternate mzXML directory"),
                        new DeltaMassArgumentDefinition("m", false, "Mass tolerance"),
                        new DirectoryToWriteArgumentDefinition("outdir", false,
                            "Output Directory (for handling multiple files)"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFiles = this.getUnnamedSeriesFileArgumentValues();


        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");
        if (inFiles.length > 1)
        {
            assertArgumentAbsent("out");
        }
        else
        {
            if (hasArgumentValue("out"))
                assertArgumentAbsent("outdir");
            else
                assertArgumentPresent("outdir");
        }

        if (hasArgumentValue("labeledresidue"))
            labeledResidue = getStringArgumentValue("labeledresidue").charAt(0);
        if (hasArgumentValue("massdiff"))
            massDiff = getFloatArgumentValue("massdiff");
        if (hasArgumentValue("forceoutput"))
            forceOutput = getBooleanArgumentValue("forceoutput");
        if (hasArgumentValue("mimicxpress"))
            mimicXpress = getBooleanArgumentValue("mimicxpress");
        if (hasArgumentValue("noesentinels"))
            noSentinels = getBooleanArgumentValue("nosentinels");
        if (hasArgumentValue("compat"))
            compatMode = getBooleanArgumentValue("compat");
        if (hasArgumentValue("debug"))
            debugMode = getBooleanArgumentValue("debug");
        stripExistingQ3 = getBooleanArgumentValue("stripoldq3");
        if (hasArgumentValue("minpeptideprophet"))
        {
            minPeptideProphet = getFloatArgumentValue("minpeptideprophet");
        }
        if (hasArgumentValue("maxfracdeltamass"))
        {
            DeltaMassArgumentDefinition.DeltaMassWithType d =
                    getDeltaMassArgumentValue("maxfracdeltamass");
            maxFracDeltaMass = d.getDeltaMass();
            maxFracDeltaMassIsPPM = d.getDeltaMassType() == AnalyzeICAT.DELTA_MASS_PPM;
        }
        alternateMzXmlDir = getStringArgumentValue("alternatemzxmldir");


        if (hasArgumentValue("n"))
        {
            parseShortLabel(getStringArgumentValue("n"));
        }
        if (hasArgumentValue("d"))
        {
            alternateMzXmlDir = getStringArgumentValue("d");
        }
        if (hasArgumentValue("m"))
        {
            // mass tolerance currently not used in Q3
        }
    }


    /**
     * Parse an XPRESS-style label definiton -nA,###
     */
    private void parseShortLabel(String label)
    {
        String[] val = label.split(",");
        if (val.length < 2 || val[0].length() != 1 || val[1].length() < 1)
            throw new RuntimeException("Could not parse residue and mass from argument " + label);
        labeledResidue = val[0].charAt(0);
        massDiff = Float.parseFloat(val[1]);
    }

    /*
     * Throw a commandline exception on bad arguments
     */
    private static void quit(String msg) throws CommandLineModuleExecutionException
    {
        throw new CommandLineModuleExecutionException(msg);
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (inFiles.length == 1)
        {
            File outputFile = outFile;
            if (outputFile == null && outDir != null)
                outputFile = new File(outDir, inFiles[0].getName());
            processFile(inFiles[0], outputFile);
        }
        else
        {
            for (File inFile : inFiles)
            {
                File outputFile = null;
                if (outDir != null)
                    outputFile = new File(outDir, inFile.getName());
                processFile(inFile, outputFile);
            }
        }
    }

    public void processFile(File inFile, File outputFile) throws CommandLineModuleExecutionException
    {
        String inFilename = inFile != null ? inFile.getPath() : null;
        String outFilename = outputFile != null ? outputFile.getPath() : null;
        
        if (null == inFilename)
            quit("Please specify a pepXML file as input");

        if (null == outFilename)
        {
            outFilename = inFilename;
            forceOutput = true;
        }

        if ('\0' == labeledResidue)
            quit("Please specifiy the amino acid that is labeled");

        if (0 == massDiff)
            quit("Please specifiy the mass difference between the heavy and light labels");

        Q3 q3 = new Q3(inFilename, labeledResidue, massDiff, outFilename);

        if (null != alternateMzXmlDir)
            q3.setAlternateMzXmlDir(alternateMzXmlDir);
        if (minPeptideProphet >= 0)
            q3.setMinPeptideProphet(minPeptideProphet);
        if (maxFracDeltaMass >= 0)
            q3.setMaxFracDeltaMass(maxFracDeltaMass, maxFracDeltaMassIsPPM);
//        if (massTolerance >= 0)
//            q3.setMassTolerance(massTolerance, massToleranceIsPPM);

        q3.setCompatMode(compatMode);
        q3.setDebugMode(debugMode);
        q3.setForceOutput(forceOutput);
        q3.setMimicXpress(mimicXpress);
        q3.setNoSentinels(noSentinels);
        q3.setStripExistingQ3(stripExistingQ3);
        try
        {
            q3.apply();
        }
        catch(XMLStreamException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        catch(IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
