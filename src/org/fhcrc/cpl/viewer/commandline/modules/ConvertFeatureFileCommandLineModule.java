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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeaturePepXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.APMLFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.BasePepXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.HardklorFeatureFileHandler;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.PepXMLFeatureFileHandler;
import org.fhcrc.cpl.viewer.feature.extraction.SpectrumResampler;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 */
public class ConvertFeatureFileCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FilterFeaturesCommandLineModule.class);

    protected File[] inFeatureFiles = null;
    protected File outFeatureFile = null;
    protected File outDir = null;
    protected MSRun run = null;
    protected int inFileFormat = FILE_FORMAT_MSINSPECT;
    protected int outFileFormat = 0;
    protected boolean dumpWindow = false;
    protected File fastaFile = null;

    protected boolean forcePeptideProphet = false;
    protected double forcePeptideProphetValue = 0;

    protected String pepXmlSearchEngine = BasePepXmlWriter.DEFAULT_SEARCH_ENGINE;

    

    protected static final int FILE_FORMAT_MSINSPECT=0;
    protected static final int FILE_FORMAT_PEPXML=1;
    protected static final int FILE_FORMAT_SPECARRAY=2;
    protected static final int FILE_FORMAT_APML=3;
    protected static final int FILE_FORMAT_HARDKLOR=4;
    protected static final int FILE_FORMAT_MULTI_MSINSPECT=5;



    protected static final String[] formatStrings = {
            "msinspect",
            "pepxml",
            "specarraytsv",
            "apml",
            "hardklor",
            "multimsinspect"
        };

    protected static final String[] formatDescriptions = {
            "msInspect tab-separated values",
            "PepXML",
            "SpecArray tab-separated values",
            "APML v2.0 (XML)",
            "Hardklor text format",
            "Multiple msInspect feature sets in a single file, with a column indicating run"
        };

    protected static final int SPECARRAY_VERSION_1_0 = 0;
    protected static final int SPECARRAY_VERSION_1_2_0 = 1;

    protected int specArrayVersion = SPECARRAY_VERSION_1_2_0;

    protected final static String[] specArrayVersionValues = {"1.0","1.2"};

    public ConvertFeatureFileCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "convertfeaturefile";
        mUsageMessage = CommandLineModule.MODULE_USAGE_AUTOMATIC;
        StringBuffer helpMessageBuf =
                new StringBuffer("This command converts between different formats of feature files.  " +
                                 "Allowed formats:\n");
        for (int i=0; i<formatStrings.length; i++)
            helpMessageBuf.append("\t" + formatStrings[i] + ": " + formatDescriptions[i] + "\n");
        mHelpMessage = helpMessageBuf.toString();
        mShortDescription = "Convert between formats of feature files";
        

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true,
                                "Input feature file(s)"),
                        new FileToWriteArgumentDefinition("out",false,"output file"),
                        new DirectoryToWriteArgumentDefinition("outdir",false,"output directory (for multiple inputs)"),

                        new EnumeratedValuesArgumentDefinition("informat", false,
                                "input file format.  To force loading of an ambiguous file as a specific file type",
                                formatStrings, "msinspect"),
                        new EnumeratedValuesArgumentDefinition("outformat", true, "output file format",formatStrings),
                        new FileToReadArgumentDefinition("mzxml",false,null),
                        new BooleanArgumentDefinition("dumpwindow",false, null),
                        new EnumeratedValuesArgumentDefinition("specarrayversion",false,
                                "specArray version, for specArray conversions (default 1.2)",
                                specArrayVersionValues),
                        new DecimalArgumentDefinition("forcepeptideprophetvalue",false,
                                "Set the PeptideProphet probability of all features to this (for pepxml output)"),
                        new FileToReadArgumentDefinition("fasta", false,
                                "FASTA filepath to include in pepXML file (for outformat=pepxml only)"),
                        new StringArgumentDefinition("searchengine", false,
                                "Search engine to store in pepXML file (for outformat=pepxml ony)", pepXmlSearchEngine)
                };
        addArgumentDefinitions(argDefs);
    }



    protected int translateFormatString(String formatString)
    {
        int result = -1;
        for (int i = 0; i < formatStrings.length; i++)
        {
            String modeString = formatStrings[i];
            if (modeString.equalsIgnoreCase(formatString))
            {
                result = i;
                break;
            }
        }
        return result;
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        if (hasArgumentValue("informat"))
            inFileFormat = translateFormatString(getStringArgumentValue("informat"));
        if (inFileFormat == FILE_FORMAT_MULTI_MSINSPECT)
            throw new ArgumentValidationException("multimsinspect is not a valid input format, only output.  " +
                    "For now, anyway.");
        outFileFormat = translateFormatString(getStringArgumentValue("outformat"));
        inFeatureFiles = this.getUnnamedSeriesFileArgumentValues();
        outFeatureFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        pepXmlSearchEngine = getStringArgumentValue("searchengine");
        if (hasArgumentValue("searchengine") && outFileFormat != FILE_FORMAT_PEPXML)
               throw new ArgumentValidationException("Argument searchengine is only for pepXML output mode");

        if (hasArgumentValue("outdir"))
            assertArgumentAbsent("out");
        else
            assertArgumentPresent("out");
        if (inFeatureFiles.length > 1)
            assertArgumentPresent("outdir");

        if (hasArgumentValue("mzxml"))
        {
            try
            {
                run = MSRun.load(getFileArgumentValue("mzxml").getAbsolutePath());
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }

        if (hasArgumentValue("dumpwindow"))
        {
            dumpWindow = getBooleanArgumentValue("dumpwindow");
            if (run == null)
                throw new ArgumentValidationException("if dumpwindow is specified, must provide an mzxml file");
        }


        if (hasArgumentValue("specarrayversion"))
            specArrayVersion = ((EnumeratedValuesArgumentDefinition)
                               getArgumentDefinition("specarrayversion")).getIndexForArgumentValue(
                                        getStringArgumentValue("specarrayversion"));

        if (inFileFormat == FILE_FORMAT_SPECARRAY)
            assertArgumentPresent("mzxml");

        if (hasArgumentValue("forcepeptideprophetvalue"))
        {
            if (!(outFileFormat == FILE_FORMAT_PEPXML))
                throw new ArgumentValidationException("Argument forcepeptideprophetvalue is only for pepXML output mode");
            forcePeptideProphet = true;
            forcePeptideProphetValue = getDoubleArgumentValue("forcepeptideprophetvalue");
        }

        if (outFileFormat != FILE_FORMAT_PEPXML)
            assertArgumentAbsent("fasta");

        fastaFile = getFileArgumentValue("fasta");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        for (File file : inFeatureFiles)
        {
            File outputFile = outFeatureFile;
            if (outputFile == null)
            {
                String outputSuffix = "";
                switch (outFileFormat)
                {
                    case FILE_FORMAT_MSINSPECT:
                    case FILE_FORMAT_MULTI_MSINSPECT:
                    case FILE_FORMAT_SPECARRAY:
                    case FILE_FORMAT_HARDKLOR:
                        outputSuffix = "tsv";
                        break;
                    case FILE_FORMAT_PEPXML:
                        outputSuffix = "pep.xml";
                        break;
                    case FILE_FORMAT_APML:
                        outputSuffix = "apml.xml";
                        break;
                }
                outputFile = CommandLineModuleUtilities.createOutputFile(file, outputSuffix, outDir);
            }
            handleFile(file, outputFile);
        }
    }

    protected void handleFile(File inputFile, File outputFile)
        throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Loading features from " + inputFile.getAbsolutePath() + "...");
        FeatureSet featureSet = null;

        //pepxml files can have multiple FeatureSets in them
        List<FeatureSet> featureSets = null;
        switch (inFileFormat)
        {
            case FILE_FORMAT_SPECARRAY:
                featureSet = loadFeatureSetFromSpecArrayTSV(inputFile, specArrayVersion);
                break;
            case FILE_FORMAT_PEPXML:
                try
                {
                    featureSets =
                            PepXMLFeatureFileHandler.getSingletonInstance().loadAllFeatureSets(inputFile);
                    //The only valid format for the output if there are multiple sets is multi_tsv.
                    //If using a different format, ignore featureSets and use the first one, featureSet
                    featureSet = featureSets.get(0);
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException(
                            "Failed to load feature sets from PepXML file",e);
                }
                break;
            case FILE_FORMAT_APML:            
            case FILE_FORMAT_MSINSPECT:
            case FILE_FORMAT_HARDKLOR:
                try
                {
                    featureSet = new FeatureSet(inputFile);
                }
                catch (Exception e)
                {
                    throw new CommandLineModuleExecutionException("Problems opening feature file");
                }
                break;                
            default:
                throw new CommandLineModuleExecutionException(
                        "Don't know how to support the specified input file format yet");
        }
        PrintWriter pw = null;
        try
        {

            switch (outFileFormat)
            {
                case FILE_FORMAT_MSINSPECT:
                    if (dumpWindow)
                        createIntensityWindows(featureSet.getFeatures(), run);
                    pw = new PrintWriter(outputFile);
                    featureSet.save(pw,dumpWindow);
                    break;
                case FILE_FORMAT_MULTI_MSINSPECT:
                    if (featureSets == null || featureSets.size() == 1)
                    {
                        if (dumpWindow)
                            createIntensityWindows(featureSet.getFeatures(), run);
                        pw = new PrintWriter(outputFile);
                        featureSet.save(pw,dumpWindow);
                    }
                    else
                    {
                        pw = new PrintWriter(outputFile);

                        //check if all featuresets have basenames.  If so, use those for
                        //run identifier.  Otherwise, use a number
                        boolean allHaveBaseNames = true;
                        for (FeatureSet fSet : featureSets)
                        {
                            String baseName = MS2ExtraInfoDef.getFeatureSetBaseName(fSet);
                            if (baseName == null || baseName.length()<1)
                                allHaveBaseNames = false;
                        }

                        for (int i=0; i<featureSets.size(); i++)
                        {
                            ApplicationContext.setMessage("Writing FeatureSet " + (i+1) + "...");
                            FeatureSet fSet = featureSets.get(i);
                            File tempFile = TempFileManager.createTempFile("fset" + i + ".tsv", this);
                            fSet.save(tempFile);

                            FileReader fr = new FileReader(tempFile);
                            BufferedReader br = new BufferedReader(fr);
                            String line = null;
                            if (i==0)
                            {
                                while ((line = br.readLine()).startsWith("#"))
                                    pw.println(line);
                                pw.println("run\t" + line);
                            }
                            else
                            {
                                while ((line = br.readLine()).startsWith("#"))
                                    continue;
                                //now at header line
                            }
                            pw.flush();
                            //next read is past the header line
                            while ((line = br.readLine()) != null)
                            {
                                String baseName = MS2ExtraInfoDef.getFeatureSetBaseName(fSet);
                                if (allHaveBaseNames)
                                    pw.println(baseName + "\t" + line);
                                else
                                    pw.println((i+1) + "\t" + line);
                                pw.flush();
                            }
                        }
                        TempFileManager.deleteTempFiles(this);
                    }
                    break;
                case FILE_FORMAT_PEPXML:
                    if (forcePeptideProphet)
                    {
                        for(Feature feature : featureSet.getFeatures())
                            MS2ExtraInfoDef.setPeptideProphet(feature, forcePeptideProphetValue);
                    }
                    FeaturePepXmlWriter pepXmlWriter =
                            new FeaturePepXmlWriter(featureSet);
                    pepXmlWriter.set_searchEngine(pepXmlSearchEngine);
                    if (fastaFile != null)
                        pepXmlWriter.setSearchDatabase(fastaFile.getAbsolutePath());
                    pepXmlWriter.write(outputFile);
                    break;
                case FILE_FORMAT_APML:
                    APMLFeatureFileHandler.getSingletonInstance().saveFeatureSet(featureSet, outputFile);
                    break;
                case FILE_FORMAT_HARDKLOR:
                    featureSet.save(outputFile,dumpWindow, HardklorFeatureFileHandler.FILE_TYPE_NAME);
                    break;
                default:
                    throw new CommandLineModuleExecutionException("Don't know how to support the specified output file format yet");
            }
            ApplicationContext.infoMessage("Successfully wrote feature file " + outputFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (pw != null)
                pw.close();
        }

    }

    protected FeatureSet loadFeatureSetFromSpecArrayTSV(File featureFile, int specArrayVersion) throws CommandLineModuleExecutionException
    {
        FeatureSet result = null;
        try
        {
            FileInputStream fis = new FileInputStream(featureFile);
            ArrayList<Feature> featureList = new ArrayList<Feature>();

            String fileLine;
            double[] timesForScans = new double[run.getScanCount()];
            for (int i=0; i<timesForScans.length; i++)
            {
                timesForScans[i] = run.getScan(i).getDoubleRetentionTime();
            }
            while ((fileLine = readLine(fis)) != null)
            {
                if (fileLine.startsWith("index") || fileLine.contains("intensity"))
                    continue;
                featureList.add(loadFeatureFromSpecArrayLine(fileLine, timesForScans, specArrayVersion));
            }
            result = new FeatureSet(featureList.toArray(new Feature[0]));
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        return result;
    }

    protected Feature loadFeatureFromSpecArrayLine(String specArrayLine, double[] timesForScans,
                                                   int specArrayVersion)
    {
        Feature result = new Feature();

        int specArrayCharge = 0;
        double specArrayTime = 0;
        double specArrayIntensity = 0;

        String[] specArrayLineArray = specArrayLine.split(" ");


        switch (specArrayVersion)
        {
            case SPECARRAY_VERSION_1_0:
                double specArrayMass = Double.parseDouble(specArrayLineArray[1]);
                specArrayCharge = Integer.parseInt(specArrayLineArray[2]);
                specArrayTime = Double.parseDouble(specArrayLineArray[11]);
                specArrayIntensity = Double.parseDouble(specArrayLineArray[9]);

                result.setMass((float) specArrayMass);
                break;
            case SPECARRAY_VERSION_1_2_0:
                double specArrayMz = Double.parseDouble(specArrayLineArray[0]);
                specArrayTime = Double.parseDouble(specArrayLineArray[1]);
                specArrayCharge = Integer.parseInt(specArrayLineArray[3]);
                specArrayIntensity = Double.parseDouble(specArrayLineArray[4]);

                result.setMz((float) specArrayMz);
                break;
        }

        result.setCharge(specArrayCharge);
        result.afterPopulate();

        float time = (float) specArrayTime * 60;
        result.setTime(time);

        float totalIntensity = (float) specArrayIntensity;
        result.setIntensity(totalIntensity);

        int scan = Arrays.binarySearch(timesForScans,time);
        if (scan < 0)
        {
            scan = -scan;
            if (Math.abs(time - timesForScans[scan-1]) < Math.abs(time - timesForScans[scan]))
                scan -= 1;
        }
        result.setScan(scan);


        return result;
    }

    protected String readLine(FileInputStream fis) throws IOException
    {
        String result = null;

        StringBuffer resultBuf = new StringBuffer();
        int charread;
        while ((charread = fis.read()) != -1 && charread != '\n')
        {
            resultBuf.append((char)charread);
        }
        if (resultBuf.length() > 0)
        {
            result = resultBuf.toString();
        }
        return result;
    }

    /**
     * Create the set of intensity windows around a given set of features that allow us to
     * dump those intensity windows to a feature file
     *
     * TODO: can we get rid of this completely?
     */
    public static void createIntensityWindows(Feature[] features, MSRun run)
    {
        Feature[] f = new Feature[features.length];
        System.arraycopy(features, 0, f, 0, f.length);
        Feature.ScanAscComparator sac = new Feature.ScanAscComparator();
        Arrays.sort(f, sac);

        MSRun.MSScan scan = null;
        float[][] spectrum = null;
        for (int i = 0; i < f.length; i++)
        {
            // If feature intensity window was already extracted, just skip.
            if (f[i].intensityLeadingPeaks == 3 && f[i].intensityTrailingPeaks == 3)
                continue;

            int n = run.getIndexForScanNum(f[i].scan);
            // This can happen if user, in error, applies feature set from one run to another.
            if (n >= run.getScanCount() || n < 0)
                continue;
            scan = run.getScan(n);
            spectrum = scan.getSpectrum();
            if (null == spectrum)
            {
                _log.error("Failed to get spectrum for scan " + f[i].scan);
                ApplicationContext.setMessage("Failed to get spectrum for scan " + f[i].scan);
                return;
            }

            f[i].intensityWindow =
                    Spectrum.Resample(spectrum, new FloatRange(f[i].mz - 3, f[i].mz + 3),
                                      SpectrumResampler.getResampleFrequency());
            f[i].intensityLeadingPeaks = 3;
            f[i].intensityTrailingPeaks = 3;
        }
        ApplicationContext.setMessage("");
    }

}
