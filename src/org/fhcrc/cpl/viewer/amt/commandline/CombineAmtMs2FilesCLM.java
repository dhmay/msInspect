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
package org.fhcrc.cpl.viewer.amt.commandline;

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.viewer.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.extraInfo.*;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.ms2.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Pair;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;


/**
 * Command linemodule for feature finding
 */
public class CombineAmtMs2FilesCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CombineAmtMs2FilesCLM.class);


    protected Protein[] fastaProteins = null;

    protected File outFile;
    protected File outDir;

    protected File[] inFiles;
    protected File inDir;
    protected File inDir2;

    protected boolean guessProteins = false;
    protected File fastaFile = null;


    protected int outFormat=FORMAT_PEPXML;


    protected static final int FORMAT_PEPXML = 0;
    protected static final int FORMAT_FEATURE = 1;


    protected boolean showCharts = false;

    protected boolean runRefreshParser = false;

    protected boolean restrictCharge = true;



    protected final static String[] formatStrings =
            {
                    "pepxml",
                    "feature"
            };

    protected final static String[] formatExplanations =
            {
                    "PepXML",
                    "feature"
            };


    public CombineAmtMs2FilesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "combineamtms2";
        mShortDescription =
                "Combine MS/MS results (pepxml or feature files) with AMT results (pepxml or feature files)";
        mHelpMessage = "Combine MS/MS results (pepxml or feature files) with AMT results (pepxml or feature files)";

        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesArgumentDefinition(ArgumentDefinitionFactory.FILE_TO_READ,false, null),
                        createFileToWriteArgumentDefinition("out",false, null),
                        createDirectoryToReadArgumentDefinition("outdir",false, null),
                        createDirectoryToReadArgumentDefinition("indir", false, "Directory of input files"),
                        createDirectoryToReadArgumentDefinition("indir2", false, "Directory of input files 2"),
                        createEnumeratedArgumentDefinition("outformat", false, "Output format", formatStrings, "pepxml"),
                        createBooleanArgumentDefinition("guessproteins", false,
                                "Guess proteins for peptides? (Requires fasta)", guessProteins),
                        createFileToReadArgumentDefinition("fasta", false, "Fasta file"),
                        createBooleanArgumentDefinition("refreshparser", false,
                                "Run RefreshParser? (RefreshParser must be on path", runRefreshParser),
                        createBooleanArgumentDefinition("restrictcharge", false,
                                "Cap feature charge in output files at 5? (the maximum allowed by ProteinProphet)"),
                };
        addArgumentDefinitions(argDefs);

    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        if (hasUnnamedSeriesArgumentValue())
        {
            assertArgumentAbsent("indir");
            inFiles = getUnnamedSeriesFileArgumentValues();
        }
        else
        {
            assertArgumentPresent("indir");
            inDir = getFileArgumentValue("indir");
//            File[] allFilesInDir = getFileArgumentValue("indir").listFiles();
//            List<File> inFilesList = new ArrayList<File>();
//            for (File file : allFilesInDir)
//                if (!file.isDirectory())
//                    inFilesList.add(file);
//            inFiles = inFilesList.toArray(new File[inFilesList.size()]);
        }


        inDir2 = getFileArgumentValue("indir2");

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (hasArgumentValue("outformat"))
            outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(getStringArgumentValue("outformat"));

        guessProteins = getBooleanArgumentValue("guessproteins");
        if (hasArgumentValue("fasta"))
            fastaFile = getFileArgumentValue("fasta");        
        if (guessProteins)
        {
            assertArgumentPresent("fasta");
            fastaProteins = ProteinUtilities.loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);
        }

        runRefreshParser = getBooleanArgumentValue("refreshparser");



    }




    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (inFiles != null)
            handleFiles(inFiles, outFile);
        else
        {
            List<Pair<File,File>> filePairs = findFilePairs();
            for (Pair<File,File> filePair : filePairs)
            {
                ApplicationContext.infoMessage("\tprocessing file pair " + filePair.first.getName() + ", " + filePair.second.getName());
                File[] fileArray = new File[] { filePair.first, filePair.second };
                String outFileName =
                        filePair.first.getName().substring(0, filePair.first.getName().lastIndexOf(".") + 1);

                switch (outFormat)
                {
                    case FORMAT_PEPXML:
                        outFileName = outFileName + "pep.xml";
                        break;
                    case FORMAT_FEATURE:
                        outFileName = outFileName + "tsv";
                        break;
                }
                handleFiles(fileArray, new File(outDir, outFileName));
            }
        }

    }

    protected void handleFiles(File[] files, File outputFile)
            throws CommandLineModuleExecutionException
    {
        FeatureSet combinedFeatureSet = createCombinedSet(files, outputFile);

        if (guessProteins)
        {
            ApplicationContext.infoMessage("\tGuessing proteins for " + outputFile.getAbsolutePath());

            ProteinUtilities.guessProteinsForFeaturePeptides(
                    combinedFeatureSet,
                    fastaProteins);

            ApplicationContext.setMessage("\tDone.");
        }




        switch (outFormat)
        {
            case FORMAT_PEPXML:
                try
                {
                    combinedFeatureSet.savePepXml(outputFile);
                    ApplicationContext.setMessage("Saved output pepXML file " + outputFile.getAbsolutePath() +
                            " with " + combinedFeatureSet.getFeatures().length + " total features");
                    Set<String> outPeptides = new HashSet<String>();
                    Set<String> outPeptidesWithRatios = new HashSet<String>();


                    for (Feature feature : combinedFeatureSet.getFeatures())
                    {
                        outPeptides.add(MS2ExtraInfoDef.getFirstPeptide(feature));
                        if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                            outPeptidesWithRatios.add(MS2ExtraInfoDef.getFirstPeptide(feature));
                    }
                    ApplicationContext.setMessage("Peptides in output file: " + outPeptides.size() +
                            ".  With ratios: " + outPeptidesWithRatios.size());

                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failed to save pepXML file " + outputFile.getAbsolutePath(),e);
                }
                break;
            case FORMAT_FEATURE:
                try
                {
                    combinedFeatureSet.save(outputFile);
                    ApplicationContext.setMessage("Saved output feature file " + outputFile.getAbsolutePath() +
                            " with " + combinedFeatureSet.getFeatures().length + " total features");
                }
                catch (IOException e)
                {
                    throw new CommandLineModuleExecutionException("Failed to save feature file " +
                            outputFile.getAbsolutePath(),e);
                }
                break;
        }

        try
        {
            if (runRefreshParser)
            {
                ApplicationContext.setMessage("\tRunning RefreshParser.  Command:");
                String cmd = "RefreshParser " + outputFile.getAbsolutePath() + " " + fastaFile.getAbsolutePath();
                ApplicationContext.setMessage("\t\t" + cmd);
                Process p = Runtime.getRuntime().exec(cmd ,null);
                int err = p.waitFor();
                _log.debug("process returned, "+err);
                if (err == 0)
                    ApplicationContext.setMessage("Successfully ran RefreshParser on " +
                            outputFile.getAbsolutePath());
                else
                    ApplicationContext.setMessage("RefreshParser failed on " +
                            outputFile.getAbsolutePath() + " with error code " + err);
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed while running RefreshParser on file " +
                    outputFile, e);
        }
    }





    /**
     * Pick pairs of files from directories by name, without extension
     * @return
     */
    protected List<Pair<File,File>> findFilePairs()
    {
        List<Pair<File,File>> result = new ArrayList<Pair<File,File>>();

        for (File file1 : inDir.listFiles())
        {
            if (file1.isDirectory())
                continue;
            String file1Name = file1.getName();
            String file1NamePrefix = file1Name.substring(0, file1Name.indexOf("."));
            for (File file2 : inDir2.listFiles())
            {
                if (file2.isDirectory())
                    continue;
                String file2NamePrefix = file2.getName().substring(0, file2.getName().indexOf("."));

                if (file1NamePrefix.equals(file2NamePrefix))
                {
                    result.add(new Pair<File,File>(file1,file2));
                    break;
                }
            }
        }
        ApplicationContext.infoMessage("Found " + result.size() + " pairs of files");
        return result;
    }



    protected FeatureSet createCombinedSet(File[] files, File outputFile) throws CommandLineModuleExecutionException
    {
        List<MS2Modification> modsForOutput = new ArrayList<MS2Modification>();
        Set<Feature> allFeatures = new HashSet<Feature>();

        //we can only keep one value for minTermini and maxCleavages, so keep the first nonzero one
        //for each
        int minTermini = 0;
        int maxCleavages = 0;

        for (int i=0; i<files.length; i++)
        {
            File file = files[i];
            ApplicationContext.infoMessage("Checking file " + i + ", " + files[i].getAbsolutePath());
            BufferedReader br = null;
            try
            {
                FeatureSet featureSet = new FeatureSet(files[i]);
                if (_log.isDebugEnabled())
                {
                    Set<String> allPeptides = new HashSet<String>();
                    Set<String> quantifiedPeptides = new HashSet<String>();
                    for (Feature feature : featureSet.getFeatures())
                    {
                        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                        if (peptide != null)
                        {
                            allPeptides.add(peptide);
                            if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                                quantifiedPeptides.add(peptide);
                        }
                    }
                    _log.debug("\tFile " + i + ": peptides: " + allPeptides.size() +
                            ", quantified: " + quantifiedPeptides.size());
                }
                if (minTermini == 0)
                    minTermini = MS2ExtraInfoDef.getFeatureSetSearchConstraintMinTermini(featureSet);
                if (maxCleavages == 0)
                    maxCleavages = MS2ExtraInfoDef.getFeatureSetSearchConstraintMaxIntCleavages(featureSet);

                MS2Modification[] featureSetMods =
                        MS2ExtraInfoDef.getFeatureSetModifications(featureSet);
                if (featureSetMods != null)
                {
                    for (MS2Modification newMod : featureSetMods)
                    {
                        boolean foundDuplicate = false;
                        MS2Modification dupMod = null;
                        for (MS2Modification existingMod : modsForOutput)
                        {
                            if (existingMod.getAminoAcid().equals(newMod.getAminoAcid()) &&
                                    (existingMod.getVariable() && newMod.getVariable() ||
                                            !existingMod.getVariable() && !newMod.getVariable()) &&
                                    Math.abs(existingMod.getMassDiff() - newMod.getMassDiff()) < .1)
                            {
                                foundDuplicate = true;
                                dupMod = existingMod;
                                break;
                            }
                        }
                        if (foundDuplicate)
                        {
                            _log.debug("Duplicate mods found.  Copy 1: " + dupMod +
                                    ", Copy 2: " + newMod);
                        }
                        else
                        {
                            _log.debug("Adding modification for output: " + newMod);
                            modsForOutput.add(newMod);
                        }
                    }
                }

                for (Feature feature : featureSet.getFeatures())
                    allFeatures.add(feature);
            }
            catch (Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
            finally
            {
                if (br != null) try { br.close();} catch (Exception e) {}
            }
            ApplicationContext.setMessage("Loaded features from file " + file.getName());
        }

        List<Feature> outFeatures = new ArrayList<Feature>();
        outFeatures.addAll(allFeatures);

        FeatureSet outFeatureSet = new FeatureSet(outFeatures.toArray(new Feature[outFeatures.size()]));
        MS2ExtraInfoDef.setFeatureSetSearchConstraintMaxIntCleavages(outFeatureSet, maxCleavages);
        MS2ExtraInfoDef.setFeatureSetSearchConstraintMinTermini(outFeatureSet, minTermini);
        _log.debug("Max cleavages: " + maxCleavages + ", Min Termini: " + minTermini);        

        if (modsForOutput.size() > 0)
        {
            MS2ExtraInfoDef.setFeatureSetModifications(outFeatureSet,
                    modsForOutput.toArray(new MS2Modification[modsForOutput.size()]));
        }
        _log.debug("Added " + modsForOutput.size() + " modifications to output file");

        outFeatureSet.setSourceFile(outputFile);

        if (restrictCharge)
        {
            FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
            sel.setMaxCharge(5);

            outFeatureSet = outFeatureSet.filter(sel);
        }

        if (fastaFile != null)
            MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(outFeatureSet, fastaFile.getAbsolutePath());

        return outFeatureSet;

    }




}