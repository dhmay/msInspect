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
package org.fhcrc.cpl.viewer.ms2.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.List;
import java.util.ArrayList;


/**
 * test
 */
public class GuessProteinsFromFastaCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(GuessProteinsFromFastaCLM.class);

    protected Protein[] fastaProteins = null;

    protected File[] featureFiles;
    protected FeatureSet[] featureSets;
    protected File fastaFile;

    protected File outFile = null;
    protected File outDir = null;

    protected boolean runRefreshParser = false;

    protected boolean stripHighCharge = true;

    protected static final int OUT_FORMAT_MSINSPECT = 0;
    protected static final int OUT_FORMAT_PEPXML = 1;

    protected boolean guessAllProteins = false;


    protected static String[] outFormatStrings = new String[]
            {
                    "msinspect",
                    "pepxml"
            };

    protected static String[] outFormatExplanations = new String[]
            {
                    "msinspect",
                    "pepxml"
            };

    protected int outFormat = OUT_FORMAT_PEPXML;


    public GuessProteinsFromFastaCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "guessproteinsfromfasta";

        mHelpMessage ="Populate (guess at) the protein for each identified peptide";
        mShortDescription = "Populate (guess at) the protein for each identified peptide, by "
                 + "looking at the specified FASTA and using the first protein that contains the peptide";

        CommandLineArgumentDefinition[] argDefs =
               {
                       createUnnamedSeriesFileArgumentDefinition(true, "Feature File to fix up"),
                       new FileToReadArgumentDefinition("fasta", true, "Fasta file"),
                       new FileToWriteArgumentDefinition("out", false, "output file (if not specified, alters in place"),
                       new DirectoryToWriteArgumentDefinition("outdir", false,
                               "output directory (if not specified, alters in place"),
                       new BooleanArgumentDefinition("refreshparser", false,
                               "Run RefreshParser?  RefreshParser executable must be on path.", runRefreshParser),
                       new BooleanArgumentDefinition("striphighcharge", false,
                               "Strip high-charge features from output? (for ProteinProphet)", stripHighCharge),
                       new EnumeratedValuesArgumentDefinition("outformat", false,
                               outFormatStrings, outFormatExplanations, "pepxml"),
                       new BooleanArgumentDefinition("guessallproteins", false,
                               "Guess all proteins?  If false, just guess one protein", guessAllProteins),
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        File[] featureFiles = getUnnamedSeriesFileArgumentValues();
        featureSets = new FeatureSet[featureFiles.length];

        fastaFile = getFileArgumentValue("fasta");

        outFormat = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("outformat")).getIndexForArgumentValue(getStringArgumentValue("outformat"));

        guessAllProteins = getBooleanArgumentValue("guessallproteins");

        stripHighCharge = getBooleanArgumentValue("striphighcharge");

        for (int i=0; i<featureSets.length; i++)
        {
            _log.debug("Loading feature file " + featureFiles[i].getName());
            try
            {
                featureSets[i] = new FeatureSet(featureFiles[i]);
                MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSets[i], fastaFile.getAbsolutePath());

                if (stripHighCharge)
                {
                    FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
                    sel.setMaxCharge(5);
                    featureSets[i] = featureSets[i].filter(sel);
                }
            }
            catch (Exception e)
            {
                throw new ArgumentValidationException(e);
            }
        }

        fastaProteins = ProteinUtilities.loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);

        runRefreshParser = getBooleanArgumentValue("refreshparser");

        if (outFormat != OUT_FORMAT_PEPXML && runRefreshParser)
            throw new ArgumentValidationException("Can't run refreshparser if not outputting pepXML");

        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            if (guessAllProteins)
                ProteinUtilities.guessAllProteinsForFeaturePeptides(featureSets, fastaFile, fastaProteins);
            else
                ProteinUtilities.guessProteinsForFeaturePeptides(featureSets, fastaFile, fastaProteins);

            List<File> outputFileList = new ArrayList<File>();

            for (FeatureSet featureSet : featureSets)
            {
                ApplicationContext.setMessage("Processing file " +
                        featureSet.getSourceFile().getAbsolutePath());

                File outputFile = null;
                if (outFile != null)
                {
                    outputFile = outFile;
                }
                else if (outDir != null)
                {
                    outputFile = new File(outDir, featureSet.getSourceFile().getName());
                }
                else
                    outputFile = featureSet.getSourceFile();

                switch(outFormat)
                {
                    case OUT_FORMAT_PEPXML:
                        featureSet.savePepXml(outputFile);
                        break;
                    default:
                        featureSet.save(outputFile);
                        break;
                }
                outputFileList.add(outputFile);
                ApplicationContext.setMessage("Saved file " + outputFile);
            }


            if (runRefreshParser)
            {
                ApplicationContext.setMessage("Running RefreshParser on all files...");
                
                for (File file : outputFileList)
                {
                    String cmd = "RefreshParser " + file.getAbsolutePath() + " " + fastaFile.getAbsolutePath();
                    _log.debug("Running RefreshParser: " + cmd);
                    Process p = Runtime.getRuntime().exec(cmd ,null);
                    int err = p.waitFor();
                    _log.debug("process returned, "+err);
                    if (err == 0)
                        ApplicationContext.setMessage("Successfully ran RefreshParser on " +
                                file.getAbsolutePath());
                    else
                        ApplicationContext.setMessage("RefreshParser failed on " +
                                file.getAbsolutePath() + " with error code " + err);
                }
            }
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
