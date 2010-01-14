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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.viewer.align.BucketedPeptideArray;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;


/**
 * Command linemodule for feature finding
 */
public class ConsensusFeatureFileCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ConsensusFeatureFileCLM.class);

    protected File[] featureFiles;
    protected int minRunsPerFeature = 2;
    protected File outFile;
    protected boolean shouldRequireSamePeptide = false;

       public ConsensusFeatureFileCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "consensusfeaturefile";

        mShortDescription = "";
        mHelpMessage = "";

        CommandLineArgumentDefinition[] argDefs =
               {
                       new FileToWriteArgumentDefinition("out",false, "output file"),
                       createUnnamedSeriesFileArgumentDefinition(true, "input feature files"),
                       new IntegerArgumentDefinition("minfeatureruns", false,
                               "minimun number of runs a feature must appear in",
                               minRunsPerFeature),
                       new BooleanArgumentDefinition("samepeptide", false, "Require all features in row to have same peptide?", shouldRequireSamePeptide)
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        outFile = getFileArgumentValue("out");
        minRunsPerFeature = getIntegerArgumentValue("minfeatureruns");
        featureFiles = getUnnamedSeriesFileArgumentValues();
        shouldRequireSamePeptide = getBooleanArgumentValue("samepeptide");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        File tempArrayFile = new File(outFile.getParent(), "temp.array.tsv");
        File tempDetailsFile = new File(outFile.getParent(), "temp.array.details.tsv");


        List<File> featureFileList = new ArrayList<File>(featureFiles.length);
        for (File featureFile : featureFiles)
            featureFileList.add(featureFile);
        BucketedPeptideArray arr = new BucketedPeptideArray(featureFileList, new FeatureSet.FeatureSelector());

        arr.setOutFileName(tempArrayFile.getAbsolutePath());
        arr.setAlign(true);
        arr.run(true);

//        arr.getAligner().getWarpingMaps()

        try
        {
            PeptideArrayAnalyzer peptideArrayAnalyzer =
                    new PeptideArrayAnalyzer(tempArrayFile);
            FeatureSet consensusFeatureSet =
                    peptideArrayAnalyzer.createConsensusFeatureSet( tempDetailsFile,
                            minRunsPerFeature, PeptideArrayAnalyzer.CONSENSUS_INTENSITY_MODE_MEDIAN,
                            shouldRequireSamePeptide);
            ApplicationContext.infoMessage("Created consensus feature set with " +
                    consensusFeatureSet.getFeatures().length + " features");
            consensusFeatureSet.save(outFile);
            ApplicationContext.infoMessage("Saved consensus feature set to file " +
                    outFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }



    }



}
