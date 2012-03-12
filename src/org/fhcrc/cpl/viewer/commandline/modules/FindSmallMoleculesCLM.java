/*
 * Copyright (c) 2003-2011 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.viewer.align.commandline.PeptideArrayCommandLineModule;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategySmallMolecule;
import org.fhcrc.cpl.viewer.feature.extraction.strategy.FeatureStrategySmallMoleculeNeg;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Command Module for feature finding.
 *
 * Not much of the actual work of feature finding is done directly in this class.  This class serves mainly to
 * parse user arguments and control the flow of work.  This design philosophy is used throughout Command Modules
 * in the msInspect platform
 *
 * dhmay changing 2009/01/02: adding default filtering of minimum 2 peaks and maximum KL 3.0, as well as a
 * "nofilter" parameter that turns this off
 */
public class FindSmallMoleculesCLM extends FindPeptidesCommandLineModule
        implements CommandLineModule
{
     //used for multiple levels of log messages
    protected static Logger _log = Logger.getLogger(FindSmallMoleculesCLM.class);

    //more appropriate default than the smaller default for peptides
    public static final int DEFAULT_SCAN_WINDOW=10000;

    public static final int DEFAULT_STITCH_ELUTION_SECONDS = 600;

    protected int stitchElutionSeconds = DEFAULT_STITCH_ELUTION_SECONDS;

    public static final int SM_MOL_MINPEAKS = 1;
    public static final float SM_MOL_MAXMZ = 1000;
    public static final float SM_MOL_MINMZ = 100;
    public static final float SM_MOL_MAXKL = 5;
    public static final int SM_MOL_MAXCHARGE = 2;


    public void init() {
        super.init();

        //override default from superclass
        scanWindowSize = FindSmallMoleculesCLM.DEFAULT_SCAN_WINDOW;

        //command name
        mCommandName = "findsmallmolecules";

        //A longer help message
        mHelpMessage = "find small-molecule features, in positive ion mode. This command adjusts several of the " +
                "default parameters from findpeptides and uses the appropriate small-molecule-finding feature strategy " +
                "for positive or negative ion mode";

        //A short (single-sentence) description of this command
        mShortDescription = "Find small-molecule features in an mzXML file, from a positive ion mode run";

        addArgumentDefinition(new BooleanArgumentDefinition("posionmode", false,
                "Positive ion mode? (false: negative ion mode)", true));
        addArgumentDefinition(new IntegerArgumentDefinition("stitchseconds", false,
                "Number of seconds over which to stitch together serially-eluting features with the same mass",
                stitchElutionSeconds));

    }

    public void assignArgumentValues()
            throws ArgumentValidationException {

        super.assignArgumentValues();

        if (!hasArgumentValue("scanwindow"))  {
            scanWindowSize = DEFAULT_SCAN_WINDOW;
        }
        if (!hasArgumentValue("minpeaks"))
            featureSelector.setMinPeaks(SM_MOL_MINPEAKS);
        if (!hasArgumentValue("maxmz"))
            featureSelector.setMaxMz(SM_MOL_MAXMZ);
        if (!hasArgumentValue("minmz"))
            featureSelector.setMinMz(SM_MOL_MINMZ);
        if (!hasArgumentValue("maxkl"))
            featureSelector.setMaxKL(SM_MOL_MAXKL);
        if (!hasArgumentValue("maxcharge"))
            featureSelector.setMaxCharge(SM_MOL_MAXCHARGE);
        //This is hacky. Properly I should remove that argument definition, but superclass code will fail if it's
        //not there
        if (hasArgumentValue("strategy"))
            throw new ArgumentValidationException("'strategy' argument not allowed for this command. " +
                    "Apologies for the confusion.");

        if (getBooleanArgumentValue("posionmode"))
            featureStrategyClass = FeatureStrategySmallMolecule.class;
        else featureStrategyClass = FeatureStrategySmallMoleculeNeg.class;

        stitchElutionSeconds = getIntegerArgumentValue("stitchseconds");
    }

    /**
     * workflow
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        System.err.println("Finding features with strategy " + featureStrategyClass.getName() + "...");

        super.execute();
        System.err.println("Features found. Post-processing...");

        List<File> origOutputFiles = new ArrayList<File>();

        if (outFile != null)
            origOutputFiles.add(outFile);
        else {
            for (File file : mzXmlFiles) {
                origOutputFiles.add(calcOutputFile(file));
            }
        }

        for (File origOutputFile : origOutputFiles) {
            File arrFile = TempFileManager.createTempFile(origOutputFile.getName() + ".peparray.tsv", this);
            PeptideArrayCommandLineModule pepCLM = new PeptideArrayCommandLineModule();
            pepCLM.setInputFiles(new File[] { arrFile });
            pepCLM.setOutFile(arrFile);

            Map<String, String> arrayCLMArgMap = new HashMap<String, String>();
            arrayCLMArgMap.put("align", "false");
            //hardcoded tolerance but doesn't matter -- this is just identity
            arrayCLMArgMap.put("masswindow", "1");
            arrayCLMArgMap.put("masstype", "ppm");

            arrayCLMArgMap.put("elutionwindow", "" + stitchElutionSeconds);
            arrayCLMArgMap.put("elutionmode", "time");
            arrayCLMArgMap.put("normalize", "false");
            arrayCLMArgMap.put(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                    origOutputFile.getAbsolutePath());
            arrayCLMArgMap.put("out",arrFile.getAbsolutePath());


            try {
                pepCLM.digestArguments(arrayCLMArgMap);
            } catch (Exception e) {
                throw new CommandLineModuleExecutionException(e);
            }

            ApplicationContext.infoMessage("Building ion array " + arrFile.getAbsolutePath() + "...");
            pepCLM.execute();
            ApplicationContext.infoMessage("Built ion array " + arrFile.getAbsolutePath() );

            ConsensusFeatureFileCLM consensusCLM = new ConsensusFeatureFileCLM();
            Map<String, String> consensusCLMArgMap = new HashMap<String, String>();

            consensusCLMArgMap.put("peparray", arrFile.getAbsolutePath());
            consensusCLMArgMap.put("out", origOutputFile.getAbsolutePath());
            consensusCLMArgMap.put("intensitymode", "max");
            consensusCLMArgMap.put("minfeatureruns", "1");

            try {
                consensusCLM.digestArguments(consensusCLMArgMap);
            } catch (Exception e) {
                throw new CommandLineModuleExecutionException(e);
            }

            ApplicationContext.infoMessage("Building stitched feature file " + origOutputFile.getAbsolutePath() + "...");
            consensusCLM.execute();
            ApplicationContext.infoMessage("Built stitched feature file " + origOutputFile.getAbsolutePath());
        }
    }
}
