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

import org.fhcrc.cpl.viewer.commandline.modules.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.commandline.arguments.*;
import org.fhcrc.cpl.viewer.amt.*;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;
import org.apache.log4j.Logger;

import java.io.File;


/**
 * Command linemodule for refining AMT databases
 */
public class AmtDatabaseManagerCommandLineModule extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AmtDatabaseManagerCommandLineModule.class);

    protected AmtDatabase amtDatabase = null;
    protected File outFile = null;

    boolean fromAcrylamideToNot = true;


    protected int mode=-1;


    protected static final int REMOVE_OUTLIER_OBSERVATIONS_MODE=0;
    protected static final int REMOVE_PREDICTED_H_OUTLIERS_MODE=1;
    protected static final int REMOVE_FEW_OBSERVATIONS_MODE=2;
    protected static final int REMOVE_RUNS_WITHOUT_MASS_MATCHES_MODE=3;
    protected static final int ALIGN_ALL_RUNS_MODE=4;
    protected static final int ADJUST_ACRYLAMIDE_MODE=5;
    protected static final int REMOVE_PEPTIDES_WITH_RESIDUE_MODE=6;


    protected float predictedHOutlierDeviationMultipleCutoff = .001f;


    protected final static String[] modeStrings = {
            "removeoutlierobservations",
            "removepredictedhoutliers",
            "removefewobservations",
            "removerunswithoutmassmatches",
            "alignallruns",
            "adjustacrylamide",
            "removepeptideswithresidue",
    };

    protected static final String[] modeExplanations =
            {
                    "Remove all individual observations that are at least 3 standard deviations from the median for that peptide, for peptides with at least three observations",
                    "Remove all peptides with only one observation, for which that observation is at least 2 standard deviations away from the prediction",
                    "Remove all peptides with fewer than minobservations observations",
                    "Make mass matches between each run's entries and an MS1 feature file.  Remove all runs that don't mass-match at least minmassmatchpercent percent of peptides",
                    "Nonlinearly align all runs in the database to a single run, starting with the run with the " +
                            "most peptide overlap with other runs in the database.  This is an extremely important step.",
                    "Adjust all Cysteine-bearing observations to take the H contribution of acrylamide into account.  Direction of adjustment depends on the 'fromacryltonot' parameter",
                    "Remove all peptides containing a given residue"
            };

    protected int minObservations = 2;
    protected int minMassMatchPercent = 20;
    protected int maxEntriesInMassMatchedDatabase = Integer.MAX_VALUE;
    protected int maxRunsInMassMatchedDatabase = Integer.MAX_VALUE;

    protected FeatureSet ms1FeatureSet = null;
    protected boolean showCharts = false;


    protected String residueToRemove = null;

    //todo: parameterize
    protected float massMatchDeltaMass = AmtDatabaseMatcher.DEFAULT_MASS_MATCH_DELTA_MASS;
    protected int massMatchDeltaMassType = AmtDatabaseMatcher.DEFAULT_MASS_MATCH_DELTA_MASS_TYPE;

    public AmtDatabaseManagerCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "manageamt";
        mShortDescription = "Tools for managing an AMT database";
        mHelpMessage = "Refine an AMT database by removing outlier observations or peptides, " +
                "or nonlinearly aligning all runs to each other";
        CommandLineArgumentDefinition[] basicArgDefs =
                {
                        createEnumeratedArgumentDefinition("mode",true,modeStrings, modeExplanations),
                        createFileToWriteArgumentDefinition("out", false, null),
                        createUnnamedFileArgumentDefinition(true, "AMT database file"),
                        createIntegerArgumentDefinition("minobservations", false,
                                "Minimum number of observations for features kept in the database",
                                minObservations),
                        createBooleanArgumentDefinition("showcharts", false,
                                "Show charts?", showCharts),
                        createStringArgumentDefinition("residue", false,
                                "Residue (for 'removepeptideswithresidue' mode)")
                };
        addArgumentDefinitions(basicArgDefs);

        CommandLineArgumentDefinition[] advancedArgDefs =
                {
                        createIntegerArgumentDefinition("minmassmatchpercent", false,
                                "Minimum percent of peptides mass-matched to MS1, per run " +
                                        "(removerunswithoutmassmatches mode only)",
                                minMassMatchPercent),
                        createIntegerArgumentDefinition("maxentries", false,
                                "Maximum DB entries (removerunswithoutmassmatches mode only)",
                                maxEntriesInMassMatchedDatabase),
                        createIntegerArgumentDefinition("maxruns", false,
                                "Maximum DB runs (removerunswithoutmassmatches mode only)",
                                maxRunsInMassMatchedDatabase),
                        createFeatureFileArgumentDefinition("ms1features", false,
                                "MS1 features (removerunswithoutmassmatches mode only)"),
                        createBooleanArgumentDefinition("fromacryltonot", false,
                                "For mode adjustacrylamide.  If true, adjusts all observations to _remove_ the effect " +
                                        " of acrylamide.  If false, adjusts observations to _add_ the effect.")
                };
        addArgumentDefinitions(advancedArgDefs, true);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mode = ((EnumeratedValuesArgumentDefinition) getArgumentDefinition("mode")).getIndexForArgumentValue(getStringArgumentValue("mode"));

        File dbFile = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
        try
        {
            AmtXmlReader amtXmlReader = new AmtXmlReader(dbFile);
            amtDatabase = amtXmlReader.getDatabase();
            _log.info("Loaded AMT database: " + amtDatabase);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(e);
        }

        if (mode == REMOVE_FEW_OBSERVATIONS_MODE)
        {
            assertArgumentPresent("minobservations");
            minObservations = getIntegerArgumentValue("minobservations");
        }
        else
            assertArgumentAbsent("minobservations");

        if (mode == REMOVE_RUNS_WITHOUT_MASS_MATCHES_MODE)
        {
            assertArgumentPresent("minmassmatchpercent");
            minMassMatchPercent = getIntegerArgumentValue("minmassmatchpercent");
            assertArgumentPresent("ms1features");
            ms1FeatureSet = getFeatureSetArgumentValue("ms1features");
            assertArgumentPresent("maxentries");
            maxEntriesInMassMatchedDatabase = getIntegerArgumentValue("maxentries");
            maxRunsInMassMatchedDatabase = getIntegerArgumentValue("maxruns");                       
        }
        else
        {
            assertArgumentAbsent("minmassmatchpercent");
            assertArgumentAbsent("ms1features");
        }

        if (mode == ADJUST_ACRYLAMIDE_MODE)
        {
            assertArgumentPresent("fromacryltonot");
            fromAcrylamideToNot = getBooleanArgumentValue("fromacryltonot");
        }
        else
        {
            assertArgumentAbsent("fromacryltonot");
        }

        residueToRemove = getStringArgumentValue("residue");
        if (mode == REMOVE_PEPTIDES_WITH_RESIDUE_MODE)
            assertArgumentPresent("residue","mode");

        outFile = getFileArgumentValue("out");

        showCharts = getBooleanArgumentValue("showcharts");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Read AMT Database with " + amtDatabase.numEntries() + " entries.");
        switch (mode)
        {
            case REMOVE_OUTLIER_OBSERVATIONS_MODE:
            {
                int numObsBefore = 0;
                for (AmtRunEntry runEntry : amtDatabase.getRuns())
                    numObsBefore += amtDatabase.getObservationsForRun(runEntry).length;

                //todo: parameterize
                AmtDatabaseManager.removeHydrophobicityOutliers(amtDatabase, 3);

                int numObsAfter = 0;
                for (AmtRunEntry runEntry : amtDatabase.getRuns())
                    numObsAfter += amtDatabase.getObservationsForRun(runEntry).length;
                ApplicationContext.infoMessage("\nRemoved " + (numObsBefore - numObsAfter) +
                        " observations (out of " + numObsBefore + ")");


                break;
            }
            case REMOVE_PREDICTED_H_OUTLIERS_MODE:
            {
                //todo: parameterize
                double maxHDiff = amtDatabase.calculateMeanDifferenceFromPredictedHydro() +
                        (predictedHOutlierDeviationMultipleCutoff * amtDatabase.calculateStandardDeviationDifferenceFromPredictedHydro());
                ApplicationContext.infoMessage("Removing entries with one observation and observed hydro > " +
                        maxHDiff + " different from predicted");
                AmtDatabaseManager.removePredictedHOutliers(amtDatabase, (float) maxHDiff);


                break;
            }
            case REMOVE_FEW_OBSERVATIONS_MODE:
            {
                ApplicationContext.infoMessage("Removing entries with less than " +
                        minObservations + " observations");
                for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
                {
                    if (peptideEntry.getNumObservations() < minObservations)
                        amtDatabase.removeEntry(peptideEntry.getPeptideSequence());
                }
                break;
            }
            case REMOVE_RUNS_WITHOUT_MASS_MATCHES_MODE:
            {
                amtDatabase =
                        AmtDatabaseManager.removeRunsWithoutMassMatches(
                                amtDatabase, ms1FeatureSet.getFeatures(),
                                minMassMatchPercent, massMatchDeltaMass,
                                massMatchDeltaMassType,
                                maxEntriesInMassMatchedDatabase,
                                maxRunsInMassMatchedDatabase,
                                AmtDatabaseMatcherCLM.defaultMS2ModificationsForMatching,
                                showCharts);
                break;
            }
            case ALIGN_ALL_RUNS_MODE:
            {
                amtDatabase =
                        AmtDatabaseManager.alignAllRunsUsingCommonPeptides(
                                amtDatabase, 10, showCharts);
//                ApplicationContext.infoMessage("Repeating...");
//                amtDatabase =
//                        AmtDatabaseManager.alignAllRunsUsingCommonPeptides(amtDatabase, 10, showCharts);
                break;
            }
            case ADJUST_ACRYLAMIDE_MODE:
            {
                amtDatabase = AmtDatabaseManager.adjustEntriesForAcrylamide(
                        amtDatabase, fromAcrylamideToNot, showCharts);

                break;
            }
            case REMOVE_PEPTIDES_WITH_RESIDUE_MODE:
            {
                ApplicationContext.infoMessage("Removing all peptides containing '" + residueToRemove + "'...");
                for (AmtPeptideEntry peptideEntry : amtDatabase.getEntries())
                {
                    if (peptideEntry.getPeptideSequence().contains(residueToRemove))
                        amtDatabase.removeEntry(peptideEntry.getPeptideSequence());
                }
            }
        }
        if (outFile != null && amtDatabase != null)
            writeAmtDatabase(amtDatabase, outFile);
    }





    /**
     *
     * @param amtDatabase
     * @param outAmtXmlFile
     */
    protected static void writeAmtDatabase(AmtDatabase amtDatabase,
                                           File outAmtXmlFile)
    {
        try
        {
            AmtXmlWriter amtXmlWriter = new AmtXmlWriter(amtDatabase);
            amtXmlWriter.write(outAmtXmlFile);
            ApplicationContext.infoMessage("Wrote " +
                    amtDatabase.numEntries() + " entries to amtxml file " +
                    outAmtXmlFile.getAbsolutePath());
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            ApplicationContext.infoMessage("Error writing amt file " +
                    outAmtXmlFile.getAbsolutePath());
        }
    }


}
