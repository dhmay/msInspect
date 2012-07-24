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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.PeptideUtilities;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.amt.AmtPeptideEntry;
import org.fhcrc.cpl.viewer.amt.AmtDatabaseFeatureSetGenerator;
import org.apache.log4j.Logger;



import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;


/**
 * Calculates mz values for peptides for a targeted MS/MS experiment
 */
public class CalcPeptideTargetMzsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CalcPeptideTargetMzsCLM.class);


    MS2Modification[] mods = new MS2Modification[] {
            ModificationListArgumentDefinition.safeCreateModificationFromString("C57.021"),
            ModificationListArgumentDefinition.safeCreateModificationFromString("M16V")
    };

    protected File outFile = null;

    List<String> peptides = new ArrayList<String>();

    int minCharge = 2;
    int maxCharge = 3;

    protected float massToAdd = 0f;


    public CalcPeptideTargetMzsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "calctargetmzs";
        mShortDescription =
                "Calculate mz (or mass) values for targeted MS/MS of specified peptides";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        new StringListArgumentDefinition("peptides",false, "comma-separated list of peptides"),
                        new FileToReadArgumentDefinition("peptidefile", false, "peptides in a file, one per line"),
                        new FileToWriteArgumentDefinition("out", true, null),
                        new IntegerArgumentDefinition("mincharge", false, "Minimum charge", minCharge),
                        new IntegerArgumentDefinition("maxcharge", false, "Maximum charge", maxCharge),
                        new ModificationListArgumentDefinition("mods", false,
                                "Static modifications (default C57.021, M16V)"),
                        new DecimalArgumentDefinition("masstoadd", false,
                                "A mass to add to all peptide mass values (e.g., for decoy database)", massToAdd),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        peptides = getStringListArgumentValue("peptides");
        if (peptides==null)
            assertArgumentPresent("peptidefile");
        else
            assertArgumentAbsent("peptidefile");

        if (hasArgumentValue("peptidefile")) {
            try {
                peptides = new ArrayList<String>();
                BufferedReader br = new BufferedReader(new FileReader(getFileArgumentValue("peptidefile")));
                String peptide = null;
                while ((peptide = br.readLine()) != null) {
                     peptides.add(peptide);
                }
            } catch (IOException e) {
                throw new ArgumentValidationException("Failed to read peptide file " +
                        getFileArgumentValue("peptidefile").getAbsolutePath());
            }
            _log.info("Loaded " + peptides.size() + " peptides from file.");
        }

        outFile = getFileArgumentValue("out");
        minCharge = getIntegerArgumentValue("mincharge");
        maxCharge = getIntegerArgumentValue("maxcharge");

        massToAdd = getFloatArgumentValue("masstoadd");
        if (massToAdd > 0)
            _log.info("WARNING: adding " + massToAdd + " to all peptide masses.");

        if (hasArgumentValue("mods"))
            mods = (MS2Modification[]) getArgumentValue("mods");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<MS2Modification> modsList = new ArrayList<MS2Modification>(Arrays.asList(mods));
        try
        {
            PrintWriter outPW = new PrintWriter(outFile);
            outPW.println("peptide\tmodpeptide\tcharge\tmass\tmz");
            for (String peptide : peptides)
            {
                for (Pair<List<ModifiedAminoAcid>[], Float> modsAndMass :
                        PeptideUtilities.calculatePeptideMassesForMods(peptide, modsList))
                {
                    float mass = modsAndMass.second + massToAdd;
                    for (int charge=minCharge; charge<=maxCharge; charge++)
                    {
                        float mz = Feature.convertMassToMz(mass, charge);
                        outPW.println(peptide + "\t" +
                                MS2ExtraInfoDef.createModifiedSequenceString(peptide, modsAndMass.first)
                           + "\t" + charge + "\t" + mass + "\t" + mz);
                        outPW.flush();
                    }
                }
            }

            outPW.flush();
            outPW.close();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
