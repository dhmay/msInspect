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
package org.fhcrc.cpl.viewer.metabologna.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MetaboliteDatabaseMatcher;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAdd2HMod;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAddWaterMod;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.metabologna.UnknownHAdditionMod;

import java.io.*;
import java.util.*;

/**
 */
public class MatchMetMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MatchMetMassesCLM.class);

    protected File[] inputFiles;
    protected MetaboliteDatabaseMatcher metDBMatcher;
    protected File outFile;
    protected File outDir;


    protected float massTolerancePPM = 2f;

    protected List<ChemicalCompound> databaseCompoundsByMass = null;

    protected TabLoader loader;

    protected boolean shouldCollapseByMax = false;

    protected boolean shouldFixFeatureMasses = false;

    protected boolean shouldUseBaseMod = false;

    public MatchMetMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matchmetmasses";
        mShortDescription = "Match metabolite masses";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedSeriesFileArgumentDefinition(true, "Input masses spreadsheet file"),
                        new FileToReadArgumentDefinition("metdb", true, "Metabolite database file"),
                        new FileToWriteArgumentDefinition("out", false, "out"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "outdir"),
                        new BooleanArgumentDefinition("fixfeaturemasses", false, "Should fix feature masses so that " +
                                "they are mz * z? (different from peptide mass)", shouldFixFeatureMasses),
                        new DecimalArgumentDefinition("deltappm", false, "Delta mass (ppm)", massTolerancePPM),
                        new BooleanArgumentDefinition("usebasemod", false, "Use base mass, i.e., [M]?", shouldUseBaseMod),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inputFiles = this.getUnnamedSeriesFileArgumentValues();

        shouldFixFeatureMasses = getBooleanArgumentValue("fixfeaturemasses");
        shouldUseBaseMod = getBooleanArgumentValue("usebasemod");

        metDBMatcher = new MetaboliteDatabaseMatcher();
        List<ChemicalModification> mods = new ArrayList<ChemicalModification>();

        //todo: paramterize
        //mods.add(new ReduceDoubleBondAdd2HMod());
        //mods.add(new ReduceDoubleBondAddWaterMod());
        mods.add(new UnknownHAdditionMod());

        metDBMatcher.setChemicalModifications(mods);
        try
        {
            databaseCompoundsByMass = ChemicalCompound.loadCompoundsFromFile(getFileArgumentValue("metdb"), 2);
            Collections.sort(databaseCompoundsByMass, new ChemicalCompound.ComparatorMassAsc());

        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load metabolite database",e);
        }
        metDBMatcher.setDatabaseCompounds(databaseCompoundsByMass);
        metDBMatcher.setUseUnmodifiedAdduct(shouldUseBaseMod);

        massTolerancePPM = getFloatArgumentValue("deltappm");

        outFile = getFileArgumentValue("out");
        outDir= getFileArgumentValue("outdir");

        if (inputFiles.length > 1 && outDir == null)
            throw new ArgumentValidationException("multiple input files, no output directory specified");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        for (File file : inputFiles)
        {
            File outputFile = outFile;
            if (outFile == null)
                outputFile = CommandLineModuleUtilities.createOutputFile(file, ".matched.tsv", outDir);
            processFile(file, outputFile);
        }
    }

    public void processFile(File massesFile, File outputFile) throws CommandLineModuleExecutionException
    {
        Map<String,Object>[] rows = null;
        try
        {
            loader = new TabLoader(massesFile);
            rows = (Map<String,Object>[]) loader.load();
            for (int i=0; i<rows.length; i++)
            {
                if (rows[i].get("mass") == null)
                    throw new CommandLineModuleExecutionException("ERROR! Missing mass values in masses file, row " + (i+1));
                if (shouldFixFeatureMasses) {
                    //fix masses
                    double oldMass = Double.parseDouble(rows[i].get("mass").toString());
                    if (rows[i].get("mz") == null || rows[i].get("charge") == null)
                        throw new CommandLineModuleExecutionException("ERROR! Missing mass values in masses file, row " + (i+1));
                     rows[i].put("mass", Double.parseDouble(rows[i].get("mz").toString()) *
                             Integer.parseInt(rows[i].get("charge").toString()));
//System.err.println("Change: " + oldMass + " -> " + Double.parseDouble(rows[i].get("mass").toString()));
                }
            }
            ApplicationContext.infoMessage("Loaded " + rows.length + " rows from masses file");
//for (ChemicalCompound comp : databaseCompoundsByMass) System.err.println(comp.getMass());
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException(e);
        }

        Arrays.sort(rows, new Comparator<Map<String, Object>>()
        {
            public int compare(Map<String, Object> o1, Map<String, Object> o2)
            {
                double mass1 = Double.parseDouble(o1.get("mass").toString());
                double mass2 = Double.parseDouble(o2.get("mass").toString());

                return mass1 == mass2 ? 0 : mass1 < mass2 ? -1 : 1;
            }
        });

        List<Feature> featuresForMasses = new ArrayList<Feature>();
        for (Map<String, Object> row : rows)
        {
            float mass = (float) Double.parseDouble(row.get("mass").toString());
            Feature feature = new Feature(0, mass, 200);
            feature.setMass(mass);
            feature.setMz(mass);
            feature.setCharge(1);
            featuresForMasses.add(feature);
        }

        ApplicationContext.infoMessage("Matching " + rows.length + " masses...");
        Map<Feature, Map<ChemicalFormula, List<Adduct>>> matchingResult =
                metDBMatcher.massMatchFull(featuresForMasses.toArray(new Feature[featuresForMasses.size()]), massTolerancePPM, 1);
        ApplicationContext.infoMessage("Done matching");

        StringBuffer headerLineBuf = new StringBuffer("id\t");
        try
        {
            PrintWriter outPW = new PrintWriter(outputFile);
            for (int i=0; i<loader.getColumns().length; i++)
            {
                TabLoader.ColumnDescriptor column = loader.getColumns()[i];
                if (i>0)
                    headerLineBuf.append("\t");
                headerLineBuf.append(column.name);
            }
            headerLineBuf.append("\tformula\tcompound\tiontype\tdeltamass\tSMILES");
            outPW.println(headerLineBuf);
            outPW.flush();
            int matchId=0;
            for (int i=0; i<rows.length; i++)
            {
                Feature feature = featuresForMasses.get(i);
                Map<ChemicalFormula, List<Adduct>> featureMatchingResult = matchingResult.get(feature);
                if (featureMatchingResult == null)
                    continue;
                matchId++;

                int matchIdForOutput = matchId;
                if (rows[i].containsKey("id"))
                    matchIdForOutput = (Integer) rows[i].get("id");
                StringBuffer lineBufOrigCols = new StringBuffer("" + matchIdForOutput);
                Map<String, Object> row = rows[i];
                for (int j=0; j<loader.getColumns().length; j++)
                {
                    TabLoader.ColumnDescriptor column = loader.getColumns()[j];
                    lineBufOrigCols.append("\t");
                    Object val = row.get(column.name);
                    if (val != null)
                        lineBufOrigCols.append(val);
                }

                for (ChemicalFormula formula : featureMatchingResult.keySet())
                {
                    for (Adduct adduct : featureMatchingResult.get(formula))
                    {
                        outPW.println(lineBufOrigCols.toString() + "\t" + formula + "\t" +
                                      adduct.getCompound().getName() + "\t" + adduct.getIonTypeString() + "\t" +
                                Rounder.round(feature.getMass() - adduct.getCommonestIsotopeMass(), 5) + "\t" +
                                      ChemCalcs.createSMILESString(adduct.getMolecule()));
                        outPW.flush();
                    }
                }   
            }
            ApplicationContext.infoMessage("Matched " + matchId + " out of " + rows.length + " rows");
            outPW.close();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed writing output file",e);
        }
    }
}


