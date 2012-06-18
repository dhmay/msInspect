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
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * Module for matching the masses in a 'peptide' array to a metabolite database. Uses the mass in the middle
 * of each cell.
 */
public class MatchArrayMetMassesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MatchArrayMetMassesCLM.class);

    protected File arrayFile;
    protected MetaboliteDatabaseMatcher metDBMatcher;
    protected File outFile;

    protected float massTolerancePPM = 2f;

    protected List<ChemicalCompound> databaseCompoundsByMass = null;

    protected boolean shouldFixArrayMasses = true;

    protected boolean shouldUseBaseMod = false;

    protected String[] annotationColumnNames =
            new String[] { "class","Accession_code","cas_number","kegg_id","subclass","species","pathway1","pathway2"};

    Map<String,Map<String,Object>> annotations = null;


    public MatchArrayMetMassesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "matcharrmetmasses";
        mShortDescription = "Match array metabolite masses to a metabolite database";
        mHelpMessage = "Metabolite database should be a .tsv file with (at least) the following columns " +
                "(all but the first 3 of which may be empty): " + 
          "name, formula, smiles, class, Accession_code, cas_number, kegg_id, subclass, species, pathway1, pathway2";
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Input array file"),
                        new FileToReadArgumentDefinition("metdb", true, "Metabolite database file"),
                        new FileToWriteArgumentDefinition("out", false, "out"),
                        new DirectoryToWriteArgumentDefinition("outdir", false, "outdir"),
                        new BooleanArgumentDefinition("fixarraymasses", false, "Should fix array masses so that " +
                                "they are mz * z? (this is a 'correction' for arrays whose masses represent " +
                                "M+H ion masses)", shouldFixArrayMasses),
                        new DecimalArgumentDefinition("deltappm", false,
                                "Delta mass (ppm) between the center of the array cell and the mass of the metabolite",
                                massTolerancePPM),
                        new BooleanArgumentDefinition("usebasemod", false, "Include base mass, i.e., [M]?",
                                shouldUseBaseMod),
                        new StringListArgumentDefinition("mods", true, "modifications to use on every compound"),
                        new BooleanArgumentDefinition("showcharts", false, "show charts?", false),
                };
        addArgumentDefinitions(argDefs);
    }

    public static ChemicalModification createModForString(String modString) {
        if (modString.equals("M+H"))
            return new UnknownHAdditionMod();
        if (modString.equals("M+Na"))
            return new UnknownFormulaAdditionMod(new ChemicalFormula("Na1"), "M+Na","UnknownNaAdditionMod");
        if (modString.equals("2M+H"))
            return new DoubleFormulaPlusSomethingMod(new ChemicalFormula("H"), "2M+H","DoubleMassPlusHMod");
        if (modString.equals("2M+Na"))
            return new DoubleFormulaPlusSomethingMod(new ChemicalFormula("Na"), "2M+Na","DoubleMassPlusNaMod");
        return null;
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        arrayFile = this.getUnnamedFileArgumentValue();

        shouldFixArrayMasses = getBooleanArgumentValue("fixarraymasses");
        if (shouldFixArrayMasses)
            System.err.println("'Fixing' feature masses by making them 1Da higher than would be the case in proteomics.");
        shouldUseBaseMod = getBooleanArgumentValue("usebasemod");

        metDBMatcher = new MetaboliteDatabaseMatcher();
        List<ChemicalModification> mods = new ArrayList<ChemicalModification>();


        for (String modString : getStringListArgumentValue("mods")) {
            mods.add(createModForString(modString));
        }

        metDBMatcher.setChemicalModifications(mods);
        try
        {
            databaseCompoundsByMass = ChemicalCompound.loadCompoundsFromFile(getFileArgumentValue("metdb"), 2);

            Collections.sort(databaseCompoundsByMass, new ChemicalCompound.ComparatorMassAsc());
            annotations = loadAnnotations(getFileArgumentValue("metdb"));
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to load metabolite database",e);
        }
        metDBMatcher.setDatabaseCompounds(databaseCompoundsByMass);
        metDBMatcher.setUseUnmodifiedAdduct(shouldUseBaseMod);




        massTolerancePPM = getFloatArgumentValue("deltappm");

        outFile = getFileArgumentValue("out");

        List<ChemicalModification> chemicalModsToAdd = new ArrayList<ChemicalModification>();

        for (String modString : getStringListArgumentValue("mods")) {
            ChemicalModification mod = createModForString(modString);
            if (mod == null)
                throw new ArgumentValidationException("Failed to process modification string: " + modString);
            chemicalModsToAdd.add(mod);
            ApplicationContext.infoMessage("Adding mod: " + mod.getSymbol());
        }
        metDBMatcher.setChemicalModifications(chemicalModsToAdd);
        metDBMatcher.setShowCharts(getBooleanArgumentValue("showcharts"));
    }

    public void execute() throws CommandLineModuleExecutionException
    {
//        PeptideArrayAnalyzer arrayAnalyzer = null;
        Map<String, Object>[] arrayRows;
        List<Feature> rowFeatures = new ArrayList<Feature>();
        TabLoader loader ;
        String headerLine;
        try {
//            arrayAnalyzer = new PeptideArrayAnalyzer(arrayFile);
//            Map<Integer, Map<String, List<Feature>>> detailsRowMap = arrayAnalyzer.loadDetailsRowMapsById();
             loader = new TabLoader(arrayFile);

            StringBuffer headerLineBuf = new StringBuffer("");
            try {
                for (TabLoader.ColumnDescriptor column : loader.getColumns())
                    headerLineBuf.append(column.name + "\t");
                headerLineBuf.append("formula\tmatchcount\tcompound\tppm\thmdb_id\tcas_number\tkegg_id\tiontype\tdeltamass\tSMILES\tclass\tsubclass\tspecies\tpathway1\tpathway2");
            } catch (IOException e) {
                throw new CommandLineModuleExecutionException(e);
            }
            headerLine = headerLineBuf.toString();


            arrayRows = (Map<String, Object>[]) loader.load();
            for (Map<String, Object> row : arrayRows) {
                int id = Integer.parseInt(row.get("id").toString());
//                Map<String, List<Feature>> runFeaturesMap = detailsRowMap.get(id);
                int charge = (Integer) row.get("charge");
                double minMass = (Double) row.get("minMass");
                double maxMass = (Double) row.get("maxMass");
                if (shouldFixArrayMasses) {
                    minMass += (charge * Spectrum.HYDROGEN_ION_MASS);
                    maxMass += (charge * Spectrum.HYDROGEN_ION_MASS);
                    row.put("minMass", minMass);
                    row.put("maxMass", minMass);
                }
                double meanMass = (minMass + maxMass) / 2;

                Feature feature = new Feature(0, (float) meanMass, 200);
                feature.setMass((float) meanMass);
                feature.setMz((float) meanMass / charge);
                feature.setCharge(charge);
                rowFeatures.add(feature);
            }
        }
        catch (IOException e) {
            throw new CommandLineModuleExecutionException(e);
        }


        ApplicationContext.infoMessage("Matching " + rowFeatures.size() + " masses...");
        Map<Feature, Map<ChemicalFormula, List<Adduct>>> matchingResult =
                metDBMatcher.massMatchFull(rowFeatures.toArray(new Feature[rowFeatures.size()]), massTolerancePPM, 1);
        ApplicationContext.infoMessage("Done matching");

        try
        {
            PrintWriter outPW = new PrintWriter(outFile);
            outPW.println(headerLine);
            outPW.flush();

            int matchCount = 0;
            ApplicationContext.infoMessage("Writing output array...");
            for (int i=0; i<arrayRows.length; i++)
            {
                int matchCountThisRow = 0;
                StringBuffer lineBuf = new StringBuffer();
                for (TabLoader.ColumnDescriptor column : loader.getColumns()) {
                    Object val = arrayRows[i].get(column.name);
                    if (val != null)
                        lineBuf.append(val.toString());
                    lineBuf.append( "\t");
                }
                Feature feature = rowFeatures.get(i);
                Map<ChemicalFormula, List<Adduct>> featureMatchingResult = matchingResult.get(feature);
                                    StringBuffer formulasBuf = new StringBuffer();
                StringBuffer namesBuf = new StringBuffer();
                StringBuffer ionTypesBuf = new StringBuffer();
                StringBuffer massDiffsBuf = new StringBuffer();
                StringBuffer smilesBuf = new StringBuffer();
                StringBuffer classBuf = new StringBuffer();
                StringBuffer subclassBuf = new StringBuffer();
                StringBuffer speciesBuf = new StringBuffer();
                StringBuffer pathway1Buf = new StringBuffer();
                StringBuffer pathway2Buf = new StringBuffer();
                StringBuffer hmdbIdBuf = new StringBuffer();
                StringBuffer casNumberBuf = new StringBuffer();
                StringBuffer ppmBuf = new StringBuffer();
                StringBuffer keggIdBuf = new StringBuffer();



                List<StringBuffer> allBufs = Arrays.asList(
                        new StringBuffer[] {formulasBuf, namesBuf, ppmBuf, hmdbIdBuf, casNumberBuf, keggIdBuf,
                                ionTypesBuf,massDiffsBuf,smilesBuf,classBuf,subclassBuf,
                                speciesBuf,pathway1Buf,pathway2Buf});


                if (featureMatchingResult != null) {
                    matchCount++;
                    boolean first = true;
                    for (ChemicalFormula formula : featureMatchingResult.keySet())
                    {
                        for (Adduct adduct : featureMatchingResult.get(formula))
                        {
                            matchCountThisRow++;
                            if (!first) {
                                for (StringBuffer buf : allBufs)
                                    buf.append(";");
                            }
                            formulasBuf.append(formula);
                            namesBuf.append(adduct.getCompound().getName());
                            ionTypesBuf.append(adduct.getIonTypeString());
                            massDiffsBuf.append(Rounder.round(feature.getMass() - adduct.getCommonestIsotopeMass(), 5));
                            if (adduct.getMolecule() != null)
                                smilesBuf.append(ChemCalcs.createSMILESString(adduct.getMolecule()));
                            String name = adduct.getCompound().getName();
                            classBuf.append(annotations.get("class").get(name));
                            subclassBuf.append(annotations.get("subclass").get(name));
                            speciesBuf.append(annotations.get("species").get(name));
                            pathway1Buf.append(annotations.get("pathway1").get(name));
                            pathway2Buf.append(annotations.get("pathway2").get(name));
                            hmdbIdBuf.append(annotations.get("Accession_code").get(name));
                            casNumberBuf.append(annotations.get("cas_number").get(name));
                            keggIdBuf.append(annotations.get("kegg_id").get(name));
                            float ppmDiff = (float) MassUtilities.convertDaToPPM(
                                    feature.getMz() - (float) adduct.getCommonestIsotopeMass(),
                                    (float) adduct.getCommonestIsotopeMass());
                            ppmBuf.append(Rounder.round(ppmDiff,1));
                            first = false;
                        }
                    }
                }
                lineBuf.append("\t" + matchCountThisRow);
                boolean firstBuf = true;
                for (StringBuffer buf : allBufs) {
                    if (!firstBuf)
                        lineBuf.append("\t" + buf.toString());
                    firstBuf = false;
                }
                outPW.println(lineBuf.toString());
                outPW.flush();
            }
            ApplicationContext.infoMessage("Matched " + matchCount + " out of " + arrayRows.length + " array rows");
            outPW.close();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed writing output file",e);
        }
    }



    protected Map<String,Map<String,Object>> loadAnnotations(File file)
            throws IOException
    {

        ApplicationContext.infoMessage("Loading annotations from file " + file.getAbsolutePath());
        TabLoader loader = new TabLoader(new FileReader(file),true);

ApplicationContext.infoMessage("Annotation file columns:");
for (TabLoader.ColumnDescriptor col : loader.getColumns()) ApplicationContext.infoMessage("\t" + col.name);
        //a separate map for each annotation column
        Map<String,Map<String,Object>> result = new HashMap<String,Map<String,Object>>();
//System.err.println("Loading annotations...");
//Set<String> allElements = new HashSet<String>();
        for (Map row : (Map[]) loader.load())
        {
//try{
//allElements.addAll(ChemCalcs.emp2atomCount((String)row.get("formula")).keySet());
//}catch(Exception e){}

//System.err.println("ROW");
//for (Object key1 : row.keySet()) System.err.println("\t" + key1.toString());
            Object key = row.get("name");
            if (key != null)
            {
//System.err.println("CHECK");
                String keyString = key.toString();

                    for (String annotationColumnName : annotationColumnNames)
                    {
//System.err.println("Checking " + annotationColumnName + "...");
                        if (row.get(annotationColumnName) == null)
                            continue;
                        Map<String, Object> mapThisAnnotColumn = result.get(annotationColumnName);
                        if (mapThisAnnotColumn == null)
                        {
                            mapThisAnnotColumn = new HashMap<String, Object>();
                            result.put(annotationColumnName, mapThisAnnotColumn);
                        }
                        mapThisAnnotColumn.put(keyString,row.get(annotationColumnName));
//System.err.println("\t found " +row.get(annotationColumnName) );
                    }

            }
        }
        ApplicationContext.infoMessage("Loaded annotation values for " + result.size() + " columns:");
        for (String annot : result.keySet())
        {
            ApplicationContext.infoMessage("\t" + annot + ": " + result.get(annot).size());
        }
//for (String element : allElements) System.err.println(element);
        return result;
    }

}


