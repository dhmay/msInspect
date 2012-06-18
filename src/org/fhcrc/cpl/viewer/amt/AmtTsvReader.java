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
package org.fhcrc.cpl.viewer.amt;

import java.io.*;
import java.util.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.text.ParseException;

import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.statistics.MatrixUtil;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.gui.chart.ScatterPlotDialog;

/**
 * A class for writing AMT TSV files.
 */
public class AmtTsvReader
{
    static Logger _log = Logger.getLogger(AmtTsvReader.class);

    protected static DateFormat dateParser =
            new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");

    /**
     * Read an AMT database from a tsv file
     * @param file
     * @throws IOException
     */
    public static AmtDatabase read(File file) throws IOException
    {
        TabLoader tabLoader = new TabLoader(file, AmtTsvRow.class);
        AmtTsvRow[] tsvRows = (AmtTsvRow[]) tabLoader.load();

        Map<Integer, List<AmtTsvRow>> runNumberTsvRowMap =
                new HashMap<Integer, List<AmtTsvRow>>();

        //should we use time values, or scan numbers?  Times are preferable,
        //if available
        boolean useTimesInsteadOfScans = true;

        //figure out the hydrophobicity algorithm and version that was used
        String hydroAlgorithmString =
                (String) tabLoader.getComments().get("HydrophobicityAlgorithm");
        String[] hydroAlgPieces = hydroAlgorithmString.split(",");
        String hydroAlgorithmName = hydroAlgPieces[0].toLowerCase();
        double hydroAlgorithmVersion = Double.parseDouble(hydroAlgPieces[1]);

        //check that the specified hydrophobicity algorithm and version
        //are known to us.  If not, can't normalize, die.
        if (!HydrophobicityNormalizer.algorithmVersionKnown(hydroAlgorithmName,
                                                            hydroAlgorithmVersion))
            throw new RuntimeException("Unknown hydrophobicity algorithm and version: " + hydroAlgorithmName + ", version " + hydroAlgorithmVersion);

        //first get all the hydrophobicity values and times/scans for all the rows
        for (AmtTsvRow tsvRow : tsvRows)
        {
            //normalize hydrophobicity
            tsvRow.setH(HydrophobicityNormalizer.normalize(tsvRow.getH(),
                    hydroAlgorithmName,hydroAlgorithmVersion));

            //if any row doesn't have retention time, use scan instead
            //TODO:  for some reason, tabloader sets time to 0 if it's not specified.
            //TODO: that isn't all that cool.  Rows with time 0 will throw everything off
            if (tsvRow.getRetentiontime() <= 0)
                useTimesInsteadOfScans = false;

            //strip tryptic ends from all rows
            tsvRow.setPeptide(stripTrypticEnds(tsvRow.getPeptide()));

            //Associate the run numbers with lists of rows
            if (!runNumberTsvRowMap.containsKey(tsvRow.getRun()))
            {
                runNumberTsvRowMap.put(tsvRow.getRun(), new ArrayList<AmtTsvRow>());
            }
            runNumberTsvRowMap.get(tsvRow.getRun()).add(tsvRow);
        }
//System.err.println("***read 2, rows:" + tsvRows.length);

        //map run numbers to runs, so we can get at them easily later
        Map<Integer, AmtRunEntry> runNumberRunMap = new HashMap<Integer, AmtRunEntry>();

        //map variable mod symbols to their associated modifications, within each run
        Map<AmtRunEntry, Map<String, MS2Modification>> runSymbolModMap =
            new HashMap<AmtRunEntry, Map<String, MS2Modification>>();

        //create the actual run entries
        for (int runNumber : runNumberTsvRowMap.keySet())
        {
            AmtRunEntry runEntry = createRunEntry(runNumberTsvRowMap.get(runNumber),
                                                  (String) tabLoader.getComments().get("Run" + runNumber),
                                                  useTimesInsteadOfScans);
            runNumberRunMap.put(runNumber, runEntry);
        }
//System.err.println("***read 3, runs: " + runSymbolModMap.keySet().size());

        AmtDatabase result = new AmtDatabase();

        result.setHydrophobicityAlgorithmName(hydroAlgorithmName);
        result.setHydrophobicityAlgorithmVersion(hydroAlgorithmVersion);

        //addRunEntry implicitly overrides duplicate modifications found in the run
        for (AmtRunEntry runEntry : runNumberRunMap.values())
            result.addRunEntry(runEntry);

        //map the symbols of modifications to their modifications
        for (AmtRunEntry runEntry : runNumberRunMap.values())
        {
            Map<String, MS2Modification> thisRunSymbolModMap =
                    new HashMap<String, MS2Modification>();
            for (MS2Modification mod : runEntry.getModifications())
            {
                if (mod.getVariable())
                    thisRunSymbolModMap.put(mod.getSymbol(), mod);
//System.err.println("   adding mod for char " + mod.getSymbol());
            }
            runSymbolModMap.put(runEntry, thisRunSymbolModMap);
        }

        //add the individual observations
        for (AmtTsvRow tsvRow : tsvRows)
        {
            String strippedPeptide = stripPeptide(tsvRow.getPeptide());
            AmtRunEntry run = runNumberRunMap.get(tsvRow.getRun());
            List<MS2Modification>[] modifications =
                    identifyAllModifications(strippedPeptide, tsvRow.getPeptide(), run,
                                             runSymbolModMap.get(run));
            result.addObservation(strippedPeptide, modifications, 0,
                   run.convertTimeToHydrophobicity(useTimesInsteadOfScans ?
                           tsvRow.getRetentiontime() : tsvRow.getScan()), run,
                           AmtPeptideEntry.AmtPeptideObservation.SPECTRAL_COUNT_UNKNOWN,
                    useTimesInsteadOfScans ?
                           tsvRow.getRetentiontime() : tsvRow.getScan());
        }
//testing code
double[] origXValues = new double[tsvRows.length];
double[] origYValues = new double[tsvRows.length];
for (int i=0; i<tsvRows.length; i++)
{
    origXValues[i] = tsvRows[i].getH();
    origYValues[i] = tsvRows[i].getScan();
}

List<Double> dbXValuesList = new ArrayList<Double>();
List<Double> dbYValuesList = new ArrayList<Double>();
for (AmtPeptideEntry entry : result.getEntries())
    for (AmtPeptideEntry.AmtPeptideObservation obs : entry.getObservations())
    {
        dbXValuesList.add(entry.getPredictedHydrophobicity());
        dbYValuesList.add(obs.getTimeInRun());
//        dbYValuesList.add(obs.getRunEntry().recoverTimeForHydrophobicity(obs.getObservedHydrophobicity()));

    }
double[] dbXValues = new double[dbXValuesList.size()];
double[] dbYValues = new double[dbXValuesList.size()];
for (int i=0; i<dbXValues.length; i++)
{
    dbXValues[i] = dbXValuesList.get(i);
    dbYValues[i] = dbYValuesList.get(i);
}

ScatterPlotDialog chartDialog = new ScatterPlotDialog(origXValues, origYValues, "from tsv");
chartDialog.addData(dbXValues, dbYValues, "from db");
chartDialog.setVisible(true);
        return result;
    }

    /**
     * Strip tryptic ends off a raw peptide string
     * @param peptideString
     * @return
     */
    protected static String stripTrypticEnds(String peptideString)
    {
        if (peptideString.charAt(1) == '.' &&
            peptideString.charAt(peptideString.length()-2) == '.')
        {
            peptideString = peptideString.substring(2, peptideString.length()-2);
        }
        return peptideString;
    }

    /**
     * Return a version of the peptide string that contains only the residues
     * @param peptideWithMods
     * @return
     */
    protected static String stripPeptide(String peptideWithMods)
    {
        StringBuffer resultBuf = new StringBuffer();

        for (byte charByte : peptideWithMods.getBytes())
        {
            char thisChar = (char) charByte;
            if (Character.isLetter(thisChar))
                resultBuf.append(thisChar);
        }
        return resultBuf.toString();
    }

    /**
     * Figure out what modifications exist in this peptide string
     * @param strippedPeptide
     * @param peptideWithMods
     * @param run
     * @param thisRunSymbolModMap
     * @return
     */
    protected static List<MS2Modification>[]
            identifyAllModifications(String strippedPeptide,
                                     String peptideWithMods,
                                     AmtRunEntry run,
                                     Map<String, MS2Modification> thisRunSymbolModMap)
    {
//System.err.println("identifyAllMods, pep=" + strippedPeptide + ", wmods=" + peptideWithMods);
        List<MS2Modification>[] result =
                (List<MS2Modification>[]) new List[strippedPeptide.length()];


        int peptideSequenceIndex = -1;
        for (byte charByte : peptideWithMods.getBytes())
        {
            String thisChar = Character.toString((char) charByte);
            if (Character.isLetter(thisChar.charAt(0)))
            {
                peptideSequenceIndex++;

                MS2Modification mod = run.getStaticMod(thisChar);
                if (mod != null)
                {
                    if (result[peptideSequenceIndex] == null)
                        result[peptideSequenceIndex] = new ArrayList<MS2Modification>();
                    result[peptideSequenceIndex].add(mod);
                }
            }
            else
            {
                if (result[peptideSequenceIndex] == null)
                    result[peptideSequenceIndex] = new ArrayList<MS2Modification>();
                MS2Modification mod = thisRunSymbolModMap.get(thisChar);
                if (mod == null)
                    throw new RuntimeException("AmtTsvReader:  row contained a mod with character " + thisChar + ", which is not defined at run level");
                result[peptideSequenceIndex].add(mod);
            }

        }
        return result;
    }

    /**
     * Create a run entry from a comment row in the tsv file
     * @param amtTsvRowList
     * @param commentLineForRun
     * @param useTimesInsteadOfScans
     * @return
     */
    protected static AmtRunEntry createRunEntry(List<AmtTsvRow> amtTsvRowList,
                                                String commentLineForRun,
                                                boolean useTimesInsteadOfScans)
    {
        double[] elutionValues = new double[amtTsvRowList.size()];
        double[] hydrophobicities = new double[amtTsvRowList.size()];

        for (int i=0; i<amtTsvRowList.size(); i++)
        {
            AmtTsvRow tsvRow = amtTsvRowList.get(i);
            hydrophobicities[i] = tsvRow.getH();
            elutionValues[i] = useTimesInsteadOfScans? tsvRow.getRetentiontime() :
                                                       tsvRow.getScan();
        }
//System.err.println("useTimes? " + useTimesInsteadOfScans);

        double[] timeToHydroArray =
                MatrixUtil.linearRegression(elutionValues, hydrophobicities);

//#Run92=Desktop/haloICAT2_30 (mascot)|haloICAT2_30.pep.xml|2006-07-09 08:37:56.234|C'=442.225;C"=450.275

        String[] fields = commentLineForRun.split("\\|");
        //TODO: do something with description
        String description = fields[0];
        String pepXmlFilename = fields[1];
        String dateString = fields[2];
        String modificationsString = fields[3];

        //modifications
        String[] modificationStrings = modificationsString.split(";");
        MS2Modification[] modifications = null;
        if (modificationStrings != null && modificationStrings.length>0 &&
            !("".equals(modificationStrings[0])))
        {
            modifications = new MS2Modification[modificationStrings.length];
            for (int i=0; i<modificationStrings.length; i++)
            {
                String modString = modificationStrings[i];

                String residueWithPossibleMod =
                        modString.substring(0, modString.indexOf('='));
                String residue = residueWithPossibleMod.substring(0,1);
//System.err.println("mod****  " + modString + "***" + residueWithPossibleMod + "***" + residue);
                boolean variable = false;
                String symbol = null;
                if (residueWithPossibleMod.length() > 1)
                {
                    variable = true;
                    symbol = residueWithPossibleMod.substring(1,2);
                }
                double mass = Double.parseDouble(modString.substring(modString.indexOf('=') + 1));

                MS2Modification mod = new MS2Modification();
                mod.setAminoAcid(residue);
                mod.setMass((float) mass);
                mod.setMassDiff((float) (mass -
                        PeptideGenerator.AMINO_ACID_MONOISOTOPIC_MASSES[residue.charAt(0)]));
                mod.setVariable(variable);
                mod.setSymbol(symbol);
//System.err.println("mod: " + mod.getAminoAcid() + ", " + mod.getMassDiff() + ", var?"+mod.getVariable() + ", sym="+mod.getSymbol());
                modifications[i] = mod;
            }
        }

        Date date = new Date();
        try
        {
            dateParser.parse(dateString);
        }
        catch (ParseException e)
        {
            _log.error("error parsing date format (" + dateString +
                       ") for TSV run, using current date");
        }
        //slope, int, slope, int
        AmtRunEntry run =new AmtRunEntry(timeToHydroArray, modifications, date);
        run.setPepXmlFilename(pepXmlFilename);
        return run;
    }

    /**
     * Used by TabLoader to represent one row.  Rows have the format:
     * run     fraction        scan    retentiontime   peptide h
     */
    protected static class AmtTsvRow
    {
        protected int run;
        protected int fraction;
        protected int scan;
        protected double retentiontime;
        protected String peptide;
        protected double h;

        public AmtTsvRow()
        {
            run=-1; fraction=-1; scan=-1; retentiontime=-1;
            peptide="";
            h=-1;
        }


        public int getRun()
        {
            return run;
        }

        public void setRun(int run)
        {
            this.run = run;
        }

        public int getFraction()
        {
            return fraction;
        }

        public void setFraction(int fraction)
        {
            this.fraction = fraction;
        }

        public int getScan()
        {
            return scan;
        }

        public void setScan(int scan)
        {
            this.scan = scan;
        }

        public double getRetentiontime()
        {
            return retentiontime;
        }

        public void setRetentiontime(double retentiontime)
        {
            this.retentiontime = retentiontime;
        }

        public String getPeptide()
        {
            return peptide;
        }

        public void setPeptide(String peptide)
        {
            this.peptide = peptide;
        }

        public double getH()
        {
            return h;
        }

        public void setH(double h)
        {
            this.h = h;
        }
    }
}
