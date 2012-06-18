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
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.viewer.amt.AmtUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;
import java.util.List;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PickTargetedMS2CandidatesCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(PickTargetedMS2CandidatesCLM.class);

    protected File[] protXmlFiles;
    protected File protListFile;
    protected File fastaFile;
    protected File outFile;
    protected MS2Modification[] modifications = null;

    protected float minProteinProphet = 0.9f;
    protected float minPeptideProphet = 0.9f;

    protected int maxPeptidesPerProtein = 3;

    protected String residueToExclude = null;

    protected float minMZ = 400;
    protected float maxMZ = 1800;

    public static final double PROTON_MASS = 1.0072766;


    public PickTargetedMS2CandidatesCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "picktargetedms2candidates";
        mShortDescription = "picktargetedms2candidates";
        mHelpMessage = "picktargetedms2candidates";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "protxml file(s)"),
                        new FileToReadArgumentDefinition("protlistfile", true, "Protein list file"),
                        new FileToReadArgumentDefinition("fasta", true, "FASTA file"),
                        new FileToWriteArgumentDefinition("out", true, "output file"),


                        new DecimalArgumentDefinition("minproteinprophet", false, "Minimum ProteinProphet probability",
                                minProteinProphet),
                        new DecimalArgumentDefinition("minpprophet", false, "Minimum PeptideProphet probability",
                                minPeptideProphet),
                        new IntegerArgumentDefinition("maxpeptidesperprotein", false, "Maximum peptides per protein",
                                maxPeptidesPerProtein),
                        new ModificationListArgumentDefinition("modifications", true, "Modifications"),
                        new StringArgumentDefinition("excluderesidue", false, "Residue to exclude"),
                        new DecimalArgumentDefinition("minmz", false, "Minimum m/z for targets",
                                minMZ),
                        new DecimalArgumentDefinition("maxmz", false, "Maximum m/z for targets",
                                maxMZ),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        protXmlFiles = getUnnamedSeriesFileArgumentValues();
        protListFile = getFileArgumentValue("protlistfile");
        fastaFile = getFileArgumentValue("fasta");
        outFile = getFileArgumentValue("out");

        minProteinProphet = getFloatArgumentValue("minproteinprophet");
        minPeptideProphet = getFloatArgumentValue("minpprophet");

        maxPeptidesPerProtein = getIntegerArgumentValue("maxpeptidesperprotein");
        modifications = getModificationListArgumentValue("modifications");
        residueToExclude = getStringArgumentValue("excluderesidue");

        minMZ = getFloatArgumentValue("minmz");
        maxMZ = getFloatArgumentValue("maxmz");

    }


    public void execute() throws CommandLineModuleExecutionException
    {
        PrintWriter pw = null;
        FileReader fileReader = null;
        BufferedReader br = null;
        try
        {
            fileReader = new FileReader(protListFile);
            br = new BufferedReader(fileReader);

            List<String> proteinList = new ArrayList<String>();

            String line = null;
            while ((line = br.readLine()) != null)
                proteinList.add(line);

            Map<String, List<ProtXmlReader.Peptide>> proteinPeptideMap =
                    new HashMap<String, List<ProtXmlReader.Peptide>>();

            for (File protXmlFile : protXmlFiles)
                processProtXmlFile(protXmlFile, proteinList, proteinPeptideMap);

            ApplicationContext.infoMessage("Found peptide evidence for " + proteinPeptideMap.size() + " proteins");

            Set<String> allPeptides = new HashSet<String>();
            for (List<ProtXmlReader.Peptide> peptideList : proteinPeptideMap.values())
            {
                for (ProtXmlReader.Peptide peptide : peptideList)
                    allPeptides.add(peptide.getPeptideSequence());
            }

            ApplicationContext.infoMessage("Loading proteins for all relevant peptides (" +
                    allPeptides.size() + ")...");
            Map<String, Set<String>> allPeptideProteinMap =
                    ProteinUtilities.findFastaProteinsForPeptides(allPeptides, fastaFile);
            ApplicationContext.infoMessage("Done loading proteins from fasta, found proteins for " + allPeptideProteinMap.size() + " peptides");


            pw = new PrintWriter(outFile);
            pw.println("protein\tpeptide\tcharge\tmass\tmz");

            Set<String> allKeptPeptides = new HashSet<String>();

            for (String protein : proteinPeptideMap.keySet())
            {

                Map<String, Integer> peptideCounts = new HashMap<String, Integer>();
                Set<String> uniquePeptides = new HashSet<String>();
                Set<String> nonUniquePeptides = new HashSet<String>();


                for (ProtXmlReader.Peptide peptide : proteinPeptideMap.get(protein))
                {
                    String peptideSequence = peptide.getPeptideSequence();
                    Integer countThisPeptide = peptideCounts.get(peptideSequence);
                    if (countThisPeptide == null)
                    {
                        countThisPeptide = 0;
                    }
                    peptideCounts.put(peptideSequence,countThisPeptide+peptide.getInstances());

                    Set<String> allProteinsThisPeptide = allPeptideProteinMap.get(peptideSequence);
                    if  (allProteinsThisPeptide != null && allProteinsThisPeptide.size() == 1)
                        uniquePeptides.add(peptideSequence);
                    else
                        nonUniquePeptides.add(peptideSequence);
                }

                List<String> uniquePeptidesList = new ArrayList<String>(uniquePeptides);
                PeptideCountsComparatorDesc comp = new PeptideCountsComparatorDesc(peptideCounts);
                Collections.sort(uniquePeptidesList,comp);

                List<String> chosenPeptidesThisProtein = new ArrayList<String>();
                for (String peptide : uniquePeptidesList)
                {
                    chosenPeptidesThisProtein.add(peptide);
                    if (chosenPeptidesThisProtein.size() == maxPeptidesPerProtein)
                        break;
                }

                if (chosenPeptidesThisProtein.size() < maxPeptidesPerProtein)
                {
                    List<String> nonUniquePeptidesList = new ArrayList<String>(nonUniquePeptides);
                    Collections.sort(nonUniquePeptidesList,comp);

                    for (String peptide : nonUniquePeptidesList)
                    {
                        chosenPeptidesThisProtein.add(peptide);
                        if (chosenPeptidesThisProtein.size() == maxPeptidesPerProtein)
                            break;
                    }
                }

                Map<String,Set<Integer>> peptideChargeStates = new HashMap<String,Set<Integer>>();

                for (ProtXmlReader.Peptide peptide : proteinPeptideMap.get(protein))
                {
                    String peptideSequence = peptide.getPeptideSequence();
                    if (!chosenPeptidesThisProtein.contains(peptideSequence))
                        continue;
                    int charge = peptide.getCharge();
                    Set<Integer> chargeStatesThisPeptide = peptideChargeStates.get(peptideSequence);
                    if (chargeStatesThisPeptide == null)
                    {
                        chargeStatesThisPeptide = new HashSet<Integer>();
                        peptideChargeStates.put(peptideSequence, chargeStatesThisPeptide);
                    }
                    chargeStatesThisPeptide.add(charge);
                }

                int linesThisProtein = 0;
                for (String peptideSequence : peptideChargeStates.keySet())
                {
                    if (allKeptPeptides.contains(peptideSequence))
                        continue;
                    allKeptPeptides.add(peptideSequence);
                    float neutralPeptideMass = MassUtilities.calcModifiedPeptideNeutralMass(peptideSequence, modifications);
                    for (int charge : peptideChargeStates.get(peptideSequence))
                    {
                        double mz = ((double)neutralPeptideMass / (double)charge) + PROTON_MASS;
                        if (mz >= minMZ && mz <= maxMZ)
                        {
                            pw.println(protein + "\t" + peptideSequence + "\t" + charge + "\t" + neutralPeptideMass + "\t" + Rounder.round(mz,6));
                            pw.flush();
                            linesThisProtein++;
                        }
                    }
                }



                ApplicationContext.infoMessage("Protein " + protein + ", peptides: " + chosenPeptidesThisProtein.size() + ", lines: " + linesThisProtein);
            }

            ApplicationContext.infoMessage("Wrote file " + outFile.getAbsolutePath() + ", total peptides: " + allKeptPeptides.size());

        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            try
            {
                if (pw != null)
                    pw.close();
                if (br != null)
                    br.close();
                if (fileReader != null)
                    fileReader.close();
            }
            catch (Exception e) {}
        }
    }

    protected static class PeptideCountsComparatorDesc implements Comparator<String>
    {
        protected Map<String,Integer> peptideCountsMap;

        public PeptideCountsComparatorDesc(Map<String,Integer> peptideCountsMap)
        {
            this.peptideCountsMap = peptideCountsMap;
        }

        public int compare(String o1, String o2)
        {
            int counts1 = peptideCountsMap.get(o1);
            int counts2 = peptideCountsMap.get(o2);

            if (counts1 > counts2)
                return -1;
            else if (counts1 == counts2)
                return 0;
            else return 1;
        }
    }

    /**
     * do the actual work
     */
    public void processProtXmlFile(File protXmlFile, List<String> proteinList,
                                    Map<String, List<ProtXmlReader.Peptide>> cumProteinPeptideMap)
            throws CommandLineModuleExecutionException
    {
        try
        {
            Map<String, List<ProtXmlReader.Peptide>> curProteinPeptideMap =
                    ProteinUtilities.loadProteinPeptideMapFromProtXML(protXmlFile, minProteinProphet);
            int numProteinsFound = 0;
            for (String protein : proteinList)
            {
                if (curProteinPeptideMap.containsKey(protein))
                {
                    List<ProtXmlReader.Peptide> cumEntry = cumProteinPeptideMap.get(protein);
                    if (cumEntry == null)
                    {
                        cumEntry = new ArrayList<ProtXmlReader.Peptide>();
                        cumProteinPeptideMap.put(protein, cumEntry);
                    }
                    for (ProtXmlReader.Peptide peptide : curProteinPeptideMap.get(protein))
                        if ((residueToExclude == null || !peptide.getPeptideSequence().contains(residueToExclude)) &&
                                peptide.getNspAdjustedProbability() >= minPeptideProphet )
                            cumEntry.add(peptide);

                    numProteinsFound++;
                }
            }

            ApplicationContext.infoMessage("Found " + numProteinsFound + " in file " + protXmlFile);


        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
