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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToReadArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;


/**
 * Command line module for analyzing a hydrophobicity algorithm.
 * This implementation analyzes Krokhin's v3 algorithm.  To analyze another algorithm,
 * make that algorithm available as a Java routine and call it from
 * calculateHydrophobicityWithAlgorithm 
 */
public class HydrophobicityAlgorithmAnalyzerCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(HydrophobicityAlgorithmAnalyzerCLM.class);

    protected File fastaFile;

    public HydrophobicityAlgorithmAnalyzerCLM()
    {
        init();
    }

    /**
     * This method does the actual calling of the hydrophobicity algorithm.
     *
     * To change the algorithm analyzed by this class, edit this method. 
     * @param peptideSequence
     * @return
     */
    protected double calculateHydrophobicityWithAlgorithm(String peptideSequence)
    {
        Protein fakeProtein = new Protein("fake",peptideSequence.getBytes());
        Peptide peptide = new Peptide(fakeProtein,0,fakeProtein.getBytes().length);
        return peptide.getHydrophobicity3();
    }


    protected void init()
    {
        mCommandName = "analyzehydroalg";
        mShortDescription = "Template command for analyzing a hydrophobicity algorithm";
        mHelpMessage = "This command analyzes the Krokhin peptide hydrophobicity prediction algorithm (v3).  It runs the algorithm on all tryptic peptides from a given FASTA file, and returns the mean and standard deviation.\n" +
                "It also serves as a template for analyzing other algorithms.";


        CommandLineArgumentDefinition[] argDefs =
               {
                new FileToReadArgumentDefinition("fasta",true,"fasta file containing database to digest")
               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        fastaFile = getFileArgumentValue("fasta");

    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {

        List<Protein> proteins = ProteinUtilities.loadProteinsFromFasta(fastaFile);
        Set<String> peptideStrings = new HashSet<String>();
        ApplicationContext.infoMessage("FASTA file contains " + proteins.size() + " proteins.  Digesting...");
        Map<String,Double> peptideHydrophobicityScores = new HashMap<String,Double>();
        int i=0;
        for (Protein protein : proteins)
        {
            i++;
            if (i > 0 && (i % (proteins.size()/5) == 0))
                System.err.println(""+ (i*100/proteins.size()) + " % complete");
            PeptideGenerator peptideGenerator = new PeptideGenerator();
            //use tryptic digest
            peptideGenerator.setDigest(PeptideGenerator.DIGEST_TRYPTIC);
            //allow one missed cleavage
            peptideGenerator.setMaxMissedCleavages(0);
            //set minimum residues
            peptideGenerator.setMinResidues(6 - 1);

            Peptide[] thisProteinPeptides = peptideGenerator.digestProtein(protein);
            for (Peptide peptide : thisProteinPeptides)
                peptideStrings.add(new String(peptide.getChars()));
        }
        Date startTime = new Date();
        i=0;
        int numPeptides = peptideStrings.size();
        ApplicationContext.infoMessage("Found " + numPeptides + " peptides in FASTA file");
        for (String peptideString : peptideStrings)
        {
//if (numPeptides > 25000)
//{
//    System.err.println("Hit peptide limit");
//    break;
//}
            i++;
            if (i > 0 && (i % (numPeptides/20) == 0))
                ApplicationContext.infoMessage("***** "+ (i*100/numPeptides) + " % complete *****");

            peptideHydrophobicityScores.put(peptideString,
                    calculateHydrophobicityWithAlgorithm(peptideString));
        }
        long deltaMS = (new Date().getTime() - startTime.getTime());
        ApplicationContext.infoMessage("Elapsed time: " + ( deltaMS/ 1000) + " seconds, " + deltaMS + " milliseconds");

        ApplicationContext.infoMessage("database contains " + peptideHydrophobicityScores.size() + " distinct peptides");

        double[] hydroValues = new double[peptideHydrophobicityScores.size()];
        i=0;
        for (Double hydro : peptideHydrophobicityScores.values())
            hydroValues[i++] = hydro;

        double meanHydro = BasicStatistics.mean(hydroValues);
        double stddevHydro = BasicStatistics.standardDeviation(hydroValues);

        ApplicationContext.infoMessage("Mean hydrophobicity: " + meanHydro);
        ApplicationContext.infoMessage("Hydrophobicity standard deviation: " + stddevHydro);
    }
}
