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
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;

import java.io.File;
import java.util.Set;
import java.util.HashSet;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class CompareFastasCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CompareFastasCLM.class);

    protected File fastaFile1 = null;
    protected File fastaFile2 = null;

    protected File outFile = null;

    int minPeptideLength = 1;

    public CompareFastasCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "comparefastas";
        mShortDescription = "compare fastas";
        mHelpMessage = "compare fastas";
        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("fasta1",true, "FASTA file 1"),
                        new FileToReadArgumentDefinition("fasta2",true, "FASTA file 2"),
                        new FileToWriteArgumentDefinition("out",true,"output file"),
                        new IntegerArgumentDefinition("minpeptidelength", false, "minimum peptide length to consider", minPeptideLength),

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        fastaFile1 = getFileArgumentValue("fasta1");
        fastaFile2 = getFileArgumentValue("fasta2");
        outFile = getFileArgumentValue("out");
        minPeptideLength = getIntegerArgumentValue("minpeptidelength");
            }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        Set<String> peptides1 = loadPeptidesFromFasta(fastaFile1);
        Set<String> peptides2 = loadPeptidesFromFasta(fastaFile2);

        Set<String> peptidesUnique1 = new HashSet<String>();
        Set<String> peptidesUnique2 = new HashSet<String>();
        Set<String> commonPeptides = new HashSet<String>();
        Set<String> allPeptides = new HashSet<String>();


        for (String peptide1 : peptides1)
        {

            allPeptides.add(peptide1);
            if (peptides2.contains(peptide1))
                commonPeptides.add(peptide1);
            else
                peptidesUnique1.add(peptide1);
        }
        for (String peptide2 : peptides2)
        {
            allPeptides.add(peptide2);

            if (!peptides1.contains(peptide2))
                peptidesUnique2.add(peptide2);
        }

        ApplicationContext.infoMessage("Total peptides: " + allPeptides.size());
        ApplicationContext.infoMessage("Unique peptides in fasta 1: " + peptidesUnique1.size() + " (" +
                (100 * peptidesUnique1.size() / allPeptides.size()) + "%)");
        ApplicationContext.infoMessage("Unique peptides in fasta 2: " + peptidesUnique2.size()+ " (" +
                (100 * peptidesUnique2.size() / allPeptides.size()) + "%)");
        ApplicationContext.infoMessage("Common peptides: " + commonPeptides.size()+ " (" +
                (100 * commonPeptides.size() / allPeptides.size()) + "%)");

        

    }

    protected Set<String> loadPeptidesFromFasta(File fastaFile)
    {

        Protein[] fastaProteins = ProteinUtilities.loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);
        PeptideGenerator pg = new PeptideGenerator();

        Set<String> result = new HashSet<String>();

        for (Protein protein : fastaProteins)
        {
            Peptide[] peptidesThisProtein = pg.digestProtein(protein);
            for (Peptide peptideThisProtein : peptidesThisProtein)
            {
                if (peptideThisProtein.getChars().length >= minPeptideLength)
                    result.add(new String(peptideThisProtein.getChars()));
            }
        }
        return result;
    }








}
