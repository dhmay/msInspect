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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FastaFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.Protein;


import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Creates a forward-reverse fasta file from an input fasta file
 */
public class RestrictFastaProteinsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(RestrictFastaProteinsCLM.class);


    protected File outFile = null;

    protected Protein[] fastaProteins = null;

    protected boolean keepProteinsOnList = true;

    protected List<String> proteinIDList = new ArrayList<String>();


    public RestrictFastaProteinsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "restrictfasta";
        mShortDescription = "Remove proteins not on (or on) a list from a fasta file";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        new FastaFileArgumentDefinition(
                                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,true, null),
                        new FileToWriteArgumentDefinition("out", true, null),
                        new FileToReadArgumentDefinition("proteinfile", true, "Protein list file, one line per protein ID"),
                        new BooleanArgumentDefinition("keeponlist", false, "Keep the proteins on the list? (if false, strip them)", keepProteinsOnList),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        fastaProteins = (Protein[]) this.getUnnamedArgumentValue();
        outFile = getFileArgumentValue("out");
        keepProteinsOnList = getBooleanArgumentValue("keeponlist");
        File proteinFile = getFileArgumentValue("proteinfile");

        try {
            BufferedReader br = new BufferedReader(new FileReader(proteinFile));
            String line = null;
            while((line = br.readLine()) != null) {
                proteinIDList.add(line);
            }
            if (proteinIDList.isEmpty())
                throw new ArgumentValidationException("No proteins in file.");
            ApplicationContext.infoMessage("Loaded " + proteinIDList.size() + " protein IDs from file.");
        } catch (Exception e) {
            throw new ArgumentValidationException("Failed to process proteinfile",e);
        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        List<Protein> proteinsToKeep = new ArrayList<Protein>();
        Set<String> idsLeft = new HashSet<String>(proteinIDList);

        try
        {
            PrintWriter outPW = new PrintWriter(outFile);
            int count = 0;
            for (Protein protein : fastaProteins)
            {
                //This if handles both cases we want -- keeping and it's there, tossing and it's not
                if (proteinIDList.contains(protein.getLookup()) == keepProteinsOnList) {
                    printProtein(protein, outPW);
                    count++;
                    if (keepProteinsOnList)
                        idsLeft.remove(protein.getLookup());
                }
            }
            ApplicationContext.infoMessage("Wrote " + count + " proteins.");
            if (keepProteinsOnList && !idsLeft.isEmpty()) {
                ApplicationContext.infoMessage("Specified IDs not encountered: ");
                for (String protein : idsLeft) System.err.println("\t" + protein);
            }
            outPW.flush();
            outPW.close();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

    /**
     * Print the protein in FASTA format, reversing sequence and adjusting header if reverse is specified
     * @param protein
     * @param outPW
     */
    protected void printProtein(Protein protein, PrintWriter outPW)
    {
        outPW.print(">");
        outPW.println(protein.getHeader());
        String forwardSequence = protein.getSequenceAsString();
        for (int i=0; i<forwardSequence.length(); i++)
        {
            int indexToPrint = i;
            outPW.print(forwardSequence.charAt(indexToPrint));
            if ((i%80 == 79 && i > 0) ||
                    i == forwardSequence.length()-1)
                outPW.println();
        }
    }

}
