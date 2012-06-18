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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToWriteArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FastaFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.Protein;


import java.io.*;


/**
 * Creates a forward-reverse fasta file from an input fasta file
 */
public class ReverseFastaCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ReverseFastaCLM.class);


    protected File outFile = null;

    protected Protein[] forwardProteins = null;



    public ReverseFastaCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "reversefasta";
        mShortDescription = "";
        mHelpMessage = "append a reversed version of a fasta file to a fasta file";
        CommandLineArgumentDefinition[] argDefs =
                {
                        new FastaFileArgumentDefinition(
                                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,true, null),
                        new FileToWriteArgumentDefinition("out", true, null)
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        forwardProteins = (Protein[]) this.getUnnamedArgumentValue();
        outFile = getFileArgumentValue("out");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            PrintWriter outPW = new PrintWriter(outFile);

            for (Protein protein : forwardProteins)
            {
                printProtein(protein, false, outPW);
            }

            for (Protein protein : forwardProteins)
            {
                printProtein(protein, true, outPW);
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
     * @param reverse
     * @param outPW
     */
    protected void printProtein(Protein protein, boolean reverse, PrintWriter outPW)
    {
        outPW.print(">");
        if (reverse)
            outPW.print("rev_");
        outPW.println(protein.getHeader());
        String forwardSequence = protein.getSequenceAsString();
        for (int i=0; i<forwardSequence.length(); i++)
        {
            int indexToPrint = i;
            if (reverse)
                indexToPrint = forwardSequence.length() - i - 1;
            outPW.print(forwardSequence.charAt(indexToPrint));
            if ((i%80 == 79 && i > 0) ||
                    i == forwardSequence.length()-1)
                outPW.println();
        }
    }

}
