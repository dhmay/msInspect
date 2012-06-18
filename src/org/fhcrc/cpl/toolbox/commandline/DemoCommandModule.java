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
package org.fhcrc.cpl.toolbox.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

/**
 * Demo CommandLineModule.  This module accepts one argument, the path to a PepXML file.  It reads every
 * search_result in the PepXML file (across all fractions, if there are multiple) and histograms the lengths
 * of all the peptide sequences.  PepXML parsing depends on PepXmlLoader, which in turn depends on JRAP, from
 * the Institute for Systems Biology.
 *
 * This Demo module implements a main() method so that it can be invoked as a Java application,  The main()
 * method calls a method in CommandLineModuleUtilities that parses the raw arguments from the command line and
 * places them in a Map that is processed by the CommandLineModule.
 *
 * This module is extremely simple.  More complicated modules can be found in the package
 * org.fhcrc.cpl.viewer.commandline.modules. The module FindPeptidesCommandLineModule is well-documented for
 * tutorial purposes.
 */
public class DemoCommandModule extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    //This is explicitly needed by PepXmlLoader, but it is always good practice for classes to define a
    //Logger that can be used for multiple levels of log messages
    protected static Logger _log = Logger.getLogger(DemoCommandModule.class);

    //Input file specified by the user
    protected File inFile;

    /**
     * No-arg constructor just calls the initializer
     */
    public DemoCommandModule()
    {
        init();
    }

    /**
     * Initialize the Command Module
     */
    protected void init()
    {
        //Name of the command
        mCommandName = "demo";

        //A help message for auto-documentation purposes
        mHelpMessage =
                "Demonstration module.  Reads every " +
                " search_result in the pepXML file (across all fractions, if there are multiple) " +
                "and histograms the lengths of all the peptide sequences.";

        //A short description of the command, for auto-documentation purposes
        mShortDescription = "Perform mass calibration on mzXML spectra";

        //Only one argument is defined.  This "unnamed" argument only requires a value, not a name.
        //I.e., the user simply enters the filepath, not "--file=<filepath>".  Only one unnamed argument
        //is allowed in a Command Module, for obvious reasons, but an unnamed series of arguments can be used instead
        CommandLineArgumentDefinition[] argDefs =
            {
                    createUnnamedFileArgumentDefinition(true, "Input PepXML file"),
            };
        //Add all the argument definitions we just created as basic arguments
        addArgumentDefinitions(argDefs);
    }

    /**
     * Argument processing -- in this case, just grab the specified File and assign it to a member variable
     * @throws ArgumentValidationException
     */
    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFile = this.getUnnamedFileArgumentValue();
    }

    /**
     * Called after arguments are processed.  Does the actual work
     * @throws CommandLineModuleExecutionException
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        //PepXmlLoader is a wrapper on top of JRAP.  Open the PepXML file with PepXmlLoader
        PepXmlLoader loader;
        try
        {
            loader = new PepXmlLoader(inFile, _log);
        }
        catch (FileNotFoundException e)
        {
            throw new CommandLineModuleExecutionException("File " + inFile.getAbsolutePath() + " not found");
        }
        catch (XMLStreamException e2)
        {
            throw new CommandLineModuleExecutionException("Failed to parse file " + inFile.getAbsolutePath(), e2);
        }

        //Maintain a list of peptide sequence lengths
        List<Float> sequenceLengths = new ArrayList<Float>();

        //Loop through all fractions in the PepXML file
        PepXmlLoader.FractionIterator fractionIterator = loader.getFractionIterator();
        while (fractionIterator.hasNext())
        {
            PepXmlLoader.PepXmlFraction fraction = fractionIterator.next();
            //Loop on all Peptides (search_results) in the fraction
            PepXmlLoader.PeptideIterator peptideIterator = fraction.getPeptideIterator();
            while (peptideIterator.hasNext())
            {
                PepXmlLoader.PepXmlPeptide peptide = peptideIterator.next();
                //Get the peptide's sequence length and add it to our list
                sequenceLengths.add((float) peptide.getPeptide().length());
            }
        }

        //Create a histogram of the sequence lengths
        PanelWithHistogram pwh = new PanelWithHistogram(sequenceLengths, "Peptide Sequence Lengths");
        //Display the histogram in a tab.  There's also a displayInDialog() method, if you know your module won't
        //ever need to display more than one chart.  Most modules tend to end up displaying multiple charts, so
        //the tab display comes in very handy
        pwh.displayInTab();
    }

    /**
     * This method instantiates our demo module, parses the arguments, and runs the module.  This main method
     * could be copied verbatim into a new class.
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception
    {
        //Instantiate the module and parse the arguments
        CommandLineModule module = new DemoCommandModule();
        try
        {
            Map<String,String> argNameValueMap = CommandLineModuleUtilities.parseRawArguments(module, args);
            module.digestArguments(argNameValueMap);
        }
        catch (ArgumentValidationException e)
        {
            ApplicationContext.infoMessage("Failure while parsing arguments:");
            if (e.shouldShowStackTrace())
                ApplicationContext.errorMessage(e.getMessage(), e);
            else
                ApplicationContext.infoMessage(e.getMessage());
            return;
        }
        catch (IllegalArgumentException iae)
        {
            ApplicationContext.infoMessage(iae.getMessage());
            return;
        }

        //Execute the module
        try
        {
            ApplicationContext.infoMessage("Reading file...");
            module.execute();
            ApplicationContext.infoMessage("Done.");            
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("ERROR: Failed to histogram sequence lengths",e);
        }
    }

}
