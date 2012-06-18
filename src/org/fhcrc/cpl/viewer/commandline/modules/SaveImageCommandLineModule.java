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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.gui.MSImageComponent;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;


/**
 * Command linemodule for feature finding
 */
public class SaveImageCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SaveImageCommandLineModule.class);

    protected File[] mzxmlFiles = null;
    protected File outFile = null;
    protected File outDir = null;


    protected int maxWidth = Integer.MAX_VALUE;
    protected int maxHeight = Integer.MAX_VALUE;

    protected boolean includeTIC = false;

//    protected FeatureSet featureSet;

    public SaveImageCommandLineModule()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "saveimage";
        mShortDescription = "Create an image file containing a visualization of a run";
        mHelpMessage = "This command loads an mzXml file and writes out an image in which " +
                "the X axis represents time, the Y axis represents M/Z, and the shade of each " +
                "pixel represents intensity.";

        CommandLineArgumentDefinition[] argDefs =
               {
                    createUnnamedSeriesFileArgumentDefinition(
                            true, "Input mzXml file(s)"),
                    new FileToWriteArgumentDefinition("out", false, "Output image file"),
                    new DirectoryToWriteArgumentDefinition("outdir", false, "Output image directory"),

                    new IntegerArgumentDefinition("maxwidth", false, "Maximum width of output image", maxWidth),
                    new IntegerArgumentDefinition("maxheight", false, "Maximum height of output image", maxHeight),
                    new BooleanArgumentDefinition("includetic", false, "Include Total Ion Chromatogram?", includeTIC),
//                       new FeatureFileArgumentDefinition("features", false, "Features to display"),

               };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        mzxmlFiles = getUnnamedSeriesFileArgumentValues();
        outFile = getFileArgumentValue("out");
        outDir = getFileArgumentValue("outdir");

        if (mzxmlFiles.length > 1)
            assertArgumentPresent("outdir");

        if (hasArgumentValue("outdir") && hasArgumentValue("out"))
            throw new ArgumentValidationException("Can't specify both out and outdir arguments");


        maxWidth = getIntegerArgumentValue("maxwidth");
        maxHeight = getIntegerArgumentValue("maxheight");

        includeTIC = getBooleanArgumentValue("includetic");

//        featureSet = getFeatureSetArgumentValue("features");
    }

    protected String createOutFileName(File inputFile)
    {
        String fileName = inputFile.getName();
        if (fileName.toLowerCase().endsWith(".mzxml"))
            fileName = fileName.substring(0, fileName.length()-".mzxml".length());
        fileName = fileName + ".png";
        return fileName;
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        if (mzxmlFiles.length == 1)
        {
            if (outDir == null)
                outDir = mzxmlFiles[0].getParentFile();
            handleFile(mzxmlFiles[0], outFile == null ?
                    new File(outDir, createOutFileName(mzxmlFiles[0])) : outFile);
        }
        else
        {
            for (File mzxmlFile : mzxmlFiles)
            {
                ApplicationContext.setMessage("Processing file " + mzxmlFile.getAbsolutePath());
                handleFile(mzxmlFile, new File(outDir, createOutFileName(mzxmlFile)));
            }
        }
    }

    public void handleFile(File mzxmlFile, File outputFile) throws CommandLineModuleExecutionException
    {
        try
        {
            MSRun run = MSRun.load(mzxmlFile.getPath());

            MSImageComponent comp = new MSImageComponent(run.getImage(MSImageComponent.getPrefColorScheme()));
            comp.setRun(run);

//this doesn't work.  Too bad.
//TODO: make it work
//            if (featureSet != null)
//            {
//                List listForProperty = new ArrayList<FeatureSet>(1);
//                listForProperty.add(featureSet);
//                ApplicationContext.setProperty(SharedProperties.FEATURE_RANGES,
//                                               listForProperty);
//                comp.app_propertyChange(new PropertyChangeEvent(Application.getInstance(), SharedProperties.FEATURE_RANGES, null, listForProperty));
//            }

            comp.setRun(run);

            comp.saveImage(outputFile, maxWidth, maxHeight, includeTIC);
            ApplicationContext.infoMessage("Saved image file " + outputFile.getAbsolutePath());
        }
        catch (IOException x)
        {
            ApplicationContext.errorMessage("error loading run", x);
        }
    }

}
