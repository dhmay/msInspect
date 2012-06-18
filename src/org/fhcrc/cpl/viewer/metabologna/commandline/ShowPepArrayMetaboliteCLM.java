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
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.chem.Adduct;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.viewer.metabologna.MoleculeRenderer2D;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAdd2HMod;
import org.fhcrc.cpl.viewer.metabologna.ReduceDoubleBondAddWaterMod;
import org.fhcrc.cpl.viewer.align.PeptideArrayAnalyzer;
import org.fhcrc.cpl.viewer.gui.FeatureViewerFrame;
import org.apache.log4j.Logger;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.exception.CDKException;
import java.util.List;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.HashMap;

/**
 */
public class ShowPepArrayMetaboliteCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ShowPepArrayMetaboliteCLM.class);

    protected PeptideArrayAnalyzer arrayAnalyzer;
    protected int rowId;

    protected int width = 400;
    protected int height = 400;

    protected boolean showHydrogens = false;
    protected boolean showCarbons = false;

    protected File mzXmlDir;

    protected Map<String, MSRun> runNameRunMap = new HashMap<String, MSRun>();


    public ShowPepArrayMetaboliteCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "showpeparraymetabolite";
        mShortDescription = "showpeparraymetabolite";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        this.createUnnamedFileArgumentDefinition(true, "Peptide array with p-values, q-values"),
                        new IntegerArgumentDefinition("id", true, "Row ID to display"),
                        new DirectoryToReadArgumentDefinition("mzxmldir", true, "mzXML directory"),
                        new FileToReadArgumentDefinition("detailsfile", false, "Peptide array details file"),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        File pepArrayFile = this.getUnnamedFileArgumentValue();

        try
        {
            ApplicationContext.infoMessage("Loading array...");
            arrayAnalyzer = new PeptideArrayAnalyzer(pepArrayFile);
            if (hasArgumentValue("detailsfile"))
                arrayAnalyzer.setDetailsFile(getFileArgumentValue("detailsfile"));
            ApplicationContext.infoMessage("Loading details...");
            arrayAnalyzer.loadDetailsMap();
            ApplicationContext.infoMessage("Done loading array");
        }
        catch (IOException e)
        {
            throw new ArgumentValidationException("Failed to open array file",e);
        }

        mzXmlDir = getFileArgumentValue("mzxmldir");

        rowId = getIntegerArgumentValue("id");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        displayMetabolite(arrayAnalyzer.getRowMapWithID(rowId));
    }

    /**
     * @param compoundName
     * @param smilesString
     * @throws CDKException
     */
    protected void showCompound(String compoundName, String smilesString)
            throws CDKException
    {
        if (compoundName == null)
            compoundName = "dummy";
        ChemicalCompound compound = ChemicalCompound.createFromSmiles(compoundName, smilesString);
        if (showHydrogens)
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(compound.getCdkMolecule());

        ApplicationContext.infoMessage("Formula: " + compound.getFormula());
        ApplicationContext.infoMessage("Mass: " + compound.getCommonestIsotopeMass());

        System.err.println("\t" + compound.getFormula());

        MoleculeRenderer2D renderer = new MoleculeRenderer2D();
        renderer.setWidth(width);
        renderer.setHeight(height);
        renderer.setShouldShowHydrogens(showHydrogens);
        renderer.setShouldShowCarbons(showCarbons);

        Image image = renderer.renderMolecule(compound.getCdkMolecule());
        PanelWithBlindImageChart chart = new PanelWithBlindImageChart((BufferedImage) image, compoundName);
        chart.displayInTab();
    }

    protected void displayFeature(Feature feature, String fileString)
             throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Displaying feature " + feature);
        MSRun run = runNameRunMap.get(fileString);
        if (run == null)
        {
            try
            {
                File mzXmlFile = CommandLineModuleUtilities.findFileLikeFile(new File(fileString), mzXmlDir, ".mzXML");
                ApplicationContext.infoMessage("Loading run from file " + mzXmlFile.getAbsolutePath());
                run = MSRun.load(mzXmlFile.getAbsolutePath());
            }
            catch (IOException e)
            {
                throw new  CommandLineModuleExecutionException("Failed to load run from file " + fileString,e);
            }
        }
        FeatureViewerFrame viewerFrame = new FeatureViewerFrame(run);
        viewerFrame.setTitle("Feature in run " + fileString);
        viewerFrame.displayFeature(feature);        
        viewerFrame.setVisible(true);
    }

    protected void displayMetabolite(Map<String, Object> rowMap)  throws CommandLineModuleExecutionException
    {
        Map<String, List<Feature>> fileFeaturesMapThisRow = arrayAnalyzer.createFileFeatureMapForRow(rowId);
        for (String fileString : fileFeaturesMapThisRow.keySet())
        {
            ApplicationContext.infoMessage("Displaying features from run " + fileString);
            for (Feature feature : fileFeaturesMapThisRow.get(fileString))
                displayFeature(feature, fileString);
        }

        if (rowMap.containsKey("SMILES"))
        {
            try
            {
                //todo: this is hokey, renaming the compound based on ion type
                showCompound(rowMap.get("compound") + ": " + rowMap.get("iontype"), (String) rowMap.get("SMILES"));
            }
            catch (CDKException e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
        }
        
    }
}


