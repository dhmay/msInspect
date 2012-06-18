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
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.filehandler.TabWriter;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;


import java.io.File;
import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class ProteinFractionsSpreadsheetCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(ProteinFractionsSpreadsheetCLM.class);

    protected File[] featureFiles;

    protected File outFile;
    protected File protXmlFile;

    double minPeptideProphet = 0.0;
    double minProteinProphet = 0.0;

    protected boolean groupLevel = false;



    public ProteinFractionsSpreadsheetCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "proteinfractionsspreadsheet";
        mShortDescription = "Create a spreadsheet assigning proteins to fractions";
        mHelpMessage = "asdfasdf";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true, "MS2 feature files"),
                        new FileToReadArgumentDefinition("protxml",true,"ProtXML File"),
                        new DecimalArgumentDefinition("minpprophet",false,"Minimum peptideprophet",minPeptideProphet),
                        new DecimalArgumentDefinition("minproteinprophet",false,"Minimum proteinprophet",minProteinProphet),
                        new BooleanArgumentDefinition("grouplevel", false, "Group-level? (default is accesion-number level)", groupLevel),
                        new FileToWriteArgumentDefinition("out",true, null),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureFiles = this.getUnnamedSeriesFileArgumentValues();
        protXmlFile = getFileArgumentValue("protxml");

        minPeptideProphet = getDoubleArgumentValue("minpprophet");
        minProteinProphet = getDoubleArgumentValue("minproteinprophet");

        outFile = getFileArgumentValue("out");
        groupLevel = getBooleanArgumentValue("grouplevel");
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {
            Map peptideProteinThingyMap = null;

            if (groupLevel)
            {
                Map<String, Set<String>> peptideProteinMap =
                        ProteinUtilities.loadPeptideProteinMapFromProtXML(protXmlFile,minProteinProphet);
                peptideProteinThingyMap = peptideProteinMap;
            }
            else
            {
                Map<String, Set<Integer>> peptideProteinGroupMap =
                        ProteinUtilities.loadPeptideProteinGroupMapFromProtXML(protXmlFile, minProteinProphet);
                peptideProteinThingyMap = peptideProteinGroupMap;
            }

            Map<Object, Set<File>> proteinFileMap = new HashMap<Object, Set<File>>();

            List<String> columnNamesList = new ArrayList<String>();
            columnNamesList.add("protein");

            List<Float> numProteinsPerFile = new ArrayList<Float>();
            List<Float> numQuantProteinsPerFile = new ArrayList<Float>();
            for (File featureFile : featureFiles)
            {

                columnNamesList.add(featureFile.getName());

                Set proteinsThisFile = new HashSet();
                Set quantProteinsThisFile = new HashSet();

                FeatureSet featureSet = new FeatureSet(featureFile);
                for (Feature feature : featureSet.getFeatures())
                {
                    if (minPeptideProphet > 0 &&
                            (MS2ExtraInfoDef.getPeptideProphet(feature) < minPeptideProphet))
                        continue;
                    String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
                    if (peptide == null)
                        continue;
                    Set proteinsThisPeptide = (Set) peptideProteinThingyMap.get(peptide);
                    if (proteinsThisPeptide != null)
                    {
                        proteinsThisFile.addAll(proteinsThisPeptide);
                        if (IsotopicLabelExtraInfoDef.hasRatio(feature))
                            quantProteinsThisFile.addAll(proteinsThisPeptide);
                    }
                }
                for (Object protein : proteinsThisFile)
                {
                    String key = protein.toString();
                    Set<File> filesThisProtein = proteinFileMap.get(key);
                    if (filesThisProtein == null)
                    {
                        filesThisProtein = new HashSet<File>();
                        proteinFileMap.put(key,filesThisProtein);
                    }
                    filesThisProtein.add(featureFile);
                }
                numProteinsPerFile.add((float) proteinsThisFile.size());
                numQuantProteinsPerFile.add((float) quantProteinsThisFile.size());

            }

            ApplicationContext.infoMessage("Mean proteins per file: " + BasicStatistics.mean(numProteinsPerFile));
            ApplicationContext.infoMessage("Mean QUANT proteins per file: " + BasicStatistics.mean(numQuantProteinsPerFile));

            List<Float> filesPerProtein = new ArrayList<Float>();

            TabWriter tw = new TabWriter(columnNamesList.toArray(new String[columnNamesList.size()]));
            tw.setOutFile(outFile);

            for (Object protein :  proteinFileMap.keySet())
            {
                Map<String,Object> row = new HashMap<String,Object>();
                row.put("protein",protein);
                filesPerProtein.add((float)proteinFileMap.get(protein).size());
                for (File file : proteinFileMap.get(protein))
                {
                    row.put(file.getName(), "X");
                }
                tw.addRow(row);
                
            }
            tw.write();

            ApplicationContext.infoMessage("Mean files per protein: " + BasicStatistics.mean(filesPerProtein));
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }

}
