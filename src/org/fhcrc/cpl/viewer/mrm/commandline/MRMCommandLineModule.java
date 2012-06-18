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
package org.fhcrc.cpl.viewer.mrm.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;

import org.fhcrc.cpl.viewer.mrm.BasicElutionCurveStrategy;
import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.io.File;
import java.awt.*;


/**
 * Commandline module for MRM viewer
 */
public class MRMCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(MRMCommandLineModule.class);

    protected File file;
    protected File outFile;
    protected boolean readSIMs = false;
    protected boolean traceAllFragments = true;
    protected Class peakStrategyClass = BasicElutionCurveStrategy.class;
    protected boolean syncLH = true;

    float meanDaughterMzTolerance = 0.1f;
    float meanPrecursorDiscoveryMzTolerance = 0.01f;
    float SICtolerance = 1.0f;
    float minAreaCutoff = 0.0f;
    float minPeakCutoff = 0.0f;
    public MRMCommandLineModule()
    {
        init();
    }

    public static final String strategies[] = {
       "BasicElutionCurveStrategy",
       //"NoiseSubtractionElutionCurveStrategy",
       "BasicLowIntensityElutionCurveStrategy",
       "ThermoElutionCurveStrategy"     
    };

    protected void init()
    {
        mCommandName = "mrm";
        mShortDescription = "MRMer is used to display, quantify, and edit MRM scans";
        mHelpMessage = "MRMer is a visual editor and analysis tool for MRM assays and assay development. " +
                       "It is written in Java and must be used with Java version 1.6 or higher.<br>" +
                       "\n<br>When you run MRMer from the command line, or if you are creating your "+
                       "own batch or script file, there are a number of options you can place on the " +
                       "java command.  Typically the command will look something like:<br><br>\n\n&nbsp;" +
                       "<tt>java -Xmx500m -jar viewerApp_v"+ TextProvider.getText("MRMER_VERSION")+".jar --mrm "+
                       "</tt><br><br>\n\n followed by one or more of the options below, followed by the mzXML " +
                       "file you wish to analyze.";
        CommandLineArgumentDefinition[] argDefs =
               {
                    createUnnamedFileArgumentDefinition(true,
                            "input mzXML file containing SRM/MRM scans"),
                    new BooleanArgumentDefinition("SELECTED_ION_MONITORING",false,
                            "Set \"SELECTED_ION_MONITORING\" to TRUE if MS1 SIM scans are present and you want to use them exclusively for precursor chromatagrams",readSIMs),
                    new BooleanArgumentDefinition("SYNCLH",false,
                            "Set \"SYNCHLH\" to TRUE to synchronize heavy and light elution regions (applicable in AQUA/SILAC analyses; transition.tsv file must be present)",syncLH),
                    new BooleanArgumentDefinition("TRACE_ALL_FRAGMENTS",false,
                            "Set \"TRACE_ALL_FRAGMENTS\" to TRUE to draw pale elution curves over all product ions, instead of gray spikes",traceAllFragments),
                    new DecimalArgumentDefinition("PRECURSOR_TOLERANCE",false,
                            "Use \"PRECURSOR_TOLERANCE=nnnn\" to set the tolerance for discriminating precursor mass sets",meanPrecursorDiscoveryMzTolerance),
                    new DecimalArgumentDefinition("PRODUCT_TOLERANCE",false,
                            "Use \"PRODUCT_TOLERANCED=nnnn\" to set the tolerance for discriminating between the masses of different product ions",meanDaughterMzTolerance),
                    new EnumeratedValuesArgumentDefinition ("PEAKSTRATEGY", false,"Use peakstrategy to set the elution curve discovery and AUC determination algorithms"
                            ,strategies,"BasicElutionCurveStrategy"),
                    new DecimalArgumentDefinition("AUC_CUTOFF",false,
                            "Use \"AUC_CUTOFF=nnnnn\" to define a minimum peak area to accept",minAreaCutoff),
                    new DecimalArgumentDefinition("PEAK_HEIGHT_CUTOFF",false,
                            "Use \"PEAK_HEIGHT_CUTOFF=nnnnn\" to define minimum peak height (within best curve) to accept",minPeakCutoff),
                    new DecimalArgumentDefinition("SIC_TOLERANCE",false,
                            "Use \"SIC_TOLERANCE=nnnn\" to set the tolerance around the precursor ion for MS1 single ion chromatograms",SICtolerance)
//new FileToWriteArgumentDefinition("out",false, null)
               };
        ((BaseArgumentDefinitionImpl)argDefs[0]).setDisplayName("mzXML file");
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        file = getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
//        outFile = getFileArgumentValue("out");
        readSIMs = getBooleanArgumentValue("SELECTED_ION_MONITORING");
        traceAllFragments = getBooleanArgumentValue("TRACE_ALL_FRAGMENTS");
        meanPrecursorDiscoveryMzTolerance = getFloatArgumentValue("PRECURSOR_TOLERANCE");
        SICtolerance = getFloatArgumentValue("SIC_TOLERANCE");
        meanDaughterMzTolerance = getFloatArgumentValue("PRODUCT_TOLERANCE");
        syncLH = getBooleanArgumentValue("SYNCLH");
        minPeakCutoff = getFloatArgumentValue("PEAK_HEIGHT_CUTOFF");
        minAreaCutoff = getFloatArgumentValue("AUC_CUTOFF");
        String strategy = getStringArgumentValue("peakstrategy");
        if (strategy != null)
        {
            String newSchoolStrategyClassName = strategy;
            //todo: package prefix ought not to be constant
            if (!newSchoolStrategyClassName.contains("."))
                newSchoolStrategyClassName =
                        "org.fhcrc.cpl.viewer.mrm." +
                                strategy;

            try
            {
                peakStrategyClass = Class.forName(newSchoolStrategyClassName);
            }
            catch (ClassNotFoundException e)
            {
               throw new ArgumentValidationException("Could not load class: " + strategy);
            }

        }
    }


    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        try
        {  
            MRMDialog mrmDialog = new MRMDialog(file,meanPrecursorDiscoveryMzTolerance,meanDaughterMzTolerance,SICtolerance,peakStrategyClass,traceAllFragments,syncLH,minPeakCutoff,minAreaCutoff);
            if(getBooleanArgumentValue("SELECTED_ION_MONITORING"))
            {
               mrmDialog.topGraphLabel.setText("SIM Intensities");
               mrmDialog._sim = true;
            }
            //mrmDialog.setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
            mrmDialog.setVisible(true);
           
            mrmDialog.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

//            Thread.currentThread().yield();

              Thread.currentThread().join();
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
    }


}
