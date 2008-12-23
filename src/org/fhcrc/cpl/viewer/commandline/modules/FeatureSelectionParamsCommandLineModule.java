/* 
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
import org.fhcrc.cpl.toolbox.commandline.arguments.DecimalArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.IntegerArgumentDefinition;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

/**
 * Command line module ancestor for modules that need feature selection parameters
 */
public abstract class FeatureSelectionParamsCommandLineModule extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(FeatureSelectionParamsCommandLineModule.class);

    //this FeatureSelector will hold the results of parsing all the filter arguments.
    //If you don't do anything with this, then you haven't used this superclass correctly.
    protected FeatureSet.FeatureSelector featureSelector = null;

    protected void init()
    {
        CommandLineArgumentDefinition[] argDefs =
            {
                    new DecimalArgumentDefinition("minmz", false, "Minimum M/Z value"),
                    new DecimalArgumentDefinition("maxmz", false, "Maximum M/Z value"),
                    new DecimalArgumentDefinition("minmass", false, "Minimum mass"),
                    new DecimalArgumentDefinition("maxmass", false, "Maximum mass"),
                    new IntegerArgumentDefinition("minpeaks", false, "Minimum number of peaks"),
                    new IntegerArgumentDefinition("maxpeaks", false, "Maximum number of peaks"),
                    new IntegerArgumentDefinition("mincharge", false, "Minimum charge"),
                    new IntegerArgumentDefinition("maxcharge", false, "Maximum charge"),
                    new DecimalArgumentDefinition("maxkl", false, "Maximum K/L score"),
                    new DecimalArgumentDefinition("minintensity", false, "Minimum intensity"),
                    new DecimalArgumentDefinition("mintotalintensity", false, "Minimum total intensity"),
                    new DecimalArgumentDefinition("mintime", false, "Minimum time"),
                    new DecimalArgumentDefinition("maxtime", false, "Maximum time"),
                    new IntegerArgumentDefinition("scanFirst", false, "Minimum scan number"),
                    new IntegerArgumentDefinition("scanLast", false, "Maximum scan number"),
                    new IntegerArgumentDefinition("minscans", false, "Minimum number of scans covered"),
                    new DecimalArgumentDefinition("minpprophet", false, "Minimum PeptideProphet score"),
                    new IntegerArgumentDefinition("maxmassdeviationppm", false,
                            "Maximum deviation from nearest theoretical mass cluster, in PPM"),
                    new DecimalArgumentDefinition("maxsumsquaresdist", false, "Maximum sum-squares distance score"),


            };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        featureSelector = new FeatureSet.FeatureSelector();

        if (hasArgumentValue("minMz"))
            featureSelector.setMinMz((float) getDoubleArgumentValue("minMz"));
        if (hasArgumentValue("maxMz"))
            featureSelector.setMaxMz((float) getDoubleArgumentValue("maxMz"));
        if (hasArgumentValue("minMass"))
            featureSelector.setMinMass((float) getDoubleArgumentValue("minMass"));
        if (hasArgumentValue("maxMass"))
            featureSelector.setMaxMass((float) getDoubleArgumentValue("maxMass"));
        if (hasArgumentValue("minPeaks"))
            featureSelector.setMinPeaks(getIntegerArgumentValue("minPeaks"));
        if (hasArgumentValue("maxPeaks"))
            featureSelector.setMaxPeaks(getIntegerArgumentValue("maxPeaks"));
        if (hasArgumentValue("minCharge"))
            featureSelector.setMinCharge(getIntegerArgumentValue("minCharge"));
        if (hasArgumentValue("maxCharge"))
            featureSelector.setMaxCharge(getIntegerArgumentValue("maxCharge"));
        if (hasArgumentValue("maxKL"))
            featureSelector.setMaxKL((float) getDoubleArgumentValue("maxKL"));
        if (hasArgumentValue("minIntensity"))
            featureSelector.setMinIntensity((float) getDoubleArgumentValue("minIntensity"));
        if (hasArgumentValue("minTotalIntensity"))
            featureSelector.setMinTotalIntensity((float) getDoubleArgumentValue("minTotalIntensity"));
        if (hasArgumentValue("minTime"))
            featureSelector.setMinTime((float) getDoubleArgumentValue("minTime"));
        if (hasArgumentValue("maxTime"))
            featureSelector.setMaxTime((float) getDoubleArgumentValue("maxTime"));
        if (hasArgumentValue("scanFirst"))
            featureSelector.setScanFirst(getIntegerArgumentValue("scanFirst"));
        if (hasArgumentValue("scanLast"))
            featureSelector.setScanLast(getIntegerArgumentValue("scanLast"));
        if (hasArgumentValue("minScans"))
            featureSelector.setMinScans(getIntegerArgumentValue("minScans"));
        if (hasArgumentValue("minpprophet"))
            featureSelector.setMinPProphet((float) getDoubleArgumentValue("minpprophet"));
        if (hasArgumentValue("maxmassdeviationppm"))
            featureSelector.setMaxMassDeviationPPM(getIntegerArgumentValue("maxmassdeviationppm"));
        if (hasArgumentValue("maxsumsquaresdist"))
            featureSelector.setMaxSumSquaresDist((float) getDoubleArgumentValue("maxsumsquaresdist"));
    }

}
