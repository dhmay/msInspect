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
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

/**
 * Command line module ancestor for modules that need feature selection parameters
 */
public abstract class FeatureSelectionParamsCommandLineModule extends BaseCommandLineModuleImpl
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
                    createDecimalArgumentDefinition("minmz", false, "Minimum M/Z value"),
                    createDecimalArgumentDefinition("maxmz", false, "Maximum M/Z value"),
                    createDecimalArgumentDefinition("minmass", false, "Minimum mass"),
                    createDecimalArgumentDefinition("maxmass", false, "Maximum mass"),
                    createIntegerArgumentDefinition("minpeaks", false, "Minimum number of peaks"),
                    createIntegerArgumentDefinition("maxpeaks", false, "Maximum number of peaks"),
                    createIntegerArgumentDefinition("mincharge", false, "Minimum charge"),
                    createIntegerArgumentDefinition("maxcharge", false, "Maximum charge"),
                    createDecimalArgumentDefinition("maxkl", false, "Maximum K/L score"),
                    createDecimalArgumentDefinition("minintensity", false, "Minimum intensity"),
                    createDecimalArgumentDefinition("mintotalintensity", false, "Minimum total intensity"),
                    createDecimalArgumentDefinition("mintime", false, "Minimum time"),
                    createDecimalArgumentDefinition("maxtime", false, "Maximum time"),
                    createIntegerArgumentDefinition("scanFirst", false, "Minimum scan number"),
                    createIntegerArgumentDefinition("scanLast", false, "Maximum scan number"),
                    createIntegerArgumentDefinition("minscans", false, "Minimum number of scans covered"),
                    createDecimalArgumentDefinition("minpprophet", false, "Minimum PeptideProphet score"),
                    createIntegerArgumentDefinition("maxmassdeviationppm", false,
                            "Maximum deviation from nearest theoretical mass cluster, in PPM"),
                    createDecimalArgumentDefinition("maxsumsquaresdist", false, "Maximum sum-squares distance score"),


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
