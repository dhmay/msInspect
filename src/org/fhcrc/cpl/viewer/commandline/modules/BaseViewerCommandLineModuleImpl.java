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

import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.apache.log4j.Logger;

/**
 * Base class for command line modules.  Command line modules do not have to
 * extend this class (they only need to implement CommandLineModule), but it's
 * highly recommended, because this class provides dozens of highly useful
 * convenience methods for argument management
 */
public abstract class BaseViewerCommandLineModuleImpl
    extends BaseCommandLineModuleImpl
{
    static Logger _log = Logger.getLogger(BaseViewerCommandLineModuleImpl.class);



    protected Protein[] getFastaFileArgumentValue(String argumentName)
    {
        return (Protein[]) getArgumentValue(argumentName.toLowerCase());
    }

    protected FeatureSet getFeatureSetArgumentValue(String argumentName)
    {
        return (FeatureSet) getArgumentValue(argumentName.toLowerCase());
    }

    protected MS2Modification[] getModificationListArgumentValue(String argumentName)
    {
        return (MS2Modification[]) getArgumentValue(argumentName.toLowerCase());
    }

    protected FeatureFileArgumentDefinition createUnnamedFeatureFileArgumentDefinition(
                                                                         boolean required,
                                                                         String helpText)
    {
        return new FeatureFileArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,
                required, helpText);
    }

    protected FeatureFileArgumentDefinition createUnnamedSeriesFeatureFileArgumentDefinition(
                                                                         boolean required,
                                                                         String helpText)
    {
        return new FeatureFileArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                required, helpText);
    }
}
