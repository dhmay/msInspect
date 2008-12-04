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

import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.commandline.arguments.FeatureFileArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.FastaFileArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ModificationListArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.ViewerArgumentDefinitionFactory;
import org.fhcrc.cpl.toolbox.commandline.BaseCommandLineModuleImpl;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;
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




    protected FeatureFileArgumentDefinition createFeatureFileArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (FeatureFileArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName, ViewerArgumentDefinitionFactory.FEATURE_FILE, required, helpText);
    }
    protected FastaFileArgumentDefinition createFastaFileArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (FastaFileArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ViewerArgumentDefinitionFactory.FASTA_FILE, required, helpText);
    }

    protected ModificationListArgumentDefinition createModificationListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         MS2Modification[] defaultValue)
    {
        return (ModificationListArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ViewerArgumentDefinitionFactory.MODIFICATION_LIST, required, helpText,
                        defaultValue);
    }

    protected ModificationListArgumentDefinition createModificationListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return createModificationListArgumentDefinition(argName, required,
                helpText, null);
    }

    /**
     * Cover method for ArgumentDefinitionFactory method
     * @param argDefType
     * @param required
     * @param helpText
     * @return
     */
    protected CommandLineArgumentDefinition createUnnamedArgumentDefinition(int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return ViewerArgumentDefinitionFactory.createUnnamedArgumentDefinition(argDefType, required, helpText);
    }


    protected CommandLineArgumentDefinition createUnnamedFileArgumentDefinition(
                                                                         boolean required,
                                                                         String helpText)
    {
        return ViewerArgumentDefinitionFactory.createUnnamedArgumentDefinition(
                ArgumentDefinitionFactory.FILE_TO_READ,
                required, helpText);
    }


    /**
     * Cover method for ArgumentDefinitionFactory method
     * @param required
     * @param helpText
     * @return
     */
    protected CommandLineArgumentDefinition createUnnamedSeriesFileArgumentDefinition(
                                                                         boolean required,
                                                                         String helpText)
    {
        return ViewerArgumentDefinitionFactory.createUnnamedSeriesArgumentDefinition(ArgumentDefinitionFactory.FILE_TO_READ,
                required, helpText);
    }

    /**
     * Cover method for ArgumentDefinitionFactory method
     * @param argDefType
     * @param required
     * @param helpText
     * @return
     */
    protected CommandLineArgumentDefinition createUnnamedSeriesArgumentDefinition(int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return ViewerArgumentDefinitionFactory.createUnnamedSeriesArgumentDefinition(
                argDefType, required, helpText);
    }


    /**
     * Cover method for ArgumentDefinitionFactory method, no default value
     * @param argDefType
     * @param argName
     * @param required
     * @param helpText
     * @return
     */
    protected CommandLineArgumentDefinition createArgumentDefinition(String argName,
                                                                         int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return ViewerArgumentDefinitionFactory.createArgumentDefinition(argName, argDefType, required,
                helpText);
    }

    /**
     * Cover method for ArgumentDefinitionFactory method
     * @param argDefType
     * @param argName
     * @param required
     * @param helpText
     * @return
     */
    protected CommandLineArgumentDefinition createArgumentDefinition(String argName,
                                                                         int argDefType,
                                                                         boolean required,
                                                                         String helpText,
                                                                         Object defaultValue)
    {
        return ViewerArgumentDefinitionFactory.createArgumentDefinition(argName, argDefType, required,
                helpText, defaultValue);
    }

}
