/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

package org.fhcrc.cpl.viewer.commandline.arguments;

import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;

/**
 * The only handling done here is for argument types not already handled in ArgumentDefinitionFactory
 */
public class ViewerArgumentDefinitionFactory
{
    public static final int FASTA_FILE = ArgumentDefinitionFactory.MAX_TYPE_CONSTANT_VALUE + 1;
    public static final int FEATURE_FILE = ArgumentDefinitionFactory.MAX_TYPE_CONSTANT_VALUE + 2;
    public static final int DELTA_MASS = ArgumentDefinitionFactory.MAX_TYPE_CONSTANT_VALUE + 3;
    public static final int MODIFICATION_LIST = ArgumentDefinitionFactory.MAX_TYPE_CONSTANT_VALUE + 4;

    /**
     * Convenience method for creating unnamed arguments
     * @param argDefType
     * @param required
     * @param helpText
     * @return
     */
    public static CommandLineArgumentDefinition createUnnamedArgumentDefinition(int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return createArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT,
                                        argDefType, required, helpText);
    }


    /**
     * Convenience method for creating unnamed series arguments
     * @param argDefType
     * @param required
     * @param helpText
     * @return
     */
    public static CommandLineArgumentDefinition createUnnamedSeriesArgumentDefinition(int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return createArgumentDefinition(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                                        argDefType, required, helpText);
    }


    /**
     * Create an argument definition with a null default value
     * @param argName
     * @param argDefType
     * @param required
     * @param helpText
     * @return
     */
    public static CommandLineArgumentDefinition createArgumentDefinition(String argName,
                                                                         int argDefType,
                                                                         boolean required,
                                                                         String helpText)
    {
        return createArgumentDefinition(argName, argDefType, required, helpText, null);
    }

    /**
     * Create a CommandLineArgumentDefinition of the appropriate type
     * @param argDefType
     * @param argName
     * @param required
     * @param helpText
     * @return
     */
    public static CommandLineArgumentDefinition createArgumentDefinition(String argName,
                                                                         int argDefType,
                                                                         boolean required,
                                                                         String helpText,
                                                                         Object defaultValue)
    {

        if (argDefType <= ArgumentDefinitionFactory.MAX_TYPE_CONSTANT_VALUE)
            return ArgumentDefinitionFactory.createArgumentDefinition(argName, argDefType, required,
                    helpText, defaultValue);

        CommandLineArgumentDefinition result;

        switch (argDefType)
        {
            case FASTA_FILE:
                result = new FastaFileArgumentDefinition(argName);
                break;
            case FEATURE_FILE:
                result = new FeatureFileArgumentDefinition(argName);
                break;
            case DELTA_MASS:
                result = new DeltaMassArgumentDefinition(argName);
                break;
            case MODIFICATION_LIST:
                result = new ModificationListArgumentDefinition(argName);
                break;
            default:
                throw new IllegalArgumentException("Unknown argument type requested, name: " + argName + ", type: " +
                                                   argDefType);
        }

        result.setRequired(required);
        if (helpText !=  null)
            result.setHelpText(helpText);
        if (defaultValue != null)
            result.setDefaultValue(defaultValue);
        return result;
    }


}
