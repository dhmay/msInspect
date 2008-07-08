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

public class ArgumentDefinitionFactory
{
    public static final int INTEGER = 0;
    public static final int STRING = 1;
    public static final int DECIMAL = 2;
    public static final int FILE_TO_READ = 3;
    public static final int FILE_TO_WRITE = 4;
    public static final int ENUMERATED = 5;
    public static final int FASTA_FILE = 6;
    public static final int FEATURE_FILE = 7;
    public static final int BOOLEAN = 8;
    public static final int DIRECTORY_TO_READ = 9;
    public static final int DELTA_MASS = 10;
    public static final int MODIFICATION_LIST = 11;
    public static final int DECIMAL_LIST = 12;
    public static final int FILE_TO_READ_LIST = 13;    
    public static final int OTHER = 14;

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
        CommandLineArgumentDefinition result;

        switch (argDefType)
        {
            case INTEGER:
                result = new IntegerArgumentDefinition(argName);
                break;
            case STRING:
                result = new StringArgumentDefinition(argName);
                break;
            case DECIMAL:
                result = new DecimalArgumentDefinition(argName);
                break;
            case FILE_TO_READ:
                result = new FileToReadArgumentDefinition(argName);
                break;
            case FILE_TO_WRITE:
                result = new FileToWriteArgumentDefinition(argName);
                break;
            case ENUMERATED:
                result = new EnumeratedValuesArgumentDefinition(argName);
                break;
            case FASTA_FILE:
                result = new FastaFileArgumentDefinition(argName);
                break;
            case FEATURE_FILE:
                result = new FeatureFileArgumentDefinition(argName);
                break;
            case BOOLEAN:
                result = new BooleanArgumentDefinition(argName);
                break;
            case DIRECTORY_TO_READ:
                result = new DirectoryToReadArgumentDefinition(argName);
                break;
            case DELTA_MASS:
                result = new DeltaMassArgumentDefinition(argName);
                break;
            case MODIFICATION_LIST:
                result = new ModificationListArgumentDefinition(argName);
                break;
            case DECIMAL_LIST:
                result = new DecimalListArgumentDefinition(argName);
                break;
            case FILE_TO_READ_LIST:
                result = new FileToReadListArgumentDefinition(argName);
                break;
            default:
                throw new RuntimeException("Unknown argument type requested, name: " + argName + ", type: " + argDefType);

        }

        result.setRequired(required);
        if (helpText !=  null)
            result.setHelpText(helpText);
        if (defaultValue != null)
            result.setDefaultValue(defaultValue);
        return result;
    }


}
