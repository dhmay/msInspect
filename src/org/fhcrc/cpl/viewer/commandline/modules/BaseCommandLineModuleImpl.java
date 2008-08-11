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

package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.PrintWriter;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Base class for command line modules.  Command line modules do not have to
 * extend this class (they only need to implement CommandLineModule), but it's
 * highly recommended, because this class provides dozens of highly useful
 * convenience methods for argument management
 */
public abstract class BaseCommandLineModuleImpl
    implements CommandLineModule
{
    static Logger _log = Logger.getLogger(BaseCommandLineModuleImpl.class);

    //store the command for a given module
    protected String mCommandName;
    protected String mUsageMessage = CommandLineModule.MODULE_USAGE_AUTOMATIC;
    protected String mHelpMessage="";
    protected String mShortDescription="";

    //A map of all argument definitions.  Key = argument name, value = argument definition
    protected Map<String,CommandLineArgumentDefinition> mArgumentDefs;

    //A map of all populated argument names/values.  Null until populated
    protected Map<String, Object> mArgumentValues;

    //A map of all arguments passed to the module (including any extraneous arguments)
    protected Map<String, String> mArgumentValueStrings;



    public String toString()
    {
        return "Commandline module for command " + mCommandName;
    }

    protected String makeHtmlSafe(String plainString)
    {
        if (plainString == null)
            return null;
        return plainString.replaceAll("<","&lt;").replaceAll(">","&gt;");
    }


    /**
     *
     * @return a String telling the user how to invoke this command
     */
    public String getUsage()
    {
        if (null == mUsageMessage ||
            !(CommandLineModule.MODULE_USAGE_AUTOMATIC.equals(mUsageMessage)))
            return mUsageMessage;

        StringBuffer generatedUsage = new StringBuffer();
        generatedUsage.append(getCommandName() + " ");

        String appendix = "";
        for (CommandLineArgumentDefinition definition : getArgumentDefinitionsSortedForDisplay())
        {
            if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT.equals(
                    definition.getArgumentName()))
            {
                appendix = definition.getValueDescriptor();
                continue;
            }
            else if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT.equals(
                    definition.getArgumentName()))
            {
                appendix = definition.getValueDescriptor() + " " +
                        definition.getValueDescriptor() + " ...";
                continue;
            }
            String thisArgUsage = "--" + definition.getArgumentName() + "=";
            String valueString = definition.getValueDescriptor();

           
            thisArgUsage = thisArgUsage + valueString;
            if (!definition.isRequired())
                thisArgUsage = "[" + thisArgUsage + "]";
            generatedUsage.append(thisArgUsage);
            generatedUsage.append(" ");
        }

        generatedUsage.append(appendix);
        return generatedUsage.toString();
    }

    /**
     * Returns the full help message to be displayed by --help <command>
     * @return
     */
    public String getFullHelp()
    {
        StringBuffer result = new StringBuffer();
        result.append(getHelpMessage());

        StringBuffer argumentHelp = new StringBuffer();
        for (CommandLineArgumentDefinition definition : getArgumentDefinitionsSortedForDisplay())
        {
            String helpText = definition.getHelpText();
            if (null == helpText)
                helpText = "(no details provided)";

            String argumentName = definition.getArgumentDisplayName();

            if (definition.hasDefaultValue())
                helpText = helpText + " (default " + definition.getDefaultValueAsString() + ")";
            String requiredString = "";
            if (definition.isRequired())
                requiredString = "*";
            String thisArgHelp = "\t" + requiredString + argumentName + ":\t" +
                    helpText + "\n";
            argumentHelp.append(thisArgHelp);
        }
        result.append("\n\nArgument Details:  ('*' indicates a required parameter)\n" + argumentHelp);


        return result.toString();
    }

    /**
     * Does this module have any 'advanced' arguments?
     * @return
     */
    public boolean hasAdvancedArguments()
    {
        try
        {
            for (CommandLineArgumentDefinition argDef : mArgumentDefs.values())
                if (argDef.isAdvanced())
                    return true;
        }
        catch (Exception e)
        {

        }
        return false;
    }

    /**
     * Returns an HTML fragment containing full help information for this module
     * @return
     */
    public String getHtmlHelpFragment()
    {
        StringBuffer result = new StringBuffer("<a name=\"" + this.getCommandName() + "\"></a>\n<p>\n");
        result.append("<H2>" + getCommandName() + "</H2>\n");
        result.append(getHelpMessage() + "\n<p>\n");
        result.append("<h3>Usage:</h3>\n<p>\n" + makeHtmlSafe("--" + getUsage()) + "<p>");
        String argumentsTitle = "Arguments:";
        if (hasAdvancedArguments())
            argumentsTitle = "Basic " + argumentsTitle;
        result.append("<h3>" + argumentsTitle + "</h3>\n");

        //basic Arguments
        result.append(createArgsTableHTML(getBasicArgumentDefinitions()));

        if (hasAdvancedArguments())
        {
            result.append("\n<p><h3>Advanced Arguments:</h3>\n");
            result.append(createArgsTableHTML(getAdvancedArgumentDefinitions()));            
        }

        result.append("\n");
        return result.toString();
    }

    protected String createArgsTableHTML(CommandLineArgumentDefinition[] argDefs)
    {
        StringBuffer result = new StringBuffer();

        result.append("<table border=\"1\">\n\t<tr><th>Argument</th><th>Usage</th><th>Default</th><th>Description</th></tr>\n");
        for (CommandLineArgumentDefinition definition : argDefs)
        {
            String helpText = definition.getHelpText();
            if (null == helpText || helpText.length() < 1)
                helpText = "(no details provided)";


            String argumentName = definition.getArgumentName();
            if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT.equals(argumentName))
            {
                argumentName = "(unnamed)";
            }
            if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT.equals(argumentName))
            {
                argumentName = "(unnamed ...)";
            }

            if (definition.isRequired())
                argumentName = "<b>" + argumentName + "</b>";

            String usage = definition.getValueDescriptor();
            if (usage == null || usage.length()<1)
                usage = "&nbsp;";


            String defaultValue = null;
            if (definition.hasDefaultValue())
                defaultValue = definition.getDefaultValueAsString();
            String thisArgHelp = "<tr><td>" + argumentName + "</td><td>" +
                    makeHtmlSafe(usage) + "</td><td>" +
                    makeHtmlSafe(defaultValue) + "</td><td>" +
                    makeHtmlSafe(helpText) + "</td></tr>\n";
            result.append(thisArgHelp);
        }
        result.append("</table>\n");

        return result.toString();
    }

    /**
     *
     * @return a String telling the user how to invoke this command; can be detailed
     */
    public String getHelpMessage()
    {
        return mHelpMessage;
    }

    /**
     *
     * @return a String with a quick description of the feature
     */
    public String getShortDescription()
    {
        return mShortDescription;
    }

    /**
     *
     * @return a String that tells Application to invoke this commandlinemodule
     * and not another one.  Had better not collide with another one
     */
    public String getCommandName()
    {
        return mCommandName;
    }

    /**
     * Add an argument definition
     * @param def
     */
    protected void addArgumentDefinition(CommandLineArgumentDefinition def)
    {
        if (mArgumentDefs == null)
            mArgumentDefs = new HashMap<String,CommandLineArgumentDefinition>();
        mArgumentDefs.put(def.getArgumentName().toLowerCase(),def);
    }

    /**
     * Add argument definitions
     * @param defArray
     */
    protected void addArgumentDefinitions(CommandLineArgumentDefinition[] defArray)
    {
        if (mArgumentDefs == null)
            mArgumentDefs = new HashMap<String,CommandLineArgumentDefinition>();
        for (CommandLineArgumentDefinition def : defArray)
            mArgumentDefs.put(def.getArgumentName().toLowerCase(),def);
    }

    /**
     * Add a argument definitions, setting the "advanced" state appropriately
     * @param defArray
     */
    protected void addArgumentDefinitions(CommandLineArgumentDefinition[] defArray, boolean advanced)
    {
        for (CommandLineArgumentDefinition def : defArray)
        {
             def.setAdvanced(advanced);
        }
        addArgumentDefinitions(defArray);
    }



    /**
     * Implements a method in CommandLineModule
     * @return an array of argument definitions
     */
    public CommandLineArgumentDefinition[] getBasicArgumentDefinitions()
    {
        CommandLineArgumentDefinition[] allArgDefs = getArgumentDefinitions();
        Map<String, CommandLineArgumentDefinition> defMap = new HashMap<String, CommandLineArgumentDefinition>();
        for (CommandLineArgumentDefinition argDef : allArgDefs)
        {
            if (!argDef.isAdvanced())
                defMap.put(argDef.getArgumentName(), argDef);
        }
        return sortArgDefsForDisplay(defMap);
    }

    /**
     * Implements a method in CommandLineModule
     * @return an array of argument definitions
     */
    public CommandLineArgumentDefinition[] getAdvancedArgumentDefinitions()
    {
        CommandLineArgumentDefinition[] allArgDefs = getArgumentDefinitions();
        Map<String, CommandLineArgumentDefinition> defMap = new HashMap<String, CommandLineArgumentDefinition>();
        for (CommandLineArgumentDefinition argDef : allArgDefs)
        {
            if (argDef.isAdvanced())
                defMap.put(argDef.getArgumentName(), argDef);
        }
        return sortArgDefsForDisplay(defMap);
    }



    /**
     * Implements a method in CommandLineModule
     * @return an array of argument definitions
     */
    public CommandLineArgumentDefinition[] getArgumentDefinitions()
    {
        if (mArgumentDefs == null)
            mArgumentDefs = new HashMap<String,CommandLineArgumentDefinition>();
        return
            mArgumentDefs.values().toArray(new CommandLineArgumentDefinition[mArgumentDefs.size()]);
    }

    protected CommandLineArgumentDefinition[] sortArgDefsForDisplay(Map<String, CommandLineArgumentDefinition> defMap)
    {

        if (defMap == null)
            return new CommandLineArgumentDefinition[0];
        List<String> mandatoryArgNames =
                new ArrayList<String>();
        List<String> optionalArgNames =
                new ArrayList<String>();
        CommandLineArgumentDefinition unnamedDef = null;

        for (CommandLineArgumentDefinition def : defMap.values())
        {
            if (def.getArgumentName().equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT) ||
                def.getArgumentName().equals(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
            {
                unnamedDef = def;
            }
            else if (def.isRequired())
                mandatoryArgNames.add(def.getArgumentName());
            else optionalArgNames.add(def.getArgumentName());
        }

        CommandLineArgumentDefinition[] result =
                new CommandLineArgumentDefinition[defMap.size()];
        int index = 0;

        String[] mandatoryArray =
                mandatoryArgNames.toArray(new String[mandatoryArgNames.size()]);
        String[] optionalArray =
                optionalArgNames.toArray(new String[optionalArgNames.size()]);

        Arrays.sort(mandatoryArray);
        Arrays.sort(optionalArray);

        if (unnamedDef != null)
            result[index++] = unnamedDef;
        for (String mandatoryArgName : mandatoryArray)
            result[index++] = getArgumentDefinition(mandatoryArgName);
        for (String optionalArgName : optionalArray)
            result[index++] = getArgumentDefinition(optionalArgName);
        return result;
    }

    /**
     * Return the array of argument definitions in the proper display order
     * @return
     */
    public CommandLineArgumentDefinition[] getArgumentDefinitionsSortedForDisplay()
    {
        return sortArgDefsForDisplay(mArgumentDefs);
    }

    /**
     * Implements a method in CommandLineModule
     * @param argumentName
     * @return an argument matching the specified name, or null if not found
     */
    public CommandLineArgumentDefinition getArgumentDefinition(String argumentName)
    {
        if (mArgumentDefs == null)
            mArgumentDefs = new HashMap<String,CommandLineArgumentDefinition>();
        return mArgumentDefs.get(argumentName.toLowerCase());
    }

    /**
     * process a bunch of arguments from the unnamed series definition, converting them
     * to a given argument type
     * @return
     */
    protected Object[] getUnnamedSeriesArgumentValues()
            throws ArgumentValidationException
    {
        return (Object[]) getArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
    }

    /**
     * process an arguments from the unnamed definition, converting it
     * to a given argument type
     * @return
     */
    protected Object getUnnamedArgumentValue()
            throws ArgumentValidationException
    {
        return (Object) getArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    /**
     * Convenience method for getting an array of files specified in an unnamed
     * series, since files are the things most commonly specified that way.
     * I do NOT recommend creating separate convenience methods for each type --
     * that'd be bloat
     * @return
     * @throws ArgumentValidationException
     */
    protected File[] getUnnamedSeriesFileArgumentValues()
            throws ArgumentValidationException
    {
        Object[] objectArray = getUnnamedSeriesArgumentValues();
        if (objectArray == null)
            return null;
        File[] result = new File[objectArray.length];
        for (int i=0; i<objectArray.length; i++)
        {
            result[i] = (File) objectArray[i];
            FileToReadArgumentDefinition.checkFileForReading(result[i]);            
        }
        return result;
    }
    

    /**
     * Digest raw string argument values.  Make sure they're all there and convert properly
     * @param argumentValueMap
     * @throws ArgumentValidationException
     */
    public void digestArguments(Map<String, String> argumentValueMap)
            throws ArgumentValidationException
    {
        mArgumentValueStrings = argumentValueMap;

        mArgumentValues = new HashMap<String, Object>();
        //Capture only the args we care about
        for (CommandLineArgumentDefinition argDef : getArgumentDefinitions())
        {
            String paramVal = argumentValueMap.get(argDef.getArgumentName().toLowerCase());
            //If we're missing a required argument, warn and die
            if (paramVal == null)
            {
                //if not required, then whatever.  If required, complain
                if (!argDef.isRequired())
                    continue;
                else
                {
                    if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT.equals(argDef.getArgumentName()))
                        throw new ArgumentValidationException("Missing Required unnamed parameter for command " + getCommandName());
                    else if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT.equals(argDef.getArgumentName()))
                        throw new ArgumentValidationException("Missing Required unnamed series of parameters for command " + getCommandName());
                    else
                        throw new ArgumentValidationException("Missing Required parameter '" + argDef.getArgumentName() +
                                "' for command " + getCommandName());
                }
            }
            else
            {
                if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT.equals(argDef.getArgumentName()))
                {
                    _log.debug("Unnamed series parameter value: \n" + paramVal);
                    String[] unnamedParamStrings = paramVal.split(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                    List<Object> unnamedParamList = new ArrayList<Object>();
                    _log.debug("Individual parameters:");
                    for (String unnamedParamString : unnamedParamStrings)
                    {
                        if (unnamedParamString != null &&
                            unnamedParamString.length() > 0)
                        {
                            unnamedParamList.add(argDef.convertArgumentValue(unnamedParamString));
                            _log.debug(unnamedParamList.get(unnamedParamList.size()-1));
                        }
                    }
                    mArgumentValues.put(argDef.getArgumentName(),
                                        unnamedParamList.toArray(new Object[unnamedParamList.size()]));
                }
                else
                {
                    Object argValue = argDef.convertArgumentValue(paramVal);
//System.err.println("adding arg: " + argDef.getArgumentName() + " = " + argValue);
                    mArgumentValues.put(argDef.getArgumentName(), argValue);
                }
            }
        }
        assignArgumentValues();
    }

    /**
     * Invoke the CommandLineModule, using the argument values specified in
     * the argumentValues array.  Essentially this should just call digestArguments()
     * and then execute()
     * @param argumentValues
     */
    public void invoke(Map<String, String> argumentValues)
            throws ArgumentValidationException, CommandLineModuleExecutionException
    {
        digestArguments(argumentValues);
        execute();
    }


    /**
     * Asserts the presence of an argument.  If absent, throws AVException.  This is for
     * situations in which the presence of certain arguments makes other 'optional' arguments
     * mandatory
     * @param argumentName
     * @throws ArgumentValidationException
     */
    protected void assertArgumentPresent(String argumentName)
            throws ArgumentValidationException
    {
        if (!hasArgumentValue(argumentName))
            throw new ArgumentValidationException("The argument '" + argumentName +
                                                  "' is required for the command usage you have attempted");
    }

    /**
     * Asserts the presence of an argument.  If absent, throws AVException.  This is for
     * situations in which the presence of certain arguments makes other 'optional' arguments
     * mandatory.  This version allows the user to specify which argument requires the other one
     * @throws ArgumentValidationException
     */
    protected void assertArgumentPresent(String requiredArgumentName,
                                         String requiredByArgumentName)
            throws ArgumentValidationException
    {
        if (!hasArgumentValue(requiredArgumentName))
            throw new ArgumentValidationException("The argument '" + requiredArgumentName +
                                                  "' is required if the argument '" + requiredByArgumentName + "' is provided or has the provided value.");
    }

   /**
     * Asserts the absence of an argument.  If presence, throws AVException.  This is for
     * situations in which the presence of certain arguments makes other 'optional' arguments
     * meaningless
     * @param argumentName
     * @throws ArgumentValidationException
     */
    protected void assertArgumentAbsent(String argumentName)
            throws ArgumentValidationException
    {
        if (hasArgumentValue(argumentName))
            throw new ArgumentValidationException("The argument '" + argumentName +
                                                  "' cannot be used with the other arguments you have provided");
    }

   /**
     * Asserts the absence of an argument.  If presence, throws AVException.  This is for
     * situations in which the presence of certain arguments makes other 'optional' arguments
     * meaningless.  This version allows the user to specify which argument made the arg
    * problematic
     * @param argumentName
     * @throws ArgumentValidationException
     */
    protected void assertArgumentAbsent(String argumentName, String requiredAbsentByArgument)
            throws ArgumentValidationException
    {
        if (hasArgumentValue(argumentName))
            throw new ArgumentValidationException("The argument '" + argumentName +
                    "' cannot be used if the argument '" + requiredAbsentByArgument + "' is provided or has the provided value.");
    }

    /**
     * Convenience method to tell whether an argument map contains an argument
     * @param argumentName
     * @return
     */
    protected boolean hasArgumentValue(String argumentName)
    {
        return mArgumentValues.containsKey(argumentName.toLowerCase());
    }

    protected boolean hasUnnamedSeriesArgumentValue()
    {
        return hasArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
    }

    protected boolean hasUnnamedArgumentValue()
    {
        return hasArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    //convenience methods for casting argument values

    protected String getStringArgumentValue(String argumentName)
    {
        return (String) getArgumentValue(argumentName.toLowerCase());
    }

    protected double getDoubleArgumentValue(String argumentName)
    {
        return (Double) getArgumentValue(argumentName.toLowerCase());
    }

    protected double[] getDoubleArrayArgumentValue(String argumentName)
    {
        return (double[]) getArgumentValue(argumentName.toLowerCase());
    }

    protected int getIntegerArgumentValue(String argumentName)
    {
        return (Integer) getArgumentValue(argumentName.toLowerCase());
    }

    protected boolean getBooleanArgumentValue(String argumentName)
    {
        return (Boolean) getArgumentValue(argumentName.toLowerCase());
    }

    protected File getFileArgumentValue(String argumentName)
    {
        return (File) getArgumentValue(argumentName.toLowerCase());
    }

    protected File[] getFileArrayArgumentValue(String argumentName)
    {
        return (File[]) getArgumentValue(argumentName.toLowerCase());
    }

    protected File getUnnamedFileArgumentValue()
    {
        return getFileArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    protected Protein[] getFastaFileArgumentValue(String argumentName)
    {
        return (Protein[]) getArgumentValue(argumentName.toLowerCase());
    }

    protected FeatureSet getFeatureSetArgumentValue(String argumentName)
    {
        return (FeatureSet) getArgumentValue(argumentName.toLowerCase());
    }

    protected DeltaMassArgumentDefinition.DeltaMassWithType
    getDeltaMassArgumentValue(String argumentName)
    {
        return (DeltaMassArgumentDefinition.DeltaMassWithType)
                getArgumentValue(argumentName.toLowerCase());
    }

    protected MS2Modification[] getModificationListArgumentValue(String argumentName)
    {
        return (MS2Modification[]) getArgumentValue(argumentName.toLowerCase());
    }


    protected Map<String, Object> getArgumentValues()
    {
        return mArgumentValues;
    }

    protected Object getArgumentValue(String argumentName)
    {
        Object result = mArgumentValues.get(argumentName);
        if (result == null)
        {
            if (mArgumentDefs.containsKey(argumentName))
                result = getArgumentDefinition(argumentName).getDefaultValue();
            else
                throw new RuntimeException("No such argument definition: " + argumentName);
        }
        return result;
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


    //convenience methods

    protected IntegerArgumentDefinition createIntegerArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (IntegerArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.INTEGER, required, helpText);
    }

    protected IntegerArgumentDefinition createIntegerArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         int defaultValue)
    {
        return (IntegerArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.INTEGER, required, helpText, defaultValue);
    }


    protected StringArgumentDefinition createStringArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (StringArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.STRING, required, helpText);
    }

    protected StringArgumentDefinition createStringArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         String defaultValue)
    {
        return (StringArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.STRING, required, helpText, defaultValue);
    }


    protected DecimalArgumentDefinition createDecimalArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (DecimalArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName, ArgumentDefinitionFactory.DECIMAL, required, helpText);
    }

    protected DecimalArgumentDefinition createDecimalArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         double defaultValue)
    {
        return (DecimalArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.DECIMAL, required, helpText, defaultValue);
    }

    protected DecimalListArgumentDefinition createDecimalListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (DecimalListArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName, ArgumentDefinitionFactory.DECIMAL_LIST, required, helpText);
    }

    protected DecimalListArgumentDefinition createDecimalListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         double[] defaultValue)
    {
        return (DecimalListArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.DECIMAL_LIST, required, helpText, defaultValue);
    }

    protected DecimalListArgumentDefinition createDecimalListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         float[] defaultValue)
    {
        double[] defaultValueDoubleArray = new double[defaultValue.length];
        for (int i=0; i<defaultValue.length; i++)
            defaultValueDoubleArray[i] = defaultValue[i];
        return createDecimalListArgumentDefinition(argName, required,helpText,
                                                   defaultValueDoubleArray);
    }


    protected FileToReadListArgumentDefinition createFileToReadListArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (FileToReadListArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.FILE_TO_READ_LIST, required, helpText);
    }


    protected BooleanArgumentDefinition createBooleanArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (BooleanArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.BOOLEAN, required, helpText);
    }

    protected BooleanArgumentDefinition createBooleanArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         boolean defaultValue)
    {
        return (BooleanArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.BOOLEAN, required, helpText, defaultValue);
    }

    protected FileToReadArgumentDefinition createFileToReadArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (FileToReadArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.FILE_TO_READ, required, helpText);
    }
    
    protected FileToWriteArgumentDefinition createFileToWriteArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText)
    {
        return (FileToWriteArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.FILE_TO_WRITE, required, helpText);
    }

    protected DirectoryToReadArgumentDefinition
            createDirectoryToReadArgumentDefinition(String argName,
                                                    boolean required,
                                                    String helpText)
    {
        return (DirectoryToReadArgumentDefinition) ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                ArgumentDefinitionFactory.DIRECTORY_TO_READ, required, helpText);
    }

    protected EnumeratedValuesArgumentDefinition createEnumeratedArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         String[] allowedValues,
                                                                         String defaultValue)
    {
        EnumeratedValuesArgumentDefinition result =
                        (EnumeratedValuesArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                        ArgumentDefinitionFactory.ENUMERATED, required, helpText, defaultValue);
        result.setEnumeratedValues(allowedValues);
        return result;
    }

    protected EnumeratedValuesArgumentDefinition createEnumeratedArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String helpText,
                                                                         String[] allowedValues)
    {
        return createEnumeratedArgumentDefinition(argName, required, helpText,
                allowedValues, null);
    }

    /**
     * Create an EnumeratedValues argument definition, and handily provide a short
     * explanation for each allowed value
     * @param argName
     * @param required
     * @param allowedValues
     * @param valueExplanations
     * @return
     */
    protected EnumeratedValuesArgumentDefinition createEnumeratedArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String[] allowedValues,
                                                                         String[] valueExplanations)
    {
        return createEnumeratedArgumentDefinition(argName, required, allowedValues, valueExplanations, null);
    }

    /**
     * Create an EnumeratedValues argument definition, and handily provide a short
     * explanation for each allowed value
     * @param argName
     * @param required
     * @param allowedValues
     * @param valueExplanations
     * @return
     */
    protected EnumeratedValuesArgumentDefinition createEnumeratedArgumentDefinition(String argName,
                                                                         boolean required,
                                                                         String[] allowedValues,
                                                                         String[] valueExplanations,
                                                                         String defaultValue)
    {
        StringBuffer helpTextBuf = new StringBuffer("\n");
        for (int i=0; i< allowedValues.length; i++)
        {
            helpTextBuf.append("\t\t" + allowedValues[i] + ":\t" + valueExplanations[i] + "\n");
        }
        return createEnumeratedArgumentDefinition(argName, required, helpTextBuf.toString(),
                allowedValues, defaultValue);
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

    protected DeltaMassArgumentDefinition createDeltaMassArgumentDefinition(
            String argName, boolean required, String helpText,
            DeltaMassArgumentDefinition.DeltaMassWithType defaultDeltaMass)
    {
        return (DeltaMassArgumentDefinition)
                ViewerArgumentDefinitionFactory.createArgumentDefinition(argName,
                ViewerArgumentDefinitionFactory.DELTA_MASS, required, helpText,
                defaultDeltaMass);
    }

    protected DeltaMassArgumentDefinition createDeltaMassArgumentDefinition(
            String argName, boolean required, String helpText)
    {
        return createDeltaMassArgumentDefinition(argName, required, helpText, null);
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


    /*
    * dhmay adding to centralize PrintWriter creation
    * Utility method: given an output filename, checks if null.
    * If not, returns a PrintWriter
    * on that file. If so, returns a PrintWriter on System.err
    *
    * This is purely a convenience method
    *
    * @param outFileName the output filename specified on the command line
    * @return an appropriate printwriter
    */
    protected static PrintWriter getPrintWriter(File outFile)
            throws java.io.FileNotFoundException
    {
      PrintWriter pw = null;
      if (null != outFile)
      {
        try
        {
          pw = new PrintWriter(new FileOutputStream(outFile));
        }
        catch (java.io.FileNotFoundException e)
        {
          ApplicationContext.infoMessage("Error creating PrintWriter for file " +
                                          outFile.getAbsolutePath());
          throw e;
        }
      }
      else
        pw = new PrintWriter(System.out);

      return pw;
    }


    /**
     * Return the original argument name-value pairs that were passed into this module via digestArguments()
     * @return
     */
    public Map<String, String> getArgumentValueStrings()
    {
        return mArgumentValueStrings;
    }
}
