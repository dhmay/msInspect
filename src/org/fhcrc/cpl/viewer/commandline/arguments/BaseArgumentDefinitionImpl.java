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

public abstract class BaseArgumentDefinitionImpl implements CommandLineArgumentDefinition
{
    protected String mArgumentName = null;
    protected String mArgumentHelpText = null;

    protected boolean mRequired = true;
    protected Object mDefaultValue = null;

    protected int mDataType = ArgumentDefinitionFactory.OTHER;

    protected String mDisplayName = null;

    public BaseArgumentDefinitionImpl()
    {
    }

    public BaseArgumentDefinitionImpl(String argumentName)
    {
        mArgumentName = argumentName.toLowerCase();
    }

    public BaseArgumentDefinitionImpl(String argumentName, String helpText)
    {
        this(argumentName);
        mArgumentHelpText = helpText.toLowerCase();
    }

    /**
     *
     * @return the name of the argument
     */
    public String getArgumentName()
    {
        return mArgumentName;
    }

    /**
     *
     * @return the name of the argument for display purposes.  For named arguments, this should really be
     * the same as getArgumentName().  But for unnamed arguments or unnamed series arguments, this gives the
     * developer a chance to use a different name in the help.  Be careful, though, not to make this return
     * value too long or too confusing -- the user should still be aware that this is an unnamed argument
     */
    public String getArgumentDisplayName()
    {
        if (mDisplayName != null)
            return mDisplayName;
        
        String result = getArgumentName();
        if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT.equals(result))
        {
            result = "(unnamed)";
        }
        if (CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT.equals(result))
        {
            result = "(unnamed ...)";
        }
        return result;
    }

    /**
     * Set the name of the argument for display purposes
     * @param mDisplayName
     */
    public void setArgumentDisplayName(String mDisplayName)
    {
        this.mDisplayName = mDisplayName;
    }

    /**
     * Convert the argument value to the type of object that the subclass likes
     * to work with.  As part of conversion, perform validation
     *
     * @param argumentValue
     * @return
     * @throws ArgumentValidationException if the argument doesn't validate
     */
    public abstract Object convertArgumentValue(String argumentValue)
            throws ArgumentValidationException;


    public boolean isRequired()
    {
        return mRequired;
    }

    public void setRequired(boolean required)
    {
        this.mRequired = required;
    }

    public String getHelpText()
    {
        return mArgumentHelpText;
    }

    public void setHelpText(String mArgumentHelpText)
    {
        this.mArgumentHelpText = mArgumentHelpText;
    }

    public String getValueDescriptor()
    {
        return "<value>";
    }

    /**
     * Get the default value for this argument, or null if there is none.  However, since
     * null is a valid default in some cases, hasDefaultValue() should be used to determine
     * whether there is a useful default value
     * @return
     */
    public Object getDefaultValue()
    {
        return mDefaultValue;
    }

    /**
     * Get the default value, as a String, for the help text.  In most cases, this will
     * be equivalent to getDefaultValue().toString();
     * @return
     */
    public String getDefaultValueAsString()
    {
        return valueToString(getDefaultValue());
    }

    /**
     * Return a String representing this value
     * @param argValue
     * @return
     */
    public String valueToString(Object argValue)
    {
        return argValue.toString();
    }

    /**
     * Set the default value for this argument
     * @param defaultValue
     */
    public void setDefaultValue(Object defaultValue)
    {
        mDefaultValue = defaultValue;
    }

    /**
     * Does this argument definition have a default value?
     * Default implementation is to check getDefaultValue() for nullness.  Subclasses
     * can override usefully if null is a valid default
     * @return
     */
    public boolean hasDefaultValue()
    {
        return (getDefaultValue() != null);
    }
        
    public int getDataType()
    {
        return mDataType;
    }


    public String getDisplayName()
    {
        return mDisplayName;
    }

    public void setDisplayName(String mDisplayName)
    {
        this.mDisplayName = mDisplayName;
    }
}
