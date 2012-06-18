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

package org.fhcrc.cpl.toolbox.commandline.arguments;

import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;

import javax.swing.*;
import java.awt.*;

public abstract class BaseArgumentDefinitionImpl implements CommandLineArgumentDefinition
{
    protected String mArgumentName = null;
    protected String mArgumentHelpText = null;

    protected boolean mRequired = true;
    protected Object mDefaultValue = null;


    protected String mDisplayName = null;

    protected boolean mIsAdvanced = false;

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
        mArgumentHelpText = helpText == null ? "" : helpText;
    }

    public BaseArgumentDefinitionImpl(String argumentName,
                                      boolean required,
                                      String helpText)
    {
        this (argumentName, helpText);
        setRequired(required);
    }

    /**
     * Set all the fundamental information about this argument definition, except basic/advanced status
     * @param argumentName
     * @param required
     * @param helpText
     * @param defaultValue
     */
    public BaseArgumentDefinitionImpl(String argumentName,
                                      boolean required,
                                      String helpText,
                                      Object defaultValue)
    {
        this (argumentName, required, helpText);
        setDefaultValue(defaultValue);
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
     * @throws org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException if the argument doesn't validate
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

    public String getDisplayName()
    {
        return mDisplayName;
    }

    public void setDisplayName(String mDisplayName)
    {
        this.mDisplayName = mDisplayName;
    }

    public boolean isAdvanced()
    {
        return mIsAdvanced;
    }

    public void setAdvanced(boolean mIsAdvanced)
    {
        this.mIsAdvanced = mIsAdvanced;
    }

    public JComponent addComponentsForGUISeries(Container parent, JDialog parentDialog, String defaultValue)
    {
        //if a series of files to read (most likely case), handle super-specially.  Else, a big text field
        JTextField argTextField = new JTextField();
        argTextField.setPreferredSize(new Dimension(220, 20));
        argTextField.setMinimumSize(new Dimension(220, 20));

        boolean fieldHasValue = (defaultValue != null && defaultValue.length() > 0);

        if (fieldHasValue)
        {
            defaultValue = defaultValue.replaceAll(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR, " ");
            argTextField.setText(defaultValue);
        }

        GridBagConstraints textFieldGBC = new GridBagConstraints();
        textFieldGBC.anchor = GridBagConstraints.LINE_START;
        textFieldGBC.gridwidth=GridBagConstraints.RELATIVE;
        textFieldGBC.insets = new Insets(0, 0, 0, 0);
        parent.add(argTextField, textFieldGBC);

        return argTextField;
    }

    public JComponent addComponentsForGUI(Container parent, JDialog parentDialog, String defaultValue)
    {
        JPanel fieldPanel = new JPanel();

        JTextField argTextField = new JTextField();
        argTextField.setPreferredSize(new Dimension(120, 20));
        argTextField.setMinimumSize(new Dimension(120, 20));

        if (defaultValue != null && defaultValue.length() > 0)
            argTextField.setText(defaultValue);

        GridBagConstraints argComponentGBC = new GridBagConstraints();
        argComponentGBC.anchor = GridBagConstraints.LINE_START;
        argComponentGBC.gridwidth = GridBagConstraints.REMAINDER;
        argComponentGBC.insets = new Insets(0,0,0,0);

        fieldPanel.add(argTextField, argComponentGBC);

        parent.add(fieldPanel, argComponentGBC);

        return argTextField;
    }

    public String getValueFromGUIComponent(JComponent component)
    {
        return ((JTextField) component).getText();
    }


}
