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

import javax.swing.*;
import java.awt.*;

public class BooleanArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    protected static final String[] trueStrings =
            {
                    "y","yes","t","true"
            };
    protected static final String[] falseStrings =
            {
                    "n","no","f","false"
            };

    public BooleanArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }

    public BooleanArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public BooleanArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public BooleanArgumentDefinition(String argumentName, boolean required, String help, boolean defaultValue)
    {
        super(argumentName, required, help, defaultValue);
    }

    /**
     * Try to parse the argument as a Double
     * @param argumentValue
     * @return the argument as a Double
     */
    public Object convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        for (String trueString : trueStrings)
            if (trueString.equalsIgnoreCase(argumentValue))
                return true;
        for (String falseString : falseStrings)
            if (falseString.equalsIgnoreCase(argumentValue))
                return false;
        throw new ArgumentValidationException("Failed to parse true-false value from argument: " + argumentValue);
    }

    public String getValueDescriptor()
    {
        return "<true | false>";
    }

    public String getValueFromGUIComponent(JComponent component)
    {
        return (String) ((JComboBox) component).getSelectedItem();
    }


    public JComponent addComponentsForGUI(Container parent, JDialog parentDialog, String defaultValue)
    {
        JPanel fieldPanel = new JPanel();

        JComboBox booleanComboBox = new JComboBox();
        booleanComboBox.addItem("true");
        booleanComboBox.addItem("false");

        if (defaultValue != null && defaultValue.length() > 0)
            booleanComboBox.setSelectedItem(defaultValue);

        GridBagConstraints argComponentGBC = new GridBagConstraints();
        argComponentGBC.anchor = GridBagConstraints.LINE_START;
        argComponentGBC.gridwidth = GridBagConstraints.REMAINDER;
        argComponentGBC.insets = new Insets(5,0,0,0);        

        fieldPanel.add(booleanComboBox, argComponentGBC);

        parent.add(fieldPanel, argComponentGBC);

        return booleanComboBox;
    }


}
