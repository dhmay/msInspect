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

public class DecimalArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public DecimalArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }
    public DecimalArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);

    }

    public DecimalArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public DecimalArgumentDefinition(String argumentName, boolean required, String help, double defaultValue)
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
        Double result;
        try
        {
            result = Double.parseDouble(argumentValue);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Failed to parse decimal value from argument: " + argumentValue);
        }
        return result;
    }

    public String getValueDescriptor()
    {
        return "<decimal>";
    }

    /**
     * Same as base method, but resize the text field
     * @param parent
     * @param parentDialog
     * @param defaultValue
     * @return
     */
    public JComponent addComponentsForGUI(Container parent, JDialog parentDialog, String defaultValue)
    {
        JTextField textField = (JTextField) super.addComponentsForGUI(parent, parentDialog, defaultValue);
        textField.setPreferredSize(new Dimension(70, 20));
        textField.setMinimumSize(new Dimension(70, 20));
        return textField;
    }


}
