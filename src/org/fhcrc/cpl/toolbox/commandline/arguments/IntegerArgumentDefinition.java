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

public class IntegerArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public IntegerArgumentDefinition(String argumentName)
    {
        super(argumentName);

    }
    public IntegerArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public IntegerArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public IntegerArgumentDefinition(String argumentName, boolean required, String help, int defaultValue)
    {
        super(argumentName, required, help, defaultValue);
    }




    /**
     * Any String is valid, no-op
     * @param argumentValue
     * @return the argument as an Integer
     */
    public Object convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        Integer result;
        try
        {
            result = Integer.parseInt(argumentValue);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Failed to parse integer value from argument: " + argumentValue);
        }
        return result;
    }

    public String getValueDescriptor()
    {
        return "<integer>";
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
        textField.setPreferredSize(new Dimension(50, 20));
        textField.setMinimumSize(new Dimension(50, 20));
        return textField;
    }
}
