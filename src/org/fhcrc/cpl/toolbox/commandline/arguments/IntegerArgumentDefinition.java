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

package org.fhcrc.cpl.toolbox.commandline.arguments;

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;

public class IntegerArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public IntegerArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.INTEGER;

    }
    public IntegerArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.INTEGER;

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
}
