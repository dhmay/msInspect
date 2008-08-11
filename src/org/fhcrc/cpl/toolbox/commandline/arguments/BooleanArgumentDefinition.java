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
        mDataType = ArgumentDefinitionFactory.BOOLEAN;
    }

    public BooleanArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.BOOLEAN;
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
}
