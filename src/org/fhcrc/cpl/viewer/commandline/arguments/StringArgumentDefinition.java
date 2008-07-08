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

import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.BaseArgumentDefinitionImpl;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;

public class StringArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public StringArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.STRING;

    }
    public StringArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.STRING;

    }

    /**
     * Any String is valid, no-op
     * @param argumentValue
     * @return the argument exactly as it came in, as a String
     */
    public Object convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        return argumentValue;
    }

    public String getValueDescriptor()
    {
        return "<value>";
    }
}
