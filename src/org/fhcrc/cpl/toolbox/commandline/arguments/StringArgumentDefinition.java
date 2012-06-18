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

public class StringArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public StringArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }
    public StringArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);

    }

    public StringArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public StringArgumentDefinition(String argumentName, boolean required, String help, String defaultValue)
    {
        super(argumentName, required, help, defaultValue);
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
