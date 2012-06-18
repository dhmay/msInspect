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

import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * A list of Strings separated by separatorString.  Returned as List<String>
 * Note: strings are trimmed of whitespace on ends
 */
public class StringListArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public static final String DEFAULT_SEPARATOR_STRING = ",";

    protected String separatorString = DEFAULT_SEPARATOR_STRING;

    protected static Logger _log = Logger.getLogger(StringListArgumentDefinition.class);


    public StringListArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }
    public StringListArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);

    }

    public StringListArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public StringListArgumentDefinition(String argumentName, boolean required, String help, String defaultValue)
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
        List<String> result = new ArrayList<String>();
        if (argumentValue == null)
            return result;
        String[] rawStrings = argumentValue.split(",");
        _log.debug("Parsing string list argument, value: " + argumentValue);
        for (String arg : rawStrings)
        {
            arg = arg.trim();
            if (arg.length() > 0)
            {
                _log.debug("\tParsed arg " + arg);
                result.add(arg);
            }
        }
        return result;
    }

    public String getSeparatorString()
    {
        return separatorString;
    }

    public void setSeparatorString(String separatorString)
    {
        this.separatorString = separatorString;
    }

    public String getValueDescriptor()
    {
        return "<value" + separatorString + "value...>";
    }
}
