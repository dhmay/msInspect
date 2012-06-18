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

import java.util.List;
import java.util.ArrayList;

public class DecimalListArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public DecimalListArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }
    public DecimalListArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public DecimalListArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public DecimalListArgumentDefinition(String argumentName, boolean required, String help, double[] defaultValue)
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
        if (argumentValue == null)
            return null;
        double[] result;

        try
        {
            List<Double> resultList = new ArrayList<Double>();
            String[] chunks = argumentValue.split(",");
            for (String chunk : chunks)
            {
                if (chunk != null && chunk.length() > 0)
                {
                    Double bucketsize = new Double(chunk);
                    resultList.add(bucketsize);
                }
            }
            result = new double[resultList.size()];
            for (int i=0; i<resultList.size(); i++)
            {
                result[i] = resultList.get(i);
            }
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Error while parsing parameter string " + argumentValue,e);
        }

        return result;
    }

    public String getValueDescriptor()
    {
        return "<decimal>,<decimal>,...";
    }
}
