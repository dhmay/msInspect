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

/**
 * Validates arguments against an explicitly enumerated list of String values.  Can be done
 * with or without case-sensitivity
 */
public class EnumeratedValuesArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{

    protected String[] mEnumeratedValues;
    protected boolean mCaseSensitive = false;


    public EnumeratedValuesArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.ENUMERATED;

    }

    /**
     * By default, validation is done without case sensitivity
     * @param argumentName
     * @param enumeratedValues
     */
    public EnumeratedValuesArgumentDefinition(String argumentName,
                                              String[] enumeratedValues)
    {
        this(argumentName);
        mEnumeratedValues = enumeratedValues;

    }
    public EnumeratedValuesArgumentDefinition(String argumentName, String[] enumeratedValues, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.ENUMERATED;

        mEnumeratedValues = enumeratedValues;
    }

    public EnumeratedValuesArgumentDefinition(String argumentName,
                                              String[] enumeratedValues,
                                              boolean caseSensitive)
    {
        this(argumentName);
        mEnumeratedValues = enumeratedValues;
        mCaseSensitive = caseSensitive;
    }

    /**
     * Try to match the argument against the set of allowed values, either with or without case-sensitivity
     * @param argumentValue
     * @return the argument as a String
     */
    public String convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        return mEnumeratedValues[getIndexForArgumentValue(argumentValue)];
    }

    /**
     * Get the index into the enumerated values array for the element that matches
     * the passed-in argumentValue.  If nothing matches, throw an exception
     * @param argumentValue
     * @return
     */
    public int getIndexForArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        int result = -1;
        for (int i=0; i<mEnumeratedValues.length; i++)
        {
            String testValue = mEnumeratedValues[i];
            if (mCaseSensitive)
            {
                if (testValue.equals(argumentValue))
                {
                    result = i;
                    break;
                }
            }
            else
            {
                if (testValue.equalsIgnoreCase(argumentValue))
                {
                    result = i;
                    break;
                }
            }
        }

        if (result < 0)
        {
            StringBuffer exceptionSB =
                    new StringBuffer("Argument value " + argumentValue +
                                     " is not an allowed value for parameter \"" +
                                     getArgumentName() + "\".\nAllowed values are:\n");
            for (String allowedValue : mEnumeratedValues)
                exceptionSB.append("\t" + allowedValue + "\n");
            throw new ArgumentValidationException(exceptionSB.toString());
        }

        return result;
    }


    public String[] getEnumeratedValues()
    {
        return mEnumeratedValues;
    }

    public void setEnumeratedValues(String[] mEnumeratedValues)
    {
        this.mEnumeratedValues = mEnumeratedValues;
    }

    public String getValueDescriptor()
    {
        String valueString = "<";
        for (int i=0; i<getEnumeratedValues().length; i++)
        {
            valueString = valueString + getEnumeratedValues()[i];
            if (i < (getEnumeratedValues().length - 1))
                valueString = valueString + " | ";
        }
        valueString = valueString + ">";

        return valueString;
    }
}
