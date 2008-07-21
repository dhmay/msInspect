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

import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import java.util.List;
import java.util.ArrayList;

/**
 * Treats an argument as a comma-separated list of modifications.
 * Each mod should have the format:
 *
 * <residue><massdiff>[V]
 *
 * where V (case-insensitive) is an optional character indicating variable (default is static)
 *
 * for example:
 *
 * M16.5V
 *
 * means a variable modification on Methionine with deltamass 16.5
 */
public class ModificationListArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public ModificationListArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.MODIFICATION_LIST;        
    }

    /**
     * Try to match the argument against the set of allowed values, either with or without case-sensitivity
     * @param argumentValue
     * @return the argument as a String
     */
    public MS2Modification[] convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        if (argumentValue.equalsIgnoreCase("none"))
            return new MS2Modification[0];
        List<MS2Modification> modList = new ArrayList<MS2Modification>();
        try
        {
            String[] modStrings = argumentValue.split(",");
            for (String modString : modStrings)
            {
                modList.add(createModificationFromString(modString));
            }
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(TextProvider.getText("COULDNT_PARSE_ARGUMENT_ARGUMENT", argumentValue));
        }
        return modList.toArray(new MS2Modification[modList.size()]);
    }

    public static MS2Modification safeCreateModificationFromString(String modString)
    {
        try 
        {
            return createModificationFromString(modString);
        }
        catch (ArgumentValidationException e)
        {
            ApplicationContext.errorMessage("Failed to create modification for string: " + modString, e);
        }
        return null;
    }

    public static MS2Modification createModificationFromString(String modString)
            throws ArgumentValidationException
    {
        Character residue = modString.charAt(0);
        if (!Character.isLetter(residue))
            throw new ArgumentValidationException("Illegal residue character '" + residue + "'");
        boolean variable = false;
        if (modString.toLowerCase().endsWith("v"))
            variable = true;
        double massDiff = Double.parseDouble(modString.substring(1,
                variable? modString.length()-1 : modString.length()));
        MS2Modification result = new MS2Modification();
        result.setAminoAcid(("" + residue).toUpperCase());
        result.setVariable(variable);
        result.setMassDiff((float)massDiff);
        return result;
    }

    public String getValueDescriptor()
    {
        return "<residue><massdiff>[V][,...]";
    }

    /**
     * Need to override default handling, since we're working with an array, here.
     * @return
     */
    public String getDefaultValueAsString()
    {
        StringBuffer resultBuf = new StringBuffer();
        MS2Modification[] defaultValue = (MS2Modification[]) getDefaultValue();
        if (defaultValue == null || defaultValue.length < 1)
            return "";
        for (int i=0; i<defaultValue.length; i++)
        {
            MS2Modification mod = defaultValue[i];
            if (i>0)
                resultBuf.append(",");
            resultBuf.append(mod.getAminoAcid() + mod.getMassDiff());
            if (mod.getVariable())
                resultBuf.append("V");
        }
        return resultBuf.toString();
    }


}