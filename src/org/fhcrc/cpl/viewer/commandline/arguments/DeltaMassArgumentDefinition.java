/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.viewer.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.BaseArgumentDefinitionImpl;

/**
 * Treats an argument as a Delta Mass -- e.g., 0.1da, 5ppm.  Keeps track of both
 * the mass value and the mass type (absolute or PPM)
 */
public class DeltaMassArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    protected int defaultDeltaMassType = AnalyzeICAT.DELTA_MASS_ABSOLUTE;



    public DeltaMassArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ViewerArgumentDefinitionFactory.DELTA_MASS;
    }

    /**
     * Try to match the argument against the set of allowed values, either with or without case-sensitivity
     * @param argumentValue
     * @return the argument as a String
     */
    public DeltaMassWithType convertArgumentValue(String argumentValue)
            throws ArgumentValidationException
    {
        try
        {
            int deltaMassType = defaultDeltaMassType;
            String paramFloatString = argumentValue;
            if (argumentValue.toLowerCase().endsWith("da"))
            {
                deltaMassType = AnalyzeICAT.DELTA_MASS_ABSOLUTE;
                paramFloatString = argumentValue.substring(0, argumentValue.toLowerCase().indexOf("da"));
            }
            if (argumentValue.toLowerCase().endsWith("ppm"))
            {
                deltaMassType = AnalyzeICAT.DELTA_MASS_PPM;
                paramFloatString = argumentValue.substring(0, argumentValue.toLowerCase().indexOf("ppm"));
            }
            return new DeltaMassWithType(Float.parseFloat(paramFloatString),
                    deltaMassType);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException(TextProvider.getText("COULDNT_PARSE_ARGUMENT_ARGUMENT", argumentValue));
        }
    }

    public static class DeltaMassWithType
    {
        protected float mDeltaMass = 0;                                                      
        protected int mDeltaMassType = AnalyzeICAT.DELTA_MASS_ABSOLUTE;

        public DeltaMassWithType(float deltaMass, int deltaMassType)
        {
            mDeltaMass = deltaMass;
            mDeltaMassType = deltaMassType;
        }

        public float getDeltaMass()
        {
            return mDeltaMass;
        }

        public int getDeltaMassType()
        {
            return mDeltaMassType;
        }

        public String toString()
        {
            String suffix = "";
            switch (mDeltaMassType)
            {
                case AnalyzeICAT.DELTA_MASS_ABSOLUTE:
                    suffix = "da";
                    break;
                case AnalyzeICAT.DELTA_MASS_PPM:
                    suffix = "ppm";
                    break;
            }

            return mDeltaMass + suffix;
        }
    }

    public String getValueDescriptor()
    {
        return "<mass value>da|ppm";
    }


    public int getDefaultDeltaMassType()
    {
        return defaultDeltaMassType;
    }

    public void setDefaultDeltaMassType(int defaultDeltaMassType)
    {
        this.defaultDeltaMassType = defaultDeltaMassType;
    }
}
