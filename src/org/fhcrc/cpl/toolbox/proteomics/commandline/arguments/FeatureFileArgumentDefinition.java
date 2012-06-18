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

package org.fhcrc.cpl.toolbox.proteomics.commandline.arguments;

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.FileToReadArgumentDefinition;

import java.io.File;

public class FeatureFileArgumentDefinition extends FileToReadArgumentDefinition
        implements CommandLineArgumentDefinition
{
    public FeatureFileArgumentDefinition(String argumentName)
    {
        super(argumentName);

    }
    public FeatureFileArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public FeatureFileArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    /**
     * Create a file for the filepath and then try to load it as a FeatureSet.
     * Note:  this loading happens at individual argument validation time, before the module
     * has a chance to view all the arguments together.  This can be pretty bad from a usability
     * perspective:  if the arguments as a whole don't mesh, you've wasted a bunch of time loading
     * feature files.  Only use this when you're pretty sure the user's gonna get the arguments
     * right most of the time.
     * @param filePath
     * @return a FeatureSet containing the features loaded from the file.
     */
    public Object convertArgumentValue(String filePath)
            throws ArgumentValidationException
    {
        File inputFile = (File) super.convertArgumentValue(filePath);
        FeatureSet result = null;
        try
        {
            result = new FeatureSet(inputFile);
        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Error loading feature file " + filePath, true);
        }
        if (result.getLoadStatus() != FeatureSet.FEATURESET_LOAD_SUCCESS)
        {
            throw new ArgumentValidationException("Failed to load features from file.  Message: " +
                                                  result.getLoadStatusMessage());
        }
        return result;
    }
}
