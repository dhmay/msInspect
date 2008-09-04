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

import java.io.File;

public class FileToReadArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public FileToReadArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.FILE_TO_READ;

    }
    public FileToReadArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.FILE_TO_READ;
        
    }

    /**
     * Create a file for the filepath and make sure we can read it
     * @param filePath
     * @return the argument as a File
     */
    public Object convertArgumentValue(String filePath)
            throws ArgumentValidationException
    {
        File fileToRead = new File(filePath);
        checkFileForReading(fileToRead);
        return fileToRead;
    }

    public static void checkFileForReading(File file)
            throws ArgumentValidationException
    {
        if (!file.exists())
            throw new ArgumentValidationException("File " + file.getAbsolutePath() + " does not exist.");
        if (file.isDirectory())
            throw new ArgumentValidationException(file.getAbsolutePath() + " is a directory.");
        if (!file.isFile())
            throw new ArgumentValidationException(file.getAbsolutePath() + " is not a file.");
        if (!file.canRead())
            throw new ArgumentValidationException("Unable to read file " + file.getAbsolutePath());
    }

    public String getValueDescriptor()
    {
        return "<filepath>";
    }
}