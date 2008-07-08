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

import java.io.File;

public class FileToWriteArgumentDefinition extends BaseArgumentDefinitionImpl
        implements CommandLineArgumentDefinition
{
    public FileToWriteArgumentDefinition(String argumentName)
    {
        super(argumentName);
        mDataType = ArgumentDefinitionFactory.FILE_TO_WRITE;

    }
    public FileToWriteArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
        mDataType = ArgumentDefinitionFactory.FILE_TO_WRITE;

    }

    /**
     * Create a file for the filepath and make sure we can read it
     * @param filePath
     * @return the argument as a File
     */
    public Object convertArgumentValue(String filePath)
            throws ArgumentValidationException
    {
        File fileToWrite = new File(filePath);
        boolean success = false;
        String failureMessage = "Unable to write file " + filePath;
        if (fileToWrite.exists())
        {
            if (!fileToWrite.isDirectory())
            {
                if (fileToWrite.canWrite())
                    success=true;
            }
            else
                failureMessage = "File " + filePath + " is a directory";
        }
        else
        {
            //annoyingly, the parentFile() of the original fileToWrite will be null if the user
            //specifies a relative path.  So to get the parent file, we have to create a new File
            //using the absolute path, then ask IT for its parent
            File parentFile = new File(fileToWrite.getAbsolutePath()).getParentFile();
            if (parentFile != null)
            {
                if (parentFile.canWrite())
                    success = true;
                else
                    failureMessage = "Can't write to directory " + parentFile.getAbsolutePath();
            }
            else
                failureMessage = "Can't determine parent directory of nonexistent file " + filePath;
        }
        if (!success)
            throw new ArgumentValidationException(failureMessage);      
        return fileToWrite;
    }

    public String getValueDescriptor()
    {
        return "<filepath>";
    }
}
