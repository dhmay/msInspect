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

import javax.swing.*;
import java.io.File;
import java.awt.*;

public class FileToWriteArgumentDefinition extends FileArgumentDefinition
        implements CommandLineArgumentDefinition
{


    public FileToWriteArgumentDefinition(String argumentName)
    {
        super(argumentName);

    }
    public FileToWriteArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);

    }

    public FileToWriteArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    public FileToWriteArgumentDefinition(String argumentName, boolean required, String help, int fileType)
    {
        super(argumentName, required, help, fileType);
    }

    /**
     * Create a file for the filepath and make sure we can read it
     * @param filePath
     * @return the argument as a File
     */
    public Object convertArgumentValue(String filePath)
            throws ArgumentValidationException
    {
        File fileToWrite = createFileFromRawPath(filePath);
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

    /**
     * Same as base method, but resize the text field
     * @param parent
     * @param parentDialog
     * @param defaultValue
     * @return
     */
    public JComponent addComponentsForGUI(Container parent, JDialog parentDialog, String defaultValue)
    {
        return addComponentsForGUI(parent, parentDialog, defaultValue, false, false);
    }

    public JComponent addComponentsForGUISeries(Container parent, JDialog parentDialog, String defaultValue,
                                                boolean isDir)
    {
        return super.addComponentsForGUISeries(parent, parentDialog, defaultValue, false);
    }    
}
