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

public class DirectoryToReadArgumentDefinition extends FileArgumentDefinition
        implements CommandLineArgumentDefinition
{
    public DirectoryToReadArgumentDefinition(String argumentName)
    {
        super(argumentName);
    }
    public DirectoryToReadArgumentDefinition(String argumentName, String help)
    {
        super(argumentName, help);
    }

    public DirectoryToReadArgumentDefinition(String argumentName, boolean required, String help)
    {
        super(argumentName, required, help);
    }

    /**
     * Create a file for the filepath and make sure we can read it
     * @param filePath
     * @return the argument as a File
     */
    public Object convertArgumentValue(String filePath)
            throws ArgumentValidationException
    {
        File fileToRead = createFileFromRawPath(filePath);
        if (!fileToRead.exists())
            throw new ArgumentValidationException("Directory " + filePath + " does not exist.");
        if (!fileToRead.isDirectory())
            throw new ArgumentValidationException(filePath + " is not a directory.");
        if (!fileToRead.canRead())
            throw new ArgumentValidationException("Unable to read directory " + filePath);
        return fileToRead;
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
        return addComponentsForGUI(parent, parentDialog, defaultValue, false, true);
    }

    public JComponent addComponentsForGUISeries(Container parent, JDialog parentDialog, String defaultValue,
                                                boolean isDir)
    {
        return super.addComponentsForGUISeries(parent, parentDialog, defaultValue, true);
    }
}
