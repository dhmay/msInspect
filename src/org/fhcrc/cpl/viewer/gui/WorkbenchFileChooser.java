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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.Application;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.prefs.Preferences;

/**
 * User: mbellew
 * Date: Sep 7, 2004
 * Time: 2:12:16 PM
 */
public class WorkbenchFileChooser extends JFileChooser
	{
	static File lastDirectory = null;
	static Preferences prefs = Preferences.userNodeForPackage(Application.class);

    public WorkbenchFileChooser()
    {
        super();
        setLocale(Localizer.getLocale());
        //translate as much text as JFileChooser allows
        setApproveButtonText(TextProvider.getText("OPEN"));
        setDialogTitle(TextProvider.getText("OPEN"));
        UIManager.put("FileChooser.acceptAllFileFilterText",
                      TextProvider.getText("ALL_FILES"));
        UIManager.put("FileChooser.lookInLabelText",
                      TextProvider.getText("LOOK_IN"));
        UIManager.put("FileChooser.cancelButtonText",
                      TextProvider.getText("CANCEL"));
        UIManager.put("FileChooser.fileNameLabelText",
                      TextProvider.getText("FILE_NAME_COLON"));
        UIManager.put("FileChooser.filesOfTypeLabelText",
                      TextProvider.getText("FILES_OF_TYPE_COLON"));        
    }

    public int showOpenDialog(Component component) throws HeadlessException
    {
        return showOpenDialog(component, null);
    }

    public int showOpenDialog(Component component, File directory) throws HeadlessException
    {
        boolean directorySet = false;
        if (directory != null)
        {
            if (directory.exists())
            {
                setCurrentDirectory(directory);
                directorySet = true;
            }
        }

        if (!directorySet)
        {
            try
            {
                if (null == lastDirectory)
                    lastDirectory = new File(prefs.get("lastDirectory", "."));
                if (lastDirectory.exists())
                    this.setCurrentDirectory(lastDirectory);
            }
            catch (Exception x)
            {
            }
        }


        int ret = super.showOpenDialog(component);
        if (APPROVE_OPTION == ret)
        {
            File f = getSelectedFile();
            lastDirectory = f.getParentFile();
            prefs.put("lastDirectory", lastDirectory.getPath());
        }
        return ret;
    }

    /**
     * utility method for choosing a file.  Check for file existence.  Display a
     * generic error message and return null if file doesn't exist
     * @param windowTitleText
     * @return
     */
    public static File chooseExistingFile(String windowTitleText)
        {
        WorkbenchFileChooser chooser = new WorkbenchFileChooser();
        if (windowTitleText != null)
            chooser.setDialogTitle(windowTitleText);
        int chooserStatus = chooser.showOpenDialog(ApplicationContext.getFrame());
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return null;            
        File file = chooser.getSelectedFile();

        if (file == null)
            return null;
        if (!file.exists())
            {
            ApplicationContext.infoMessage(TextProvider.getText("FILE_FILE_DOES_NOT_EXIST",
                    file.getAbsolutePath()));
            return null;
            }
        return file;
        }
    }
