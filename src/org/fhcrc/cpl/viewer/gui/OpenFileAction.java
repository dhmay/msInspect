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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.Application;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.event.ActionEvent;
import java.io.File;


/**
 * User: mbellew
 * Date: May 21, 2004
 * Time: 6:20:03 PM
 */
public class OpenFileAction extends AbstractAction
{
    static JFileChooser chooser;

    public OpenFileAction(JFileChooser chooser)
    {
        super("Open...");
        if (null == chooser)
            chooser = new WorkbenchFileChooser();
        this.chooser = chooser;
    }


    public void actionPerformed(ActionEvent evt)
    {
        JFrame frame = ApplicationContext.getFrame();
        //this file filter is optional.  User can filter on * instead
        chooser.setFileFilter(new FileFilter()
        {
            public boolean accept(File f)
            {
                if (f.isDirectory())
                    return true;
                return f.getName().toLowerCase().endsWith(".mzxml");
            }

            public String getDescription()
            {
                return "*.mzxml";
            }
        }
        );
        int chooserStatus = chooser.showOpenDialog(frame);
        //if user didn't hit OK, ignore
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
            return;
        final File file = chooser.getSelectedFile();
        if (null == file)
            return;

        VerifyFile:
        {
            //dhmay removing explicit check for .mzxml extension
            if (!file.exists() || !file.canRead())
                break VerifyFile;
            Application.getInstance().OpenFile(file.getPath());
            return;
        }

        //JOptionPane.showMessageDialog(frame, "Could not open file.", "Open File", JOptionPane.INFORMATION_MESSAGE, null);
        ApplicationContext.errorMessage("Could not open file: " + file.getPath(), null);
    }
}
