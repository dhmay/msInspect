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
package org.fhcrc.cpl.toolbox.proteomics.feature.filehandler;

import java.io.File;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;

/**
 * base class for featurefile handlers
 */
public abstract class BaseFeatureSetFileHandler
{
    public static final String FILE_TYPE_NAME = "";

    public boolean dumpWindow = false;

    /**
     * Return a string identifying this type of file
     * @return
     */
    public String getFileTypeName()
    {
        return FILE_TYPE_NAME;
    }

    public void setDumpWindow(boolean dumpWindow)
    {
        this.dumpWindow = dumpWindow;
    }


    public static boolean isXMLFile(File inputFile)
            throws IOException
    {
        BufferedReader reader = null;
        boolean result = false;
        try
        {
            reader = new BufferedReader(new FileReader(inputFile));
            String line = reader.readLine();
            if (line.startsWith("<?xml"))
                result = true;
            else if (line.startsWith("<!DOCTYPE"))
                result = true;
        }
        finally
        {
            try
            {
                if (reader != null)
                    reader.close();
            }
            catch (Exception ex)
            {
                //tough luck, presume not XML
            }
        }
        return result;
    }
}
