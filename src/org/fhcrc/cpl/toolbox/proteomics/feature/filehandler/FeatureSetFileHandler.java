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

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * An interface to be implemented by classes that will each handle a particular
 * type of feature file. 
 */
public interface FeatureSetFileHandler
{
    /**
     * Controls whether to dump a window of intensity values around each feature
     * @param dumpWindow
     */
    public void setDumpWindow(boolean dumpWindow);

    /**
     * Return a string identifying this type of file
     * @return
     */
    public String getFileTypeName();

    /**
     * Load a FeatureSet
     * @param featureFile
     * @return
     * @throws IOException
     */
    public FeatureSet loadFeatureSet(File featureFile)
        throws IOException;

    /**
     * Save a FeatureSet
     * @param outputFile
     * @throws IOException
     */
    public void saveFeatureSet(FeatureSet featureSet, File outputFile)
        throws IOException;

    public void saveFeatureSet(FeatureSet featureSet, PrintWriter out);

    /**
     * Can this type of file handler handle this specific file?
     * Implementation is up to the handler, but this should be as low-cost as possible
     * @param file
     * @return
     * @throws IOException
     */
    public boolean canHandleFile(File file)
        throws IOException;
}
