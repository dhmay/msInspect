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

package org.fhcrc.cpl.viewer.commandline;

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import java.io.File;
import java.io.IOException;

/**
 * Utilities related to commandline modules in the msInspect application
 */
public class ViewerCommandModuleUtilities
{
    /**
     * Given a feature file, find another feature file with a name like that one, in directory "directory"
     * @param baseFeatureFile
     * @param directory
     * @return
     * @throws IOException if not found or can't open
     */
    public static FeatureSet findCorrespondingFeatureSet(File baseFeatureFile, File directory)
            throws IOException
    {
        String baseName = baseFeatureFile.getName().substring(0,
                        baseFeatureFile.getName().indexOf("."));
        String pepXmlFilename =
                (baseName + ".pep.xml");
        String tsvFilename =
                (baseName + ".tsv");
        String regularTsvFilename =
                (baseName + ".peptides.tsv");
        String filteredTsvFilename =
                (baseName + ".filtered.tsv");

        File resultFile = null;
        for (String potentialFilename : directory.list())
        {
            if (potentialFilename.equalsIgnoreCase(pepXmlFilename) ||
                    potentialFilename.equalsIgnoreCase(tsvFilename) ||
                    potentialFilename.equalsIgnoreCase(regularTsvFilename) ||
                    potentialFilename.equalsIgnoreCase(filteredTsvFilename) )
            {
                resultFile = new File(directory.getAbsolutePath() + File.separatorChar +
                        potentialFilename);
                break;
            }
        }
        if (resultFile == null)
            throw new IOException("No corresponding feature file for file " + baseFeatureFile.getAbsolutePath());

        FeatureSet result = null;
        try
        {
            result = new FeatureSet(resultFile);
            ApplicationContext.setMessage("Located feature file " +
                    resultFile.getAbsolutePath() +
                    " with " + result.getFeatures().length + " features");
        }
        catch (Exception e)
        {
            throw new IOException("Failed to load feature file " + resultFile.getAbsolutePath() +
                    ": " + e.getMessage());
        }

        return result;
    }

    /**
     * Given a feature file, find its mzXML file in mzXmlDir
     * @param featureFile
     * @param mzXmlDir
     * @return
     * @throws IOException if not found or can't open
     */
    public static File findCorrespondingMzXmlFile(File featureFile, File mzXmlDir)
            throws IOException
    {
        String featureFilename = featureFile.getName();
        String mzXmlFileName =
                (featureFilename.substring(0, featureFilename.indexOf(".")) +
                        ".mzXML");
        File result = null;
        for (File potentialMzXmlFile : mzXmlDir.listFiles())
        {
            String potentialMzXmlFilename = potentialMzXmlFile.getName();
            if (potentialMzXmlFilename.equalsIgnoreCase(mzXmlFileName))
            {
                result = potentialMzXmlFile;
                ApplicationContext.setMessage("Located mzXML file " +
                        potentialMzXmlFile.getAbsolutePath());
            }
        }
        if (result == null)
            throw new IOException("No corresponding mzXML file for feature file " + featureFile.getAbsolutePath());
        return result;
    }
}
