/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.toolbox.normalize;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.RInterface;
import org.fhcrc.cpl.toolbox.TempFileManager;

import java.io.*;
import java.util.List;

/**
 * Apply Pei's normalization code to a peptide array
 */
public class Normalizer
{
    private static Logger _log = Logger.getLogger(Normalizer.class);

    private static final String normalizationScriptName = "normalize_array.R";
    private static final String intensityFileName = "ArrayIntensities.tsv";
    private static final String normalizedFileName = "NormalizedIntensities.tsv";

    /**
     * Returns true on success, false otherwise
     */
    public static boolean normalize(List<float[]> rows)
    {
        StringBuffer fakeTempFileCaller =
                new StringBuffer("fake StringBuffer to identify files that should be cleaned up for Normalizer");

        File normalizationScript;
        File intensityFile;
        File normalizedFile;

        try
        {
            InputStream in = Normalizer.class.getResourceAsStream(normalizationScriptName);
            normalizationScript = TempFileManager.createTempFile(normalizationScriptName, fakeTempFileCaller);
            OutputStream out = new FileOutputStream(normalizationScript);
            byte[] buf = new byte[1024];
            int len;
            while ((len = in.read(buf)) > 0)
                out.write(buf, 0, len);
            in.close();
            out.close();
        }
        catch (IOException e)
        {
            ApplicationContext.errorMessage("Error writing temp files", e);
            TempFileManager.deleteTempFiles(fakeTempFileCaller);
            return false;
        }

        try
        {
            intensityFile = TempFileManager.createTempFile(intensityFileName, fakeTempFileCaller);
            writeIntensityFile(rows, intensityFile);

            RInterface.runRScript(normalizationScript, fakeTempFileCaller);
            normalizedFile = TempFileManager.createTempFile(normalizedFileName, fakeTempFileCaller);
            readIntensityFile(rows, normalizedFile);
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("Error running normalization in R. Make sure R is on your path before launching msInspect.  R temporary files can be viewed in the directory " + TempFileManager.getTmpDir(), e);
            return false;
        }

        TempFileManager.deleteTempFiles(fakeTempFileCaller);

        return true;
    }

    /**
     * Write a tsv file containing *only* the intensities for the array
     */
    private static void writeIntensityFile(List rows, File intensityFile) throws IOException
    {
        PrintWriter pw = new PrintWriter(new FileOutputStream(intensityFile));

        float[] intensities = (float[])rows.get(1);

        for (int i = 1; i <= intensities.length; i++)
            pw.print((i > 1 ? "\t" : "" ) + "intensity" + i);
        pw.print("\n");

        for (int r = 0; r < rows.size(); r++)
        {
            intensities = (float[])rows.get(r);
            for (int i = 0; i < intensities.length; i++)
                pw.print((i > 0 ? "\t" : "" ) + intensities[i]);
            pw.print("\n");
        }
        pw.close();
    }

    /**
     * Read a tsv file containing the normalized intensities for the array.
     * Assumes that the normalization script wrote everything in the same
     * order as it read them in.
     */
    private static void readIntensityFile(List<float[]> rows, File normalizedFile)
            throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(normalizedFile));

        String header = br.readLine();

        for (int r = 0; r < rows.size(); r++)
        {
            String line = br.readLine();
            String[] tokens = line.split("\t");
            float[] intensities = rows.get(r);
            if (tokens.length != intensities.length)
                _log.error("Length mismatch on reading normalized intensities " + line);
            for (int i = 0; i < tokens.length; i++)
                intensities[i] = Float.parseFloat(tokens[i]);
        }
        br.close();
    }

}
