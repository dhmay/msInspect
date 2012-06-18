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
package org.fhcrc.cpl.viewer.align;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.statistics.RInterface;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;

import java.io.*;

/**
 * Spline-based alignment algorithm
 *
 *
 */
public class SplineAligner extends Aligner
{
    private static Logger _log = Logger.getLogger(SplineAligner.class);

    // 20 is pretty darn high.  Lower values may
    //be more appropriate in most cases.  20 is for backward compatibility.
    public static final int DEFAULT_DEGREES_OF_FREEDOM = 20;


    //This controls, technically speaking, how squiggly the alignment gets.
    //Higher numbers -> more squiggly.
    protected int degreesOfFreedom = DEFAULT_DEGREES_OF_FREEDOM;

    protected boolean copiedRCode = false;


    public SplineAligner()
    {
    }

    /**
     * Cover method. In the absence of filename information, use a filename that's
     * likely to be unique
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public double[] alignPairs(Pair<Double,Double>[] pairs, int maxValueToWarp)
    {
        return alignPairs(pairs, maxValueToWarp, "length" + pairs.length);
    }
    
    /**
     * Generic method for creating a warping based on pairs of values.
     * Create temp files with names based on tempFileNameStart, in case the user
     * wants to peruse them
     * @param pairs
     * @param maxValueToWarp
     * @return
     */
    public double[] alignPairs(Pair<Double,Double>[] pairs, double maxValueToWarp,
                               String tempFileNameStart)
    {
        _log.debug("alignPairs.  Pairs: " + pairs.length + ", max to warp: " +
                   maxValueToWarp + ", DoF: " + degreesOfFreedom);
        PrintWriter pw = null;
        double[] result = null;
        if (pairs.length < 4)
            throw new RuntimeException("SplineAligner.alignPairs: at least 4 pairs are necessary for alignment, only " +
                    pairs.length + " were provided");
        File pairsFile = null;
        try
        {
            if (!copiedRCode)
            {
                copyRCode();
                copiedRCode = true;
            }
            String pairFileName = tempFileNameStart + ".pairs.tsv";
            pairsFile = writePairsFile(pairs, pairFileName);

            File outFile = TempFileManager.createTempFile(tempFileNameStart + ".align_out.tsv", this);

            File rScriptFile = TempFileManager.createTempFile(tempFileNameStart + ".align.R", this);

            pw = new PrintWriter(new FileOutputStream(rScriptFile));
            pw.println("source('warp_scans.R')");

            String line = "map_file_scans_onto_file('" +
                    pairFileName + "', '" +
                    outFile.getName() + "', maxScan=" + maxValueToWarp + ", df=" +
                    degreesOfFreedom + ");";
            pw.println(line);
            pw.flush();
            pw.close();

            RInterface.runRScript(rScriptFile, this);

            result = parseWarpingMap(outFile, (int) maxValueToWarp);

//            if (buildCharts)
//                createChart(pairs, result);
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Failed alignment: " + x.getMessage(), x);
            TempFileManager.unmarkFileForDeletion(pairsFile, this);
            TempFileManager.deleteTempFiles(this);
            return null;
        }
        finally
        {
            if (pw != null)
                pw.close();
        }

        return result;
    }

    /**
     * copy warp_scans.R to the temp directory
     * @throws IOException
     * @throws FileNotFoundException
     */
    protected void copyRCode() throws IOException, FileNotFoundException
    {
        InputStream in = SplineAligner.class.getResourceAsStream("warp_scans.R");
        File alignFile = TempFileManager.createTempFile("warp_scans.R", this);
        OutputStream out = new FileOutputStream(alignFile);
        byte[] buf = new byte[1024];
        int len;
        while ((len = in.read(buf)) > 0)
            out.write(buf, 0, len);
        in.close();
        out.close();
	    if (!alignFile.exists())
           throw new FileNotFoundException("Failed to copy R code to temp directory");
    }

    /**
     * Generic code for parsing a map of integers to doubles
     * @param mapFile
     * @return
     * @throws IOException
     * @throws FileNotFoundException
     */
    public double[] parseWarpingMap(File mapFile, int maxValueToWarp)
            throws IOException, FileNotFoundException
    {
        double[] warpingMap = new double[maxValueToWarp + 1];

        BufferedReader br =
                new BufferedReader(new InputStreamReader(new FileInputStream(mapFile)));
        while (true)
        {
            String line = br.readLine();
            if (line == null) break;
            if (line.startsWith("source_scan"))
                continue;
            String[] sourceAndDest = line.split("\t");
            warpingMap[Integer.parseInt(sourceAndDest[0])] =
                    Double.parseDouble(sourceAndDest[1]);
        }
        return warpingMap;
    }

    /**
     * Write a file containing pairs
     * @param fileName
     * @throws FileNotFoundException
     */
    protected File writePairsFile(Pair<Double,Double>[] pairs, String fileName)
            throws FileNotFoundException
    {
        File currentPairsFile =
                    TempFileManager.createTempFile(fileName,
                                                   this);

        PrintWriter currentPairsPW = new PrintWriter(currentPairsFile);
        currentPairsPW.println("source\tdest");
        for (Pair<Double,Double> pair : pairs)
            currentPairsPW.println(pair.first + "\t" + pair.second);

        currentPairsPW.flush();
        _log.debug("wrote feature pairs to file " + currentPairsFile.getAbsolutePath());
        return currentPairsFile;
    }

    public int getDegreesOfFreedom()
    {
        return degreesOfFreedom;
    }

    public void setDegreesOfFreedom(int degreesOfFreedom)
    {
        this.degreesOfFreedom = degreesOfFreedom;
    }

}
