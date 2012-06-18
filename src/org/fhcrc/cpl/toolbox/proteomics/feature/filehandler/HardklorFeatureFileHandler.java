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
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * File handler for Hardklor feature files.   To be specific, this handles the output format of kronik,
 * the post-processor for Hardklor output that stitches single-scan features into multi-scan features
 */
public class HardklorFeatureFileHandler extends BaseFeatureSetFileHandler
    implements FeatureSetFileHandler
{
    static Logger _log = Logger.getLogger(HardklorFeatureFileHandler.class);


    public static final String FILE_TYPE_NAME = "HARDKLOR";

    protected static HardklorFeatureFileHandler singletonInstance = null;

    public static HardklorFeatureFileHandler getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new HardklorFeatureFileHandler();
        return singletonInstance;
    }

    /**
     * Load a FeatureSet
     * @param file
     * @return
     * @throws IOException
     */
    public FeatureSet loadFeatureSet(File file)
            throws IOException
    {
        FeatureSet result = new FeatureSet();
        TabLoader loader = new TabLoader(file);
        TabLoader.ColumnDescriptor[] cols = loader.getColumns();

        if (cols.length > 0)
        {
            TabLoader.TabLoaderIterator iterator = loader.iterator();

            List<Feature>featureList = new ArrayList<Feature>();
            while (iterator.hasNext())
            {
                Map<String, Object> rowMap = (Map<String,Object>) iterator.next();
                int firstScan = (Integer) rowMap.get("First Scan");

                Float mz = ((Double) rowMap.get("Base Isotope Peak")).floatValue();
                Float intensity = ((Double) rowMap.get("Best Intensity")).floatValue();                

                Feature feature = new Feature(firstScan, mz, intensity);
                feature.setScanFirst(firstScan);
                feature.setScanLast((Integer) rowMap.get("Last Scan"));
                feature.setScanCount((Integer) rowMap.get("Num of Scans"));
                feature.setCharge((Integer) rowMap.get("Charge"));
                //ignore Hardklor mass because it's one less than base peak mass
                feature.updateMass();
                feature.setTotalIntensity(((Double) rowMap.get("Summed Intensity")).floatValue());
                //ignore Hardklor first and last time because we can't store those
                feature.setTime(((Double) rowMap.get("Best RTime")).floatValue());              
                
                featureList.add(feature);

            }

            result.setFeatures(featureList.toArray(new Feature[featureList.size()]));
        }

        return result;
    }

    public void saveFeatureSet(FeatureSet featureSet, File outFile)
            throws IOException
    {
       PrintWriter out = null;
        assert null != featureSet.getFeatures();
        try
        {
            out = new PrintWriter(new FileOutputStream(outFile));
            saveFeatureSet(featureSet, out);
        }
        catch(IOException e)
        {
            throw e;
        }
        finally
        {
            if (out != null)
                out.close();
        }
    }

    /**
     * Save a FeatureSet
     * @param featureSet
     * @param out
     */
    public void saveFeatureSet(FeatureSet featureSet, PrintWriter out)
    {
        assert null != featureSet.getFeatures();

        String header = "File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\t" +
                "Best Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications";
        out.println(header);
        out.flush();

        Feature[] features = featureSet.getFeatures();
        for (Feature feature : features)
        {
            out.println(featureSet.getSourceFile().getAbsolutePath() + "\t" +
                        feature.getScanFirst() + "\t" + feature.getScanLast() + "\t" + feature.getScanCount() + "\t" +
                        feature.getCharge() + "\t" + (feature.getMass() -1) + "\t" + feature.getMz() + "\t" +
                        feature.getIntensity() + "\t" + feature.getTotalIntensity() + "\t" +
//we don't store start and end time
                        feature.getTime() + "\t" + feature.getTime() + "\t" + feature.getTime() + "\t" +
                        0 + "\t" + "_");
        }
        out.flush();
    }

    // TODO: use StringUtils instead
    private static String join(float[] f, String delim)
    {
        if ( null == f || 0 >= f.length )
            return null;
        StringBuffer sb = new StringBuffer();
        sb.append("" + f[0]);
        for(int i = 1; i < f.length; i++)
        {
            if ( null != delim)
                sb.append(delim);
            sb.append("" + f[i]);
        }
        return sb.toString();
    }

    /**
     * Can this type of file handler handle this specific file?
     * Implementation is up to the handler, but this should be as low-cost as possible
     * @param file
     * @return
     * @throws IOException
     */
    public boolean canHandleFile(File file)
        throws IOException
    {
        if (isXMLFile(file))
            return false;

        //Check out the first line and look for columns typical of Kronik output
        FileInputStream fis = new FileInputStream(file);
        int headerLength = 250;
        if (fis.available() < headerLength)
            return false;

        byte[] headerChars = new byte[headerLength];
        if (fis.read(headerChars) < headerLength)
            return false;

        String header = new String(headerChars);

        //check for several things that should be in Kronik files and not in other types
        return (header.contains("File") && header.contains("First Scan") && header.contains("Best Correlation") &&
                header.contains("Monoisotopic Mass"));
    }
}
