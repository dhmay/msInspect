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
package org.fhcrc.cpl.viewer;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.feature.FeatureStrategyCentroided;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;

import org.apache.log4j.Logger;

public class DumpWindow2D
{
    private static Logger _log = Logger.getLogger(DumpWindow2D.class);

    private static Comparator compareFeaturesByScanAndMzAsc = new Comparator(){
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature)o1;
            Feature f2 = (Feature)o2;
            if ( f1.scan != f2.scan )
                return f1.scan > f2.scan ? 1 : -1;
            return f1.mz == f2.scan ? 0 : ( f1.mz > f2.mz ? 1 : -1 );
        }
    };

    public static void dump(File runFile, File featureFile, int leadingScans, int trailingScans, int leadingPeaks, int trailingPeaks, File outFile) throws Exception
    {
        if (!runFile.exists())
        {
            ApplicationContext.errorMessage("Couldn't find mzXML file " + runFile.getPath(), null);
            return;
        }
        if (!featureFile.exists())
        {
            ApplicationContext.errorMessage("Couldn't find feature file" + featureFile.getPath(), null);
            return;
        }

        PrintWriter out = null;
        try
        {
            out = new PrintWriter(new FileWriter(outFile));

            MSRun run = MSRun.load(runFile.getPath());

            if (null == run)
            {
                ApplicationContext.errorMessage("Couldn't load mzXML file" + runFile.getPath(), null);
                return;
            }

            FeatureSet fs = new FeatureSet(featureFile);
            Feature[] features = fs.getFeatures();

            Arrays.sort(features, compareFeaturesByScanAndMzAsc);

            // Check this for LTQ-FT; may just come from mzXML header...
            boolean centroided = run.getHeaderInfo().getDataProcessing().getCentroided() == 1;

            _log.debug("Data is " + (centroided ? "" : "not ") + "centroided.");

            out.println("FeatureID\tLabel\tPeptide\tProtein\tCharge\tm/z\tMass\tScan\tScan\tTime\tWindow");

            for (int i = 0; i < features.length; i++)
            {
                // Skip charge zero features
                if ( features[i].charge == 0 )
                    continue;

                int nScans = leadingScans + trailingScans;
                MSRun.MSScan[] scans = new MSRun.MSScan[nScans];

                // search back for some number of MS1 scans
                // int startIndex = run.getIndexForFeature(features[i]) - 1;
                int startIndex = run.getIndexForScanNum(features[i].scan) - 1;

                // The acrylquant.R code currently has trouble if you are
                // within the last 20 scans; don't output such features
                if ( startIndex < 10 || run.getScanCount() - startIndex < 20 )
                {
                    _log.warn("Skipping feature (" + features[i].scan + ", " + features[i].mass + ", " + features[i].charge + ") " + startIndex + " " + run.getScanCount());
                    continue;
                }

                int s = leadingScans - 1;
                while ( s >= 0 && startIndex >= 0 )
                {
                    scans[s] = run.getScan(startIndex);
                    if ( scans[s].getMsLevel() == 1 )
                        s--;
                    startIndex--;
                }

                // now search forward for some number of MS1 scans
                int endIndex = run.getIndexForFeature(features[i]);
                s = leadingScans;
                while ( s < nScans && endIndex < run.getScanCount() )
                {
                    scans[s] = run.getScan(endIndex);
                    if ( scans[s].getMsLevel() == 1 )
                        s++;
                    endIndex++;
                }

                nScans = s; // get the number of scans actually read

                float startMz = features[i].mz - leadingPeaks;
                float endMz =   features[i].mz + trailingPeaks;

                // HACK: We expect the description to contain feature ID, label type, etc.
                String desc = features[i].getDescription();
                String[] tokens = desc.split(" ");
                desc = join("\t", tokens);

                for (s = 0; s < nScans; s++)
                {
                    float[][] spectrum = scans[s].getSpectrum();
                    if ( centroided )
                        spectrum = FeatureStrategyCentroided.cleanSpectrum(spectrum);
                    float[] signal = Spectrum.Resample(spectrum, new FloatRange(startMz, endMz), 36);
                    StringBuffer sb = new StringBuffer();

                    sb.append(desc + "\t");
                    sb.append(features[i].charge + "\t"); // charge
                    sb.append(features[i].mz + "\t"); // mz
                    sb.append(features[i].mass + "\t"); // mass
                    sb.append(features[i].scan + "\t"); // scan
                    sb.append(scans[s].getNum() + "\t"); // scan
                    sb.append(scans[s].getRetentionTime() + "\t"); // retention time
                    for (int m = 0; m < signal.length; m++)
                        sb.append(signal[m] + "\t");
                    out.println(sb.toString());
                }
            }
        }
        finally
        {
            if (null != out)
                out.close();
        }
    }

    private static String join(String delim, String[] s)
    {
        if ( null == s )
            return null;

        if ( s.length == 0 )
            return "";

        StringBuffer sb = new StringBuffer(s[0]);
        for (int i = 1; i < s.length; i++)
        {
            if ( null != delim )
                sb.append(delim);
            sb.append(s[i]);
        }
        return sb.toString();
    }

}
