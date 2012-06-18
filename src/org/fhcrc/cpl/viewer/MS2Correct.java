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

import org.apache.commons.collections.map.MultiValueMap;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

/**
 * User: migra
 * Date: Jun 7, 2005
 * Time: 3:31:49 PM
 */
public class MS2Correct
{
    public static void correct(File mzXMLFile, File ms1File, FeatureSet.FeatureSelector sel, File outFile) throws Exception
    {
        if (!mzXMLFile.exists())
        {
            ApplicationContext.errorMessage("Couldn't find file " + mzXMLFile.getPath(), null);
            return;
        }
        if (!ms1File.exists())
        {
            ApplicationContext.errorMessage("Couldn't find feature file" + ms1File.getPath(), null);
            return;
        }
        File ms2File = new File(mzXMLFile.getPath() + ".ms2.tsv");
        if (!ms2File.exists())
        {
            File inspectFile = new File(mzXMLFile.getPath() + ".inspect");
            if (inspectFile.exists())
                inspectFile.delete();
            ApplicationContext.setMessage("Reindexing to find ms2 scans");
        }

        PrintWriter out = null;
        try
        {
            out = new PrintWriter(new FileWriter(outFile));
            out.println("# ms1FeatureFile=" + ms1File.getPath());
            out.println("# minPeaks=" + sel.getMinPeaks());
            out.println("# minScans=" + sel.getMinScans());
            out.println("# maxKL=" + sel.getMaxKL());
            out.println("# maxCharge=" + sel.getMaxCharge());
            out.println("# minCharge=" + sel.getMinCharge());

            MSRun run = MSRun.load(mzXMLFile.getPath());
            FeatureSet fsMS2 = new FeatureSet(ms2File);
            Feature[] ms2Features = fsMS2.getFeatures();

            FeatureSet fs = new FeatureSet(ms1File);
            //sel.setMinPeaks(2);
            //sel.setMinCharge(1);
            //fs = fs.filter(sel);
            Feature[] features = fs.getFeatures(sel);
            correct(run, ms2Features, features, out);
        }
        finally
        {
            if (null != out)
                out.close();
        }
    }

    public static void correct(MSRun run, Feature[] ms2Features, Feature[] ms1Features, PrintWriter out)
    {
        Map assignMap = new MultiValueMap();
        Tree2D ms2Tree = new Tree2D();
        for (int i = 1; i < ms2Features.length; i++)
        {
            Feature feature = ms2Features[i];
            ms2Tree.add(feature.scan, feature.mz, feature);
        }

        Arrays.sort(ms1Features, new Feature.IntensityDescComparator());
        for (Feature feature : ms1Features)
        {
            int minScanIndex = run.getIndexForScanNum(feature.getScanFirst());
            if (minScanIndex > 0)
                minScanIndex--;
            int maxScanIndex = run.getIndexForScanNum(feature.getScanLast());
            if (maxScanIndex < run.getScanCount() - 1)
                maxScanIndex++;


            int minScan = run.getScan(minScanIndex).getNum();
            int maxScan = run.getScan(maxScanIndex).getNum() + 1;
            double minMz = feature.mz - .25 / (feature.charge == 0 ? 1 : feature.charge);
            double maxMz = feature.mz + (feature.peaks + .25)/ (feature.charge == 0 ? 1 : feature.charge);
            List<Feature> list = ms2Tree.getPoints(minScan, (float) minMz, maxScan, (float) maxMz);
            if (null != list)
                for (Feature ms2Feature : list)
                {
                    assignMap.put(ms2Feature, feature);
                }

        }

        out.println ("scan\tmz\tcharge\tintensity\tpeaks\tkl\tmatch");
        //Write them out in the same order they came in...
        for (Feature ms2Feature : ms2Features)
        {
            Collection<Feature> col = (Collection) assignMap.get(ms2Feature);
            if (null != col && col.size() > 0)
            {
                int match = 0;
                for (Feature feature : col)
                {
                    out.println (ms2Feature.scan + "\t" + feature.mz + "\t" + + feature.charge + "\t" + feature.intensity + "\t" + feature.peaks + "\t" + feature.kl + "\t" +  match);
                    match++;
                }
            }
        }
    }
}
