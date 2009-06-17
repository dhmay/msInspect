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
package org.fhcrc.cpl.viewer.align;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureGrouper;

import java.awt.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * User: migra
 * Date: Mar 9, 2005
 * Time: 11:06:52 AM
 */
public class BucketedPeptideArray implements Runnable
{
    private static Logger _log = Logger.getLogger(BucketedPeptideArray.class);

    private java.util.List _sets;
    private FeatureSet.FeatureSelector _sel;
    boolean _align = true;
    boolean _normalize = false;
    int _deconvoluteScanWindow = 6;
    double _deconvoluteMassWindow = .2;
    int _scanBucket = 50;
    double _massBucket = .2;
    boolean _sumDeconvolutedIntensities = true;

    protected Aligner _aligner = null;

    //different bucket sizes used for optimization
    protected int[] _scanBuckets = {20, 30, 50, 75, 100, 150, 200, 300, 400};
    protected double[] _massBuckets = {.025, .05, .1, .15, .2};

    protected FeatureGrouper _featureGrouper;

    protected Aligner.FeaturePairSelector _featurePairSelector =
            Aligner.DEFAULT_FEATURE_PAIR_SELECTOR;

    private String _outFileName;
    protected boolean _shouldWriteUndeconvolutedDetails;

    int _conflictResolver = FeatureGrouper.DEFAULT_CONFLICT_RESOLVER;

    public BucketedPeptideArray(List<?> sets)
    {
        _sets = sets;
        _featureGrouper = new FeatureGrouper();
        _aligner = new SplineAligner();
    }

    public BucketedPeptideArray(List<?> sets, FeatureSet.FeatureSelector sel)
    {
        this(sets);
        _sel = sel;
    }

    public BucketedPeptideArray(List<?> sets, FeatureSet.FeatureSelector sel, int scanBucket, double massBucket,
                                String outFileName, boolean align)
    {
        this(sets,sel);
        _scanBucket = scanBucket;
        _massBucket = massBucket;
        _outFileName = outFileName;
        _align = align;
    }

    public void run()
    {
        run(false);
    }

    public void run(boolean optimize)
    {
        run(optimize, FeatureGrouper.BUCKET_EVALUATION_MODE_ONE_FROM_EACH);
    }

    public void run(boolean optimize, int optimizationMode)
    {
        run(optimize, optimizationMode, false);
    }

    public void run(boolean optimize, int optimizationMode, boolean showCharts)
    {




        String detailsFileName = null;

        PrintWriter out = null;
        if (_sel == null)
        {
            _log.debug("*************CREATING DUMMY SELECTOR");
            _sel = new FeatureSet.FeatureSelector();
        }

        try
        {
            ApplicationContext.setMessage("Loading");
            java.util.List<FeatureSet> sourceFeatureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < _sets.size(); i++)
            {
                FeatureSet fs;
                if (_sets.get(i) instanceof File)
                    fs = new FeatureSet((File) _sets.get(i), Color.RED);
                else if (_sets.get(i) instanceof FeatureSet)
                    fs = (FeatureSet) _sets.get(i);
                else if (_sets.get(i) instanceof String)
                    fs = new FeatureSet(new File((String) _sets.get(i)), Color.RED);
                else
                {
                    ApplicationContext.errorMessage("Couldn't load feature set due to bad object: " + _sets.get(i).toString(), null);
                    return;
                }

                sourceFeatureSets.add(fs);
            }

            //Filter
            ApplicationContext.setMessage("Filtering");
            List<FeatureSet> featureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < sourceFeatureSets.size(); i++)
            {
                FeatureSet fs = (sourceFeatureSets.get(i)).filter(_sel);
                featureSets.add(fs);
            }

            //Align
            ApplicationContext.setMessage("Aligning");
            if (_align)
            {
                _aligner.setFeaturePairSelector(_featurePairSelector);
                featureSets = _aligner.alignFeatureSets(featureSets, showCharts);
            }

            //Deconvolute
            ApplicationContext.setMessage("Deconvoluting");
            java.util.List<FeatureSet> deconvoluted = new ArrayList<FeatureSet>();
            for (FeatureSet fs : featureSets)
            {
                deconvoluted.add(fs.deconvolute(_deconvoluteScanWindow, _deconvoluteMassWindow,
                        _sumDeconvolutedIntensities));
            }



            _featureGrouper.setGroupByMass(true);
            _featureGrouper.setConflictResolver(_conflictResolver);
            for (int i = 0; i < deconvoluted.size(); i++)
            {
                FeatureSet fs = (FeatureSet) deconvoluted.get(i);
                _featureGrouper.addSet(fs);
            }

            //dhmay adding for in-line optimization
            if (optimize)
            {
                StringBuffer massBucketsToPrint = new StringBuffer();
                for (double massBucket : _massBuckets)
                    massBucketsToPrint.append(massBucket + ", ");
                StringBuffer scanBucketsToPrint = new StringBuffer();
                for (double scanBucket : _scanBuckets)
                    scanBucketsToPrint.append(scanBucket + ", ");
                ApplicationContext.setMessage("Optimizing: mass buckets " + massBucketsToPrint +
                        " scan buckets " + scanBucketsToPrint);

                Pair<Double, Integer> bestBuckets =
                        _featureGrouper.calculateBestBuckets(_massBuckets, _scanBuckets, optimizationMode);
                _massBucket = bestBuckets.first;
                _scanBucket = bestBuckets.second;
                ApplicationContext.setMessage("Using mass and scan buckets " +
                        _massBucket + " and " + _scanBucket);
            }

            _featureGrouper.split2D(_massBucket, _scanBucket);

            if (_outFileName != null)
                out = new PrintWriter(new FileOutputStream(_outFileName));
            else
                out = new PrintWriter(System.out);

            writeHeader(out);
            _featureGrouper.writePeptideArray(out, _normalize); // ???? check for successful write?
            out.flush();

            //Write details file if outfilename != null
            if (_outFileName != null)
            {
                out.close();
                detailsFileName = calcDetailsFilepath(_outFileName, false);
                out = new PrintWriter(new FileOutputStream(detailsFileName));
                writeHeader(out);
                //todo: figure out dynamically if we need to write out MS2 extrainfo                
                _featureGrouper.writeArrayDetails(out, false, true);
                out.flush();

                if (_shouldWriteUndeconvolutedDetails)
                {
                    out.close();
                    String undeconvolutedDetailsFilename = calcDetailsFilepath(_outFileName, true);

                    out = new PrintWriter(new FileOutputStream(undeconvolutedDetailsFilename));
                    writeHeader(out);
                    out.println("#Undeconvoluted features");
                    //todo: figure out dynamically if we need to write out MS2 extrainfo
                    _featureGrouper.writeArrayDetails(out, true, true);
                    out.flush();
                }
            }
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Error building peptide array", x);
        }
        finally
        {
            ApplicationContext.setMessage("Peptide array complete." + (_outFileName == null ? "" : " See " + _outFileName + " and " + detailsFileName));
            if (out != null && _outFileName == null)
            {
                out.close();
//System.err.println("Closing out");
            }
        }
    }

    /**
     * Figure out the filepath for the details file related to an array file
     * @param arrayFilepath
     * @param undeconvoluted
     * @return
     */
    public static String calcDetailsFilepath(String arrayFilepath, boolean undeconvoluted)
    {
        int dotPos = arrayFilepath.lastIndexOf('.');
        String detailsFileName = null;
        String newFilenamePart = (undeconvoluted ? ".details.nodecon" : ".details");
        if (dotPos < 0)
            detailsFileName = arrayFilepath + newFilenamePart + ".tsv";
        else
            detailsFileName = arrayFilepath.substring(0, dotPos) + newFilenamePart + arrayFilepath.substring(dotPos);
        return detailsFileName;
    }

    /**
     * Look up the tag for the given feature set. Tags a are supplied as key/value pairs,
     * where the key is the prefix of the last element of the file path (everything from
     * the first "." on is stripped off).
     */
    private static String lookupTag(FeatureSet fs, Map<String, String> tags)
    {
        if (null == fs || null == fs.getSourceFile())
            return null;

        String filename = fs.getSourceFile().getName();

        if (tags.containsKey(filename))
            return tags.get(filename);

        return null;
    }

    /**
     * Group feature sets by tag. 
     * Align within group to generate a representative FeatureSet for each.
     * Then align the representatives.
     */
    public void alignByGroup(Map<String,String> tags, boolean strict)
    {
        String detailsFileName = null;

        PrintWriter out = null;
        if (_sel == null)
            _sel = new FeatureSet.FeatureSelector();

        try
        {
            ApplicationContext.setMessage("Loading");
            List<FeatureSet> sourceFeatureSets = new ArrayList<FeatureSet>();
            for (int i = 0; i < _sets.size(); i++)
            {
                FeatureSet fs;
                if (_sets.get(i) instanceof File)
                    fs = new FeatureSet((File) _sets.get(i), Color.RED);
                else if (_sets.get(i) instanceof FeatureSet)
                    fs = (FeatureSet) _sets.get(i);
                else if (_sets.get(i) instanceof String)
                    fs = new FeatureSet(new File((String) _sets.get(i)), Color.RED);
                else
                {
                    ApplicationContext.errorMessage("Couldn't load feature set due to bad object: " + _sets.get(i).toString(), null);
                    return;
                }

                sourceFeatureSets.add(fs);
            }

            // ???? At this point, if there is only one tag, we could revert to regular alignment

            // Filter all source feature sets and add each to the appropriate tagged list
            ApplicationContext.setMessage("Filtering");

            HashMap<String, List<FeatureSet>> taggedFeatureSets = new HashMap<String, List<FeatureSet>>();
            ArrayList<String> tagList = new ArrayList<String>(); // Try to keep order of big alignment consistent
            for (int i = 0; i < sourceFeatureSets.size(); i++)
            {
                FeatureSet fs = ((FeatureSet) sourceFeatureSets.get(i)).filter(_sel);
                fs.setTag(lookupTag(fs, tags));
                if (!taggedFeatureSets.containsKey(fs.getTag()))
                    taggedFeatureSets.put(fs.getTag(), new ArrayList<FeatureSet>());
                taggedFeatureSets.get(fs.getTag()).add(fs);
                if (!tagList.contains(fs.getTag()))
                    tagList.add(fs.getTag());
            }

            // Loop through the tags and align all associated FeatureSets
            ArrayList<FeatureSet> filteredFeatureSets = new ArrayList<FeatureSet>();
            for (String tag : tagList)
            {
                List<FeatureSet> featureSets = taggedFeatureSets.get(tag);

                // If only one FeatureSet in this group just add it to the list for the next round
                if (featureSets.size() == 1)
                {
                    filteredFeatureSets.add(featureSets.get(0));
                }
                else
                {
                    //Align
                    ApplicationContext.setMessage("Aligning group " + tag);
                    _aligner.setFeaturePairSelector(_featurePairSelector);
                    featureSets = _aligner.alignFeatureSets(featureSets, false);
            
                    //Deconvolve
                    List deconvoluted = new ArrayList();
                    for (int i = 0; i < featureSets.size(); i++)
                    {
                        FeatureSet fs = featureSets.get(i);
                        deconvoluted.add(fs.deconvolute(_deconvoluteScanWindow, _deconvoluteMassWindow,
                                _sumDeconvolutedIntensities));
                    }
            
                    FeatureGrouper grouper = new FeatureGrouper();
            
                    grouper.setGroupByMass(true);
                    grouper.setConflictResolver(_conflictResolver);
                    for (int i = 0; i < deconvoluted.size(); i++)
                    {
                        FeatureSet fs = (FeatureSet) deconvoluted.get(i);
                        grouper.addSet(fs);
                    }
            
                    ApplicationContext.setMessage("Creating array " + tag);
                    grouper.split2D(_massBucket, _scanBucket);

                    // If strict, a feature must match across all runs to be included.
                    // Otherwise ~75% of the runs is good enough
                    int minAligned = strict ? featureSets.size() :(int) (featureSets.size() * 3.0/4.0 + 0.5);
                    for (FeatureSet fs : grouper.filterByGroupedAlignment(minAligned))
                        filteredFeatureSets.add(fs);
                }
            }

            // Now create a new BucketedPeptideArray to align the filtered feature sets.
            // Note that we've already applied any other filters during the within group
            // alignment, so there should be no need for a selector here.
            BucketedPeptideArray arr = new BucketedPeptideArray(filteredFeatureSets);
            arr.setScanBucket(_scanBucket);
            arr.setMassBucket(_massBucket);
            arr.setOutFileName(_outFileName);
            arr.setAlign(true);
            arr.setNormalize(_normalize);
            arr.setConflictResolver(_conflictResolver);
            arr.run();
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Error building peptide array", x);
        }
    }

    private void writeHeader(PrintWriter out)
    {
        String revision = (String) ApplicationContext.getProperty("REVISION");
        out.println("# Revision: " + revision);
        out.print("# Files:");
        for (int i = 0; i < _sets.size(); i++)
            out.print(" " + (i + 1) + ": " + getSetName(_sets.get(i)));
        out.println();
        out.print("# Params: " + _sel.toString());
        out.println(" --massWindow=" + _massBucket + " --scanWindow=" + _scanBucket + (_normalize ? " --normalize" : ""));
        // ???? If normalize is true, it would be good to write the individual scales here
    }

    private String getSetName(Object set)
    {
        if (set instanceof File)
            return ((File) set).getName();
        else if (set instanceof FeatureSet)
            return ((FeatureSet) set).getSourceFile() != null ? ((FeatureSet) set).getSourceFile().getName() : set.toString();
        else
            return set.toString();
    }

    public java.util.List<?> getFiles()
    {
        return _sets;
    }

    public void setFiles(java.util.List<?> files)
    {
        this._sets = files;
    }

    public FeatureSet.FeatureSelector get_sel()
    {
        return _sel;
    }

    public void set_sel(FeatureSet.FeatureSelector _sel)
    {
        this._sel = _sel;
    }

    public boolean isAlign()
    {
        return _align;
    }

    public void setAlign(boolean align)
    {
        this._align = align;
    }

    public boolean getNormalize()
    {
        return _normalize;
    }

    public void setNormalize(boolean normalize)
    {
        this._normalize = normalize;
    }

    public int getDeconvoluteScanWindow()
    {
        return _deconvoluteScanWindow;
    }

    public void setDeconvoluteScanWindow(int deconvoluteScanWindow)
    {
        this._deconvoluteScanWindow = deconvoluteScanWindow;
    }

    public double getDeconvoluteMassWindow()
    {
        return _deconvoluteMassWindow;
    }

    public void setDeconvoluteMassWindow(double deconvoluteMassWindow)
    {
        this._deconvoluteMassWindow = deconvoluteMassWindow;
    }

    public double getMassBucket()
    {
        return _massBucket;
    }

    public void setMassBucket(double massBucket)
    {
        this._massBucket = massBucket;
    }

    public int getScanBucket()
    {
        return _scanBucket;
    }

    public void setScanBucket(int scanBucket)
    {
        this._scanBucket = scanBucket;
    }

    public String getOutFileName()
    {
        return _outFileName;
    }

    public void setOutFileName(String outFileName)
    {
        this._outFileName = outFileName;
    }

    public int[] getScanBuckets()
    {
        return _scanBuckets;
    }

    public void setScanBuckets(int[] _scanBuckets)
    {
        this._scanBuckets = _scanBuckets;
    }

    public double[] getMassBuckets()
    {
        return _massBuckets;
    }

    public void setMassBuckets(double[] _massBuckets)
    {
        this._massBuckets = _massBuckets;
    }

    public int getConflictResolver()
    {
        return _conflictResolver;
    }

    public void setConflictResolver(int conflictResolver)
    {
        this._conflictResolver = conflictResolver;
    }

    public FeatureGrouper getFeatureGrouper()
    {
        return _featureGrouper;
    }

    public void setFeatureGrouper(FeatureGrouper _featureGrouper)
    {
        this._featureGrouper = _featureGrouper;
    }

    public Aligner.FeaturePairSelector getFeaturePairSelector()
    {
        return _featurePairSelector;
    }

    public void setFeaturePairSelector(Aligner.FeaturePairSelector featurePairSelector)
    {
        this._featurePairSelector = featurePairSelector;
    }

    public Aligner getAligner()
    {
        return _aligner;
    }

    public void setAligner(Aligner _aligner)
    {
        this._aligner = _aligner;
    }

    public boolean isSumDeconvolutedIntensities()
    {
        return _sumDeconvolutedIntensities;
    }

    public void setSumDeconvolutedIntensities(boolean sumDeconvolutedIntensities)
    {
        _sumDeconvolutedIntensities = sumDeconvolutedIntensities;
    }

    public boolean isShouldWriteUndeconvolutedDetails()
    {
        return _shouldWriteUndeconvolutedDetails;
    }

    public void setShouldWriteUndeconvolutedDetails(boolean _shouldWriteUndeconvolutedDetails)
    {
        this._shouldWriteUndeconvolutedDetails = _shouldWriteUndeconvolutedDetails;
    }
}
