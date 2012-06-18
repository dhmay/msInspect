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
package org.fhcrc.cpl.toolbox.proteomics;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithRPerspectivePlot;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * _Highly_ generic 2-dimensional clusterer.  Adapted from FeatureGrouper, which was specific
 * to clustering Features by mass and scan.
 *
 * Arrays of Clusterables are added to a tree structure, which divides them up first by
 * dimension1 and then by dimension2.  Each Clusterable is aware of its value in both dimensions
 * and (as far as Clusterer2D is concerned) nothing else.
 *
 * Any calling class will have to implement its own Clusterer class in order for clustering to have
 * any meaning.
 *
 * So far so good.  Here's the weird bit:  we need to be able to cluster Features by mass, and
 * the max bucket size may be in terms of Daltons or in terms of PPM.  If it's in terms of Daltons,
 * great.  If ppm, then the canSplit() method has some very funky behavior.  By default it behaves as
 * you might expect, though.
 */
public class Clusterer2D
{
    private static Logger _log = Logger.getLogger(Clusterer2D.class);

    private List<Clusterable[]> _clusterableArrays = new ArrayList<Clusterable[]>();
    public Node _root;
    BucketSummary[] summaries;

    //determines how to determine whether to split buckets in each dimension
    protected ClusterDimensionSplitCalculator _dimensionSplitCalculator =
            new DefaultClusterDimensionSplitCalculator();


    //keeps track of whether each dimension should really be treated as int values (not doubles).
    // For formatting of output only!
    protected boolean _dimension1IsInt = false;
    protected boolean _dimension2IsInt = false;

    public Clusterer2D()
    {
    }

    public Clusterer2D(Clusterable[] clusterableArray)
    {
        this();
        this.addSet(clusterableArray);
    }




    public Pair<Double, Double> calculateBestBuckets(double[] dimension1Buckets,
                                                     double[] dimension2Buckets, boolean showCharts)
    {
        double bestDimension1Bucket = 0;
        double bestDimension2Bucket = 0;

        //maximum number of perfect matches we could possibly have.  Used in order to stop searching
        //if we find that many.  Not likely to encounter this in real scenarios, maybe should remove?
        int maxPossiblePerfectMatches = Integer.MAX_VALUE;
        for (Clusterable[] clusterableArray : _clusterableArrays)
            if (clusterableArray.length < maxPossiblePerfectMatches)
                maxPossiblePerfectMatches = clusterableArray.length;
        
        double[][] perfectMatches = new double[dimension1Buckets.length][dimension2Buckets.length];
        int bestNumPerfectMatches = -1;

        //evaluate all combinations of mass and hydrophobicity buckets
        for (int iMass = 0; iMass < dimension1Buckets.length; iMass++)
        {
            for (int iElution = 0; iElution < dimension2Buckets.length;
                 iElution++)
            {
                double dimension2Bucketsize = dimension2Buckets[iElution];
                double dimension1Bucketsize = dimension1Buckets[iMass];
                if (dimension2Bucketsize <= 0 || dimension1Bucketsize <= 0)
                    continue;
                split2D(dimension1Bucketsize, dimension2Bucketsize);
                perfectMatches[iMass][iElution] = rowsWithOneFromEach();
                _log.debug("Evaluating bucket size: mass " + dimension1Bucketsize + ", elution " + dimension2Bucketsize);
                _log.debug("    Bucket score: " + perfectMatches[iMass][iElution]);
                if (perfectMatches[iMass][iElution] >
                    bestNumPerfectMatches)
                {
                    bestDimension1Bucket = dimension1Bucketsize;
                    bestDimension2Bucket = dimension2Bucketsize;
                    bestNumPerfectMatches = (int) perfectMatches[iMass][iElution];

                    //if we've got as many perfect matches as we can get, stop evaluating buckets
                    if (bestNumPerfectMatches == maxPossiblePerfectMatches)
                    {
                        _log.debug("Found best matches early, saved " + (dimension1Buckets.length-iMass-1) + " mass evals");
                        break;
                    }
                }
            }
        }

        if (showCharts)
        {
            //don't want to require R, so catching exceptions
            try
            {
                PanelWithRPerspectivePlot perspPlot = new PanelWithRPerspectivePlot();
                perspPlot.setName("Perfect Buckets");
                perspPlot.setAxisRVariableNames("dim1","dim2","perfect_buckets");                
                perspPlot.plot(dimension1Buckets, dimension2Buckets, perfectMatches);
                perspPlot.displayInTab();
            }
            catch (Exception e)
            {
                _log.debug("Error showing plot");
            }
        }
        _log.debug("Best Mass/Elution buckets: " + bestDimension1Bucket + ", " +  bestDimension2Bucket);
        ApplicationContext.setMessage("Num perfect buckets: " + bestNumPerfectMatches);

        return new Pair<Double, Double>(bestDimension1Bucket, bestDimension2Bucket);
    }



    /**
     * A single entry in the clustering tree, aware of its value in only one dimension
     */
    public static class TreeEntry implements Comparable
    {
        public Clusterable parent;
        public int iSet;
        public double value;

        public int compareTo(Object o)
        {
            TreeEntry treeEntry = (TreeEntry) o;
            if (value < treeEntry.value)
                return 0;
            else if (value > treeEntry.value)
                return 1;
            return 0;
        }
    }

    public void setDimensionSplitCalculator(ClusterDimensionSplitCalculator dimensionSplitCalculator)
    {
        this._dimensionSplitCalculator = dimensionSplitCalculator;
    }

    /**
     * abstract class for determining whether to split in each dimension.  See usage in canSplit()
     */
    public abstract static class ClusterDimensionSplitCalculator
    {
        public abstract double calculateDimension1ForSplit(double referenceValue, double dimensionValue);
        public abstract double calculateDimension2ForSplit(double referenceValue, double dimensionValue);
    }

    /**
     * Default split calculator is very simple!
     */
    public static class DefaultClusterDimensionSplitCalculator
        extends ClusterDimensionSplitCalculator
    {
        public double calculateDimension1ForSplit(double referenceValue, double dimensionValue)
        {
            return dimensionValue;
        }
        public double calculateDimension2ForSplit(double referenceValue, double dimensionValue)
        {
            return dimensionValue;
        }
    }

    /**
     * defines values in two dimensions.  Individual implementations will likely have some other
     * object that they're aware of that they grab these values from
     */
    public interface Clusterable
    {
        public abstract double getDimension1Value();
        public abstract double getDimension2Value();
    }



    /**
     * A node in the clustering tree, aware of its children
     */
    private class Node
    {
        TreeEntry[] entries;
        int start;
        int length;
        Node left, right;
        Node dimension2; //used to divide up in a different dimension

        Node(TreeEntry[] entries, int start, int length)
        {
            assert length > 0;
            this.entries = entries;
            this.start = start;
            this.length = length;
        }

        double getMin()
        {
            return entries[start].value;
        }

        double getMax()
        {
            return entries[start + length - 1].value;
        }

        TreeEntry[] getEntries()
        {
            TreeEntry[] newArray = new TreeEntry[length];
            System.arraycopy(entries, start, newArray, 0, length);
            return newArray;
        }

        void addEntries(TreeEntry[] addEntries)
        {
            TreeEntry[] newEntries = new TreeEntry[length + addEntries.length];
            if (null != this.entries)
                System.arraycopy(this.entries, start, newEntries, 0, length);

            System.arraycopy(addEntries, 0, newEntries, length, addEntries.length);
            this.entries = newEntries;

            start = 0;
            length = this.entries.length;
            Arrays.sort(this.entries);
        }

        /**
         * In the case of splitting Features by mass, we might have to convert
         * from PPM to Da before figuring out whether we can split.  This means we need to
         * take into account getMin() as well as maxBucket.  So I've abstracted the calculation of
         * the allowed bucket width.
         * @param maxBucket
         * @param dimension the dimension we're considering -- can have different behavior for both
         * @return
         */
        boolean canSplit(double maxBucket, int dimension)
        {
            return (length > 1) &&
                   (getMax() - getMin() >
                    (dimension == 2 ? _dimensionSplitCalculator.calculateDimension2ForSplit(getMin(), maxBucket)
                                    : _dimensionSplitCalculator.calculateDimension1ForSplit(getMin(), maxBucket)));
        }

        void reSplit(double maxBucket, int dimension)
        {

            if (canSplit(maxBucket, dimension))
            {
                double maxDistance = 0;
                int bestIndex = -1;
                //TODO: Can pre-compute all of this...
                for (int i = start; i < start + length - 1; i++)
                {
                    double distance = entries[i + 1].value - entries[i].value;
                    if (distance > maxDistance)
                    {
                        maxDistance = distance;
                        bestIndex = i;
                    }
                }

                assert bestIndex >= start;
                left = new Node(entries, start, bestIndex + 1 - start);
                assert start + length - bestIndex - 1 > 0;
                right = new Node(entries, bestIndex + 1, start + length - bestIndex - 1);
                left.reSplit(maxBucket, dimension);
                right.reSplit(maxBucket, dimension);
            }
            else
                left = right = null;
        }

        List<Node> appendLeafNodes(List<Node> nodes)
        {
            if (nodes == null)
                nodes = new ArrayList<Node>();

            if (null == left)
                nodes.add(this);
            else
            {
                left.appendLeafNodes(nodes);
                right.appendLeafNodes(nodes);
            }

            return nodes;
        }

        List<Node> getLeafNodes()
        {
            return appendLeafNodes(null);
        }
    }

    public void addSet(Clusterable[] clusterableArray)
    {
        summaries = null;
        int index = _clusterableArrays.size();
        _clusterableArrays.add(clusterableArray);
        TreeEntry[] entries = new TreeEntry[clusterableArray.length];

        for (int i = 0; i < entries.length; i++)
        {
            Clusterable parent = clusterableArray[i];
            TreeEntry treeEntry = new TreeEntry();
            treeEntry.value = parent.getDimension1Value();
            treeEntry.parent = parent;
            treeEntry.iSet = index;
            entries[i] = treeEntry;
        }

        Arrays.sort(entries);
        if (null == _root)
            _root = new Node(entries, 0, entries.length);
        else
            _root.addEntries(entries);
   }

    public void split2D(double maxDimension1Bucket, double maxDimension2Bucket)
    {
        summaries = null;
        _root.reSplit(maxDimension1Bucket, 1);
        List<Node> leaves = _root.getLeafNodes();
        for (Node node : leaves)
        {
            TreeEntry[] massEntries = node.getEntries();
            TreeEntry[] hydrophobicityEntries = new TreeEntry[massEntries.length];
            for (int j = 0; j < massEntries.length; j++)
            {
                TreeEntry dimension1TreeEntry = massEntries[j];
                TreeEntry hydrophobicityTreeEntry = new TreeEntry();
                hydrophobicityTreeEntry.parent = dimension1TreeEntry.parent;
                hydrophobicityTreeEntry.iSet = dimension1TreeEntry.iSet;
                hydrophobicityTreeEntry.value = dimension1TreeEntry.parent.getDimension2Value();

                hydrophobicityEntries[j] = hydrophobicityTreeEntry;
            }
            Arrays.sort(hydrophobicityEntries);
            node.dimension2 = new Node(hydrophobicityEntries, 0, hydrophobicityEntries.length);
            node.dimension2.reSplit(maxDimension2Bucket, 2);
        }
    }


    public BucketSummary[] summarize()
    {
        if (summaries != null)
            return summaries;

        List<BucketSummary> summaryList = new ArrayList<BucketSummary>();
        List<Node> dimension1Leaves = _root.appendLeafNodes(null);
        for (Node dimension1Leaf : dimension1Leaves)
        {
            List<Node> dimension2Leaves = dimension1Leaf.dimension2.appendLeafNodes(null);
            for (Node dimension2Leaf : dimension2Leaves)
            {
                BucketSummary summary = new BucketSummary(dimension1Leaf, dimension2Leaf);
                summaryList.add(summary);
            }
        }

        summaries = summaryList.toArray(new BucketSummary[summaryList.size()]);
        return summaries;
    }

    public int[] histogramBucketCounts()
    {
        BucketSummary[] bucketSummaries = summarize();
        int maxCount = 0;
        for (int i = 0; i < bucketSummaries.length; i++)
            if (bucketSummaries[i].featureCount > maxCount)
                maxCount = bucketSummaries[i].featureCount;

        int[] histogram = new int[maxCount + 1];
        for (int i = 0; i < bucketSummaries.length; i++)
            histogram[bucketSummaries[i].featureCount]++;

        return histogram;
    }

    public int[] histogramSetCounts()
    {
        BucketSummary[] bucketSummaries = summarize();
        int maxCount = 0;
        for (int i = 0; i < bucketSummaries.length; i++)
            if (bucketSummaries[i].setCount > maxCount)
                maxCount = bucketSummaries[i].setCount;

        int[] histogram = new int[maxCount + 1];
        for (int i = 0; i < bucketSummaries.length; i++)
            histogram[bucketSummaries[i].setCount]++;

        return histogram;
    }

    public int numBuckets()
    {
        BucketSummary[] bucketSummaries = summarize();
        return bucketSummaries.length;
    }

    public int rowsWithOneFromEach()
    {
        BucketSummary[] bucketSummaries = summarize();
        return rowsWithOneFromEach(bucketSummaries);
    }

    /**
     * Calculates how many buckets have exactly 1 feature from each feature set.
     * Useful for testing align
     *
     * @return
     */
    public int rowsWithOneFromEach(BucketSummary[] bucketSummaries)
    {
        int count = 0;
        int numSets = _clusterableArrays.size();

        for (int i = 0; i < bucketSummaries.length; i++)
            if (bucketSummaries[i].featureCount == numSets && bucketSummaries[i].setCount == numSets)
                count++;

//TODO: Figure out what to do about this.  It would be nice to have the minimum number of
//TODO: featuresets with one from each configurable.        
//        for (int i = 0; i < bucketSummaries.length; i++)
//            if (bucketSummaries[i].featureCount > 50 &&
//                    bucketSummaries[i].featureCount == bucketSummaries[i].setCount)
//                count++;


        return count;
    }

    public class BucketSummary
    {
        private Node _dimension2Leaf; //Innermost node

        public double minDimension1;
        public double maxDimension1;
        public double minDimension2;
        public double maxDimension2;
        public int featureCount = 0;
        public int setCount = 0;
        public Map<Integer, List<TreeEntry>> setIndexTreeEntryListMap =
                new HashMap<Integer, List<TreeEntry>>();


        public BucketSummary(Node dimension1Leaf, Node dimension2Leaf)
        {
            this._dimension2Leaf = dimension2Leaf;
            this.minDimension1 = dimension1Leaf.getMin();
            this.maxDimension1 = dimension1Leaf.getMax();
            this.minDimension2 = dimension2Leaf.getMin();
            this.maxDimension2 = dimension2Leaf.getMax();
            this.featureCount = dimension2Leaf.length;
            int[] setCounts = new int[_clusterableArrays.size()];

            for (int i = dimension2Leaf.start; i < dimension2Leaf.start + dimension2Leaf.length; i++)
            {
                TreeEntry treeEntry = dimension2Leaf.entries[i];
                setCounts[treeEntry.iSet]++;
                getTreeEntryListForSetIndex(treeEntry.iSet).add(treeEntry);
            }

            for (int i = 0; i < _clusterableArrays.size(); i++)
                if (setCounts[i] > 0)
                {
                    setCount++;
                }
        }

        /**
         * No time is wasted here:  If there's only one entry, return it
         * @param valueToMatch
         * @param dimension
         * @param treeEntryList
         * @return
         */
        protected TreeEntry pickClosestTreeEntry(double valueToMatch, int dimension,
                                                 List<TreeEntry> treeEntryList)
        {
            TreeEntry bestTreeEntry = treeEntryList.get(0);

            if (treeEntryList.size() > 1)
            {
                double bestTreeEntryDifference =
                        Math.abs((dimension == 1 ? bestTreeEntry.parent.getDimension1Value() :
                                bestTreeEntry.parent.getDimension2Value()) -
                                valueToMatch);
                for (int i = 1; i < treeEntryList.size(); i++)
                {
                    TreeEntry treeEntry = treeEntryList.get(i);
//System.err.println("Possibility: " + (dimension == 1 ? treeEntry.parent.getDimension1Value() : treeEntry.parent.getDimension2Value()));

                    double treeEntryDifference =
                            Math.abs((dimension == 1 ? bestTreeEntry.parent.getDimension1Value() :
                                    bestTreeEntry.parent.getDimension2Value()) -
                                    valueToMatch);
                    if (treeEntryDifference < bestTreeEntryDifference)
                    {
                        bestTreeEntry = treeEntry;
                        bestTreeEntryDifference = treeEntryDifference;
                    }
                }
            }
            return bestTreeEntry;
        }

    /**
     * Resolve pairing conflicts for features in the same bucket.  If there's more than one
     * feature from each set in a bucket, it gets hairy.
     * Time is wasted if there's already just one for each set, because I don't want to
     * waste time checking for that here.  If that's the case, don't call this method.
     *
     * The way one feature from each set is picked is as follows:
     *  -Along the interestingDimension, find the average value in the bucket
     *  -For each set, find each feature's deviance from that value and pick the feature with the smallest
     *
     * TODO: The picking in this method represents the most questionable choices in all the clustering code.
     * We need to revisit:
     *  -is it reasonable to do this picking based on just one dimension?
     *  -is it reasonable to use the center of the bucket in order to pick the "closest" feature?
     *      -if so, when calculating the center of the bucket, should we exclude features from the current
     * set (this would take more time)
     *
     * @return A list of Clusterables, one representing each set, in order of increasing
     * set index
     */
    public Clusterable[] pickOneFromEachSet(int interestingDimension)
    {
        double interestingDimensionSum = 0;

        for (TreeEntry treeEntry : _dimension2Leaf.entries)
        {
            interestingDimensionSum +=
                    (interestingDimension == 1? treeEntry.parent.getDimension1Value() :
                                                treeEntry.parent.getDimension2Value());
        }
        double interestingDimensionMean =
                interestingDimensionSum / (double) _dimension2Leaf.entries.length;
//System.err.println("mean: " + interestingDimensionMean);

        Map<Integer, List<TreeEntry>> interestingDimensionSetIndexTreeEntryListMap =
                (interestingDimension == 2 ? setIndexTreeEntryListMap :
                                             setIndexTreeEntryListMap);


        Clusterable[] result = new Clusterable[_clusterableArrays.size()];

        for (int i=0; i<_clusterableArrays.size(); i++)
        {
            List<TreeEntry> setTreeEntryList =
                    interestingDimensionSetIndexTreeEntryListMap.get(i);
            if (setTreeEntryList == null || setTreeEntryList.size() == 0)
                result[i] = null;
            else
            {
                result[i] = pickClosestTreeEntry(interestingDimensionMean,
                                                 interestingDimension, setTreeEntryList).parent;
            }
        }

        return result;
    }

        public List<Clusterable> getParentList()
        {
            List<Clusterable> result = new ArrayList<Clusterable>();
            //dhmay fixing 20100512.  This used to iterate through _dimension2Leaf.entries,
            //which is <sigh> completely different
            for (TreeEntry treeEntry : _dimension2Leaf.getEntries())
            {
                result.add(treeEntry.parent);
            }
            return result;
        }

        public List<Clusterable> getParentListForSetIndex(int setIndex)
        {
            //Dimension doesn't matter, because both tree lists will have the same parents
            List<TreeEntry> treeEntryList = getTreeEntryListForSetIndex(setIndex);
            List<Clusterable> result = new ArrayList<Clusterable>(treeEntryList.size());
            for (TreeEntry treeEntry : treeEntryList)
            {
                result.add(treeEntry.parent);
            }
            return result;
        }

        public List<TreeEntry> getTreeEntryListForSetIndex(int setIndex)
        {
            List<TreeEntry> result = setIndexTreeEntryListMap.get(setIndex);
            if (result == null)
            {
                result = new ArrayList<TreeEntry>();
                setIndexTreeEntryListMap.put(setIndex,result);
            }
            return result;
        }

        public Object[] objects()
        {
            Object[] objects = new Feature[_dimension2Leaf.length];
            for (int i = 0; i < _dimension2Leaf.length; i++)
            {
                TreeEntry treeEntry = _dimension2Leaf.entries[i + _dimension2Leaf.start];
                objects[i] = treeEntry.parent;
            }

            return objects;
        }

        public TreeEntry[] entries()
        {
            TreeEntry[] entries = new TreeEntry[_dimension2Leaf.length];
            for (int i = 0; i < _dimension2Leaf.length; i++)
            {
                entries[i] = _dimension2Leaf.entries[i + _dimension2Leaf.start];
            }

            return entries;
        }

        public String arrayRowDetail()
        {
            StringBuffer sb = new StringBuffer(arrayRow());
            for (int i = _dimension2Leaf.start; i < _dimension2Leaf.start + _dimension2Leaf.length; i++)
            {
                TreeEntry treeEntry = _dimension2Leaf.entries[i];
                sb.append("\t");
                sb.append("" + i);
                sb.append("\t");
                sb.append(treeEntry.parent.toString());
            }

            return sb.toString();
        }

        public String arrayRow()
        {
            StringBuffer sb = new StringBuffer(toString());

            return sb.toString();
        }

        public String toString()
        {
            StringBuffer dimension1String = new StringBuffer(minDimension1 + "\t" + maxDimension1);
            StringBuffer dimension2String = new StringBuffer(minDimension2 + "\t" + maxDimension2);

            if (dimension1IsInt())
                dimension1String = new StringBuffer((int) minDimension1 + "\t" + (int) maxDimension1);
            if (dimension2IsInt())
                dimension2String = new StringBuffer((int) minDimension2 + "\t" + (int) maxDimension2);

            return dimension1String + "\t" + dimension2String + "\t" + featureCount + "\t" + setCount;
        }
    }


    //getters and setters

    public boolean dimension1IsInt()
    {
        return _dimension1IsInt;
    }

    public void setDimension1IsInt(boolean _dimension1IsInt)
    {
        this._dimension1IsInt = _dimension1IsInt;
    }

    public boolean dimension2IsInt()
    {
        return _dimension2IsInt;
    }

    public void setDimension2IsInt(boolean _dimension2IsInt)
    {
        this._dimension2IsInt = _dimension2IsInt;
    }

    public List<Clusterable[]> getClusterableArrays()
    {
        return _clusterableArrays;
    }

    public void setClusterableArrays(List<Clusterable[]> clusterableArrays)
    {
        _clusterableArrays = clusterableArrays;
    }

    public int countAllEntries()
    {
        if (_root == null)
            return 0;
        return _root.entries.length;
    }
}
