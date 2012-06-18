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
package org.fhcrc.cpl.toolbox.datastructure;

import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Tree node that's aware of spatial information and has up to four children,
 * one in each quadrant.
 *
 * Each node keeps track of its theoretical boundaries as well as the bounding
 * box around the actual points in the node (which may be smaller)
 *
 * Only leaf nodes will have data values
 */
public class SpatialQuadrantTreeNode<T extends List<XYDataPoint>>
        extends ParentAwareTreeNode<T>
{
    protected static Logger _log = Logger.getLogger(SpatialQuadrantTreeNode.class);

    //values for scaffolding printlns
    protected static float minXScaff = 450f;
    protected static float maxXScaff = 500f;
    protected static float minYScaff = -600f;
    protected static float maxYScaff = -500f;
    protected static boolean isInScaffRange(double x, double y)
    { return (x>=minXScaff && x<=maxXScaff && y>=minYScaff && y<=maxYScaff); }

    protected double minDataX = Double.MAX_VALUE;
    protected double maxDataX = Double.MIN_VALUE;
    protected double minDataY = Double.MAX_VALUE;
    protected double maxDataY = Double.MIN_VALUE;

    protected double minBoundaryX = Double.MAX_VALUE;
    protected double maxBoundaryX = Double.MIN_VALUE;
    protected double minBoundaryY = Double.MAX_VALUE;
    protected double maxBoundaryY = Double.MIN_VALUE;

    protected SpatialQuadrantTreeNode<T> upperLeftChild = null;
    protected SpatialQuadrantTreeNode<T> upperRightChild = null;
    protected SpatialQuadrantTreeNode<T> lowerLeftChild = null;
    protected SpatialQuadrantTreeNode<T> lowerRightChild = null;

    protected int depth = 0;

    public static final DataPointXComparatorAsc dataPointXComparator =
            new DataPointXComparatorAsc();

    public SpatialQuadrantTreeNode(double minBoundaryX, double maxBoundaryX,
                                   double minBoundaryY, double maxBoundaryY,
                                   double minDataX, double maxDataX,
                                   double minDataY, double maxDataY)
    {
        this.minBoundaryX = minBoundaryX;
        this.maxBoundaryX = maxBoundaryX;
        this.minBoundaryY = minBoundaryY;
        this.maxBoundaryY = maxBoundaryY;

        this.minDataX = minDataX;
        this.maxDataX = maxDataX;
        this.minDataY = minDataY;
        this.maxDataY = maxDataY;
    }


    public List<SpatialQuadrantTreeNode<T>> findAllLeaves()
    {
        List<SpatialQuadrantTreeNode<T>> result =  new ArrayList<SpatialQuadrantTreeNode<T>>();
        if (isLeafNode())
        {
            result.add(this);
        }
        else
        {
            for (TreeNode child : getChildren())
            {
                result.addAll(((SpatialQuadrantTreeNode<T>) child).findAllLeaves());
            }
        }
        return result;
    }

    /**
     * Recursively build the quadrant tree
     * @param dataPoints
     * @return
     */
    public static SpatialQuadrantTreeNode buildTree(
            double minBoundaryX, double maxBoundaryX,
            double minBoundaryY, double maxBoundaryY,
            List<XYDataPoint> dataPoints, int currentDepth, int maxDepth)
    {
//System.err.println(minBoundaryX + ", " + maxBoundaryX  + ", " +minBoundaryY + ", " + maxBoundaryY);
        Collections.sort(dataPoints, dataPointXComparator);

        int numPoints = dataPoints.size();

        double minDataX = dataPoints.get(0).getX();
        double maxDataX = dataPoints.get(numPoints-1).getX();
//System.err.println("minDataX=" + minDataX + ", maxDataX=" + maxDataX);
        double minDataY = Double.MAX_VALUE;
        double maxDataY = Double.MIN_VALUE;
        for (XYDataPoint point : dataPoints)
        {
            if (point.getY() < minDataY)
                minDataY = point.getY();
            if (point.getY() > maxDataY)
                maxDataY = point.getY();
        }

        SpatialQuadrantTreeNode<List<XYDataPoint>> newNode =
                new SpatialQuadrantTreeNode<List<XYDataPoint>>(
                        minBoundaryX, maxBoundaryX, minBoundaryY, maxBoundaryY,
                        minDataX, maxDataX, minDataY, maxDataY);
        newNode.setDepth(0);

        if (dataPoints.size() == 1 || currentDepth >= maxDepth)
        {
            newNode.setData(dataPoints);
            return newNode;
        }

        double middleX = BasicStatistics.mean(new double[] {maxBoundaryX,minBoundaryX});
        double middleY = BasicStatistics.mean(new double[] {maxBoundaryY,minBoundaryY});

        List<XYDataPoint> lowerLeftPoints = new ArrayList<XYDataPoint>(0);
        List<XYDataPoint> lowerRightPoints = new ArrayList<XYDataPoint>(0);
        List<XYDataPoint> upperLeftPoints = new ArrayList<XYDataPoint>(0);
        List<XYDataPoint> upperRightPoints = new ArrayList<XYDataPoint>(0);

        boolean rightHalf = false;
        for (XYDataPoint point : dataPoints)
        {
            if (!rightHalf && point.getX() > middleX)
                rightHalf = true;
            if (point.getY() < middleY)
            {
                if (rightHalf)
                    lowerRightPoints.add(point);
                else
                    lowerLeftPoints.add(point);
            }
            else
            {
                {
                    if (rightHalf)
                        upperRightPoints.add(point);
                    else
                        upperLeftPoints.add(point);
                }
            }
        }

        if (lowerLeftPoints.size() > 0)
        {
            newNode.setLowerLeftChild(buildTree(
                    minBoundaryX, middleX, minBoundaryY, middleY,
                    lowerLeftPoints, currentDepth+1, maxDepth));
        }
        if (lowerRightPoints.size() > 0)
        {
            newNode.setLowerRightChild(buildTree(
                    middleX, maxBoundaryX, minBoundaryY, middleY,
                    lowerRightPoints, currentDepth+1, maxDepth));
        }
        if (upperLeftPoints.size() > 0)
        {
            newNode.setUpperLeftChild(buildTree(
                    minBoundaryX, middleX, middleY, maxBoundaryY,
                    upperLeftPoints, currentDepth+1, maxDepth));
        }
        if (upperRightPoints.size() > 0)
        {
            newNode.setUpperRightChild(buildTree(
                    middleX, maxBoundaryX, middleY, maxBoundaryY,
                    upperRightPoints, currentDepth+1, maxDepth));
        }

        return newNode;
    }


    /**
     * Search the tree to find the node that contains the specified point.  May or may not be a leaf.
     * @param treeRoot MUST contain the originPoint.  This is an assumption.
     * @param originPoint
     * @return
     */
    public static SpatialQuadrantTreeNode findNodeContainingPoint(
            SpatialQuadrantTreeNode<List<XYDataPoint>> treeRoot, XYDataPoint originPoint)
    {
        SpatialQuadrantTreeNode<List<XYDataPoint>> appropriateChild = null;
        boolean isRightOfMiddle = originPoint.getX() >= treeRoot.getMiddleXByBoundaries();
        boolean isUpOfMiddle = originPoint.getY() >= treeRoot.getMiddleYByBoundaries();
        if (isRightOfMiddle)
        {
            if (isUpOfMiddle)
                appropriateChild = treeRoot.getUpperRightChild();
            else
                appropriateChild = treeRoot.getLowerRightChild();
        }
        else
        {
            if (isUpOfMiddle)
                appropriateChild = treeRoot.getUpperLeftChild();
            else
                appropriateChild = treeRoot.getLowerLeftChild();
        }

        if (appropriateChild == null)
            return treeRoot;
        if (appropriateChild.getMinDataX() <= originPoint.getX() &&
                appropriateChild.getMaxDataX() >= originPoint.getX() &&
                appropriateChild.getMinDataY() <= originPoint.getY() &&
                appropriateChild.getMaxDataY() >= originPoint.getY())
            return findNodeContainingPoint(appropriateChild, originPoint);
        return treeRoot;
    }

    public static List<XYDataPoint> findClosestDataPoints(
            SpatialQuadrantTreeNode<List<XYDataPoint>> treeRoot,
            XYDataPoint originPoint,
            int numPointsToFind)
    {
         return findClosestDataPoints(treeRoot, originPoint, null,
                 new ArrayList<XYDataPoint>(), 0f, numPointsToFind);
    }

    public static List<XYDataPoint> findClosestDataPoints(
            SpatialQuadrantTreeNode<List<XYDataPoint>> treeRoot,
            XYDataPoint originPoint,
            SpatialQuadrantTreeNode<List<XYDataPoint>> containingNode,
            int numPointsToFind)
    {
         return findClosestDataPoints(treeRoot, originPoint, containingNode,
                 new ArrayList<XYDataPoint>(), 0f, numPointsToFind);
    }


    /**
     * Find the closest nodes to a given starting point (which may not have a node in the tree)
     * @param treeRoot
     * @param originPoint
     * @param numPointsToFind
     * @return a list of nodes, sorted by ascending distance
     */
    public static List<XYDataPoint> findClosestDataPoints(
            SpatialQuadrantTreeNode<List<XYDataPoint>> treeRoot,
            XYDataPoint originPoint,
            SpatialQuadrantTreeNode<List<XYDataPoint>> containingNode,
            List<XYDataPoint> seedPoints,
            float farthestSeedPointDistance,
            int numPointsToFind)
    {
        if (containingNode == null)
            containingNode = findNodeContainingPoint(treeRoot, originPoint);

        float farthestDistanceSoFar = 0f;
        if (seedPoints.size() > 0)
           farthestDistanceSoFar = farthestSeedPointDistance;

//        List<SpatialQuadrantTreeNode<XYDataPoint>> allLeafNodesInContainer =
//                containerNode.findAllLeaves();

//if (!containerNode.isLeafNode())
//    System.err.println("Node containing " + originPoint + " is not a leaf node!");

        List<XYDataPoint> result =
                containingNode.recursivelyFindClosestDataPoints(
                        originPoint, numPointsToFind,
                        seedPoints, farthestDistanceSoFar, null);

        DistanceToInterestingPointComparatorAsc distComp =
                new DistanceToInterestingPointComparatorAsc(originPoint);
        Collections.sort(result, distComp);



        //if too many points, sort them by distance and get rid of the excess
//System.err.println("  result size " + result.size());
        if (result.size() > numPointsToFind)
        {
            while (result.size() > numPointsToFind)
            {
//System.err.println("      removed " + SpatialQuadrantTreeNode.calculateEuclidianDistance(originPoint.getX(), originPoint.getY(), result.get(result.size()-1).getX(), result.get(result.size()-1).getY()));
                result.remove(result.size()-1);
            }
        }

if (_log.isDebugEnabled() && isInScaffRange(originPoint.getX(), originPoint.getY() ))
{
    StringBuffer sb = new StringBuffer();
    for (XYDataPoint point : result)
        sb.append(point + "DIST" +calculateEuclidianDistance(point.getX(), point.getY(),originPoint.getX(), originPoint.getY() ));
    System.err.println(sb);
}        

        return result;
    }





    /**
     * Find the closest nodes to a given starting node
     * @param originNode
     * @param numPointsToFind
     * @return a list of nodes, sorted by ascending distance
     */
    public static List<XYDataPoint> findClosestDataPointsToPointInNode(
            SpatialQuadrantTreeNode<List<XYDataPoint>> originNode,
            XYDataPoint originPoint,
            int numPointsToFind)
    {
        if (originNode.data == null)
            throw new IllegalArgumentException("findClosestDataPointsToPointInNode was called on non-leaf node");

        List<XYDataPoint> result =
                originNode.recursivelyFindClosestDataPoints(
                        originPoint, numPointsToFind,
                        new ArrayList<XYDataPoint>(),
                        0.0, null);

        DistanceToInterestingPointComparatorAsc distComp =
                new DistanceToInterestingPointComparatorAsc(originPoint);
        Collections.sort(result, distComp);

        //if too many points, sort them by distance and get rid of the excess
        if (result.size() > numPointsToFind)
        {
            while (result.size() > numPointsToFind)
                result.remove(result.size()-1);
        }

        return result;
    }

    public String toString()
    {
        return "Node, Boundaries: minX=" + minBoundaryX + ", maxX=" + maxBoundaryX +
                      ", minY=" + minBoundaryY + ", maxY=" + maxBoundaryY + ", leaf? " + isLeafNode();
    }


    public List<XYDataPoint> recursivelyFindClosestDataPoints(
            XYDataPoint originPoint,
            int numPointsToFind,
            List<XYDataPoint> identifiedPointsList,
            double farthestDistanceSoFar,
            SpatialQuadrantTreeNode<T> referringNode)
    {                                                       
        //Order the children by their distance from the origin point,
        //toss in the parent, with the closest possible sibling distance to the origin,
        //then search each one if it makes sense to
        Map<Double, SpatialQuadrantTreeNode> minDistanceNodeMap =
                new HashMap<Double, SpatialQuadrantTreeNode>();
        List<Double> distanceList = new ArrayList<Double>();
        if (lowerLeftChild != null && lowerLeftChild != referringNode)
        {
            double distance = lowerLeftChild.calculateMinPointDistanceToPoint(originPoint);
            minDistanceNodeMap.put(distance, lowerLeftChild);
            distanceList.add(distance);
        }
        if (lowerRightChild != null && lowerRightChild != referringNode)
        {
            double distance = lowerRightChild.calculateMinPointDistanceToPoint(originPoint);
            minDistanceNodeMap.put(distance, lowerRightChild);
            distanceList.add(distance);
        }
        if (upperLeftChild != null && upperLeftChild != referringNode)
        {
            double distance = upperLeftChild.calculateMinPointDistanceToPoint(originPoint);
            minDistanceNodeMap.put(distance, upperLeftChild);
            distanceList.add(distance);
        }
        if (upperRightChild != null && upperRightChild != referringNode)
        {
            double distance = upperRightChild.calculateMinPointDistanceToPoint(originPoint);
            minDistanceNodeMap.put(distance, upperRightChild);
            distanceList.add(distance);
        }
        if (parentNode != null && parentNode != referringNode)
        {
            double minDistanceToSibling =
                    calculateDistanceFromPointInBoundariesToBoundary(originPoint.getX(),
                            originPoint.getY());
//if (originPoint.getX() > 26 && originPoint.getX() < 26.1 &&
//    originPoint.getY() > -19.1 && originPoint.getY() < -19.0)
//{
//    System.err.println("Checking parent. Boundaries this box=" + getMinBoundaryX() + "," + getMaxBoundaryX() + ", " + getMinBoundaryY() + "," + getMaxBoundaryY() + ".  Distance: " + minDistanceToSibling);
//}
            minDistanceNodeMap.put(minDistanceToSibling,
                    (SpatialQuadrantTreeNode) parentNode);
            distanceList.add(minDistanceToSibling);
        }

        Collections.sort(distanceList);

        for (double distance : distanceList)
        {
            if (distance < farthestDistanceSoFar ||
               (identifiedPointsList.size() < numPointsToFind))
            {
                SpatialQuadrantTreeNode nodeToCheck = minDistanceNodeMap.get(distance);
//if (nodeToCheck == parentNode) System.err.println("Checking parent, " + nodeToCheck + ", distance " + distance);
//if (nodeToCheck == lowerLeftChild) System.err.println("Checking lowerLeft, " + nodeToCheck + ", distance " + distance);
//if (nodeToCheck == lowerRightChild) System.err.println("Checking lowerRight, " + nodeToCheck + ", distance " + distance);
//if (nodeToCheck == upperLeftChild) System.err.println("Checking upperLeft, " + nodeToCheck + ", distance " + distance);
//if (nodeToCheck == upperRightChild) System.err.println("Checking upperRight, " + nodeToCheck + ", distance " + distance);

                int sizeBefore = identifiedPointsList.size();
                nodeToCheck.recursivelyFindClosestDataPoints(originPoint, numPointsToFind,
                    identifiedPointsList, farthestDistanceSoFar, this);
                if (identifiedPointsList.size() > sizeBefore)
                {
                    XYDataPoint farthestPointSoFar =
                            identifiedPointsList.get(identifiedPointsList.size()-1);
                    farthestDistanceSoFar =
                            calculateEuclidianDistance(originPoint.getX(), originPoint.getY(),
                                                      farthestPointSoFar.getX(), farthestPointSoFar.getY());
                }


            }
        }

        //Leaf-node action.  Always guarantee that farthest point is last
        if (data != null)
        {
            for (XYDataPoint dataPoint : data)
            {
                if (identifiedPointsList.contains(dataPoint))
                    continue;

                double distanceToOrigin =
                        calculateEuclidianDistance(originPoint.getX(), originPoint.getY(),
                                dataPoint.getX(), dataPoint.getY());

                //add this node if it's closer than the farthest one we've found, or
                //if we haven't found enough nodes yet
                boolean addedPoint = false;
                if (distanceToOrigin < farthestDistanceSoFar)
                {
                    //add to the beginning of the list, so the farthest point remains last
                    addedPoint=true;
                    identifiedPointsList.add(0, dataPoint);
//                    _log.debug("Added datapoint " + data + " with distance " + distanceToOrigin);
                }
                else
                {
                    if  (identifiedPointsList.size() < numPointsToFind)
                    {
                        //add to end of list
                        addedPoint=true;
                        identifiedPointsList.add(dataPoint);
                        farthestDistanceSoFar = distanceToOrigin;
//                        _log.debug("Added datapoint " + data + " with distance " + distanceToOrigin);
                    }
                }
//if (addedPoint && originPoint.getX() > 26 && originPoint.getX() < 26.1 &&
//    originPoint.getY() > -19.1 && originPoint.getY() < -19.0)
//{
//    System.err.println("Adding point " + dataPoint + ", distance=" + distanceToOrigin);
//}
            }
        }


        
        return identifiedPointsList;
    }


    protected static class DataPointXComparatorAsc implements Comparator<XYDataPoint>
    {
        public int compare(XYDataPoint point1, XYDataPoint point2)
        {
            return Double.compare(point1.getX(), point2.getX());
        }
    }


    protected static class DataPointYComparatorAsc implements Comparator<XYDataPoint>
    {
        public int compare(XYDataPoint point1, XYDataPoint point2)
        {
            return Double.compare(point1.getY(), point2.getY());
        }
    }

    public static class DistanceToInterestingPointComparatorAsc
            implements Comparator<XYDataPoint>
    {
        XYDataPoint interestingPoint = null;

        public DistanceToInterestingPointComparatorAsc(XYDataPoint originPoint)
        {
            this.interestingPoint = originPoint;
        }

        public int compare(XYDataPoint point1, XYDataPoint point2)
        {
            return Double.compare(calculateEuclidianDistance(point1.getX(),
                                          point1.getY(),
                                          interestingPoint.getX(), interestingPoint.getY()),
                                 calculateEuclidianDistance(point2.getX(),
                                          point2.getY(),
                                          interestingPoint.getX(), interestingPoint.getY()));
        }
    }

    public double calculateMinPointDistanceToPoint(XYDataPoint dataPoint)
    {
        return calculateMinPointDistanceToPoint(dataPoint.getX(), dataPoint.getY());
    }

    /**
     * Calculate the _minimum_ distance that could possibly be between
     * any point under this node and the target point.
     *
     * The actual distance between the closest point under this node and the
     * target point may be higher, in fact.  We're just calculating the closest distance
     * between the bounding box that surrounds all points under this node
     * (rather than the boundaries of the node)
     *
     * If the point is out of the boundaries in both dimensions, this is simply
     * the distance to the closest corner.  If it's inside in one dimension, this is
     * the distance to the closest spot on the closest side.  Inside in /both/
     * dimensions works out fine in the second category (comes out to 0). 
     * @param pointX
     * @param pointY
     * @return
     */
    public double calculateMinPointDistanceToPoint(double pointX, double pointY)
    {
        double thisXValForCalculation;
        double thisYValForCalculation;

        if (pointX < minDataX)
            thisXValForCalculation = minDataX;
        else if (pointX > maxDataX)
            thisXValForCalculation = maxDataX;
        else
            thisXValForCalculation = pointX;

        if (pointY < minDataY)
            thisYValForCalculation = minDataY;
        else if (pointY > maxDataY)
            thisYValForCalculation = maxDataY;
        else
            thisYValForCalculation = pointY;

        return calculateEuclidianDistance(thisXValForCalculation,
                                          thisYValForCalculation,
                                          pointX, pointY);
    }

    public double getMiddleXByBoundaries()
    {
        return (maxBoundaryX + minBoundaryX) / 2;
    }

    public double getMiddleYByBoundaries()
    {
        return (maxBoundaryY + minBoundaryY) / 2;
    }

    /**
     * calc the maximum distance between a given point and any point inside
     * the boundaries of this node.
     *
     * This will always be calculated between the point and a corner of the boundaries
     *
     * Assumption: point is within the boundaries
     * @param pointX
     * @param pointY
     * @return
     */
    public double calculateDistanceFromPointInBoundariesToBoundary(double pointX, double pointY)
    {
        double closestXBoundary = maxBoundaryX;
        double closestYBoundary = maxBoundaryY;

        if (pointX < getMiddleXByBoundaries())
            closestXBoundary= minBoundaryX;
        if (pointY < getMiddleYByBoundaries())
            closestYBoundary= minBoundaryY;

        return Math.min(Math.abs(pointX-closestXBoundary), Math.abs(pointY-closestYBoundary));
    }


    public static double calculateEuclidianDistance(double x1, double y1, double x2, double y2)
    {
        return Math.sqrt(Math.pow(x1-x2,2) + Math.pow(y1-y2,2));
    }




    public SpatialQuadrantTreeNode<T> getUpperLeftChild()
    {
        return upperLeftChild;
    }

    public void setUpperLeftChild(SpatialQuadrantTreeNode<T> child)
    {
        if (getUpperLeftChild() != null)
            super.removeChild(this.upperLeftChild);
        super.addChild(child);
        this.upperLeftChild = child;
        child.setDepth(depth+1);
        updateDataBoundaries(child);
    }

    public SpatialQuadrantTreeNode<T> getUpperRightChild()
    {
        return upperRightChild;
    }

    public void setUpperRightChild(SpatialQuadrantTreeNode<T> child)
    {
        if (getUpperRightChild() != null)
            super.removeChild(this.upperRightChild);
        super.addChild(child);
        this.upperRightChild = child;
        child.setDepth(depth+1);
        updateDataBoundaries(child);
    }

    public SpatialQuadrantTreeNode<T> getLowerLeftChild()
    {
        return lowerLeftChild;
    }

    public void setLowerLeftChild(SpatialQuadrantTreeNode<T> child)
    {
        if (getLowerLeftChild() != null)
            super.removeChild(this.lowerLeftChild);
        super.addChild(child);
        this.lowerLeftChild = child;
        child.setDepth(depth+1);
        updateDataBoundaries(child);
    }

    public SpatialQuadrantTreeNode<T> getLowerRightChild()
    {
        return lowerRightChild;
    }

    public void setLowerRightChild(SpatialQuadrantTreeNode<T> child)
    {
        if (getLowerRightChild() != null)
            super.removeChild(this.lowerRightChild);
        super.addChild(child);
        this.lowerRightChild = child;
        child.setDepth(depth+1);
        updateDataBoundaries(child);
    }

    protected void updateDataBoundaries(SpatialQuadrantTreeNode<T> newChild)
    {

        if (newChild.getMinDataX() < minDataX)
            minDataX = newChild.getMinDataX();
        if (newChild.getMaxDataX() > maxDataX)
            maxDataX = newChild.getMaxDataX();
        if (newChild.getMinDataY() < minDataY)
            minDataY = newChild.getMinDataY();
        if (newChild.getMaxDataY() > maxDataY)
            maxDataY = newChild.getMaxDataY();
    }

    public void addChild(TreeNode<T> child)
    {
        throw new RuntimeException("This method (SpatialQuadrantTreeNode.addChild()) " +
                " should not be called directly");
    }

    public double getMinDataX()
    {
        return minDataX;
    }

    public double getMaxDataX()
    {
        return maxDataX;
    }

    public double getMinDataY()
    {
        return minDataY;
    }

    public double getMaxDataY()
    {
        return maxDataY;
    }


    public double getMinBoundaryX()
    {
        return minBoundaryX;
    }

    public void setMinBoundaryX(double minBoundaryX)
    {
        this.minBoundaryX = minBoundaryX;
    }

    public double getMaxBoundaryX()
    {
        return maxBoundaryX;
    }

    public void setMaxBoundaryX(double maxBoundaryX)
    {
        this.maxBoundaryX = maxBoundaryX;
    }

    public double getMinBoundaryY()
    {
        return minBoundaryY;
    }

    public void setMinBoundaryY(double minBoundaryY)
    {
        this.minBoundaryY = minBoundaryY;
    }

    public double getMaxBoundaryY()
    {
        return maxBoundaryY;
    }

    public void setMaxBoundaryY(double maxBoundaryY)
    {
        this.maxBoundaryY = maxBoundaryY;
    }


    public int getDepth()
    {
        return depth;
    }

    public void setDepth(int depth)
    {
        this.depth = depth;
    }
}
