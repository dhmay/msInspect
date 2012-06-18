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

import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;


/**
 * An in-order (left child < parent < right child) binary tree implementation
 *
 * This is a Tree that only allows two children, uses a Comparator to ensure that the
 * left child is somehow "less than or equal to" the parent and the right child is "greater".
 * Keeps track of which child is left and which is right.
 *
 * Be careful with this one... I haven't overridden methods in the superclass that might
 * cause trouble.  Also, there should be another layer of indirection here -- a BinaryTree class
 * that this descends from.  Maybe someday, if we need it.
 *
 */
public class InOrderBinaryTree<T> extends Tree<T>
{
    protected Comparator<T> comparator = null;

    public InOrderBinaryTree(Comparator<T> comparator)
    {
        this.comparator = comparator;
    }

    public void addNode(BinaryTreeNode<T> newNode)
    {
        if (getRootElement() == null)
            setRootElement(newNode);
        else
        {
            addNode((BinaryTreeNode) getRootElement(), newNode);
        }
    }

    /**
     * by convention, if comparator returns equality, stick the new node to the left
     * @param parentNode
     * @param newNode
     */
    public void addNode(BinaryTreeNode<T> parentNode, BinaryTreeNode<T> newNode)
    {
        int result = comparator.compare(parentNode.getData(), newNode.getData());
        if (result >= 0)
        {
            if (parentNode.hasLeftChild())
                addNode(parentNode.getLeftChild(), newNode);
            else
                parentNode.setLeftChild(newNode);
        }
        else
        {
            if (parentNode.hasRightChild())
                addNode(parentNode.getRightChild(), newNode);
            else
                parentNode.setRightChild(newNode);
        }
    }

    /**
     * Add a bunch of nodes at once.  This should really be in Tree, but I didn't want to
     * mess with Tree
     * @param newNodes
     */
    public void populateTree(List<BinaryTreeNode<T>> newNodes)
    {
        for (BinaryTreeNode<T> newNode : newNodes)
            addNode(newNode);
    }

    /**
     * Cover method for the root node
     * @param minValue
     * @param maxValue
     * @return
     */
    public List<BinaryTreeNode<T>> getNodesInRange(T minValue, T maxValue)
    {
        return getNodesInRange((BinaryTreeNode<T>) getRootElement(), minValue, maxValue);
    }

    /**
     * Recursive.  Given a parent node, return all the nodes of the subtree rooted at that
     * parent that fall within the range provided.
     * @param parent
     * @param minValue
     * @param maxValue
     * @return
     */
    public List<BinaryTreeNode<T>> getNodesInRange(BinaryTreeNode<T> parent, T minValue, T maxValue)
    {
        boolean greaterThanMin = false;
        boolean lessThanMax = false;

        List<BinaryTreeNode<T>> result = new ArrayList<BinaryTreeNode<T>>();

        if (comparator.compare(parent.getData(), minValue) >= 0)
            greaterThanMin = true;
        if (comparator.compare(parent.getData(), maxValue) <= 0)
            lessThanMax = true;

        if (lessThanMax && parent.hasRightChild())
                result.addAll(getNodesInRange(parent.getRightChild(), minValue, maxValue));
        if (greaterThanMin && parent.hasLeftChild())
                result.addAll(getNodesInRange(parent.getLeftChild(), minValue, maxValue));          

        if (greaterThanMin && lessThanMax)
            result.add(parent);
        return result;
    }

    /**
     * Convenience method to look at the nodes from getNodesInRange() and spit out the
     * data values
     * @param minValue
     * @param maxValue
     * @return
     */
    public List<T> getDataValuesInRange(T minValue, T maxValue)
    {
        List<BinaryTreeNode<T>> nodesInRange = getNodesInRange(minValue, maxValue);

        List<T> result = new ArrayList<T>(nodesInRange.size());
        for (BinaryTreeNode<T> nodeInRange : nodesInRange)
            result.add(nodeInRange.getData());
        return result;
    }
}

