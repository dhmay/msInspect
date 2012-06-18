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


/**
 * Binary tree node
 */
public class BinaryTreeNode<T> extends TreeNode<T>
{
    protected BinaryTreeNode leftChild = null;
    protected BinaryTreeNode rightChild = null;


    public void addChild(TreeNode<T> child)
    {
        if (this.getNumberOfChildren() >= 2)
        {
            throw new IllegalArgumentException("Attempted to add a 3rd child to a BinaryTreeNode");
        }

        super.addChild(child);
    }

    public void insertChildAt(int index, TreeNode<T> child) throws IndexOutOfBoundsException
    {
        if (this.getNumberOfChildren() >= 2)
        {
            throw new IllegalArgumentException("Attempted to add a 3rd child to a BinaryTreeNode.");
        }

        super.insertChildAt(index, child);
    }

    public void setLeftChild(TreeNode<T> newNode)
    {
        this.addChild(newNode);
        leftChild = (BinaryTreeNode<T>) newNode;
    }

    public void setRightChild(TreeNode<T> newNode)
    {
        this.addChild(newNode);
        rightChild = (BinaryTreeNode<T>) newNode;
    }

    public boolean hasLeftChild()
    {
        return leftChild != null;
    }

    public boolean hasRightChild()
    {
        return rightChild != null;
    }


    public BinaryTreeNode<T> getLeftChild()
    {
        return leftChild;
    }

    public BinaryTreeNode<T> getRightChild()
    {
        return rightChild;
    }


}
