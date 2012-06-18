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
 * Tree node that's aware of its parent
 */
public class ParentAwareTreeNode<T> extends TreeNode<T>
{
    protected ParentAwareTreeNode parentNode = null;

    public void addChild(ParentAwareTreeNode<T> child)
    {
        super.addChild(child);
        child.setParentNode(this);
    }

    public ParentAwareTreeNode getParentNode()
    {
        return parentNode;
    }

    public void setParentNode(ParentAwareTreeNode parentNode)
    {
        this.parentNode = parentNode;
    }

    public boolean isRoot()
    {
        return parentNode == null;
    }
}
