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

package org.fhcrc.cpl.viewer.amt;

import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;

/**
 * Creates instances of particular FeatureSetMatchers 
 */
public class AmtFeatureSetMatcherFactory
{
    public static final String AMT_PACKAGE = "org.fhcrc.cpl.viewer.amt";

    /**
     * Given a partial class name, tacks on "FeatureSetMatcher", puts it in the right package, and
     * looks for it.
     * This method is not forgiving.  Capitalization matters, etc.
     * @param partialClassName
     * @return
     * @throws InstantiationException
     */
    public static FeatureSetMatcher createInstanceFromPartialClassName(String partialClassName)
            throws InstantiationException
    {
        String className =  AMT_PACKAGE + "." + partialClassName + "FeatureSetMatcher";
        try
        {
            return createInstance(Class.forName(className));

        }
        catch (ClassNotFoundException e)
        {
            throw new InstantiationException("Class " + className + " not found");
        }
    }

    /**
     * Catch different kinds of exceptions and throw InstantiationExceptions, because I'm lazy
     * farther up in the code
     * @param classToInstantiate
     * @return
     * @throws InstantiationException
     */
    public static FeatureSetMatcher createInstance(Class classToInstantiate)
            throws InstantiationException
    {
        try
        {
            return (FeatureSetMatcher) classToInstantiate.newInstance();
        }
        catch (ClassCastException e)
        {
            throw new InstantiationException("Class " + classToInstantiate.getName() +
                                             " does not implement interface FeatureSetMatcher");
        }
        catch (IllegalAccessException e)
        {
            throw new InstantiationException("IllegalAccessException while trying to instantiate class " +
                                             classToInstantiate);
        }
    }

    /**
     * Instantiate the class and then call the init() method
     * @param classToInstantiate
     * @param deltaMass
     * @param deltaMassType
     * @param deltaHydrophobicity
     * @return
     * @throws InstantiationException
     */
    public static FeatureSetMatcher createInstance(Class classToInstantiate,
                                               float deltaMass, int deltaMassType,
                                               float deltaHydrophobicity)
            throws InstantiationException
    {
        FeatureSetMatcher result = createInstance(classToInstantiate);
        result.init(deltaMass, deltaMassType, deltaHydrophobicity);

        return result;
    }
}
