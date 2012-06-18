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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;

import java.util.Map;
import java.util.Set;
import java.util.Collection;


/**
 * This is a subclass of Feature that implements Map.  All the Map methods
 * put values in and take values out of the properties Map.
 *
 * This class exists to allow TabLoader to populate arbitrary properties of
 * a Feature.  Those properties may be known at development time, and they
 * may not.  This is a separate class so that Feature isn't cluttered with
 * all these Map methods, which really don't make a lot of sense for Feature.
 */
public class FeatureAsMap extends Feature implements Map
{
    protected static Logger _log = Logger.getLogger(FeatureAsMap.class);

    /**
     * A cover on setProperty().  But before calling setProperty(), it tries to determine
     * if this property name is known and if the value is a String.  If both of those are
     * true, then, before calling setProperty, convert the property to the appropriate
     * class.
     * @param propertyName
     * @param propertyValue
     * @return
     */
    public Object put(Object propertyName, Object propertyValue)
    {
        if (!(propertyName instanceof String))
            throw new RuntimeException("FeatureAsMap: tried to set a non-String property");
        Object result = get(propertyName);
//System.err.println("*****put, propname=" + propertyName + ", value class: " + propertyValue.getClass().getTextCode());

        FeatureExtraInformationDef infoDef =
                FeatureExtraInformationDef.getInfoTypeForColumn((String) propertyName);
        try
        {
            if (infoDef != null)
            {
                propertyValue =
                        infoDef.convertStringValue((String) propertyName,
                                propertyValue.toString());
            }
        }
        catch (Exception e)
        {
            _log.error("Error setting property " + propertyName, e);
        }

        setProperty((String) propertyName, propertyValue);

        return result;
    }

    public Object get(Object propertyName)
    {
        if (!(propertyName instanceof String))
            throw new RuntimeException("FeatureAsMap: tried to set a non-String property");
        return getProperty((String) propertyName);
    }





    //no-brainer methods to implement Map

    public void putAll(Map namesAndValues)
    {
        for (Object propertyName : namesAndValues.keySet())
            put(propertyName, namesAndValues.get(propertyName));
    }

    public Object remove(Object propertyName)
    {
        return getPropertyMap().remove(propertyName);
    }

    public void clear()
    {
        getPropertyMap().clear();
    }

    public Set keySet()
    {
        return getPropertyMap().keySet();
    }

    public Collection values()
    {
        return getPropertyMap().values();
    }

    public Set entrySet()
    {
        return getPropertyMap().entrySet();
    }

    public int size()
    {
        return getPropertyMap().size();
    }

    public boolean isEmpty()
    {
        return getPropertyMap().isEmpty();
    }

    public boolean containsKey(Object o)
    {
        return getPropertyMap().containsKey(o);
    }

    public boolean containsValue(Object o)
    {
        return getPropertyMap().containsValue(o);
    }

}
