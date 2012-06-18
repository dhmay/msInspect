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
package org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.AmtExtraInfoDef;

import javax.swing.*;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Constructor;

/**
 * Feature has a HashMap<String, Object> of properties.  Anything that isn't defined
 * in a field of Feature goes there.
 *
 * Also, Properties at the FeatureSet level can be managed by FeatureExtraInformationDefs.
 * For instance, only MS2 FeatureSets care about the "modifications" and "search_database_path"
 * Properties, so those properties are handled explicitly by MS2ExtraInfoDef.  That handling
 * is somewhat less rigid than the handling for Feature-level properties, in that you don't
 * need to declare the Class of the property ahead of time, just ways to convert the property
 * to and from String form.
 *
 * Groups of properties are tied together in instances of the FeatureExtraInformationDef
 * class.  Not the actual values, but the definitions of them: for each property, a
 * column name and datatype (Class).
 *
 * When FeatureSet loads a file, if it finds any
 * column from a given ExtraInformationDef, it considers itself to possess that set of
 * data -- i.e., all of those columns will be written out when the FeatureSet is
 * saved to a file.  If two column names from different ExtraInformationDefs are the same,
 * an exception will be thrown, since that would cause a collision in the hashtable.
 *
 * When reading data from and writing to files, FeatureExtraInformationDef uses a
 * single-String constructor to instantiate an object of the appropriate type and the
 * toString() method to write it out.
 *
 * The vanilla way to reference these data is with Feature.getProperty() and
 * Feature.setProperty().  But that's dangerous, because there's no compile-time validation
 * that you got the column name right, or the class of the property.  So
 * FeatureExtraInformationDef can be subclassed to provide static methods to do this
 * validation.  For instance: MS2ExtraInfoDef.getPeptideProphet(feature)
 *
 * Subclasses can also do their own thing to read in and write out data, instead of using
 * single-String constructors and toString(), in case that ever becomes necessary.
 *
 * FeatureExtraInformationDefs may also provide popup menus that appear when a feature
 * containing this extra information type is right-clicked from the main image pane
 * (see createPopupMenuItems(Feature feature))
 */
public class FeatureExtraInformationDef
{
    protected static final String MULTI_VALUE_LIST_SEPARATOR = ";";


    static Logger _log = Logger.getLogger(FeatureExtraInformationDef.class);

    //static initialization stuff

    public static final FeatureExtraInformationDef intensityWindowsInformationDef =
            new FeatureExtraInformationDef(
                    "INTENSITY_WINDOWS",
                    new String[]{"window"},
                    new Class[]{String.class}
            );

    public static final FeatureExtraInformationDef cidExtraInfoDef =
            new FeatureExtraInformationDef(
                    "CID",
                    new String[]{"cidscan"},
                    new Class[]{Integer.class}
            );

    //all known extra info types.  In theory this could contain custom info types
    //not envisioned at development time
    protected static List<FeatureExtraInformationDef> _knownExtraInfoTypes;

    //Initialize the standard info types.  In a static block for clarity
    //TODO: find a good way to initialize custom info types
    static
    {
        _standardExtraInformationTypes =
                    new FeatureExtraInformationDef[]{
                            IsotopicLabelExtraInfoDef.getSingletonInstance(),
                            MS2ExtraInfoDef.getSingletonInstance(),
                            TimeExtraInfoDef.getSingletonInstance(),
                            AmtExtraInfoDef.getSingletonInstance(),
                            intensityWindowsInformationDef,
                            cidExtraInfoDef
                    };
        for (FeatureExtraInformationDef def : getStandardExtraInformationTypes())
        {
            addKnownExtraInfoType(def);
        }
    }



    //instance variables

    
    //each one of these needs a name, primarily so that it can appear in an error
    //string if anything goes wrong
    protected String textCode;
    protected Map<String, Class<?>> columnNameDatatypeMap;
    //unfortunately we have to store the column names in an array, too, to preserve ordering
    protected String[] columnNames;
    //feature set properties
    protected String[] featureSetPropertyNames = null;

    public FeatureExtraInformationDef()
    {
        init();
    }

    public FeatureExtraInformationDef(String name, String[] columnNames, Class<?>[] dataTypes)
    {
        this();
        setTextCode(name);
        if (columnNames.length != dataTypes.length)
            throw new RuntimeException("FeatureSetExtraInformation: not all column datatypes defined");

        for (int i = 0; i < columnNames.length; i++)
            columnNameDatatypeMap.put(columnNames[i], dataTypes[i]);
        this.columnNames = columnNames;
    }

    public FeatureExtraInformationDef(String name, String[] columnNames, Class[] dataTypes,
                                      String[] featureSetPropertyNames)
    {
        this(name, columnNames, dataTypes);
        this.featureSetPropertyNames = featureSetPropertyNames;
    }

    public String getTextCode()
    {
        return textCode;
    }

    public void setTextCode(String textCode)
    {
        this.textCode = textCode;
    }

    protected Object getFeatureSetProperty(FeatureSet featureSet, String propertyName)
    {
        return featureSet.getProperty(createFeatureSetPropertyName(propertyName));
    }

    protected void setFeatureSetProperty(FeatureSet featureSet, String propertyName,
                                         Object propertyValue)
    {
        _log.debug("Setting featureset property " + propertyName + " to " + propertyValue);
        featureSet.setProperty(createFeatureSetPropertyName(propertyName), propertyValue);
    }

    public String createFeatureSetPropertyName(String propertyName)
    {
        return textCode + ":" + propertyName;
    }

    public String stripPrefixFromFeatureSetPropertyName(String propertyName)
    {
        return propertyName.substring(propertyName.indexOf(":") + 1);
    }

    public boolean isThisTypeOfFeatureSetProperty(String propertyName)
    {
        return propertyName.startsWith(textCode + ":");
    }

    /**
     * Doing the call to TextProvider here so that custom subclasses can
     * do whatever they want instead
     * @return
     */
    public String getTranslatedText()
    {
        return TextProvider.getText(getTextCode());
    }

    protected void init()
    {
        columnNameDatatypeMap = new HashMap<String, Class<?>>();
    }

    /**
     * Return the array that preserves original column order
     * @return
     */
    public String[] getColumnNames()
    {
        return columnNames;
    }

    public Class<?> getDatatypeForColumnName(String columnName)
    {
        return columnNameDatatypeMap.get(columnName);
    }

    public String toString()
    {
        StringBuffer resultBuf = new StringBuffer("FeatureExtraInformation, columns: ");

        for (String columnName : getColumnNames())
            resultBuf.append(columnName + " (" +
                    getDatatypeForColumnName(columnName).getName() + "), ");

        return resultBuf.toString();
    }

    /**
     * Convert a String value to the appropriate datatype for the column.
     * In this implementation, simply calls a one-arg constructor.
     * This can be overridden in a subclass to do something more interesting,
     * if necessary.
     *
     * Many things can conceivably go wrong here.  The column might not be
     * defined in the hashmap (NullPointerException).  The class might not
     * have a one-String constructor.  That constructor might not be public.
     * The constructor might fail on the provided String.  Any of these
     * things will throw an exception
     *
     * @param columnName
     * @param value
     * @return
     * @throws InstantiationException
     * @throws IllegalAccessException
     * @throws InvocationTargetException
     * @throws NoSuchMethodException
     */
    public Object convertStringValue(String columnName, String value)
            throws InstantiationException, IllegalAccessException,
            InvocationTargetException, NoSuchMethodException
    {
        Class dataTypeClass = getDatatypeForColumnName(columnName);
        //if we don't need to convert, don't convert
        if (String.class.isAssignableFrom(dataTypeClass))
            return value;

        Constructor<String> constructor =
                dataTypeClass.getConstructor(String.class);
        return constructor.newInstance(value);
    }

    /**
     * For output to feature file.  This is used so that the Class containing
     * the value doesn't necessarily have to implement a toString() method that does
     * what we want it to.
     *
     * In the default implementation, however, that's what happens.
     * @param columnName
     * @param value
     * @return
     */
    public String convertToString(String columnName, Object value)
    {
        if (value == null)
            return null;
        return value.toString();
    }

    /**
     * Same as convertToString, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public String convertFeatureSetPropertyToString(String propertyName, Object value)
    {
        if (value == null)
            return null;
        return value.toString();
    }

    /**
     * Save as convertStringValue, but for feature set properties
     * @param propertyName
     * @param value
     * @return
     */
    public Object convertFeatureSetPropertyStringValue(String propertyName, String value)
    {
        return value;
    }


    protected static FeatureExtraInformationDef[] _standardExtraInformationTypes;
    protected static Map<String, Class> _knownColumnDataClassMap;

    public static Map<String, Class> getKnownColumnDataClassMap()
    {
        if (_knownColumnDataClassMap == null)
            _knownColumnDataClassMap = new HashMap<String, Class>();
        return _knownColumnDataClassMap;
    }

    public static Class getKnownDataClassForColumn(String columnName)
    {
        return getKnownColumnDataClassMap().get(columnName);
    }

    public static List<FeatureExtraInformationDef> getKnownExtraInfoTypes()
    {
        if (_knownExtraInfoTypes == null)
        {
            _knownExtraInfoTypes =
                    new ArrayList<FeatureExtraInformationDef>(_standardExtraInformationTypes.length);

        }
        return _knownExtraInfoTypes;
    }


    protected static Map<String, FeatureExtraInformationDef> _columnInfoTypeMap;

    public static Map<String, FeatureExtraInformationDef> getColumnInfoTypeMap()
    {
        if (_columnInfoTypeMap == null)
        {
            _columnInfoTypeMap =
                    new HashMap<String, FeatureExtraInformationDef>();
        }
        return _columnInfoTypeMap;
    }

    /**
     * This method is how you make the system aware of custom information types
     * @param infoType
     */
    public static void addKnownExtraInfoType(FeatureExtraInformationDef infoType)
    {
        getKnownExtraInfoTypes().add(infoType);
        for (String columnName : infoType.getColumnNames())
        {
            if (getColumnInfoTypeMap().containsKey(columnName))
                _log.info("WARNING!! double-defining extra info type with column " +
                          columnName + "; new info type is " + infoType.getTextCode() +
                          ", will override older one");
            getColumnInfoTypeMap().put(columnName, infoType);
            getKnownColumnDataClassMap().put(columnName, infoType.getDatatypeForColumnName(columnName));
        }
    }

    public static FeatureExtraInformationDef getInfoTypeForColumn(String columnName)
    {
        return getColumnInfoTypeMap().get(columnName);
    }

    public static FeatureExtraInformationDef[] getStandardExtraInformationTypes()
    {
        return _standardExtraInformationTypes;
    }

    /**
     * This method should, in the subclasses, add popup menu items that should be
     * shown when the feature is right-clicked from the main image pane
     * @param feature
     * @return
     */
    public List<JMenuItem> createPopupMenuItems(Feature feature)
    {
        return new ArrayList<JMenuItem>();
    }

    public static String convertStringArrayToString(String[] stringArray)
    {
        return convertStringArrayToString(stringArray, MULTI_VALUE_LIST_SEPARATOR);
    }

    public static String convertStringArrayToString(String[] stringArray, String separatorString)
    {
        if (stringArray == null || stringArray.length == 0)
            return "";

        //check to see if every value is null.  If so, return empty string.
        //This is for situations where we've added a bunch of peptides
        //without adding proteins.  We keep adding items to the list, to preserve
        //association between proteins and peptides, but there's no need to write
        //out the separators
        boolean hasNonNullValue = false;
        for (String value : stringArray)
        {
            if (value != null) hasNonNullValue = true;
        }
        if (!hasNonNullValue)
            return "";        
        StringBuffer result = new StringBuffer(stringArray[0]);
        for (int i = 1; i < stringArray.length; i++)
            result.append(separatorString + stringArray[i]);
        return result.toString();
    }

    public static String convertIntListToString(List<Integer> intList)
    {
        if (intList == null)
            return "";
        List<String> intListAsString = new ArrayList<String>(intList.size());
        for (Integer val : intList)
            intListAsString.add("" + val);
        return convertStringListToString(intListAsString);
    }

    public static String convertStringListToString(List<String> stringList)
    {
        return convertStringListToString(stringList, MULTI_VALUE_LIST_SEPARATOR);
    }

    /**
     * Convert a list of peptides or proteins into a String for storage in a feature file
     *
     * @param stringList
     * @return
     */
    public static String convertStringListToString(List<String> stringList, String listSeparator)
    {
        if (stringList == null)
            return "";
        return convertStringArrayToString(stringList.toArray(new String[stringList.size()]),
                                          listSeparator);
    }

    public static List<String>
            parseFormulaListString(String stringListString)
    {
        return parseStringListString(stringListString, MULTI_VALUE_LIST_SEPARATOR);
    }

    public static List<String>
            parseStringListString(String stringListString)
    {
        return parseStringListString(stringListString, MULTI_VALUE_LIST_SEPARATOR);
    }

    public static List<Integer> parseIntListString(String intListString)
    {
        if (intListString == null)
            return new ArrayList<Integer>();
        List<String> stringList = parseStringListString(intListString);
        List<Integer> intList = new ArrayList<Integer>(stringList.size());
        for (String string : stringList)
            intList.add(Integer.parseInt(string));
        return intList;
    }

    /**
     * Parse the value of a peptide or protein list column in a feature file
     * @param stringListString
     * @return
     */
    public static List<String>
            parseStringListString(String stringListString, String separatorString)
    {
        List<String> result = new ArrayList<String>();
        String[] stringArray =
                stringListString.split(separatorString);
        for (String string : stringArray)
        {
            result.add(string);
        }
        return result;
    }
}
