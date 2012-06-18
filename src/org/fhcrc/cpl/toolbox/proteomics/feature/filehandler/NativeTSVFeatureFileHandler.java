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
package org.fhcrc.cpl.toolbox.proteomics.feature.filehandler;

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureAsMap;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Arrays;

/**
 * File handler for native msInspect feature files
 */
public class NativeTSVFeatureFileHandler extends BaseFeatureSetFileHandler
    implements FeatureSetFileHandler
{
    static Logger _log = Logger.getLogger(NativeTSVFeatureFileHandler.class);


    public static final String FILE_TYPE_NAME = "NATIVE_TSV";

    protected static NativeTSVFeatureFileHandler singletonInstance = null;

    public static NativeTSVFeatureFileHandler getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new NativeTSVFeatureFileHandler();
        return singletonInstance;
    }

    /**
     * Load a FeatureSet
     * @param file
     * @return
     * @throws IOException
     */
    public FeatureSet loadFeatureSet(File file)
            throws IOException
    {
        FeatureSet result = new FeatureSet();
        TabLoader loader = new TabLoader(file, FeatureAsMap.class);
        TabLoader.ColumnDescriptor[] cols = loader.getColumns();

        Feature[] features = null;

        //change the columns named "start" and "end" to scanFirst and scanLast,
        //and make sure the setters get called.
        //Also, register any known extra information types that we find for
        //these columns
        for (TabLoader.ColumnDescriptor col : cols)
        {
            if (col.name.equals("start"))
            {
                col.name = "scanFirst";
                col.load = true;
            }
            else if (col.name.equals("end"))
            {
                col.name = "scanLast";
                col.load = true;
            }
            //if this column indicates a known extra info type, register it
            FeatureExtraInformationDef infoDef =
                    FeatureExtraInformationDef.getInfoTypeForColumn(col.name);
            if (infoDef != null)
            {
                if (!result.hasExtraInformationType(infoDef))
                {
                    result.addExtraInformationType(infoDef);
                }
            }
        }
        if (cols.length > 0)
        {
            features = (Feature[]) loader.load();
            for (Feature feature : features)
                feature.afterPopulate();
            result.setFeatures(features);
        }

        for (FeatureExtraInformationDef infoType : result.getExtraInformationTypes())
            _log.debug("Discovered extra information type " + infoType.getTextCode());

        Map<String,String> commentsMap = (Map<String,String>) loader.getComments();
        for (String commentKey : commentsMap.keySet())
        {
            String propertyValueString = commentsMap.get(commentKey);
            _log.debug("Loading property " + commentKey + " = " + propertyValueString);
            Object propertyValue = propertyValueString;
            //if this property is known about by an extra info def, let the def handle it
            Set<FeatureExtraInformationDef> infoDefsToTry = new HashSet<FeatureExtraInformationDef>();
            infoDefsToTry.addAll(result.getExtraInformationTypes());
            for (FeatureExtraInformationDef standardInfoDef : FeatureExtraInformationDef.getStandardExtraInformationTypes())
            {
                infoDefsToTry.add(standardInfoDef);
            }
            for (FeatureExtraInformationDef extraInfoDef : infoDefsToTry)
            {
                if (extraInfoDef.isThisTypeOfFeatureSetProperty(commentKey))
                {
                    propertyValue =
                            extraInfoDef.convertFeatureSetPropertyStringValue(
                                    extraInfoDef.stripPrefixFromFeatureSetPropertyName(commentKey),
                                    propertyValueString);
                    _log.debug("\tproperty claimed by extrainfodef " + extraInfoDef.getTextCode() +
                            ", class " + propertyValue.getClass().getName());
                    break;
                }
            }
            result.setProperty(commentKey, propertyValue);
        }

        return result;
    }

    public void saveFeatureSet(FeatureSet featureSet, File outFile)
            throws IOException
    {
       PrintWriter out = null;
        assert null != featureSet.getFeatures();
        try
        {
            out = new PrintWriter(new FileOutputStream(outFile));
            saveFeatureSet(featureSet, out);
        }
        catch(IOException e)
        {
            throw e;
        }
        finally
        {
            if (out != null)
                out.close();
        }
    }

    /**
     * Save a FeatureSet
     * @param featureSet
     * @param out
     */
    public void saveFeatureSet(FeatureSet featureSet, PrintWriter out)
    {
        assert null != featureSet.getFeatures();
        Set<?> keySet = featureSet.getProperties().keySet();
        String[] keys = new String[keySet.size()];
        int i = 0;
        for (Object key : keySet)
            keys[i++] = key.toString();

        Arrays.sort(keys);
        for (i = 0; i < keys.length; i++)
        {
            Object propertyValue = featureSet.getProperties().get(keys[i]);

            //if this property is known about by an extra info def, let the def handle it
            Set<FeatureExtraInformationDef> infoDefsToTry = new HashSet<FeatureExtraInformationDef>();
            infoDefsToTry.addAll(featureSet.getExtraInformationTypes());
            for (FeatureExtraInformationDef standardInfoDef :
                    FeatureExtraInformationDef.getStandardExtraInformationTypes())
            {
                infoDefsToTry.add(standardInfoDef);
            }
            for (FeatureExtraInformationDef extraInfoDef : infoDefsToTry)
            {
                if (extraInfoDef.isThisTypeOfFeatureSetProperty(keys[i]))
                {
                    propertyValue =
                            extraInfoDef.convertFeatureSetPropertyToString(
                                    extraInfoDef.stripPrefixFromFeatureSetPropertyName(keys[i]),
                                    propertyValue);
                    break;
                }
            }
            out.println("# " + keys[i] + "=" + propertyValue);
        }

        String header = Feature.getFeatureHeader(featureSet.getExtraInformationTypesArray());
        if (dumpWindow)
            header += "\twindow";
        out.println(header);
        out.flush();

        Feature[] features = featureSet.getFeatures();
        for (i = 0; i < features.length; i++)
        {
            out.print(features[i].toString(featureSet.getExtraInformationTypesArray()));
            if (dumpWindow && null != features[i].intensityWindow)
                out.print("\t" + join(features[i].intensityWindow, ","));
            out.println();
        }
        out.flush();
    }

    // TODO: use StringUtils instead
    private static String join(float[] f, String delim)
    {
        if ( null == f || 0 >= f.length )
            return null;
        StringBuffer sb = new StringBuffer();
        sb.append("" + f[0]);
        for(int i = 1; i < f.length; i++)
        {
            if ( null != delim)
                sb.append(delim);
            sb.append("" + f[i]);
        }
        return sb.toString();
    }

    /**
     * Can this type of file handler handle this specific file?
     * Implementation is up to the handler, but this should be as low-cost as possible
     * @param file
     * @return
     * @throws IOException
     */
    public boolean canHandleFile(File file)
        throws IOException
    {
        if (isXMLFile(file))
            return false;

        //if it's not an XML file, give it a whirl, since this is the default handler
        return true;
    }
}
