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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.viewer.Application;

import javax.swing.table.AbstractTableModel;
import java.util.*;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: Feb 1, 2005
 * Time: 2:29:14 PM
 */
public class FeatureTableModel extends AbstractTableModel
    {
    Feature[] _features;


    Pair[] _properties =
        {
        new Pair("Scan",Integer.class),
        new Pair("MZ", Float.class),
        new Pair("Intensity", Float.class),
        new Pair("Charge",Integer.class),
        new Pair("KL", Float.class),
        new Pair("Peaks",Integer.class),
//		new Pair("Background", Float.class),
//		new Pair("Median", Float.class),
        new Pair("Mass", Float.class),
        new Pair("Description", String.class)
        };


    public FeatureTableModel(Feature[] features)
        {
        _features = features;
        }


    public int getColumnCount()
        {
        return _properties.length;
        }


    public String getColumnName(int column)
        {
        return (String)_properties[column].first;
        }

    //TODO: Put editability in the _properties array
    private static Set editableCols = new HashSet();
    static
        {
        editableCols.add("MZ");
        editableCols.add("Scan");
        editableCols.add("Charge");
        editableCols.add("Peaks");
        editableCols.add("Description");
        editableCols.add("Intensity");
        }
    public boolean isCellEditable(int rowIndex, int columnIndex)
        {
        Pair prop = _properties[columnIndex];
        String val = (String) prop.first;
        return editableCols.contains(val);
        }

    public Class getColumnClass(int column)
        {
        return (Class)_properties[column].second;
        }

    public int getRowCount()
        {
        return null == _features ? 0 : _features.length;
        }

    public Object getValueAt(int rowIndex, int columnIndex)
        {
        Feature f = _features[rowIndex];
        String prop = (String)_properties[columnIndex].first;
        switch (prop.charAt(0))
            {
        case 'S':
            return new Integer(f.getScan());
        case 'M':
            if ("MZ".equals(prop))
                return new Float(f.getMz());
            if ("Median".equals(prop))
                return new Float(f.getMedian());
            if ("Mass".equals(prop))
                return new Float(f.getMass());
            break;
        case 'I':
            return new Float(f.getIntensity());
        case 'B':
            return new Float(f.getBackground());
        case 'C':
            return new Integer(f.getCharge());
        case 'K':
            return new Float(f.getKl());
        case 'P':
            return new Integer(f.getPeaks());
        case 'D': //Description
            return f.getDescription();
            }
        throw new IllegalArgumentException();
        }

    //CONSIDER: Use bean utils
    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
        {
        Feature f = _features[rowIndex];
        String prop = (String)_properties[columnIndex].first;

        saveOrigValue(prop, f, rowIndex, columnIndex);
        switch (prop.charAt(0))
            {
        case 'S':
            f.setScan(((Integer) aValue).intValue());
            break;
        case 'M':
            if ("MZ".equals(prop))
                {
                f.setMz(((Float) aValue).floatValue());
                f.updateMass();
                }
            break;
        case 'C':
            f.setCharge(((Integer) aValue).intValue());
            f.updateMass();
            break;
        case 'P':
            f.setPeaks(((Integer) aValue).intValue());
            break;
        case 'I':
            f.setIntensity(((Float) aValue).floatValue());
            break;
        case 'D':
            f.setDescription((String) aValue);
            break;
            }

        List featureSets = (List) Application.getInstance().getProperty(SharedProperties.FEATURE_SETS);
        FeatureSelectionFrame.FeatureSelectionDialog.getInstance().updateFeatureSets(featureSets, false);
        return;
        }

    private void saveOrigValue(String propName, Feature f, int rowIndex, int colIndex)
        {
        if (propName == "Description")
            return;

        Object oldVal = getValueAt(rowIndex, colIndex);
        String oldDesc = f.getDescription();
        if (null == oldDesc)
            oldDesc = "";
        String match = "orig" + propName + "=";
        if (oldDesc.indexOf(match) < 0)
            {
            if (oldDesc.length() > 0)
                oldDesc += "; ";

            f.setDescription(oldDesc + "orig" + propName + "=" + oldVal.toString());
            }
        }

    }
