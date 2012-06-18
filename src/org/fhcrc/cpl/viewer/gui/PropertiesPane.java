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

import org.fhcrc.cpl.toolbox.datastructure.BoundMap;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithRPerspectivePlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBlindImageChart;
import org.fhcrc.cpl.toolbox.gui.chart.MultiChartDisplayPanel;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.quant.gui.PanelWithSpectrumChart;
import org.fhcrc.cpl.viewer.util.SharedProperties;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.util.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.io.File;


/**
 * User: mbellew
 * Date: Jun 9, 2005
 * Time: 8:00:44 PM
 *
 * dhmay adding functionality to display charts for the active feature
 */
public class PropertiesPane
{
    static JScrollPane _component;
    DefaultTableModel _model;
    JTable _table;

    static JDialog _popup;
    static PropertiesPane _popupPane;

    protected RightMousePopupMenu _rightMousePopupMenu = new RightMousePopupMenu();


    protected Feature selectedFeature = null;


    //parameters for displaying area around feature
    public static final int FEATURECHART_SCAN_PADDING = 5;
    public static final int FEATURECHART_MIN_PEAKS = 4;
    public static final float FEATURECHART_MZ_PADDING_UP = 0.5f;
    public static final float FEATURECHART_MZ_PADDING_DOWN = 1.5f;



    public static synchronized void popup(Object objectWithProperties)
    {
        if (null == _popup)
        {
            PropertiesPane pp = new PropertiesPane();
            JDialog d = new JDialog(ApplicationContext.getFrame(),
                    "Properties");
            d.setContentPane(pp.getComponent());
            d.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            d.setSize(400, 400);
            _popup = d;
            _popupPane = pp;
        }
        _popupPane.setProperties(objectWithProperties);
        _popup.setVisible(true);
    }


    public Container getComponent()
    {
        return _component;
    }


    public PropertiesPane()
    {
        _model = new DefaultTableModel(0, 2);
        _table = new JTable(_model);
        _table.getColumnModel().getColumn(0).setHeaderValue(TextProvider.getText("PROPERTY_LOWERCASE"));
        _table.getColumnModel().getColumn(1).setHeaderValue(TextProvider.getText("VALUE_LOWERCASE"));
        _component = new JScrollPane(_table);
    }


    public PropertiesPane(Object objectWithProperties)
    {
        this();
        setProperties(objectWithProperties);
    }


    public void setProperties(Object objectWithProperties)
    {
        //dhmay adding special handling for features
        if (Feature.class.isAssignableFrom(objectWithProperties.getClass()))
        {
            Feature feature = (Feature) objectWithProperties;

            //blank out existing properties
            clearProperties();

            //write out standard properties
            addPropertyToModel(TextProvider.getText("SCAN"),feature.getScan());
            addPropertyToModel(TextProvider.getText("TIME"),feature.getTime());
            addPropertyToModel(TextProvider.getText("MZ"),feature.getMz());
            addPropertyToModel(TextProvider.getText("MASS"),feature.getMass());
            addPropertyToModel(TextProvider.getText("CHARGE"),feature.getCharge());
            addPropertyToModel(TextProvider.getText("CHARGE_STATES"),feature.getChargeStates());
            addPropertyToModel(TextProvider.getText("PEAKS"),feature.getPeaks());
            addPropertyToModel(TextProvider.getText("INTENSITY"),feature.getIntensity());
            addPropertyToModel(TextProvider.getText("TOTAL_INTENSITY"),feature.getTotalIntensity());
            addPropertyToModel(TextProvider.getText("SCAN_COUNT"),feature.getScanCount());
            addPropertyToModel(TextProvider.getText("START_SCAN"),feature.getScanFirst());
            addPropertyToModel(TextProvider.getText("END_SCAN"),feature.getScanLast());
            addPropertyToModel(TextProvider.getText("KL"),feature.getKl());
            addPropertyToModel(TextProvider.getText("DESCRIPTION"),feature.getDescription());

            //write out rows for extra information types
            for (FeatureExtraInformationDef extraInformationType :
                    feature.determineExtraInformationTypes())
            {
                addPropertyToModel("<html><b>" + extraInformationType.getTextCode() + "</b></html>",
                                   null);
                for (String columnName : extraInformationType.getColumnNames())
                    addPropertyToModel(columnName, feature.getProperty(columnName));
            }
            selectedFeature = feature;
            _table.setComponentPopupMenu(_rightMousePopupMenu);
            _component.setComponentPopupMenu(_rightMousePopupMenu);
        }
        else if (objectWithProperties instanceof Map)
        {
            setProperties((Map<String, ?>) objectWithProperties, null);
            _table.setComponentPopupMenu(null);
            _component.setComponentPopupMenu(null);

        }
        else
        {
            setProperties(new BoundMap(objectWithProperties), null);
            _table.setComponentPopupMenu(null);
            _component.setComponentPopupMenu(null);
            
        }
    }

    protected class RightMousePopupMenu extends JPopupMenu
    {

        public RightMousePopupMenu()
        {
            JMenuItem showZoomMenuItem = new JMenuItem("Show Region");
            showZoomMenuItem.addActionListener(new ActionListener()
            {
                /**
                 * Show some charts describing the region around the selected feature
                 * @param event
                 */
                public void actionPerformed(ActionEvent event)
                {
                    if (selectedFeature == null)
                        return;
                    ApplicationContext.setMessage("Building chart...");
                    PanelWithSpectrumChart pwsc =
                            new PanelWithSpectrumChart((MSRun)ApplicationContext.getProperty(SharedProperties.MS_RUN),
                                    selectedFeature.getScanFirst()-FEATURECHART_SCAN_PADDING,
                                    selectedFeature.getScanLast()+FEATURECHART_SCAN_PADDING,
                                    selectedFeature.getMz() - FEATURECHART_MZ_PADDING_DOWN,
                                    selectedFeature.getMz() +
                                            Math.max(FEATURECHART_MIN_PEAKS, selectedFeature.getPeaks()) +
                                            FEATURECHART_MZ_PADDING_UP,
                                    selectedFeature.getScanFirst(), selectedFeature.getScanLast(),
                                    selectedFeature.getScanFirst(), selectedFeature.getScanLast(),
                                    selectedFeature.getMz(), 0, selectedFeature.getCharge());
                    pwsc.setGenerate3DChart(true);
                    pwsc.generateCharts();
                    ApplicationContext.setMessage("");

                    if (MultiChartDisplayPanel.getSingletonInstance().getNumCharts() > 0)
                        MultiChartDisplayPanel.createNewSingletonInstance();

                    pwsc.getContourPlot().displayInTab("3D");
                    pwsc.getIntensitySumChart().displayInTab("Intensity Sum");
                    try
                    {
                        File imageFile = TempFileManager.createTempFile("tempchart.png","asdf");
                        pwsc.saveChartToImageFile(imageFile);
                        new PanelWithBlindImageChart(imageFile, "Heatmap").displayInTab();
                    }
                    catch (IOException e) {}
                    TempFileManager.deleteTempFiles("asdf");
                }
            }
            );
            add(showZoomMenuItem);
        }
    }

    protected void addPropertyToModel(Object propertyName, Object propertyValue)
    {
        int numRows = _model.getRowCount();
        _model.setRowCount(numRows + 1);
        _model.setValueAt(propertyName, numRows, 0);
        _model.setValueAt(propertyValue, numRows, 1);
    }

    protected void clearProperties()
    {
        for (int i=_model.getRowCount()-1; i>=0; i--)
        {
            _model.setValueAt(null, i, 0);
            _model.setValueAt(null, i, 1);
        }
        _model.setRowCount(0);
    }


    public <E extends Map.Entry<?, ?>> void setProperties(Collection<E> c)
    {
        clearProperties();
        for (Object o : c)
        {
            Map.Entry entry = (Map.Entry) o;
            addPropertyToModel(entry.getKey(), entry.getValue());
        }
    }


    public void setProperties(Map<String, ?> m, java.util.List<String> keys)
    {
        if (null == keys)
        {
            keys = new ArrayList<String>();
            keys.addAll(m.keySet());
            Collections.sort(keys, new Comparator<String>()
            {
                public int compare(String a, String b)
                {
                    return (a.compareToIgnoreCase(b));
                }
            });
        }

        for (String key : keys)
        {
            Object value = m.get(key);
            addPropertyToModel(key, value);
        }
    }
}
