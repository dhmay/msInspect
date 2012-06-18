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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.systemsbiology.jrap.stax.DataProcessingInfo;
import org.systemsbiology.jrap.stax.MSInstrumentInfo;
import org.systemsbiology.jrap.stax.MZXMLFileInfo;
import org.systemsbiology.jrap.stax.SoftwareInfo;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyEvent;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: Feb 28, 2005
 * Time: 3:09:11 PM
 */
public class MSRunInfoAction extends AbstractAction
	{
	JDialog _dialog = null;

	public MSRunInfoAction(String name)
		{
		super(name);
		}

	public void actionPerformed(ActionEvent e)
		{
		MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
		if (null == run)
			return;

		showProperties(run);
		//showPopup(run);
		}


	public void showPopup(MSRun run)
		{
		TableModel model = new RunTableModel(run);
		JTable table = new JTable(model);
		table.setEnabled(false);
		JScrollPane sp = new JScrollPane(table);
		JDialog dialog = new JDialog(ApplicationContext.getFrame(), "File Properties -- " + run.getFileName(), true);
		dialog.getContentPane().add(sp);

		ListenerHelper helper = new ListenerHelper(this);
		helper.addListener(table, "dialog_keyPressed");
		helper.addListener(dialog, "dialog_keyPressed");
		helper.addListener(dialog, "dialog_componentHidden");

		dialog.setSize(400,300);
		_dialog = dialog;
		dialog.setVisible(true);
		}


	public void dialog_keyPressed(KeyEvent e)
		{
		if (e.getKeyCode() == KeyEvent.VK_ESCAPE && null != _dialog)
			{
			_dialog.dispose();
			e.consume();
			return;
			}
		}


	public void dialog_componentHidden(ComponentEvent e)
		{
        _dialog = null;
		}


	static class RunTableModel extends DefaultTableModel
		{
        RunTableModel(MSRun run)
	        {
	        super(new Object[] {"Property", "Value"}, 0);

	        MZXMLFileInfo fileInfo = run.getHeaderInfo();
	        if (null == fileInfo)
		        return;

	        MSInstrumentInfo inst = fileInfo.getInstrumentInfo();
	        if (null != inst)
		        {
		        //add("<html><b>Instrument<b></html>","");
				add("Manufacturer", inst.getManufacturer());
				add("Detector", inst.getDetector());
				add("MassAnalyser", inst.getMassAnalyzer());
				add("Model", inst.getModel());
				add("Operator", valueOf(inst.getOperator()));
				add("Software", valueOf(inst.getSoftwareInfo().name));
				add("Ionization", inst.getIonization());
		        }

	        DataProcessingInfo data = fileInfo.getDataProcessing();
	        if (null != data) //  && (data.getCentroided() != -1 || data.getChargeDeconvoluted() != -1 || data.getDeisotoped() != -1))
		        {
		        //add("<html><b>Data Processing<b></html>","");
				add("Centroided", String.valueOf(data.getCentroided()));
				add("Charge Deconvoluted", String.valueOf(data.getChargeDeconvoluted()));
				add("Deisotoped", String.valueOf(data.getDeisotoped()));
				add("Intensity Cutoff", String.valueOf(data.getIntensityCutoff()));
		        add("Spot Integration", String.valueOf(data.getSpotIntegration()));
		        List<SoftwareInfo> soft = data.getSoftwareUsed();
		        if (null != soft)
			        {
			        String key = "Software";
			        for (int i=0 ; i<soft.size() ; i++)
				        {
				        add(key, valueOf(soft.get(i).name));
				        key = "";
				        }
			        }
		        }
	        }

		void add(String key, String value)
			{
			int row = getRowCount();
			setRowCount(row+1);
			setValueAt(key, row, 0);
			setValueAt(value, row, 1);
			}
		}


	void showProperties(MSRun run)
		{
		if (ApplicationContext.getFrame() instanceof WorkbenchFrame)
			((WorkbenchFrame)ApplicationContext.getFrame()).showPropertiesPane();
		ApplicationContext.setProperty(SharedProperties.SELECTED, run);
		}


	static String valueOf(Object o)
		{
		return null == o ? "" : String.valueOf(o);
		}
	}
