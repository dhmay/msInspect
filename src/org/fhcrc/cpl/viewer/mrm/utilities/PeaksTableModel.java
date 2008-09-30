/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.fhcrc.cpl.viewer.mrm.Utils;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import javax.swing.table.AbstractTableModel;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Oct 18, 2007
 * Time: 12:51:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class PeaksTableModel extends AbstractTableModel {
    private static String columnNames[] =  {
       MRMDialog.peaksData.Accept.toString(),
       MRMDialog.peaksData.Peptide.toString(),
       MRMDialog.peaksData.Precursor.toString(),
       MRMDialog.peaksData.Daughter.toString(),
       MRMDialog.peaksData.CoStart.toString(),
       MRMDialog.peaksData.CoEnd.toString(),
       MRMDialog.peaksData.CoDelta.toString(),
       MRMDialog.peaksData.AUC.toString(),
       MRMDialog.peaksData.MaxPeak.toString(),
       MRMDialog.peaksData.MidTime.toString(),
       MRMDialog.peaksData.Quality.toString(), 
       MRMDialog.peaksData.Label.toString(),
       MRMDialog.peaksData.Code.toString(),
       MRMDialog.peaksData.LHRatio.toString(),     
       MRMDialog.peaksData.Comment.toString()
    };

     public Object[][] data = new Object[1][columnNames.length];

     public int getColumnCount() {
         return columnNames.length;
     }

     public int getRowCount() {
         return data.length;
     }

     public String getColumnName(int col) {
         return columnNames[col];
     }

     public Object getValueAt(int row, int col) {
         return data[row][col];
     }

     public Class getColumnClass(int c) {
         for(MRMDialog.peaksData pd: MRMDialog.peaksData.values()) {
             if(pd.ordinal() == c)
                 return pd.colClass;
         }
         return String.class;
     }

     public boolean isCellEditable(int row, int col) {
         if (col == MRMDialog.peaksData.Accept.colno) {
              return data[row][MRMDialog.peaksData.Daughter.colno] != null;
         }
         if (col != MRMDialog.peaksData.Precursor.colno && col != MRMDialog.peaksData.Daughter.colno && col != MRMDialog.peaksData.CoDelta.colno) {
             return true;
         } else {
             return false;
         }
     }

     public void setValueAt(Object value, int row, int col) {
         data[row][col] = value;
         fireTableCellUpdated(row, col);
     }

     public boolean saveSaveModelAsTSV(File out){
        BufferedWriter br = null;
        try {
            br = new BufferedWriter(new FileWriter(out));
            for(int i = 0; i<columnNames.length; i++) {
                if(columnNames[i].equalsIgnoreCase("Quality") && Utils.qualColIsEmpty()) continue;
                br.write(columnNames[i]);
                if(i < (columnNames.length-1)) {
                    br.write("\t");
                } else {
                    br.write("\n");
                }
            }
            for(int rows = 0; rows < data.length; rows++) {
                for(MRMDialog.peaksData pd: MRMDialog.peaksData.values()) {
                   if(data[rows][pd.colno] != null){
                      if(pd == MRMDialog.peaksData.Quality && Utils.qualColIsEmpty()) continue;
                      br.write(pd.SaveFileFormat(data[rows][pd.colno]));
                   }
                   if(pd.ordinal() < (columnNames.length-1)) {
                        br.write("\t");
                   } else {
                        br.write("\n");
                   }
                }
            }
            br.flush();
            br.close();
        } catch (Exception e) {
           ApplicationContext.infoMessage("Can't write table to file: "+e);
           return false;
        }
        return true;
     }

}

