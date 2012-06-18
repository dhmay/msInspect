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
package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.gui.MRMDialog;
import org.fhcrc.cpl.viewer.mrm.Utils;
import org.fhcrc.cpl.viewer.mrm.MRMTransition;
import org.fhcrc.cpl.viewer.mrm.MRMDaughter;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import javax.swing.table.AbstractTableModel;
import java.io.*;

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

     public boolean restoreModelFromTSV(File inf) {
         BufferedReader br = null;
         try {
            //check to see that the headers are correct
            br = new BufferedReader(new FileReader(inf));
            String line = br.readLine();
            String headers[] = line.split("\\t");
            int i = 0;
            for(String currentHeader: columnNames) {
                if(!currentHeader.equals(headers[i])){
                    if(currentHeader.equalsIgnoreCase("Quality") && Utils.qualColIsEmpty()) continue;
                    throw new Exception("Headers do not match");
                } else i++;
            }
            //Check that TSV file and model have same transitions
            for(i = 0; i<getRowCount(); i++) {
                line = br.readLine();
                String tokens[] = line.split("\\t");
                //if we're processing a Precursor line
                if(data[i][MRMDialog.peaksData.Daughter.colno] == null) {
                    //is TSV also a precursor line?
                    if(!tokens[3].equals(""))throw new Exception("Transitions are not the same in this file as in current data");
                    float curPrecursor = ((MRMTransition)data[i][MRMDialog.peaksData.Precursor.colno]).getPrecursorMz();
                    float curTSVPrecursor = Float.parseFloat(tokens[2]);
                    if(Math.abs(curPrecursor-curTSVPrecursor) > 0.00001F)
                       throw new Exception("TSV is not from the same set of transitions");
                } else { //product line
                    float curProduct = ((MRMDaughter)data[i][MRMDialog.peaksData.Daughter.colno]).getMeanMz();
                    float curTSVProduct = Float.parseFloat(tokens[3]);
                    if(Math.abs(curProduct-curTSVProduct) > 0.00001F)
                        throw new Exception("TSV is not from the same set of transitions");
                    
                }
            }
            //At this point we assume the table and the TSV file have the same
            //transitions, otherwise an exception would have been thrown.  So
            //its "safe" to update the table from the TSV file.
            br.close();
            br = new BufferedReader(new FileReader(inf));
            line = br.readLine();
            boolean qualityInTSVfile = line.contains("Quality");
            for(i = 0; i<getRowCount(); i++) {
               line = br.readLine();
               String tokens[] = line.split("\\t");
               //insert phony quality line if it doesn't exist
               if(!qualityInTSVfile && tokens.length >= 10) {
                   String newTokens[] = new String[tokens.length+1];
                   for(int k=0; k<=9; k++) newTokens[k] = tokens[k];
                   newTokens[10] = "-1";
                   for(int k = 11; k < tokens.length+1; k++) newTokens[k] = tokens[k-1];
                   tokens = newTokens;
               }

               //if we're processing a Precursor line
               MRMDialog.peaksTable.setRowSelectionInterval(i,i); 

                  for(int j=0; j<tokens.length;j++) {
                      switch(j){
                         //accept
                         case 0: if(tokens[0].equals("")) {
                                    setValueAt(null,i,MRMDialog.peaksData.Accept.colno);
                                 } else {
                                    setValueAt(new Boolean(tokens[0]),i,MRMDialog.peaksData.Accept.colno);
                                 }
                                 break;
                         //peptide
                         case 1: setValueAt(tokens[j],i,MRMDialog.peaksData.Peptide.colno);
                                 break;
                         //precursor and product masses
                         case 2:
                         case 3: break;
                         //Coelution start time
                         case 4: setValueAt(Float.parseFloat(tokens[4]),i,MRMDialog.peaksData.CoStart.colno);
                                 break;
                         //Coelution end time
                         case 5: setValueAt(Float.parseFloat(tokens[5]),i,MRMDialog.peaksData.CoEnd.colno);
                                 break;
                         //Coelution difference
                         case 6: setValueAt(Float.parseFloat(tokens[6]),i,MRMDialog.peaksData.CoDelta.colno);
                                 break;
                         //AUC  (is it the real current AUC?)
                         case 7: setValueAt(Float.parseFloat(tokens[7]),i,MRMDialog.peaksData.AUC.colno);
                                 break;
                         //MaxPeak (is it the real current MAXPEAK?)
                         case 8: setValueAt(Float.parseFloat(tokens[8]),i,MRMDialog.peaksData.MaxPeak.colno);
                                 break;
                         //MidTime (remember to move mid-time pointer)
                         case 9: setValueAt(Float.parseFloat(tokens[9]),i,MRMDialog.peaksData.MidTime.colno);
                                 break;
                         //Quality
                         case 10:setValueAt(Float.parseFloat(tokens[10]),i,MRMDialog.peaksData.Quality.colno);
                                 break;
                         //Label
                         case 11: setValueAt(tokens[11],i,MRMDialog.peaksData.Label.colno);
                                 break;
                         //Code
                         case 12: if(tokens[12].equals("")) {
                                     setValueAt(null,i,MRMDialog.peaksData.Code.colno);
                                  } else {
                                     setValueAt(Integer.parseInt(tokens[12]),i,MRMDialog.peaksData.Code.colno);
                                  }   
                                 break;
                         //LHRatio
                         case 13: if(tokens[13].equals("")) {
                                     setValueAt(null,i,MRMDialog.peaksData.LHRatio.colno);
                                  } else {
                                     setValueAt(Float.parseFloat(tokens[13]),i,MRMDialog.peaksData.LHRatio.colno);
                                  } 
                                  break;
                         //Comment
                         case 14: setValueAt(tokens[14],i,MRMDialog.peaksData.Comment.colno);
                                  break;
                         default:  ApplicationContext.infoMessage("Too many columns in TSV data!");

                      }
                  }
            }

         
         } catch (Exception e) {
            ApplicationContext.infoMessage("File is not compatible with program or data: "+e);
            e.printStackTrace();
            return false;
         } finally {
             try{br.close();}catch (Exception ee) {}
         }
         MRMDialog.peaksTable.setRowSelectionInterval(1,1); 
         return true;
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

