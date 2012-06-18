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
package org.fhcrc.cpl.toolbox.gui;


import javax.swing.*;
import javax.swing.table.TableModel;
import javax.swing.table.TableColumnModel;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;


/**
 * Utilities for Swing objects.  Possibly break this up if it gets big
 */
public class SwingUtils<T,V>
{
    /**
     * Save the contents of a table to a TSV file
     * Note:  uses toString() on the header cells as well as the data cells.  If you've got funny columns,
     * expect funny behavior
     * @param table
     * @param outFile
     * @throws IOException
     */
    public static void SaveTableAsTSV(JTable table, File outFile)
            throws IOException
    {
        PrintWriter outPW = new PrintWriter(outFile);        

        TableModel tableModel = table.getModel();
        TableColumnModel columnModel = table.getColumnModel();

        StringBuffer headerLineBuf = new StringBuffer();
        for (int i=0; i<columnModel.getColumnCount(); i++)
        {
            if (i>0)
                headerLineBuf.append("\t");
            headerLineBuf.append(columnModel.getColumn(i).getHeaderValue().toString());
        }
        outPW.println(headerLineBuf.toString());
        outPW.flush();
        for (int i=0; i<tableModel.getRowCount(); i++)
        {
            StringBuffer lineBuf = new StringBuffer();
            for (int j=0; j<tableModel.getColumnCount(); j++)
            {
                if (j>0)
                    lineBuf.append("\t");
                lineBuf.append(tableModel.getValueAt(i,j).toString());
            }
            outPW.println(lineBuf.toString());
            outPW.flush();
        }
        outPW.close();
    }


}
