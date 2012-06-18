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
package org.fhcrc.cpl.toolbox.filehandler;

import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * for writing tab files
 */
public  class TabWriter
{
    protected String[] columnNames;
    protected List<Map<String,Object>> rows;
    protected File outFile;



    public TabWriter(String[] columnNames)
    {
        this.columnNames = columnNames;
        rows = new ArrayList<Map<String,Object>>();
    }

    public TabWriter(String[] columnNames, Map<String,Object>[] rows, File outFile)
    {
        this(columnNames);

        this.rows = new ArrayList<Map<String,Object>>(rows.length);
        for (Map<String,Object> row : rows)
            this.rows.add(row);
        this.outFile = outFile;
    }

    public TabWriter(String[] columnNames, List<Map<String,Object>> rows, File outFile)
    {
        this(columnNames);
        this.rows = rows;
        this.outFile = outFile;
    }

    public void addRow(Map<String,Object> newRow)
    {
        rows.add(newRow);
    }


    public void write()
            throws IOException
    {
        PrintWriter outPW = null;
        try
        {
            outPW = new PrintWriter(outFile);
            StringBuffer headerLine = new StringBuffer();
            for (int i=0; i<columnNames.length; i++)
            {
                if (i>0)
                    headerLine.append("\t");
                headerLine.append(columnNames[i]);

            }
            outPW.println(headerLine.toString());
            outPW.flush();

            for (Map<String,Object> row : rows)
            {
                StringBuffer line = new StringBuffer();
                for (int i=0; i<columnNames.length; i++)
                {
                    if (i>0)
                        line.append("\t");
                    if (row.containsKey(columnNames[i]))
                        line.append(row.get(columnNames[i]).toString());
                }
                outPW.println(line.toString());
                outPW.flush();
            }
        }
        catch (IOException e)
        {
            throw e;
        }
        finally
        {
            if (outPW != null)
                outPW.close();
        }
    }


    public List<Map<String, Object>> getRows()
    {
        return rows;
    }

    public void setRows(List<Map<String, Object>> rows)
    {
        this.rows = rows;
    }

    public String[] getColumnNames()
    {
        return columnNames;
    }

    public void setColumnNames(String[] columnNames)
    {
        this.columnNames = columnNames;
    }

    public File getOutFile()
    {
        return outFile;
    }

    public void setOutFile(File outFile)
    {
        this.outFile = outFile;
    }
}