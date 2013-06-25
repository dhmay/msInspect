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
package org.fhcrc.cpl.viewer.commandline.modules;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class CollapseSpreadsheetRowsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CollapseSpreadsheetRowsCLM.class);

    protected String collapseColumnName = null;
    protected File inFile;
    protected File outFile;
    protected String separatorString = ";";
    protected List<String> otherColsToCollapse = new ArrayList<String>();
    protected boolean shouldCollapseToUnique = false;


    public CollapseSpreadsheetRowsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "collapsespreadsheetrows";
        mShortDescription = "collapse multiple rows in a spreadsheet with the same value for collapsecolumn";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true,"input spreadsheet"),
                        new StringArgumentDefinition("collapsecolumn", true, "column to split on"),
                        new FileToWriteArgumentDefinition("out", false, "output file"),
                        new StringArgumentDefinition("separator", false,
                                "String to split on", separatorString),
                        new StringListArgumentDefinition("othercolstocollapse", false,
                                "other columns to collapse (rather than keep one value). Any column NOT on this list " +
                                        "must have only one unique value, or an error is thrown."),
                        new BooleanArgumentDefinition("collapsetounique", false, 
                                "for other columns, collapse into unique values.", false),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFile = getUnnamedFileArgumentValue();
        collapseColumnName = getStringArgumentValue("collapsecolumn");
        separatorString = getStringArgumentValue("separator");
        outFile = getFileArgumentValue("out");
        otherColsToCollapse = getStringListArgumentValue("othercolstocollapse");

        if (otherColsToCollapse == null) 
            otherColsToCollapse = new ArrayList<String>();
        if (!otherColsToCollapse.isEmpty()) {
            ApplicationContext.infoMessage("collapsing other columns:");
            for (String col : otherColsToCollapse)
                ApplicationContext.infoMessage("\t" + col);
        }
        shouldCollapseToUnique = getBooleanArgumentValue("collapsetounique");
    }
    
    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PrintWriter outPW = null;
        try
        {
            outPW = new PrintWriter(outFile);

            TabLoader loader = new TabLoader(inFile);
            List<TabLoader.ColumnDescriptor> columns =
                    new ArrayList<TabLoader.ColumnDescriptor>();
            StringBuffer headerLine = new StringBuffer();
            boolean firstColumn = true;
            for (TabLoader.ColumnDescriptor column : loader.getColumns())
            {
                columns.add(column);
                if (firstColumn)
                    firstColumn = false;
                else
                    headerLine.append("\t");
                headerLine.append(column.name);
            }
            outPW.println(headerLine.toString());
            outPW.flush();

            Map<String, List<Map>> keyRowsMap = new HashMap<String, List<Map>>();

            Map[] rowsAsMaps = (Map[])loader.load();

            int lineNum = 0;
            for (Map row : rowsAsMaps)
            {
                lineNum++;
                Object collapseVal = row.get(collapseColumnName);
                if (collapseVal == null) {
                    throw new CommandLineModuleExecutionException("Null collapse value for line " + lineNum);
                }
                String key = collapseVal.toString();

                List<Map> rows = keyRowsMap.get(key);
                if (rows == null) {
                    rows = new ArrayList<Map>();
                    keyRowsMap.put(key, rows);
                }
                rows.add(row);
            }

            for (String key : keyRowsMap.keySet()) {
                StringBuffer lineBuf = new StringBuffer();
                List<Map> rows = keyRowsMap.get(key);
                boolean firstCol = true;
                for (TabLoader.ColumnDescriptor column : loader.getColumns()) {
                    String valString = "";

                    if (key.equals(collapseColumnName))
                        valString = key;
                    else {
                        if (shouldCollapseToUnique) {
                            Set<String> uniqueValsThisCol = new HashSet<String>();
                            for (Map row : rows) {
                                String val = "";
                                if (row.containsKey(column.name) && row.get(column.name) != null)
                                    val = row.get(column.name).toString();
                                uniqueValsThisCol.add(val);
                            }
                            boolean firstVal = true;
                            for (String val : uniqueValsThisCol) {
                                if (!firstVal)
                                    valString = valString + separatorString;
                                valString = valString + val;
                                firstVal = false;
                            }
                        } else {
                            List<String> allValsThisCol = new ArrayList<String>();
                            for (Map row : rows) {
                                String val = "";
                                if (row.containsKey(column.name) && row.get(column.name) != null)
                                    val = row.get(column.name).toString();
                                allValsThisCol.add(val);
                            }
    
                            if (otherColsToCollapse.contains(column.name)) {
                                boolean firstVal = true;
                                for (String val : allValsThisCol) {
                                    if (!firstVal)
                                        valString = valString + separatorString;
                                    valString = valString + val;
                                    firstVal = false;
                                }
                            }
                            else {
                                Set<String> uniqueVals = new HashSet<String> (allValsThisCol);
                                if (uniqueVals.size() > 1 && uniqueVals.contains("")) {
                                    uniqueVals.remove("");
                                }
                                if (uniqueVals.size() > 1)       {
                                    boolean firstVal = true;
                                    for (String val : uniqueVals) {
                                        if (!firstVal)
                                            valString = valString + separatorString;
                                            valString = valString + val;
                                    }
    //                                ApplicationContext.infoMessage("Vals:");
    //                                for (String val : uniqueVals)
    //                                    ApplicationContext.infoMessage("\t" + val);
    //                                throw new CommandLineModuleExecutionException("Key " + key + ", col " + column.name +
    //                                        ", vals: " + uniqueVals.size());
                                }
                                else {
                                    valString = uniqueVals.iterator().next();
                                }
                            }
                        }
                    }
    
                    lineBuf.append(firstCol ? "" : "\t");
                    lineBuf.append(valString);

                    firstCol = false;
                }
                outPW.println(lineBuf.toString());
                outPW.flush();

            }


        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException(e);
        }
        finally
        {
            if (outPW != null)
                outPW.close();
        }
    }

}
