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
public class DuplicateSpreadsheetRowsCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(DuplicateSpreadsheetRowsCLM.class);

    protected String splitColumnName = null;
    protected File inFile;
    protected File outFile;
    protected String stringToSplitOn = ";";
    protected List<String> otherColsToSplit = new ArrayList<String>();


    public DuplicateSpreadsheetRowsCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "duplicatespreadsheetrows";
        mShortDescription = "Create multiple rows for rows with multiple values in " +
                "splitcolumn separated by stringtospliton, duplicating all other values";
        mHelpMessage = "duplicatespreadsheetrows";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true,"input spreadsheet"),
                        new StringArgumentDefinition("splitcolumn", true, "column to split on"),
                        new FileToWriteArgumentDefinition("out", false, "output file"),
                        new StringArgumentDefinition("stringtospliton", false,
                                "String to split on", stringToSplitOn),
                        new StringListArgumentDefinition("othercolstosplit", false,
                                "other columns to split (rather that duplicate).  All of these must have same cardinality as splitcolumn, per row."),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFile = getUnnamedFileArgumentValue();
        splitColumnName = getStringArgumentValue("splitcolumn");
        stringToSplitOn = getStringArgumentValue("stringtospliton");
        outFile = getFileArgumentValue("out");
        otherColsToSplit = getStringListArgumentValue("othercolstosplit");

        if (otherColsToSplit != null && !otherColsToSplit.isEmpty()) {
            ApplicationContext.infoMessage("Splitting other columns:");
            for (String col : otherColsToSplit)
                ApplicationContext.infoMessage("\t" + col);
        }
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

            Map[] rowsAsMaps = (Map[])loader.load();

            for (Map row : rowsAsMaps)
            {
                String splitColumnValue = (String) row.get(splitColumnName);
                String[] splitColumnVals = null;
                if (splitColumnValue == null)
                    splitColumnVals = new String[] {null};
                else
                    splitColumnVals = splitColumnValue.split(stringToSplitOn);

                Map<String, String[]> otherColumnVals = new HashMap<String, String[]>();
                if (splitColumnValue != null && otherColsToSplit != null && !otherColsToSplit.isEmpty()) {
                    for (String otherColumn : otherColsToSplit) {
                        String otherColString = (String) row.get(otherColumn);
                        if (otherColString == null) {
                            ApplicationContext.infoMessage("Other column " +
                                    otherColumn + " had null value for split column value " + splitColumnValue);
                            otherColumnVals.put(otherColumn, new String[splitColumnVals.length]);
                        }
                        else {
                            String[] otherVals = otherColString.split(stringToSplitOn);
                            if (otherVals.length != splitColumnVals.length) {
                                ApplicationContext.infoMessage("Other column " +
                                        otherColumn + " had wrong number of values (" + otherColumnVals +
                                        ") for split column value " + splitColumnValue);
                                otherColumnVals.put(otherColumn, new String[splitColumnVals.length]);

                            }
                            else
                                otherColumnVals.put(otherColumn, otherVals);
                        }
                    }
                }

                for (int i=0; i<splitColumnVals.length; i++)
                {
                    String splitColumnVal = splitColumnVals[i];
                    firstColumn = true;

                    StringBuffer line = new StringBuffer();

                    for (TabLoader.ColumnDescriptor column : columns)
                    {
                        if (firstColumn)
                            firstColumn = false;
                        else line.append("\t");
                        if (column.name.equals(splitColumnName))
                            line.append(splitColumnVal == null ? "" : splitColumnVal);
                        else if (otherColumnVals.containsKey(column.name))
                            line.append(otherColumnVals.get(column.name)[i]);
                        else
                            line.append(null == row.get(column.name) ? "" : row.get(column.name));
                    }
                    outPW.println(line);

                }
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
