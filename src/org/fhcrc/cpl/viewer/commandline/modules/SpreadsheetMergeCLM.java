/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;


import java.io.*;
import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class SpreadsheetMergeCLM extends BaseCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SpreadsheetMergeCLM.class);

    protected String mergeColumnName = null;
    protected File[] inFiles;
    protected File outFile;
    protected File compareOutFile;
    protected File outUnique2File;
    protected String file2ColumnName = null;

    protected boolean keepAllFile1Values = false;



    protected String plotColumnName =  null;

    protected boolean plotLog = false;

    public SpreadsheetMergeCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "spreadsheetmerge";
        mShortDescription = "merge spreadsheets";
        mHelpMessage = "merge spreadsheets";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true,"input spreadsheets"),
                        createStringArgumentDefinition("mergecolumn", true, "column to merge on"),
                        createFileToWriteArgumentDefinition("out", false, "output file"),
                        createStringArgumentDefinition("plotcolumn", false, "column to plot, one vs. the other"),
                        createStringArgumentDefinition("file2column", false,
                                "column to add from the second file.  If not specified, all columns added"),
                        createBooleanArgumentDefinition("plotlog", false, "Plot in log scale", false),
                        createFileToWriteArgumentDefinition("compareout", false,
                                "output file for comparing values of plotcolumn"),
                        createFileToWriteArgumentDefinition("outunique2file", false,
                                "output file for rows unique to the second spreadsheet"),                        

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFiles = getUnnamedSeriesFileArgumentValues();
        if (inFiles.length < 2)
            throw new ArgumentValidationException("Must specify at least two input files");
        mergeColumnName = getStringArgumentValue("mergecolumn");
        plotColumnName = getStringArgumentValue("plotcolumn");
        file2ColumnName = getStringArgumentValue("file2column");
        if (file2ColumnName != null)
        {
            ApplicationContext.infoMessage("File 2 column specified.  Will keep ALL rows from file 1, " +
                    "annotating those with file2column with the appropriate vaelu");
            keepAllFile1Values = true;
        }
        else
            ApplicationContext.infoMessage("Will keep only those rows with mergecolumn common to all files");

        plotLog = getBooleanArgumentValue("plotlog");
        outFile = getFileArgumentValue("out");
        compareOutFile = getFileArgumentValue("compareout");
        outUnique2File = getFileArgumentValue("outunique2file");

    }

    protected Map<String,Map> loadRowsFromFile(TabLoader loader)
            throws IOException
    {
        Map[] rowsAsMaps = (Map[])loader.load();

        Map<String,Map> result = new HashMap<String,Map>();

        try
        {
            for (Map row : rowsAsMaps)
            {
                String key = (String) row.get(mergeColumnName);
                if (key != null)
                    result.put(key,row);
            }
        }
        catch (ClassCastException e)
        {
            throw new IOException("Error: All mergecolumn entries must be Strings");
        }
        return result;
    }

    String createFileLine(String key, List<TabLoader.ColumnDescriptor>[] columnsAllFiles,
                          Map[] rowMapsAllFiles)
    {
        StringBuffer resultBuf = new StringBuffer(key);
        for (int i=0; i<columnsAllFiles.length; i++)
        {
            Map rowMapThisFile = rowMapsAllFiles[i];
            List<TabLoader.ColumnDescriptor> columnsThisFile = columnsAllFiles[i];
            for (TabLoader.ColumnDescriptor column : columnsThisFile)
            {
                if (!mergeColumnName.equals(column.name))
                {
                    resultBuf.append("\t");
                    if (rowMapThisFile != null && rowMapThisFile.get(column.name) != null)
                        resultBuf.append(rowMapThisFile.get(column.name));
                }
            }
        }
        return resultBuf.toString();
    }

    /**
     * do the actual work
     */
    public void execute() throws CommandLineModuleExecutionException
    {
        PrintWriter outPW = null;
        try
        {
            if (outFile != null)
                outPW = new PrintWriter(outFile);

            List<String> combinedColumns = new ArrayList<String>();
            combinedColumns.add(mergeColumnName);

            StringBuffer headerLine = new StringBuffer(mergeColumnName);
            TabLoader[] tabLoaders = new TabLoader[inFiles.length];
            List<TabLoader.ColumnDescriptor>[] columnsAllFiles =
                    new List[inFiles.length];

            for (int i=0; i<inFiles.length; i++)
            {
                File inFile = inFiles[i];
                System.err.println("Loading file " + inFile.getName());
                TabLoader loader = new TabLoader(inFile);
                List<TabLoader.ColumnDescriptor> columnsThisFile =
                        new ArrayList<TabLoader.ColumnDescriptor>();
                if (i==1 && file2ColumnName != null)
                {

                    for (TabLoader.ColumnDescriptor column : loader.getColumns())
                    {
                        if (column.name.equals(file2ColumnName))
                            columnsThisFile.add(column);
                    }
                    headerLine.append("\t" + file2ColumnName);
                }
                else
                {
                    for (TabLoader.ColumnDescriptor column : loader.getColumns())
                    {
                        columnsThisFile.add(column);
                        if (!mergeColumnName.equals(column.name))
                            headerLine.append("\t" + column.name);
                    }
                }
                columnsAllFiles[i] = columnsThisFile;
                tabLoaders[i] = loader;
            }

            if (outFile != null)
            {
                outPW.println(headerLine.toString());
                outPW.flush();
            }

            Map<String,Map>[] rowMaps = new Map[inFiles.length];
            for (int i=0; i<inFiles.length; i++)
            {
                rowMaps[i] = loadRowsFromFile(tabLoaders[i]);
            }

            Set<String> keysInAllFiles = new HashSet<String>();

            for (String key : rowMaps[0].keySet())
            {
                boolean notFoundSomewhere = false;
                for (int i=1; i<rowMaps.length; i++)
                {
                    if (!rowMaps[i].containsKey(key))
                    {
                        notFoundSomewhere = true;
                        continue;
                    }
                }
                if (!notFoundSomewhere)
                    keysInAllFiles.add(key);

            }
            ApplicationContext.infoMessage("Rows in common: " + keysInAllFiles.size());

            Set<String> keysToWrite = keysInAllFiles;
            if (keepAllFile1Values)
                keysToWrite = rowMaps[0].keySet();


            for (String key : keysToWrite)
            {
                    Map[] mapsAllFiles = new Map[rowMaps.length];
                    for (int j=0; j<rowMaps.length; j++)
                        mapsAllFiles[j] = rowMaps[j].get(key);

                    String line = createFileLine(key, columnsAllFiles,
                            mapsAllFiles);
                    if (outFile != null)
                    {
                        outPW.println(line);
                        outPW.flush();
                    }
            }


            if (outUnique2File != null)
            {
                PrintWriter unique2OutWriter = new PrintWriter(outUnique2File);
                StringBuffer headerLineBuf = new StringBuffer(mergeColumnName);
                for (TabLoader.ColumnDescriptor column : columnsAllFiles[1])
                {
                    if (!mergeColumnName.equals(column.name))
                        headerLineBuf.append("\t" + column.name);
                }
                unique2OutWriter.println(headerLineBuf);
                for (String key : rowMaps[1].keySet())
                {
                    List<TabLoader.ColumnDescriptor>[] colArray = new List[] { columnsAllFiles[1] };
                    Map[] colMap = new Map[] { rowMaps[1].get(key) };

                    String line = createFileLine(key, colArray,
                            colMap);
                    if (outUnique2File != null)
                    {
                        unique2OutWriter.println(line);
                        unique2OutWriter.flush();
                    }
                }
                unique2OutWriter.close();
                ApplicationContext.infoMessage("Wrote lines unique to file 2 in " + outUnique2File.getAbsolutePath());
            }

            //first two files only
            if (plotColumnName != null)
            {
                List<Float> values1 = new ArrayList<Float>();
                List<Float> values2 = new ArrayList<Float>();
                List<String> commonKeys = new ArrayList<String>();


                Map<String,Map> rowMaps1 = rowMaps[0];
                Map<String,Map> rowMaps2 = rowMaps[1];


                for (String key : rowMaps1.keySet())
                {
                    if (rowMaps2.containsKey(key))
                    {
                        Object o1 = rowMaps1.get(key).get(plotColumnName);
                        Object o2 = rowMaps2.get(key).get(plotColumnName);

                        if (o1 == null || o2 == null)
                            continue;

                        Double value1;
                        Double value2;



                        try
                        {
                            value1 = (Double) o1;
                            value2 = (Double) o2;
                        }
                        catch(ClassCastException e)
                        {
                            try
                            {
                                  value1 = (double) ((Integer) o1);
                                  value2 = (double) ((Integer) o2);
                            }
                            catch(ClassCastException e2)
                            {
                                try
                                {
                                value1 = (double) ((Float) o1);
                                value2 = (double) ((Float) o2);
                                }
                                catch(ClassCastException e3)
                                {
                                    try
                                    {
                                        value1 = (double) (Float.parseFloat ((String) o1));
                                        value2 = (double) (Float.parseFloat  ((String) o2));
                                    }
                                    catch (ClassCastException e4)
                                    {
                                        ApplicationContext.infoMessage("Crap!  Can't process value " +
                                                rowMaps1.get(key).get(plotColumnName) + " or " + rowMaps2.get(key).get(plotColumnName));
                                        throw new CommandLineModuleExecutionException(e3);
                                    }
                                }
                            }
                        }
                        if (value1 != null && value2 != null)
                        {
                            float displayValue1 = value1.floatValue();
                            float displayValue2 =  value2.floatValue();
                            if (plotLog)
                            {
                                if (displayValue1 == 0)
                                    displayValue1 += 0.000001;
                                if (displayValue2 == 0)
                                    displayValue2 += 0.000001;
                                displayValue1 = (float) Math.log(displayValue1);
                                displayValue2 = (float) Math.log(displayValue2);
                            }
                            if (!Float.isInfinite(displayValue1) && !Float.isInfinite(displayValue2) &&
                                    !Float.isNaN(displayValue1) && !Float.isNaN(displayValue2))
                            {
                                values1.add(displayValue1);
                                values2.add(displayValue2);
                                commonKeys.add(key);
                            }
                        }
                    }
                }
                ApplicationContext.infoMessage("Rows in common and plottable: " + values1.size());
                PanelWithScatterPlot pwsp = new PanelWithScatterPlot(values1, values2, plotColumnName);
                pwsp.setAxisLabels("File 1","File 2");
                pwsp.displayDialog("Corresponding values");

                if (compareOutFile != null)
                {
                    PrintWriter compareOutWriter = new PrintWriter(compareOutFile);
                    compareOutWriter.println(mergeColumnName + "\t" + plotColumnName + "_1\t" +
                            plotColumnName + "_2");
                    for (int i=0; i<values1.size(); i++)
                    {
                        compareOutWriter.println(commonKeys.get(i) + "\t" + values1.get(i) +
                                "\t" + values2.get(i));
                        compareOutWriter.flush();
                    }
                    compareOutWriter.close();
                }
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