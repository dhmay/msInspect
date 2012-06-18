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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;
import org.apache.commons.lang.StringUtils;


import java.io.*;
import java.util.*;


/**
 * This is an unwieldy jumble of code to help merge together two spreadsheets based on shared values in the
 * 'mergecolumn' column.  In several different ways.
 */
public class SpreadsheetMergeCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(SpreadsheetMergeCLM.class);

    protected String mergeColumnName = null;
    protected File[] inFiles;
    protected File outFile;
    protected File compareOutFile;
    protected File outUnique2File;
    protected String file2ColumnName = null;
    protected Set<String> valuesToTrack = null;
    protected boolean keepAllFile1Values = false;
    protected boolean keepAllValuesAllFiles = false;

    protected boolean multipleMergeColumnValuesFirstFile = false;

    protected String newColumnName = null;
    protected String presenceAnnotation = null;

    protected String plotColumnName =  null;
    protected String plotColumnNameFile2 =  null;


    protected boolean plotLog = false;

    protected boolean isCaseSensitive = true;

    public SpreadsheetMergeCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "spreadsheetmerge"; 
        mShortDescription = "merge spreadsheets";
        mHelpMessage = "This is for merging and comparing spreadsheets based on the value in some column " +
                "('mergecolumn').";
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedSeriesFileArgumentDefinition(true,"input spreadsheets"),
                        new StringArgumentDefinition("mergecolumn", true, "column to merge on"),
                        new FileToWriteArgumentDefinition("out", false, "output file"),
                        new StringArgumentDefinition("plotcolumn", false, "column to plot, one vs. the other"),
                        new StringArgumentDefinition("file2plotcolumn", false, "column to plot from file 2, if different from file 1"),
                        new StringArgumentDefinition("file2column", false,
                                "column to add from the second file.  If not specified, all columns added"),
                        new BooleanArgumentDefinition("plotlog", false, "Plot in log scale", false),
                        new StringArgumentDefinition("presenceannotation", false,
                                "Rather than adding columns from second file, add this string as the value for all " +
                                        "matching mergecolumn rows, as a new column with name 'newcolname'"),
                        new StringArgumentDefinition("newcolname", false,
                                "New column name, for presence annotations"),
                        new FileToWriteArgumentDefinition("compareout", false,
                                "output file for comparing values of plotcolumn"),
                        new FileToWriteArgumentDefinition("outunique2file", false,
                                "output file for rows unique to the second spreadsheet"),
                        new BooleanArgumentDefinition("keepallfile1values", false,
                                "Keep all values from the first file, even if they don't occur in other files?",
                                keepAllFile1Values),
                        new BooleanArgumentDefinition("keepallvaluesallfiles", false,
                                "Keep all values from the all files, even if they don't occur in other files?",
                                keepAllValuesAllFiles),
                        new BooleanArgumentDefinition("multiplemergecolvaluesfirstfile", false,
                                "check for multiple merge-column values in the first file, separated by ';'",
                                multipleMergeColumnValuesFirstFile),
                        new FileToReadArgumentDefinition("trackvaluesfile", false,
                                "File containing values (one per line, everything before whitespace) to track") ,
                        new BooleanArgumentDefinition("casesensitive", false, "Case-sensitive match on merge column? " +
                                "(if false, all merge column values will be lowercased)",isCaseSensitive)
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
        plotColumnNameFile2 = getStringArgumentValue("file2plotcolumn");
        if (plotColumnName != null && plotColumnNameFile2 == null)
            plotColumnNameFile2 = plotColumnName;

        isCaseSensitive = getBooleanArgumentValue("casesensitive");

        file2ColumnName = getStringArgumentValue("file2column");

        if (hasArgumentValue("newcolname"))
        {
            assertArgumentPresent("presenceannotation", "newcolname");
            assertArgumentAbsent("file2column", "newcolname");
            newColumnName = getStringArgumentValue("newcolname");
            presenceAnnotation = getStringArgumentValue("presenceannotation");
        }
        else
            assertArgumentAbsent("presenceannotation");

        File trackValuesFile = getFileArgumentValue("trackvaluesfile");
        if (trackValuesFile != null)
            valuesToTrack = new HashSet<String>(readOneStringPerLine(trackValuesFile));


        multipleMergeColumnValuesFirstFile = getBooleanArgumentValue("multiplemergecolvaluesfirstfile");

        keepAllFile1Values = getBooleanArgumentValue("keepallfile1values");
        keepAllValuesAllFiles = getBooleanArgumentValue("keepallvaluesallfiles");
        if (!keepAllFile1Values && hasArgumentValue("keepallfile1values") && keepAllValuesAllFiles)
            throw new ArgumentValidationException("Keep-values arguments don't make sense");

        if (file2ColumnName != null)
        {
            ApplicationContext.infoMessage("File 2 column specified.  Will keep ALL rows from file 1, " +
                    "annotating those with file2column with the appropriate value");
            keepAllFile1Values = true;
        }
        else
            ApplicationContext.infoMessage("Will keep only those rows with mergecolumn common to all files");

        plotLog = getBooleanArgumentValue("plotlog");
        outFile = getFileArgumentValue("out");
        compareOutFile = getFileArgumentValue("compareout");
        outUnique2File = getFileArgumentValue("outunique2file");
    }

    protected Map<String,Map> mapRowsByMergeCol(Map[] rowsAsMaps)
            throws IOException
    {


        Map<String,Map> result = new HashMap<String,Map>();

        for (Map row : rowsAsMaps)
        {
            Object key = row.get(mergeColumnName);
            if (key != null)
            {
////special-purpose junk
//                try{
//                    float newIntensity = Float.parseFloat(row.get("intensity").toString());
//                    float oldIntensity =  Float.parseFloat(result.get(key.toString()).get("intensity").toString());
//                    if (oldIntensity > newIntensity)
//                    {
//                        row = result.get(key.toString());
//                    }
//                }
//                catch (Exception e) {}

                String keyString = key.toString();
                if (!isCaseSensitive)
                    keyString = keyString.toLowerCase();
                result.put(keyString,row);
            }
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
                //forcing reading the first file line as a header line, since we're dependent on column names
                TabLoader loader = new TabLoader(new FileReader(inFile), true);
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
                else if (i==1 && newColumnName != null)
                {
                    headerLine.append("\t" + newColumnName);
                    columnsThisFile.add(new TabLoader.ColumnDescriptor(newColumnName, String.class));
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
                List<String> columnNamesThisFile = new ArrayList<String>();
                for (TabLoader.ColumnDescriptor column : columnsThisFile)
                    columnNamesThisFile.add(column.name);
                if (!columnNamesThisFile.contains(mergeColumnName))
                    throw new CommandLineModuleExecutionException("File " + inFile.getAbsolutePath() +
                            " does not contain column " + mergeColumnName);
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
                rowMaps[i] = mapRowsByMergeCol((Map[])new TabLoader(inFiles[i]).load());
                //Replace map with  presence annotations if that's what we're doing
                if (i == 1 && newColumnName != null)
                {
                    Map<String, Map> annotRowMap = new HashMap<String, Map>();
                    for (String key : rowMaps[i].keySet())
                    {
                        Map<String, String> keyMap = new HashMap<String, String>();
                        keyMap.put(newColumnName, presenceAnnotation);
                        annotRowMap.put(key, keyMap);
                    }
                    rowMaps[i] = annotRowMap;
                }
                ApplicationContext.infoMessage("Loaded " + rowMaps[i].size() + " rows from file " + (i+1));
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
                        break;
                    }
                }
                if (!notFoundSomewhere)
                    keysInAllFiles.add(key);
            }
            ApplicationContext.infoMessage("Rows in common: " + keysInAllFiles.size());

            Set<String> keysToWrite = new HashSet<String>();

            if (keepAllValuesAllFiles)
            {
                for (Map<String,Map> rowMap : rowMaps)
                    keysToWrite.addAll(rowMap.keySet());
            }
            else if (keepAllFile1Values)
                keysToWrite = rowMaps[0].keySet();
            else
            {
                keysToWrite = keysInAllFiles;
            }


            for (String key : keysToWrite)
            {
//if (key.equals("IPI00115660")) System.err.println("***!");
                Map[] mapsAllFiles = new Map[rowMaps.length];
                for (int j=0; j<rowMaps.length; j++)
                {
                    mapsAllFiles[j] = rowMaps[j].get(key);
//if (key.equals("IPI00115660")) System.err.println("\t" + j + ", " + mapsAllFiles[j].get("geommean_defghi"));                    
                    //if we don't find it, and we're supposed to split up keys by ";", do so
                    if (mapsAllFiles[j] == null && j > 0 && multipleMergeColumnValuesFirstFile &&
                            key.contains(";"))
                    {
                        for (String partKey : key.split(";"))
                        {
                            if (rowMaps[j].containsKey(partKey))
                            {
                                mapsAllFiles[j] = rowMaps[j].get(partKey);
ApplicationContext.infoMessage("Split up multi-key " + key + ", found match for " + partKey);
                                break;
                            }
                        }
                    }
                }

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
                List<Float> plotColumnUnique2Values = new ArrayList<Float>();
                for (String key : rowMaps[1].keySet())
                {
                    if (keysInAllFiles.contains(key))
                        continue;
                    List<TabLoader.ColumnDescriptor>[] colArray = new List[] { columnsAllFiles[1] };
                    Map[] colMap = new Map[] { rowMaps[1].get(key) };

                    String line = createFileLine(key, colArray,
                            colMap);
                    if (outUnique2File != null)
                    {
                        unique2OutWriter.println(line);
                        unique2OutWriter.flush();
                    }

                    if (plotColumnName != null && colMap[0].get(plotColumnName) != null)
                    {
                        try
                        {
                            plotColumnUnique2Values.add(columnValueAsFloat(colMap[0].get(plotColumnNameFile2)));
                        }
                        catch (ClassCastException e)
                        {}
                    }

                }
                if (plotColumnName != null & !plotColumnUnique2Values.isEmpty())
                {
                    PanelWithHistogram pwh = new PanelWithHistogram(plotColumnUnique2Values, "Values unique to 2");
                    pwh.displayInTab();
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

                List<Float> trackValues1 = new ArrayList<Float>();
                List<Float> trackValues2 = new ArrayList<Float>();


                Map<String,Map> rowMaps1 = rowMaps[0];
                Map<String,Map> rowMaps2 = rowMaps[1];

                for (String key : rowMaps2.keySet())
                {
                    if (!rowMaps1.containsKey(key))
                    {
                        Object o2 = rowMaps2.get(key).get(plotColumnNameFile2);
                        if (valuesToTrack != null && valuesToTrack.contains(key))
                        {
                            System.err.println(key + "\tNA\t" + o2);
                        }
                    }
                }
                for (String key : rowMaps1.keySet())
                {
//System.err.println("Key: " + key);
                    Object o1 = rowMaps1.get(key).get(plotColumnName);

                    if (rowMaps2.containsKey(key))
                    {
//System.err.println("\t" + o1 + rowMaps2.get(key).get(plotColumnName));
                        Object o2 = rowMaps2.get(key).get(plotColumnNameFile2);
                        if (o1 == null || o2 == null)
                            continue;
//if (key.equals("IPI00115660")) System.err.println("@@@" + o1 + ", " + o2);
                        try
                        {
                            float value1 = columnValueAsFloat(o1);
                            float value2 = columnValueAsFloat(o2);
// System.err.println("Unplottable! " + value1 + ", " + value2);

                            float displayValue1 = value1;
                            float displayValue2 =  value2;
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
//System.err.println("***" + displayValue1 + ", " + displayValue2);                                
                                values1.add(displayValue1);
                                values2.add(displayValue2);
                                commonKeys.add(key);

                                if (valuesToTrack != null && valuesToTrack.contains(key))
                                {
//         System.err.println(key + "\t" + displayValue1 + "\t" + displayValue2);

                                    trackValues1.add(displayValue1);
                                    trackValues2.add(displayValue2);
                                }
                            }
                        }
                        catch (ClassCastException e)
                        {        ApplicationContext.infoMessage("Crap!  Can't process value " +
                                rowMaps1.get(key).get(plotColumnName) + " or " + rowMaps2.get(key).get(plotColumnName));
                        }

                    }
                    else
                    {
                                if (valuesToTrack != null && valuesToTrack.contains(key))
                                {
//         System.err.println(key + "\t" + o1 + "\tNA");
                                }
                    }
                }
                ApplicationContext.infoMessage("Rows in common and plottable: " + values1.size());
                ApplicationContext.infoMessage("Correlation coefficient: " + BasicStatistics.correlationCoefficient(values1, values2));
                PanelWithScatterPlot pwsp = new PanelWithScatterPlot(values1, values2, plotColumnName);
                pwsp.setAxisLabels("File 1","File 2");
                pwsp.displayInTab();

                List<Float> differences = new ArrayList<Float>();
                for (int i=0; i< values1.size(); i++)
                    differences.add(values2.get(i)-values1.get(i));
                new PanelWithHistogram(differences, "Differences").displayInTab();

                if (valuesToTrack != null && trackValues1.size() > 0)
                {
                    PanelWithScatterPlot pwsp2 = new PanelWithScatterPlot(trackValues1, trackValues2, plotColumnName + "_track");
                    pwsp2.setAxisLabels("File 1","File 2");
                    pwsp2.displayInTab();
                }

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

    protected float columnValueAsFloat(Object columnValueObject)
            throws ClassCastException
    {
        float result = 0;

        try
        {
            result = ((Double) columnValueObject).floatValue();
        }
        catch(ClassCastException e)
        {
            try
            {
                result = ((Integer) columnValueObject).floatValue();
            }
            catch(ClassCastException e2)
            {
                try
                {
                    result = ((Float) columnValueObject);
                }
                catch(ClassCastException e3)
                {
                    result = Float.parseFloat ((String) columnValueObject);
                }
            }
        }
        return result;
    }


    /**
     * Read each line as a String, stopping at first whitespace
     * @param file
     * @return
     * @throws ArgumentValidationException
     */
    protected List<String> readOneStringPerLine(File file)
            throws ArgumentValidationException
    {
        List<String> result =  new ArrayList<String>();
        try
        {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = null;
            while ((line = br.readLine()) != null)
            {
                String protein = StringUtils.strip(line);
                //if there's more than one column in the file, take first
                protein = protein.replaceFirst("\\s.*", "");

                result.add(protein);
            }

        }
        catch (Exception e)
        {
            throw new ArgumentValidationException("Failed to retrieve list from file " +
                    file + ".  Please make sure file contains only a list of identifiers, one per line");
        }
        return result;
    }

}
