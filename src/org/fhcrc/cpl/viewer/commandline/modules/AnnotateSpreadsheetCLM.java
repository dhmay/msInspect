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
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.chem.ChemCalcs;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * This is an unwieldy jumble of code to help merge together two spreadsheets based on shared values in the
 * 'mergecolumn' column.  In several different ways.
 */
public class AnnotateSpreadsheetCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(AnnotateSpreadsheetCLM.class);

    protected String mergeColumnName = null;
    protected String annotationMergeColumnName = null;
    protected File inFile;
    protected File annotationFile;
    protected File outFile;
    protected List<String> annotationColumnNames = null;
    protected List<String> annotationColumnNamesNotAlreadyPresent = null;

    protected String multiAnnotationSeparator = ";";

    protected boolean shouldCollapseByMax = false;
    protected boolean shouldMergeValueForMissing = false;


    public AnnotateSpreadsheetCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "annotatespreadsheet";
        mShortDescription = "annotate a .tsv spreadsheet with annotations from another .tsv spreadsheet";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        createUnnamedFileArgumentDefinition(true,"input spreadsheet"),
                        new StringArgumentDefinition("mergecolumn", true, "column to merge on, in the main file (and in the annotation file, unless annotationmergecolumn specified"),
                        new StringArgumentDefinition("annotationmergecolumn", false, "column name in the annotation file to merge on"),
                        new StringListArgumentDefinition("annotationcolumns", true, "column(s) to annotate with, separated by commas"),
                        new FileToWriteArgumentDefinition("out", true, "output file"),
                        new FileToReadArgumentDefinition("annotationfile", true, "Annotation file"),
                        new StringArgumentDefinition("multiannotseparator", false, "symbol that separates multiple merge keys in the annotation list",multiAnnotationSeparator),
                        new BooleanArgumentDefinition("collapsebymax", false, "If multiple keys with multiple annotations, use the maximum? (default: semicolon-separated", shouldCollapseByMax),
                        new BooleanArgumentDefinition("usemergevalueformissing", false,
                                "Iff this is true, then if the annotation value is missing, it will be populated with the merge column value",
                                shouldMergeValueForMissing),
                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        inFile = getUnnamedFileArgumentValue();
        mergeColumnName = getStringArgumentValue("mergecolumn");
        annotationMergeColumnName = getStringArgumentValue("annotationmergecolumn");
        if (annotationMergeColumnName == null)
            annotationMergeColumnName = mergeColumnName;
ApplicationContext.infoMessage("Main file merge column: " + mergeColumnName);
ApplicationContext.infoMessage("Annotation file merge column: " + annotationMergeColumnName);

        annotationColumnNames = (List<String>) getArgumentValue("annotationcolumns");
        StringBuffer buf = new StringBuffer("Annotation columns: ");
        for (String col : annotationColumnNames)
            buf.append(" " + col);
        ApplicationContext.infoMessage(buf.toString());
        annotationFile = getFileArgumentValue("annotationfile");
        outFile = getFileArgumentValue("out");
        multiAnnotationSeparator = getStringArgumentValue("multiannotseparator");
        shouldCollapseByMax = getBooleanArgumentValue("collapsebymax");
        shouldMergeValueForMissing = getBooleanArgumentValue("usemergevalueformissing");

ApplicationContext.infoMessage("Multi-annotation separator: " + multiAnnotationSeparator);        
    }

    protected Map<String,Map<String,Object>> loadAnnotations(File file)
            throws IOException
    {

        ApplicationContext.infoMessage("Loading annotations from file " + file.getAbsolutePath());
        TabLoader loader = new TabLoader(new FileReader(file),true);

ApplicationContext.infoMessage("Annotation file columns:");
for (TabLoader.ColumnDescriptor col : loader.getColumns()) ApplicationContext.infoMessage("\t" + col.name);
        //a separate map for each annotation column
        Map<String,Map<String,Object>> result = new HashMap<String,Map<String,Object>>();
//System.err.println("Loading annotations...");
//Set<String> allElements = new HashSet<String>();
        for (Map row : (Map[]) loader.load())
        {
//try{
//allElements.addAll(ChemCalcs.emp2atomCount((String)row.get("formula")).keySet());
//}catch(Exception e){}

//System.err.println("ROW");
//for (Object key1 : row.keySet()) System.err.println("\t" + key1.toString());
            Object key = row.get(annotationMergeColumnName);
            if (key != null)
            {
//System.err.println("CHECK");
                String keyString = key.toString();
                List<String> allKeys = new ArrayList<String>();
                if (multiAnnotationSeparator != null)
                {
                    String[] chunks = keyString.split(multiAnnotationSeparator);
                    for (String chunk : chunks)
                        allKeys.add(chunk);
                }
                else
                    allKeys.add(keyString);

                for (String akey : allKeys)
                {
                    for (String annotationColumnName : annotationColumnNames)
                    {
//System.err.println("Checking " + annotationColumnName + "...");
                        if (row.get(annotationColumnName) == null)
                            continue;
                        Map<String, Object> mapThisAnnotColumn = result.get(annotationColumnName);
                        if (mapThisAnnotColumn == null)
                        {
                            mapThisAnnotColumn = new HashMap<String, Object>();
                            result.put(annotationColumnName, mapThisAnnotColumn);
                        }
                        mapThisAnnotColumn.put(akey,row.get(annotationColumnName));
//System.err.println("\t found " +row.get(annotationColumnName) );                        
                    }
                }
            }
        }
        ApplicationContext.infoMessage("Loaded annotation values for " + result.size() + " columns:");
        for (String annot : result.keySet())
        {
            ApplicationContext.infoMessage("\t" + annot + ": " + result.get(annot).size());
        }
//for (String element : allElements) System.err.println(element);
        return result;
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

            Map<String,Map<String,Object>> annotationsMap = loadAnnotations(annotationFile);

            for (String annotColName : annotationColumnNames)
                if (!annotationsMap.containsKey(annotColName))
                    throw new CommandLineModuleExecutionException("Requested annotation columns " + annotColName +
                            " not present in annotation file " + annotationFile.getName());
//            ApplicationContext.infoMessage("Loaded annotations for " + annotationsMap.size() + " columns:");
//            for (String annotColumn : annotationsMap.keySet())
//                ApplicationContext.infoMessage("\t" + annotColumn);


            StringBuffer headerLine = new StringBuffer();

            System.err.println("Loading file " + inFile.getName());
            //forcing reading the first file line as a header line, since we're dependent on column names
            TabLoader loader = new TabLoader(new FileReader(inFile), true);
            TabLoader.ColumnDescriptor[] columnDescriptors = loader.getColumns();
ApplicationContext.infoMessage("Main file columns:");
for (TabLoader.ColumnDescriptor col : columnDescriptors) ApplicationContext.infoMessage("\t" + col.name);
            ApplicationContext.infoMessage("Checking annotation columns...");
            annotationColumnNamesNotAlreadyPresent = new ArrayList<String>();

            for (String annotationColumnName : annotationColumnNames)
            {
                boolean alreadyPresent = false;
                for (TabLoader.ColumnDescriptor column : columnDescriptors)
                {
                    if (annotationColumnName.equals(column.name))
                    {
                        ApplicationContext.infoMessage("\tAnnotation column " + annotationColumnName + " already present");
                        alreadyPresent = true;
                        break;
                    }
                }
                if (!alreadyPresent)
                {
                    ApplicationContext.infoMessage("\tAnnotation column " + annotationColumnName + " NOT already present");
                    annotationColumnNamesNotAlreadyPresent.add(annotationColumnName);
                }
            }

            boolean first = true;
            for (TabLoader.ColumnDescriptor column : columnDescriptors)
            {
                if (!first)
                    headerLine.append("\t");
                headerLine.append(column.name);
                first = false;
            }
            for (String annotationColumnName : annotationColumnNamesNotAlreadyPresent)
                headerLine.append("\t" + annotationColumnName);

            if (outFile != null)
            {
                outPW.println(headerLine.toString());
                outPW.flush();
            }

            int annotatedRowCount = 0;
            int fullRowCount = 0;

            for (Map rowMap : (Map[]) loader.load())
            {
                fullRowCount++;
                StringBuffer fileLineBuf = new StringBuffer();
                boolean firstCol = true;

                Map<String, Object> myRowMap = new HashMap<String, Object>();
                for (Object key : rowMap.keySet())
                    myRowMap.put(key.toString(), rowMap.get(key));


                Object mergeColumnValue = myRowMap.get(mergeColumnName);
                List<String> allKeys = new ArrayList<String>();
                boolean isThisRowAnnotated = false;

                if (mergeColumnValue != null)
                {
                    if (multiAnnotationSeparator != null)
                    {
                        String[] chunks = mergeColumnValue.toString().split(multiAnnotationSeparator);
                        for (String chunk : chunks)
                        {
                            if (chunk != null)
                                allKeys.add(chunk);

                        }
                    }
                    else
                        allKeys.add(mergeColumnValue.toString());
                }

                for (String annotationColumnName : annotationColumnNames)
                {
                    if (!allKeys.isEmpty())
                    {
                        List<String> annotationValuesAllKeys = new ArrayList<String>();
                        Map<String, Object> annotationMap = annotationsMap.get(annotationColumnName);
                        if (annotationMap == null)
                            throw new CommandLineModuleExecutionException("Missing annotation map for column " + annotationColumnName);
                        for (String key : allKeys)
                        {
                            if (annotationMap.containsKey(key))
                            {
                                annotationValuesAllKeys.add(annotationMap.get(key).toString());
                            }
                        }
                        if (!annotationValuesAllKeys.isEmpty())
                        {
                            isThisRowAnnotated=true;

                            myRowMap.put(annotationColumnName,collapseAnnotations(annotationValuesAllKeys));
                        }

                    }
                }

                for (TabLoader.ColumnDescriptor column : columnDescriptors)
                {
                    if (!firstCol)
                        fileLineBuf.append("\t");
                    if (myRowMap != null && myRowMap.get(column.name) != null)
                        fileLineBuf.append(myRowMap.get(column.name));
                    firstCol = false;
                }

                for (String annotationColumnName : annotationColumnNamesNotAlreadyPresent)
                {
                    fileLineBuf.append("\t");
                    if (myRowMap != null)
                    {
                        if (myRowMap.get(annotationColumnName) != null)
                            fileLineBuf.append(myRowMap.get(annotationColumnName));
                        else if (shouldMergeValueForMissing)
                            fileLineBuf.append(mergeColumnValue.toString());

                    }
                }

                if (isThisRowAnnotated)
                    annotatedRowCount++;
                outPW.println(fileLineBuf.toString());
                outPW.flush();
            }
            ApplicationContext.infoMessage("Annotated " + annotatedRowCount + " out of " + fullRowCount + " rows");
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("FAIL!!",e);
        }
        outPW.close();
    }

    protected String collapseAnnotations(List<String> annotationValues)
    {
        String result = null;
        if (shouldCollapseByMax)
        {
            //todo: generalize to float/double
            int max = Integer.MIN_VALUE;
            for (String value : annotationValues)
            {
                int current = Integer.parseInt(value);
                if (current > max)
                    max = current;
            }
            result = "" + max;
        }
        else
        {
            StringBuffer resultBuf = new StringBuffer();
            for (int i=0; i<annotationValues.size(); i++)
            {
                if (i>0)
                    resultBuf.append(";");
                resultBuf.append(annotationValues.get(i));

            }
            result = resultBuf.toString();
        }
        return result;
    }

}
