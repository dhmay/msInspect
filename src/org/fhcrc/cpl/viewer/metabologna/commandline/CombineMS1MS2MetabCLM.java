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
package org.fhcrc.cpl.viewer.metabologna.commandline;

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.BaseViewerCommandLineModuleImpl;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * This is an unwieldy jumble of code to help merge together two spreadsheets based on shared values in the
 * 'mergecolumn' column.  In several different ways.
 */
public class CombineMS1MS2MetabCLM extends BaseViewerCommandLineModuleImpl
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(CombineMS1MS2MetabCLM.class);

    protected File ms1File;
    protected File ms2File;
    protected File outFile;


    protected boolean shouldCollapseByMax = false;

    public CombineMS1MS2MetabCLM()
    {
        init();
    }

    protected void init()
    {
        mCommandName = "metabms1ms2";
        mShortDescription = "Combine MS1 and MS2 metabolite ID info";
        mHelpMessage = mShortDescription;
        CommandLineArgumentDefinition[] argDefs =
                {
                        new FileToReadArgumentDefinition("ms1",true,"input MS1 spreadsheet"),
                        new FileToReadArgumentDefinition("ms2",true,"input MS1 spreadsheet"),
                        new FileToWriteArgumentDefinition("out",true,"output file"),

                };
        addArgumentDefinitions(argDefs);
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        ms1File = getFileArgumentValue("ms1");
        ms2File = getFileArgumentValue("ms2");


        outFile = getFileArgumentValue("out");
    }

    public void execute() throws CommandLineModuleExecutionException
    {
        Map<Integer,List<Map<String,Object>>> ms1RowMap = loadRowsByScan(ms1File,"MS2Scans");
        Map<Integer,List<Map<String,Object>>> ms2RowMap = loadRowsByScan(ms2File,"scan");

        PrintWriter outPW = null;
        try
        {
            outPW = new PrintWriter(outFile);
        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to load file " + outFile.getAbsolutePath());
        }

        int scansInCommon = 0;
        outPW.println("scan\tid");
        List<Integer> agreedScans = new ArrayList<Integer>();
        List<Float> ms1Masses = new ArrayList<Float>();
        List<Float> massDiffs = new ArrayList<Float>();

        for (int scan : ms1RowMap.keySet())
        {
            if (ms2RowMap.containsKey(scan))
            {
                List<Map<String,Object>> ms1Rows = ms1RowMap.get(scan);
                List<Map<String,Object>> ms2Rows = ms2RowMap.get(scan);

                scansInCommon++;

                float ms1Mass = 0;
                float ms2Mass = 0;


System.err.println("*****" + scan);
                Set<String> ms1IDs = new HashSet<String>();
                for (Map<String,Object> rows : ms1Rows)
                {
System.err.println("\t" + rows.get("Compound"));
                    ms1IDs.add(stripOuterQuotes((String) rows.get("Compound")));
                    ms1Mass = Float.parseFloat(rows.get("FeatureMass").toString());
                }
                Set<String> ms2IDs = new HashSet<String>();
                for (Map<String,Object> rows : ms2Rows)
                {
  System.err.println("\t\t" + rows.get("Metlin ID Name"));

                    ms2IDs.add(stripOuterQuotes((String) rows.get("Metlin ID Name")));

                    ms2Mass = Float.parseFloat(rows.get("Metlin Mass").toString());

                }

                ms1Masses.add(ms1Mass);
                massDiffs.add(ms1Mass - ms2Mass);

                Set<String> intersectionIds = new HashSet<String>(ms1IDs);
                intersectionIds.retainAll(ms2IDs);
                if (!intersectionIds.isEmpty())
                {
                    agreedScans.add(scan);
                    for (String id : intersectionIds)
                        outPW.println(scan + "\t" + id);
                    outPW.flush();
                }
            }

        }
        outPW.close();
        System.err.println(scansInCommon + " scans in common.  Of those, " + agreedScans.size() + " agreed");
        new PanelWithScatterPlot(ms1Masses, massDiffs, "MS1 vs Mass Diff").displayInTab();
    }

    protected String stripOuterQuotes(String input)
    {
        if (input.charAt(0) == '"' && input.charAt(input.length()-1) == '"')
            return input.substring(1, input.length()-1);
        return input;
    }



    protected Map<Integer,List<Map<String,Object>>> loadRowsByScan(File file, String scanColName)
            throws CommandLineModuleExecutionException
    {
        ApplicationContext.infoMessage("Loading info from file " + file.getAbsolutePath());
        try
        {
            TabLoader loader = new TabLoader(new FileReader(file),true);
//for (TabLoader.ColumnDescriptor col : loader.getColumns()) System.err.println("\t" + col.name);
            Map<Integer,List<Map<String,Object>>> result = new HashMap<Integer,List<Map<String,Object>>>();
            for (Map row : (Map[]) loader.load())
            {
                Object scanObj = row.get(scanColName);
                if (scanObj == null)
                    continue;
                List<Integer> scans = new ArrayList<Integer>();
                if (scanObj.getClass().isAssignableFrom(String.class))
                    scans.addAll(MS2ExtraInfoDef.parseIntListString((String) scanObj));
                else
                    scans.add((Integer) scanObj);
                for (int scan : scans)
                {
System.err.println(scan);                    
                    List<Map<String,Object>> rowsThisScan = result.get(scan);
                    if (rowsThisScan == null)
                    {
                        rowsThisScan = new ArrayList<Map<String, Object>>();
                        result.put(scan, rowsThisScan);
                    }
                    rowsThisScan.add(row);
                }
            }
            ApplicationContext.setMessage("Loaded " + result.size() + " scans");
            return result;

        }
        catch (IOException e)
        {
            throw new CommandLineModuleExecutionException("Failed to load file " + file.getAbsolutePath(),e);
        }
    }

}
