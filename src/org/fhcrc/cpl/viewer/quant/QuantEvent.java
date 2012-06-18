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

package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithPeakChart;
import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.Rounder;
import org.fhcrc.cpl.viewer.quant.commandline.FlagQuantEventsCLM;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.table.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.event.*;

/**
 * Holds all the information related to a quantitative event for display.
 *
 * And does a lot of other stuff, like the inner classes for showing tables of QuantEvents.
 * TODO: break some of this stuff up
 */
public class QuantEvent
{
    protected static Logger _log = Logger.getLogger(QuantEvent.class);


    public static final float MISSING_RATIO = -1;

    //curation status values
    public static final int CURATION_STATUS_UNKNOWN = 0;
    public static final int CURATION_STATUS_GOOD = 1;
    public static final int CURATION_STATUS_BAD = 2;
    public static final int CURATION_STATUS_RATIO_ONEPEAK = 3;


    //calculated ratio using only the theoretically highest peak
    protected float ratioOnePeak;

    protected String protein;
    protected String peptide;
    protected String fraction;
    protected int charge;
    protected String modificationState;
    protected int scan;
    protected float mz;
    protected File spectrumFile;
    protected File scansFile;
    protected File file3D;
    protected File intensitySumFile;

    protected float ratio;
    protected float peptideProphet;

    protected float lightMz;
    protected float heavyMz;
    protected float lightIntensity;
    protected float heavyIntensity;
    protected int firstLightQuantScan;
    protected int lastLightQuantScan;
    protected int firstHeavyQuantScan;
    protected int lastHeavyQuantScan;
    protected List<QuantEvent> otherEvents;

    protected Feature sourceFeature;

    protected String comment;


    protected int quantCurationStatus = CURATION_STATUS_UNKNOWN;
    protected int idCurationStatus = CURATION_STATUS_UNKNOWN;

    protected List<ActionListener> quantCurationStatusListeners = new ArrayList<ActionListener>();
    protected List<ActionListener> algAssessmentStatusListeners = new ArrayList<ActionListener>();


    protected QuantEventAssessor.QuantEventAssessment algorithmicAssessment = null;

    /**
     * Create a deep copy of another QuantEvent, preserving everything except curation status
     * @param eventToCopy
     */
    public QuantEvent(QuantEvent eventToCopy)
    {
        this.charge = eventToCopy.getCharge();
        this.mz = eventToCopy.getMz();
        this.peptide = eventToCopy.getPeptide();
        this.peptideProphet = eventToCopy.getPeptideProphet();
        this.protein = eventToCopy.getProtein();
        this.lightMz = eventToCopy.lightMz;
        this.heavyMz = eventToCopy.heavyMz;
        this.lightIntensity = eventToCopy.lightIntensity;
        this.heavyIntensity = eventToCopy.heavyIntensity;
        this.firstLightQuantScan = eventToCopy.firstLightQuantScan;
        this.lastLightQuantScan = eventToCopy.lastLightQuantScan;
        this.firstHeavyQuantScan = eventToCopy.firstHeavyQuantScan;
        this.lastHeavyQuantScan = eventToCopy.lastHeavyQuantScan;
        this.ratio = eventToCopy.ratio;
        this.spectrumFile = eventToCopy.spectrumFile;
        this.scansFile = eventToCopy.scansFile;
        this.file3D = eventToCopy.file3D;
        this.intensitySumFile = eventToCopy.intensitySumFile;
        this.fraction = eventToCopy.fraction;
        this.otherEvents = eventToCopy.otherEvents;
        this.modificationState = eventToCopy.modificationState;
        this.ratioOnePeak = eventToCopy.ratioOnePeak;
        this.algorithmicAssessment = eventToCopy.algorithmicAssessment;
    }

    /**
     * Create a QuantEvent from a Feature and a fraction name.
     * HACK: use description as comment
     * @param feature
     *
     * @param fraction
     */
    public QuantEvent(Feature feature, String fraction)
    {
        sourceFeature = feature;
        String protein = MS2ExtraInfoDef.getFirstProtein(feature);
        String peptide = MS2ExtraInfoDef.getFirstPeptide(feature);
        float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(feature);
        float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(feature);
        float heavyIntensity = (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);
        float lightNeutralMass = (float) IsotopicLabelExtraInfoDef.getLightMass(feature) -
                Spectrum.HYDROGEN_ION_MASS;
        float heavyNeutralMass = (float) IsotopicLabelExtraInfoDef.getHeavyMass(feature) -
                Spectrum.HYDROGEN_ION_MASS;
        int charge = feature.getCharge();

        float lightMz = lightNeutralMass / charge + Spectrum.HYDROGEN_ION_MASS;
        float heavyMz = heavyNeutralMass / charge + Spectrum.HYDROGEN_ION_MASS;

        int firstLightQuantScan = IsotopicLabelExtraInfoDef.getLightFirstScan(feature);
        int lastLightQuantScan = IsotopicLabelExtraInfoDef.getLightLastScan(feature);
        int firstHeavyQuantScan = IsotopicLabelExtraInfoDef.getHeavyFirstScan(feature);
        int lastHeavyQuantScan = IsotopicLabelExtraInfoDef.getHeavyLastScan(feature);

        if (feature.comprised != null && feature.comprised.length > 0)
        {
            otherEvents = new ArrayList<QuantEvent>();
            for (Spectrum.Peak peak : feature.comprised)
            {
                otherEvents.add(new QuantEvent((Feature) peak, fraction));
            }
        }

        String eventComment = feature.getDescription();
        //This is a HACK.  The quant-event flagger uses a dummy search score to store the String describing why
        //an event was flagged
        if (eventComment == null)
            eventComment = MS2ExtraInfoDef.getSearchScore(feature, FlagQuantEventsCLM.REASON_DUMMY_SEARCH_SCORE_NAME);

        init(protein, peptide, fraction, charge,
                MS2ExtraInfoDef.convertModifiedAminoAcidsMapToString(MS2ExtraInfoDef.getModifiedAminoAcidsMap(feature)),
                feature.getScan(), feature.getMz(),
                null, null, null, null,
                ratio, lightMz, heavyMz, lightIntensity, heavyIntensity,
                firstLightQuantScan, lastLightQuantScan, firstHeavyQuantScan, lastHeavyQuantScan,
                otherEvents,
                (float) MS2ExtraInfoDef.getPeptideProphet(feature),
                CURATION_STATUS_UNKNOWN, CURATION_STATUS_UNKNOWN, eventComment);
    }

    /**
     * Create a QuantEvent from a Feature and a fraction name, pointing at existing chart files
     * @param feature
     * @param fraction
     * @param spectrumFile
     * @param scansFile
     * @param file3D
     * @param intensitySumFile
     */
    public QuantEvent(Feature feature, String fraction, File spectrumFile, File scansFile,
                          File file3D, File intensitySumFile)
    {
        this(feature, fraction);
        this.spectrumFile = spectrumFile;
        this.scansFile = scansFile;
        this.file3D = file3D;
        this.intensitySumFile = intensitySumFile;
    }

    /**
     * This is a hack.  We're modeling 'other' events subsumed by a QuantEvent in the tsv file
     * as lists of scan numbers and mz values.  So when we read that in, need to turn those into a list
     * of events identical to this one except scan and mz.
     * @param protein
     * @param peptide
     * @param fraction
     * @param charge
     * @param scan
     * @param mz
     * @param spectrumFile
     * @param scansFile
     * @param file3D
     * @param intensitySumFile
     * @param ratio
     * @param lightMz
     * @param heavyMz
     * @param lightIntensity
     * @param heavyIntensity
     * @param firstLightQuantScan
     * @param lastLightQuantScan
     * @param firstHeavyQuantScan
     * @param lastHeavyQuantScan
     * @param otherEventScans
     * @param otherEventMzs
     * @param peptideProphet
     * @param quantCurationStatus
     * @param idCurationStatus
     * @param comment
     */
    public QuantEvent(String protein, String peptide, String fraction, int charge, String modificationState, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D,
                          File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<Integer> otherEventScans, List<Float> otherEventMzs,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {

        init(protein, peptide, fraction, charge, modificationState, scan, mz, spectrumFile, scansFile, file3D, intensitySumFile,
                ratio,
                lightMz, heavyMz, lightIntensity, heavyIntensity, firstLightQuantScan, lastLightQuantScan,
                firstHeavyQuantScan, lastHeavyQuantScan, null,
                peptideProphet, quantCurationStatus, idCurationStatus, comment);
        otherEvents = new ArrayList<QuantEvent>();
        for (int i=0; i<otherEventScans.size(); i++)
        {
            QuantEvent otherEvent = new QuantEvent(this);
            otherEvent.otherEvents = null;
            otherEvent.scan = otherEventScans.get(i);
            otherEvent.mz = otherEventMzs.get(i);
            otherEvents.add(otherEvent);
        }

    }

    /**
     * QuantEvent constructor that takes everything stored in a QuantEvent
     * @param protein
     * @param peptide
     * @param fraction
     * @param charge
     * @param modificationState
     * @param scan
     * @param mz
     * @param spectrumFile
     * @param scansFile
     * @param file3D
     * @param intensitySumFile
     * @param ratio
     * @param lightMz
     * @param heavyMz
     * @param lightIntensity
     * @param heavyIntensity
     * @param firstLightQuantScan
     * @param lastLightQuantScan
     * @param firstHeavyQuantScan
     * @param lastHeavyQuantScan
     * @param otherEvents
     * @param peptideProphet
     * @param quantCurationStatus
     * @param idCurationStatus
     * @param comment
     */
    public QuantEvent(String protein, String peptide, String fraction, int charge, String modificationState,
                          int scan, float mz, File spectrumFile, File scansFile, File file3D,
                          File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<QuantEvent> otherEvents,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {
        init(protein, peptide, fraction, charge, modificationState, scan, mz, spectrumFile, scansFile, file3D, intensitySumFile,
                ratio,
                lightMz, heavyMz, lightIntensity, heavyIntensity, firstLightQuantScan, lastLightQuantScan,
                firstHeavyQuantScan, lastHeavyQuantScan, otherEvents,
                peptideProphet, quantCurationStatus, idCurationStatus, comment);
    }


    /**
     * Initialize all values
     * @param protein
     * @param peptide
     * @param fraction
     * @param charge
     * @param modificationState
     * @param scan
     * @param mz
     * @param spectrumFile
     * @param scansFile
     * @param file3D
     * @param intensitySumFile
     * @param ratio
     * @param lightMz
     * @param heavyMz
     * @param lightIntensity
     * @param heavyIntensity
     * @param firstLightQuantScan
     * @param lastLightQuantScan
     * @param firstHeavyQuantScan
     * @param lastHeavyQuantScan
     * @param otherEvents
     * @param peptideProphet
     * @param quantCurationStatus
     * @param idCurationStatus
     * @param comment
     */
    protected void init(String protein, String peptide, String fraction, int charge, String modificationState, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D, File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<QuantEvent> otherEvents,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {
        this.protein = protein;
        this.peptide = peptide;
        this.fraction = fraction;
        this.charge = charge;
        this.modificationState = modificationState;
        this.scan=scan;
        this.mz=mz;
        this.spectrumFile = spectrumFile;
        this.scansFile = scansFile;
        this.file3D = file3D;
        this.intensitySumFile = intensitySumFile;

        this.ratio = ratio;
        this.lightMz = lightMz;
        this.heavyMz = heavyMz;
        this.lightIntensity = lightIntensity;
        this.heavyIntensity = heavyIntensity;
        this.firstLightQuantScan = firstLightQuantScan;
        this.lastLightQuantScan = lastLightQuantScan;
        this.firstHeavyQuantScan = firstHeavyQuantScan;
        this.lastHeavyQuantScan = lastHeavyQuantScan;
        this.otherEvents = otherEvents;
        this.peptideProphet = peptideProphet;

        this.quantCurationStatus = quantCurationStatus;
        this.idCurationStatus = idCurationStatus;
        this.comment = comment;
    }

    public String toString()
    {
        return "QuantEvent: peptide=" + peptide + ", protein=" + protein + ", fraction=" + fraction +
                ", charge=" + charge + ", scan=" + scan + ", scanrange=" + this.firstLightQuantScan + "-" +
                this.lastLightQuantScan;
    }

    /**
     * Check to see if this is the same quantitative event as another event.  This is NOT the same as
     * equals() -- this checks the peptide, fraction, and charge fields for equivalence and checks
     * whether the scan number is the same or whether any otherEvents scans are the same.  If all
     * that passes, true, otherwise false.
     *
     * This is used for checking whether we have already seen this event before in a list of events.
     * @param compareEvent
     * @return
     */
    public boolean isSameEvent(QuantEvent compareEvent)
    {
        if (!peptide.equals(compareEvent.peptide) || !fraction.equals(compareEvent.fraction) ||
                charge != compareEvent.charge)
            return false;
        List<Integer> otherEventScans = new ArrayList<Integer>();
        otherEventScans.add(compareEvent.scan);
        if (compareEvent.otherEvents != null)
            for (QuantEvent compareOtherEvent : compareEvent.otherEvents)
                otherEventScans.add(compareOtherEvent.scan);
        if (otherEventScans.contains(scan))
            return true;
        if (otherEvents != null)
            for (QuantEvent otherEvent : otherEvents)
                if (otherEventScans.contains(otherEvent.getScan()))
                    return true;
        return false;

    }


    public String createOutputRowHtml(String outChartsRelativeDirPath, boolean showProteinColumn,
                                  boolean show3DColumn)
    {
        return createOutputRow(outChartsRelativeDirPath, true, showProteinColumn, show3DColumn);
    }

    public String createOutputRowTsv(boolean showProteinColumn,
                                  boolean show3DColumn)
    {
        return createOutputRow(null, false, showProteinColumn, show3DColumn);
    }

    /**
     * Get a Map that can be used for display/writing of all fields
     * @return
     */
    public Map<String, String> getNameValueMapNoCharts()
    {
        Map<String, String> result = new HashMap<String, String>();

        if (protein != null)
            result.put("Protein", protein);
        result.put("Peptide",    peptide);
        result.put("Fraction", fraction);
        result.put("Charge",  "" + charge);
        result.put("Modifications",  modificationState);         
        result.put("Scan",  "" + scan);
        result.put("Mz",  "" + mz);
        result.put("Ratio",  "" + ratio);
        result.put("LightMz", "" + lightMz);
        result.put("HeavyMz", "" + heavyMz);
        result.put("LightInt", "" + lightIntensity);
        result.put("HeavyInt", "" + heavyIntensity);
        result.put("LightFirstScan", "" + firstLightQuantScan);
        result.put("LightLastScan", "" + lastLightQuantScan);
        result.put("HeavyFirstScan",  "" + firstHeavyQuantScan);
        result.put("HeavyLastScan", "" + lastHeavyQuantScan);
        result.put("OtherEventScans", convertOtherEventScansToString());
        result.put("OtherEventMZs", convertOtherEventMZsToString());
        result.put("PProphet", "" + peptideProphet);

        return result;
    }

    protected String convertOtherEventScansToString()
    {
        List<String> allScansAsStrings = new ArrayList<String>();
        if (otherEvents != null)
        {
            for (QuantEvent otherEvent : otherEvents)
                allScansAsStrings.add("" + otherEvent.getScan());
        }
        return MS2ExtraInfoDef.convertStringListToString(allScansAsStrings);
    }

    protected String convertOtherEventMZsToString()
    {
        List<String> allMzsAsStrings = new ArrayList<String>();
        if (otherEvents != null)
        {
            for (QuantEvent otherEvent : otherEvents)
                allMzsAsStrings.add("" + otherEvent.getMz());
        }
        return MS2ExtraInfoDef.convertStringListToString(allMzsAsStrings);
    }

    /**
     * Create a line of text representing this QuantEvent, either for an HTML page or for a TSV file
     * TODO: get rid of HTML generation entirely?
     * @param outChartsRelativeDirPath
     * @param isHtml
     * @param showProteinColumn
     * @param show3DColumn
     * @return
     */
    protected String createOutputRow(String outChartsRelativeDirPath, boolean isHtml, boolean showProteinColumn,
                                  boolean show3DColumn)
    {
        String spectrumFileString = spectrumFile == null ? "" : spectrumFile.getAbsolutePath();
        String scansFileString = scansFile == null ? "" : scansFile.getAbsolutePath();
        String fileString3D =  null == file3D ? null : file3D.getAbsolutePath();
        String intensitySumFileString = null == intensitySumFile ? null : intensitySumFile.getAbsolutePath();


        if (isHtml)
        {
            spectrumFileString =  spectrumFile == null ? "" :HtmlGenerator.createLink(outChartsRelativeDirPath + spectrumFile.getName(),"Spectrum");
            scansFileString =  scansFile == null ? "" :HtmlGenerator.createLink(outChartsRelativeDirPath + scansFile.getName(), "Scans");
            if (show3DColumn)
                fileString3D = null == file3D ? null :HtmlGenerator.createLink(outChartsRelativeDirPath + file3D.getName(), "3D");
            intensitySumFileString =
                    null == intensitySumFile ? null :HtmlGenerator.createLink(outChartsRelativeDirPath + intensitySumFile.getName(), "intSumChart");

        }

        List<String> stringValuesForRow = new ArrayList<String>();
        if (showProteinColumn)
            stringValuesForRow.add(protein);
        stringValuesForRow.add(peptide);
        stringValuesForRow.add(fraction);
        stringValuesForRow.add( "" + charge);
        stringValuesForRow.add(modificationState);
        stringValuesForRow.add( "" + scan);
        stringValuesForRow.add( "" + mz);
        stringValuesForRow.add( "" + peptideProphet);
        stringValuesForRow.add(spectrumFileString);
        stringValuesForRow.add(scansFileString);
        stringValuesForRow.add(intensitySumFileString);           
        if (show3DColumn)
            stringValuesForRow.add(fileString3D);
        stringValuesForRow.add("" + ratio);
        stringValuesForRow.add("" + lightMz);
        stringValuesForRow.add("" + heavyMz);
        stringValuesForRow.add("" + lightIntensity);
        stringValuesForRow.add("" + heavyIntensity);
        stringValuesForRow.add("" + firstLightQuantScan);
        stringValuesForRow.add("" + lastLightQuantScan);
        stringValuesForRow.add( "" + firstHeavyQuantScan);
        stringValuesForRow.add("" + lastHeavyQuantScan);

        //dhmay adding 20090723
        stringValuesForRow.add("" + ratioOnePeak);

        stringValuesForRow.add("" + convertOtherEventScansToString());
        stringValuesForRow.add("" + convertOtherEventMZsToString());

        stringValuesForRow.add(convertCurationStatusToString(quantCurationStatus));
        stringValuesForRow.add(convertCurationStatusToString(idCurationStatus));

        stringValuesForRow.add(null == comment ? "" : comment);        

        String assessmentStatus = "";
        String assessmentDesc = "";
        if (algorithmicAssessment != null)
        {
            assessmentStatus = QuantEventAssessor.flagReasonCodes[algorithmicAssessment.getStatus()];
//System.err.println("**" + algorithmicAssessment.getStatus() + " --- " + assessmentStatus);          
            assessmentDesc = algorithmicAssessment.getExplanation();
        }
        stringValuesForRow.add("" + assessmentStatus);
        stringValuesForRow.add("" + assessmentDesc);

        //dhmay adding 20091123
        char[] bitmapChars = new char[QuantEventAssessor.flagReasonCodes.length];
        for (int i=0; i<algorithmicAssessment.getFlagBitmap().length; i++)
             bitmapChars[i] = algorithmicAssessment.getFlagBitmap()[i] ? '1' : '0';
        stringValuesForRow.add(new String(bitmapChars));


        String result = null;
        if (isHtml)
            result = HtmlGenerator.createTableRow(stringValuesForRow);
        else
        {
            StringBuffer resultBuf = new StringBuffer();
            boolean firstCol = true;
            for (String value : stringValuesForRow)
            {
                if (!firstCol)
                    resultBuf.append("\t");
                resultBuf.append(value);
                firstCol = false;
            }
            result = resultBuf.toString();
        }
        return result;
    }

    //all the columns in the properties table
    public static final String[] dataColumnNames = new String[]
            {
                    "Protein",
                    "Peptide",
                    "Fraction",
                    "Charge",
                    "Modifications",
                    "Scan",
                    "Mz",
                    "PProphet",
                    "Spectrum",
                    "Scans",
                    "intSumChart",
                    "3D",
                    "Ratio",
                    "LightMz",
                    "HeavyMz",
                    "LightInt",
                    "HeavyInt",
                    "LightFirstScan",
                    "LightLastScan",
                    "HeavyFirstScan",
                    "HeavyLastScan",
        //dhmay adding 20090723
                    "RatioOnePeak",
                    "OtherEventScans",
                    "OtherEventMZs",
                    "QuantCuration",
                    "IDCuration",
                    "Comment",
                    "AssessmentStatus",
                    "AssessmentDesc",
                    //dhmay adding 20091123
                    "AssessmentComponents"
            };

    protected boolean[] columnsAreFiles = new boolean[]
            {
                    false, false,false,false, false, false,true,true,true,
                    false,false,false,false,false,false,false,false,false,false, false, false
            };

    /**
     * Load quantitation events from a file and put them into a nested map by fraction and then by scan
     * @param eventFile
     * @return
     * @throws IOException
     */
    public static Map<String, Map<Integer, QuantEvent>> loadFractionScanQuantEventMap(File eventFile)
            throws IOException
    {
        List<QuantEvent> quantEvents = loadQuantEvents(eventFile);
        Map<String, Map<Integer, QuantEvent>> result = new HashMap<String, Map<Integer, QuantEvent>>();
        for (QuantEvent event : quantEvents)
        {
            String fraction = event.getFraction();
            Map<Integer, QuantEvent> scanEventMap = result.get(fraction);
            if (scanEventMap == null)
            {
                scanEventMap = new HashMap<Integer, QuantEvent>();
                result.put(fraction, scanEventMap);
            }
            scanEventMap.put(event.getScan(), event);
        }
        return result;
    }


    /**
     * Load quantitation events from a TSV file
     * @param eventFile
     * @return
     */
    public static List<QuantEvent> loadQuantEvents(File eventFile)
            throws IOException
    {
        TabLoader loader;

        _log.debug("Loading tab-delimited event file " + eventFile.getAbsolutePath());
        loader = new TabLoader(eventFile);
        _log.debug("Loaded tab file");

        List<QuantEvent> result = new ArrayList<QuantEvent>();

        Map[] rowsAsMaps = (Map[])loader.load();
        //nothing inthe file.  Return empty list
        if (rowsAsMaps == null || rowsAsMaps.length == 0 ||
                (rowsAsMaps.length == 1 && !rowsAsMaps[0].containsKey("Peptide")))
        return result;
        for (Map row : rowsAsMaps)
        {
            String protein = null;
            if (row.containsKey("Protein"))
                protein = row.get("Protein").toString();
            String peptide = row.get("Peptide").toString();
//System.err.println(peptide);
            String fraction = row.get("Fraction").toString();
//System.err.println(row.get("ID"));
//            System.err.println(row.get("image_url"));
//            System.err.println(row.get("algratio"));
//            System.err.println(row.get("singlepeakratio"));
//            System.err.println(row.get("evalstatus"));
//            System.err.println(row.get("evalnotes"));
//            System.err.println(row.get("Peptide"));

            int charge = Integer.parseInt(row.get("Charge").toString());
            String modificationState = "";
            if (row.containsKey("Modifications") && row.get("Modifications") != null)
                modificationState = row.get("Modifications").toString();
            int scan = Integer.parseInt(row.get("Scan").toString());
            float mz = 0;
            if (row.containsKey("Mz"))
                mz = Float.parseFloat(row.get("Mz").toString());

            //dhmay adding 20090723
            float ratioOnePeak = MISSING_RATIO;
            if (row.containsKey("RatioOnePeak"))
                ratioOnePeak = Float.parseFloat(row.get("RatioOnePeak").toString());

            File spectrumFile = null;
            if (row.get("Spectrum") != null)
                spectrumFile = new File(row.get("Spectrum").toString());
            File scansFile = null;
            if (row.get("Scans") != null)
                scansFile = new File(row.get("Scans").toString());
            File file3D = null;
            if (row.get("3D") != null)
                file3D = new File(row.get("3D").toString());
            File intensitySumFile = null;
            if (row.get("intSumChart") != null)
                intensitySumFile = new File(row.get("intSumChart").toString());
            Float ratio  = Float.parseFloat(row.get("Ratio").toString());
            Float lightMz = 0f;
            if (row.containsKey("LightMz"))
                lightMz = Float.parseFloat(row.get("LightMz").toString());
            Float heavyMz  = 0f;
            if (row.containsKey("HeavyMz"))
                heavyMz = Float.parseFloat(row.get("HeavyMz").toString());
            Float lightIntensity  = Float.parseFloat(row.get("LightInt").toString());
            Float heavyIntensity  = Float.parseFloat(row.get("HeavyInt").toString());
            int firstLightQuantScan = Integer.parseInt(row.get("LightFirstScan").toString());
            int lastLightQuantScan = Integer.parseInt(row.get("LightLastScan").toString());
            int firstHeavyQuantScan = Integer.parseInt(row.get("HeavyFirstScan").toString());
            int lastHeavyQuantScan = Integer.parseInt(row.get("HeavyLastScan").toString());
            int quantCurationStatus = QuantEvent.CURATION_STATUS_UNKNOWN;
            try
            {
                quantCurationStatus = parseCurationStatusString(row.get("QuantCuration").toString());
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Warning: problem loading curation status",e);
            }
            int idCurationStatus = QuantEvent.CURATION_STATUS_UNKNOWN;
            try
            {
                idCurationStatus = parseCurationStatusString(row.get("IDCuration").toString());
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Warning: problem loading curation status",e);
            }

            List<Integer> otherEventScans = new ArrayList<Integer>();
            List<Float> otherEventMZs = new ArrayList<Float>();
            if (row.get("OtherEventScans") != null && row.get("OtherEventMZs") != null)
            {
                otherEventScans =  MS2ExtraInfoDef.parseIntListString(row.get("OtherEventScans").toString());
                List<String> otherEventMZsStrings =
                        MS2ExtraInfoDef.parseStringListString(row.get("OtherEventMZs").toString());
                for (String mzString : otherEventMZsStrings)
                    otherEventMZs.add(Float.parseFloat(mzString));
            }

            float peptideProphet = 0;
            if (row.get("PProphet") != null)
                peptideProphet = Float.parseFloat(row.get("PProphet").toString());
            String comment = null;
            if (row.get("Comment") != null)
                comment = row.get("Comment").toString();

            QuantEvent quantEvent = new QuantEvent(protein,  peptide,  fraction,
                    charge, modificationState, scan, mz,
                    spectrumFile, scansFile, file3D, intensitySumFile,
                    ratio,  lightMz, heavyMz,
                    lightIntensity,  heavyIntensity,
                    firstLightQuantScan,  lastLightQuantScan,
                    firstHeavyQuantScan,  lastHeavyQuantScan,
                    otherEventScans, otherEventMZs, peptideProphet,
                    quantCurationStatus, idCurationStatus, comment);
            quantEvent.setRatioOnePeak(ratioOnePeak);

            if (row.get("AssessmentStatus") != null && row.get("AssessmentStatus").toString().length() > 0)
            {
                String assessmentDesc = "";
                if (row.get("AssessmentDesc") != null)
                    assessmentDesc = row.get("AssessmentDesc").toString();
                QuantEventAssessor.QuantEventAssessment assessment =
                        new QuantEventAssessor.QuantEventAssessment(
                                QuantEventAssessor.parseAssessmentCodeString(row.get("AssessmentStatus").toString()),
                                assessmentDesc);
                if (row.get("AssessmentComponents") != null)
                {
                    String bitmapString = row.get("AssessmentComponents").toString();
                    if (bitmapString.length() == QuantEventAssessor.flagReasonCodes.length)
                    {
                        boolean[] flagBitmap = new boolean[QuantEventAssessor.flagReasonCodes.length];
                        for (int i=0; i<bitmapString.length(); i++)
                            flagBitmap[i] = (bitmapString.charAt(i) == '1');
                        assessment.setFlagBitmap(flagBitmap);
                    }
                }
                quantEvent.setAlgorithmicAssessment(assessment);
            }

            result.add(quantEvent);
        }

        return result;
    }

    public static String convertCurationStatusToString(int curationStatus)
    {
        switch (curationStatus)
        {
            case CURATION_STATUS_GOOD:
                return "Good";
            case CURATION_STATUS_BAD:
                return "Bad";
            case CURATION_STATUS_RATIO_ONEPEAK:
                return "OnePeakRatio";
            default:
                return "Unknown";
        }
    }

    public static int parseCurationStatusString(String curationStatusString)
    {
        if ("Good".equals(curationStatusString))
            return CURATION_STATUS_GOOD;
        else if ("Bad".equals(curationStatusString))
            return CURATION_STATUS_BAD;
        else if ("OnePeakRatio".equals(curationStatusString))
            return CURATION_STATUS_RATIO_ONEPEAK;
        else return CURATION_STATUS_UNKNOWN;
    }

    public static void saveQuantEventsToTSV(Collection<QuantEvent> quantEvents,
            File outTsvFile, boolean showProteinColumn, boolean show3DPlots)
            throws IOException
    {
        PrintWriter outTsvPW = new PrintWriter(outTsvFile);
        _log.debug("QuantEvent.saveQuantEventsToTSV");        
        writeHeader(null, outTsvPW, showProteinColumn, show3DPlots);
        for (QuantEvent quantEvent : quantEvents)
        {
            outTsvPW.println(quantEvent.createOutputRow(null, false, showProteinColumn, show3DPlots));
            outTsvPW.flush();
        }
        writeFooterAndClose(null, outTsvPW);
    }
    

    public static void writeHeader(PrintWriter outHtmlPW, PrintWriter outTsvPW,
                            boolean showProteinColumn, boolean show3DPlots)
    {
        _log.debug("QuantEvent.writeHeader");
        if (outHtmlPW != null)
        {
            outHtmlPW.println(HtmlGenerator.createDocumentBeginning("Quantitative Events"));
            outHtmlPW.print("<table border=\"1\"><tr>");
            for (String columnHeader : dataColumnNames)
            {
                if (shouldShowColumn(columnHeader, showProteinColumn, show3DPlots))
                    outHtmlPW.print("<th>" + columnHeader + "</th>");
            }
            outHtmlPW.println("</tr>");
            outHtmlPW.flush();
        }
        if (outTsvPW != null)
        {
            boolean firstColumnShown = true;
            for (int i=0; i < dataColumnNames.length; i++)
            {
                String columnHeader = dataColumnNames[i];
                if (shouldShowColumn(columnHeader, showProteinColumn, show3DPlots))
                {
                    if (!firstColumnShown)
                        outTsvPW.print("\t");
                    outTsvPW.print(columnHeader);
                    firstColumnShown = false;
                }
            }
            outTsvPW.println();
            outTsvPW.flush();
        }
    }
    



    public static boolean shouldShowColumn(String columnName, boolean showProteinColumn, boolean show3DPlots)
    {
        if (!showProteinColumn && "Protein".equals(columnName))
            return false;
        if (!show3DPlots && "3D".equals(columnName))
            return false;
        return true;
    }

    public static void writeFooterAndClose(PrintWriter outHtmlPW, PrintWriter outTsvPW)
    {
        if (outHtmlPW != null)
        {
            outHtmlPW.println("</table>");
            outHtmlPW.flush();
            outHtmlPW.close();
        }
        if (outTsvPW != null)
            outTsvPW.close();
    }
    
    public void addQuantCurationStatusListener(ActionListener listener)
    {
        quantCurationStatusListeners.add(listener);
    }

    public void addAlgAssessmentStatusListener(ActionListener listener)
    {
        algAssessmentStatusListeners.add(listener);
    }

    /**
     *
     * @return true for light, false for heavy
     */
    public boolean idIsLight()
    {
        return mz == lightMz;
    }


    public String getProtein()
    {
        return protein;
    }

    public void setProtein(String protein)
    {
        this.protein = protein;
    }

    public String getPeptide()
    {
        return peptide;
    }

    public void setPeptide(String peptide)
    {
        this.peptide = peptide;
    }

    public String getFraction()
    {
        return fraction;
    }

    public void setFraction(String fraction)
    {
        this.fraction = fraction;
    }

    public int getCharge()
    {
        return charge;
    }

    public void setCharge(int charge)
    {
        this.charge = charge;
    }

    public int getScan()
    {
        return scan;
    }

    public void setScan(int scan)
    {
        this.scan = scan;
    }

    public File getSpectrumFile()
    {
        return spectrumFile;
    }

    public void setSpectrumFile(File spectrumFile)
    {
        this.spectrumFile = spectrumFile;
    }

    public File getScansFile()
    {
        return scansFile;
    }

    public void setScansFile(File scansFile)
    {
        this.scansFile = scansFile;
    }

    public File getFile3D()
    {
        return file3D;
    }

    public void setFile3D(File file3D)
    {
        this.file3D = file3D;
    }

    public float getRatio()
    {
        return ratio;
    }

    public void setRatio(float ratio)
    {
        this.ratio = ratio;
    }

    public float getLightIntensity()
    {
        return lightIntensity;
    }

    public void setLightIntensity(float lightIntensity)
    {
        this.lightIntensity = lightIntensity;
    }

    public float getHeavyIntensity()
    {
        return heavyIntensity;
    }

    public void setHeavyIntensity(float heavyIntensity)
    {
        this.heavyIntensity = heavyIntensity;
    }

    public int getFirstLightQuantScan()
    {
        return firstLightQuantScan;
    }

    public void setFirstLightQuantScan(int firstLightQuantScan)
    {
        this.firstLightQuantScan = firstLightQuantScan;
    }

    public int getLastLightQuantScan()
    {
        return lastLightQuantScan;
    }

    public void setLastLightQuantScan(int lastLightQuantScan)
    {
        this.lastLightQuantScan = lastLightQuantScan;
    }

    public int getFirstHeavyQuantScan()
    {
        return firstHeavyQuantScan;
    }

    public void setFirstHeavyQuantScan(int firstHeavyQuantScan)
    {
        this.firstHeavyQuantScan = firstHeavyQuantScan;
    }

    public int getLastHeavyQuantScan()
    {
        return lastHeavyQuantScan;
    }

    public void setLastHeavyQuantScan(int lastHeavyQuantScan)
    {
        this.lastHeavyQuantScan = lastHeavyQuantScan;
    }

    public int getQuantCurationStatus()
    {
        return quantCurationStatus;
    }

    /**
     * In addition to changing the status, updates any listeners
     * @param quantCurationStatus
     */
    public void setQuantCurationStatus(int quantCurationStatus)
    {
        this.quantCurationStatus = quantCurationStatus;
        for (ActionListener listener: quantCurationStatusListeners)
        {
            listener.actionPerformed(new ActionEvent(this, 0, "" + quantCurationStatus));
        }
    }

    public float getLightMz()
    {
        return lightMz;
    }

    public void setLightMz(float lightMz)
    {
        this.lightMz = lightMz;
    }

    public float getHeavyMz()
    {
        return heavyMz;
    }

    public void setHeavyMz(float heavyMz)
    {
        this.heavyMz = heavyMz;
    }

    public float getPeptideProphet()
    {
        return peptideProphet;
    }

    public void setPeptideProphet(float peptideProphet)
    {
        this.peptideProphet = peptideProphet;
    }

    public float getMz()
    {
        return mz;
    }

    public void setMz(float mz)
    {
        this.mz = mz;
    }

    public File getIntensitySumFile()
    {
        return intensitySumFile;
    }

    public void setIntensitySumFile(File intensitySumFile)
    {
        this.intensitySumFile = intensitySumFile;
    }

    public int getIdCurationStatus()
    {
        return idCurationStatus;
    }

    public void setIdCurationStatus(int idCurationStatus)
    {
        this.idCurationStatus = idCurationStatus;
    }

    public String getComment()
    {
        return comment;
    }

    public void setComment(String comment)
    {
        this.comment = comment;
    }

    public float calcLightNeutralMass()
    {
        return (lightMz - Spectrum.HYDROGEN_ION_MASS) * charge;
    }

    public float calcHeavyNeutralMass()
    {
        return (heavyMz - Spectrum.HYDROGEN_ION_MASS) * charge;
    }

    public String getModificationState()
    {
        return modificationState;
    }

    public void setModificationState(String modificationState)
    {
        this.modificationState = modificationState;
    }

    public float getRatioOnePeak()
    {
        return ratioOnePeak;
    }

    public void setRatioOnePeak(float ratioOnePeak)
    {
        this.ratioOnePeak = ratioOnePeak;
    }

    public QuantEventAssessor.QuantEventAssessment getAlgorithmicAssessment()
    {
        return algorithmicAssessment;
    }

    /**
     * Notifies any algorithmic assessment listeners
     * NOTE!  This will not be triggered if the algorithmicAssessment itself changes!
     * @param algorithmicAssessment
     */
    public void setAlgorithmicAssessment(QuantEventAssessor.QuantEventAssessment algorithmicAssessment)
    {
        this.algorithmicAssessment = algorithmicAssessment;
        for (ActionListener listener: algAssessmentStatusListeners)
        {
            listener.actionPerformed(new ActionEvent(this, 0, "" + algorithmicAssessment.getStatus()));
        }
    }

    public static class ScanAscComparator implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            if (o1.getScan() > o2.getScan())
                return 1;
            if (o1.getScan() < o2.getScan())
                return -1;
            return 0;
        }
    }

    public static class RatioAscComparator implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            if (o1.getRatio() > o2.getRatio())
                return 1;
            if (o1.getRatio() < o2.getRatio())
                return -1;
            return 0;
        }
    }

    public static class PeptideSequenceAscComparator implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            return o1.getPeptide().compareTo(o2.getPeptide());
        }
    }

    public static class FractionAscComparator implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            return o1.getFraction().compareTo(o2.getFraction());
        }
    }

    /**
     *  sort by peptide, then fraction, then charge, then modifications.
     * This is somewhat special-purpose, for ProteinQuantSummaryFrame, maybe should be moved there
     */
    public static class ProteinPeptideFractionChargeModificationsRatioAscComparator
            implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            float diff = o1.getProtein().compareTo(o2.getProtein());
            if (diff == 0)
                diff = o1.getPeptide().compareTo(o2.getPeptide());
            if (diff == 0)
                diff = o1.getFraction().compareTo(o2.getFraction());
            if (diff == 0)
                diff = o1.getCharge() - o2.getCharge();
            if (diff == 0)
                diff = o1.getModificationState().compareTo(o2.getModificationState());
            if (diff == 0)
                diff = o1.getRatio() - o2.getRatio();
            if (diff > 0)
                return 1;
            else if (diff < 0)
                return -1;
            return 0;
        }
    }

    /**
     *  sort by peptide, then fraction, then charge, then modifications.
     * This is somewhat special-purpose, for ProteinQuantSummaryFrame, maybe should be moved there
     */
    public static class PeptideSequenceAscFractionAscChargeModificationsAscRatioAscComparator
            implements Comparator<QuantEvent>
    {
        public int compare(QuantEvent o1, QuantEvent o2)
        {
            float diff = o1.getPeptide().compareTo(o2.getPeptide());
            if (diff == 0)
                diff = o1.getFraction().compareTo(o2.getFraction());
            if (diff == 0)
                diff = o1.getCharge() - o2.getCharge();
            if (diff == 0)
                diff = o1.getModificationState().compareTo(o2.getModificationState());            
            if (diff == 0)
                diff = o1.getRatio() - o2.getRatio();
            if (diff > 0)
                return 1;
            else if (diff < 0)
                return -1;
            return 0;
        }
    }

    public Feature getSourceFeature()
    {
        return sourceFeature;
    }

    public void setSourceFeature(Feature sourceFeature)
    {
        this.sourceFeature = sourceFeature;
    }

    public List<QuantEvent> getOtherEvents()
    {
        return otherEvents;
    }

    public void setOtherEvents(List<QuantEvent> otherEvents)
    {
        this.otherEvents = otherEvents;
    }

    public static class QuantEventPropertiesTable extends JTable
    {
        DefaultTableModel model = new DefaultTableModel(0, 2)
            {
                //all cells uneditable
                public boolean isCellEditable(int row, int column)
                {
                    return false;
                }

                public Class getColumnClass(int columnIndex)
                {
                    return String.class;
                }
            };

        public QuantEventPropertiesTable()
        { 
            setModel(model);
            getColumnModel().getColumn(0).setHeaderValue(TextProvider.getText("PROPERTY_LOWERCASE"));
            getColumnModel().getColumn(1).setHeaderValue(TextProvider.getText("VALUE_LOWERCASE"));


        }

        public QuantEventPropertiesTable(QuantEvent quantEvent)
        {
            this();
            displayQuantEvent(quantEvent);
        }

        public void displayQuantEvent(QuantEvent quantEvent)
        {
            clearProperties();
            Map<String, String> propMap = quantEvent.getNameValueMapNoCharts();
            for (String propName : QuantEvent.dataColumnNames)
            {
                if (propMap.containsKey(propName))
                    addPropertyToModel(propName, propMap.get(propName));
            }

            addPropertyToModel("Light Mass", "" + quantEvent.calcLightNeutralMass());
            addPropertyToModel("Heavy Mass", "" + quantEvent.calcHeavyNeutralMass());

            String lightOrHeavyID = "Light";
            if (Math.abs(quantEvent.getHeavyMz() - quantEvent.getMz()) < 0.25)
                lightOrHeavyID = "Heavy";
            //this is to accommodate legacy files without MZ
            if (quantEvent.getMz() == 0)
                lightOrHeavyID = "Unknown";
            addPropertyToModel("ID Light/Heavy", lightOrHeavyID);
        }

        //show tooltip with contents of cells
        public Component prepareRenderer(TableCellRenderer
                renderer,
                                         int rowIndex, int vColIndex)
        {
            Component c = super.prepareRenderer(renderer, rowIndex, vColIndex);
            if (c instanceof JComponent)
            {
                JComponent jc = (JComponent)c;
                jc.setToolTipText((String)getValueAt(rowIndex, vColIndex));
            }
            return c;
        }

        /**
         * Remove all properties from table
         */
        public void clearProperties()
        {
            for (int i=model.getRowCount()-1; i>=0; i--)
            {
                model.setValueAt(null, i, 0);
                model.setValueAt(null, i, 1);
            }
            model.setRowCount(0);
        }
        /**
         * Update the properties table model, adding a property value
         * @param propertyName
         * @param propertyValue
         */
        public void addPropertyToModel(Object propertyName, Object propertyValue)
        {
            int numRows = model.getRowCount();
            model.setRowCount(numRows + 1);
            model.setValueAt(propertyName, numRows, 0);
            model.setValueAt(propertyValue, numRows, 1);
        }
    }


}
