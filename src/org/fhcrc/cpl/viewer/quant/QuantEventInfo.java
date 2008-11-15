package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;

import javax.swing.*;
import javax.swing.table.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.event.MouseListener;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.ItemEvent;

/**
     * Holds all the information related to a quantitative event for display
 */
public class QuantEventInfo
{
    public static final int CURATION_STATUS_UNKNOWN = 0;
    public static final int CURATION_STATUS_GOOD = 1;
    public static final int CURATION_STATUS_BAD = 2;
    
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
    protected List<QuantEventInfo> otherEvents;

    protected Feature sourceFeature;

    protected String comment;


    protected int quantCurationStatus = CURATION_STATUS_UNKNOWN;
    protected int idCurationStatus = CURATION_STATUS_UNKNOWN;

    public QuantEventInfo(QuantEventInfo eventToCopy)
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
    }

    public QuantEventInfo(Feature feature, String fraction)
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
            otherEvents = new ArrayList<QuantEventInfo>();
            for (Spectrum.Peak peak : feature.comprised)
            {
                otherEvents.add(new QuantEventInfo((Feature) peak, fraction));
            }
        }

        init(protein, peptide, fraction, charge,
                MS2ExtraInfoDef.convertModifiedAminoAcidsMapToString(MS2ExtraInfoDef.getModifiedAminoAcidsMap(feature)),
                feature.getScan(), feature.getMz(),
                null, null, null, null,
                ratio, lightMz, heavyMz, lightIntensity, heavyIntensity,
                firstLightQuantScan, lastLightQuantScan, firstHeavyQuantScan, lastHeavyQuantScan,
                otherEvents,
                (float) MS2ExtraInfoDef.getPeptideProphet(feature),
                CURATION_STATUS_UNKNOWN, CURATION_STATUS_UNKNOWN, null);
    }

    public QuantEventInfo(Feature feature, String fraction, File spectrumFile, File scansFile,
                          File file3D, File intensitySumFile)
    {
        this(feature, fraction);
        this.spectrumFile = spectrumFile;
        this.scansFile = scansFile;
        this.file3D = file3D;
        this.intensitySumFile = intensitySumFile;
    }

    /**
     * This is a hack.  We're modeling other events in the tsv file as lists of scan numbers and
     * mz values.  So when we read that in, need to turn those into a list of events identical to
     * this one except scan and mz.
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
    public QuantEventInfo(String protein, String peptide, String fraction, int charge, String modificationState, int scan,
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
        otherEvents = new ArrayList<QuantEventInfo>();
        for (int i=0; i<otherEventScans.size(); i++)
        {
            QuantEventInfo otherEvent = new QuantEventInfo(this);
            otherEvent.otherEvents = null;
            otherEvent.scan = otherEventScans.get(i);
            otherEvent.mz = otherEventMzs.get(i);
            otherEvents.add(otherEvent);
        }

    }

    public QuantEventInfo(String protein, String peptide, String fraction, int charge, String modificationState,
                          int scan, float mz, File spectrumFile, File scansFile, File file3D,
                          File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<QuantEventInfo> otherEvents,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {
        init(protein, peptide, fraction, charge, modificationState, scan, mz, spectrumFile, scansFile, file3D, intensitySumFile,
                ratio,
                lightMz, heavyMz, lightIntensity, heavyIntensity, firstLightQuantScan, lastLightQuantScan,
                firstHeavyQuantScan, lastHeavyQuantScan, otherEvents,
                peptideProphet, quantCurationStatus, idCurationStatus, comment);
    }


    protected void init(String protein, String peptide, String fraction, int charge, String modificationState, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D, File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<QuantEventInfo> otherEvents,
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
            for (QuantEventInfo otherEvent : otherEvents)
                allScansAsStrings.add("" + otherEvent.getScan());
        }
        return MS2ExtraInfoDef.convertStringListToString(allScansAsStrings);
    }

    protected String convertOtherEventMZsToString()
    {
        List<String> allMzsAsStrings = new ArrayList<String>();
        if (otherEvents != null)
        {
            for (QuantEventInfo otherEvent : otherEvents)
                allMzsAsStrings.add("" + otherEvent.getMz());
        }
        return MS2ExtraInfoDef.convertStringListToString(allMzsAsStrings);
    }

    protected String createOutputRow(String outChartsRelativeDirPath, boolean isHtml, boolean showProteinColumn,
                                  boolean show3DColumn)
    {
        String spectrumFileString = spectrumFile.getAbsolutePath();
        String scansFileString = scansFile.getAbsolutePath();
        String fileString3D =  null == file3D ? null : file3D.getAbsolutePath();
        String intensitySumFileString = null == intensitySumFile ? null : intensitySumFile.getAbsolutePath();


        if (isHtml)
        {
            spectrumFileString = HtmlGenerator.createLink(outChartsRelativeDirPath + spectrumFile.getName(),"Spectrum");
            scansFileString = HtmlGenerator.createLink(outChartsRelativeDirPath + scansFile.getName(), "Scans");
            if (show3DColumn)
                fileString3D = HtmlGenerator.createLink(outChartsRelativeDirPath + file3D.getName(), "3D");
            intensitySumFileString =
                    HtmlGenerator.createLink(outChartsRelativeDirPath + intensitySumFile.getName(), "intSumChart");

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

        stringValuesForRow.add("" + convertOtherEventScansToString());
        stringValuesForRow.add("" + convertOtherEventMZsToString());

        stringValuesForRow.add(convertCurationStatusToString(quantCurationStatus));
        stringValuesForRow.add(convertCurationStatusToString(idCurationStatus));

        stringValuesForRow.add(null == comment ? "" : comment);


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
                    "OtherEventScans",
                    "OtherEventMZs",
                    "QuantCuration",
                    "IDCuration",
                    "Comment",
            };

    protected boolean[] columnsAreFiles = new boolean[]
            {
                    false, false,false,false, false, false,true,true,true,
                    false,false,false,false,false,false,false,false,false, false
            };


    /**
     * @param eventFile
     * @return
     */
    public static List<QuantEventInfo> loadQuantEvents(File eventFile)
            throws IOException
    {
        TabLoader loader;

        loader = new TabLoader(eventFile);

        List<QuantEventInfo> result = new ArrayList<QuantEventInfo>();

        Map[] rowsAsMaps = (Map[])loader.load();

        for (Map row : rowsAsMaps)
        {
            String protein = null;
            if (row.containsKey("Protein"))
                protein = row.get("Protein").toString();
            String peptide = row.get("Peptide").toString();
            String fraction = row.get("Fraction").toString();
            int charge = Integer.parseInt(row.get("Charge").toString());
            String modificationState = "";
            if (row.containsKey("Modifications"))
                modificationState = row.get("Modifications").toString();
            int scan = Integer.parseInt(row.get("Scan").toString());
            float mz = 0;
            if (row.containsKey("Mz"))
                mz = Float.parseFloat(row.get("Mz").toString());

            File spectrumFile = new File(row.get("Spectrum").toString());
            File scansFile = new File(row.get("Scans").toString());
            File file3D = null;
            if (row.containsKey("3D"))
                file3D = new File(row.get("3D").toString());
            File intensitySumFile = null;
            if (row.containsKey("intSumChart"))
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
            int quantCurationStatus = QuantEventInfo.CURATION_STATUS_UNKNOWN;
            try
            {
                quantCurationStatus = parseCurationStatusString(row.get("QuantCuration").toString());
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Warning: problem loading curation status",e);
            }
            int idCurationStatus = QuantEventInfo.CURATION_STATUS_UNKNOWN;
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

            QuantEventInfo quantEvent = new QuantEventInfo(protein,  peptide,  fraction,
                    charge, modificationState, scan, mz,
                    spectrumFile, scansFile, file3D, intensitySumFile,
                    ratio,  lightMz, heavyMz,
                    lightIntensity,  heavyIntensity,
                    firstLightQuantScan,  lastLightQuantScan,
                    firstHeavyQuantScan,  lastHeavyQuantScan,
                    otherEventScans, otherEventMZs, peptideProphet,
                    quantCurationStatus, idCurationStatus, comment);
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
        else return CURATION_STATUS_UNKNOWN;
    }

    public static void saveQuantEventsToTSV(Collection<QuantEventInfo> quantEvents,
            File outTsvFile, boolean showProteinColumn, boolean show3DPlots)
            throws IOException
    {
        PrintWriter outTsvPW = new PrintWriter(outTsvFile);
        writeHeader(null, outTsvPW, showProteinColumn, show3DPlots);
        for (QuantEventInfo quantEvent : quantEvents)
        {
            outTsvPW.println(quantEvent.createOutputRow(null, false, showProteinColumn, show3DPlots));
            outTsvPW.flush();
        }
        writeFooterAndClose(null, outTsvPW);
    }

    public static void writeHeader(PrintWriter outHtmlPW, PrintWriter outTsvPW,
                            boolean showProteinColumn, boolean show3DPlots)
    {

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

    protected static boolean shouldShowColumn(String columnName, boolean showProteinColumn, boolean show3DPlots)
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

    public void setQuantCurationStatus(int quantCurationStatus)
    {
        this.quantCurationStatus = quantCurationStatus;
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

    public static class ScanAscComparator implements Comparator<QuantEventInfo>
    {
        public int compare(QuantEventInfo o1, QuantEventInfo o2)
        {
            if (o1.getScan() > o2.getScan())
                return 1;
            if (o1.getScan() < o2.getScan())
                return -1;
            return 0;
        }
    }

    public static class RatioAscComparator implements Comparator<QuantEventInfo>
    {
        public int compare(QuantEventInfo o1, QuantEventInfo o2)
        {
            if (o1.getRatio() > o2.getRatio())
                return 1;
            if (o1.getRatio() < o2.getRatio())
                return -1;
            return 0;
        }
    }

    public static class PeptideSequenceAscComparator implements Comparator<QuantEventInfo>
    {
        public int compare(QuantEventInfo o1, QuantEventInfo o2)
        {
            return o1.getPeptide().compareTo(o2.getPeptide());
        }
    }

    public static class PeptideSequenceAscFractionAscChargeModificationsAscRatioAscComparator implements Comparator<QuantEventInfo>
    {
        public int compare(QuantEventInfo o1, QuantEventInfo o2)
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

    public List<QuantEventInfo> getOtherEvents()
    {
        return otherEvents;
    }

    public void setOtherEvents(List<QuantEventInfo> otherEvents)
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

        public QuantEventPropertiesTable(QuantEventInfo quantEvent)
        {
            this();
            displayQuantEvent(quantEvent);
        }

        public void displayQuantEvent(QuantEventInfo quantEvent)
        {
            clearProperties();
            Map<String, String> propMap = quantEvent.getNameValueMapNoCharts();
            for (String propName : QuantEventInfo.dataColumnNames)
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

    /**
     *
     */
    public static class QuantEventsSummaryTable extends JTable
    {
        protected List<Integer> shadedTableRows = new ArrayList<Integer>();
        protected List<QuantEventInfo> quantEvents = new ArrayList<QuantEventInfo>();

        DefaultTableModel model = new DefaultTableModel(0, 9)
        {
            //all cells uneditable
            public boolean isCellEditable(int row, int column)
            {
                if (column == 0)
                    return true;
                return false;
            }

            public Class getColumnClass(int columnIndex)
            {
                switch (columnIndex)
                {
                    case 0:
                        return Boolean.class;
                    case 6:
                        return JSlider.class;
                    default:
                        return String.class;
                }
            }
        };

        /**
         * Hide the checkbox column.  There's no undoing this
         */
        public void hideSelectionColumn()
        {
            this.removeColumn(getColumnModel().getColumn(0));
        }

        /**
         * Hide the checkbox column.  There's no undoing this
         */
        public void hideProteinColumn()
        {
            this.removeColumn(getColumnModel().getColumn(1));
        }

        public QuantEventsSummaryTable()
        {
            setModel(model);



            TableColumn checkboxColumn = getColumnModel().getColumn(0);
            checkboxColumn.setHeaderRenderer(new CheckBoxHeader(new SelectAllListener()));
//        checkboxColumn.setCellEditor(new DefaultCellEditor(new JCheckBox()));
            checkboxColumn.setPreferredWidth(20);
            checkboxColumn.setMaxWidth(20);

            TableColumn proteinColumn = getColumnModel().getColumn(1);
            proteinColumn.setHeaderValue("Protein");
            proteinColumn.setPreferredWidth(90);

            TableColumn peptideColumn = getColumnModel().getColumn(2);
            peptideColumn.setHeaderValue("Peptide");
            peptideColumn.setPreferredWidth(170);
            peptideColumn.setMinWidth(140);

            getColumnModel().getColumn(3).setHeaderValue("Charge");
            getColumnModel().getColumn(4).setHeaderValue("Probability");
            getColumnModel().getColumn(5).setHeaderValue("Ratio");
            getColumnModel().getColumn(6).setHeaderValue("Light");
            getColumnModel().getColumn(7).setHeaderValue("Heavy");
            
            TableColumn logRatioSliderColumn = getColumnModel().getColumn(8);
            logRatioSliderColumn.setHeaderValue("LogRatio");
            JSliderRenderer sliderRenderer = new JSliderRenderer();
            logRatioSliderColumn.setCellRenderer(sliderRenderer);
            logRatioSliderColumn.setPreferredWidth(280);
            logRatioSliderColumn.setMinWidth(100);

            getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        }

        protected Color altRowColor = new Color(235, 235, 235);
        /**
         * Shades alternate peptides in different colors.
         */
        public Component prepareRenderer(TableCellRenderer renderer, int row, int column)
        {
            Component c = super.prepareRenderer(renderer, row, column);
            if (isCellSelected(row, column) == false)
            {
                Color rowColor = UIManager.getColor("Table.background");
                if (shadedTableRows.contains(row))
                    rowColor = altRowColor;
                c.setBackground(rowColor);
                c.setForeground(UIManager.getColor("Table.foreground"));
            } else
            {
                c.setBackground(UIManager.getColor("Table.selectionBackground"));
                c.setForeground(UIManager.getColor("Table.selectionForeground"));
            }
            return c;
        }

        /**
         * Remove all properties from table
         */
        public void clearProperties()
        {
            while (model.getRowCount() > 0)
            {
                model.removeRow(0);
            }
        }

        public void displayEvents(List<QuantEventInfo> quantEvents)
        {
            clearProperties();
            this.quantEvents = quantEvents;
            shadedTableRows = new ArrayList<Integer>();
            boolean shaded = true;
            String previousPeptide = "";
            for (QuantEventInfo quantEvent : quantEvents)
            {
                if (!previousPeptide.equals(quantEvent.getPeptide()))
                {
                    shaded = !shaded;
                    previousPeptide = quantEvent.getPeptide();
                }

                int numRows = model.getRowCount();

                if (shaded)
                    shadedTableRows.add(numRows);                
                model.setRowCount(numRows + 1);
//            JCheckBox thisEventCheckBox = new JCheckBox();
                model.setValueAt(new Boolean(false), numRows, 0);
                model.setValueAt(quantEvent.getProtein(), numRows, 1);                
                model.setValueAt(quantEvent.getPeptide(), numRows, 2);
                model.setValueAt("" + quantEvent.getCharge(), numRows, 3);
                model.setValueAt("" + quantEvent.getPeptideProphet(), numRows, 4);
                model.setValueAt("" + quantEvent.getRatio(), numRows, 5);
                model.setValueAt("" + quantEvent.getLightIntensity(), numRows, 6);
                model.setValueAt("" + quantEvent.getHeavyIntensity(), numRows, 7);

                float ratioBound = 10f;
                float logRatioBounded =
                        (float) Math.log(Math.min(ratioBound, Math.max(1.0f / ratioBound, quantEvent.getRatio())));
                int logRatioIntegerizedHundredScale =
                        (int) (logRatioBounded * 100 / (2 * Math.log(ratioBound))) + 50;
                model.setValueAt(logRatioIntegerizedHundredScale, numRows, 8);
            }
        }

        public List<QuantEventInfo> getSelectedEvents()
        {
            List<QuantEventInfo> selectedQuantEvents = new ArrayList<QuantEventInfo>();
            for (int i=0; i<model.getRowCount(); i++)
            {
                Boolean isSelected = (Boolean) model.getValueAt(i, 0);
                if (isSelected)
                    selectedQuantEvents.add(quantEvents.get(i));
            }
            return selectedQuantEvents;
        }

        public class JSliderRenderer implements TableCellRenderer
        {
            protected JSlider slider = new JSlider();

            public JSliderRenderer()
            {
                slider.setMinimum(0);
                slider.setMaximum(100);
                slider.setPaintLabels(false);
                slider.setPaintTicks(false);
                slider.setMajorTickSpacing(25);
                slider.setPreferredSize(new Dimension(280, 15));
                slider.setPreferredSize(new Dimension(100, 15));
                slider.setToolTipText("Log ratio, bounded at 0.1 and 10");
            }

            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                           boolean hasFocus, int row, int column)
            {
                Integer val = (Integer)value;
                slider.setValue(val.intValue());
                return slider;
            }
        }


        class CheckBoxHeader extends JCheckBox
                implements TableCellRenderer, MouseListener
        {
            protected CheckBoxHeader rendererComponent;
            protected int column;
            protected boolean mousePressed = false;
            public CheckBoxHeader(ItemListener itemListener) {
                rendererComponent = this;
                rendererComponent.addItemListener(itemListener);
            }
            public Component getTableCellRendererComponent(
                    JTable table, Object value,
                    boolean isSelected, boolean hasFocus, int row, int column) {
                if (table != null) {
                    JTableHeader header = table.getTableHeader();
                    if (header != null) {
                        rendererComponent.setForeground(header.getForeground());
                        rendererComponent.setBackground(header.getBackground());
                        rendererComponent.setFont(header.getFont());
                        header.addMouseListener(rendererComponent);
                    }
                }
                setColumn(column);
                rendererComponent.setText("Check All");
                setBorder(UIManager.getBorder("TableHeader.cellBorder"));
                return rendererComponent;
            }
            protected void setColumn(int column) {
                this.column = column;
            }
            public int getColumn() {
                return column;
            }
            protected void handleClickEvent(MouseEvent e) {
                if (mousePressed) {
                    mousePressed=false;
                    JTableHeader header = (JTableHeader)(e.getSource());
                    JTable tableView = header.getTable();
                    TableColumnModel columnModel = tableView.getColumnModel();
                    int viewColumn = columnModel.getColumnIndexAtX(e.getX());
                    int column = tableView.convertColumnIndexToModel(viewColumn);

                    if (viewColumn == this.column && e.getClickCount() == 1 && column != -1) {
                        doClick();
                    }
                }
            }
            public void mouseClicked(MouseEvent e) {
                handleClickEvent(e);
                ((JTableHeader)e.getSource()).repaint();
            }
            public void mousePressed(MouseEvent e) {
                mousePressed = true;
            }
            public void mouseReleased(MouseEvent e) {
            }
            public void mouseEntered(MouseEvent e) {
            }
            public void mouseExited(MouseEvent e) {
            }
        }

        class SelectAllListener implements ItemListener
        {
            public void itemStateChanged(ItemEvent e) {
                Object source = e.getSource();
                if (source instanceof AbstractButton == false) return;
                boolean checked = e.getStateChange() == ItemEvent.SELECTED;
                for(int x = 0, y = getRowCount(); x < y; x++)
                {
                    setValueAt(new Boolean(checked),x,0);
                }
            }
        }

    }
}
