package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.HtmlGenerator;
import org.fhcrc.cpl.toolbox.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

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
    protected List<Integer> otherEventScans;
    protected List<Float> otherEventMZs;

    protected String comment;


    protected int quantCurationStatus = CURATION_STATUS_UNKNOWN;
    protected int idCurationStatus = CURATION_STATUS_UNKNOWN;


    public QuantEventInfo(Feature feature, String fraction)
    {
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

        init(protein, peptide, fraction, charge, feature.getScan(), feature.getMz(),
                null, null, null, null,
                ratio, lightMz, heavyMz, lightIntensity, heavyIntensity,
                firstLightQuantScan, lastLightQuantScan, firstHeavyQuantScan, lastHeavyQuantScan, feature.comprised,
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

    public QuantEventInfo(String protein, String peptide, String fraction, int charge, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D,
                          File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<Integer> otherEventScans, List<Float> otherEventMZs,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {
        init(protein, peptide, fraction, charge, scan, mz, spectrumFile, scansFile, file3D, intensitySumFile,
                ratio,
                lightMz, heavyMz, lightIntensity, heavyIntensity, firstLightQuantScan, lastLightQuantScan,
                firstHeavyQuantScan, lastHeavyQuantScan, otherEventScans, otherEventMZs,
                peptideProphet, quantCurationStatus, idCurationStatus, comment);
    }

    protected void init(String protein, String peptide, String fraction, int charge, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D, File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          Spectrum.Peak[] otherFeaturesAsPeaks,
                          float peptideProphet, int curationStatus, int idCurationStatus,
                          String comment)
    {
        List<Integer> otherEventScans = new ArrayList<Integer>();
        List<Float> otherEventMZs = new ArrayList<Float>();

        if (otherFeaturesAsPeaks != null)
        {
            for (Spectrum.Peak otherFeatureAsPeak : otherFeaturesAsPeaks)
            {
                otherEventScans.add(otherFeatureAsPeak.getScan());
                otherEventMZs.add(otherFeatureAsPeak.getMz());
            }
        }

        init(protein, peptide, fraction, charge, scan, mz, spectrumFile, scansFile, file3D, intensitySumFile,
                ratio,
                lightMz, heavyMz, lightIntensity, heavyIntensity, firstLightQuantScan, lastLightQuantScan,
                firstHeavyQuantScan, lastHeavyQuantScan, otherEventScans, otherEventMZs,
                peptideProphet, curationStatus, idCurationStatus, comment);
    }

    protected void init(String protein, String peptide, String fraction, int charge, int scan,
                          float mz, File spectrumFile, File scansFile, File file3D, File intensitySumFile,
                          float ratio, float lightMz, float heavyMz,
                          float lightIntensity, float heavyIntensity,
                          int firstLightQuantScan, int lastLightQuantScan,
                          int firstHeavyQuantScan, int lastHeavyQuantScan,
                          List<Integer> otherEventScans, List<Float> otherEventMZs,
                          float peptideProphet, int quantCurationStatus, int idCurationStatus,
                          String comment)
    {
        this.protein = protein;
        this.peptide = peptide;
        this.fraction = fraction;
        this.charge = charge;
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
        this.otherEventScans = otherEventScans;
        this.otherEventMZs = otherEventMZs;
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
        for (int scan : otherEventScans)
            allScansAsStrings.add("" + scan);
        return MS2ExtraInfoDef.convertStringListToString(allScansAsStrings);
    }

    protected String convertOtherEventMZsToString()
    {
        List<String> allMZsAsStrings = new ArrayList<String>();
        for (float mz : otherEventMZs)
            allMZsAsStrings.add("" + mz);
        return MS2ExtraInfoDef.convertStringListToString(allMZsAsStrings);
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
                    charge,  scan, mz,
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

    public List<Integer> getOtherEventScans()
    {
        return otherEventScans;
    }

    public void setOtherEventScans(List<Integer> otherEventScans)
    {
        this.otherEventScans = otherEventScans;
    }

    public List<Float> getOtherEventMZs()
    {
        return otherEventMZs;
    }

    public void setOtherEventMZs(List<Float> otherEventMZs)
    {
        this.otherEventMZs = otherEventMZs;
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
    
}
