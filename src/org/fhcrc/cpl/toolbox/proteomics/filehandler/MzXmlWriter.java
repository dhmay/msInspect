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
package org.fhcrc.cpl.toolbox.proteomics.filehandler;

import java.math.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import net.sourceforge.sashimi.schemaRevision.mzXML21.*;
import org.w3c.dom.Node;
import org.w3c.dom.Attr;
import org.apache.xmlbeans.XmlOptions;
import org.apache.xmlbeans.GDuration;
import org.systemsbiology.jrap.stax.SoftwareInfo;
import org.systemsbiology.jrap.stax.MZXMLFileInfo;
import org.systemsbiology.jrap.stax.ParentFile;
import org.systemsbiology.jrap.stax.MSInstrumentInfo;
import org.fhcrc.cpl.toolbox.proteomics.MassCalibrationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.datastructure.Pair;


/**
 * A restrictive wrapper for writing mzXml files based on a single run.  We take advantage of XmlBeans to build
 * the structure of the file, and to build individual spectrum_queries (for features)
 * and search_summaries (for modifications), but we stitch the XmlBeans XML output for features together by
 * hand, writing out to a file as we go, so that we don't have to hold the whole structure
 * in memory.
 * A special wrinkle of mzXML files is that the file must contain offsets for finding individual
 * scans, and even for finding the index that holds those offsets.  For that reason we have to
 * maintain a _currentFilePosition variable that tells us how many bytes we've written so far.
 *
 * Also, XmlBeans likes namespaces.  A lot.  I can't figure out how to get XmlBeans to stop prefixing
 * every tag with <mzx:>.  This causes problems for jrap.  So I'm manually stripping out mzx: every time
 * I write anything to the file.  This could conceivably cause problems.
 *
 * For extra special fun, we can calibrate the spectra masses as we write out the file.  That is,
 * for every m/z value we write out, perform a specified linear transformation on it.  We do
 * NOT calibrate MS/MS precursor masses here.  That can be done elsewhere.  The spectra can't
 * be done elsewhere, since they're only soft-referenced in Java, and if we try to do them
 * in place in memory the changes will get lost.
 *
 * dhmay, 2008-12-31: changing the scan-writing behavior to preserve the numbering of scans as
 * received.  Previously, scans were re-numbered starting with 1, because our mzXML-parsing software
 * had trouble with non-sequential scans.  This has been addressed.
 */
public class MzXmlWriter
{
    static Logger _log = Logger.getLogger(MzXmlWriter.class);

    //doc representation
    protected MsRunDocument _xmlBeansMsRunDoc = null;
    //Index representation
    protected MzXMLDocument.MzXML.Index _xmlBeansIndex = null;
    //MsRun representation
    protected MsRunDocument.MsRun _xmlBeansMsRun = null;

    //MSRun that we want to write out
    MSRun _run = null;

    //Strings of xml representing the structure before and after the scan content
    protected String _documentPrefix = null;
    protected String _documentPostscript = null;

    //encapsulates printing options for all fragments
    protected XmlOptions _optionsForPrinting = null;

    //current position in the output file
    protected long _currentFilePosition = 0;

    //attribute name indicating precursor mz has been corrected
    public static final String PRECURSORMZ_ATTR_MSINSPECT_CORRECTED = "msInspect_corrected";

    public static final double UNSET_WAVELENGTH_OR_OFFSET = 999999999;
    public static final double UNSET_INTENSITY_SCALE = -1;

    //for mass recalibration
    protected Pair<Integer, Pair<Double, Double>>[] _massCalibrationParameters = null;

    //Should we calibrate precursor masses?
    protected boolean _shouldCalibratePrecursorMz = false;

    //should we calibrate the spectra themselves?  If not, and if calibration
    //parameters are set, only calibrate precursor masses
    protected boolean _shouldCalibrateSpectra = false;

    //for intensity scaling
    protected double intensityScaleFactor = UNSET_INTENSITY_SCALE;

    //optionally exclude MS1 scans
    protected boolean excludeMS1Scans = false;

    /**
     * Constructor creates the XmlBeans representing the shell of a mzXml document, and
     * creates the "prefix" and "postscript" strings representing that shell
     */
    public MzXmlWriter()
    {
        //initialize basic document structure and printing options

        _xmlBeansMsRunDoc = MsRunDocument.Factory.newInstance();

        //Construct generic document structure
        _xmlBeansMsRun = _xmlBeansMsRunDoc.addNewMsRun();

        //oddly, there's no Index element in the xsd under MsRun, even though that's where we want it
        MzXMLDocument fakeMzXmlDoc = MzXMLDocument.Factory.newInstance();
        _xmlBeansIndex = fakeMzXmlDoc.addNewMzXML().addNewIndex();

        //set printing options for xml fragments
        _optionsForPrinting = new XmlOptions();
        _optionsForPrinting.setSaveOuter();
        _optionsForPrinting.setSavePrettyPrint();
        _optionsForPrinting.setSavePrettyPrintOffset(0);
    }

    /**
     * Create doc structure, populate features and modifications
     * @param run
     */
    public MzXmlWriter(MSRun run)
    {
        this();
        setRun(run);
    }

    /**
     * This should only be called after the run is set.  This creates the whole XML shell around
     * the scans and index, which thankfully go right next to each other
     */
    protected void createDocumentShellXml()
    {
        //add a sentinel node that tells us where to split the document to insert scans,
        Node runNode = _xmlBeansMsRun.getDomNode();
        Node scanLocationNode = runNode.getOwnerDocument().createElement("SENTINEL_SCAN_LOCATION");
        runNode.appendChild(scanLocationNode);

        //create and break up the xml that defines the document structure
        String documentShell = _xmlBeansMsRunDoc.xmlText(_optionsForPrinting);
        String[] halves = documentShell.split("<SENTINEL_SCAN_LOCATION[^\\/]*\\/>");
        if (halves.length != 2)
        {
            _log.error("Failed to create document shell for writing");
            return;
        }

        _documentPrefix = removeNamespaceColon("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" + halves[0]);
        _documentPostscript = removeNamespaceColon(halves[1]);

        //remove our dummy node
        runNode.removeChild(scanLocationNode);
    }

    /**
     * Does all the processing to create an XmlBeans structure that
     * represents the run.  Doesn't do anything with the scans.
     * if restricted is true, the scan count is set to the passed-in numScans,
     * otherwise it's gathered from the run
     */
    public void buildDocStructure(int firstScanNum, int lastScanNum, boolean restricted)
    {
        MZXMLFileInfo fileInfo = _run.getHeaderInfo();

        int numMs1ScansToWrite = _run.getScanCount();
        int numMs2ScansToWrite = 0;
        int numMs3ScansToWrite = 0;
        if (_run.getMS2Scans() != null)
            numMs2ScansToWrite = _run.getMS2Scans().length;
        if (_run.getMS3Scans() != null)
            numMs3ScansToWrite = _run.getMS3Scans().length;

        if (restricted)
        {
            //since we're given scan numbers, not indexes, we don't immediately know how many
            //scans we're actually writing
            int firstScanIndex = _run.getIndexForScanNum(firstScanNum);
            int lastScanIndex = _run.getIndexForScanNum(lastScanNum);
            numMs1ScansToWrite = lastScanIndex-firstScanIndex+1;
            if (excludeMS1Scans)
                numMs1ScansToWrite = 0;

            numMs2ScansToWrite = 0;
            if (_run.getMS2Scans() != null)
            {
                for (int i=0; i<_run.getMS2Scans().length; i++)
                {
                    int currentMS2ScanNum = _run.getMS2Scans()[i].getNum();
                    if (currentMS2ScanNum >= firstScanNum)
                    {
                        if (currentMS2ScanNum <= lastScanNum)
                        {
                            numMs2ScansToWrite++;
                        }
                        else
                            break;
                    }
                }
            }
            numMs3ScansToWrite = 0;
            if (_run.getMS3Scans() != null)
            {
                for (int i=0; i<_run.getMS3Scans().length; i++)
                {
                    int currentMS3ScanNum = _run.getMS3Scans()[i].getNum();
                    if (currentMS3ScanNum >= firstScanNum)
                    {
                        if (currentMS3ScanNum <= lastScanNum)
                        {
                            numMs3ScansToWrite++;
                        }
                        else
                            break;
                    }
                }
            }
        }
        _xmlBeansMsRun.setScanCount(BigInteger.valueOf(numMs1ScansToWrite + numMs2ScansToWrite +
                                                       numMs3ScansToWrite));

        MsRunDocument.MsRun.DataProcessing dataProcessing = _xmlBeansMsRun.addNewDataProcessing();
        List<SoftwareInfo> jrapSoftwareInfoList = fileInfo.getDataProcessing().getSoftwareUsed();
        for (int i=0; i<jrapSoftwareInfoList.size(); i++)
        {
            SoftwareInfo jrapSoftwareInfo = jrapSoftwareInfoList.get(i);
            SoftwareDocument.Software xmlBeansSoftware = dataProcessing.addNewSoftware();
            xmlBeansSoftware.setName(jrapSoftwareInfo.name);
            xmlBeansSoftware.setType(SoftwareDocument.Software.Type.Enum.forString(jrapSoftwareInfo.type));
            xmlBeansSoftware.setVersion(jrapSoftwareInfo.version);
        }

        List<ParentFile> jrapParentFileArray = fileInfo.getParentFiles();
        for (int i=0; i<jrapParentFileArray.size(); i++)
        {
            ParentFile jrapParentFile = jrapParentFileArray.get(i);
            MsRunDocument.MsRun.ParentFile xmlBeansParentFile = _xmlBeansMsRun.addNewParentFile();
            //is this right?! ParentFile doesn't seem to have a Name
            xmlBeansParentFile.setFileName(jrapParentFile.getURI());
            xmlBeansParentFile.setFileSha1(jrapParentFile.getSha1());
            xmlBeansParentFile.setFileType(MsRunDocument.MsRun.ParentFile.FileType.Enum.forString(jrapParentFile.getType()));
        }

        //for some reason we don't always seem to capture this, so null-checking
        MSInstrumentInfo jrapInstrumentInfo = fileInfo.getInstrumentInfo();
        if (jrapInstrumentInfo != null)
        {
            MsRunDocument.MsRun.MsInstrument xmlBeansMsInstrument = _xmlBeansMsRun.addNewMsInstrument();
            MsRunDocument.MsRun.MsInstrument.MsManufacturer xmlBeansManufacturer =
                    xmlBeansMsInstrument.addNewMsManufacturer();
            xmlBeansManufacturer.setCategory("msManufacturer");            
            xmlBeansManufacturer.setValue(jrapInstrumentInfo.getManufacturer());
            OntologyEntryType xmlBeansMsModel = xmlBeansMsInstrument.addNewMsModel();
            xmlBeansMsModel.setValue(jrapInstrumentInfo.getModel());
            SoftwareInfo jrapInstrumentSoftwareInfo = jrapInstrumentInfo.getSoftwareInfo();
            if (jrapInstrumentSoftwareInfo != null)
            {
                SoftwareDocument.Software xmlBeansSoftware = dataProcessing.addNewSoftware();
                xmlBeansSoftware.setName(jrapInstrumentSoftwareInfo.name);
                xmlBeansSoftware.setType(SoftwareDocument.Software.Type.Enum.forString(jrapInstrumentSoftwareInfo.type));
                xmlBeansSoftware.setVersion(jrapInstrumentSoftwareInfo.version);
                xmlBeansMsInstrument.setSoftware(xmlBeansSoftware);
            }
        }

        //TODO: I'm not at all sure that this is correct, but these items are required
        //and they're not stored on the run
        String runStartTime = "PT0.0S";
        String runEndTime = "PT0.0S";

        try
        {
            MSRun.MSScan firstScan = _run.getScan(0);
            if (restricted)
            {
                int firstScanIndex = _run.getIndexForScanNum(firstScanNum);
                if (firstScanIndex < 0)
                    firstScanIndex = -firstScanIndex;
                firstScan = _run.getScan(firstScanIndex);

            }
            runStartTime = firstScan.getRetentionTime();

            MSRun.MSScan lastScan = _run.getScan(_run.getScanCount()-1);
            //if restricting scans, then the last scan isn't the last scan of the run
            int lastScanIndex = _run.getIndexForScanNum(lastScanNum);
            if (lastScanIndex < 0)
                lastScanIndex = -lastScanIndex;
            if (restricted)
                lastScan = _run.getScan(lastScanIndex);
            runEndTime = lastScan.getRetentionTime();
        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Warning:  Unable to get first and last scan information from run.  Defaulting run start/end times to 0.");
        }
        _xmlBeansMsRun.setStartTime(new GDuration(runStartTime));
        _xmlBeansMsRun.setEndTime(new GDuration(runEndTime));

        createDocumentShellXml();
    }


    /**
     * setter for the run.
     * @param run
     */
    public void setRun(MSRun run)
    {
        _run = run;
    }

    /**
     * Write out the index
     * @param pw
     */
    public void writeIndex(PrintWriter pw)
    {
        long indexOffset = _currentFilePosition;

        String indexFragment = removeNamespaceColon(_xmlBeansIndex.xmlText(_optionsForPrinting));
        pw.print(indexFragment);

        _currentFilePosition += indexFragment.length();

        String indexOffsetFragment = removeNamespaceColon("<mzx:indexOffset>" + indexOffset + "</mzx:indexOffset>");
        pw.print(indexOffsetFragment);

        _currentFilePosition += indexOffsetFragment.length();

        pw.flush();
    }

    /**
     * Return the index of the next MS2 scan _after_ index oldMs2ScanIndex that falls into
     * the given range (if restrict is true)
     * @param oldMs2ScanIndex
     * @param ms2Scans
     * @param firstScan
     * @param lastScan
     * @param restrict
     * @return The index of the next MS2 scan for writing.  If none, return -1 
     */
    protected int queueNextMs2ScanIndex(int oldMs2ScanIndex, MSRun.MSScan[] ms2Scans,
                                        int firstScan, int lastScan,
                                        boolean restrict)
    {
        int result = -1;
        if (ms2Scans == null || ms2Scans.length==0)
            return result;
        int currentScanIndex = oldMs2ScanIndex+1;
        while (currentScanIndex < ms2Scans.length)
        {
            if (!restrict ||
                (ms2Scans[currentScanIndex].getNum() >= firstScan &&
                 ms2Scans[currentScanIndex].getNum() <= lastScan))
            {
                result = currentScanIndex;
                break;
            }
            currentScanIndex++;
        }
        return result;
    }

    /**
     * Return the index of the next MS3 scan _after_ index oldMs3ScanIndex that falls into
     * the given range (if restrict is true)
     * This could really be folded into queueNextMs2ScanIndex, but it's just easier not to
     * @param oldMs3ScanIndex
     * @param ms3Scans
     * @param firstScan
     * @param lastScan
     * @param restrict
     * @return The index of the next MS2 scan for writing.  If none, return -1
     */
    protected int queueNextMs3ScanIndex(int oldMs3ScanIndex, MSRun.MSScan[] ms3Scans,
                                        int firstScan, int lastScan,
                                        boolean restrict)
    {
        int result = -1;
        if (ms3Scans == null || ms3Scans.length==0)
            return result;
        int currentScanIndex = oldMs3ScanIndex+1;
        while (currentScanIndex < ms3Scans.length)
        {
            if (!restrict ||
                (ms3Scans[currentScanIndex].getNum() >= firstScan &&
                 ms3Scans[currentScanIndex].getNum() <= lastScan))
            {
                result = currentScanIndex;
                break;
            }
            currentScanIndex++;
        }
        return result;
    }


    /**
     * Write out scans immediately.  Either write out all scans, or just a subregion, depending
     * on the value of the restrict argument.  If restrict == false, ignore the int and float args
     *
     * @param pw
     * @param firstScan
     * @param lastScan
     * @param lowMz
     * @param highMz
     * @param restrict
     */
    public void writeScans(PrintWriter pw, int firstScan, int lastScan, float lowMz, float highMz,
                              boolean restrict)
    {
        _log.debug("writeScans start");

        if (_run == null)
            return;
        MSRun.MSScan[] ms2Scans = _run.getMS2Scans();
        MSRun.MSScan[] ms3Scans = _run.getMS3Scans();

        int nextMs2ScanIndex = queueNextMs2ScanIndex(-1, ms2Scans, firstScan, lastScan,
                                                              restrict);
        int nextMs3ScanIndex = queueNextMs3ScanIndex(-1, ms3Scans, firstScan, lastScan,
                                                              restrict);

        int ms1ScanCount = _run.getScanCount();

        _log.debug("MS1 scans: " + ms1ScanCount);
        if (ms2Scans !=null) _log.debug("MS2 scans: " + ms2Scans.length);
        if (ms3Scans !=null) _log.debug("MS3 scans: " + ms3Scans.length);


        int nextCalibChangeIndex = 0;
        double massCalibrationWavelength = UNSET_WAVELENGTH_OR_OFFSET;
        double massCalibrationOffset = UNSET_WAVELENGTH_OR_OFFSET;

        int nextCalibChangeScan = 0;
        if (_massCalibrationParameters != null &&
                (_shouldCalibrateSpectra || _shouldCalibratePrecursorMz))
        {
            nextCalibChangeScan = _massCalibrationParameters[0].first;
        }
        for (int i=0; i<ms1ScanCount; i++)
        {
            if (i % (ms1ScanCount / 100) == 0)
                ApplicationContext.setMessage((i * 100 / ms1ScanCount) + " % complete");
            MSRun.MSScan ms1Scan = _run.getScan(i);

            if (_massCalibrationParameters != null &&
                (_shouldCalibrateSpectra || _shouldCalibratePrecursorMz))
            {
                if (nextCalibChangeScan <= ms1Scan.getNum())
                {
                    massCalibrationWavelength = _massCalibrationParameters[nextCalibChangeIndex].second.first;
                    massCalibrationOffset = _massCalibrationParameters[nextCalibChangeIndex].second.second;

                    if (_massCalibrationParameters.length >= nextCalibChangeIndex+2)
                    {
                        nextCalibChangeIndex++;
                        nextCalibChangeScan = _massCalibrationParameters[nextCalibChangeIndex].first;
                    }
                }
            }

//System.err.println("*****MS1 scan: " + ms1Scan.getNum());

            //Write out all MS2 scans that should be written out _before_ this ms1 scan is
            //written.
            //check that the next ms2 scan exists,
            //and that its number is less than the next ms1 scan's number,
            //and that it's number is in the range we want to write
            while (nextMs2ScanIndex > -1 &&
                   ms1Scan.getNum() > ms2Scans[nextMs2ScanIndex].getNum())
            {
                MSRun.MSScan ms2Scan = ms2Scans[nextMs2ScanIndex];
//System.err.println("  *****MS2 scan: " + ms2Scan.getNum());
                //Write out all MS3 scans that should be written out _before_ this ms1 scan is
                //written.
                //check that the next ms3 scan exists,
                //and that its number is less than the next ms2 scan's number,
                //and that it's number is in the range we want to write
                while (nextMs3ScanIndex > -1 &&
                        ms2Scan.getNum() > ms3Scans[nextMs3ScanIndex].getNum())
                {
                    writeScan(ms3Scans[nextMs3ScanIndex], ms3Scans[nextMs3ScanIndex].getNum(), pw, lowMz,
                              highMz, restrict, massCalibrationWavelength, massCalibrationOffset);
                    nextMs3ScanIndex = queueNextMs3ScanIndex(nextMs3ScanIndex, ms3Scans,
                                                             firstScan, lastScan, restrict);
                }

                writeScan(ms2Scan, ms2Scan.getNum(), pw, lowMz, highMz, restrict,
                          massCalibrationWavelength, massCalibrationOffset);
                nextMs2ScanIndex = queueNextMs2ScanIndex(nextMs2ScanIndex, ms2Scans,
                                                         firstScan, lastScan, restrict);
            }

            if (restrict)
            {
                //if we're restricting the scan and mz window, need to check whether this scan
                //is in the window.  Since restriction is based on scan number, not index,
                //need to check the scan to get the scan number
                int scanNum = ms1Scan.getNum();
                if (scanNum < firstScan)
                    continue;
                //scan numbers are monotonically increasing with scan index, so if we're out of
                //range already, we're done
                if (scanNum > lastScan)
                    break;
            }


            //if we get here, this is a scan we want to write
            if (!excludeMS1Scans)
                writeScan(ms1Scan, ms1Scan.getNum(), pw, lowMz, highMz, restrict,
                      massCalibrationWavelength, massCalibrationOffset);

        }

        //write out any remaining MS2 scans with higher indexes than the ms1 scans
        //we're writing, and any interstitial MS3 scans.
        while (nextMs2ScanIndex > -1)
        {
            MSRun.MSScan ms2Scan = ms2Scans[nextMs2ScanIndex];
            while (nextMs3ScanIndex > -1 &&
                    ms2Scan.getNum() > ms3Scans[nextMs3ScanIndex].getNum())
            {
                writeScan(ms3Scans[nextMs3ScanIndex], ms3Scans[nextMs3ScanIndex].getNum(), pw, lowMz, highMz,
                          restrict,
                          massCalibrationWavelength, massCalibrationOffset);
                nextMs3ScanIndex = queueNextMs3ScanIndex(nextMs3ScanIndex, ms3Scans,
                        firstScan, lastScan, restrict);
            }
            //no recalibration for MS/MS scans
            writeScan(ms2Scan, ms2Scan.getNum(), pw, lowMz, highMz, restrict);
            nextMs2ScanIndex = queueNextMs2ScanIndex(nextMs2ScanIndex, ms2Scans,
                                                     firstScan, lastScan, restrict);
        }

        //write out any remaining MS3 scans after all ms1 and ms2 scans.
        //No calibration for MS3 scans
        while (nextMs3ScanIndex > -1)
        {
            writeScan(ms3Scans[nextMs3ScanIndex], ms3Scans[nextMs3ScanIndex].getNum(), pw,
                    lowMz, highMz, restrict);
            nextMs3ScanIndex = queueNextMs3ScanIndex(nextMs3ScanIndex, ms3Scans,
                                                     firstScan, lastScan, restrict);
        }

        ApplicationContext.setMessage("100% complete");
    }


    protected void writeScan(MSRun.MSScan scan, int scanNumber, PrintWriter pw,
                             float lowMz, float highMz,
                             boolean restrict)
    {
        writeScan(scan, scanNumber, pw, lowMz, highMz, restrict,
                  UNSET_WAVELENGTH_OR_OFFSET, UNSET_WAVELENGTH_OR_OFFSET);
    }

    /**
     * Add a MsScan representing the passed-in MSScan to the MsRun.  If restrict==false, then only
     * write out a subregion of the scan's mz values
     * Write out the XML for that search result
     * @param scan
     * @param scanNumber
     * @param pw
     * @param lowMz
     * @param highMz
     * @param restrict
     * @param wavelengthForCalibration
     * @param offsetForCalibration
     */
    public void writeScan(MSRun.MSScan scan, int scanNumber, PrintWriter pw, float lowMz, float highMz,
                             boolean restrict, double wavelengthForCalibration, double offsetForCalibration)
    {
        ScanDocument.Scan xmlBeansScan = _xmlBeansMsRun.addNewScan();

        //scan information

        //this often seems to be null
        if (scan.getScanType() != null)
        {
            //some mzXml files have been seen to have bad scanType values, i.e., values not
            //in the list of allowed values, which includes Full, SRM, CRM, etc.
            //In those cases, just trapping the exception and refusing to write it out.
            try
            {
                xmlBeansScan.setScanType(ScanDocument.Scan.ScanType.Enum.forString(
                                               scan.getScanType()));
            }
            catch (Exception e)
            {
                _log.debug("Error setting scan type, " + e.getMessage());
            }
        }
        xmlBeansScan.setNum(BigInteger.valueOf(scanNumber));
        xmlBeansScan.setMsLevel(BigInteger.valueOf(scan.getMsLevel()));

//Hack for "removing" scans with no peaks
//TODO: make this hack ACTUALLY, OPTIONALLY remove those scans        
//if (scan.getPeaksCount() == 0 || scan.getSpectrum().length == 0)
//{
//    xmlBeansScan.setMsLevel(BigInteger.valueOf(4));
//    System.err.println("setting mslevel of scan " + scan.getNum() + " to 4");
//}
        xmlBeansScan.setRetentionTime(new GDuration(scan.getRetentionTime()));

        //These will be overridden if we're restricting
        float basePeakMz = scan.getBasePeakMz();
        float basePeakIntensity = scan.getBasePeakIntensity();

        //peaks
        ScanDocument.Scan.Peaks xmlBeansPeaks = xmlBeansScan.addNewPeaks();

        float[][] spectrum = null;
        try
        {
            //getSpectrum() can fail with an NPE if a scan is empty.
System.err.println("BEFORE GETSPEC");            
            spectrum = scan.getSpectrum();
  System.err.println("AFTER GETSPEC!");            

        }
        catch (NullPointerException e)
        {
            spectrum = new float[2][0];
System.err.println("CAUGHT NPE!");            
        }

        float[] mzSpectrum = spectrum[0];
        float[] intensitySpectrum = spectrum[1];

        byte[] spectrumArray = null;
        int peaksCount = scan.getPeaksCount();

        //get the low and high mz values from the scan.  If we're restricting, these will be changed
        float scanLowMz = scan.getLowMz();
        float scanHighMz = scan.getHighMz();

        //These are the MZ range imposed on the scan by the machine.  In the case of restriction,
        //should be further constrained by the window we're restricting to
        float scanStartMz = scan.getStartMz();
        float scanEndMz = scan.getEndMz();

        Float[] totIonCurrent = new Float[1];
        totIonCurrent[0] = scan.getTotIonCurrent();

        float precursorMz = scan.getPrecursorMz();

        //scale intensities, if specified
        if (intensityScaleFactor != UNSET_INTENSITY_SCALE)
        {
            for (int i=0; i<intensitySpectrum.length; i++)
                intensitySpectrum[i] *= intensitySpectrum.length;
        }

        //calibrate masses, if specified
        if (wavelengthForCalibration != UNSET_WAVELENGTH_OR_OFFSET &&
                offsetForCalibration != UNSET_WAVELENGTH_OR_OFFSET)
        {
            if (_shouldCalibrateSpectra)
            {
                for (int i=0; i<mzSpectrum.length; i++)
                    mzSpectrum[i] +=
                            mzSpectrum[i] *
                                    (MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH -
                                            wavelengthForCalibration) -
                                    offsetForCalibration;
            }
            if (_shouldCalibratePrecursorMz)
            {
                precursorMz += precursorMz * (MassCalibrationUtilities.DEFAULT_THEORETICAL_MASS_WAVELENGTH -
                                            wavelengthForCalibration) -
                                    offsetForCalibration;
            }
        }


        if (xmlBeansScan.getMsLevel().intValue() == 2)
        {
            ScanDocument.Scan.PrecursorMz xmlBeansPrecursorMz =
                    xmlBeansScan.addNewPrecursorMz();
            xmlBeansPrecursorMz.setPrecursorCharge(BigInteger.valueOf(scan.getPrecursorCharge()));
            xmlBeansPrecursorMz.setPrecursorScanNum(BigInteger.valueOf(scan.getPrecursorScanNum()));
            xmlBeansPrecursorMz.setFloatValue(precursorMz);

            if (scan.isPrecursorMzCorrected())
            {
                Attr correctedAttribute = _xmlBeansMsRun.getDomNode().getOwnerDocument().createAttribute(
                        PRECURSORMZ_ATTR_MSINSPECT_CORRECTED);
                correctedAttribute.setValue("true");
                xmlBeansPrecursorMz.getDomNode().getAttributes().setNamedItem(correctedAttribute);
//                xmlBeansPrecursorMz.getDomNode().appendChild(correctedAttribute);
            }

        }

        if (restrict)
        {
            //this will hold the lowest and highest mz values that we find in the scan
            Pair<Float, Float> lowAndHighFoundMz =
                    new Pair<Float, Float>(lowMz, highMz);
            //this will hold the mz value of the peak with the highest intensity
            Pair<Float, Float> basePeakMzIntensity =
                    new Pair<Float, Float>(0f, 0f);
            spectrumArray = createRestrictedByteArrayForSpectrum(mzSpectrum, intensitySpectrum,
                                                                 lowMz, highMz, lowAndHighFoundMz,
                                                                 basePeakMzIntensity, totIonCurrent);
            //if we're restricting, need to recalculate peaksCount.  every 8 bytes in the
            //spectrum array represent one peak
            peaksCount = spectrumArray.length / 8;

            //record the actual found lowest and highest MZ values
            scanStartMz = Math.max(lowMz,scanLowMz);
            scanEndMz = Math.min(highMz,scanHighMz);

            //record the base (highest-intensity) peak mz and intensity
            basePeakMz = basePeakMzIntensity.first;
            basePeakIntensity = basePeakMzIntensity.second;

            //Set the scan _range_ appropriately.  The most restrictive window defined by
            //the scan's low and high values and the low and high values we're restricting to here
            scanLowMz = lowAndHighFoundMz.first;
            scanHighMz = lowAndHighFoundMz.second;
        }
        else
        {
            spectrumArray = createByteArrayForSpectrum(mzSpectrum, intensitySpectrum,
                                                       scan.getPeaksCount());
        }

        xmlBeansScan.setLowMz(scanLowMz);
        xmlBeansScan.setHighMz(scanHighMz);
        xmlBeansScan.setStartMz(scanStartMz);
        xmlBeansScan.setEndMz(scanEndMz);

        xmlBeansScan.setTotIonCurrent(totIonCurrent[0]);
        xmlBeansScan.setBasePeakMz(basePeakMz);
        xmlBeansScan.setBasePeakIntensity(basePeakIntensity);

        xmlBeansScan.setPeaksCount(BigInteger.valueOf(peaksCount));
        xmlBeansPeaks.setPrecision(BigInteger.valueOf(scan.getPrecision()));
        xmlBeansPeaks.setByteArrayValue(spectrumArray);


        try
        {
            String fragment = removeNamespaceColon(_xmlBeansMsRun.getScanArray(0).xmlText(_optionsForPrinting));
            MzXMLDocument.MzXML.Index.Offset scanOffset = _xmlBeansIndex.addNewOffset();
            scanOffset.setId(BigInteger.valueOf(scanNumber));
            scanOffset.setLongValue(_currentFilePosition);
            pw.print(fragment);
            _currentFilePosition += fragment.length();
            pw.flush();
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }

        _xmlBeansMsRun.removeScan(0);
    }

    /**
     * Creates a base-64 byte array, with alternating mz values and intensities, to represent
     * the two float arrays we store in MSScan.
     * If we got REALLY parsimonious wrt memory, we could change things so that we
     * don't even store this whole thing all at once.  That'd get hairy, though.
     * @return
     */
    protected byte[] createByteArrayForSpectrum(float[] mzSpectrum, float[]intensitySpectrum,
                                                int peaksCount)
    {
        byte[] result = new byte[peaksCount * 2 * 4];
        for (int i=0; i<mzSpectrum.length; i++)
        {
            //populate 4 bytes representing the mz
            fillByteArrayForFloat(mzSpectrum[i], result, 8*i);
            //populate 4 bytes representing the intensity
            fillByteArrayForFloat(intensitySpectrum[i], result, 8*i+4);
        }

        return result;
    }

    /**
     * Separate handling for restricted spectra, since we have to build an ArrayList of bytes
     * because we don't know in itially how many bytes we're writing out.
     * Populate lowAndHighFoundMz with the actual lowest and highest values we find
     * @param mzSpectrum
     * @param intensitySpectrum
     * @param lowMz
     * @param highMz
     * @param outLowAndHighFoundMz will contain the lowest and highest found mz values
     * @param outBasePeakMzIntensity will contain the mz and intensity of the base peak
     * @param outTotIonCurrent will contain the sum of all intensities as its first element (hack)
     * @return
     */
    protected byte[] createRestrictedByteArrayForSpectrum(float[] mzSpectrum,
                                                          float[]intensitySpectrum,
                                                          float lowMz, float highMz,
                                                          Pair outLowAndHighFoundMz,
                                                          Pair outBasePeakMzIntensity,
                                                          Float[] outTotIonCurrent)
    {
        ArrayList<Byte> resultList = new ArrayList<Byte>();

        //initialize the lowest and highest MZ values so that they're sure to get overwritten
        float lowestMz = highMz;
        float highestMz = lowMz;

        //initialize the base peak (peak with highest intensity) values so that they're sure to get
        //overwritten
        float basePeakMz = 0f;
        float basePeakIntensity = 0f;

        float totIonCurrent = 0f;

        byte[] byteHolder = new byte[4];

        for (int i=0; i<mzSpectrum.length; i++)
        {
            if (mzSpectrum[i] < lowMz)
                continue;
            if (mzSpectrum[i] > highMz)
                break;

            //populate 4 bytes representing the mz
            fillByteArrayForFloat(mzSpectrum[i], byteHolder, 0);
            resultList.add(byteHolder[0]);
            resultList.add(byteHolder[1]);
            resultList.add(byteHolder[2]);
            resultList.add(byteHolder[3]);
            //populate 4 bytes representing the intensity
            fillByteArrayForFloat(intensitySpectrum[i], byteHolder, 0);
            resultList.add(byteHolder[0]);
            resultList.add(byteHolder[1]);
            resultList.add(byteHolder[2]);
            resultList.add(byteHolder[3]);

            if (mzSpectrum[i] < lowestMz)
                lowestMz = mzSpectrum[i];
            if (mzSpectrum[i] > highestMz)
                highestMz = mzSpectrum[i];

            if (intensitySpectrum[i] > basePeakIntensity)
            {
                basePeakIntensity = intensitySpectrum[i];
                basePeakMz = mzSpectrum[i];
            }
            totIonCurrent += intensitySpectrum[i];
        }

        //TODO: This is probably very memory-inefficient
        byte[] result = new byte[resultList.size()];
        for (int i=0; i<resultList.size(); i++)
            result[i] = resultList.get(i);
        outLowAndHighFoundMz.first = Float.valueOf(lowestMz);
        outLowAndHighFoundMz.second = Float.valueOf(highestMz);

        outBasePeakMzIntensity.first = Float.valueOf(basePeakMz);
        outBasePeakMzIntensity.second = Float.valueOf(basePeakIntensity);

        outTotIonCurrent[0]=totIonCurrent;
        return result;
    }

    /**
     * Given a float, populate the given byte array with four bytes representing the float.
     * The first byte will go in at the offset given
     * @param floatVal
     * @param bytes
     */
    protected void fillByteArrayForFloat(float floatVal, byte[] bytes, int offset)
    {
        assert(bytes.length >= (offset + 4));

        int intBits = Float.floatToIntBits(floatVal);
        bytes[offset] = (byte)(intBits >> 24);
        bytes[offset + 1] = (byte)(intBits >> 16);
        bytes[offset + 2] = (byte)(intBits >> 8);
        bytes[offset + 3] = (byte)(intBits);
    }


    /**
     * Remove all occurrences of the namespace prefix from an xml fragment.
     * TODO: Boy, I'd sure like XmlBeans to do this for me.
     * @param input
     * @return
     */
    protected String removeNamespaceColon(String input)
    {
        return input.replaceAll("mzx:","");
    }



    /**
     * Write out the full document, with all modifications and features, to a file
     * @param file
     * @throws IOException
     */
    public void write(File file) throws IOException
    {
       _write(file, 0,0,0,0,false);
    }


    /**
     * Write out just a subregion
     * @param file
     * @param firstScan
     * @param lastScan
     * @param lowMz
     * @param highMz
     */
    public void writeSubregion(File file, int firstScan, int lastScan, float lowMz, float highMz)
            throws IOException
    {
        _write(file,firstScan,lastScan,lowMz,highMz,true);
    }

    /**
     * Workhorse writing method.  Can write just a subregion, or the whole file
     * @param file
     * @param firstScanNum
     * @param lastScanNum
     * @param lowMz
     * @param highMz
     */
    protected void _write(File file, int firstScanNum, int lastScanNum, float lowMz, float highMz,
                          boolean restrict)
            throws IOException
    {
        //build the XML structure of the document
        buildDocStructure(firstScanNum, lastScanNum, restrict);
        PrintWriter pw = new PrintWriter(file);
        _log.debug("Writing document start");
        printDocPrefix(pw);

        writeScans(pw, firstScanNum, lastScanNum, lowMz, highMz, restrict);
        _log.debug("Writing index");
        writeIndex(pw);
        _log.debug("Finishing document");
        printDocPostscript(pw);

        pw.flush();
        _log.debug("Done.");
    }

    public void printDocPrefix(PrintWriter pw)
    {
        _currentFilePosition = 0;
        pw.print(_documentPrefix);
        _currentFilePosition += _documentPrefix.length();
    }

    public void printDocPostscript(PrintWriter pw)
    {
        pw.print(_documentPostscript);
        _currentFilePosition += _documentPostscript.length();
    }

    public void setMassCalibrationParameters(Pair<Integer,Pair<Double,Double>>[] newParameters)
    {
        _massCalibrationParameters = newParameters;
        _shouldCalibrateSpectra = true;
    }

    public void setMassCalibrationWavelengthOffset(double massCalibrationWavelength,
                                                   double massCalibrationOffset)
    {
        _massCalibrationParameters = (Pair<Integer, Pair<Double,Double>>[]) new Pair[1];
        _massCalibrationParameters[0] =
                new Pair<Integer, Pair<Double,Double>>(0,
                        new Pair<Double,Double>(massCalibrationWavelength,
                                                massCalibrationOffset));
    }

    public boolean shouldCalibrateSpectra()
    {
        return _shouldCalibrateSpectra;
    }

    public void setShouldCalibrateSpectra(boolean _shouldCalibrateSpectra)
    {
        this._shouldCalibrateSpectra = _shouldCalibrateSpectra;
    }


    public boolean isShouldCalibratePrecursorMz()
    {
        return _shouldCalibratePrecursorMz;
    }

    public void setShouldCalibratePrecursorMasses(boolean _shouldCalibratePrecursorMasses)
    {
        this._shouldCalibratePrecursorMz = _shouldCalibratePrecursorMasses;
    }

    public double getIntensityScaleFactor()
    {
        return intensityScaleFactor;
    }

    public void setIntensityScaleFactor(double intensityScaleFactor)
    {
        this.intensityScaleFactor = intensityScaleFactor;
    }


    public boolean isExcludeMS1Scans()
    {
        return excludeMS1Scans;
    }

    public void setExcludeMS1Scans(boolean excludeMS1Scans)
    {
        this.excludeMS1Scans = excludeMS1Scans;
    }
}
