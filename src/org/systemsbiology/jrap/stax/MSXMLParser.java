/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) MSXMLParser.java * Author: * Ning Zhang
 * nzhang@systemsbiology.org
 * ****************************************************************************** * * *
 * This software is provided ``AS IS'' and any express or implied * *
 * warranties, including, but not limited to, the implied warranties of * *
 * merchantability and fitness for a particular purpose, are disclaimed. * * In
 * no event shall the authors or the Institute for Systems Biology * * liable
 * for any direct, indirect, incidental, special, exemplary, or * *
 * consequential damages (including, but not limited to, procurement of * *
 * substitute goods or services; loss of use, data, or profits; or * * business
 * interruption) however caused and on any theory of liability, * * whether in
 * contract, strict liability, or tort (including negligence * * or otherwise)
 * arising in any way out of the use of this software, even * * if advised of
 * the possibility of such damage. * * *
 * ******************************************************************************/

package org.systemsbiology.jrap.stax;

import javax.xml.stream.XMLStreamReader;
import java.io.*;
import java.util.*;



/**
 * A generic utility class for reading an MSXML file in a random access fashion
 * and utilizing a stored scan index for fast reads.
 *
 * tholzman 200911xx
 * -Inserting changes for compatibility with sequential scan iterators 
 *
 */
public final class MSXMLParser
{
    /** The file we are in charge of reading */
    protected String fileName = null;

    /** The indexes */
    protected Map<Integer, Long> offsets;
    protected int maxScan;
    protected long chrogramIndex;

    protected boolean isXML = false;
    protected boolean isML = false;

    /* TAH Nov 2009 */
    int currentScanIndex;
    EndPatternStringIterator epsi = null;
    public void setEpsi(EndPatternStringIterator e) {
       this.epsi = e;
    }
    public EndPatternStringIterator getEpsi(){
        return epsi;
    }

    public static boolean isMzXML(String fn) {
        return fn.indexOf("mzXML") != -1;
    }

    private void commonInits(String fileName) {
        if (isMzXML(fileName))
            isXML = true;
	    else
	        isML = true;
        this.fileName = fileName;
    }

    private void sequentialInits() throws IOException {
       String leftPat="<msRun", rightPat=">", attr = "scanCount";
       String nextLeftPat="<scan",nextRightPat="</peaks>";
       if(isML){
           leftPat = "<spectrumList"; rightPat = ">"; attr = "count";
           nextLeftPat="<spectrum"; nextRightPat="</spectrum>";
       }
       setEpsi(new EndPatternStringIterator(leftPat,rightPat,fileName));
       XMLStreamReader xmlSR = epsi.xmlsrNext();
       try {xmlSR.next();} catch (Exception e) {throw new IOException(e);};
       maxScan = Integer.parseInt(xmlSR.getAttributeValue(null,attr));
       offsets = new HashMap<Integer,Long>();
       getEpsi().setLeftPatStr(nextLeftPat);
       getEpsi().setRightPatStr(nextRightPat);
       currentScanIndex = 0;
    }

    private void randomInits() {
       //using IndexParser get indexes
       IndexParser indexParser = new IndexParser(fileName);
       indexParser.parseIndexes();
       offsets = indexParser.getOffsetMap();
       maxScan = indexParser.getMaxScan();
       chrogramIndex = indexParser.getChrogramIndex();
    }

    public MSXMLParser(String fn, boolean isSequential) throws IOException {
       commonInits(fn);
       if(isSequential) {
          sequentialInits();
       } else {
          randomInits();
       }
    }

    public MSXMLParser(String fileName)  {
       commonInits(fileName);
       randomInits();
    }

    /* end TAH */

    /**this gives back the file header (info before scan)
     *@return the file header info (MZXMLFileInfo)
     */
    public MZXMLFileInfo rapFileHeader()
   {
	FileHeaderParser fileParser = new FileHeaderParser(fileName);
	fileParser.parseFileHeader();
	return (fileParser.getInfo());
    }

    /**
     *@return a scan header object without peaks information.
     * dhmay changing 20091021 to set the scanOffset on the returned scanHeader.  This was earlier behavior that
     * was removed by Ning.  Replacing, because it increases efficiency quite a bit for calling code.
     */
    public ScanHeader rapHeader(int scanNumber)
    {
        FileInputStream fileIN = null;
        long scanOffset = -1;
        try
        {
            fileIN = new FileInputStream(fileName);
            scanOffset = getScanOffset(scanNumber);
            if (scanOffset == -1)
            {
                return null;
            }

            fileIN.skip(scanOffset);
        }
        catch (Exception e)
        {
            System.out.println("File exception:" + e);
            e.printStackTrace();
        }
        ScanHeader scanHeader = null;
        if(isXML)
        {
            ScanAndHeaderParser headerParser = new ScanAndHeaderParser();
            headerParser.setIsScan(false);
            headerParser.setFileInputStream(fileIN);
            headerParser.parseScanAndHeader();

            closeFile(fileIN);
            scanHeader = headerParser.getHeader();
        }
        else
        {
            MLScanAndHeaderParser headerParser = new MLScanAndHeaderParser();
            headerParser.setIsScan(false);
            headerParser.setFileInputStream(fileIN);
            headerParser.parseMLScanAndHeader();
            closeFile(fileIN);
            scanHeader = headerParser.getHeader();
        }
        scanHeader.setScanOffset(scanOffset);
        return scanHeader;
    }

    /* TAH Nov 2009 */


    public ScanHeader nextHeader()
    {
       ScanHeader scanHeader = null;
       StringBuffer curScanInfo = epsi.next();
       if(curScanInfo == null || curScanInfo.length() == 0 || curScanInfo.charAt(0) != '<') return null;
       currentScanIndex++;
       if(isXML) {
           ScanAndHeaderParser headerParser = new ScanAndHeaderParser();
           headerParser.setIsScan(false);
           try{headerParser.parseScanAndHeader(epsi.xmlsrCur());}
           catch (Exception e){};
           scanHeader = headerParser.getHeader();
       } else {
           MLScanAndHeaderParser headerParser = new MLScanAndHeaderParser();
           headerParser.setIsScan(false);
           try{headerParser.parseMLScanAndHeader(epsi.xmlsrCur());}
           catch (Exception e){};
           scanHeader = headerParser.getHeader();
       }
       offsets.put(currentScanIndex,epsi.getFilePos());
       scanHeader.setScanOffset(epsi.getFilePos());
       return scanHeader;
    }
    /* end TAH      */

    private void closeFile(FileInputStream fileIN) {
        if(fileIN != null) {
            try {
                fileIN.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Read a particular scan from a MSXML file and return a generic Scan object
     * with it's data. Note: scanNumbers are 1-based, so scanNumber must be at
     * least 1 and be not greater than getScanCount() + 1
     *@return a scan object. It has all the infomation in a scanheader object and also
     * peaks information that doesn't included in scanHeader object.
     */
    public Scan rap(int scanNumber)
    {
	FileInputStream fileIN = null;
	try
	    {
		fileIN = new FileInputStream(fileName);
		long scanOffset = getScanOffset(scanNumber);
        if (scanOffset == -1)
		{
			return null;
		}

		fileIN.skip(scanOffset);
	    } catch (Exception e)
	    {
		System.out.println("File exception:" + e);
		e.printStackTrace();
	    }

	if(isXML)
	    {
		ScanAndHeaderParser scanParser = new ScanAndHeaderParser();
		scanParser.setIsScan(true);
		scanParser.setFileInputStream(fileIN);
		scanParser.parseScanAndHeader();
		closeFile(fileIN);
        return ( scanParser.getScan());
	    }
	else
	    {
		MLScanAndHeaderParser scanParser = new MLScanAndHeaderParser();
		scanParser.setIsScan(true);
		scanParser.setFileInputStream(fileIN);
		scanParser.parseMLScanAndHeader();
		
		closeFile(fileIN);
		return (scanParser.getScan());
	    }
    }

    /**
     * Get the total number of scans in the mzXMLfile handled by this parser.
     *
     * @return The number of scans.
     */
    public int getScanCount()    /* TAH Nov 2009 */
    {
        if(epsi != null) { //sequential scan, index scan hasn't been run
           return maxScan;
        } else {
           return offsets.size();
        }
    }         /* end TAH */

    public int getMaxScanNumber()
    {
		return maxScan;
    }


    /**
     *get scan offset, scan number is 1 based.
     */

    public long getScanOffset(int scanNumber)
    {
	if (scanNumber > 0 && offsets.containsKey(scanNumber))
        {
            return ((offsets.get(scanNumber)).longValue());
        } else
        {
            return (-1);
        }
    }
}