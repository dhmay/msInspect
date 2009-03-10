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

import java.io.*;
import java.util.*;



/**
 * A generic utility class for reading an MSXML file in a random access fashion
 * and utilizing a stored scan index for fast reads.
 *
 */
public final class MSXMLParser
{
    /** The file we are in charge of reading */
    protected String fileName = null;

    protected boolean debug = false;

   
    /** The indexes */
    protected Map<Integer, Long> offsets;
    protected int maxScan;
    protected long chrogramIndex;


    protected boolean isXML = false;
    protected boolean isML = false;

    public MSXMLParser(String fileName) {
	this.fileName = fileName;

	if(fileName.indexOf("mzXML") != -1)
	    isXML = true;
	else
	    isML = true;

	//using IndexParser get indexes
	IndexParser indexParser = new IndexParser(fileName);
	indexParser.parseIndexes();
	offsets = indexParser.getOffsetMap();
    if (debug)
    {
        System.err.println("Offsets map created with " + offsets.size() + " entries");
//        for (int offset : offsets.keySet())
//            System.err.println("\t" + offset + "=" + offsets.get(offset));
    }
    maxScan = indexParser.getMaxScan();
	chrogramIndex = indexParser.getChrogramIndex();
    }

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
     */
    public ScanHeader rapHeader(int scanNumber)
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
		ScanAndHeaderParser headerParser = new ScanAndHeaderParser();
		headerParser.setIsScan(false);
		headerParser.setFileInputStream(fileIN);
		headerParser.parseScanAndHeader();

		return (headerParser.getHeader());
	    }
	else
	    {
//System.err.println("***rapHeader, scan " + scanNumber);
        MLScanAndHeaderParser headerParser = new MLScanAndHeaderParser();
		headerParser.setIsScan(false);
		headerParser.setFileInputStream(fileIN);
		headerParser.parseMLScanAndHeader();
        return (headerParser.getHeader());
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
//System.err.println("Rap 1");
	try
	    {
		fileIN = new FileInputStream(fileName);
		long scanOffset = getScanOffset(scanNumber);
		if (scanOffset == -1)
		{
			return null;
		}
//  System.err.println("Rap to offset " + scanOffset);

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

		return ( scanParser.getScan());
	    }
	else
	    {
        if (debug)
            System.err.println("parsing mzML scan spectra for scan " + scanNumber);
        MLScanAndHeaderParser scanParser = new MLScanAndHeaderParser();
		scanParser.setIsScan(true);
		scanParser.setFileInputStream(fileIN);
		scanParser.parseMLScanAndHeader();
		return (scanParser.getScan());
	    }
    }

    /**
     * Get the total number of scans in the mzXMLfile handled by this parser.
     *
     * @return The number of scans.
     */
    public int getScanCount()
    {
		return offsets.size();
    }

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
