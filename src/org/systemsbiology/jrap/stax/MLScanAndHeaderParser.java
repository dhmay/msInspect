/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) MLScanHeaderParser.java * Author: * Ning Zhang
 * nzhang@systemsbiology.org
 * ****************************************************************************** 
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
 
//support mzML parsing
package org.systemsbiology.jrap.stax;
import java.io.*;
import javax.xml.stream.*;
import javax.xml.stream.events.*;
import java.util.*;
import java.nio.*;
import java.util.zip.*;

/**
 * dhmay changing 2009/03/09, because mzML 1.1 changes the way scan IDs are stored.  They are now stored in
 * the "id" attribute of "spectrum", which is being used to contain multiple name-value pairs; the
 * name of the name-value pair containing the scan number is "scan", so I'm knocking off everything but that pair.
 * Also changing to cobble together Scan.massIntensityList from its components, which was missed earlier.
 * Also calling tmpScanHeader.setRetentionTime(), which was previously not set.
 */
public class MLScanAndHeaderParser{

    public ScanHeader tmpScanHeader;
    public Scan tmpScan;

    
    FileInputStream fileIN = null;

    boolean isScan = false;
 

    public void setIsScan(boolean isScan)
    {
	this.isScan = isScan;
    }

 
   
    public void setFileInputStream(FileInputStream in)
    {
	try{
	    this.fileIN = in;
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
   }

    public ScanHeader getHeader()
    {
	return tmpScanHeader;
    }

    public Scan getScan()
    {
	return tmpScan;
    }

    /**
     * mzML 1.1 changes the way scan IDs are stored.  They are now stored in
     * the "id" attribute of "spectrum", which is being used to contain multiple name-value pairs; the
     * name of the name-value pair containing the scan number is "scan", so I'm knocking off everything but that pair.
     * @param idString
     * @return
     */
    protected int parseScanNumberFromSpectrumIdField(String idString)
    {
        if (idString.contains("scan="))
            idString = idString.substring(idString.indexOf("scan=") + "scan=".length());
        if (idString.contains(" "))
            idString = idString.substring(0, idString.indexOf(" "));
        return Integer.parseInt(idString);
    }

    public void parseMLScanAndHeader()
    {
	try{
	    XMLInputFactory inputFactory = XMLInputFactory.newInstance();
	    XMLStreamReader xmlSR = inputFactory.createXMLStreamReader(fileIN,"ISO-8859-1");

	    boolean inSpectrum = false;
	    boolean inPeaks = false;
	    boolean isPeaks = false;
	    String elementName = null;
	    String attriName = null;
	    String attriValue = null;
	  
	    StringBuffer peaksBuffer = null;
	    int count = 0;

	    while(xmlSR.hasNext())
		{
		    int event = xmlSR.next();
		    if(event == xmlSR.START_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("spectrum"))
				{
				    inSpectrum = true;
				    count=0;
				    tmpScanHeader = new ScanHeader();
                    //dhmay changing 2009/03/09.  mzML 1.1 changes the way scan IDs are stored
                    tmpScanHeader.setNum(parseScanNumberFromSpectrumIdField(getStringValue(xmlSR, "id")));
				    tmpScanHeader.setPeaksCount(getIntValue(xmlSR, "defaultArrayLength"));
				    
				}
			    if(elementName.equals("cvParam"))
				{
				    attriName = xmlSR.getAttributeValue(null,"name");
				    if(inSpectrum)
					{
					    
					    if(attriName.equals("ms level"))
						tmpScanHeader.setMsLevel(getIntValue(xmlSR, "value"));
					    if(attriName.equals("centroid mass spectrum"))
						tmpScanHeader.setCentroided(1);
					    if(attriName.equals("base peak m/z"))
						tmpScanHeader.setBasePeakMz(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("base peak intensity"))
						tmpScanHeader.setBasePeakIntensity(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("total ion current"))
						tmpScanHeader.setTotIonCurrent(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("lowest m/z value"))
						tmpScanHeader.setStartMz(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("highest m/z value"))
						tmpScanHeader.setEndMz(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("scan m/z lower limit"))
						tmpScanHeader.setLowMz(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("scan m/z upper limit"))
						tmpScanHeader.setHighMz(getFloatValue(xmlSR, "value"));
					    if(attriName.equals("filter string"))
						tmpScanHeader.setFilterLine(getStringValue(xmlSR,"value"));
					    if(attriName.equals("full scan"))
						tmpScanHeader.setScanType("full scan");
					    if(attriName.equals("positive scan"))
						tmpScanHeader.setPolarity("+");
					    if(attriName.equals("scan time"))
						{
						    String timeType = xmlSR.getAttributeValue(null,"unitName");
						    double rt = Double.parseDouble(xmlSR.getAttributeValue(null, "value"));
						    if(timeType.equals("minute"))
							    tmpScanHeader.setRT(rt*60);
						    else
                                tmpScanHeader.setRT(rt);

                            //dhmay adding for backward compatibility.  Probably this should be rewired so that
                            //getDoubleRetentionTime just accesses the rt variable, but I don't want to sort out
                            //that tangle
                            tmpScanHeader.setRetentionTime("PT" + rt + "S");
                        }
					    //precursor
					    if(attriName.equals("m/z"))
						tmpScanHeader.setPrecursorMz(getFloatValue(xmlSR,"value"));
					    if(attriName.equals("intensity"))
						tmpScanHeader.setPrecursorIntensity(getFloatValue(xmlSR,"value"));
					    if(attriName.equals("charge"))
						tmpScanHeader.setPrecursorCharge(getIntValue(xmlSR,"value"));
					    if(attriName.equals("collision energy"))
						tmpScanHeader.setCollisionEnergy(getFloatValue(xmlSR,"value"));
					}
				    if(inPeaks)
					{
					    if(attriName.equals("64-bit float"))
						{
						    if(count == 1)
							tmpScanHeader.setMassPrecision(64);	
						    if(count == 2)
							tmpScanHeader.setIntenPrecision(64);
						}
					    if(attriName.equals("32-bit float"))
						{
						    if(count == 1)
							tmpScanHeader.setMassPrecision(32);
						    if(count == 2)
							tmpScanHeader.setIntenPrecision(32);
						}
					    if(attriName.equals("no compression"))
						{
						    if(count == 1)
							{
							    tmpScanHeader.setMassCompressionType("None");
							}
						    if(count == 2)
							{
							    tmpScanHeader.setIntenCompressionType("None");
							}
						}
					    if(attriName.indexOf("zlib") != -1)
						{
						    if(count == 1)
							{
							    tmpScanHeader.setMassCompressionType("zlib");
							}
						    if(count == 2)
							{
							    tmpScanHeader.setIntenCompressionType("zlib");
							}
						}
					}
					    
						
				}

			    if(elementName.equals("binaryDataArrayList"))
				{
				    if(isScan)
					{
					    tmpScan = new Scan();
					    tmpScan.setHeader(tmpScanHeader);
					}
				    else
					throw new XMLStreamException("ScanHeaderEndFoundException");
				}
				    
			    if(elementName.equals("binaryDataArray"))
				{
				    inPeaks = true;
				    count++;
				    //System.out.println("count "+count);
				    
				    if(count == 1)
					tmpScanHeader.setMassCompressedLen(getIntValue(xmlSR, "encodedLength"));
				    if(count == 2)
                    tmpScanHeader.setIntenCompressedLen(getIntValue(xmlSR, "encodedLength"));

				   
				}
			    if(elementName.equals("binary"))
				{
				    inPeaks = false;
				    isPeaks = true;
				    peaksBuffer = new StringBuffer();
				}
			}
		    if(event == xmlSR.CHARACTERS)
			{
			    if(isPeaks)
				peaksBuffer.append(xmlSR.getText());
			}
		    if(event ==xmlSR.END_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("spectrumDescription"))
				{
				    inSpectrum = false;
				}
			   
			    if(elementName.equals("binary"))
				{
//System.err.println("End binary, count=" + count + ", peaks=" + peaksBuffer.toString());
                    getPeaks(peaksBuffer.toString(),count);
				    isPeaks = false;
				    peaksBuffer = null;
				}
			    if(elementName.equals("binaryDataArrayList"))
				{
				    throw new XMLStreamException("ScanEndFoundException");
				}  
			}
		}
    }
	catch(Exception e)
	    {
		String exception1=e.getMessage();
		if(!exception1.equals("ScanHeaderEndFoundException"))
		    {
			if(!exception1.equals("ScanEndFoundException"))
			    e.printStackTrace();
		    }
	    }
    }

    public String getStringValue(XMLStreamReader xmlSR, String name)
    {
	String value="";
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value="";
	    else
		value=xmlSR.getAttributeValue(null,name);
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public int getIntValue(XMLStreamReader xmlSR, String name)
    {
	int value=-1;
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value = -1;
	    else
		value = Integer.parseInt(xmlSR.getAttributeValue(null,name));
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public float getFloatValue(XMLStreamReader xmlSR, String name)
    {
	float value=-1f;
	try{
	    if(xmlSR.getAttributeValue(null,name) == null)
		value= -1f;
	    else
		value=Float.parseFloat(xmlSR.getAttributeValue(null,name));
	}
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	return value;
    }

    public void getPeaks(String peakData, int count)
    {
	int precision = -1;
	if(count == 1)
	    {
		precision = tmpScanHeader.getMassPrecision();
	    }
	if(count == 2)
	    {
		precision = tmpScanHeader.getIntenPrecision();
	    }
	//support non-zlib
	byte[] peakArray = peakData.getBytes();
	byte[] outPeakArray = peakArray;
	int outpos = Base64.decode(peakArray,0,peakArray.length,outPeakArray);

      
	double[] doubleMassList = null;
       
	double[] doubleIntenList = null;

	ByteBuffer peakBuffer = null;
	//check if it's compressed
	byte[] result=null;
	int unCompLen = outpos;
	String compressType = "None";
	if(count == 1)
	    compressType = tmpScanHeader.getMassCompressionType();
	if(count == 2)
	    compressType = tmpScanHeader.getIntenCompressionType();
//System.err.println("***Compress type: " + compressType);
	if(compressType.equals("zlib"))
	    {
		try{
		    Inflater decompresser = new Inflater();
		    decompresser.setInput(outPeakArray, 0, outpos);
		    unCompLen = (tmpScanHeader.getPeaksCount())*(precision/4);
		    result = new byte[unCompLen];
		    decompresser.inflate(result);
		    decompresser.end();

		}
		catch(DataFormatException e)
		    {
			e.printStackTrace();
		    }
		
		peakBuffer = ByteBuffer.wrap(result);
		peakBuffer.order(ByteOrder.LITTLE_ENDIAN);
		
	    }
	else
	    {
		peakBuffer = ByteBuffer.wrap(outPeakArray,0,outpos);
		peakBuffer.order(ByteOrder.LITTLE_ENDIAN);
	    }
//System.err.println("Creating list, precision=" + precision + ", count=" + count);
	if(precision == 64)
	    {
		if(count == 1)
		    {
			doubleMassList = new double[unCompLen/8];
			 int i=0;
			 while(peakBuffer.hasRemaining())
			     {
				 doubleMassList[i] = peakBuffer.getDouble();
				 i++;
			     }
			 tmpScan.setDoubleMassList(doubleMassList);
			 //System.out.println("massList size "+tmpScan.getDoubleMassList().length);
			 
		    }
	
		if(count == 2)
		    {
			doubleIntenList = new double[unCompLen/8];
			 int i=0;
			 while(peakBuffer.hasRemaining())
			     {
				 doubleIntenList[i] = peakBuffer.getDouble();
				 i++;
			     }
			 tmpScan.setDoubleIntensityList(doubleIntenList);
		    }
	    }
	else
	    {
			
		if(count == 1)
		    {
			doubleMassList = new double[unCompLen/4];
			 int i=0;
			 while(peakBuffer.hasRemaining())
			     {
				 doubleMassList[i] = (double)peakBuffer.getFloat();
				 i++;
			     }
			 tmpScan.setDoubleMassList(doubleMassList);
		    }
		if(count == 2)
		    {
			doubleIntenList = new double[unCompLen/4];
			int i=0;
			 while(peakBuffer.hasRemaining())
			     {
				 doubleIntenList[i] = (double)peakBuffer.getFloat();
				 i++;
			     }
			 tmpScan.setDoubleIntensityList(doubleIntenList);
			 //System.out.println("intenList size "+tmpScan.getFloatIntensityList().length);
		    }
	    }
        //dhmay fixing up the massIntensityList, 2009/03/09.  This seems to have been missed initially
        if (count == 2)
        {
//System.err.println("****Setting mass-int list");
            double[][] massIntensityList = new double[2][];
            massIntensityList[0] = tmpScan.getDoubleMassList();
            massIntensityList[1] = tmpScan.getDoubleIntensityList();
            tmpScan.setMassIntensityList(massIntensityList);
        }
//System.err.println("getPeaks done, masslist length=" + (tmpScan.getDoubleMassList() == null ? "null" : tmpScan.getDoubleMassList().length) + ", inten=" + (tmpScan.getDoubleIntensityList() == null? "null" : tmpScan.getDoubleIntensityList().length));
    }
 
    
}
