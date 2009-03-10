/***************************************************************************** 
* --------------------------------------------------------------------------- 
 * File: * @(#) IndexParser.java * Author: * Ning Zhang
 * nzhang@systemsbiology.org
 * ****************************************************************************** * * *
 * This software is provided ``AS IS'' and any express or implied 
 * warranties, including, but not limited to, the implied warranties of 
 * merchantability and fitness for a particular purpose, are disclaimed.  In
 * no event shall the authors or the Institute for Systems Biology  liable
 * for any direct, indirect, incidental, special, exemplary, or
 * consequential damages (including, but not limited to, procurement of
 * substitute goods or services; loss of use, data, or profits; or business
 * interruption) however caused and on any theory of liability, whether in
 * contract, strict liability, or tort (including negligence or otherwise)
 * arising in any way out of the use of this software, even if advised of
 * the possibility of such damage. 
 * ******************************************************************************/

package org.systemsbiology.jrap.stax;

import javax.xml.stream.*;
import javax.xml.stream.events.*;
import java.io.*;
import java.util.*;

public class IndexParser{

    String inputMZXMLfile;
    Map<Integer, Long> offsetMap = new HashMap<Integer, Long>(10000);
    
    //for mzML
    long chrogramIndex=-1;

    int maxScan = -1;
    int currentScan = -1;

    //Is the file mzXML or mzML?
    boolean isXML = false;
    boolean isML = false;

    public IndexParser(String inputMZXMLfile)
    {
	this.inputMZXMLfile = inputMZXMLfile;
	
	//determine whether the file is mzXML or mzML?
	if(inputMZXMLfile.indexOf("mzXML") != -1)
	    isXML = true;
	else
	    isML = true;
    }

    public Map<Integer, Long> getOffsetMap()
    {
	System.out.println("offset size "+offsetMap.size());
	return offsetMap;
	
    }

    public long getChrogramIndex()
    {
	System.out.println("chrogramIndex "+chrogramIndex);
	return chrogramIndex;
	
    }

    public int getMaxScan()
    {
	System.out.println("maxScan "+maxScan);
        return maxScan;
    }

    private long getIndexPosition()
    {
	FileInputStream fileIN = null;
	long indexPosition = -1;
	File tmpXML = null;

	//mzXML and mzML have different element name
	String indexName = null;
	if(isXML)
	    indexName = "indexOffset";
	else
	    indexName = "indexListOffset";

	try
	    {
		tmpXML = new File(inputMZXMLfile);
		//System.out.println(inputMZXMLfile +" length is "+ tmpXML.length());

		fileIN = new FileInputStream(inputMZXMLfile);
		fileIN.skip(tmpXML.length() - 500);
		byte[] bytes = new byte[500];
		int bytesRead = fileIN.read(bytes);
		String footer = new String(bytes, 0, bytesRead);
		int offset;
		if ((offset = footer.indexOf("<"+indexName+">")) == -1)
		    {
			System.err.println("<"+indexName+">"+" not found!!!");
		    }
		footer = footer.substring(offset + indexName.length()+2);
		int endIndex = footer.indexOf("</"+indexName+">");
		if (endIndex == -1)
		    {
			System.err.println("</"+indexName+"> not found!!!");
		    }

		footer = footer.substring(0, endIndex);
		indexPosition = Long.parseLong(footer);
	    
		fileIN.close();
		System.out.println("indexPosition is "+indexPosition);
	    } catch (Exception e)
	    {
		System.out.println("exception:" + e);
		e.printStackTrace();
	    }
	
	return indexPosition;
    }

    public void parseIndexes() 
    {
	try{
	    long indexPos = getIndexPosition();
	    //System.out.println("indexPos "+indexPos);

	    FileInputStream fileIN = null;
	    fileIN = new FileInputStream(inputMZXMLfile);
	    fileIN.skip(indexPos);
	    
	    XMLInputFactory inputFactory = XMLInputFactory.newInstance();
	    XMLStreamReader xmlSR = inputFactory.createXMLStreamReader(fileIN);

	    String elementName = null;
	    StringBuffer indexBuffer = new StringBuffer();

	    boolean inOffSet = false;

	    //for mzML
	    boolean inSpec = false;
	    boolean inChrogram = false;

	    while(xmlSR.hasNext())
		{
		    int event = xmlSR.next();
		    if(event == xmlSR.START_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    //System.out.println("elementName "+elementName);

			    if(elementName.equals("index"))
				{
				    if(isML)
					{
					    if((xmlSR.getAttributeValue(null,"name")).equals("spectrum"))
						inSpec = true;
					    if((xmlSR.getAttributeValue(null,"name")).equals("chromatogram"))
						{
						    inSpec = false;
						    inChrogram = true;
						}
					}
				}
	 
			    if(elementName.equals("offset"))
				{
				    if(indexBuffer.length()>0)
					indexBuffer.delete(0,indexBuffer.capacity());
				    inOffSet = true;
				    if(isXML)
					currentScan = Integer.parseInt(xmlSR.getAttributeValue(null, "id"));
				    else if(inSpec)
					currentScan = Integer.parseInt(xmlSR.getAttributeValue(null, "nativeID"));
				}
			}
		    if(event == xmlSR.CHARACTERS)
			{
			    if(inOffSet)
				indexBuffer.append(xmlSR.getText());
			}
		    if(event == xmlSR.END_ELEMENT)
			{
			    elementName = xmlSR.getLocalName();
			    if(elementName.equals("offset"))
				{
				    long offset = Long.valueOf(indexBuffer.toString());
				    if(inChrogram)
					chrogramIndex = offset;
				    else
					{
					    assert currentScan != -1 : "Did not find the scan number associated with offset " + indexBuffer.toString();
					    //System.out.println("index "+indexBuffer.toString());
					    offsetMap.put(currentScan, offset);
					    maxScan = currentScan;
					    currentScan = -1;
					    indexBuffer.delete(0,indexBuffer.capacity());
					    inOffSet = false;
					}
				}
			    if(elementName.equals("index"))
				{
				    
				    if(isXML)
					throw new XMLStreamException("IndexEndFoundException");
				    else if(inChrogram)
					throw new XMLStreamException("IndexEndFoundException");

				    
				    
				}
			}
		}
	    
	}
	catch(Exception e)
	    {
		if(!(e.getMessage()).equals("IndexEndFoundException"))
		    e.printStackTrace(System.err);
	    }
    }

    
}
