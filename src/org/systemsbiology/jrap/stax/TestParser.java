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

/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) TestParser.java * Author: * Robert M. Hubley
 * rhubley@systemsbiology.org
 * modified by Ning Zhang 
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
 * ******************************************************************************
 * 
 * ChangeLog
 * 
 * $Log: TestParser.java,v $ Revision 1.1.1.1 2003/04/09 00:02:54 ppatrick
 * Initial import.
 * 
 *  
 ******************************************************************************/
package org.systemsbiology.jrap.stax;



/**
 * A demo program which uses the MSXMLParser to read in scans from a file you
 * specify on the command line.
 *  
 */
public class TestParser
{

	public static void main(String[] asArg)
	{
		Scan myScan = null;
		ScanHeader myHeader = null;
		MZXMLFileInfo fileInfo = null;

		if (asArg.length > 1)
		{
		    MSXMLParser myParser = new MSXMLParser(asArg[0]);
			System.out.println("There are total "+myParser.getScanCount()+" scans.");
			
			System.out.println("=============================================");
			System.out.println("Writing out FileHeader info!");
			fileInfo = myParser.rapFileHeader();
			System.out.println(fileInfo);

			
			
			System.out.println("=============================================");
			System.out.println("Writing out ScanHeader info!");
			for (int i = 1; i < asArg.length; i++)
			{
				myHeader = myParser.rapHeader(Integer.parseInt(asArg[i]));
				System.out.println(myHeader);
			}

			System.out.println("==============================================");
			System.out.println("Writing out Scan info!");
			for (int i = 1; i < asArg.length; i++)
			{
			    myScan = myParser.rap(Integer.parseInt(asArg[i]));
			    System.out.println(myScan);
			}
			
		       
		} else
		{
			System.out.println(
				"Invalid number of arguments: TestParser fileName scanNumber1 <scanNumber2>..");
		}
	}

}
