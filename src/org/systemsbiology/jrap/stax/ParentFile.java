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
 * File: * @(#) ParentFile.java * Author: * Mathijs Vogelzang
 * m_v@dds.nl
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
 * 10-05-2004 Added this header
 * 
 * Created on May 21, 2004
 * 
 * on Jan. 22, 2008, changed the int type to a String type
 *  
 ******************************************************************************/
package org.systemsbiology.jrap.stax;

/**
 * The ParentFile class contains information about parent files
 * of an mzXML file.
 * 
 * @author Mathijs
 */
public class ParentFile
{
    //public final static int TYPE_RAW = 1, TYPE_PROCESSED = 2;
	
    protected String URI, sha1,type;
       
	
	public ParentFile(String URI, String type, String sha1)
	{
		this.URI = URI;
		this.sha1 = sha1;
		this.type = type;
	}
	
	/**
	 * Get the URI of this file.
	 */
	public String getURI()
	{
		return URI;
	}
	
	/**
	 * Get the sha1-sum of this file.
	 */
	public String getSha1()
	{
		return sha1;
	}
	
	/**
	 * Return the type of parent file. 
	 * 
	 * This value is either TYPE_RAW or TYPE_PROCESSED.
	 * 
	 * @return the type of parent file.
	 */
	public String getType()
	{
		return type;
	}

    public String toString()
    {
	return ("URI "+URI+" sha1 "+sha1+" type "+type);
    }
	
}
