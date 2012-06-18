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
 * File: * @(#) SoftwareInfo.java * Author: * Mathijs Vogelzang
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
 ******************************************************************************/
package org.systemsbiology.jrap.stax;

/**
 * SoftwareInfo represents data that is available on software that has processed
 * a mzXML file.
 * 
 * @author Mathijs Vogelzang
 */
public class SoftwareInfo
{
	public String type, name, version;

    public SoftwareInfo(String type, String name, String version)
    {
	this.type = type;
	this.name = name;
	this.version = version;
    }

    public void setType(String type)
    {
	this.type = type;
    }

    public void setName(String name)
    {
	this.name = name;
    }
    public void setVersion(String version)
    {
	this.version = version;
    }

    public String toString()
    {
	return ("type "+type+" name "+name+" version "+version);
    }
}
