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

import org.fhcrc.cpl.toolbox.proteomics.Protein;

import java.io.File;
import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.text.Format;
import java.text.DecimalFormat;

/**
 * User: adam
 * Date: Jan 12, 2008
 * Time: 7:43:36 PM
 */
public class FastaValidator
{
    private Set<String> _proteinNames = new HashSet<String>(1000);
    private List<String> _errors = new ArrayList<String>();
    private File _fastaFile;


    public FastaValidator(File fastaFile)
    {
        _fastaFile = fastaFile;
    }


    public List<String> getErrors()
    {
        return _errors;
    }


    public List<String> validate()
    {
        Format lineFormat = DecimalFormat.getIntegerInstance();
        FastaLoader curLoader = new FastaLoader(_fastaFile);

        for (FastaLoader.ProteinIterator proteinIterator = curLoader.iterator(); proteinIterator.hasNext();)
        {
            Protein protein = proteinIterator.next();
            String lookup = protein.getLookup();

            if (_proteinNames.contains(lookup))
            {
                _errors.add("Line " + lineFormat.format(proteinIterator.getLastHeaderLine()) + ": " + lookup + " is a duplicate protein name");
            }
            else
            {
                _proteinNames.add(lookup);
            }
        }

        return _errors;
    }
}
