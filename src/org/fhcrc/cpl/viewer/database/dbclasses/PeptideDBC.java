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
package org.fhcrc.cpl.viewer.database.dbclasses;

import org.apache.log4j.Logger;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class PeptideDBC
{
    protected static Logger _log = Logger.getLogger(PeptideDBC.class);

    protected int id;
    protected String sequence = null;
    protected Set<ProteinDBC> proteins = null;
    protected Set<PeptideHitDBC> peptideHits = null;


    public PeptideDBC()
    {
        proteins = new HashSet<ProteinDBC>();
    }

    public PeptideDBC(String sequence)
    {
        this();
        this.sequence = sequence;
    }

    public String toString()
    {
        return "Peptide: " + sequence;
    }

    public String getSequence()
    {
        return sequence;
    }

    public void setSequence(String sequence)
    {
        this.sequence = sequence;
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }

    public Set<ProteinDBC> getProteins()
    {
        return proteins;
    }

    public void setProteins(Set<ProteinDBC> proteins)
    {
        this.proteins = proteins;
    }

    public void addProtein(ProteinDBC protein)
    {
        if (proteins == null)
            proteins = new HashSet<ProteinDBC>();
        proteins.add(protein);
    }


    public Set<PeptideHitDBC> getPeptideHits()
    {
        return peptideHits;
    }

    public void setPeptideHits(Set<PeptideHitDBC> peptideHits)
    {
        this.peptideHits = peptideHits;
    }
}
