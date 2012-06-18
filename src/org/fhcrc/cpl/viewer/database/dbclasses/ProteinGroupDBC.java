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
import java.util.Set;
import java.util.HashSet;


/**
 * Models a Protein Group
 */
public class ProteinGroupDBC
{
    protected static Logger _log = Logger.getLogger(ProteinGroupDBC.class);

    protected int id;
    protected int number;
    protected ExperimentRunDBC experimentRun;
    protected ProteinSearchDBC proteinSearch;
    protected Set<ProteinDBC> proteins = null;
    protected Set<PeptideDBC> peptides = null;

    public ProteinGroupDBC()
    {
        proteins = new HashSet<ProteinDBC>();
        peptides = new HashSet<PeptideDBC>();
    }

    public String toString()
    {
        return "Protein Group: " + number;
    }

    public int getNumber()
    {
        return number;
    }

    public void setNumber(int number)
    {
        this.number = number;
    }

    public ExperimentRunDBC getExperimentRun()
    {
        return experimentRun;
    }

    public void setExperimentRun(ExperimentRunDBC experimentRun)
    {
        this.experimentRun = experimentRun;
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
        proteins.add(protein);
    }

    public Set<PeptideDBC> getPeptides()
    {
        return peptides;
    }

    public void setPeptides(Set<PeptideDBC> peptides)
    {
        this.peptides = peptides;
    }

    public void addPeptide(PeptideDBC peptide)
    {
        peptides.add(peptide);
    }

    public ProteinSearchDBC getProteinSearch()
    {
        return proteinSearch;
    }

    public void setProteinSearch(ProteinSearchDBC proteinSearch)
    {
        this.proteinSearch = proteinSearch;
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }      
}
