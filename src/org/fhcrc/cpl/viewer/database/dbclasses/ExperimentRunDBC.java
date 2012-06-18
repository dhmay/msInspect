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

import java.util.*;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class ExperimentRunDBC
{
    protected static Logger _log = Logger.getLogger(ExperimentRunDBC.class);

    protected int id;
    protected String name = null;
    protected Date date = null;

    protected Set<ProteinSearchDBC> proteinSearches = null;
    protected Set<PeptideHitDBC> peptideHits = null;


    public ExperimentRunDBC()
    {
        proteinSearches = new HashSet<ProteinSearchDBC>();
    }

    public ExperimentRunDBC(String name)
    {
        this();
        this.name = name;
    }

    public String toString()
    {
        return "ExperimentRun: " + name;
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }

    public Date getDate()
    {
        return date;
    }

    public void setDate(Date date)
    {
        this.date = date;
    }

    public Set<ProteinSearchDBC> getProteinSearches()
    {
        return proteinSearches;
    }

    public void setProteinSearches(Set<ProteinSearchDBC> proteinSearches)
    {
        this.proteinSearches = proteinSearches;
    }

    public void addProteinSearch(ProteinSearchDBC proteinSearch)
    {
        proteinSearches.add(proteinSearch);
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
