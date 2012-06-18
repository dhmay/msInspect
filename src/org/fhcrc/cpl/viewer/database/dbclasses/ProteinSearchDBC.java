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
import java.util.Date;


/**
 * Models a Protein Search
 */
public class ProteinSearchDBC
{
    protected static Logger _log = Logger.getLogger(ProteinSearchDBC.class);

    protected int id;
    protected String name;
    protected Date date;
    protected Set<ProteinGroupDBC> proteinGroups = null;
    protected Set<ExperimentRunDBC> experimentRuns = null;


    public ProteinSearchDBC()
    {
        proteinGroups = new HashSet<ProteinGroupDBC>();
    }

    public String toString()
    {
        return "Protein Search";
    }

    public Set<ProteinGroupDBC> getProteinGroups()
    {
        return proteinGroups;
    }

    public void setProteinGroups(Set<ProteinGroupDBC> proteinGroups)
    {
        this.proteinGroups = proteinGroups;
    }

    public Set<ExperimentRunDBC> getExperimentRuns()
    {
        return experimentRuns;
    }

    public void setExperimentRuns(Set<ExperimentRunDBC> experimentRuns)
    {
        this.experimentRuns = experimentRuns;
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }


    public Date getDate()
    {
        return date;
    }

    public void setDate(Date date)
    {
        this.date = date;
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
