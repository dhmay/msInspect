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
public class ProteinHitDBC
{
    protected static Logger _log = Logger.getLogger(ProteinHitDBC.class);

    protected int id;
    protected ProteinDBC protein = null;
    protected float proteinProphet;
    protected Set<PeptideHitDBC> peptideHits = null;
    protected ProteinSearchDBC proteinSearch = null;

    public ProteinHitDBC()
    {
        peptideHits = new HashSet<PeptideHitDBC>();
    }

    public String toString()
    {
        return "Peptide Hit: " + protein + ", pprophet=" + proteinProphet;
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }


    public ProteinDBC getProtein()
    {
        return protein;
    }

    public void setProtein(ProteinDBC protein)
    {
        this.protein = protein;
    }

    public Set<PeptideHitDBC> getPeptideHits()
    {
        return peptideHits;
    }

    public void setPeptideHits(Set<PeptideHitDBC> peptideHits)
    {
        this.peptideHits = peptideHits;
    }

    public void addPeptideHit(PeptideHitDBC peptideHit)
    {
        peptideHits.add(peptideHit);
    }

    public float getProteinProphet()
    {
        return proteinProphet;
    }

    public void setProteinProphet(float proteinProphet)
    {
        this.proteinProphet = proteinProphet;
    }

    public ProteinSearchDBC getProteinSearch()
    {
        return proteinSearch;
    }

    public void setProteinSearch(ProteinSearchDBC proteinSearch)
    {
        this.proteinSearch = proteinSearch;
    }
}
