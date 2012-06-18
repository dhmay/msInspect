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
public class PeptideHitDBC
{
    protected static Logger _log = Logger.getLogger(PeptideHitDBC.class);

    protected int id;
    protected PeptideDBC peptide = null;
    protected float peptideProphet;
    protected Set<ProteinHitDBC> proteinHits = null;
    protected ExperimentRunDBC experimentRunDBC = null;

    public PeptideHitDBC()
    {
        proteinHits = new HashSet<ProteinHitDBC>();
    }

    public String toString()
    {
        return "Peptide Hit: " + peptide + ", pprophet=" + peptideProphet;
    }

    public String getSequence()
    {
        return peptide.getSequence();
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }

    public Set<ProteinHitDBC> getProteinHits()
    {
        return proteinHits;
    }

    public void setProteinHits(Set<ProteinHitDBC> proteinHits)
    {
        this.proteinHits = proteinHits;
    }

    public void addProteinHit(ProteinHitDBC proteinHit)
    {
        if (proteinHits == null)
            proteinHits = new HashSet<ProteinHitDBC>();
        proteinHits.add(proteinHit);
    }

    public float getPeptideProphet()
    {
        return peptideProphet;
    }

    public void setPeptideProphet(float peptideProphet)
    {
        this.peptideProphet = peptideProphet;
    }


    public ExperimentRunDBC getExperimentRun()
    {
        return experimentRunDBC;
    }

    public void setExperimentRun(ExperimentRunDBC experimentRunDBC)
    {
        this.experimentRunDBC = experimentRunDBC;
    }


    public PeptideDBC getPeptide()
    {
        return peptide;
    }

    public void setPeptide(PeptideDBC peptide)
    {
        this.peptide = peptide;
    }
}
