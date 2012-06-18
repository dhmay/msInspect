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
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;

import java.util.Set;
import java.util.HashSet;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class ProteinDBC
{
    protected static Logger _log = Logger.getLogger(ProteinDBC.class);

    protected int id;
    protected Protein protein = null;
    protected Set<GeneDBC> genes = null;
    protected Set<ProteinGroupDBC> proteinGroups = null;
    protected Set<PeptideDBC> peptides = null;
    protected BioSequenceDBC bioSequence = null;


    public ProteinDBC()
    {
        proteinGroups = new HashSet<ProteinGroupDBC>();
        peptides = new HashSet<PeptideDBC>();
        genes = new HashSet<GeneDBC>();

        //interim protein until one is set
        protein = new Protein("","".getBytes());
    }

    public ProteinDBC(Protein toolsProtein)
    {
        this();
        protein = toolsProtein;
    }

    public ProteinDBC(String header)
    {
        this(new Protein(header, new byte[0]));
    }

    /**
     * TODO:  is this correct?
     * @param xmlProtein
     */
    public ProteinDBC(ProtXmlReader.Protein xmlProtein)
    {
        this(new Protein(xmlProtein.getProteinName(), null));
    }

    public String toString()
    {
        if (protein == null)
            return "Unknown Protein";
        return "Protein: " + protein.getLookup();
    }

    public String getHeader()
    {
        String proteinHeader = protein.getHeader();
        if ("".equals(proteinHeader))
            return null;
        return proteinHeader;
    }

    public void setHeader(String header)
    {
        protein.setHeader(header);
    }

    public String getLookup()
    {
        String proteinLookup = protein.getLookup();
        if ("".equals(proteinLookup))
            return null;
        return proteinLookup;
    }

    public void setLookup(String lookup)
    {
        protein.setLookup(lookup);
    }

    public int getId()
    {
        return id;
    }

    public void setId(int id)
    {
        this.id = id;
    }

    public Set<GeneDBC> getGenes()
    {
        return genes;
    }

    public void setGenes(Set<GeneDBC> genes)
    {
        this.genes = genes;
    }

    public void addGene(GeneDBC gene)
    {
        genes.add(gene);
    }

    public void setProteinGroups(Set<ProteinGroupDBC> proteinGroups)
    {
        this.proteinGroups = proteinGroups;
    }

    public void addProteinGroup(ProteinGroupDBC proteinGroup)
    {
        proteinGroups.add(proteinGroup);
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

    public Set<ProteinGroupDBC> getProteinGroups()
    {
        return proteinGroups;
    }

    public Protein getProtein()
    {
        return protein;
    }

    public void setProtein(Protein protein)
    {
        this.protein = protein;
    }

    public BioSequenceDBC getBioSequence()
    {
        return bioSequence;
    }

    public void setBioSequence(BioSequenceDBC bioSequence)
    {
        this.bioSequence = bioSequence;
    }
}
