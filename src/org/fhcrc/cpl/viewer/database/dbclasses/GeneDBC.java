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
public class GeneDBC
{
    protected static Logger _log = Logger.getLogger(GeneDBC.class);

    protected int id;
    protected String symbol = null;
    protected String name = null;
    protected Set<ProteinDBC> proteins = null;

    public GeneDBC()
    {
        proteins = new HashSet<ProteinDBC>();
    }

    public GeneDBC(String symbol, String name)
    {
        this();
        this.symbol = symbol;
        this.name = name;
    }

    public String toString()
    {
        return symbol + "\t" + name;
    }

    public String getSymbol()
    {
        return symbol;
    }

    public void setSymbol(String symbol)
    {
        this.symbol = symbol;
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
        if (!protein.getGenes().contains(this))
            protein.addGene(this);
    }
}
