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
package org.fhcrc.cpl.viewer.database.loader;

import org.fhcrc.cpl.viewer.database.HibernateManager;
import org.fhcrc.cpl.viewer.database.dbclasses.ProteinDBC;
import org.fhcrc.cpl.viewer.database.dbclasses.GeneDBC;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import org.hibernate.*;

import java.util.*;
import java.io.*;

/**
 * Utility class for loading fasta files
 */
public class ProteinGeneMapDBLoader
{
    protected static Logger _log = Logger.getLogger(ProteinGeneMapDBLoader.class);

    public ProteinGeneMapDBLoader()
    {
    }

    public void load(File proteinGeneMapFile)
            throws IOException
    {
        Session session = HibernateManager.getInstance().openSession();

        session.beginTransaction();

        Set<Object> objectsToSave = new HashSet<Object>();
        Set<Object> objectsToUpdate = new HashSet<Object>();


        Map<String,GeneDBC> geneSymbolGeneMap = new HashMap<String,GeneDBC>();

        Map<String, String[]> ipiGeneSymbolMap =
                readIPIGeneMap(proteinGeneMapFile);

        Map<String, GeneDBC> dbSymbolGeneMap =
                new HashMap<String, GeneDBC>();
        for (String[] geneSymbols : ipiGeneSymbolMap.values())
        {
            for (String geneSymbol : geneSymbols)
            {
                GeneDBC geneDBC = geneSymbolGeneMap.get(geneSymbol);
                if (geneDBC == null)
                {
                    geneDBC = new GeneDBC(geneSymbol, null);
                    geneSymbolGeneMap.put(geneSymbol, geneDBC);
                }
                if (!dbSymbolGeneMap.containsKey(geneSymbol))
                    dbSymbolGeneMap.put(geneSymbol, geneDBC);
                objectsToSave.add(geneDBC);
            }
        }

        ApplicationContext.infoMessage("About to query proteins from database");

        List<ProteinDBC> proteinsFromDB =
                session.createQuery("from ProteinDBC").list();
        ApplicationContext.infoMessage("Queried " + proteinsFromDB.size() + " proteins");
        
        Map<String,ProteinDBC> proteinsByIPI = new HashMap<String,ProteinDBC>();
        for (ProteinDBC protein : proteinsFromDB)
        {
            proteinsByIPI.put(protein.getLookup(), protein);
        }

        ApplicationContext.infoMessage("Created IPI protein map");

        for (String ipi : ipiGeneSymbolMap.keySet())
        {
            ProteinDBC protein = proteinsByIPI.get(ipi);
            if (protein != null)
            {
                for (String geneSymbol : ipiGeneSymbolMap.get(ipi))
                    protein.addGene(dbSymbolGeneMap.get(geneSymbol));
                objectsToUpdate.add(protein);
            }
        }

        ApplicationContext.infoMessage("About to save " + objectsToSave.size() +
                " objects to database");

        session = HibernateManager.getInstance().openSession();
        Transaction transaction = session.beginTransaction();
        for (Object objectToSave : objectsToSave)
        {
            session.save(objectToSave);
        }
        for (Object objectToUpdate : objectsToUpdate)
        {
            session.update(objectToUpdate);
        }
//        session.flush();

        transaction.commit();
        session.close();

        ApplicationContext.infoMessage("Done with DB commit");
    }


    protected Map<String,String[]> readIPIGeneMap(File proteinGeneMapFile)
            throws IOException
    {
        Map<String,String[]> ipiGeneMap = new HashMap<String,String[]>();

        BufferedReader br = new BufferedReader(new FileReader(proteinGeneMapFile));
        String line;
        while ((line = br.readLine()) != null)
        {
            String[] words = line.split("\t");
            String ipi = words[0];
            if (words.length < 2 || words[1] == null || words[1].length() < 1)
                continue;
            String geneListString = words[1];

            String[] geneArray = geneListString.split("//");

            ipiGeneMap.put(ipi,geneArray);
        }

        ApplicationContext.infoMessage("Loaded " + ipiGeneMap.size() + " genes from lookup file");

        return ipiGeneMap;
    }
}
