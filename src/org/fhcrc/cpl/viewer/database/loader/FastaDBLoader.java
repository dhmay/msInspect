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
import org.fhcrc.cpl.viewer.database.dbclasses.BioSequenceDBC;
import org.fhcrc.cpl.toolbox.proteomics.ProteinUtilities;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.apache.log4j.Logger;

import org.hibernate.*;
import org.fhcrc.cpl.toolbox.proteomics.Protein;

import java.util.*;
import java.io.File;

/**
 * Utility class for loading fasta files
 */
public class FastaDBLoader
{
    protected static Logger _log = Logger.getLogger(FastaDBLoader.class);

    protected boolean loadProteinSequences = false;

    public FastaDBLoader()
    {
    }

    public void load(File fastaFile)
    {
        Protein[] proteinsInFasta =
                ProteinUtilities.loadProteinsFromFasta(fastaFile).toArray(new Protein[0]);
        load(proteinsInFasta);
    }

    public void load(Protein[] proteinsInFasta)
    {
        Session session = HibernateManager.getInstance().openSession();

        session.beginTransaction();

        Set<Object> objectsToSave = new HashSet<Object>();

        Map<String, ProteinDBC> dbIPIProteinMap = new HashMap<String, ProteinDBC>();
        for (Protein protein : proteinsInFasta)
        {
            ProteinDBC proteinDBC = new ProteinDBC(protein);
            dbIPIProteinMap.put(protein.getLookup(), proteinDBC);
            objectsToSave.add(proteinDBC);
        }

        if (loadProteinSequences)
        {
            for (Protein protein : proteinsInFasta)
            {
                BioSequenceDBC bioSequence = new BioSequenceDBC();
                bioSequence.setSequence(protein.getSequenceAsString());
                ProteinDBC dbProtein = dbIPIProteinMap.get(protein.getLookup());
                dbProtein.setBioSequence(bioSequence);

                objectsToSave.add(bioSequence);
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
        session.flush();

        transaction.commit();

        ApplicationContext.infoMessage("Done with DB commit");
    }


    public boolean isLoadProteinSequences()
    {
        return loadProteinSequences;
    }

    public void setLoadProteinSequences(boolean loadProteinSequences)
    {
        this.loadProteinSequences = loadProteinSequences;
    }
}
