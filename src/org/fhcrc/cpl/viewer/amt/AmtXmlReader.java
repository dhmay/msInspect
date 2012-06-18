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
package org.fhcrc.cpl.viewer.amt;

import java.io.*;
import java.util.List;
import java.util.ArrayList;

import org.w3c.dom.Node;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.filehandler.Stax2DomBuilder;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.proteomics.schemaRevision.amtXml10.*;

import javax.xml.stream.XMLStreamException;

/**
* A restrictive wrapper for reading AmtXml files.  The method is a compromise between
 * the elegance of DOM/XmlBeans and the efficiency of stax.  We use Stax2DomBuilder to pull DOM subtrees
 * out of the pepxml file for bits that we care about, namely modifications and features
 */
public class AmtXmlReader
{
    static Logger _log = Logger.getLogger(AmtXmlReader.class);

    //stores the database read from a file
    protected AmtDatabase mAmtDatabase = null;

    /**
     * Read in a amtxml file
     * @param file
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public AmtXmlReader(File file) throws FileNotFoundException, XMLStreamException
    {
        read(file);
    }

    /** 
     * Read in a file, extract modifications and features
     * @param file
     */
    public void read(File file) throws FileNotFoundException, XMLStreamException
    {
        //this is not the most efficient way to do it.  If performance is bad, reimplement.
        //On the other hand, it's simple, and it doesn't waste any memory to speak of
        _log.debug("Reading metadata");
        extractMetadataAndModifications(file);
        _log.debug("Extracting runs");        
        extractRuns(file);
        _log.debug("Extracting entries");
        extractEntries(file);
    }

    /**
     * TODO: get hydrophobicity calculator and version
     * @param inputFile
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    protected void extractMetadataAndModifications(File inputFile)
            throws FileNotFoundException, XMLStreamException
    {
        if (mAmtDatabase == null)
            mAmtDatabase = new AmtDatabase();
        mAmtDatabase.setAmtDBSourceFile(inputFile);

        Stax2DomBuilder builder = new Stax2DomBuilder(inputFile);

        //these XmlBeans objects are necessary throway structure to get XmlBeans to parse the xml
        //we pull out with Stax2DomBuilder
        AmtDatabaseDocument xmlBeansAmtDatabaseDocument =
                AmtDatabaseDocument.Factory.newInstance();
        AmtDatabaseDocument.AmtDatabase xmlBeansAmtDatabase =
                xmlBeansAmtDatabaseDocument.addNewAmtDatabase();
        //extract top-level stuff

        mAmtDatabase.setHydrophobicityAlgorithmName(xmlBeansAmtDatabase.getHydrophobicityAlgorithm());
        mAmtDatabase.setHydrophobicityAlgorithmVersion(xmlBeansAmtDatabase.getHydrophobicityAlgVersion());


        while (true)
        {
            Node treeRoot = builder.findTreeForName(xmlBeansAmtDatabase.getDomNode().getOwnerDocument(), "aminoacid_modification","amt_database");
            if (treeRoot == null)
                break;
            xmlBeansAmtDatabase.getDomNode().appendChild(treeRoot);
            AmtDatabaseDocument.AmtDatabase.AminoacidModification xmlBeansModification =
                    xmlBeansAmtDatabase.getAminoacidModificationArray(0);

            MS2Modification dbMod = new MS2Modification();
            dbMod.setAminoAcid(xmlBeansModification.getResidue());
            dbMod.setVariable(xmlBeansModification.getVariableFlag());
            dbMod.setMassDiff(xmlBeansModification.getMassDifference().floatValue());
            mAmtDatabase.addAminoacidModification(dbMod);
            //remove the dummy node
            xmlBeansAmtDatabase.getDomNode().removeChild(treeRoot);            
        }
    }

    /**
     * extracts interesting stuff at the amt_database level, and all runs
     * @param inputFile
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    protected void extractRuns(File inputFile) throws FileNotFoundException, XMLStreamException
    {
        Stax2DomBuilder builder = new Stax2DomBuilder(inputFile);

        //these XmlBeans objects are necessary throway structure to get XmlBeans to parse the xml
        //we pull out with Stax2DomBuilder
        AmtDatabaseDocument xmlBeansAmtDatabaseDocument =
                AmtDatabaseDocument.Factory.newInstance();
        AmtDatabaseDocument.AmtDatabase xmlBeansAmtDatabase =
                xmlBeansAmtDatabaseDocument.addNewAmtDatabase();

        while (true)
        {
            Node treeRoot = builder.findTreeForName(xmlBeansAmtDatabase.getDomNode().getOwnerDocument(), "run","amt_database");
            if (treeRoot == null)
                break;
            xmlBeansAmtDatabase.getDomNode().appendChild(treeRoot);
            AmtDatabaseDocument.AmtDatabase.Run xmlBeansRun =
                    xmlBeansAmtDatabase.getRunArray(0);
            List<MS2Modification> ms2ModList = new ArrayList<MS2Modification>();

            for (AmtDatabaseDocument.AmtDatabase.Run.AminoacidModUsage modUsage :
                    xmlBeansRun.getAminoacidModUsageArray())
            {
                ms2ModList.add(mAmtDatabase.getAminoacidModificationBySequence(modUsage.getModificationId()));
            }

            AmtDatabaseDocument.AmtDatabase.Run.TimeHydroMappingCoeff[] mappingCoeffs
                    = xmlBeansRun.getTimeHydroMappingCoeffArray();

            double[] coeffs = new double[mappingCoeffs.length];
            for (int i=0; i<coeffs.length; i++)
            {
                coeffs[mappingCoeffs[i].getDegree()] = mappingCoeffs[i].getValue().doubleValue();
            }
            AmtRunEntry runEntry =
                    new AmtRunEntry(coeffs,
                            ms2ModList.toArray(new MS2Modification[0]),
                            xmlBeansRun.getTimeAdded().getTime());

            mAmtDatabase.addRunEntry(runEntry);

            //check that the sequence worked out ok
            assert(mAmtDatabase.getSequenceForRun(runEntry) == xmlBeansRun.getRunId());

            if (xmlBeansRun.getMzxmlFilename() != null &&
                !("".equals(xmlBeansRun.getMzxmlFilename())))
                runEntry.setMzXmlFilename(xmlBeansRun.getMzxmlFilename());
            if (xmlBeansRun.getPepxmlFilename() != null &&
                !("".equals(xmlBeansRun.getPepxmlFilename())))
                runEntry.setPepXmlFilename(xmlBeansRun.getPepxmlFilename());
            if (xmlBeansRun.getLSID() != null &&
                !("".equals(xmlBeansRun.getLSID())))
                runEntry.setLSID(xmlBeansRun.getLSID());
            if (xmlBeansRun.getMinPeptideProphet() != null)
                runEntry.setMinPeptideProphet(xmlBeansRun.getMinPeptideProphet().doubleValue());
            if (xmlBeansRun.getTimeAnalyzed() != null)
                runEntry.setTimeAnalyzed(xmlBeansRun.getTimeAnalyzed().getTime());

            //remove the dummy node
            xmlBeansAmtDatabase.getDomNode().removeChild(treeRoot);
        }
    }


    /**
     * @param inputFile
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    protected void extractEntries(File inputFile) throws FileNotFoundException, XMLStreamException
    {
        if (mAmtDatabase == null)
            mAmtDatabase = new AmtDatabase();

        Stax2DomBuilder builder = new Stax2DomBuilder(inputFile);

        //these XmlBeans objects are necessary throway structure to get XmlBeans to parse the xml
        //we pull out with Stax2DomBuilder
        AmtDatabaseDocument xmlBeansAmtDatabaseDocument = AmtDatabaseDocument.Factory.newInstance();
        AmtDatabaseDocument.AmtDatabase xmlBeansAmtDatabase = xmlBeansAmtDatabaseDocument.addNewAmtDatabase();

        //Because we have potentially a Very Large Number of entries, it's necessary to add them
        //and remove them one-by-one, so as not to fill up the available memory.
        while (true)
        {
            Node treeRoot = builder.findTreeForName(xmlBeansAmtDatabase.getDomNode().getOwnerDocument(), "peptide_entry","amt_database");
            if (treeRoot == null)
                break;
            xmlBeansAmtDatabase.getDomNode().appendChild(treeRoot);
            AmtDatabaseDocument.AmtDatabase.PeptideEntry xmlBeansPeptideEntry =
                    xmlBeansAmtDatabase.getPeptideEntryArray(0);

            AmtPeptideEntry entry = new AmtPeptideEntry();
            entry.setPeptideSequence(xmlBeansPeptideEntry.getPeptideSequence());
//            entry.setMedianObservedHydrophobicity(xmlBeansPeptideEntry.getMedianObservedHydrophobicity().doubleValue());
            _log.debug("extractEntries: peptide " + entry.getPeptideSequence());

            entry.setPredictedHydrophobicity(xmlBeansPeptideEntry.getCalculatedHydrophobicity().doubleValue());
            for (AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry xmlBeansModState :
                    xmlBeansPeptideEntry.getModificationStateEntryArray())
            {
                _log.debug("extractEntries: mod state");

                List<MS2Modification>[] ms2Modifications = null;
                AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry.AminoacidModInstance[]
                    xmlBeansModifications = xmlBeansModState.getAminoacidModInstanceArray();
                if (xmlBeansModifications != null && xmlBeansModifications.length > 0)
                {
                    ms2Modifications =
                            (List<MS2Modification>[]) new List[xmlBeansPeptideEntry.getPeptideSequence().length()];

                    //remember, position = index + 1
                    for (AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry.AminoacidModInstance xmlBeansModification :
                            xmlBeansModifications)
                    {
                        if (ms2Modifications[xmlBeansModification.getPosition()] == null)
                            ms2Modifications[xmlBeansModification.getPosition()] =
                                    new ArrayList<MS2Modification>();
                        ms2Modifications[xmlBeansModification.getPosition()].add(
                                mAmtDatabase.getAminoacidModificationBySequence(xmlBeansModification.getModificationId()));
                    }
                }
                _log.debug("extractEntries: mod state 2");

                AmtPeptideEntry.AmtPeptideModificationStateEntry modState =
                        entry.addModificationStateEntry(
                                xmlBeansModState.getModifiedSequence(),
                                xmlBeansModState.getModifiedMass().doubleValue(),
                                ms2Modifications);
                _log.debug("extractEntries: mod state 3");

                for (AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry.Observation xmlBeansObservation :
                        xmlBeansModState.getObservationArray())
                {
                    AmtPeptideEntry.AmtPeptideObservation observation =
                            AmtPeptideEntry.AmtPeptideObservation.createObservation(
                                xmlBeansObservation.getObservedHydrophobicity().doubleValue(),
                                xmlBeansObservation.getPeptideProphet().doubleValue(),
                                mAmtDatabase.getRunBySequence(xmlBeansObservation.getRunId()),
                                xmlBeansObservation.getTimeInRun().doubleValue());
                    if (xmlBeansObservation.isSetSpectralCount())
                        observation.setSpectralCount(xmlBeansObservation.getSpectralCount());
                    modState.addObservationNoRecalc(observation);
                }

                _log.debug("extractEntries: mod state 4");

                modState.recalculateStats();
            }
            entry.recalculateStats();
            mAmtDatabase.addObservationsFromEntry(entry);
            xmlBeansAmtDatabase.removePeptideEntry(0);

        }
    }

    /**
     * Accessor for loaded database
     * @return
     */
    public AmtDatabase getDatabase()
    {
        return mAmtDatabase;
    }
}
