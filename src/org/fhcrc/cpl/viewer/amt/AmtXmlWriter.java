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
import java.util.Calendar;
import java.util.List;
import java.math.BigDecimal;

import org.apache.log4j.Logger;

import org.fhcrc.proteomics.schemaRevision.amtXml10.*;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.w3c.dom.Node;
import org.apache.xmlbeans.XmlOptions;

/**
 * A restrictive wrapper for writing AmtXml files.  We take advantage of XmlBeans to build
 * the structure of the amtxmlxml file, and to build individual peptide_entries,
 * but we stitch the XmlBeans XML output for features together by
 * hand, writing out to a file as we go, so that we don't have to hold the whole structure
 * in memory
 */
public class AmtXmlWriter
{
    static Logger _log = Logger.getLogger(AmtXmlWriter.class);

    //document shell structure
    AmtDatabaseDocument mXmlBeansAmtDatabaseDoc = null;
    AmtDatabaseDocument.AmtDatabase mXmlBeansAmtDatabase = null;

    //all entries to be written
    protected AmtDatabase mAmtDatabase;

    //Strings of xml representing the structure before and after the feature content
    protected String _documentPrefix = null;
    protected String _documentPostscript = null;

    //encapsulates printing options for all fragments
    protected XmlOptions _optionsForPrinting = null;


    /**
     * Constructor creates the XmlBeans representing the shell of a AmtXml document, and
     * creates the "prefix" and "postscript" strings representing that shell
     */
    public AmtXmlWriter()
    {
        init();
        generateStubs();
    }

    protected void init()
    {
        //Construct generic document structure
        mXmlBeansAmtDatabaseDoc = AmtDatabaseDocument.Factory.newInstance();
        mXmlBeansAmtDatabase = mXmlBeansAmtDatabaseDoc.addNewAmtDatabase();


    }

    protected void generateStubs()
    {
        //add a sentinel node that tells us where to split the document to insert entries
        Node amtDatabaseNode = mXmlBeansAmtDatabase.getDomNode();
        Node dummyEntryLocationNode =
                amtDatabaseNode.getOwnerDocument().createElement("SENTINEL_FEATURE_LOCATION");
        amtDatabaseNode.appendChild(dummyEntryLocationNode);

        //set printing options for xml fragments
        _optionsForPrinting = new XmlOptions();
        _optionsForPrinting.setSaveOuter();
        _optionsForPrinting.setSavePrettyPrint();
        _optionsForPrinting.setSavePrettyPrintOffset(0);

        //create and break up the xml that defines the document structure
        String documentShell = mXmlBeansAmtDatabaseDoc.xmlText(_optionsForPrinting);
        String[] halves = documentShell.split("<SENTINEL_FEATURE_LOCATION[^\\/]*\\/>");
        if (halves.length != 2)
        {
            _log.error("Failed to create document shell for writing");
            return;
        }

        _documentPrefix = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" + halves[0];
        _documentPostscript = halves[1];

        //remove our dummy node
        amtDatabaseNode.removeChild(dummyEntryLocationNode);
    }


    public AmtXmlWriter(AmtDatabase amtDatabase)
    {
        setAmtDatabase(amtDatabase);
    }

    public void setAmtDatabase(AmtDatabase database)
    {
        init();
        mAmtDatabase = database;
        mXmlBeansAmtDatabase.setHydrophobicityAlgorithm(database.getHydrophobicityAlgorithmName());
        mXmlBeansAmtDatabase.setHydrophobicityAlgVersion(database.getHydrophobicityAlgorithmVersion());        

        generateStubs();
    }

    protected void writeEntries(PrintWriter pw)
    {
        AmtPeptideEntry[] entries = mAmtDatabase.getEntries();
        for (int i=0; i<entries.length; i++)
            writeEntry(entries[i], pw);
    }

    /**
     * Write out a single peptide entry
     * @param entry
     * @param pw
     */
    protected void writeEntry(AmtPeptideEntry entry, PrintWriter pw)
    {
        AmtDatabaseDocument.AmtDatabase.PeptideEntry xmlBeansEntry =
                mXmlBeansAmtDatabase.addNewPeptideEntry();
        xmlBeansEntry.setPeptideSequence(entry.getPeptideSequence());
        xmlBeansEntry.setCalculatedHydrophobicity(BigDecimal.valueOf(entry.getPredictedHydrophobicity()));
        xmlBeansEntry.setMedianObservedHydrophobicity(BigDecimal.valueOf(entry.getMedianObservedHydrophobicity()));
        xmlBeansEntry.setMedianPeptideProphet(BigDecimal.valueOf(entry.getMedianPeptideProphet()));

        for (AmtPeptideEntry.AmtPeptideModificationStateEntry modState :
                entry.getModificationStateEntries())
        {
            AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry xmlBeansModEntry =
                    xmlBeansEntry.addNewModificationStateEntry();
            xmlBeansModEntry.setModifiedSequence(modState.getModifiedSequence());
            xmlBeansModEntry.setModifiedMass(BigDecimal.valueOf(modState.getModifiedMass()));
            xmlBeansModEntry.setMedianObservedHydrophobicity(BigDecimal.valueOf(modState.getMedianObservedHydrophobicity()));
            xmlBeansModEntry.setMedianPeptideProphet(BigDecimal.valueOf(modState.getMedianPeptideProphet()));
            List<MS2Modification>[] modifications = modState.getModifications();

            if (modifications != null)
            {
                for (int i=0; i<modifications.length; i++)
                {
                    List<MS2Modification> modificationList = modifications[i];
                    if (modificationList != null)
                    {
                        for (MS2Modification modification : modificationList)
                        {
                            int modId = mAmtDatabase.getSequenceForAminoacidModification(modification);


                            AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry.AminoacidModInstance
                                    xmlBeansModInstance = xmlBeansModEntry.addNewAminoacidModInstance();
                            xmlBeansModInstance.setModificationId(modId);
                            xmlBeansModInstance.setPosition(i);
                        }
                    }
                }
            }

            for (AmtPeptideEntry.AmtPeptideObservation observation :
                    modState.getObservations())
            {
                AmtDatabaseDocument.AmtDatabase.PeptideEntry.ModificationStateEntry.Observation
                        xmlBeansObservation =
                        xmlBeansModEntry.addNewObservation();
                xmlBeansObservation.setObservedHydrophobicity(BigDecimal.valueOf(observation.getObservedHydrophobicity()));
                xmlBeansObservation.setPeptideProphet(BigDecimal.valueOf(observation.getPeptideProphet()));
                //TODO: getSequenceForRun could return -1, in which case the run isn't found.  That should throw an exception
                xmlBeansObservation.setRunId(mAmtDatabase.getSequenceForRun(observation.getRunEntry()));
                xmlBeansObservation.setTimeInRun(BigDecimal.valueOf(observation.getTimeInRun()));
                if (observation.hasSpectralCount())
                    xmlBeansObservation.setSpectralCount(observation.getSpectralCount());
            }

        }

        try
        {
            String fragment = mXmlBeansAmtDatabase.getPeptideEntryArray(0).xmlText(_optionsForPrinting);
            pw.print(fragment);
            pw.flush();
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        mXmlBeansAmtDatabase.removePeptideEntry(0);
    }

    protected void writeModifications(PrintWriter pw)
    {
        if (mAmtDatabase.numAminoacidModifications() == 0)
            return;
        for (MS2Modification mod : mAmtDatabase.getAminoacidModifications())
            writeModification(mod, pw);
    }

    protected void writeModification(MS2Modification ms2Modification, PrintWriter pw)
    {
        AmtDatabaseDocument.AmtDatabase.AminoacidModification xmlBeansMod =
                mXmlBeansAmtDatabase.addNewAminoacidModification();
        xmlBeansMod.setResidue(ms2Modification.getAminoAcid());
        xmlBeansMod.setMassDifference(BigDecimal.valueOf(ms2Modification.getMassDiff()));
        xmlBeansMod.setVariableFlag(ms2Modification.getVariable());
        xmlBeansMod.setModificationId(mAmtDatabase.getSequenceForAminoacidModification(ms2Modification));
        try
        {
            String fragment = mXmlBeansAmtDatabase.getAminoacidModificationArray(0).xmlText(_optionsForPrinting);
            pw.print(fragment);
            pw.flush();
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        mXmlBeansAmtDatabase.removeAminoacidModification(0);
    }

    /**
     * Write out all the runs
     * @param pw
     */
    protected void writeRuns(PrintWriter pw)
    {
        for (AmtRunEntry run : mAmtDatabase.getRuns())
        {
            writeRun(run, pw);
        }
    }

    /**
     * Write out a single run
     * @param pw
     */
    protected void writeRun(AmtRunEntry run, PrintWriter pw)
    {
        AmtDatabaseDocument.AmtDatabase.Run xmlBeansRun =
            mXmlBeansAmtDatabase.addNewRun();
        //set required attributes
        xmlBeansRun.setRunId(mAmtDatabase.getSequenceForRun(run));
        double[] coeffs = run.getTimeHydMapCoefficients();
        for (int i=0; i<coeffs.length; i++)
        {
            AmtDatabaseDocument.AmtDatabase.Run.TimeHydroMappingCoeff xmlBeansCoeff =
                    xmlBeansRun.addNewTimeHydroMappingCoeff();
            xmlBeansCoeff.setDegree(i);
            xmlBeansCoeff.setValue(BigDecimal.valueOf(coeffs[i]));
        }        Calendar timeAddedCalendar = Calendar.getInstance();
        timeAddedCalendar.setTime(run.getTimeAdded());
        xmlBeansRun.setTimeAdded(timeAddedCalendar);

        //set modifications, if any
        MS2Modification[] modifications = run.getModifications();
        if (modifications != null && modifications.length > 0)
        {
            for (MS2Modification mod : modifications)
            {
                AmtDatabaseDocument.AmtDatabase.Run.AminoacidModUsage modUsage =
                        xmlBeansRun.addNewAminoacidModUsage();
                modUsage.setModificationId(mAmtDatabase.getSequenceForAminoacidModification(mod));
            }
        }

        //set optional attributes if present
        if (run.getMzXmlFilename() != null)
            xmlBeansRun.setMzxmlFilename(run.getMzXmlFilename());
        if (run.getPepXmlFilename() != null)
            xmlBeansRun.setPepxmlFilename(run.getPepXmlFilename());
        if (run.getLSID() != null)
            xmlBeansRun.setLSID(run.getLSID());
        if (run.getMinPeptideProphet() > 0.0)
            xmlBeansRun.setMinPeptideProphet(BigDecimal.valueOf(run.getMinPeptideProphet()));
        if (run.getTimeAnalyzed() != null)
        {
            Calendar timeAnalyzedCalendar = Calendar.getInstance();
            timeAnalyzedCalendar.setTime(run.getTimeAnalyzed());
            xmlBeansRun.setTimeAnalyzed(timeAnalyzedCalendar);
        }

        try
        {
            String fragment = mXmlBeansAmtDatabase.getRunArray(0).xmlText(_optionsForPrinting);
            pw.print(fragment);
            pw.flush();
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
        }
        mXmlBeansAmtDatabase.removeRun(0);        
    }


    /**
     * Write out the full document, with all modifications and features, to a file
     * @param file
     * @throws IOException
     */
    public void write(File file) throws IOException
    {
        PrintWriter pw = new PrintWriter(file);
        pw.print(_documentPrefix);
        writeModifications(pw);
        pw.println("");
        writeRuns(pw);
        pw.println("");        
        writeEntries(pw);
        pw.print(_documentPostscript);
        pw.flush();
    }
}
