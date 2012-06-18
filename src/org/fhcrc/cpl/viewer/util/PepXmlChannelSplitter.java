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
package org.fhcrc.cpl.viewer.util;

import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLEventRewriter;

import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.ArrayList;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.stream.events.Attribute;

/**
 * Split the "heavy" or "light" channel of a pepXML file into a new pepXML
 * file. Strips out analysis results as well (since some, like PeptideProphet,
 * will no longer be valid for a subset).
 *
 * TODO: emits one blank line for each excluded analysis/spectrum block
 */
public class PepXmlChannelSplitter
{
    private static Logger _log = Logger.getLogger(PepXmlChannelSplitter.class);

    public static void split(String inFilename,
                             char labeledResidue,
                             boolean preserveLabeled,
                             String outFilename)
        throws IOException, XMLStreamException
    {
        PepXmlSplitRewriter splitter = new PepXmlSplitRewriter(inFilename,
                                                               labeledResidue,
                                                               preserveLabeled,
                                                               outFilename);
        try
        {
            splitter.rewrite();
        }
        finally
        {
            splitter.close();
        }
    }


    /**
     * 
     */
    static class PepXmlSplitRewriter extends SimpleXMLEventRewriter
    {
        private static final float EPSILON = .001f;

        String inFilename;
        String outFilename;
        char labeledResidue;
        boolean preserveLabeled;
        float labeledMass;

        ArrayList<XMLEvent> spectrumQueryEvents;
        boolean insideAnalysisBlock = false;
        boolean insideSpectrumQuery = false;
        int totalResidueCount;
        int labeledResidueCount;
        String currentPeptide;

        public PepXmlSplitRewriter(String inFilename, char labeledResidue, boolean preserveLabeled, String outFilename)
        {
            super(inFilename, outFilename);
            this.labeledResidue = labeledResidue;
            this.preserveLabeled = preserveLabeled;

            spectrumQueryEvents = new ArrayList<XMLEvent>();
            totalResidueCount = 0;
            labeledResidueCount = 0;
            insideAnalysisBlock = false;
            insideSpectrumQuery = false;
        }

        public void handleStartElement(StartElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();
            // ???? final QName ANALYSIS_SUMMARY = new QName("analysis_summary");
            // ???? ANALYSIS_SUMMARY.equals(qname)

            if ("analysis_summary".equals(qname.getLocalPart()))
            {
                insideAnalysisBlock = true;
            }
            else if ("analysis_result".equals(qname.getLocalPart()))
            {
                insideAnalysisBlock = true;
            }
            else if ("spectrum_query".equals(qname.getLocalPart()))
            {
                insideSpectrumQuery = true;
            }
            else if ("aminoacid_modification".equals(qname.getLocalPart()))
            {
                processModificationDefinition(event);
            }
            else if ("search_hit".equals(qname.getLocalPart()))
            {
                processSearchHit(event);
            }
            else if ("mod_aminoacid_mass".equals(qname.getLocalPart()))
            {
                processModificationMass(event);
            }

            conditionalAdd(event);
        }

        public void handleEndElement(EndElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();

            conditionalAdd(event);

            if ("spectrum_query".equals(qname.getLocalPart()))
            {
                insideSpectrumQuery = false;
                processSpectrumQueryEvents();
                spectrumQueryEvents.clear();
            }
            else if ("analysis_summary".equals(qname.getLocalPart()))
            {
                insideAnalysisBlock = false;
            }
            else if ("analysis_result".equals(qname.getLocalPart()))
            {
                insideAnalysisBlock = false;
            }
        }

        public void handleDefault(XMLEvent event)
            throws XMLStreamException
        {
            conditionalAdd(event);
        }

        private boolean conditionalAdd(XMLEvent event)
            throws XMLStreamException
        {
            
            if (insideAnalysisBlock)
            {
                return false;
            }

            if (insideSpectrumQuery)
            {
                spectrumQueryEvents.add(event);
                return false;
            }

            add(event);
            return false;
        }

        private void processSpectrumQueryEvents()
            throws XMLStreamException
        {
            if (totalResidueCount == 0)
            {
                // Unlabeled; could come from heavy or light
                addSpectrumQueryEvents();
            }
            else if (labeledResidueCount == 0)
            {
                // Light
                if (!preserveLabeled)
                {
                    addSpectrumQueryEvents();
                }
            }
            else if (labeledResidueCount == totalResidueCount)
            {
                // Heavy
                if (preserveLabeled)
                {
                    addSpectrumQueryEvents();
                }
            }
            else
            {
                // Partial!
//                _log.warn("Skipping partially labeled peptide " + currentPeptide);
            }
        }

        private void addSpectrumQueryEvents()
            throws XMLStreamException
        {
            for (XMLEvent e : spectrumQueryEvents)
            {
                conditionalAdd(e);
            }
        }

        private void processModificationDefinition(StartElement event)
        {
            final QName AMINOACID = new QName("aminoacid");
            final QName VARIABLE = new QName("variable");
            final QName MASS = new QName("mass");
            String residue = getAttribute(event, AMINOACID);
            if (residue != null && residue.charAt(0) == labeledResidue)
            {
                String variable = getAttribute(event, VARIABLE);
                String mass = getAttribute(event, MASS);
                if ("Y".equals(variable) && mass != null)
                {
                    labeledMass = Float.parseFloat(mass);
                    _log.info("Using " + labeledMass + "@" + labeledResidue);
                }
            }
        }

        /**
         * Count the number of potential labeled residues in each peptide hit
         */
        private void processSearchHit(StartElement event)
        {
            final QName PEPTIDE = new QName("peptide");
            currentPeptide = getAttribute(event, PEPTIDE);

            totalResidueCount = 0;
            labeledResidueCount = 0;
            for (int i = 0; i < currentPeptide.length(); i++)
            {
                if (currentPeptide.charAt(i) == labeledResidue)
                {
                    totalResidueCount++;
                }
            }
        }

        /**
         * Count the number of labeled residues in each peptide hit
         */
        private void processModificationMass(StartElement event)
        {
            final QName MASS = new QName("mass");
            String mass = getAttribute(event, MASS);
            if (mass != null && Math.abs(Float.parseFloat(mass) - labeledMass) < EPSILON)
            {
                labeledResidueCount++;
            }
        }

        private String getAttribute(StartElement event, QName name)
        {
            Attribute attr = event.getAttributeByName(name);
            if (attr == null)
            {
                return null;
            }
            return attr.getValue();
        }
    }

    /**
     *
     */
    public static void main(String[] av)
    {
        if (av.length < 3)
        {
            System.err.println("PepXmlChannelSplitter [--light|--heavy] inFile outFile");
            System.exit(1);
        }

        boolean preserveHeavy = "--heavy".equals(av[0]);

        try
        {
            PepXmlChannelSplitter.split(av[1], 'C', preserveHeavy, av[2]);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

}
