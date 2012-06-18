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
package org.fhcrc.cpl.viewer.qa;

import org.fhcrc.cpl.toolbox.filehandler.TabLoader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 * test
 */
public class QAUtilities
{
    protected static Logger _log = Logger.getLogger(QAUtilities.class);

    public static void createAllProtText(File allProtXmlFile, File protGeneFile, File outFile,
                                         float minProteinProphet)
            throws IOException, XMLStreamException
    {

        Map<String,List<String>> ipiGeneListMap = QAUtilities.loadIpiGeneListMap(protGeneFile);

        ApplicationContext.setMessage("Gene map loaded. " + ipiGeneListMap.size() + " proteins with genes");



        String headerLine = "*Group\tGroup_Probability\tProtein_Probability\tL2H_Mean\tL2H_StdDev\tRatio_Peps\tH2L_Mean\tH2L_StdDev\tnum_Indistinguishable_Proteins\tIndistinguishable_Proteins\tnum_Genes\tGene_Symbol\tPeptides";

        PrintWriter pw = new PrintWriter(outFile);
        pw.println(headerLine);

        ProtXmlReader protXmlReader = new ProtXmlReader(allProtXmlFile);
        ProtXmlReader.ProteinGroupIterator groupIterator = protXmlReader.iterator();

        while (groupIterator.hasNext())
        {
            ProteinGroup proteinGroup = groupIterator.next();

            List<ProtXmlReader.Protein> proteins = proteinGroup.getProteins();

            for (int i=0; i<proteins.size(); i++)
            {
                ProtXmlReader.Protein protein = proteins.get(i);

                if (minProteinProphet > 0 && protein.getProbability() < minProteinProphet)
                    continue;

                //For proteins after the first, "group" number is actually group number with _# after it, where _# is
                //the count
                StringBuffer lineBuf = new StringBuffer("" + proteinGroup.getGroupNumber());
                if (i>0)
                    lineBuf.append("_" + i);


                lineBuf.append("\t" + proteinGroup.getGroupProbability());
                lineBuf.append("\t" + protein.getProbability());

                ProtXmlReader.QuantitationRatio ratio = protein.getQuantitationRatio();
                if (ratio == null)
                    lineBuf.append("\t-666\t-666\t-666\t-666\t-666");
                else
                {
                    lineBuf.append("\t" +  ratio.getRatioMean());
                    lineBuf.append("\t" + ratio.getRatioStandardDev());
                    lineBuf.append("\t" + ratio.getRatioNumberPeptides());
                    lineBuf.append("\t" +  ratio.getHeavy2lightRatioMean());
                    lineBuf.append("\t" + ratio.getHeavy2lightRatioStandardDev());
                }

                StringBuffer indistProteinNamesBuf = new StringBuffer("");
                List<String> indistProteinNames = protein.getIndistinguishableProteinNames();
                if (indistProteinNames == null)
                    indistProteinNames = new ArrayList<String>(1);
                indistProteinNames.add(protein.getProteinName());
                for (int j=0; j<indistProteinNames.size(); j++)
                {
                    if (j>0)
                        indistProteinNamesBuf.append(";");
                    indistProteinNamesBuf.append(indistProteinNames.get(j));

                }
                lineBuf.append("\t" + indistProteinNames.size());
                lineBuf.append("\t" + indistProteinNamesBuf);


                //num_genes, gene_symbol
                Set<String> genesAllIndistProteins = new HashSet<String>();
                for (String indistProteinName : indistProteinNames)
                {
                    List<String> geneList = ipiGeneListMap.get(indistProteinName);
                    if (geneList != null)
                        genesAllIndistProteins.addAll(geneList);
                }
                lineBuf.append("\t" + genesAllIndistProteins.size());
                StringBuffer geneSymbolBuf = new StringBuffer("");
                boolean firstGene = true;
                for (String gene : genesAllIndistProteins)
                {
                    if (!firstGene)
                        geneSymbolBuf.append(";");
                    geneSymbolBuf.append(gene);
                    firstGene = false;
                }
                if (geneSymbolBuf.length() == 0)
                    geneSymbolBuf.append("NA");
                lineBuf.append("\t" + geneSymbolBuf);


                StringBuffer peptidesBuf = new StringBuffer("");
                List<ProtXmlReader.Peptide> peptides = protein.getPeptides();
                int numWrittenPeptides = 0;
                Set<String> peptidesThisProtein = new HashSet<String>();
                for (ProtXmlReader.Peptide peptide : peptides)
                {
                    if (peptide.isContributingEvidence() && peptide.isNondegenerateEvidence())
                    {
                        if (numWrittenPeptides>0)
                            peptidesBuf.append(";");
                        peptidesBuf.append(peptide.getPeptideSequence());
                        peptidesThisProtein.add(peptide.getPeptideSequence());
                        numWrittenPeptides++;
                    }
                }

                if (peptidesBuf.length() == 0)
                    peptidesBuf.append("NA");
                lineBuf.append("\t" + peptidesBuf);



                pw.println(lineBuf);
                pw.flush();
            }
        }
        pw.flush();
        pw.close();
    }

    public static Map<String,List<String>> loadGeneIpiListMap(File protGeneFile)
            throws IOException
    {
        Map<String,List<String>> geneIpiListMap = new HashMap<String,List<String>>();
        ApplicationContext.setMessage("Loading gene mapping file...");

        TabLoader tabLoader = new TabLoader(protGeneFile);
        tabLoader.setColumns(new TabLoader.ColumnDescriptor[] {
                new TabLoader.ColumnDescriptor("protein", String.class),
                new TabLoader.ColumnDescriptor("genes", String.class)});

        TabLoader.TabLoaderIterator iter = tabLoader.iterator();
        while (iter.hasNext())
        {
            Map rowMap = (Map) iter.next();
            String genes = (String) rowMap.get("genes");
            if (genes == null)
                continue;
            String[] geneArray = genes.toString().split("//");
            for (String gene : geneArray)
            {
                List<String> proteins = geneIpiListMap.get(gene);
                if (proteins == null)
                {
                    proteins = new ArrayList<String>();
                    geneIpiListMap.put(gene, proteins);
                }
                String protein = (String) rowMap.get("protein");

                proteins.add(protein);
            }
        }
        return geneIpiListMap;
    }



    public static Map<String,List<String>> loadIpiGeneListMap(File protGeneFile)
            throws IOException
    {
        Map<String,List<String>> ipiGeneListMap = new HashMap<String,List<String>>();
        ApplicationContext.setMessage("Loading gene mapping file...");

        TabLoader tabLoader = new TabLoader(protGeneFile);
        tabLoader.setColumns(new TabLoader.ColumnDescriptor[] {
                new TabLoader.ColumnDescriptor("protein", String.class),
                new TabLoader.ColumnDescriptor("genes", String.class)});

        TabLoader.TabLoaderIterator iter = tabLoader.iterator();
        while (iter.hasNext())
        {
            Map rowMap = (Map) iter.next();
            String protein = (String) rowMap.get("protein");
            if (protein != null)
            {
                Object genesStringObj = rowMap.get("genes");
                if (genesStringObj == null)
                    continue;

                String[] geneArray = genesStringObj.toString().split("//");
                List<String> geneList = new ArrayList<String>();
                for (String gene : geneArray)
                {
                    if (gene != null && gene.length() > 0)
                        geneList.add(gene);
                }
                if (!geneList.isEmpty())
                    ipiGeneListMap.put(protein, geneList);
            }
        }
        return ipiGeneListMap;
    }

    



    


}
