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
package org.fhcrc.cpl.viewer.ms2;

import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;


import java.io.*;
import java.util.*;


/**
 * This class doesn't do parsing of ProtXML files itself.  It uses GeneMappingReader for that.
 *
 * This is for utilities to work with the output of GeneMappingReader
 */
public class GeneMappingUtilities
{
    protected static Logger _log = Logger.getLogger(GeneMappingUtilities.class);

    public static Map<String,List<String>> loadIPIGeneMap(File proteinGeneMapFile)
            throws CommandLineModuleExecutionException
    {
        BufferedReader br = null;
        Map<String,List<String>> ipiGeneMap = new HashMap<String,List<String>>();

        try
        {
            br = new BufferedReader(new FileReader(proteinGeneMapFile));
            String line = null;
            while ((line = br.readLine()) != null)
            {
                String[] words = line.split("\t");
                String ipi = words[0];
                if (words.length < 2 || words[1] == null || words[1].length() < 1)
                    continue;
                String geneListString = words[1];

                String[] geneArray = geneListString.split("//");

                List<String> geneList = new ArrayList<String>();
                for (String gene : geneArray)
                    geneList.add(gene);

                ipiGeneMap.put(ipi, geneList);
            }

            ApplicationContext.setMessage("Loaded " + ipiGeneMap.size() + " genes from lookup file");
        }
        catch (Exception e)
        {
            throw new CommandLineModuleExecutionException("Failed to load gene file: ", e);
        }
        finally
        {
            if (br != null)
                try {br.close();} catch (Exception e) {}
        }
        return ipiGeneMap;
    }

    public static class InfoForGene
    {
        protected  List<Double> ratios;
        protected  List<Double> intensities1;
        protected  List<Double> intensities2;


        protected Set<String> proteins;
        protected List<String> peptides;
        protected String symbol;

        public InfoForGene(String symbol, List<String> peptides, List<Double> ratios,
                           List<Double> intensities1, List<Double> intensities2,
                           Set<String> proteins)
        {
            this.symbol = symbol;
            this.ratios = ratios;
            this.proteins = proteins;
            this.peptides = peptides;
            this.intensities1 = intensities1;
            this.intensities2 = intensities2;
        }

        public InfoForGene(String symbol)
        {
            this(symbol, new ArrayList<String>(), new ArrayList<Double>(),
                    new ArrayList<Double>(), new ArrayList<Double>(),
                    new HashSet<String>()
                    );
        }


        public List<Double> getRatios()
        {
            return ratios;
        }

        public void setRatios(List<Double> ratios)
        {
            this.ratios = ratios;
        }


        public List<Double> getIntensities1()
        {
            return intensities1;
        }

        public void setIntensities1(List<Double> intensities1)
        {
            this.intensities1 = intensities1;
        }

        public List<Double> getIntensities2()
        {
            return intensities2;
        }

        public void setIntensities2(List<Double> intensities2)
        {
            this.intensities2 = intensities2;
        }

        public Set<String> getProteins()
        {
            return proteins;
        }

        public void setProteins(Set<String> proteins)
        {
            this.proteins = proteins;
        }

        public String getSymbol()
        {
            return symbol;
        }

        public void setSymbol(String symbol)
        {
            this.symbol = symbol;
        }


        public List<String> getPeptides()
        {
            return peptides;
        }

        public void setPeptides(List<String> peptides)
        {
            this.peptides = peptides;
        }

        public String toString()
        {
            List<Double> ratiosWithoutZeroes = new ArrayList<Double>();
            boolean hasNumeratorNonzero = false;

//System.err.println("ratios: " + getRatios().size() + ", i1: " + intensities1.size() + ", i2: " + intensities2.size());

            for (int i=0; i<getRatios().size(); i++)
            {
//System.err.println("\ti1: " + intensities1.get(i));
                if (intensities1.get(i) > 0)
                {
                    hasNumeratorNonzero = true;
                    if (intensities2.get(i) > 0)
                        ratiosWithoutZeroes.add(ratios.get(i));
                }

            }
            if (ratiosWithoutZeroes.size() == 0)
            {
                if (hasNumeratorNonzero)
                    ratiosWithoutZeroes.add(20.0);
                else
                    ratiosWithoutZeroes.add(0.0);
            }

            int numNumeratorPeptides = 0;
            int numDenomenatorPeptides = 0;
            int numNumerDenomPeptides = 0;
            for (int i=0; i<ratios.size(); i++)
            {
                boolean numer = false;
                boolean denom = false;
                if (intensities1.get(i) > 0)
                    numer = true;
                if (intensities2.get(i) > 0)
                    denom = true;
                if (numer)
                {
                    numNumeratorPeptides++;
                    if (denom) numNumerDenomPeptides++;
                }
                if (denom) numDenomenatorPeptides++;
            }

            double[] ratiosForMean = new double[ratios.size()];
            for(int i=0; i<ratios.size(); i++)
                if (ratios.get(i) == 0)
                    ratiosForMean[i] = 0.0001;
                else
                    ratiosForMean[i] = ratios.get(i);

            double meanRatio = BasicStatistics.geometricMean(ratiosForMean);

            return getSymbol() + "\t" +
                    meanRatio + "\t" +
                    BasicStatistics.geometricMean(ratiosWithoutZeroes) + "\t" +
                    numNumeratorPeptides + "\t" + numDenomenatorPeptides + "\t" +
                    numNumerDenomPeptides + "\t" +
                    BasicStatistics.standardDeviation(ratiosWithoutZeroes) + "\t" +
                    getProteins().size() + "\t" +
                    (BasicStatistics.mean(getIntensities1()) +
                            BasicStatistics.mean(getIntensities2()));
        }

        public static void writeGeneRatioFile(Collection<InfoForGene> infosForGenes, File outFile)
                throws CommandLineModuleExecutionException
        {
            PrintWriter outPW = null;

            try
            {
                outPW = new PrintWriter(outFile);
                outPW.println("gene\tratio\tratio_no_zeroes\tnum_numerator_peptides\tnum_denominator_peptides\tnum_num_denom_peptides\trationozero_std_dev\tnum_proteins\tmean_sum_int1int2");

                for (InfoForGene infoForGene : infosForGenes)
                {
                    outPW.println(infoForGene.toString());
                    outPW.flush();
                }
            }
            catch(Exception e)
            {
                throw new CommandLineModuleExecutionException(e);
            }
            finally
            {
                if (outPW != null)
                    outPW.close();
            }
        }
    }
}
