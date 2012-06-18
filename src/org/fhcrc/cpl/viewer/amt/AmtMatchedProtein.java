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

import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.proteomics.Protein;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: dhmay
 * Date: Jun 14, 2007
 * Time: 10:15:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class AmtMatchedProtein
{
    protected String identifier;
    protected int numMatches;
    protected int numUnique;
    protected float meanIntensity;
    protected String header;

    public Map<String, AmtMatchedPeptide> peptideInfoMap;

    public AmtMatchedProtein(String identifier, String header)
    {
        this.identifier = identifier;
        this.header = header;
    }

    public AmtMatchedProtein(Protein protein)
    {
        this(protein.getLookup(), protein.getHeader());
    }


    public AmtMatchedProtein(String identifier, int numMatches,
                             int numUnique, float meanIntensity, String header)
    {
        this(identifier, header);
        this.numMatches = numMatches;
        this.numUnique = numUnique;
        this.meanIntensity = meanIntensity;
    }

    public static String getFileHeaderString()
    {
        return "protein\tpeptide\tintensity\tunique";
    }

    public void addPeptideInfo(AmtMatchedPeptide newPeptideInfo)
    {
        if (peptideInfoMap == null)
            peptideInfoMap = new HashMap<String, AmtMatchedPeptide>(1);
        if (getPeptideSequences().contains(newPeptideInfo.getSequence()))
        {
            AmtMatchedPeptide thisProteinPeptideInfo =
                    getPeptideInfo(newPeptideInfo.getSequence());
            thisProteinPeptideInfo.intensity =
                    (thisProteinPeptideInfo.intensity +
                            newPeptideInfo.intensity)
                            /2;
        }
        else
        {
            peptideInfoMap.put(newPeptideInfo.sequence, newPeptideInfo);
            numMatches++;
            if (newPeptideInfo.uniqueThisProtein)
                numUnique++;
        }
    }

    public Set<String> getPeptideSequences()
    {
        return peptideInfoMap.keySet();
    }

    public AmtMatchedPeptide getPeptideInfo(String sequence)
    {
        return peptideInfoMap.get(sequence);
    }

    public void addMoreInfo(AmtMatchedProtein otherProteinInfo)
    {
        for (String otherProteinPeptide : otherProteinInfo.getPeptideSequences())
        {
            if (getPeptideSequences().contains(otherProteinPeptide))
            {
                AmtMatchedPeptide thisProteinPeptideInfo =
                        getPeptideInfo(otherProteinPeptide);
                thisProteinPeptideInfo.intensity =
                        (thisProteinPeptideInfo.intensity +
                                otherProteinInfo.getPeptideInfo(otherProteinPeptide).intensity)
                                /2;
            }
            else
            {
                addPeptideInfo(otherProteinInfo.getPeptideInfo(otherProteinPeptide));
            }
        }
        recalculateMeanIntensity();
    }

    public void recalculateMeanIntensity()
    {
        float[] peptideIntensities = new float[getPeptideSequences().size()];
        int i=0;
        for (AmtMatchedPeptide peptideInfo : peptideInfoMap.values())
        {
            peptideIntensities[i++] = peptideInfo.intensity;
        }
        meanIntensity = BasicStatistics.mean(peptideIntensities);
    }

    public String oldToString()
    {
        StringBuffer resultBuf = new StringBuffer();

        resultBuf.append(">" + identifier + "\t" +
                numMatches + "\t" + numUnique + "\t" +
                meanIntensity +"\t" +
                0 + "\t" +
                header + "\n");
        for (AmtMatchedPeptide peptideInfo : peptideInfoMap.values())
        {
            resultBuf.append(peptideInfo.toString());
        }
        return resultBuf.toString();
    }

    public String toString(String prefixEachLine)
    {
        StringBuffer resultBuf= new StringBuffer();
        for (AmtMatchedPeptide peptideInfo : peptideInfoMap.values())
        {
            resultBuf.append(prefixEachLine + peptideInfo.toString(identifier));
        }
        return resultBuf.toString();
    }

    public String toString()
    {
        return toString("");
    }


    public String getIdentifier()
    {
        return identifier;
    }

    public void setIdentifier(String identifier)
    {
        this.identifier = identifier;
    }

    public int getNumMatches()
    {
        return numMatches;
    }

    public void setNumMatches(int numMatches)
    {
        this.numMatches = numMatches;
    }

    public int getNumUniqueMatches()
    {
        return numUnique;
    }

    public void setNumUniqueMatches(int numUnique)
    {
        this.numUnique = numUnique;
    }

    public float getMeanIntensity()
    {
        return meanIntensity;
    }

    public void setMeanIntensity(float meanIntensity)
    {
        this.meanIntensity = meanIntensity;
    }

    public String getHeader()
    {
        return header;
    }

    public void setHeader(String header)
    {
        this.header = header;
    }

    public Map<String, AmtMatchedPeptide> getPeptideInfoMap()
    {
        return peptideInfoMap;
    }

    public void setPeptideInfoMap(Map<String, AmtMatchedPeptide> peptideInfoMap)
    {
        this.peptideInfoMap = peptideInfoMap;
    }

    public static Map<String, AmtMatchedProtein> loadMatchedProteins(File file)
            throws IOException
    {

        Map<String, AmtMatchedProtein> result =
                new HashMap<String, AmtMatchedProtein>();
        BufferedReader br = null;
        try
        {
            br = new BufferedReader(new FileReader(file));

            //skip header line
            String line = br.readLine();

            line = br.readLine();
            String currentProteinID = "";
            AmtMatchedProtein currentMatchedProtein = null;
            while (line != null)
            {
                String[] lineFields = line.split("\t");
                String proteinName = lineFields[0];

                if (!proteinName.equals(currentProteinID))
                {
                    currentMatchedProtein = new AmtMatchedProtein(currentProteinID,"");
                    result.put(currentProteinID, currentMatchedProtein);
                    currentProteinID = proteinName;
                }

                String peptideSequence = lineFields[1];
                float intensity = Float.parseFloat(lineFields[2]);
                boolean unique = Boolean.parseBoolean(lineFields[3]);

                AmtMatchedPeptide newMatchedPeptide =
                        new AmtMatchedPeptide(peptideSequence, intensity, unique);

                currentMatchedProtein.addPeptideInfo(newMatchedPeptide);

                line = br.readLine();
            }
        }
        catch (IOException e)
        {
            throw e;
        }
        finally
        {
            if (br != null)
            {
                try { br.close();} catch (Exception e) {}
            }
        }
        return result;
    }



    public static class AmtMatchedPeptide
    {
        protected String sequence;
        protected float intensity;
        protected boolean uniqueThisProtein;

        public AmtMatchedPeptide(String sequence, float intensity, boolean uniqueThisProtein)
        {
            this.sequence = sequence;
            this.intensity = intensity;
            this.uniqueThisProtein = uniqueThisProtein;
        }

        public String toString(String proteinId)
        {
            return proteinId + "\t" + sequence + "\t" + intensity + "\t" + uniqueThisProtein + "\n";
        }

        public String oldToString()
        {
            return sequence + "\t" + intensity + "\t" + uniqueThisProtein + "\n";
        }

        public String getSequence()
        {
            return sequence;
        }

        public void setSequence(String sequence)
        {
            this.sequence = sequence;
        }

        public float getIntensity()
        {
            return intensity;
        }

        public void setIntensity(float intensity)
        {
            this.intensity = intensity;
        }

        public boolean isUniqueThisProtein()
        {
            return uniqueThisProtein;
        }

        public void setUniqueThisProtein(boolean uniqueThisProtein)
        {
            this.uniqueThisProtein = uniqueThisProtein;
        }
    }
}
