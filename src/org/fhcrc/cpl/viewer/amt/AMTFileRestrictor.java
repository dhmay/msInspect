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

import java.util.*;
import java.io.*;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProtXmlReader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.ProteinGroup;
import org.apache.log4j.Logger;

public class AMTFileRestrictor
{
    private static Logger _log = Logger.getLogger(AMTFileRestrictor.class);


    public static float DEFAULT_MIN_PROTEIN_PROPHET_FOR_RESTRICTION = .95f;
    /**
     * Load all proteins with a minimum protein prophet score from a protxml file
     * @return
     */
    public static ArrayList<Protein> restrictProteinListWithProtXml(List<Protein> fullProteinList,
                                                              File protXmlFile,
                                                              float minProbability,
                                                              int minPeptidesPerProtein)
    {
        ArrayList<Protein> proteinsInProtXml = new ArrayList<Protein>();
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);

        try
        {
            Iterator<ProteinGroup> iterator = protXmlReader.iterator();

            while (iterator.hasNext())
            {
                ProteinGroup group = iterator.next();
                if (group.getGroupProbability() < minProbability)
                    continue;
                for (int i = 0; i < group.getProteins().size(); i++)
                {
                    // have to check all proteins in the group
                    ProtXmlReader.Protein protXmlReaderProtein = group.getProteins().get(i);
                    for (int j = 0; j < fullProteinList.size(); j++)
                    {
                        Protein currentProtein = fullProteinList.get(j);

                        // protXml protein_name should be the whole start of the protein header
                        // up to ' ' or carriage return
                        if (currentProtein.getHeader().startsWith(protXmlReaderProtein.getProteinName()) &&
                            (currentProtein.getHeader().length() == protXmlReaderProtein.getProteinName().length() ||
                             currentProtein.getHeader().charAt(protXmlReaderProtein.getProteinName().length()) == ' '))
                        {
                            if (protXmlReaderProtein.getTotalNumberPeptides() >= minPeptidesPerProtein)
                                proteinsInProtXml.add(currentProtein);
                            break;
                        }
                    }
                }
            }
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            return null;
        }
        return proteinsInProtXml;
    }

    /**
     * remove all features that match features that were used to support a given list
     * of proteins in a protxml file
     * @return
     */
    public static ArrayList<Feature> removeFeaturesWithProteinsFromProtXml(
                                                              Feature[] inputFeatures,
                                                              List<Protein> proteinList,
                                                              File protXmlFile)
    {
        //our output feature list starts out with everything
        ArrayList<Feature> outputFeatureList = new ArrayList<Feature>(inputFeatures.length);
//ArrayList<Feature> matchedFeatureList = new ArrayList<Feature>(inputFeatures.length);
        for (Feature feature : inputFeatures)
            outputFeatureList.add(feature);
        ArrayList<Feature> identifiedMatchingFeatures = identifyCorrespondingProteinFeaturesInProtXml(
                inputFeatures, proteinList, protXmlFile);
        for (Feature matchingFeature : identifiedMatchingFeatures)
        {
            outputFeatureList.remove(matchingFeature);
//matchedFeatureList.add(matchingFeature);
        }
/*
//HACK!!!
FeaturePepXmlWriter pepXmlWriter = new FeaturePepXmlWriter(matchedFeatureList.toArray(new Feature[0]),
  null);
try {pepXmlWriter.write(new File("/home/dhmay/temp/outunmatched.pep.xml"));} catch (Exception e){}
//end HACK!!!
*/
        return outputFeatureList;

    }


    /**
     * identify all features that match peptides that were used to support a given list
     * of proteins in a protxml file
     * @return
     */
    public static ArrayList<Feature> identifyCorrespondingProteinFeaturesInProtXml(
            Feature[] inputFeatures,
            List<Protein> proteinList,
            File protXmlFile)
    {
        ArrayList<Feature> result = new ArrayList<Feature>();
        ProtXmlReader protXmlReader = new ProtXmlReader(protXmlFile);
int[] proteinBuckets = new int[30];
int[] peptideBuckets = new int[30];
HashMap <String,Integer> peptideHitMap = new HashMap<String,Integer>();
        try
        {
            Iterator<ProteinGroup> iterator = protXmlReader.iterator();

            while (iterator.hasNext())
            {
                ProteinGroup group = iterator.next();

                for (int i = 0; i < group.getProteins().size(); i++)
                {
                    // have to check all proteins in the group
                    ProtXmlReader.Protein protXmlReaderProtein = group.getProteins().get(i);
                    for (Protein currentProtein: proteinList)
                    {
                        // protXml protein_name should be the whole start of the protein header
                        // up to ' ' or carriage return
                        if (currentProtein.getHeader().startsWith(protXmlReaderProtein.getProteinName()) &&
                            (currentProtein.getHeader().length() == protXmlReaderProtein.getProteinName().length() ||
                             currentProtein.getHeader().charAt(protXmlReaderProtein.getProteinName().length()) == ' '))
                        {
int numPeptidesForProtein=0;
                            //found one of our proteins!  identify & remove all the peptides
                            List<ProtXmlReader.Peptide> peptideList =
                                    protXmlReaderProtein.getPeptides();
                            for (ProtXmlReader.Peptide peptide : peptideList)
                            {
                                if (!peptide.isContributingEvidence())
                                    continue;
                                //for loop on inputfeatures, not outputfeatures, so messing
                                //with the list doesn't mess up the loop
                                for (Feature feature : inputFeatures)
                                {
                                    if (peptide.getPeptideSequence().equalsIgnoreCase(MS2ExtraInfoDef.getFirstPeptide(feature)) &&
                                        peptide.getCharge() == feature.getCharge() &&
                                        peptide.getCalcNeutralPepMass() == feature.getMass())
                                    {
numPeptidesForProtein++;
Integer whee = peptideHitMap.get(peptide.getPeptideSequence());
if (whee == null) peptideHitMap.put(peptide.getPeptideSequence(), new Integer(1));
else peptideHitMap.put(peptide.getPeptideSequence(),new Integer(whee+1));
                                        result.add(feature);
      //System.err.println("Removed " + feature.getPeptide());
                                    }
                                }
                            }
proteinBuckets[numPeptidesForProtein]++;

                        }
                    }
                }
            }
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            return null;
        }
int a=0;
for (String whee : peptideHitMap.keySet())
{
    a++;
    peptideBuckets[peptideHitMap.get(whee)]++;
}
System.err.println("peptides: " + a);
//for (String peptide : peptideHitMap.keySet()) System.err.println(peptide);
for (int i=0; i<20; i++)
    System.err.println("Peptides: " + i + ", proteins: " + proteinBuckets[i]);
for (int i=0; i<20; i++)
    System.err.println("Proteins: " + i + ", peptides: " + peptideBuckets[i]);
        return result;

    }
}
