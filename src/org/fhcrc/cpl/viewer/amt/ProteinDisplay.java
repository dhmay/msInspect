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

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.BrowserController;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;

import java.util.*;
import java.net.URL;

/**
 * Methods for HTML display of protein sequences.  Also, methods for creating URLs to access
 * online resources for protein information
 */
public class ProteinDisplay
{
    // constants for html protein display
    protected static int DISPLAY_TYPE_BASE=0;
    protected static int DISPLAY_TYPE_MATCHED_MS1=1;
    protected static int DISPLAY_TYPE_MATCHED_MS2=2;
    protected static int DISPLAY_TYPE_MATCHED_MS1_AND_MS2=3;
    protected static int DISPLAY_TYPE_SELECTED=4;

    //display modes
    public static final int DISPLAY_MATCHED_FEATURES_MODE=0;
    public static final int DISPLAY_UNMATCHED_MS2_FEATURES_MODE=1;
    public static final int DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE=2;


    private static Logger _log = Logger.getLogger(ProteinDisplay.class);

    /**
     * For a given protein and feature list, create an HTML string to display
     * the protein
     * @return
     */
    public static String getProteinSequenceHtml(Protein protein, ArrayList<Feature> ms1MatchedFeatures,
                                                ArrayList<Feature> ms2MatchedFeatures,
                                                ArrayList<Feature> ms1AndMS2MatchedFeatures,
                                                ArrayList<Feature> selectedFeatures,
                                                int featureDisplayMode)
    {
        String proteinSequence = protein.getSequenceAsString();
        int[] characterTypes = new int[proteinSequence.length()];
        Arrays.fill(characterTypes, DISPLAY_TYPE_BASE);
        boolean[] underlines = new boolean[proteinSequence.length()];
        for (int i=0; i<underlines.length; i++)
            underlines[i] = false;

        if (featureDisplayMode == DISPLAY_MATCHED_FEATURES_MODE)
        {
            highlightFeatures(characterTypes, ms1MatchedFeatures, proteinSequence,
                              DISPLAY_TYPE_MATCHED_MS1);
            underlineFeatures(underlines, ms1MatchedFeatures, proteinSequence);
            if (ms2MatchedFeatures != null)
                highlightFeatures(characterTypes, ms2MatchedFeatures, proteinSequence, DISPLAY_TYPE_MATCHED_MS2);
            if (ms1AndMS2MatchedFeatures != null)
                highlightFeatures(characterTypes, ms1AndMS2MatchedFeatures, proteinSequence, DISPLAY_TYPE_MATCHED_MS1_AND_MS2);
        }

        highlightFeatures(characterTypes, selectedFeatures, proteinSequence, DISPLAY_TYPE_SELECTED);

        ArrayList<String> snippetsArray = createHtmlSnippets(characterTypes, underlines, proteinSequence);
        return createHtmlString(snippetsArray);
    }

    /**
     * Return html strings for the ms1MatchedFeatures that highlight them appropriately
     * @param displayedFeatures
     * @param ms1AndMS2MatchedFeatures
     * @param selectedFeatures
     * @param featureDisplayMode
     * @return
     */
    public static String[] getHtmlForFeatures(  ArrayList<Feature> displayedFeatures,
                                                ArrayList<Feature> ms1AndMS2MatchedFeatures,
                                                ArrayList<Feature> selectedFeatures,
                                                int featureDisplayMode)
    {
        if (displayedFeatures == null) return null;
        String[] result = new String[displayedFeatures.size()];
        int[] displayTypes = new int[displayedFeatures.size()];
        for (int i = 0; i < displayedFeatures.size(); i++)
        {
            displayTypes[i] = DISPLAY_TYPE_BASE;
            Feature currentFeature = displayedFeatures.get(i);
            if (featureDisplayMode == DISPLAY_MATCHED_FEATURES_MODE)
            {
                displayTypes[i] = DISPLAY_TYPE_MATCHED_MS1;
                if (peptideIsPartOf(MS2ExtraInfoDef.getFirstPeptide(currentFeature),ms1AndMS2MatchedFeatures))
                {
                    displayTypes[i] = DISPLAY_TYPE_MATCHED_MS1_AND_MS2;
                }
            }
            else if (featureDisplayMode == DISPLAY_UNMATCHED_MS2_FEATURES_MODE)
            {

            }
            else if (featureDisplayMode == DISPLAY_UNMATCHED_PROTEIN_PEPTIDES_MODE)
            {

            }

            if (selectedFeatures != null && selectedFeatures.contains(currentFeature))
                displayTypes[i] = DISPLAY_TYPE_SELECTED;
        }
        ArrayList<String> throwaway = new ArrayList<String>(1);
        throwaway.add("");
        for (int i=0; i<displayedFeatures.size(); i++)
        {
            String snippet = createHtmlSnippet(MS2ExtraInfoDef.getFirstPeptide(displayedFeatures.get(i)),
                                               displayTypes[i],
                                               true);
            throwaway.set(0,snippet);
            result[i] = createHtmlString(throwaway);
        }
        return result;
    }


    protected static boolean peptideIsPartOf(String peptide, ArrayList<Feature> featureList)
    {
        if (featureList == null || peptide == null)
            return false;
        for (int i=0; i<featureList.size(); i++)
        {
            if (peptide.equals(MS2ExtraInfoDef.getFirstPeptide(featureList.get(i))))
                return true;
        }
        return false;
    }

    /**
     * Calculate the percent coverage of a protein by the peptides in a featureset.
     * The only reason this is in ProteinDisplay is because it makes use of methods I'd
     * already written for displaying the protein sequence
     * @param features
     * @param proteinSequence
     * @return
     */
    public static double calculatePercentCovered(ArrayList<Feature> features,
                                                 String proteinSequence)
    {
        int[] characterTypes = new int[proteinSequence.length()];
        Arrays.fill(characterTypes, DISPLAY_TYPE_BASE);

        highlightFeatures(characterTypes, features, proteinSequence, DISPLAY_TYPE_SELECTED);

        int numResiduesCovered = 0;
        for (int characterType : characterTypes)
        {
            if (characterType == DISPLAY_TYPE_SELECTED)
                numResiduesCovered++;
        }
        //don't forget to convert to a percent
        return (100.0 * (double) numResiduesCovered / (double) proteinSequence.length());
    }

    /**
     * Find a set of peptides from a set of features and mark them, obliterating whatever
     * marking they had already.  So, note:  last highlight takes all.
     * @param characterTypes
     * @param features
     * @param proteinSequence
     * @param highlightStyle
     */
    protected static void highlightFeatures(int[] characterTypes, ArrayList<Feature> features,
                                            String proteinSequence, int highlightStyle)
    {
        if (features == null)
            return;

        for (int i=0; i<features.size(); i++)
        {
            String currentPeptideString = MS2ExtraInfoDef.getFirstPeptide(features.get(i));
            int currentPeptideLength = currentPeptideString.length();
            int pepIndex = proteinSequence.indexOf(currentPeptideString);
            while (pepIndex > 0)
            {
                for (int j=pepIndex; (j-pepIndex)<currentPeptideLength; j++)
                    characterTypes[j] = highlightStyle;

                int restIndex = pepIndex + currentPeptideLength;
                pepIndex = proteinSequence.indexOf(currentPeptideString, restIndex);
            }
        }
    }

    /**
     * Indicate that certain features' residues should be underlined
     * @param underlines
     * @param features
     * @param proteinSequence
     */
    protected static void underlineFeatures(boolean[] underlines, ArrayList<Feature> features,
                                            String proteinSequence)
    {
        if (features == null)
            return;

        for (int i=0; i<features.size(); i++)
        {
            String currentPeptideString = MS2ExtraInfoDef.getFirstPeptide(features.get(i));
            int currentPeptideLength = currentPeptideString.length();
            int pepIndex = proteinSequence.indexOf(currentPeptideString);
            while (pepIndex > 0)
            {
                for (int j=pepIndex; (j-pepIndex)<currentPeptideLength; j++)
                    underlines[j] = true;

                int restIndex = pepIndex + currentPeptideLength;
                pepIndex = proteinSequence.indexOf(currentPeptideString, restIndex);
            }
        }
    }


    /**
     * From an array of html snippets, create an html string for the whole protein.  This can
     * include extra formatting, since each snippet is independent
     * @param snippetsArray
     * @return
     */
    protected static String createHtmlString(ArrayList<String> snippetsArray)
    {
        StringBuffer resultSB = new StringBuffer();
        resultSB.append("<html><pre>");
        for (int i=0; i<snippetsArray.size(); i++)
        {
           resultSB.append(snippetsArray.get(i));
           //insert speces every 10th residue
           int ioffset = i+1;
           if (ioffset > 0 && (ioffset % 10 == 0))
            if (ioffset % 60 == 0)
                resultSB.append("\n");
            else
                resultSB.append(' ');
        }
       resultSB.append("</pre></html>");
        return resultSB.toString();
    }

    /**
     * Create an arraylist of html snippets for each residue.  Styles for each character are controlled
     * by the characterTypes and undlines arrays
     * @param characterTypes
     * @param proteinSequence
     * @return
     */
    protected static ArrayList<String> createHtmlSnippets(int[] characterTypes,
                                                          boolean[] underlines,
                                                          String proteinSequence)
    {

        ArrayList<String> result = new ArrayList<String>(proteinSequence.length());

        for (int i=0; i<characterTypes.length; i++)
        {
            int displayType = characterTypes[i];
            char residue = proteinSequence.charAt(i);
            String snippet = createHtmlSnippet("" + residue,displayType,underlines[i]);
            result.add(snippet);
        }
        return result;
    }

    /**
     * Create an html snippet for a series of residues (may be just one residue)
     * @param residues
     * @param displayType
     * @param shouldUnderline
     * @return
     */
    public static String createHtmlSnippet(String residues, int displayType,
                                           boolean shouldUnderline)
    {
        String snippet = "";
        if (displayType == DISPLAY_TYPE_BASE)
            snippet = "" + residues;
        else
        {
            boolean shouldBold = false;
            String colorString = "#000000";
            if (displayType == DISPLAY_TYPE_MATCHED_MS1)
            {
                //yellow
                colorString = "#D0D000";
            }
            else if (displayType == DISPLAY_TYPE_MATCHED_MS2)
            {
                //red
                colorString = "#FF0000";
            }
            else if (displayType == DISPLAY_TYPE_MATCHED_MS1_AND_MS2)
            {
                //Orange, bold
                colorString = "#FFA000";
                shouldBold = true;
            }
            else if (displayType == DISPLAY_TYPE_SELECTED)
            {
                //blue
                colorString = "#0000FF";
            }
            snippet = "<font color=\"" + colorString + "\">" + residues + "</font>";
            if (shouldBold)
                snippet = "<b>" + snippet + "</b>";
            if (shouldUnderline)
                snippet = "<u>" + snippet + "</u>";
        }
        return snippet;
    }




    //base URL of the online NCBI Entrez viewer
    protected static final String NCBIEntrezViewerBaseUrl="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi";

    /**
     * Create a URL that can be used to query the NCBI database
     * @param queryValue
     * @return
     */
    public static URL createNCBIEntrezProteinURL(String queryValue)
    {
        String query = "db=protein&val=" + queryValue;
        URL result = null;
        try
        {
            result = new URL(NCBIEntrezViewerBaseUrl + "?" + query);
        }
        catch (Exception e){}
        return result;
    }

    /**
     * Go through all the lookups associated with this protein, in descending preference order
     * based on the likelihood of finding the data at NCBI
     * @param protein
     */
    public static void openNCBIBrowserWindow(Protein protein)
    {
        Map<String, Set<String>> identifiers = protein.getIdentifierMap();

        String queryValue=null;

        if (identifiers.containsKey("SwissProt"))
            queryValue = identifiers.get("SwissProt").iterator().next();
        else if (identifiers.containsKey("SwissProtAccn"))
            queryValue = identifiers.get("SwissProtAccn").iterator().next();
        else if (identifiers.containsKey("IPI"))
            queryValue = identifiers.get("IPI").iterator().next();
        else if (identifiers.containsKey("ENSEMBL"))
            queryValue = identifiers.get("ENSEMBL").iterator().next();

        if (queryValue != null)
        {
            try
            {
                BrowserController.navigate(createNCBIEntrezProteinURL(queryValue));
            }
            catch (Exception e) {}
        }
    }
}

