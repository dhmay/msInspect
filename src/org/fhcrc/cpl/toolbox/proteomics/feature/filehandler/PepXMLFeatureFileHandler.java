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
package org.fhcrc.cpl.toolbox.proteomics.feature.filehandler;

import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeaturePepXmlWriter;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLStreamReader;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.QuantitationUtilities;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.*;
import org.apache.log4j.Logger;

import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * File handler for native msInspect feature files
 */
public class PepXMLFeatureFileHandler extends BaseFeatureSetFileHandler
    implements FeatureSetFileHandler
{
    static Logger _log = Logger.getLogger(PepXMLFeatureFileHandler.class);

    protected int firstSpectrumQueryIndex = 1;

    public static final int QUANT_ALGORITHM_Q3 = 0;
    public static final int QUANT_ALGORITHM_XPRESS = 1;

    public static final String FILE_TYPE_NAME = "PEPXML";

    protected static PepXMLFeatureFileHandler singletonInstance = null;

    protected String _searchEngine = BasePepXmlWriter.SEARCH_ENGINE_XTANDEM_COMET;


    public static PepXMLFeatureFileHandler getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new PepXMLFeatureFileHandler();
        return singletonInstance;
    }

    /**
     * Load all FeatureSets
     * @param file
     * @return
     * @throws IOException
     */
    public List<FeatureSet> loadAllFeatureSets(File file)
            throws IOException
    {
        PepXMLFeatureSetIterator fsi = new PepXMLFeatureSetIterator(file);
        List<FeatureSet> result = new ArrayList<FeatureSet>();
        while (fsi.hasNext())
            result.add(fsi.next());
        return result;
    }

    /**
     * Load first FeatureSet in a file
     * @param file
     * @return
     * @throws IOException
     */
    public FeatureSet loadFeatureSet(File file)
            throws IOException
    {
        PepXMLFeatureSetIterator fsi = new PepXMLFeatureSetIterator(file);
        return fsi.next();
    }


    /**
     * Get an iterator on the FeatureSets in a pepXML file
     */
    public static class PepXMLFeatureSetIterator implements Iterator<FeatureSet>
    {
        protected PepXmlLoader pepXmlLoader;
        PepXmlLoader.FractionIterator fractionIterator;
        protected File sourceFile;

        public PepXMLFeatureSetIterator(File file)
                throws IOException
        {
            try
            {
                sourceFile = file;
                pepXmlLoader = new PepXmlLoader(sourceFile, _log);
                _log.debug("Instantiated PepXmlLoader");
                fractionIterator = pepXmlLoader.getFractionIterator();
                _log.debug("Got FractionIterator");
            }
            catch (XMLStreamException xmlse)
            {
                xmlse.printStackTrace(System.err);
                throw new IOException("XML Stream Failure while loading XML file.  Message: " + xmlse.getMessage());
            }
        }

        public void remove()
        {
            //do nothing
        }

        public FeatureSet next()
        {
            _log.debug("Accessing next fraction");
            FeatureSet newFeatureSet = new FeatureSet();
            newFeatureSet.setSourceFile(sourceFile);
            newFeatureSet.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());

            PepXmlLoader.PepXmlFraction fraction = fractionIterator.next();
            PepXMLFeatureFileHandler.getSingletonInstance().setFeatureSetPropertiesFromFraction(
                    fraction, newFeatureSet);
            Feature[] features =
                    PepXMLFeatureFileHandler.getSingletonInstance().getFeaturesFromPepXmlFraction(
                            fraction, pepXmlLoader, newFeatureSet).toArray(new Feature[0]);
            newFeatureSet.setFeatures(features);
            return newFeatureSet;
        }

        public boolean hasNext()
        {
            return fractionIterator.hasNext();
        }
    }

    /**
     * nterm and cterm mods handled same way as aminoacid
     * @param fraction
     * @param featureSet
     */
    protected void setFeatureSetPropertiesFromFraction(PepXmlLoader.PepXmlFraction fraction,
                                                       FeatureSet featureSet)
    {
        MS2Modification[] modifications =
                fraction.getModifications().toArray(new MS2Modification[0]);
        if (modifications != null && modifications.length > 0)
        {
            for (MS2Modification ms2Mod : modifications)
            {
                if (!ms2Mod.getAminoAcid().equals("n") && !ms2Mod.getAminoAcid().equals("c"))
                    MS2ExtraInfoDef.updateMS2ModMassOrDiff(ms2Mod);
            }
        }
        MS2ExtraInfoDef.correctMS2ModMasses(modifications);
        _log.debug("\tloaded " + (modifications == null ? 0 : modifications.length + " MS2 modifications"));
        MS2ExtraInfoDef.setFeatureSetModifications(featureSet, modifications);

        String databasePath = fraction.getDatabaseLocalPath();
        if (databasePath != null)
            MS2ExtraInfoDef.setFeatureSetSearchDatabasePath(featureSet, databasePath);
        int maxCleavages = fraction.getSearchConstraintMaxInternalCleavages();
        if (maxCleavages > 0)
            MS2ExtraInfoDef.setFeatureSetSearchConstraintMaxIntCleavages(featureSet, maxCleavages);
        int minTermini = fraction.getSearchConstraintMaxInternalCleavages();
        if (minTermini > 0)
            MS2ExtraInfoDef.setFeatureSetSearchConstraintMinTermini(featureSet, minTermini);

        String baseName = fraction.getDataBasename();
        if (baseName != null)
            MS2ExtraInfoDef.setFeatureSetBaseName(featureSet, baseName);
        _log.debug("\tmin termini=" + minTermini + ", max cleavages=" + maxCleavages);

    }

    public FeatureSet createFeatureSetFromPepXMLFraction(PepXmlLoader.PepXmlFraction fraction,
                                                       PepXmlLoader pepXmlLoader)
    {
        FeatureSet result = new FeatureSet();
        List<Feature> featureList =
                getFeaturesFromPepXmlFraction(fraction, pepXmlLoader, result);
        result.setFeatures(featureList.toArray(new Feature[featureList.size()]));
        setFeatureSetPropertiesFromFraction(fraction, result);        
        return result;
    }

    /**
     * get a list of all features contained in a PepXmlFraction.
     * Exposing this for a tool that just pulls one fraction from a pepXml file.
     *
     * Adds the feature information types inferred from the feature list to resultFeatureSet
     * @param fraction
     * @return
     */
    public List<Feature> getFeaturesFromPepXmlFraction(PepXmlLoader.PepXmlFraction fraction,
                                                       PepXmlLoader pepXmlLoader,
                                                       FeatureSet resultFeatureSet)
    {
        _log.debug("getFeaturesFromPepXmlFraction 1");
        boolean hasQuant = false;
        int quantAlgorithmType = -1;
        AnalyzeICAT.IsotopicLabel quantIsotopicLabel = null;

        List<RelativeQuantAnalysisSummary> quantSummaryList =
                pepXmlLoader.getQuantSummaries();
        if (quantSummaryList != null && quantSummaryList.size() == 1)
        {
            _log.debug("Has quantitation summary.  Determining type...");
            try
            {
                RelativeQuantAnalysisSummary quantSummary = quantSummaryList.get(0);
                String massDiffValueString = quantSummary.getMassDiff();
                hasQuant=true;
                resultFeatureSet.addExtraInformationType(IsotopicLabelExtraInfoDef.getSingletonInstance());
                String quantAlgorithmName = quantSummary.getAnalysisAlgorithm();
                if (quantAlgorithmName.contains(Q3AnalysisSummary.analysisType))
                {
                    quantAlgorithmType = QUANT_ALGORITHM_Q3;
                    IsotopicLabelExtraInfoDef.setFeatureSetAlgorithm(resultFeatureSet,
                            QuantitationUtilities.ALGORITHM_Q3);
                }
                else
                if (quantAlgorithmName.contains(XPressAnalysisSummary.analysisType))
                {
                    quantAlgorithmType = QUANT_ALGORITHM_XPRESS;
                    IsotopicLabelExtraInfoDef.setFeatureSetAlgorithm(resultFeatureSet,
                            QuantitationUtilities.ALGORITHM_XPRESS);
                }
                _log.debug("Successfully got basic quantitation information.  Algorithm: " +
                        (quantAlgorithmType==QUANT_ALGORITHM_Q3 ? "Q3" : "XPress"));

                double quantMassDiff = 0;
                //Q3 seems to stick the labeled residue in the mass diff string
                try
                {
                    quantMassDiff = Float.parseFloat(massDiffValueString);
                }
                catch (Exception e)
                {
                    massDiffValueString = massDiffValueString.substring(massDiffValueString.indexOf(',') + 1);
                    quantMassDiff = Float.parseFloat(massDiffValueString);
                }
                char quantLabeledResidue = quantSummary.getLabeledResidues().charAt(0);
                if (quantLabeledResidue == 'n')
                {
                    //TODO: handle n-terminal labels better
                    ApplicationContext.setMessage("Warning: n-terminal label declared.  We don't handle that very well yet");
                    quantLabeledResidue = ' ';
                }
                double lightMass = PeptideGenerator.getMasses(true)[quantLabeledResidue];

                quantIsotopicLabel =
                        new AnalyzeICAT.IsotopicLabel((float) lightMass,
                                (float) (lightMass + quantMassDiff),
                                quantLabeledResidue, 3);
            }
            catch (Exception e)
            {
                if (hasQuant)
                    ApplicationContext.setMessage("WARNING: quantitation information loaded, but couldn't get details of label.  Error: " +
                            e.getMessage());
                else
                    ApplicationContext.setMessage("WARNING: Error loading pepXML quantitation information.  Quantitation will be unavailable for this file. Error: " +
                            e.getMessage());
                e.printStackTrace(System.err);
            }
        }

        PepXmlLoader.PeptideIterator pi = fraction.getPeptideIterator();


        List<Feature> featuresList = new ArrayList<Feature>();

        while (pi.hasNext())
        {

            PepXmlLoader.PepXmlPeptide peptide = pi.next();

            Feature feature = createFeatureFromPepXmlPeptide(peptide, hasQuant, quantAlgorithmType, quantIsotopicLabel);
            if (feature != null)
                featuresList.add(feature);
        }
        return featuresList;
    }

    public Feature createFeatureFromPepXmlPeptide(PepXmlLoader.PepXmlPeptide peptide, boolean hasQuant,
                                                  int quantAlgorithmType, AnalyzeICAT.IsotopicLabel quantIsotopicLabel)
    {
        ModifiedAminoAcid[] modArray = peptide.getModifiedAminoAcids();
        List<ModifiedAminoAcid>[] modificationListArray = null;
        if (modArray !=null)
        {
            //necessary cast
            modificationListArray =
                    (List<ModifiedAminoAcid>[])
                            new ArrayList[modArray.length];

            for (int i=0; i<modArray.length; i++)
                if (modArray[i] != null)
                {
                    modificationListArray[i] =
                            new ArrayList<ModifiedAminoAcid>();
                    modificationListArray[i].add(modArray[i]);
                }

        }


        List<String> peptideList = new ArrayList<String>(1);
        peptideList.add(peptide.getTrimmedPeptide());

        List<String> proteinList = new ArrayList<String>();
        proteinList.add(peptide.getProtein());
        proteinList.addAll(peptide.getAlternativeProteins());

        List<Integer> altProteinNTTs = peptide.getAlternativeProteinNTTs();

        Feature currentFeature = MS2ExtraInfoDef.createMS2Feature(
                peptide.getScan(),
                (float) peptide.getCalculatedNeutralMass(),
                peptide.getCharge(),
                peptideList,
                proteinList,
                modificationListArray);

        if (peptide.getNTerminalModMass() != 0)
            MS2ExtraInfoDef.setNtermModMass(currentFeature, peptide.getNTerminalModMass());
        if (peptide.getCTerminalModMass() != 0)
            MS2ExtraInfoDef.setNtermModMass(currentFeature, peptide.getCTerminalModMass());        

        MS2ExtraInfoDef.setAltProteinNTTs(currentFeature, altProteinNTTs);

        MS2ExtraInfoDef.setDeltaMass(currentFeature, peptide.getDeltaMass());
        MS2ExtraInfoDef.setSearchScores(currentFeature, peptide.getScores());
        String prevAA = peptide.getPrevAA();
        String nextAA = peptide.getNextAA();
        if (prevAA != null)
            MS2ExtraInfoDef.setPrevAminoAcid(currentFeature, prevAA.charAt(0));
        if (peptide.getNextAA() != null)
            MS2ExtraInfoDef.setNextAminoAcid(currentFeature, nextAA.charAt(0));

//System.err.println("prevaa: " + MS2ExtraInfoDef.getPrevAminoAcid(currentFeature) + ", nextaa: " + MS2ExtraInfoDef.getNextAminoAcid(currentFeature));



        int numTrypticEnds = peptide.getNumTolTerm();
        //Try to load NTT from pepXML file, but if it wasn't there, try to recalculate it.  -1 is a sentinel for
        //"wasn't there"
        if (numTrypticEnds < 0)
        {
            //WARNING WARNING WARNING!!!
            //This behavior is trypsin-specific.  If another enzyme is used, number of enzymatic
            //ends will be set incorrectly.
            String peptideSequence = peptide.getTrimmedPeptide();
            //check for start of protein sequence or trypsin digestion at start of peptide  (remember proline)
            if ((prevAA != null && prevAA.startsWith("-")) ||
                    ((prevAA != null && (prevAA.startsWith("K") || prevAA.startsWith("R"))) && !peptideSequence.startsWith("P")))
                numTrypticEnds++;
            //check for end of protein sequence or trypsin digestion at end of peptide
            if ((nextAA != null && nextAA.startsWith("-")) ||
                    ((nextAA != null && !nextAA.startsWith("P")) &&(peptide.getTrimmedPeptide().endsWith("K") ||
                            peptide.getTrimmedPeptide().endsWith("R"))))
                numTrypticEnds++;
        }
        MS2ExtraInfoDef.setNumEnzymaticEnds(currentFeature, numTrypticEnds);


        //If retention time is available, grab it.  If not, time will have its
        //default value (0).  It's possible to have a mix of set and unset values
        //in the same FeatureSet.
        if (peptide.getRetentionTime() != null)
            currentFeature.setTime(peptide.getRetentionTime().floatValue());

//System.err.println("Delta mass: " + peptide.getDeltaMass() + ", " +MS2ExtraInfoDef.getDeltaMass(currentFeature));
        //analysis results:
        //peptideProphet: pprophet
        PeptideProphetHandler.PeptideProphetResult ppar = peptide.getPeptideProphetResult();
        if (null != ppar)
        {
            MS2ExtraInfoDef.setPeptideProphet(currentFeature,ppar.getProbability());
            MS2ExtraInfoDef.setFval(currentFeature, ppar.getProphetFval());
            String allNttProb = ppar.getAllNttProb();
            if (allNttProb != null)
                MS2ExtraInfoDef.setAllNttProb(currentFeature, allNttProb);
        }

        if (hasQuant)
        {
            boolean foundResult = false;
            switch(quantAlgorithmType)
            {
                case QUANT_ALGORITHM_Q3:
                    Q3Handler.Q3Result q3ar = peptide.getQ3Result();
                    if (null != q3ar)
                    {                     
                        IsotopicLabelExtraInfoDef.setRatio(currentFeature,
                                q3ar.getDecimalRatio());
                        IsotopicLabelExtraInfoDef.setHeavyIntensity(currentFeature,
                                q3ar.getHeavyArea());
                        IsotopicLabelExtraInfoDef.setLightIntensity(currentFeature,
                                q3ar.getLightArea());
                        currentFeature.setIntensity(q3ar.getLightArea());
                        currentFeature.setTotalIntensity(q3ar.getLightArea());
                        IsotopicLabelExtraInfoDef.setLightMass(currentFeature,
                                q3ar.getLightMass());
                        IsotopicLabelExtraInfoDef.setHeavyMass(currentFeature,
                                q3ar.getHeavyMass());
                        IsotopicLabelExtraInfoDef.setLightFirstScan(currentFeature, q3ar.getLightFirstscan());
                        IsotopicLabelExtraInfoDef.setLightLastScan(currentFeature, q3ar.getLightLastscan());
                        IsotopicLabelExtraInfoDef.setHeavyFirstScan(currentFeature, q3ar.getHeavyFirstscan());
                        IsotopicLabelExtraInfoDef.setHeavyLastScan(currentFeature, q3ar.getHeavyLastscan());

                        foundResult = true;
                    }
                    break;
                case QUANT_ALGORITHM_XPRESS:
                    XPressHandler.XPressResult xpar = peptide.getXPressResult();
                    if (null!= xpar)
                    {
                        //This is a hack.  XPress returns an intensity that's the sum of
                        //all the intensity values for one peak over its elution profile.
                        //Really that doesn't correspond exactly to our notion of intensity,
                        //or of total intensity.  But I'm setting both.
                        currentFeature.setTotalIntensity(xpar.getLightArea());
                        currentFeature.setIntensity(xpar.getLightArea());
                        IsotopicLabelExtraInfoDef.setRatio(currentFeature,
                                xpar.getDecimalRatio());
                        IsotopicLabelExtraInfoDef.setHeavyIntensity(currentFeature,
                                xpar.getHeavyArea());
                        IsotopicLabelExtraInfoDef.setLightIntensity(currentFeature,
                                xpar.getLightArea());
                        IsotopicLabelExtraInfoDef.setLightFirstScan(currentFeature, xpar.getLightFirstscan());
                        IsotopicLabelExtraInfoDef.setLightLastScan(currentFeature, xpar.getLightLastscan());                        
                        IsotopicLabelExtraInfoDef.setHeavyFirstScan(currentFeature, xpar.getHeavyFirstscan());
                        IsotopicLabelExtraInfoDef.setHeavyLastScan(currentFeature, xpar.getHeavyLastscan());

                        IsotopicLabelExtraInfoDef.setHeavyMass(currentFeature, xpar.getHeavyMass());
                        IsotopicLabelExtraInfoDef.setLightMass(currentFeature, xpar.getLightMass());


                        foundResult = true;
                    }
                    break;
            }
            if (foundResult)
            {
                IsotopicLabelExtraInfoDef.setLabel(currentFeature, quantIsotopicLabel);
            }

        }

        //TODO: This is arbitrary.  Is it excusable??
        //pepxml files don't have intensity, unless xpress results are present.
        //Some software featfeatures, like align,
        //require intensity strictly >100 to work correctly
        if (currentFeature.getIntensity() <= 0)
            currentFeature.setIntensity(200);
        return currentFeature;
    }

    public void saveFeatureSet(FeatureSet featureSet, File outFile)
            throws IOException
    {
        FeaturePepXmlWriter pepXmlWriter = new FeaturePepXmlWriter(featureSet);
        String baseName = outFile.getName();
        if (baseName.contains("."))
            baseName = baseName.substring(0, baseName.indexOf("."));
        pepXmlWriter.setBaseName(baseName);
        pepXmlWriter.setFirstSpectrumQueryIndex(firstSpectrumQueryIndex);
        pepXmlWriter.set_searchEngine(_searchEngine);
        try
        {
            pepXmlWriter.write(outFile);
        }
        catch (Exception e)
        {
            _log.error("Failed to save pepXML",e);
        }
    }

    /**
     * 
     * @param pepXmlFiles
     * @param outFile
     * @throws IOException
     */
    public void combinePepXmlFiles(List<File> pepXmlFiles, File outFile)
            throws IOException
    {
        PrintWriter outPW = null;
        _log.debug("Combining " + pepXmlFiles.size() + " pepXML files into one file...");
        try
        {
            outPW = new PrintWriter(outFile);

            for (int i=0; i<pepXmlFiles.size(); i++)
            {
                File pepXmlFile = pepXmlFiles.get(i);
                _log.debug("\tProcessing file " + pepXmlFile.getAbsolutePath());
                FileReader fr = new FileReader(pepXmlFile);
                BufferedReader br = new BufferedReader(fr);
                String line = null;

                boolean pastHeader = false;
                boolean reachedFooter = false;
                while ((line = br.readLine()) != null)
                {
                    if (line.contains("/msms_pipeline_analysis"))
                        reachedFooter = true;
                    else if (!pastHeader && line.contains("msms_run_summary"))
                    {
                        //if no base_name given
                        if (line.contains("<msms_run_summary>"))
                        {
                            String baseName = pepXmlFile.getName();
                            if (baseName.contains(".") && baseName.length() > baseName.indexOf(".") + 1)
                                baseName = baseName.substring(baseName.indexOf(".") + 1);
                            line = "<msms_run_summary base_name=\"" + baseName + "\">";
                        }
                        pastHeader = true;
                    }

                    if ((i==0 || pastHeader) && (!reachedFooter || i==pepXmlFiles.size()-1))
                        outPW.println(line);
                    outPW.flush();
                    if (reachedFooter && i<pepXmlFiles.size()-1)
                        break;
                }
                outPW.flush();
                _log.debug("Saved pepXML file " + outFile.getAbsolutePath());
            }
        }
        catch (Exception e)
        {
            _log.error("Failed to save pepXML",e);
            throw new IOException("Failed to save pepXML file, " + e.getMessage());
        }
        finally
        {
            if (outPW != null)
                outPW.close();
        }
    }



    /**
     * Write an "all.pep.xml" file
     * @param featureSets
     * @param outFile
     * @throws IOException
     */
    public void saveFeatureSets(List<FeatureSet> featureSets, File outFile)
            throws IOException
    {
        try
        {
            List<File> tempFiles = new ArrayList<File>();
            for (int i=0; i<featureSets.size(); i++)
            {
                FeatureSet featureSet = featureSets.get(i);

                FeaturePepXmlWriter pepXmlWriter =
                        new FeaturePepXmlWriter(featureSet);
                pepXmlWriter.setFirstSpectrumQueryIndex(firstSpectrumQueryIndex);
                pepXmlWriter.set_searchEngine(_searchEngine);                
                File tempFile = TempFileManager.createTempFile("saveFeatureSet"+i+".tmp", this);
                pepXmlWriter.write(tempFile);
                tempFiles.add(tempFile);
            }
            combinePepXmlFiles(tempFiles, outFile);
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            throw new IOException(e.getMessage());
        }
        finally
        {
            TempFileManager.deleteTempFiles(this);
        }
    }

    /**
     * Save a FeatureSet
     * @param featureSet
     * @param out
     */
    public void saveFeatureSet(FeatureSet featureSet, PrintWriter out)
    {
        throw new IllegalArgumentException(
                "This version of saveFeatureSet not implemented in PepXMLFeatureFileHandler");
    }

    /**
     * Can this type of file handler handle this specific file?
     *
     * @param file
     * @return
     * @throws IOException
     */
    public boolean canHandleFile(File file)
        throws IOException
    {

        if (!isXMLFile(file))
        {
            _log.debug("canHandleFile, File is not XML");
            return false;
        }
        _log.debug("canHandleFile, File is XML...");
        FileInputStream fis = null;
        boolean result = false;
        try
        {
            fis = new FileInputStream(file);
            SimpleXMLStreamReader parser = new SimpleXMLStreamReader(fis);

            while (!parser.isStartElement())
                parser.next();
            String startElementName = parser.getLocalName();

            //check that the first element is an msms_pipeline_analysis.  I'm pretty
            //sure that this is required by the pepXML spec, but whether it is or not,
            //if we run into files where this doesn't hold true, will need to change.
            if ("msms_pipeline_analysis".equalsIgnoreCase(startElementName))
                result = true;
            else
            {
                _log.debug("canHandleFile, First element is not msms_pipeline_analysis... it's " + startElementName);
                return false;
            }
        }
        catch (XMLStreamException xse)
        {
            _log.debug("canHandleFile, throwing exception with message " + xse.getMessage());            
            throw new IOException(xse.getMessage());
        }
        finally
        {
            if (fis != null)
                fis.close();
        }
        _log.debug("canHandleFile, returning true!");
        return result;
    }


    public int getFirstSpectrumQueryIndex()
    {
        return firstSpectrumQueryIndex;
    }

    public void setFirstSpectrumQueryIndex(int firstSpectrumQueryIndex)
    {
        this.firstSpectrumQueryIndex = firstSpectrumQueryIndex;
    }

    public String getSearchEngine()
    {
        return _searchEngine;
    }

    public void setSearchEngine(String _searchEngine)
    {
        this._searchEngine = _searchEngine;
    }
}
