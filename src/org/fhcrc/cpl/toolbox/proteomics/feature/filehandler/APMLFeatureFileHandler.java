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
import org.fhcrc.cpl.toolbox.proteomics.feature.AnalyzeICAT;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.MS2ExtraInfoDef;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLStreamReader;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlException;
import org.systemsbiology.apmlparser.v2.*;
import org.systemsbiology.apmlparser.v2.datatype.*;
import org.systemsbiology.apmlparser.v2.datatype.Feature;

import javax.xml.stream.XMLStreamException;
import java.io.*;
import java.util.*;

/**
 * File handler for native msInspect feature files
 */
public class APMLFeatureFileHandler extends BaseFeatureSetFileHandler
    implements FeatureSetFileHandler
{
    static Logger _log = Logger.getLogger(APMLFeatureFileHandler.class);

    public static final String FILE_TYPE_NAME = "APML";
    

    protected int minFeaturesCount = 0;

    //TODO:  get rid of these when implementing APML 2.0
    //These are dummy values to assign to all features because APML 1.0 doesn't store them
    protected static final float DUMMY_KL_SCORE = 2;
    protected static final int DUMMY_NUM_PEAKS = 2;

    //standardize quality score names for MS1 features
    public static final String APML_QUALITY_SCORE_NAME_KL = "KL";
    public static final String APML_QUALITY_SCORE_ACCMASS = "accurate_mass_flag";
    public static final String APML_QUALITY_SCORE_SUMSQUARESDIST = "sum_squares_dist";

    //standardize cluster type descriptions
    public static final String APML_CLUSTERS_DESCRIPTION_LABELED_QUANT = "labeled_quantitation";


    //Software parameter names
    public static final String APML_SOFTWARE_PARAM_FEATURE_STRATEGY = "feature_strategy";
    public static final String APML_SOFTWARE_PARAM_QUANT_LABEL = "quantitation_label";




    //TODO: standardize these search score names with Mi-Youn
    //Search score names, as appropriate for an APML file
    public static final String APML_SEARCH_SCORE_NAME_PEPTIDE_PROPHET = "PeptideProphet";

    protected static APMLFeatureFileHandler singletonInstance = null;

    public static APMLFeatureFileHandler getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new APMLFeatureFileHandler();
        return singletonInstance;
    }

    /**
     * Load a FeatureSet
     * @param file
     * @return
     * @throws IOException
     */
    public FeatureSet loadFeatureSet(File file)
            throws IOException
    {
        _log.debug("loadFeatureSet begin");

        DefaultAPMLReaderListener apmlReaderListener = new DefaultAPMLReaderListener();

        APMLReader apmlReader = new APMLReader(apmlReaderListener);
        apmlReader.setValidate(true);
        apmlReader.setReadSingleScanPeaks(false);
        try
        {
            apmlReader.read(file);
        }
        catch (Exception e)
        {
            e.printStackTrace(System.err);
            throw new IOException("Failed to load APML from file " + file.getAbsolutePath() +
                                  ", type=" + e.getClass().getName() +
                                  ", message=" + e.getMessage());
        }

        if (apmlReaderListener.getDataType() != APMLReaderListener.DATA_TYPE_PEAK_LISTS)
        {
            throw new IOException("Tried to load a non-peak-list APML file.  Quitting");
        }

        FeatureSet result = new FeatureSet();
         

        DataProcessing apmlDataProcessing = apmlReaderListener.getDataProcessing();
        result.setProperty("date",apmlDataProcessing.getProcessingDate().toString());
        List<DataProcessing.Software> apmlSoftwares = apmlDataProcessing.getSoftwares();
        boolean hasQuantitation = false;
        AnalyzeICAT.IsotopicLabel quantLabel = null;
        if (apmlSoftwares != null)
        {
            for (DataProcessing.Software apmlSoftware : apmlSoftwares)
            {
                if ("msinspect".equals(apmlSoftware.getName()))
                {
                    String strategy =
                            apmlSoftware.getDataProcessingParam(APML_SOFTWARE_PARAM_FEATURE_STRATEGY);
                    if (strategy != null)
                         result.setProperty("algorithm",strategy);

                    String quantLabelString = apmlSoftware.getDataProcessingParam(APML_SOFTWARE_PARAM_QUANT_LABEL);
                    if (quantLabelString != null)
                    {
                        _log.debug("Found software parameter indicating quantitation");
                        hasQuantitation = true;
                        quantLabel = new AnalyzeICAT.IsotopicLabel(quantLabelString);
                        _log.debug("\tisotopic label: " + quantLabel);
                        result.addExtraInformationType(IsotopicLabelExtraInfoDef.getSingletonInstance());
                    }
                }
            }
            //TODO: other params?
        }

        List<DefaultPeakListListener> apmlPeakListListeners = apmlReaderListener.getPeakListListeners();

        if (apmlPeakListListeners.size() > 1)
        {
            ApplicationContext.infoMessage("WARNING: loading an APML file with more " +
                    "than one peak (feature) list.  All peak lists after the first " +
                    "will be ignored");
        }

        DefaultPeakListListener peakListListener = apmlPeakListListeners.get(0);
        result.setSourceFile(new File(peakListListener.getSource()));

        List<org.fhcrc.cpl.toolbox.proteomics.feature.Feature> msInspectFeatures =
                new ArrayList<org.fhcrc.cpl.toolbox.proteomics.feature.Feature>(peakListListener.getFeatures().size());
        _log.debug("loadFeatureSet loading features");

        //FeatureSet-level stuff to keep track of
        boolean hasPpids = false;
        List<MS2Modification>allMS2Mods = new ArrayList<MS2Modification>();


        Map<org.systemsbiology.apmlparser.v2.datatype.Feature, org.systemsbiology.apmlparser.v2.datatype.Feature>
                lightHeavyFeatureMap =
                new HashMap<org.systemsbiology.apmlparser.v2.datatype.Feature,
                            org.systemsbiology.apmlparser.v2.datatype.Feature>();
        List<BaseDefaultClustersListener> clustersListeners = apmlReaderListener.getClustersListeners();
        boolean hasLabeledPairs = false;
        if (clustersListeners != null && !clustersListeners.isEmpty())
        {
            for (BaseDefaultClustersListener clustersListener : clustersListeners)
            {
                if (APML_CLUSTERS_DESCRIPTION_LABELED_QUANT.equals(clustersListener.getDescription()))
                {
                    hasLabeledPairs = true;
                    DefaultFeatureClustersListener featureClustersListener =
                            (DefaultFeatureClustersListener) clustersListener;

                    for (FeatureCluster featureCluster : featureClustersListener.getFeatureClusters())
                    {
                        List<org.systemsbiology.apmlparser.v2.datatype.Feature> clusterFeatures =
                                featureCluster.getFeatures();
                        if (clusterFeatures.size() != 2)
                            throw new IOException("Labeled quantitation cluster had " + clusterFeatures.size() +
                                                  " members!");
                        if (clusterFeatures.get(0).getCoord().getMz() < clusterFeatures.get(1).getCoord().getMz())
                            lightHeavyFeatureMap.put(clusterFeatures.get(0), clusterFeatures.get(1));
                        else
                            lightHeavyFeatureMap.put(clusterFeatures.get(1), clusterFeatures.get(0));
                    }
                    break;
                }
            }
        }

        if (hasQuantitation && !hasLabeledPairs)
            throw new IOException("File has quantitation, but no labeled pairs found.  Quitting");

        //populate the feature list
        for (org.systemsbiology.apmlparser.v2.datatype.Feature apmlFeature: peakListListener.getFeatures())
        {
            //if this is isotopically-labeled data and we've already identified this as a heavy feature, pass
            if (hasLabeledPairs && lightHeavyFeatureMap.containsValue(apmlFeature))
                continue;

            org.fhcrc.cpl.toolbox.proteomics.feature.Feature msInspectFeature =
                    createMsInspectFeatureFromAPMLFeature(apmlFeature);
            msInspectFeatures.add(msInspectFeature);
            if (MS2ExtraInfoDef.getPeptideList(msInspectFeature) != null)
            {
                hasPpids = true;
                List<ModifiedAminoAcid>[] modsThisPeptide =
                        MS2ExtraInfoDef.getModifiedAminoAcids(msInspectFeature);
                if (modsThisPeptide != null)
                {
                    for (List<ModifiedAminoAcid> modList : modsThisPeptide)
                    {
                        if (modList != null)
                        {
                            for (ModifiedAminoAcid mod : modList)
                            {
                                boolean foundIt = false;
                                for (MS2Modification ms2Mod : allMS2Mods)
                                {
                                    if (ms2Mod.getAminoAcid().charAt(0) == mod.getAminoAcid() &&
                                            ms2Mod.getMassDiff() == mod.getMass())
                                    {
                                        foundIt = true;
                                        break;
                                    }
                                }
                                if (!foundIt)
                                {
                                    MS2Modification ms2Mod = new MS2Modification();
                                    ms2Mod.setAminoAcid("" + mod.getAminoAcid());
                                    ms2Mod.setMassDiff((float)mod.getMass());
                                    allMS2Mods.add(ms2Mod);
                                }
                            }
                        }
                    }
                }
            }


            if (hasQuantitation && hasLabeledPairs && lightHeavyFeatureMap.containsKey(apmlFeature))
            {
                org.systemsbiology.apmlparser.v2.datatype.Feature heavyFeature =
                        lightHeavyFeatureMap.get(apmlFeature);
                float massDiff = heavyFeature.getCoord().getMass() - apmlFeature.getCoord().getMass();
                float labelMassDiff = quantLabel.getHeavy() - quantLabel.getLight();

                int numLabels = (int) (massDiff / labelMassDiff);

                IsotopicLabelExtraInfoDef.setLabel(msInspectFeature, quantLabel);
                IsotopicLabelExtraInfoDef.setLabelCount(msInspectFeature, numLabels);

                float heavyIntensity = heavyFeature.getCoord().getApexIntensity();
                float lightIntensity = apmlFeature.getCoord().getApexIntensity();
                float ratio = lightIntensity / heavyIntensity;
                IsotopicLabelExtraInfoDef.setRatio(msInspectFeature, ratio);
                IsotopicLabelExtraInfoDef.setLightIntensity(msInspectFeature, lightIntensity);
                IsotopicLabelExtraInfoDef.setLightIntensity(msInspectFeature, heavyIntensity);

                IsotopicLabelExtraInfoDef.setHeavyMass(msInspectFeature, heavyFeature.getCoord().getMass());
            }
        }


        result.setFeatures(msInspectFeatures.toArray(
               new org.fhcrc.cpl.toolbox.proteomics.feature.Feature[msInspectFeatures.size()]));
        _log.debug("Loaded " + result.getFeatures().length + " features from APML file.");

        if (hasPpids)
        {
            result.addExtraInformationType(MS2ExtraInfoDef.getSingletonInstance());
        }

        if (allMS2Mods.size() > 0)
        {

            MS2Modification[] modifications =
                     allMS2Mods.toArray(new MS2Modification[0]);
            _log.debug("\tloaded " + (modifications == null ? 0 : modifications.length +
                       " MS2 modifications"));
            MS2ExtraInfoDef.setFeatureSetModifications(result, modifications);
        }

        return result;
    }

    /**
     * This method knows the mapping between APML feature fields and
     * msInspect feature fields
     * @param apmlFeature
     * @return
     */
    public org.fhcrc.cpl.toolbox.proteomics.feature.Feature
              createMsInspectFeatureFromAPMLFeature(
            org.systemsbiology.apmlparser.v2.datatype.Feature apmlFeature)
    {
        org.fhcrc.cpl.toolbox.proteomics.feature.Feature result =
                 new org.fhcrc.cpl.toolbox.proteomics.feature.Feature();

        //handle coordinate information
        Coordinate coord = apmlFeature.getCoord();
        result.setScan(coord.getApexScan());
        Coordinate.Range<Integer> apmlScanRange = coord.getScanRange();
        if (apmlScanRange != null)
        {
            result.setScanFirst(apmlScanRange.getMin());
            result.setScanLast(apmlScanRange.getMax());
        }
        else
        {
            result.setScanFirst(result.getScan());
            result.setScanLast(result.getScan());
        }
        //if APML scancount isn't set, default to 1 scan
        result.setScanCount(Math.max(coord.getScanCount(), 1));
        result.setMz(coord.getMz());
        result.setMass(coord.getMass());
        result.setCharge(coord.getCharge());

        //minTime and maxTime, and minMz and maxMz, are ignored
        result.setTime(coord.getRt());
        result.setIntensity(coord.getApexIntensity());
        result.setTotalIntensity(coord.getIntensity());
        result.setPeaks(apmlFeature.getNumPeaks());

        List<MultiScanPeak> multiScanPeaks = apmlFeature.getMultiScanPeaks();
        if (multiScanPeaks != null && !multiScanPeaks.isEmpty())
        {
            Spectrum.Peak[] msInspectPeaks = new Spectrum.Peak[multiScanPeaks.size()];
            for (int i=0; i<multiScanPeaks.size(); i++)
            {
                MultiScanPeak multiScanPeak = multiScanPeaks.get(i);

                org.fhcrc.cpl.toolbox.proteomics.feature.Feature peakFeature =
                        new org.fhcrc.cpl.toolbox.proteomics.feature.Feature(
                        multiScanPeak.getCoordinate().getApexScan(), multiScanPeak.getCoordinate().getMz(),
                        multiScanPeak.getCoordinate().getApexIntensity());
                peakFeature.setTotalIntensity(multiScanPeak.getCoordinate().getIntensity());
                peakFeature.setScanFirst(multiScanPeak.getCoordinate().getScanRange().getMin());
                peakFeature.setScanLast(multiScanPeak.getCoordinate().getScanRange().getMax());
                msInspectPeaks[i] = peakFeature;
            }
            result.comprised = msInspectPeaks;
            result.peaks = msInspectPeaks.length;
        }

        List<org.systemsbiology.apmlparser.v2.datatype.Feature.QualityScore> apmlQualityScores =
                apmlFeature.getQualityScores();
        if (apmlQualityScores != null && !apmlQualityScores.isEmpty())
        {
            for (org.systemsbiology.apmlparser.v2.datatype.Feature.QualityScore apmlQualityScore : apmlQualityScores)
            {
                String scoreName = apmlQualityScore.getScoreName();
                String scoreValue = apmlQualityScore.getScoreValue();
                if (APML_QUALITY_SCORE_NAME_KL.equals(scoreName))
                    result.setKl(Float.parseFloat(scoreValue));
                else if (APML_QUALITY_SCORE_ACCMASS.equals(scoreName))
                    result.setAccurateMZ(Boolean.parseBoolean(scoreValue));
                else if (APML_QUALITY_SCORE_SUMSQUARESDIST.equals(scoreName))
                    result.setSumSquaresDist(Float.parseFloat(scoreValue));
            }
        }

        if (apmlFeature.getAnnotation() != null)
            result.setDescription(apmlFeature.getAnnotation());


        if (apmlFeature.getPpids() != null && apmlFeature.getPpidsSize() > 0)
        {
            //TODO: handle multiple ppids
            PutativePeptideId ppid = apmlFeature.getPpids().get(0);

            String peptideSequence = ppid.getPeptideSequence();
            MS2ExtraInfoDef.addPeptide(result, peptideSequence);

            List<Modification> apmlMods = ppid.getModifications();
            if (apmlMods != null && !apmlMods.isEmpty())
            {
                List<ModifiedAminoAcid>[] modifiedAminoAcids = new List[peptideSequence.length()];

                for (Modification apmlMod : apmlMods)
                {
                    int position = apmlMod.getPosition();
                    List<ModifiedAminoAcid> modList = modifiedAminoAcids[position];
                    if (modList == null)
                    {
                        modList = new ArrayList<ModifiedAminoAcid>();
                        modifiedAminoAcids[position] = modList;
                    }

                    modList.add(new ModifiedAminoAcid(peptideSequence.charAt(position), apmlMod.getPosition()));
                }

                MS2ExtraInfoDef.setModifiedAminoAcids(result, modifiedAminoAcids);
            }


            //TODO: right now we don't handle multiple possible proteins per peptide
            List<String> proteinAccessionNum = ppid.getProteinAccessionNumbers();
            if (proteinAccessionNum != null && proteinAccessionNum.size() > 0)
                MS2ExtraInfoDef.addProtein(result, proteinAccessionNum.get(0));

            List<PutativePeptideId.SearchScore> ms2SearchScores = ppid.getMs2SearchScores();
            if (ms2SearchScores != null && !ms2SearchScores.isEmpty())
            {
                for (PutativePeptideId.SearchScore ms2SearchScore : ms2SearchScores)
                {
                    String searchScoreName = ms2SearchScore.getScoreName();
                    if (APML_SEARCH_SCORE_NAME_PEPTIDE_PROPHET.equals(searchScoreName))
                    {
                        //todo: what if this doesn't parse?  Should check type?
                        MS2ExtraInfoDef.setPeptideProphet(result,
                                Double.parseDouble(ms2SearchScore.getScoreValue()));
                    }
                }
            }

            List<Modification> apmlModifications = ppid.getModifications();
            if (apmlModifications != null && apmlModifications.size() > 0)
            {
                List<ModifiedAminoAcid>[] modListArray =
                        (List<ModifiedAminoAcid>[])
                                new List[ppid.getPeptideSequence().length()];
                for (Modification apmlMod : apmlModifications)
                {
                    char newModChar = ppid.getPeptideSequence().charAt(apmlMod.getPosition());

                    //TODO: we store these zero-based.  How does APML store them?
                    //TODO: assuming zero-based for now
                    List<ModifiedAminoAcid> modListThisIndex =
                            modListArray[apmlMod.getPosition()];
                    if (modListThisIndex == null)
                    {
                        modListThisIndex =
                                new ArrayList<ModifiedAminoAcid>();
                        modListArray[apmlMod.getPosition()] = modListThisIndex;
                    }

                    ModifiedAminoAcid msInspectMod =
                            new ModifiedAminoAcid(newModChar, apmlMod.getValue());
                    modListThisIndex.add(msInspectMod);
                }
                MS2ExtraInfoDef.setModifiedAminoAcids(result, modListArray);
            }

        }
        return result;
    }


    public void saveFeatureSet(FeatureSet featureSet, File outFile)
            throws IOException
    {
        APMLWriter apmlWriter = new APMLWriter();

        DataProcessing dataProcessing = new DataProcessing();
        dataProcessing.setProcessingDate(new GregorianCalendar());

        DataProcessing.Software msInspectSoftware = new DataProcessing.Software();
        msInspectSoftware.setName("msinspect");
        msInspectSoftware.setType(DataProcessing.Software.TYPE_PEAK_PICKING);
        String featureFindingAlgorithm = (String) featureSet.getProperty("algorithm");
        if (featureFindingAlgorithm != null)
            msInspectSoftware.addDataProcessingParam(APML_SOFTWARE_PARAM_FEATURE_STRATEGY, featureFindingAlgorithm);
        dataProcessing.addSoftware(msInspectSoftware);

        boolean hasQuant = featureSet.hasExtraInformationType(IsotopicLabelExtraInfoDef.getSingletonInstance());

//        dp.setType(ProcessStatus.PEAK_PICKING);
//        dp.setSoftware("msInspect");
//        //todo: what if revision is null?
//        dp.setVersion((String)featureSet.getProperty("revision"));

        org.fhcrc.cpl.toolbox.proteomics.feature.Feature[] msInspectFeatures =
                featureSet.getFeatures();


        ArrayList<org.systemsbiology.apmlparser.v2.datatype.Feature> apmlFeatureList =
                new ArrayList<org.systemsbiology.apmlparser.v2.datatype.Feature>(
                        msInspectFeatures.length);

        int currentFeatureId = 1;
        List<Cluster> quantClusters = new ArrayList<Cluster>();
        int currentClusterId = 1;
        boolean foundLabel = false;
        for (org.fhcrc.cpl.toolbox.proteomics.feature.Feature msInspectFeature : msInspectFeatures)
        {
            org.systemsbiology.apmlparser.v2.datatype.Feature apmlFeature =
                 createAPMLFeature(msInspectFeature);
            apmlFeature.setId(currentFeatureId++);
            apmlFeatureList.add(apmlFeature);

            if (msInspectFeature.comprised != null)
            {

                for (Spectrum.Peak peak : msInspectFeature.comprised)
                {
                    if (peak == null)
                        continue;
                    org.fhcrc.cpl.toolbox.proteomics.feature.Feature peakFeature =
                            (org.fhcrc.cpl.toolbox.proteomics.feature.Feature) peak;
                    Coordinate peakCoordinate = new Coordinate();
                    peakCoordinate.setMz(peakFeature.getMz());
                    peakCoordinate.setApexIntensity(peakFeature.getIntensity());
                    peakCoordinate.setIntensity(peakFeature.getTotalIntensity());
                    peakCoordinate.setScanRange(new Coordinate.Range<Integer>(
                            peakFeature.getScanFirst(), peakFeature.getScanLast()));
                    peakCoordinate.setApexScan(peakFeature.getScan());
                    float mzDiff = peak.getMz() - msInspectFeature.getMz();

                    int peakOffset = Math.round(mzDiff * msInspectFeature.charge);
//if (peakOffset < 0)
//    System.err.println("***" + peak.getMz() + ", " + msInspectFeature.getMz() + ", " + msInspectFeature.charge + " ******" + msInspectFeature);
                    MultiScanPeak multiScanPeak = new MultiScanPeak(peakCoordinate, peakOffset);
                    apmlFeature.addMultiScanPeak(multiScanPeak);
                }
            }

            if (hasQuant && IsotopicLabelExtraInfoDef.hasRatio(msInspectFeature))
            {
                float ratio = (float) IsotopicLabelExtraInfoDef.getRatio(msInspectFeature);

                //create a heavy feature, initially identical to the light one
                org.systemsbiology.apmlparser.v2.datatype.Feature heavyFeature = createAPMLFeature(msInspectFeature);
                heavyFeature.setId(currentFeatureId++);
                apmlFeatureList.add(heavyFeature);

                float lightIntensity = (float) IsotopicLabelExtraInfoDef.getLightIntensity(msInspectFeature);
                float heavyIntensity = (float) IsotopicLabelExtraInfoDef.getHeavyIntensity(msInspectFeature);
                float massDiff = IsotopicLabelExtraInfoDef.getLabel(msInspectFeature).getHeavy() -
                        IsotopicLabelExtraInfoDef.getLabel(msInspectFeature).getLight();
                heavyFeature.getCoord().setMass(apmlFeature.getCoord().getMass() + massDiff);

                float mzDiff = massDiff * IsotopicLabelExtraInfoDef.getLabelCount(msInspectFeature) /
                               msInspectFeature.getCharge();
                heavyFeature.getCoord().setMz(apmlFeature.getCoord().getMz() + mzDiff);
                if (lightIntensity != 0 && heavyIntensity != 0)
                {
                    apmlFeature.getCoord().setApexIntensity(lightIntensity);
                    heavyFeature.getCoord().setApexIntensity(heavyIntensity);
                }
                else
                {
                    heavyFeature.getCoord().setApexIntensity(msInspectFeature.getIntensity() / ratio);
                }

                FeatureCluster labelCluster = new FeatureCluster();
                labelCluster.addFeature(apmlFeature);
                labelCluster.addFeature(heavyFeature);
                labelCluster.setId(currentClusterId++);
                labelCluster.setClassification(APML_CLUSTERS_DESCRIPTION_LABELED_QUANT);
                quantClusters.add(labelCluster);

                if (!foundLabel)
                {
                    msInspectSoftware.addDataProcessingParam(APML_SOFTWARE_PARAM_QUANT_LABEL,
                            IsotopicLabelExtraInfoDef.getLabel(msInspectFeature).toString());
                    foundLabel = true;
                }

            }
        }

        if (hasQuant && !foundLabel)
            throw new IOException("File seems to have quantitation, but no labeled features were found. Quitting");

        String filePath = outFile.getAbsolutePath();
        if (featureSet.getSourceFile() != null)
            filePath = featureSet.getSourceFile().getAbsolutePath();

        try
        {
            if (hasQuant)
                apmlWriter.writePeakListFile(outFile, dataProcessing, apmlFeatureList.iterator(),
                    filePath, apmlFeatureList.size(),
                    1, APML_CLUSTERS_DESCRIPTION_LABELED_QUANT, quantClusters);
            else
                apmlWriter.writePeakListFile(outFile, dataProcessing, apmlFeatureList.iterator(),
                    filePath, apmlFeatureList.size());

        }
        catch (XMLStreamException e)
        {
            ApplicationContext.errorMessage("Failed to save APML",e);
        }
        catch (XmlException xe)
        {
            ApplicationContext.errorMessage("Failed to save APML",xe);
        }
    }

    /**
     * This method knows the mapping of fields from msInspect-style features to
     * APML-style features
     * @param msInspectFeature
     * @return
     */
    public org.systemsbiology.apmlparser.v2.datatype.Feature createAPMLFeature(
            org.fhcrc.cpl.toolbox.proteomics.feature.Feature msInspectFeature)
    {
        org.systemsbiology.apmlparser.v2.datatype.Feature result =
                new org.systemsbiology.apmlparser.v2.datatype.Feature();

        Coordinate coord = new Coordinate();
        coord.setApexScan(msInspectFeature.getScan());
        Coordinate.Range<Integer> scanRange =
                new Coordinate.Range<Integer>(msInspectFeature.getScanFirst(),
                                   msInspectFeature.getScanLast());
        coord.setScanCount(msInspectFeature.getScanCount());
        coord.setScanRange(scanRange);
        coord.setMz(msInspectFeature.getMz());
        coord.setMass(msInspectFeature.getMass());
        coord.setRt(msInspectFeature.getTime());
        coord.setApexIntensity(msInspectFeature.getIntensity());
        coord.setIntensity(msInspectFeature.getTotalIntensity());
        coord.setCharge(msInspectFeature.getCharge());
        result.setCoord(coord);


        Feature.QualityScore klQualityScore = new Feature.QualityScore(APML_QUALITY_SCORE_NAME_KL,
                                                                       Float.toString(msInspectFeature.getKl()),
                                                                       Feature.QualityScore.TYPE_DECIMAL);
        result.addQualityScore(klQualityScore);

        Feature.QualityScore accMassQualityScore =
                new Feature.QualityScore(APML_QUALITY_SCORE_ACCMASS,
                                         Boolean.toString(msInspectFeature.isAccurateMZ()),
                                         Feature.QualityScore.TYPE_BOOLEAN);
        result.addQualityScore(accMassQualityScore);

        Feature.QualityScore sumSquaresQualityScore =
                new Feature.QualityScore(APML_QUALITY_SCORE_SUMSQUARESDIST,
                                         Float.toString(msInspectFeature.getSumSquaresDist()),
                                         Feature.QualityScore.TYPE_DECIMAL);
        result.addQualityScore(sumSquaresQualityScore);


        if (msInspectFeature.getDescription() != null)
            result.setAnnotation(msInspectFeature.getDescription());

        //TODO: support multiple peptides per feature
        String peptide = MS2ExtraInfoDef.getFirstPeptide(msInspectFeature);
        if (peptide != null)
        {
            PutativePeptideId ppid = new PutativePeptideId();
            ppid.setPeptideSequence(peptide);

            if (MS2ExtraInfoDef.hasPeptideProphet(msInspectFeature))
            {
                PutativePeptideId.SearchScore searchScore =
                        new PutativePeptideId.SearchScore(APML_SEARCH_SCORE_NAME_PEPTIDE_PROPHET,
                                                          "" + MS2ExtraInfoDef.getPeptideProphet(msInspectFeature));
                ppid.addMs2SearchScore(searchScore);
            }

            List<ModifiedAminoAcid>[] modifiedAAs = MS2ExtraInfoDef.getModifiedAminoAcids(msInspectFeature);
            if (modifiedAAs != null && modifiedAAs.length > 0)
            {
                for (int i=0; i<modifiedAAs.length; i++)
                {
                    List<ModifiedAminoAcid> modsThisIndex = modifiedAAs[i];
                    if (modsThisIndex != null)
                    {
                        for (ModifiedAminoAcid mod : modsThisIndex)
                            ppid.addModification(new Modification(i, (float) mod.getMass()));
                    }
                }
            }

            ArrayList<PutativePeptideId> ppids =
                    new ArrayList<PutativePeptideId>(1);
            ppids.add(ppid);
            result.setPpids(ppids);
        }

        return result;
    }


    /**
     * Save a FeatureSet
     * @param featureSet
     * @param out
     */
    public void saveFeatureSet(FeatureSet featureSet, PrintWriter out)
    {
        throw new IllegalArgumentException(
                "This version of saveFeatureSet not implemented in APMLFeatureFileHandler");
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
            return false;
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
            if ("apml".equalsIgnoreCase(startElementName))
                result = true;
        }
        catch (XMLStreamException xse)
        {
            throw new IOException(xse.getMessage());
        }
        finally
        {
            if (fis != null)
                fis.close();
        }
        return result;
    }
     
}
