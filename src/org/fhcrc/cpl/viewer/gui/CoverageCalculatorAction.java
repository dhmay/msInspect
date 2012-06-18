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
package org.fhcrc.cpl.viewer.gui;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.Protein;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.Peptide;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.FastaLoader;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Map;

/**
 * User: migra
 * Date: Jul 13, 2004
 * Time: 1:08:05 PM
 */
public class CoverageCalculatorAction extends AbstractAction implements PropertyChangeListener
{
    public CoverageCalculatorAction()
    {
        super("Calculate Peptide Coverage");
        ApplicationContext.addPropertyChangeListener("fastaFile", this);
        ApplicationContext.addPropertyChangeListener("featureSelector", this);
        ApplicationContext.addPropertyChangeListener(SharedProperties.FEATURE_RANGES, this);
    }

    private CoverageInfo computeCoverage(Peptide[] peptides, float[] distances, float maxDistance)
    {
        int pepCount = 0;
        int byteCount = 0;
        int nextByte = 0;

        for (int i = 0; i < distances.length; i++)
        {
            if (distances[i] > maxDistance)
                continue;

            pepCount++;
            //Careful not to double-count bytes for overlapping peptides (missed cleavages)
            int offset = peptides[i].getStart() + peptides[i].getLength();
            if (offset > nextByte)
            {
                int len = Math.min(offset - nextByte, peptides[i].getLength());
                byteCount += len;

                nextByte = offset;
            }
        }

        return new CoverageInfo(byteCount, pepCount);
    }

    /**
     * Output a file with feature coverage information for
     * each protein in the current fasta-file and the current set of
     * features.
     * Tab delimted output has following columns:
     * Protein id, distance, length, coveredLength, peptides, coveredPeptides
     * <p/>
     * Search is done according to current feature & fasta settings
     */
    public void actionPerformed(ActionEvent event)
    {
        OpenFastaDialog openFastaDialog = new OpenFastaDialog();
        File fastaFile = openFastaDialog.getFastaFile();
        if (null == fastaFile)
        {
            //JOptionPane.showMessageDialog(null, "No fasta file is open");
            ApplicationContext.errorMessage("No fasta file is open.", null);
            return;
        }

        MSRun run = (MSRun) ApplicationContext.getProperty(SharedProperties.MS_RUN);
        if (null == run)
        {
            //JOptionPane.showMessageDialog(null, "No run file is open");
            ApplicationContext.errorMessage("No run file is open.", null);
            return;
        }

        java.util.List featureSets = (java.util.List) ApplicationContext.getProperty(SharedProperties.FEATURE_RANGES);

        File runFile = run.getFile();
        String outputFileName = runFile.getAbsolutePath() + ".coverage.tsv";

        File outputFile = new File(outputFileName);

        CoverageCalculator2D calculator = new CoverageCalculator2D(fastaFile, featureSets, outputFile, (float) openFastaDialog.getTolerance());
        Thread t = new Thread(calculator);
        t.start();
    }

    public void propertyChange(PropertyChangeEvent event)
    {
        Application app = Application.getInstance();
        setEnabled(null != app.getProperty(SharedProperties.FEATURE_RANGES) && null != app.getProperty("fastaFile"));
    }

    private static class CoverageInfo
    {
        int coveredLength;
        int coveredPeptides;

        public CoverageInfo(int coveredLength, int coveredPeptides)
        {
            this.coveredLength = coveredLength;
            this.coveredPeptides = coveredPeptides;
        }
    }

    private class CoverageCalculator implements Runnable
    {
        private File _fastaFile;
        private java.util.List _featureSets;
        private File _outputFile;
        private float _maxDistance;

        public CoverageCalculator(File fastaFile, java.util.List featureSets, File outputFile, float maxDistance)
        {
            _fastaFile = fastaFile;
            _featureSets = featureSets;
            _outputFile = outputFile;
            _maxDistance = maxDistance;
        }

        public void run()
        {
            PrintWriter out = null;
            FastaLoader.ProteinIterator iterator = null;

            try
            {
            out = new PrintWriter(new FileOutputStream(_outputFile));
            out.print("protein\tlength\tpeptides");
            OpenFastaDialog openFastaDialog = new OpenFastaDialog();
            FastaLoader loader = new FastaLoader(_fastaFile);
            iterator = (FastaLoader.ProteinIterator) loader.iterator();

            float[][] featureSetMassArrays = new float[_featureSets.size()][];
            for (int i = 0; i < _featureSets.size(); i++)
            {
                FeatureSet fs = (FeatureSet) _featureSets.get(i);
                Feature[] features = fs.getFeatures();
                float[] featureMasses = new float[features.length];
                for (int j = 0; j < features.length; j++)
                    featureMasses[j] = features[j].mass;

                Arrays.sort(featureMasses);
                featureSetMassArrays[i] = featureMasses;
            }
            for (int i = 0; i < _featureSets.size(); i++)
                out.print("\tcoveredLength" + i);
            for (int i = 0; i < _featureSets.size(); i++)
                 out.print("\tcoveredPeptides" + i);
            out.println();

            int nProteins = 0;
            while (iterator.hasNext())
            {
                Protein protein = (Protein) iterator.next();
                String proteinId = protein.getLookup();


                PeptideGenerator pepGen = new PeptideGenerator();
                double[] massTab = openFastaDialog.getMassTab();
                pepGen.setMassTable(massTab);
                pepGen.setMaxMissedCleavages(openFastaDialog.getMissedCleavages());

                Peptide[] peptides = pepGen.digestProtein(protein);
                float[] distances = new float[peptides.length];

                out.print(proteinId);
                out.print("\t");
                out.print(protein.getBytes().length);
                out.print("\t");
                out.print(peptides.length);

                CoverageInfo[] coverage = new CoverageInfo[featureSetMassArrays.length];
                for(int i = 0; i < featureSetMassArrays.length; i++)
                {
                    float[] featureMasses = featureSetMassArrays[i];
                    for (int j = 0; j < peptides.length; j++)
                    {
                        double pepMass = peptides[j].getMass();
                        int index = Arrays.binarySearch(featureMasses, (float) pepMass);
                        double distance = Double.MAX_VALUE;
                        if (index >= 0)
                            distance = 0; //found an exact match
                        else
                        {
                            index = -index -1;
                            if (index > 0)
                                distance = Math.abs(pepMass - featureMasses[index - 1]);
                            if (index < featureMasses.length)
                            {
                                double dist2 = Math.abs(pepMass - featureMasses[index]);
                                if (dist2 < distance)
                                    distance = dist2;
                            }
                        }

                        distances[j] = (float) distance;
                    }

                    coverage[i] = computeCoverage(peptides, distances, _maxDistance);
                }

                for (int i = 0; i < coverage.length; i++)
                {
                    out.print("\t");
                    out.print(coverage[i].coveredLength);
                }
                for (int i = 0; i < coverage.length; i++)
                {
                    out.print("\t");
                    out.print(coverage[i].coveredPeptides);
                }
                out.println();

                if (nProteins++ % 100 == 0)
                    ApplicationContext.setMessage(String.valueOf(nProteins - 1) + " proteins processed.");
            }

            iterator = null;
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage(null, x);
            }
            finally
            {
                if (null != out)
                    out.close();
                if (null != iterator)
                    iterator.close();
                ApplicationContext.setMessage("Coverage information is in file " + _outputFile.getAbsolutePath());
            }
        }

    }

    private class CoverageCalculator2D implements Runnable
        {
        private static final String SLOPE_PROPERTY = "hydroSlope";
        private static final String INTERCEPT_PROPERTY = "hydroIntercept";
        private static final String SIGMA_PROPERTY = "hydroSigma";

        private File _fastaFile;
        private java.util.List _featureSets;
        private File _outputFile;
        private float _maxDistance;
        private Tree2D[] _trees;
        private double[] _slope;
        private double[] _intercept;
        private double[] _scanWindow;

        public CoverageCalculator2D(File fastaFile, java.util.List featureSets, File outputFile, float maxDistance)
        {
            _fastaFile = fastaFile;
            _featureSets = featureSets;
            _outputFile = outputFile;
            _maxDistance = maxDistance;
            _slope = new double[featureSets.size()];
            Arrays.fill(_slope, Double.NaN);
            _intercept = new double[featureSets.size()];
            Arrays.fill(_intercept, Double.NaN);
            _scanWindow = new double[featureSets.size()];
            Arrays.fill(_scanWindow, Double.NaN);
        }

        public CoverageCalculator2D(File fastaFile, java.util.List featureSets, File outputFile, float maxDistance,
                                    double slope, double intercept, double sigma)
            {
            this(fastaFile, featureSets, outputFile, maxDistance);
            Arrays.fill(_slope, slope);
            Arrays.fill(_intercept, intercept);
            Arrays.fill(_scanWindow, sigma * 2);
            }


        public void run()
            {

            _trees = new Tree2D[_featureSets.size()];
            //Populate the trees
            for (int i = 0; i < _trees.length; i++)
                {
                Tree2D tree = new Tree2D();
                FeatureSet fs = ((FeatureSet) _featureSets.get(i));
                Map properties = fs.getProperties();
                if (properties.containsKey(SLOPE_PROPERTY))
                    _slope[i] = Double.parseDouble((String) properties.get(SLOPE_PROPERTY));
                if (properties.containsKey(INTERCEPT_PROPERTY))
                    _intercept[i] = Double.parseDouble((String) properties.get(INTERCEPT_PROPERTY));
                if (properties.containsKey(SIGMA_PROPERTY))
                    {
                    String sigma = (String ) properties.get(SIGMA_PROPERTY);
                    if (null != sigma)
                        _scanWindow[i] = 2 * Double.parseDouble(sigma);
                    }

                Feature[] features = fs.getFeatures();
                for (int j = 0; j < features.length; j++)
                    {
                    Feature feature = features[j];
                    tree.add(feature.getScan(),feature.getMass(), feature);
                    }
                _trees[i] = tree;
                }

            PrintWriter out = null;
            FastaLoader.ProteinIterator iterator = null;

            try
            {
            out = new PrintWriter(new FileOutputStream(_outputFile));
            out.print("protein\tlength\tpeptides");
            OpenFastaDialog openFastaDialog = new OpenFastaDialog();
            FastaLoader loader = new FastaLoader(_fastaFile);
            iterator = (FastaLoader.ProteinIterator) loader.iterator();


            for (int i = 0; i < _featureSets.size(); i++)
                out.print("\tcoveredLength" + i);
            for (int i = 0; i < _featureSets.size(); i++)
                 out.print("\tcoveredPeptides" + i);
            out.println();

            int nProteins = 0;
            while (iterator.hasNext())
            {
                Protein protein = (Protein) iterator.next();
                String proteinId = protein.getLookup();


                PeptideGenerator pepGen = new PeptideGenerator();
                double[] massTab = openFastaDialog.getMassTab();
                pepGen.setMassTable(massTab);
                pepGen.setMaxMissedCleavages(openFastaDialog.getMissedCleavages());

                Peptide[] peptides = pepGen.digestProtein(protein);
                float[] distances = new float[peptides.length];

                out.print(proteinId);
                out.print("\t");
                out.print(protein.getBytes().length);
                out.print("\t");
                out.print(peptides.length);

                CoverageInfo[] coverage = new CoverageInfo[_trees.length];
                for(int i = 0; i < _trees.length; i++)
                {
                    Tree2D tree = _trees[i];
                    for (int j = 0; j < peptides.length; j++)
                    {
                        float pepMass = (float) peptides[j].getMass();
                        float pepScan = 0;
                        //TODO: Make this depend on actual size of file.
                        float scanRange = 100000;
                        if (_slope[i] != Double.NaN)
                            {
                            pepScan = (float) (peptides[j].getHydrophobicity() * _slope[i] + _intercept[i]);
                            scanRange = (float) _scanWindow[i];
                            }

                        if (tree.containsPoints(pepScan - scanRange, pepMass - _maxDistance, pepScan + scanRange, pepMass + _maxDistance))
                            distances[i] = 0;
                        else
                            distances[i] = Float.MAX_VALUE;
                    }

                    coverage[i] = computeCoverage(peptides, distances, _maxDistance);
                }

                for (int i = 0; i < coverage.length; i++)
                {
                    out.print("\t");
                    out.print(coverage[i].coveredLength);
                }
                for (int i = 0; i < coverage.length; i++)
                {
                    out.print("\t");
                    out.print(coverage[i].coveredPeptides);
                }
                out.println();

                if (nProteins++ % 100 == 0)
                    ApplicationContext.setMessage(String.valueOf(nProteins - 1) + " proteins processed.");
            }

            iterator = null;
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage(null, x);
            }
            finally
            {
                if (null != out)
                    out.close();
                if (null != iterator)
                    iterator.close();
                ApplicationContext.setMessage("Coverage information is in file " + _outputFile.getAbsolutePath());
            }
        }
    }
}
