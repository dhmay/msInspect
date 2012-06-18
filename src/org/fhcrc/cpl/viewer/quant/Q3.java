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
package org.fhcrc.cpl.viewer.quant;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.proteomics.ModifiedAminoAcid;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.proteomics.PeptideGenerator;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PeptideProphetHandler;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader.FractionIterator;
// import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader.MS2TerminalModification;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader.PepXmlFraction;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader.PepXmlPeptide;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.PepXmlLoader.PeptideIterator;
import org.fhcrc.cpl.toolbox.proteomics.filehandler.RelativeQuantAnalysisResult;
import org.fhcrc.cpl.toolbox.filehandler.SimpleXMLEventRewriter;
import org.fhcrc.cpl.toolbox.proteomics.MSRun.MSScan;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;

import Jama.Matrix;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.stream.events.Attribute;

/**
 * Apply Marc Coram's q3 quantitation method to a pepXML file.
 *
 * TODO: output minPeptideProphet and maxFracDeltaMass
 *
 * TODO: Check for additional exceptions:
 * - Missing/invalid pepXML file
 * - Missing/invalid mzXML file
 * - No variable mod for indicated residue character
 */
public class Q3
{
    private static Logger _log = Logger.getLogger(Q3.class);

    private static final String VERSION = "1.22a";
    private static final String AUTHOR = "Marc Coram";
    private static final double PROTON_MASS = Spectrum.HYDROGEN_ION_MASS;
    private static final double MASS_MATCH_THRESHOLD = 0.01;
    private static final double PEAK_MATCH_THRESHOLD = 25E-6;

    private static final String outputSuffix = "_q3.pep.xml"; // Must end .pep.xml so we can import into CPAS

    private static final double[] monoMasses = PeptideGenerator.getMasses(true);

    private String pepXmlFilename;
    private String q3XmlFilename;

    private Map<Character,IsotopicLabel> labels = null;

    private String alternateMzXmlDir = null; // Check another directory for associated mzXML files

    // Filtering options
    private float minPeptideProphet;
    private boolean filterByMinPeptideProphet;

    private float maxFracDeltaMass;
    private boolean maxFracDeltaMassIsPPM;
    private boolean filterByMaxFracDeltaMass;

    // Flags to control behavior
    private boolean forceOutput = false; // Force overwrite of existing output file
    private boolean mimicXpress = false; // Output an XPRESS-like file
    private boolean noSentinels = false; // Do not use sentinel values for ill-defined ratios
    private boolean debugMode = false;   // Output extra attributes for debugging
    private boolean compatMode = true;   // Replicate behavior of the R code when center scan has too few matches
    private boolean stripExistingQ3 = false; // If there are existing Q3 analysis_results, remove them

    /**
     * Construct a new Q3 processor
     *
     * @param pepXmlFilename Path to the pepXML file to process.
     * @param labels 
     * @param q3XmlFilename Path to the output pepXML to write.
     */
    public Q3(String pepXmlFilename,
              Map<Character,IsotopicLabel> labels,
              String q3XmlFilename)
    {
        this.pepXmlFilename = pepXmlFilename;
        this.labels = labels;
        this.q3XmlFilename = q3XmlFilename;
    }

    /**
     * Construct a new Q3 processor; old-school entry point supports
     * only a single labeled residue per run
     *
     * @param pepXmlFilename Path to the pepXML file to process.
     * @param labels 
     * @param q3XmlFilename Path to the output pepXML to write.
     */
    public Q3(String pepXmlFilename,
              char labeledResidue,
              float massdiff,
              String q3XmlFilename)
    {
        this.pepXmlFilename = pepXmlFilename;
        this.q3XmlFilename = q3XmlFilename;
        this.labels = new HashMap<Character,IsotopicLabel>();
        this.labels.put(labeledResidue, new IsotopicLabel(labeledResidue, massdiff));
    }

    /**
     * Set another directory to check for mzXML files for each fraction
     *
     * @param alternateMzXmlDir Directory to check
     */
    public void setAlternateMzXmlDir(String alternateMzXmlDir)
    {
        _log.debug("Will search for mzXML files in " + alternateMzXmlDir);
        this.alternateMzXmlDir = alternateMzXmlDir;
    }

    /**
     * Filter out peptides with PeptideProphet probability lower than given threshold
     *
     * @param minPeptideProphet Lowest PeptideProphet probability to consider
     */
    public void setMinPeptideProphet(float minPeptideProphet)
    {
        this.minPeptideProphet = minPeptideProphet;
        this.filterByMinPeptideProphet = true;
    }

    /**
     * Filter out peptides with fractional delta mass greater than the given threshold
     *
     * @param maxFracDeltaMass Largest fractional delta mass to allow
     * @param isPPM If true, maxFracDeltaMass is considered to be in PPM
     */
    public void setMaxFracDeltaMass(float maxFracDeltaMass, boolean isPPM)
    {
        this.maxFracDeltaMass = maxFracDeltaMass;
        this.maxFracDeltaMassIsPPM = isPPM;
        this.filterByMaxFracDeltaMass = true;
    }

    /**
     * Control overwriting of output pepXML file
     *
     * @param forceOutput If true, will force output of pepXML file even if q3XmlFilename already exists
     */
    public void setForceOutput(boolean forceOutput)
    {
        this.forceOutput = forceOutput;
    }

    /**
     * Strip out existing Q3 analysis_results blocks and analysis summaries
     * @return
     */
    public boolean isStripExistingQ3()
    {
        return stripExistingQ3;
    }

    /**
     * Strip out existing Q3 analysis_results blocks and analysis summaries
     * @return
     */    
    public void setStripExistingQ3(boolean stripExistingQ3)
    {
        this.stripExistingQ3 = stripExistingQ3;
    }

    /**
     * Use XPRESS output format
     *
     * @param mimicXpress If true, will cause output to be written in a format compatible with XPRESS
     */
    public void setMimicXpress(boolean mimicXpress)
    {
        this.mimicXpress = mimicXpress;
    }

    /**
     * Control use of sentinel values
     *
     * @param noSentinels If true, sentinel values will <b>not<b> be used in ratios when heavy area is zero
     */
    public void setNoSentinels(boolean noSentinels)
    {
        this.noSentinels = noSentinels;
    }

    /**
     * Control output of intermediate values for debugging
     *
     * @param debugMode If true, extra debugging values will be written as attributes in
     * each <code>q3ratio</code> block
     */
    public void setDebugMode(boolean debugMode)
    {
        this.debugMode = debugMode;
    }

    /**
     * Control compatibility with original Q3 code
     *
     * @param compatMode If true, replicate behavior of original Q3 code when fewer than 3
     * isotopes where matched in the center scan between heavy and light. This behavior can
     * give a slight boost to the heavy area.
     */
    public void setCompatMode(boolean compatMode)
    {
        this.compatMode = compatMode;
    }

    /**
     * Read a pepXML file and apply q3 to each fraction
     *
     * @param pepXmlFilename Path to the pepXML file to process.
     * @param labeledResidue The single-letter code for the amino acid labeled.
     *        Multiple labeled residues or N-terminal/C-terminal mods are
     *        currently not supported.
     * @param minPeptideProphet The lowest PeptideProphet score to be
     *        considered
     * @param forceOutput Force writing output even if the output file already exists
     * @param useSentinels When heavy area is zero, use XPRESS style sentinel values in ratios
     * @param mimicXpress Use an Xpress analysis result format
     * @return
     * @throws IOExceptiona
     * @throws XMLStreamException
    public static void doQ3(String pepXmlFilename,
                            char labeledResidue,
                            float massDiff,
                            float minPeptideProphet,
                            float maxFracDeltaMass,
                            boolean forceOutput,
                            boolean mimicXpress,
                            boolean noSentinels,
                            String q3XmlFilename)
        throws XMLStreamException, IOException
    {
        Q3 q3 = new Q3(pepXmlFilename, labeledResidue, massDiff, q3XmlFilename);
        q3.setMinPeptideProphet(minPeptideProphet);
        q3.setMaxFracDeltaMass(maxFracDeltaMass, true);
        q3.setForceOutput(forceOutput);
        q3.setMimicXpress(mimicXpress);
        q3.setNoSentinels(noSentinels);
        q3.apply();
    }
     */


    /**
     * Read a pepXML file and apply q3 to each fraction
     *
     * @throws IOException
     * @throws XMLStreamException
     */
    public void apply()
        throws XMLStreamException, IOException
    {
        if (null == pepXmlFilename)
            throw new Q3RuntimeException("No pepXML filename was given");

        if (minPeptideProphet > 1.0)
            throw new Q3RuntimeException("PeptideProphet cutoff must not exceed 1.0");

/* ????
        if (labeledResidue < 'A' || labeledResidue >= 'Z')
            throw new Q3RuntimeException("Invalid labeled residue character '" + labeledResidue + "'");
*/

        if (null == q3XmlFilename)
            q3XmlFilename = getOutputFilename(pepXmlFilename);

        File q3XmlFile = new File(q3XmlFilename);
        if (!forceOutput)
        {
            if (q3XmlFile.exists())
                throw new Q3RuntimeException("Operation would overwrite existing file " + q3XmlFilename);
        }
            
        ApplicationContext.setMessage("Q3 output will be written to " + q3XmlFilename);

        String tmpFilename = q3XmlFilename + ".copy";
        File tmpFile = new File(tmpFilename);

        try
        {
            // process all fractions. Returns a list of lists
            List<List<Q3Peptide>> results = quantitate();

            // write a new pepXML file and add corresponding analysis blocks
            Q3PepXmlRewriter rewriter = new Q3PepXmlRewriter(results, pepXmlFilename, labels, mimicXpress, noSentinels, debugMode, stripExistingQ3, tmpFilename);

            rewriter.rewrite();
            rewriter.close();
            if (!tmpFile.renameTo(q3XmlFile))
            {
                if (!q3XmlFile.exists())
                    throw new Q3RuntimeException("Unable to save output to " + q3XmlFile);
                // Windows does not allow rename if destination exists (e.g. if update is in-place)
                if (!q3XmlFile.delete())
                    throw new Q3RuntimeException("Unable to save output; destination file " + q3XmlFile + " exists and could not be deleted");
                if (!tmpFile.renameTo(q3XmlFile))
                    throw new Q3RuntimeException("Unable to save output to " + q3XmlFile);
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("Error running Q3: " + e.getMessage(), e);
            tmpFile.delete();
        }
    }

    /**
     * Make up an output filename
     */ 
    private static String getOutputFilename(String pepXmlFilename)
    {
        if (null == pepXmlFilename)
            throw new IllegalArgumentException("pepXML file name can not be null");
        if (pepXmlFilename.toLowerCase().endsWith(".pep.xml"))
            return pepXmlFilename.substring(0, pepXmlFilename.length()-8) + outputSuffix;
        else if (pepXmlFilename.toLowerCase().endsWith(".xml"))
            return pepXmlFilename.substring(0, pepXmlFilename.length()-4) + outputSuffix;
        else
            return pepXmlFilename + outputSuffix;
    }

    /**
     * Read a pepXML file and run Q3 on each fraction
     *
     * @return the tag mass difference
     * @throws Q3RuntimeException if different fractions have different tag weights
     */
    private List<List<Q3Peptide>> quantitate()
        throws XMLStreamException, IOException
    {
        File pepXmlFile = new File(pepXmlFilename);
        PepXmlLoader loader = new PepXmlLoader(pepXmlFile, _log);

        ArrayList<List<Q3Peptide>> master = new ArrayList<List<Q3Peptide>>();

        try
        {
            FractionIterator fi = loader.getFractionIterator();
            int fractionId = 0;
            while (fi.hasNext())
            {
                fractionId++;

                PepXmlFraction fraction = (PepXmlFraction) fi.next();

                File mzXmlFile = findMzXmlFile(fraction, pepXmlFile);

                if (null == mzXmlFile)
                    throw new Q3RuntimeException("Could not find mzXML file associated with " + pepXmlFilename + " fraction number " + fractionId);

                // Write feature pairs for processing by Q3
                List<Q3Peptide> fractionPeps = readFraction(fraction, mzXmlFile);

                // Run Q3 on peptides from this fraction
                quantitateFraction(fractionPeps, mzXmlFile.getPath());

                master.add(fractionPeps);
            }
        }
        finally
        {
            loader.close();
        }
        return master;
    }

    /**
     * Try to locate the mzXML file associated with this fraction.
     */
    private File findMzXmlFile(PepXmlFraction fraction, File pepXmlFile)
        throws IOException
    {
        String base = fraction.getDataBasename();

        if (base == null)
        {
            base = pepXmlFile.getName();
            if (base.contains("."))
                base = base.substring(0, base.indexOf("."));
        }
        File f = null;

        // Check any explicitly specified directory first
        if (null != alternateMzXmlDir)
        {
            //dhmay changing 2008/06/25:
            //base can contain a partial or full directory structure.  That structure isn't
            //appropriate if the user has specified a directory, so remove it for this check.
            String altBase = base;
            if (altBase != null && altBase.contains(File.separator) &&
                    altBase.length() > altBase.lastIndexOf(File.separator) + 1)
            {
                _log.debug("Removing filepath from base: " + base);
                altBase = altBase.substring(altBase.lastIndexOf(File.separator) + 1);
            }
            f = checkForMzXml(alternateMzXmlDir, altBase);
            if (null != f)
                return f;
        }

        // Next check the given spectrum path
        String mzXmlFilename = fraction.getSpectrumPath();
        if (null != mzXmlFilename)
        {
            f = new File(mzXmlFilename);

            // If it exists, just return
            if (f.exists())
                return f;
        }

        // Finally, try looking in the same directory as the pepXML file
        String pepXmlDir = pepXmlFile.getParent();
        if (null == pepXmlDir) pepXmlDir = "."; // Use CWD
        f = checkForMzXml(pepXmlDir, base);
        if (null != f)
            return f;

        return null;
    }

    /**
     *
     */
    private static File checkForFile(String dir, String base, String suffix)
    {
        File f = new File(dir, base + suffix);
        _log.debug("Checking for " + f.getPath());
        if (f.exists())
            return f;
        if (base.endsWith(".pep"))
        {
            base = base.substring(0,base.length()-4);
            f = new File(dir, base + suffix);
            _log.debug("Checking for " + f.getPath());
            if (f.exists())
                return f;
        }
        return null;
    }

    /**
     *
     */
    private static File checkForMzXml(String dir, String base)
    {
        _log.debug("checkForMzXml: dir=" + dir + ", base=" + base);        
        File f = checkForFile(dir, base, ".mzXML");
        if (null != f)
            return f;
        f = checkForFile(dir, base, ".mzxml");
        if (null != f)
            return f;
        return null;
    }

    /**
     * Try to locate the mzXML file associated with this fraction. First
     * look for the full path recorded in pepXML, then in the same directory
     * as the pepXML file.
     */
    private static File findMzXmlFileOld(PepXmlFraction fraction, File pepXmlFile)
        throws IOException
    {
        String mzXmlFilename = fraction.getSpectrumPath();
        if (null != mzXmlFilename)
        {
            File f = new File(mzXmlFilename);

            // If it exists, just return
            if (f.exists())
                return f;

            // Otherwise look in the same directory as the pepXML file
            if (null != pepXmlFile)
            {
                f = new File(pepXmlFile.getParentFile(), f.getName());
                if (f.exists())
                    return f;
            }
        }

        // If not yet found, look for mzXML file with same base name
        // as input pepXML; this can only work for single-fraction files
        String pepXmlFilename = pepXmlFile.getCanonicalPath();
        if (pepXmlFilename.toLowerCase().endsWith(".pep.xml"))
            mzXmlFilename = pepXmlFilename.substring(0, pepXmlFilename.length() - 8);
        else if (pepXmlFilename.toLowerCase().endsWith(".xml"))
            mzXmlFilename = pepXmlFilename.substring(0, pepXmlFilename.length() - 4);
        else
            return null;

        File f = new File(mzXmlFilename + ".mzXML");
        if (f.exists())
            return f;

        return null;
    }

    /**
     * Given the modifications defined for a fraction, figure out
     * the light and heavy masses that we'll see for a labeled residue
     */
    private void calculateLabelMasses(PepXmlFraction fraction)
    {
        // Process amino acid mods
        List<MS2Modification> mods = fraction.getModifications();
        for (MS2Modification mod : mods)
        {
            char aa = mod.getAminoAcid().charAt(0);
            IsotopicLabel label = labels.get(aa);
            if (null != label)
            {
                float mass = mod.getMass();
                if (mod.getVariable())
                    label.setHeavyMass(mass);
                else
                    label.setLightMass(mass);
            }
        }

/**********************************************************************
   ???? Pending commit of PepXmlLoader updates...

        // Process terminal mods
        List<MS2TerminalModification> tmods = fraction.getTerminalModifications();
        for (MS2TerminalModification mod : tmods)
        {
            char aa = mod.isNTerminus() ? 'n' : 'c'; // use psuedo amino acids for termini
            IsotopicLabel label = labels.get(aa);
            if (null != label)
            {
                float mass = mod.getMass();
                if (mod.getVariable())
                    label.setHeavyMass(mass);
                else
                    label.setLightMass(mass);
            }
        }
**********************************************************************/

        // Make sure we've got a complete table of masses for the labels
        for (IsotopicLabel label : labels.values())
        {
            float light = label.getLightMass();
            float heavy = label.getHeavyMass();
            float diff = label.getMassDiff();

            _log.debug("Label Before = " + label.toString());
            
            if (heavy == 0f && light == 0f)
                throw new Q3RuntimeException("No static or variable modification found for labeled residue '" + label.getResidue() + "'");
            if (heavy == 0f && light != 0f)
                label.setHeavyMass(light + diff);
            else if (light == 0f && heavy != 0f)
                label.setLightMass(heavy - diff);

            _log.debug("Label After = " + label.toString());

            if (label.getLightMass() < 0f ||
                label.getHeavyMass() < 0f ||
                label.getHeavyMass() < label.getLightMass())
            {
                throw new Q3RuntimeException("Inconsistent light and heavy masses (" + label.getLightMass() + " and " + label.getHeavyMass() + ")");
            }

            if (Math.abs(label.getHeavyMass() - label.getLightMass() - diff) > 0.01)
                _log.warn("Specified quantitation mass diff, " + diff + ", may not be compatible with mass diff used in search " + (heavy - light));
        }
    }

    /**
     * Read one fraction of the pepXML file and return a list of Q3Peptides
     */
    private List<Q3Peptide> readFraction(PepXmlFraction fraction, File mzXmlFile)
    {
        calculateLabelMasses(fraction);

        ArrayList<Q3Peptide> results = new ArrayList<Q3Peptide>();

        PeptideIterator pi = fraction.getPeptideIterator();
        int peptideId = 0;
        while (pi.hasNext())
        {
            peptideId++; // Count each feature within a fraction so we can join up results later..
            PepXmlPeptide peptide = (PepXmlPeptide) pi.next();

            // Filter by PeptideProphet probability, if requested and available
            if (filterByMinPeptideProphet)
            {
                PeptideProphetHandler.PeptideProphetResult pepProphet = peptide.getPeptideProphetResult();
                if (null != pepProphet && pepProphet.getProbability() < minPeptideProphet)
                    continue;
            }

            // Filter by fractional delta mass
            if (filterByMaxFracDeltaMass)
            {
                double deltaMass = peptide.getDeltaMass();
                double fracDeltaMass = Math.abs(deltaMass - Math.round(deltaMass));
                if (maxFracDeltaMassIsPPM)
                    fracDeltaMass = 1000000f * fracDeltaMass / peptide.getCalculatedNeutralMass();
                if (fracDeltaMass > maxFracDeltaMass)
                    continue;
            }

            char[] trimmed = peptide.getTrimmedPeptide().toCharArray();
            ModifiedAminoAcid[] mods = peptide.getModifiedAminoAcids();

            float totalMassDiff = 0f;
            int countPotential = 0;
            int countLight = 0;
            int countHeavy = 0;

/**********************************************************************
   ???? Pending commit of PepXmlLoader updates...
            char[] termini = {'n', 'c'};

            // Look for terminal labeling first
            for (char terminus : termini)
            {
                IsotopicLabel label = labels.get(terminus);
                if (null != label)
                {
                    countPotential++;
                    totalMassDiff += label.getMassDiff();

                    float mass = (terminus == 'n' ?
                                  peptide.getNtermModificationMass() :
                                  peptide.getCtermModificationMass());
                    if (mass == 0f || closeEnough(mass, label.getLightMass()))
                        countLight++;
                    else if (closeEnough(mass, label.getHeavyMass()))
                        countHeavy++;
                }
            }
**********************************************************************/

            for (int i = 0; i < trimmed.length; i++)
            {
                IsotopicLabel label = labels.get(trimmed[i]);
                if (null != label)
                {
                    countPotential++;
                    totalMassDiff += label.getMassDiff();

                    // Unmodified is assumed light (e.g. 12C lysine in SILAC)
                    if (null == mods || null == mods[i])
                    {
                        countLight++;
                    }
                    else
                    {
                        float mass = (float) mods[i].getMass();
                        if (mass == 0f || closeEnough(mass, label.getLightMass()))
                            countLight++;
                        else if (closeEnough(mass, label.getHeavyMass()))
                            countHeavy++;
                    }
                }

            }

            // Output only peptides that contain at least one label;
            // skip partially labeled peptides
            if (countPotential > 0)
            {
                String type = "";
                if (countLight < countPotential && countHeavy < countPotential)
                    type = "Partial";
                else
                    type = (countHeavy == countPotential ? "Heavy" : "Light");

                _log.debug(peptide.getPeptide() + " : " + type + " " + totalMassDiff + "Da " + countPotential + " " + countLight + " " + countHeavy);
                if (countLight < countPotential && countHeavy < countPotential)
                {
                    _log.warn("Skipping partially labeled peptide: " + peptide.getPeptide());
                }
                else
                {
                    boolean isHeavy = (countHeavy == countPotential);
                    results.add(new Q3Peptide(peptideId, peptide, totalMassDiff, isHeavy));
                }
            }
        }
        return results;
    }

    /**
     *
     */
    private boolean closeEnough(float m1, float m2)
    {
        return Math.abs(m1 - m2) < MASS_MATCH_THRESHOLD;
    }

    /**
     *
     */
    void quantitateFraction(List<Q3Peptide> peptides, String mzXmlFilename)
        throws IOException
    {
        double cutoff = PEAK_MATCH_THRESHOLD;

        MSRun.setShowIndexBuilderProgress(false);
        MSRun run = MSRun.load(mzXmlFilename, false); // Do NOT save a .inspect file for Q3

        // loadf
        // setup0

        _log.info("Processing peptides from " + mzXmlFilename);

        for (Q3Peptide p : peptides)
        {
//            StringBuilder sb = new StringBuilder();

//            _log.info("Got " + p.toString());

            // getblock2
            int safeiso = (int) Math.round((p.mzH - p.mzL) * p.charge) - 1;
            int safeiso5 = safeiso < 5 ? safeiso : 5;

            _log.debug("peptide = " + p.toString());
            _log.debug("mzL = " + p.mzL + " mzH = " + p.mzH);
            _log.debug("safeiso5 = " + safeiso5 + " (" + safeiso + ")");

            int[] L = getScanIndexList(run, p.scan);
            double[] targ = new double[safeiso5 + 1 + 5 + 1];
            int i = 0;
            for (int j = 0; j <= safeiso5; j++)
                targ[i++] = p.mzL + j * 1.f / p.charge;
            for (int j = 0; j <= 5; j++)
                targ[i++] = p.mzH + j * 1.f / p.charge;
            double[][] D = new double[L.length][targ.length];
            for (int iL = 0; iL < L.length; iL++)
            {
                double[][] peaks = getSpectrum(run, L[iL]);

                // get [mzL - 1 to mzH + 3] into x and y
                int start = fixIndex(Arrays.binarySearch(peaks[0], p.mzL - 1));
                int end = fixIndex(Arrays.binarySearch(peaks[0], p.mzH + 3));

                if (end > peaks[0].length)
                    end = peaks[0].length - 1;
                if (start >= end)
                    continue;

                double[] x = new double[end - start];
                double[] y = new double[end - start];

                System.arraycopy(peaks[0], start, x, 0, end - start);
                System.arraycopy(peaks[1], start, y, 0, end - start);

                int[] mi = localmaxes(y);

                double[] xx = new double[mi.length];
                double[] yy = new double[mi.length];

                for (int j = 0; j < mi.length; j++)
                {
                    xx[j] = x[mi[j]];
                    yy[j] = y[mi[j]];
                }

                // b=bindtarg(xx,targ,cutoff*mzH[ii])
                int r[][] = bindtarg(xx, targ, cutoff * p.mzH);

                // D[iL,b$targi]=yy[b$xxi]
                for (int j = 0; j < r[0].length; j++)
                    D[iL][r[1][j]] = yy[r[0][j]];
            }

            double[] isoL = getiso(p.lightMass, targ.length);
            // pad isoH with safeiso5 leading zeroes
            // ???? Check this; should it be safeiso5+1?
            double[] tmp = getiso(p.heavyMass, targ.length - safeiso5);
            double[] isoH  = new double[targ.length];
            for (int j = 0; j < safeiso5; j++)
                isoH[j] = 0;
            for (int j = safeiso5; j < targ.length; j++)
                isoH[j] = tmp[j - safeiso5];

            // Lhit=(D[,1:(safeiso5+1)]>0)
            // Rhit=(D[,safeiso5+1+(1:(safeiso5+1))]>0)
            // match=Lhit&Rhit
            boolean[][] match = getMatches(D, safeiso5);

            // matchD=D*cbind(match,match,array(0,c(nrow(D),ncol(D)-2*ncol(match))))
            // we just don't bother creating the entries that get zeroed out in R

            double[][] matchD = new double[L.length][(safeiso5+1)*2];
            for (int scan = 0; scan < L.length; scan++)
            {
                for (int iso = 0; iso <= safeiso5; iso++)
                    if (match[scan][iso])
                    {
                        matchD[scan][iso] = D[scan][iso];
                        matchD[scan][iso + safeiso5 + 1] = D[scan][iso + safeiso5 + 1];
                    }
                    else
                    {
                        matchD[scan][iso] = 0;
                        matchD[scan][iso + safeiso5 + 1] = 0;
                    }
            }

            // matchct=rowSums(match)
            int[] matchct = new int[L.length];
            for (int scan = 0; scan < L.length; scan++)
            {
                matchct[scan] = 0;
                for (int iso = 0; iso <= safeiso5; iso++)
                    if (match[scan][iso])
                        matchct[scan] += 1;
            }

            // Walk from the center scan in either direction to find extent
            int center = L.length > 9 ? 9 : L.length - 1;

            int left = center;
            while (left >= 1)
            {
                if (matchct[left-1] <= 0)
                    break;
                left--;
            }

            int right = center;
            while (right < L.length - 1)
            {
                if (matchct[right+1] <= 0)
                    break;
                right++;
            }

            p.firstScan = run.getScan(L[left]).getNum();
            p.lastScan = run.getScan(L[right]).getNum();
            p.centerMatchct = matchct[center];

            // treat 10 specially as a failsafe in case there aren't a reasonable number of matches available
            // xxx=matchD[setdiff(left:right,10),] 
            // #if 10 has too few matches, use all of it; otherwise just use the matches
            // if (matchct[10]<=2)
            //   xxx=rbind(D[10,],xxx)
            // else
            //   xxx=rbind(matchD[10,],xxx)

            if (matchct[center]<=2)
            {
                for (int j = 0; j < matchD[center].length; j++)
                    matchD[center][j] = D[center][j];
            }

            double[] xxs = new double[matchD[0].length];
            for (int j = 0; j < matchD[0].length; j++)
            {
                xxs[j] = 0;
                for (int k = left; k <= right; k++)
                    xxs[j] += matchD[k][j];
            }

            double[] q2 = new double[2];
            q2[0] = 0;
            for (int j = 0; j <= safeiso5; j++)
                q2[0] += xxs[j];
            q2[1] = 0;
            for (int j = safeiso5+1; j < xxs.length; j++)
                q2[1] += xxs[j];

            if (compatMode && matchct[center] <= 2)
            {
                // Replicate the behavior of the R code when only two matches are available
                for (int k = xxs.length; k < D[center].length; k++)
                    q2[1] += D[center][k];
            }

            p.q2L = q2[0];
            p.q2H = q2[1];

            // # Make a naive correction for the overlap
            // shift=safeiso[ii]+1
            // head=sum(dpois(0:(shift-1),mass[ii]/1800))
            // tail=sum(dpois(shift:(shift+shift-1),mass[ii]/1800))
            // Xht=rbind(c(head,0),c(tail,head))
            int shift = safeiso + 1;
            double[] t2 = getiso(p.heavyMass, shift);

            double head = 0;
            for (int j = 0; j < shift; j++)
                head += t2[j];

            t2 = getiso(p.heavyMass, shift+shift);

            double tail = 0;
            for (int j = shift; j < shift + shift; j++)
                tail += t2[j];
            double[][] Xht = new double[2][2];
            Xht[0][0] = head;
            Xht[0][1] = 0;
            Xht[1][0] = tail;
            Xht[1][1] = head;

            Matrix A = new Matrix(Xht);
            Matrix B = new Matrix(q2,2);
            Matrix R = A.solve(B);
            double[][] r = R.getArray();

            p.q3L = r[0][0];
            p.q3H = r[1][0] >= 0.0 ? r[1][0] : 0.0; // Clamp corrected heavy areas to zero
        }
    }


    /**
     * Build a boolean array indicating, for each peak, whether a match was found between light and heavy
     */
    static boolean[][] getMatches(double[][] D, int safeiso5)
    {
        boolean[][] match = new boolean[D.length][safeiso5+1];
        for (int i = 0; i < D.length; i++)
            for (int j = 0; j < match[0].length; j++)
                match[i][j] = D[i][j] > 0 && D[i][j + safeiso5+1] > 0;
        return match;
    }

    /**
     * Build a match list ignoring the heavy species; used for "unlabeled" quant only. Hacktastic
    static boolean[][] fakeMatches(double[][] D, int safeiso5)
    {
        boolean[][] match = new boolean[D.length][safeiso5+1];
        for (int i = 0; i < D.length; i++)
            for (int j = 0; j < match[0].length; j++)
                match[i][j] = D[i][j] > 0;
        return match;
    }
     */

    static double[] fact = new double[]{1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0, 1.307674e+12, 2.092279e+13, 3.556874e+14, 6.402374e+15, 1.216451e+17, 2.432902e+18, 5.109094e+19, 1.124001e+21, 2.585202e+22, 6.204484e+23};

    /**
     * Get isotopes based on a simple poisson model
     *              p(x) = lambda^x exp(-lambda)/x!
     */
    static double[] getiso(double mass, int len)
    {
        double lambda = mass / 1800;
        double lambdaExp = Math.exp(-lambda);

        double[] r = new double[len];

        for (int i = 0; i < len; i++)
            if (i < fact.length)
                r[i] = (double) (Math.pow(lambda, i) / fact[i] * lambdaExp);
            else
                r[i] = 0; // Not really zero, but too small to matter

        return r;
    }

    /**
     *
     */
    public static void logArray(String desc, double[][] a)
    {
        StringBuilder sb = new StringBuilder();
        _log.debug("**** " + desc);
        for (int outer = 0; outer < a[0].length; outer += 7)
        {
            sb.setLength(0);
            for (int j = 0; j < 7; j++)
                sb.append("\t[," + (outer + j) + "]");
            _log.debug(sb.toString());

            for (int scan = 0; scan < a.length; scan++)
            {
                sb.setLength(0);
                sb.append(scan + ":\t");
                for (int iso = outer; iso < a[0].length && iso < outer + 7; iso++)
                    sb.append(a[scan][iso] + "\t");
                _log.debug(sb.toString());
            }
        }

    }

    /**
     * Get peaks for one scan and cast to a double array
     */
    static double[][] getSpectrum(MSRun run, int scan)
    {
        MSScan s = run.getScan(scan);
        float[][] t1 = s.getSpectrum();
        double[][] peaks = new double[2][t1[0].length];
        for (int j = 0; j < t1[0].length; j++)
        {
            peaks[0][j] = t1[0][j];
            peaks[1][j] = t1[1][j];
        }
        return peaks;
    }

    /**
     * A double tagged with an integer indicating which array it orginally came from
     * as well as the index relative to that original array
     */
    static class TaggedDouble
    {
        static Comparator<TaggedDouble> compareByDoubleAsc = new Comparator<TaggedDouble>()
        {
            public int compare(TaggedDouble a, TaggedDouble b)
            {
                return a.f < b.f ? -1 : (a.f == b.f ? 0 : 1);
            }
        };

        double f;
        int tag;
        int index;

        TaggedDouble(double f, int tag, int index)
        {
            this.f = f;
            this.tag = tag;
            this.index = index;
        }

        public String toString()
        {
            return "[" + f + ", " + tag + ", " + index + "]";
        }
    }

    /**
     * for each element of targ, find the closest element of xx and if it's within
     * distance tol record the match; return the xx and targ indices of all such matches
     */
    static int[][] bindtarg(double[] xx, double[] targ, double tol)
    {
        int[][] r = new int[2][];

        // Combine both arrays, preserving source with tags (1 == xx, 2 == targ), and sort
        TaggedDouble[] b = new TaggedDouble[xx.length + targ.length];
        for (int i = 0; i < xx.length; i++)
            b[i] = new TaggedDouble(xx[i], 1, i);
        for (int i = 0; i < targ.length; i++)
            b[i + xx.length] = new TaggedDouble(targ[i], 2, i);
        Arrays.sort(b, TaggedDouble.compareByDoubleAsc);

        // dbo = diff(bo)
        // hit=dbo<tol
        // ch=bido[1:(length(bido)-1)]!=bido[2:length(bido)]
        // hit=hit&ch
        double[] dbo = new double[b.length];
        boolean[] hit = new boolean[b.length];
        boolean[] ch = new boolean[b.length];
        for (int i = 1; i < b.length; i++)
        {
            dbo[i] = b[i].f - b[i - 1].f;
            ch[i] = b[i].tag != b[i - 1].tag;
            hit[i] = dbo[i] < tol && ch[i];
        }

        // rhit=c(hit[2:length(hit)],F)
        // rbet=c(dbo[2:length(dbo)]<dbo[1:(length(dbo)-1)],F)
        boolean[] rhit = new boolean[b.length];
        boolean[] rbet = new boolean[b.length];
        for (int i = 0; i < b.length - 1; i++)
        {
            rhit[i] = hit[i+1];
            rbet[i] = dbo[i+1]<dbo[i];
        }
        rhit[b.length-1] = false;
        rbet[b.length-1] = false;

        // lhit=c(F,hit[1:(length(hit)-1)])
        // lbet=c(F,dbo[2:length(dbo)]>dbo[1:(length(dbo)-1)])
        boolean[] lhit = new boolean[b.length];
        boolean[] lbet = new boolean[b.length];
        lhit[0] = false;
        lbet[0] = false;
        for (int i = 1; i < b.length; i++)
        {
            lhit[i] = hit[i-1];
            lbet[i] = dbo[i]>dbo[i-1];
        }

        // hit=hit&(!(rhit&rbet))&(!(lhit&lbet))
        int count = 0;
        int xxcount = 0;
        int targcount = 0;
        for (int i = 0; i < b.length; i++)
        {
            hit[i] = hit[i] && (!(rhit[i]&&rbet[i]) && !(lhit[i]&&lbet[i]));
            if (hit[i])
            {
                count++;
                if (1 == b[i].tag)
                    xxcount++;
                if (2 == b[i].tag)
                    targcount++;
            }
        }

        _log.debug("Found " + count + " hits (" + xxcount + " " + targcount + ")");

        r[0] = new int[count];
        r[1] = new int[count];
        int xc = 0;
        int tc = 0;
        for (int i = 0; i < b.length; i++)
            if (hit[i])
                if (1 == b[i].tag)
                {
                    r[0][xc++] = b[i].index;
                    r[1][tc++] = b[i-1].index;
                }
                else if (2 == b[i].tag)
                {
                    r[1][tc++] = b[i].index;
                    r[0][xc++] = b[i - 1].index;
                }

        return r;
    }

    /**
     *
     */
    static double[] diff(double[] x)
    {
        double[] r = new double[x.length - 1];
        for (int i = 0; i < r.length; i++)
            r[i] = x[i + 1] - x[i];
        return r;
    }

    /**
     *
     */
    static int[] localmaxes(double[] x)
    {
        double epsilon = 1E-10;
        int n = x.length;
        if (n < 2)
            return new int[0];
        double[] dx = diff(x);

        //  ismax=c(x[1]>x[2],dx[-(n-1)]>epsilon & dx[-1]< -epsilon, x[n]>x[n-1])
        boolean[] ismax = new boolean[n];
        ismax[0] = x[0] > x[1];
        int maxcount = 0;
        for (int i = 1; i < n - 1; i++)
        {
            ismax[i] = dx[i-1] > epsilon && dx[i] < -epsilon;
            if (ismax[i])
                maxcount++;
        }
        ismax[n-1] = x[n - 1] > x[n - 2];
        if (ismax[0])
            maxcount++;
        if (ismax[n-1])
            maxcount++;

        //(1:n)[ismax]
        int[] maxes = new int[maxcount];
        int mi = 0;
        for (int i = 0; i < n; i++)
            if (ismax[i])
                maxes[mi++] = i;

        return maxes;
    }

    /**
     * Convert negative indices ("insertion point") to valid index
     */
    static int fixIndex(int i)
    {
        return i >= 0 ? i : -(i + 1);
    }

    /**
     * Returns a list of MS1 scan *indicies*!
     */
    static int[] getScanIndexList(MSRun run, int scan)
    {
        int scanIndex = run.getIndexForScanNum(scan, true);
        int leftIndex = scanIndex - 10;
        int rightIndex = scanIndex + 20;
        if (leftIndex < 0)
            leftIndex = 0;
        if (rightIndex >= run.getScanCount())
            rightIndex = run.getScanCount() - 1;
        int[] list = new int[rightIndex - leftIndex + 1];
        for (int i = leftIndex; i <= rightIndex; i++)
            list[i-leftIndex] = i;
        Arrays.sort(list);
        return list;
    }

    /**
     * 
     */
    static class Q3PepXmlRewriter extends SimpleXMLEventRewriter
    {
        Map<Character,IsotopicLabel> labels = null;
        Q3Peptide[] q3peptides = null;
        GregorianCalendar now = null;
        boolean mimicXpress;
        boolean debugMode;
        int peptideId;
        double sentinelInf;
        double sentinelNaN;
        boolean stripExistingQ3;

        // Match XPRESS convention when reporting masses
        boolean reportSinglyProtonatedMasses = true;

        //keep track of whether we're in an analysis_result (for stripExistingAnalysisResults)
        boolean insideOldQ3AnalysisResult = false;
        boolean insideOldQ3AnalysisSummary = false;

        List<List<Q3Peptide>> results;
        List<Q3Peptide> currentFraction;
        Q3Peptide currentPeptide;
        
        public Q3PepXmlRewriter(List<List<Q3Peptide>> results, String inFilename, Map<Character,IsotopicLabel> labels, boolean mimicXpress, boolean noSentinels, boolean debugMode, boolean stripExistingQ3, String outFilename)

        {
            super(inFilename, outFilename);
            this.results = results;
            this.labels = labels;
            this.mimicXpress = mimicXpress;
            this.stripExistingQ3 = stripExistingQ3;
            peptideId = 0;
            now = new GregorianCalendar();
            if (noSentinels)
            {
                // Use standard IEEE special values for ill-defined ratios
                sentinelNaN = Double.NaN;
                sentinelInf = Double.POSITIVE_INFINITY;
            }
            else
            {
                // Use sentinel values for ill-defined ratios
                sentinelNaN = (double) RelativeQuantAnalysisResult.SENTINEL_NAN;
                sentinelInf = (double) RelativeQuantAnalysisResult.SENTINEL_POSITIVE_INFINITY;
            }
        }

        /**
         * overrides parent add() method.  Only add events if we're not inside an old analysis result
         * @param event
         * @throws XMLStreamException
         */
        public void add(XMLEvent event)
                throws XMLStreamException
        {
            if (!stripExistingQ3 || (!insideOldQ3AnalysisResult && !insideOldQ3AnalysisSummary))
            {
//if (event.isStartElement()) System.err.println("adding Start " + event.asStartElement().getName().getLocalPart());
//else if (event.isEndElement()) System.err.println("adding End " + event.asEndElement().getName().getLocalPart());
//else System.err.println("*********adding unknown event " + event.getEventType() + ", compare " +     XMLStreamConstants.START_ELEMENT +","+ XMLStreamConstants.END_ELEMENT+","+ XMLStreamConstants.CHARACTERS+","+ XMLStreamConstants.ATTRIBUTE+","+ XMLStreamConstants.NAMESPACE+","+ XMLStreamConstants.PROCESSING_INSTRUCTION+","+ XMLStreamConstants.COMMENT+","+XMLStreamConstants.START_DOCUMENT+","+ XMLStreamConstants.END_DOCUMENT+","+ XMLStreamConstants.DTD);

                //dhmay adding a workaround here to avoid a nullpointerexception from super.add(),
                //specifically doWriteDefaultNs.  Namespace can't be null, or wstx falls apart.  This is apparently
                //addressed in version 3.2.3, so we can take out this hack when we upgrade.

                /**** mfitzgib 090326 -
                      This work-around emits xmlns="" namespaces on each
                      start tag and doesn't seem to be needed even with
                      Woodstox 3.2.1; delete this block when we upgrade
                      to 3.2.8 and are feeling comparitively comfortable.
                if (event.isStartElement())
                {
                    SimpleStartElement newStart = new SimpleStartElement(event.asStartElement().getName().getLocalPart());

                    // getAttributes returns a non-genericized Iterator; OK
                    @SuppressWarnings("unchecked")
                    Iterator<Attribute> attIter =  event.asStartElement().getAttributes();

                    boolean alreadyHasNS = false;
                    while (attIter.hasNext())
                    {
                        Attribute attr = attIter.next();
                        if (attr.getName().getLocalPart().equals("xmlns"))
                        {
                            alreadyHasNS = true;
                            break;
                        }
                        newStart.addAttribute(attr.getName().getLocalPart(), attr.getValue());
                    }
                    if (!alreadyHasNS)
                    {
                        newStart.addAttribute("xmlns","");
                        event = newStart.getEvent();
                    }
                }
                ****/
                super.add(event);
            }
            else
            {
//                System.err.println("SKIPPING " + event.getEventType());
//                if (event.isStartElement())
//                    System.err.println("         (" + event.asStartElement().getName() + ")");
            }
        }

        public void handleDefault(XMLEvent event)
            throws XMLStreamException
        {
            if (!stripExistingQ3 || (!insideOldQ3AnalysisResult && !insideOldQ3AnalysisSummary))
            {
                super.handleDefault(event);
            }
        }

        public void handleStartElement(StartElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();

            if ("msms_run_summary".equals(qname.getLocalPart()))
            {
                currentFraction = (results.size() > 0 ? results.remove(0) : null);
                currentPeptide = nextQ3Peptide();
                peptideId = 0;
                add(event);
            }
            else if ("msms_pipeline_analysis".equals(qname.getLocalPart()))
            {
                add(event);
                addAnalysisSummary();
            }
            else if ("analysis_result".equals(qname.getLocalPart()))
            {
                Attribute analysisAttribute = event.getAttributeByName(new QName("analysis"));
                if (analysisAttribute != null && "q3".equals(analysisAttribute.getValue()))
                {
                    if (stripExistingQ3)
                        _log.debug("Found existing analysis result, skipping...");
                    insideOldQ3AnalysisResult = true;
                }
                add(event);
            }
            else if ("analysis_summary".equals(qname.getLocalPart()))
            {
                Attribute analysisAttribute = event.getAttributeByName(new QName("analysis"));
                if (analysisAttribute != null && "q3".equals(analysisAttribute.getValue()))
                {
                    if (stripExistingQ3)
                        _log.debug("Found existing Q3 analysis summary, skipping...");
                    insideOldQ3AnalysisSummary = true;
                }
                add(event);
            }
            else
            {
//System.err.println(qname.getLocalPart());
                add(event);
            }
        }

        public void handleEndElement(EndElement event)
            throws XMLStreamException
        {
            QName qname = event.getName();
            if ("search_summary".equals(qname.getLocalPart()))
            {
                add(event);
                addAnalysisTimestamp();
            }
            else if ("search_hit".equals(qname.getLocalPart()))
            {
                peptideId++;
                // if a corresponding q3 peptide is loaded, output an analysis_result block
                if (null != currentPeptide && peptideId == currentPeptide.peptideId)
                {
                    addAnalysisResult(currentPeptide);
                    currentPeptide = nextQ3Peptide();
                }
                add(event);
            }
            else if ("analysis_result".equals(qname.getLocalPart()))
            {
                add(event);
                if (insideOldQ3AnalysisResult && stripExistingQ3)
                    _log.debug("End of existing Q3 analysis result, resuming.");
                insideOldQ3AnalysisResult = false;
            }
            else if ("analysis_summary".equals(qname.getLocalPart()))
            {
                add(event);
                if (insideOldQ3AnalysisSummary && stripExistingQ3)
                    _log.debug("End of existing Q3 analysis summary, resuming.");
                insideOldQ3AnalysisSummary = false;
            }
            else
            {
                add(event);
            }
        }

        /**
         *
         */
        private Q3Peptide nextQ3Peptide()
        {
            if (null == currentFraction || currentFraction.size() <= 0)
                return null;
            return currentFraction.remove(0);
        }

        /**
         * Add a series of events for the analysis_summary block
         */
        private void addAnalysisSummary()
            throws XMLStreamException
        {
            SimpleStartElement start;
            SimpleEndElement end;

            addNewline();

            start = new SimpleStartElement("analysis_summary");
            start.addAttribute("analysis", mimicXpress ? "xpress" : "q3");
            start.addAttribute("time", now);
//            start.addAttribute("xmlns","");
            add(start);

            addNewline();

            String tagName = mimicXpress ? "xpressratio_summary" : "q3ratio_summary";
            start = new SimpleStartElement(tagName);
            start.addAttribute("version", mimicXpress ? "0.0" : Q3.VERSION);
            start.addAttribute("author", Q3.AUTHOR);

            // ???? Must ensure format of multiples is compatible with TPP
            StringBuilder residues = new StringBuilder();
            StringBuilder specs = new StringBuilder();
            boolean first = true;
            for (IsotopicLabel label : labels.values())
            {
                residues.append(label.getResidue());
                if (first)
                    first = false;
                else
                    specs.append(' ');
                specs.append(label.getResidue());
                specs.append(',');
                specs.append(String.format("%f", label.getMassDiff()));
            }

            start.addAttribute("labeled_residues", residues.toString());
            start.addAttribute("massdiff", specs.toString());

            start.addAttribute("massTol", mimicXpress ? "0.75" : ".1");

            if (mimicXpress)
            {
                start.addAttribute("same_scan_range", "Y");
                start.addAttribute("xpress_light", "0");
            }

            add(start);

            end = new SimpleEndElement(tagName);
            add(end);

            addNewline();

            end = new SimpleEndElement("analysis_summary");
            add(end);
        }

        /**
         * Add the analysis_timestamp block
         */
        private void addAnalysisTimestamp()
            throws XMLStreamException
        {
            SimpleStartElement start;
            SimpleEndElement end;

            addNewline();

            start = new SimpleStartElement("analysis_timestamp");
            start.addAttribute("analysis", mimicXpress ? "xpress" : "q3");
            start.addAttribute("time", now);
            start.addAttribute("id", "1");
            add(start);

            addNewline();

            if (mimicXpress)
            {
                start = new SimpleStartElement("xpressratio_timestamp");
                start.addAttribute("xpress_light", "0");
                add(start);

                end = new SimpleEndElement("xpressratio_timestamp");
                add(end);

                addNewline();
            }

            end = new SimpleEndElement("analysis_timestamp");
            add(end);
        }

        /**
         *
         */
        private void addAnalysisResult(Q3Peptide p)
            throws XMLStreamException
        {
            SimpleStartElement start;
            SimpleEndElement end;

            start = new SimpleStartElement("analysis_result");
            start.addAttribute("analysis", mimicXpress ? "xpress" : "q3");
            add(start);

            addNewline();

            double lightMass = p.lightMass + (reportSinglyProtonatedMasses ? PROTON_MASS : 0);
            double heavyMass = p.heavyMass + (reportSinglyProtonatedMasses ? PROTON_MASS : 0);

            start = new SimpleStartElement(mimicXpress ? "xpressratio_result" : "q3ratio_result");
            start.addAttribute("light_firstscan", p.firstScan);
            start.addAttribute("light_lastscan", p.lastScan);
            start.addAttribute("light_mass", lightMass);
            start.addAttribute("heavy_firstscan", p.firstScan);
            start.addAttribute("heavy_lastscan", p.lastScan);
            start.addAttribute("heavy_mass", heavyMass);
            start.addAttribute("light_area", p.q3L);
            start.addAttribute("heavy_area", p.q3H);
            if (mimicXpress)
            {
                start.addAttribute("ratio", getLight2HeavyString(p.q3L, p.q3H));
                start.addAttribute("heavy2light_ratio", getHeavy2LightString(p.q3L, p.q3H));
            }
            else
            {
                start.addAttribute("q2_light_area", p.q2L);
                start.addAttribute("q2_heavy_area", p.q2H);
                if (debugMode)
                    start.addAttribute("center_match_count", p.centerMatchct);
            }
            start.addAttribute("decimal_ratio", getLight2HeavyRatio(p.q3L, p.q3H));
            add(start);

            end = new SimpleEndElement(mimicXpress ? "xpressratio_result" : "q3ratio_result");
            add(end);

            addNewline();

            end = new SimpleEndElement("analysis_result");
            add(end);

            addNewline();
        }

        /**
         * 
         */ 
        private double getLight2HeavyRatio(double light, double heavy)
        {
            if (heavy != 0.0)
                return light/heavy;

            return (light == 0.0) ? sentinelNaN : sentinelInf;
        }

        /**
         * To match XPRESS default, this ratio is normalized so that the larger
         * area is 1.
         */ 
        private String getLight2HeavyString(double light, double heavy)
        {
            if (heavy == 0.0)
                return (light == 0.0) ? Double.toString(sentinelNaN) : Double.toString(sentinelInf);

            if (light < heavy)
                return String.format("%f", light/heavy) + ":1";
            return "1:" + String.format("%f", heavy/light);
        }

        /**
         * To match XPRESS defaul, this ratio is normalized so that the smaller
         * area is 1.
         */ 
        private String getHeavy2LightString(double light, double heavy)
        {
            if (light == 0.0)
                return (heavy == 0.0) ? Double.toString(sentinelNaN) : Double.toString(sentinelInf);

            if (light < heavy)
                return String.format("%f", heavy/light) + ":1";
            return "1:" + String.format("%f", light/heavy);
        }

    }

    /**
     *
     */
    static class Q3Peptide
    {
        int peptideId;
        int scan;
        int charge;
        double lightMass;
        double heavyMass;
        double mzL;
        double mzH;
        String strippedPeptide;
        String protein;

        // Calculated by Q3 ...
        double q3L;
        double q3H;
        int firstScan;
        int lastScan;

        // These are only output in debugging mode
        int centerMatchct;
        double q2L;
        double q2H;

        /**
         *
         */
        public Q3Peptide(int peptideId, PepXmlPeptide peptide, float massDiff, boolean isHeavy)
        {
            this.peptideId = peptideId;
            this.strippedPeptide = peptide.getTrimmedPeptide();
            this.scan = peptide.getScan();
            this.protein = peptide.getProtein();
            this.charge = peptide.getCharge();

            double mass = peptide.getCalculatedNeutralMass();

            if (isHeavy)
            {
                heavyMass = mass;
                lightMass = mass - massDiff;
            }
            else
            {
                lightMass = mass;
                heavyMass = mass + massDiff;
            }

            this.mzL = lightMass / this.charge + PROTON_MASS;
            this.mzH = heavyMass / this.charge + PROTON_MASS;
        }

        /**
         *
         */
        public String toString()
        {
            return "[" + peptideId + ", " +  scan + ", " + charge + "+, " + lightMass + "(" + mzL + "), " + heavyMass + "(" + mzH + "), " + strippedPeptide + "]";
        }

    }

    /**
     *
     */
    static Comparator<Q3Peptide> compareQ3PeptideIdAsc = new Comparator<Q3Peptide>()
    {
        public int compare(Q3Peptide a, Q3Peptide b)
        {
            return a.peptideId < b.peptideId ? -1 : (a.peptideId == b.peptideId ? 0 : 1);
        }
    };

    /**
     *
     */
    private static class Q3RuntimeException extends RuntimeException
    {
        public Q3RuntimeException(String msg)
        {
            super(msg);
        }

        public Q3RuntimeException(String msg, Exception e)
        {
            super(msg, e);
        }
    }

    public static class IsotopicLabel
    {
        private char residue = '\0';
        private float massdiff = 0f;
        private float lightMass = 0f;
        private float heavyMass = 0f;

        public IsotopicLabel(char residue, float massdiff)
        {
            setResidue(residue);
            setMassDiff(massdiff);
        }

        public void setResidue(char residue)
        {
            // handle X! Tandem style termini
            if (residue == '[')
                this.residue = 'n';
            else if (residue == ']')
                this.residue = 'c';
            else
                this.residue = residue;

            if (this.residue == 'n' || this.residue == 'c')
                throw new IllegalArgumentException("N-terminal and C-terminal labels are not currently supported");
        }

        public char getResidue()
        {
            return this.residue;
        }

        public void setMassDiff(float massdiff)
        {
            this.massdiff = massdiff;
        }

        public float getMassDiff()
        {
            return this.massdiff;
        }

        public void setLightMass(float lightMass)
        {
            this.lightMass = lightMass;
        }

        public float getLightMass()
        {
            return this.lightMass;
        }

        public void setHeavyMass(float heavyMass)
        {
            this.heavyMass = heavyMass;
        }

        public float getHeavyMass()
        {
            return this.heavyMass;
        }
    }

    /**
     * Parse command-line args and launch Q3
     *
     * TODO: Since we are forced to be XPRESS-compatible, the commandline framework needs
     *   to handle multiple instances of same argument before we can replace this
     *
     * --q3 --labeledResidue=C --massDiff=3.0 --minPeptideProphet=0.75 --maxFracDeltaMass=15 --forceOutput --mimicXpress --noSentinels --out=out.pep.xml in.pep.xml
     */
    public static void run(String[] args) throws Exception
    {
        String inFileName = null;
        String outFileName = null;

        String alternateMzXmlDir = null;

        HashMap<Character,IsotopicLabel> labels = new HashMap<Character,IsotopicLabel>();

        char labeledResidue = '\0';
        float massDiff = 0f;
        float massTol = 25.f; // PPM

        float minPeptideProphet = 0.f;
        boolean filterByMinPeptideProphet = false;

        float maxFracDeltaMass = 9999999.9f;
        boolean maxFracDeltaMassIsPPM = true;
        boolean filterByMaxFracDeltaMass = false;

        boolean mimicXpress = false;
        boolean noSentinels = false;
        boolean forceOutput = false;
        boolean debugMode = false;
        boolean compatMode = false;
        boolean stripExistingQ3 = false;

        for (int i = 1; i < args.length; i++)
        {
            if (args[i].equals("--mimicXpress"))
                mimicXpress = true;
            else if (args[i].equals("--noSentinels"))
                noSentinels = true;
            else if (args[i].equals("--forceOutput"))
                forceOutput = true;
            else if (args[i].equals("--stripExistingQ3"))
                stripExistingQ3 = true;
            else if (args[i].equals("--debug"))
                debugMode = true;
            else if (args[i].equals("--compat"))
                compatMode = true;
            else if (args[i].startsWith("-n"))
            {
                // Look for XPRESS args of the form -nA,###
                String[] val = args[i].substring(2).split(",");
                if (val.length < 2 || val[0].length() != 1 || val[1].length() < 1)
                    throw new RuntimeException("Could not parse residue and mass from argument " + args[i]);
                char tmpResidue = val[0].charAt(0);
                float tmpDelta = Float.parseFloat(val[1]);
                labels.put(tmpResidue, new IsotopicLabel(tmpResidue, tmpDelta));
            }
            else if (args[i].startsWith("-m"))
            {
                // Look for XPRESS style mass tolerance; ignored for now
                // massTol = Float.parseFloat(args[i].substring(2));
                massTol = 25.f;
            }
            else if (args[i].startsWith("-d"))
            {
                // Look for mzXML files in an alternate directory
                alternateMzXmlDir = args[i].substring(2);
            }
            else if (args[i].startsWith("--"))
            {
                String[] param = args[i].split("=");
                if (param.length != 2)
                    quit("Unknown param: " + args[i]);

                String paramName = param[0];
                String paramVal = param[1];

                if ("--out".equals(paramName))
                    outFileName = paramVal;
                else if ("--labeledResidue".equals(paramName))
                    labeledResidue = paramVal.charAt(0);
                else if ("--massDiff".equals(paramName))
                    massDiff = Float.parseFloat(paramVal);
                else if ("--massTol".equals(paramName))
                    massTol = 25.f;
                else if ("--minPeptideProphet".equals(paramName))
                {
                    minPeptideProphet = Float.parseFloat(paramVal);
                    filterByMinPeptideProphet = true;
                }
                else if ("--maxFracDeltaMass".equals(paramName))
                {
                    if (paramVal.toLowerCase().endsWith("da"))
                        maxFracDeltaMassIsPPM = false;
                    maxFracDeltaMass = Float.parseFloat(paramVal);
                    filterByMaxFracDeltaMass = true;
                }
                else
                    quit("Unknown parameter: " + paramName);

                if (labeledResidue != '\0' && massDiff != 0f)
                {
                    labels.put(labeledResidue, new IsotopicLabel(labeledResidue, massDiff));
                    labeledResidue = '\0';
                    massDiff = 0f;
                }
            }
            else
            {
                if (null != inFileName)
                    quit("Can only process one pepXML file at a time. Found '" + args[i] + "'");
                inFileName = args[i];
            }
        }

        if (null == inFileName)
            quit("Must specify a pepXML file");

        if (null == outFileName)
        {
            outFileName = inFileName;
            forceOutput = true;
        }

        if (labels.size() == 0)
            quit("Must specify at least one isotopic label");

        try
        {
            Q3 q3 = new Q3(inFileName, labels, outFileName);

            if (null != alternateMzXmlDir)
                q3.setAlternateMzXmlDir(alternateMzXmlDir);
            if (filterByMinPeptideProphet)
                q3.setMinPeptideProphet(minPeptideProphet);
            if (filterByMaxFracDeltaMass)
                q3.setMaxFracDeltaMass(maxFracDeltaMass, maxFracDeltaMassIsPPM);
            q3.setCompatMode(compatMode);
            q3.setDebugMode(debugMode);
            q3.setForceOutput(forceOutput);
            q3.setMimicXpress(mimicXpress);
            q3.setNoSentinels(noSentinels);
            q3.setStripExistingQ3(stripExistingQ3);
            q3.apply();
        }
        catch (Exception e)
        {
//           e.printStackTrace(System.err);
            quit("Error running Q3: " + e.getMessage());
        }

    }

    private static void quit(String err)
    {
        System.err.println(err);
//        showUsage();
        System.exit(1);
    }


    /**
     *
     */
    public static void main(String[] av)
    {
        String infile = null;
        String outfile = null;
        // file = "/home/mfitzgib/AcrylForMarc/data/L_04_IPAS0012_AX01_RP_SG04to10.pep.xml";
        // file = "./tmp/test/IP0012_AX01_RP_SG04to10.pep.xml";
        // file = "./tmp/test/IP0012_AX01_RP_SG37to41.pep.xml";
        // file = "./tmp/test/IP0022_AX01_RP_SG16to27.pep.xml";

        if (av.length < 1)
        {
            System.err.println("Please specify a pepXML filename");
        }

        infile = av[0];
        if (av.length >= 2)
            outfile = av[1];

        try
        {
            Q3 q3 = new Q3(infile, 'C', 3.01f, outfile);
            q3.setMinPeptideProphet(0.75f);
            q3.setMaxFracDeltaMass(20.f, true);
            q3.setForceOutput(false);
            q3.setMimicXpress(false);
            q3.setNoSentinels(false);
            q3.apply();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

}
