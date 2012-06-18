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
package org.fhcrc.cpl.toolbox.proteomics;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.CPUTimer;
import org.fhcrc.cpl.toolbox.proteomics.gui.IntensityPlot;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.fhcrc.cpl.toolbox.*;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithHistogram;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.datastructure.FloatArray;
import org.fhcrc.cpl.toolbox.datastructure.FloatRange;
import org.systemsbiology.jrap.stax.Base64;
import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.MZXMLFileInfo;
import org.systemsbiology.jrap.stax.ScanHeader;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.lang.ref.SoftReference;
import java.nio.ByteBuffer;
import java.nio.channels.ClosedByInterruptException;
import java.nio.channels.ClosedChannelException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.Map;


/**
 * User: mbellew
 * Date: May 20, 2004
 * Time: 9:54:13 AM
 */
public class MSRun implements Serializable
{
    //flag to control the use of the sequential parser, which doesn't require an mzXML index
    protected boolean useSequentialParser = true;
//    protected boolean showTimeCharts = false;

    private static Logger _log = Logger.getLogger(MSRun.class);
    private static final int FLOATBYTES = 4;

    //	private static final long serialVersionUID = 8280766319681981127L;
    private static final long serialVersionUID = 8280766319681981128L;
    static float IMAGE_THRESHOLD = 10;

    // source file info
    String _filename;
    long _lastModified;
    long _filelength;
    private transient File _file;

    // for retrieving spectrum data
    transient FileChannel _fileChannel;
    transient MSXMLParser parser;

    MSScan[] _scans;
    FloatRange _mzRange;
    transient float[] _templateSpectrum = null;

    // BufferedImage is not serializable
    transient BufferedImage _image;
    float[] _scanArray;
    float[] _mzArray;
    float[] _intensityArray;

    // OK, these are redundant, but it's makes getIndex() easier
    transient int[] _indexMap;
    transient double[] _timeMap;

    //dhmay adding for ms2 scan indexing
    transient int[] _ms2IndexMap;

    transient MZXMLFileInfo _fileInfo = null;
    //dhmay removing "transient" 2/8/06 -- need MS2 data for AMT, and Mark thinks
    //MS3 data might be useful in the future
    MSScan[] _scans2 = null;
    MSScan[] _scans3 = null;
    TreeMap<Integer,MSScan> _allScans = null;

    private static boolean showIndexBuilderProgress = true;

    // used by Scan.getSpectrum()
    static CPUTimer readTimer = new CPUTimer("scan read");
    static CPUTimer decodeTimer = new CPUTimer("scan decode");
    static CPUTimer toFloatTimer = new CPUTimer("scan toFloat");
    transient private byte[] encodedData = null;



    private MSRun(String path) throws IOException
    {
        File f = new File(path);
        if (!f.exists())
            throw new FileNotFoundException(path);

        _file = f;
        _filename = f.getName();
        _lastModified = f.lastModified();
        _filelength = f.length();

        // data to generate BufferedImage
        FloatArray scanArray = new FloatArray();
        FloatArray mzArray = new FloatArray();
        FloatArray intensityArray = new FloatArray();
        _allScans = new TreeMap<Integer,MSScan>();
        try
        {                                           
            ApplicationContext.setMessage("Building index file: first pass");
            if (null != ApplicationContext.getFrame())
                ApplicationContext.getFrame().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

            parser = new MSXMLParser(path, useSequentialParser);
            int count = parser.getScanCount();
            _log.debug("JRAP scan count: " + count);
            int percent = Math.max(1, count / 100);
            float min = Float.MAX_VALUE;
            float max = 0.0F;
            int precursorScan = 0;
            //_scans = new MSScan[count];
            ArrayList list = new ArrayList();
            ArrayList list2 = new ArrayList();
            ArrayList list3 = new ArrayList();

            int currentNumScansPresent = 0;
            //long start = -1L;
            //long finish = -1L;
            //dhmay switching on sequential parser
            int loopMax = useSequentialParser ? count : parser.getMaxScanNumber();

            long jrapStart = System.currentTimeMillis();
            //for bookkeeping.  Not populated if showTimeCharts is false;
            java.util.List<Float> scanLoadTimes = new ArrayList<Float>();
            java.util.List<Float> scanSizes = new ArrayList<Float>();

            for (int i = 1; i <= loopMax; i++)
            {
                if (0 == (i % percent))
                {
                    if (showIndexBuilderProgress)
                        ApplicationContext.setMessage("Building index file: " + (currentNumScansPresent * 100 / Math.max(1,count)) + "%");
                    Thread.yield();
                }
                //start = System.currentTimeMillis();
                //dhmay switching on sequential parser
                ScanHeader scan = useSequentialParser ? parser.nextHeader() : parser.rapHeader(i);
                //finish = System.currentTimeMillis();
/*
                if (showTimeCharts)
                {
                    scanLoadTimes.add((float) (finish-start));
                    scanSizes.add((float)scan.getPeaksCount());
                }
*/
                //System.out.println("get scanheader for scan "+(i-1)+" "+(finish-start)+" ms");
//System.err.println("Scan index " + i);
                if (scan == null)
                    continue;
//System.err.println("\tScan #" + scan.getNum());
                currentNumScansPresent++;
                MSScan msScan = new MSScan(scan);
                _allScans.put(i,msScan);
                // ???? Consider: might be better to accept only msLevel=1 scanType=Full|unspecified
                if ("calibration".equals(scan.getScanType()))
                    continue;
                if ("zoom".equals(scan.getScanType()))
                    continue;
//System.err.println("*2");
                // Some ABI 4800 files contain partially filled headers when a scan detects no peaks
                // allowing ms2 scans with 0 peaks for MRM analyses
                if (scan.getMsLevel() == 1 && scan.getPeaksCount() <= 0)
                    continue;

//  System.err.println("*3, level=" + scan.getMsLevel() + ", num=" + scan.getNum());

                if (scan.getMsLevel() == 1)
                    precursorScan = scan.getNum();
                else
                {
                    if (-1 == scan.getPrecursorScanNum())
                        scan.setPrecursorScanNum(precursorScan);
                    if (scan.getMsLevel() == 2)
                        list2.add(msScan);
                    else if (scan.getMsLevel() == 3)
                        list3.add(msScan);
                    continue;
                }
                int index = list.size();
                list.add(msScan);
//   System.err.println("List size: " + index + ", added");

                float[][] spectrum = msScan.getSpectrum();
//   System.err.println("*4, " + spectrum.length);

                if (scan.getMsLevel() == 1 && spectrum[0].length > 0)
                {
                    min = Math.min(min, spectrum[0][0]);
                    max = Math.max(max, spectrum[0][spectrum[0].length - 1]);
                    if (min < 0 || max > 1e5)
                    {
                        RuntimeException x = new RuntimeException("Suspect mzxml data file, scan=" + i + " mz=" + (min < 0 ? min : max));
                        ApplicationContext.errorMessage(null, x);
                        throw x;
                    }
                }
                _computeImagePoints(index, spectrum, scanArray, mzArray, intensityArray);
            }
/*
            //dhmay adding conditional showing of charts
            if (showTimeCharts)
            {
                System.err.println("Total ms in JRAP: " + (System.currentTimeMillis() - jrapStart));
                PanelWithHistogram pwh = new PanelWithHistogram(scanLoadTimes, "Scan load ms"); pwh.setBreaks(200); pwh.displayInTab();
                new PanelWithScatterPlot(scanSizes, scanLoadTimes, "scan size (x) vs load time (y)").displayInTab();
            }
*/          _log.debug("Total ms in JRAP: "+(System.currentTimeMillis() - jrapStart));

            _scans = (MSScan[])list.toArray(new MSScan[0]);
            _scans2 = (MSScan[])list2.toArray(new MSScan[0]);
            _scans3 = (MSScan[])list3.toArray(new MSScan[0]);
            _mzRange = new FloatRange(min, max);

            _scanArray = scanArray.toArray(null);
            _mzArray = mzArray.toArray(null);
            _intensityArray = intensityArray.toArray(null);

            _initMaps();
        }
        finally
        {
            if (null != ApplicationContext.getFrame())
                ApplicationContext.getFrame().setCursor(Cursor.getDefaultCursor());
            ApplicationContext.setMessage("");
        }
    }


    private void _initMaps()
    {
        _indexMap = new int[_scans.length];
        _timeMap = new double[_scans.length];
        for (int i = 0; i < _scans.length; i++)
        {
            _indexMap[i] = _scans[i].getNum();
            _timeMap[i] = _scans[i].getDoubleRetentionTime();
        }

        //dhmay adding for ms2 scan indexing
        _ms2IndexMap = new int[_scans2.length];
        for (int i = 0; i < _scans2.length; i++)
        {
            _ms2IndexMap[i] = _scans2[i].getNum();
        }
    }

    /**
     * Turn index builder progress on or off (useful for batch scripts)
     */
    public static void setShowIndexBuilderProgress(boolean b)
    {
        showIndexBuilderProgress = b;
    }

    public String getFileName()
    {
        return _filename;
    }


    /* We have a bit of a mismatch with IntensityPlot
         * reduce and convert spectrum to expected format
         */
    private void _computeImagePoints(float scan, float[][] spectrum, FloatArray scanArray, FloatArray mzArray, FloatArray intensityArray)
    {
        int len = spectrum[0].length;
        if (len == 0) return;
        int mzCurrent = (int)spectrum[0][0];
        float intensityMax = (int)spectrum[1][0];
        for (int i = 1; i < len; i++)
        {
            int mz = (int)spectrum[0][i];
            float intensity = spectrum[1][i];
            if (mz == mzCurrent)
            {
                intensityMax = intensityMax > intensity ? intensityMax : intensity;
            }
            else
            {
                if (intensityMax > IMAGE_THRESHOLD)
                {
                    scanArray.add(scan);
                    mzArray.add((float)mzCurrent);
                    intensityArray.add(intensityMax);
                }
                mzCurrent = mz;
                intensityMax = intensity;
            }
        }
        scanArray.add(scan);
        mzArray.add((float)mzCurrent);
        intensityArray.add(intensityMax);
    }


    public static String _indexName(String filename)
    {
        File f = new File(filename);
        String s = f.getName() + ".inspect";
        if (null != f.getParent())
            s = f.getParent() + "/" + s;
        return s;
    }


/*	private static String _imageName(String filename)
		{
		File f = new File(filename);
		String s = f.getTextCode() + ".png";
		if (null != f.getParent())
			s = f.getParent() + "/" + s;
		return s;
		} */


    static MSRun _loadFromIndex(String path, String indexname)
    {
        File sourceFile = new File(path);
        File indexFile = new File(indexname);
        if (!indexFile.exists() || !sourceFile.exists())
            return null;
        _log.debug("Found .inspect file " + indexFile.getAbsolutePath() + ", checking...");
        if (indexFile.lastModified() < sourceFile.lastModified())
        {
            _log.debug(".inspect file is older than mzXML file, not using.");
            return null;
        }
        FileInputStream in = null;
        try
        {
            in = new FileInputStream(indexFile);
            ObjectInputStream ois = new ObjectInputStream(in);
            Object o = ois.readObject();
            MSRun run = (MSRun)o;

            //  don't check to make it easier to rename both files together
            //if (!sourceFile.getTextCode().equals(run._filename))
            //	return null;
            if (sourceFile.lastModified() != run._lastModified)
            {
                _log.debug("Incorrect modified date for mzXML file in .inspect file");
                return null;
            }
            if (sourceFile.length() != run._filelength)
            {
                _log.debug("Bad source file length in .inspect file");
                return null;
            }

            if (null == run._scanArray || null == run._mzArray || null == run._intensityArray)
                return null;
            if (null == run._mzRange)
            {
                _log.debug("Null _mzRange in .inspect file");
                return null;
            }
            // If index file contained no scans, assume it is busted.
            // dhmay 9/19/2007 removing a check for at least one MS1 scan, which fails on MRM files
            if (null == run._scans)
            {
                _log.debug("Null _scans in .inspect file");
                return null;
            }
            // dhmay adding 2/8/06.  If _null_ ms2/3 scans (empty is ok), assume it is busted
            if (null == run._scans2)
            {
                _log.debug("Null _scans2 in .inspect file");
                return null;
            }
            if (null == run._scans3)
            {
                _log.debug("Null _scans3 in .inspect file");
                return null;
            }

            run._initMaps();
            run._file = sourceFile;
            run._filename = sourceFile.getName();

            _log.debug("Index file loaded successfully");

            return run;
        }
        catch (FileNotFoundException x)
        {
            _log.debug("FileNotFoundException on .inspect file");
        }
        catch (IOException x)
        {
            _log.debug("IOException on index file");
        }
        catch (ClassNotFoundException x)
        {
            _log.debug("ClassNotFoundException on .inspect file");
        }
        catch (ClassCastException x)
        {
            _log.debug("ClassCastException on .inspect file");
        }
        finally
        {
            if (null != in)
                try
                {
                    in.close();
                }
                catch (IOException x)
                {
                }
        }
        return null;
    }


    boolean _writeIndex(String indexname)
    {
        File f = new File(indexname);
        try
        {
            // CONSIDER: use Deflater
            if (f.exists())
                f.delete();
            f.createNewFile();
            FileOutputStream out = new FileOutputStream(indexname);
            ObjectOutputStream oos = new ObjectOutputStream(out);
            oos.writeObject(this);
            out.close();
        }
        catch (FileNotFoundException x)
        {
            ApplicationContext.infoMessage(TextProvider.getText("WARNING_FAILED_TO_WRITE_AUXILIARY_FILE_FILE", f.getAbsolutePath()));
        }
        catch (IOException x)
        {
            ApplicationContext.infoMessage(TextProvider.getText("WARNING_FAILED_TO_WRITE_AUXILIARY_FILE_FILE", f.getAbsolutePath()));
        }
        return false;
    }

    public static MSRun load(String filename) throws IOException
    {
        return load(filename, true); // Write both .inspect and .ms2.tsv files unless run was read from index file
    }

    public static MSRun load(String filename, boolean writeIndex) throws IOException
    {
        if (!(new File(filename).exists()))
            throw new FileNotFoundException();

        String indexName = _indexName(filename);
        MSRun run = _loadFromIndex(filename, indexName);
        if (null == run)
        {
            _log.debug("No valid index file found, loading from mzXML");
            run = new MSRun(filename);
            if (writeIndex)
            {
                ApplicationContext.setMessage("Writing .inspect file...");
                run._writeIndex(indexName);

                /* NOTE: should really make saving these an Action, but since we've just read the .mzxml file,
                * I'd hate to have to read it again.
                */
                FeatureSet fs = run.getTandemFeatureSet(2);
                if (null != fs)
                {
                    ApplicationContext.setMessage("Writing MS2 features...");
                    File ms2File = new File(run.getFile().getPath() + ".ms2.tsv");
                    try
                    {
                        fs.save(ms2File);
                    }
                    catch (IOException e)
                    {
                        ApplicationContext.infoMessage(
                                TextProvider.getText(
                                        "WARNING_FAILED_TO_WRITE_AUXILIARY_FILE_FILE",
                                        ms2File.getAbsolutePath()));
                    }
                }
                fs = run.getTandemFeatureSet(3);
                if (null != fs)
                {
                    ApplicationContext.setMessage("Writing MS3 features...");
                    File ms3File = new File(run.getFile().getPath() + ".ms2.tsv");
                    try
                    {
                        fs.save(new File(run.getFile().getPath() + ".ms3.tsv"));
                    }
                    catch (IOException e)
                    {
                        ApplicationContext.infoMessage(
                                TextProvider.getText(
                                        "WARNING_FAILED_TO_WRITE_AUXILIARY_FILE_FILE",
                                        ms3File.getAbsolutePath()));
                    }
                }
            }
            ApplicationContext.setMessage("");
        }
        return run;
    }


    private long _checkIO() throws IOException
    {
        if (null != _fileChannel && _fileChannel.isOpen())
        {
            try
            {
                long size = _fileChannel.size();
                return size;
            }
            catch (IOException x)
            {
            }
        }
        close();
        _initIO();
        return _fileChannel.size();
    }


    private void _initIO()
    {
        String path = _file.getPath();
        try
        {
            _fileChannel = new RandomAccessFile(path, "r").getChannel();
        }
        catch (Exception x)
        {
            close();
            ApplicationContext.errorMessage(null, x);
        }
    }


    public void close()
    {
        if (null != _fileChannel)
            try
            {
                _fileChannel.close();
            }
            catch (IOException x)
            {
            }
        _fileChannel = null;
    }


    protected void finalize()
    {
        close();
    }

    public void invalidateImage()
    {
        _image = null;
    }

    public BufferedImage getImage(String colorScheme)
    {
        if (null == _image)
        {
            try
            {
                ApplicationContext.setMessage("Building image...");
                IntensityPlot plot = new IntensityPlot();
                float median = Spectrum.MedianSampled(_intensityArray, false);
                float threshold = median / 2;
                plot.setData(FloatArray.asFloatArray(_scanArray),
                        FloatArray.asFloatArray(_mzArray),
                        FloatArray.asFloatArray(_intensityArray));

                _image = plot.plot(threshold, true, colorScheme);

                float maxTIC = 0.0F;
                for (int s = 0; s < getScanCount(); s++)
                {
                    MSScan scan = getScan(s);
                    maxTIC = Math.max(maxTIC, scan.getTotIonCurrent());
                }

                // TIC Chart
                if (true)
                {
                    Graphics g = _image.getGraphics();
                    g.setColor(Color.BLACK);
                    int height = _image.getHeight();
                    int yPrev = 0;
                    for (int s = 0; s < getScanCount(); s++)
                    {
                        MSScan scan = getScan(s);
                        int y = (int)(99 * scan.getTotIonCurrent() / maxTIC);
                        if (s > 0)
                            g.drawLine(s - 1, height - yPrev - 1, s, height - y - 1);
                        yPrev = y;
                    }
                }

                ApplicationContext.setMessage("");
            }
            catch (Exception x)
            {
                ApplicationContext.errorMessage("Error generating image", x);
            }
        }
        return _image;
    }


    public int getScanCount()
    {
        return _scans.length;
    }

    public Map<Integer,MSScan> getAllScans()
    {
        return _allScans;
    }

    public MSScan getScanByNum(int num) {
        return this.getAllScans().get(num);
    }

    public MSScan[] getScans()
    {
        return _scans;
    }

    /**
     * Get a scan by index (NOT scan number)
     * @param scanIndex
     * @return
     */
    public MSScan getScan(int scanIndex)
    {
        return _scans[scanIndex];
    }

    public MSScan[] getMS2Scans()
    {
        return _scans2;
    }

    /**
     * get an MS2 scan by index (NOT scan number)
     * @param scanIndex
     * @return
     */
    public MSScan getMS2Scan(int scanIndex)
    {
        return _scans2[scanIndex];
    }

    public MSScan[] getMS3Scans()
    {
        return _scans3;
    }


    // this is the actual range, not the range reported by scan.getLowMz() and scan.getHighMz()
    public FloatRange getMzRange()
    {
        return _mzRange;
    }


    public float[] getTemplateSpectrum()
    {
        if (null == _templateSpectrum)
        {
            // since it's easy find spectrum with most data points
            int scanMax = 0;
            int countMax = 0;
            long totalPeaks = 0;
            for (int s = 0; s < _scans.length; s++)
            {
                totalPeaks += _scans[s].getPeaksCount();
                if (_scans[s].getPeaksCount() > countMax)
                {
                    countMax = _scans[s].getPeaksCount();
                    scanMax = s;
                }
            }
            _log.info(this.getFileName() + " " + totalPeaks + " total peaks");

            float[][] spectrum = _scans[scanMax].getSpectrum();
            _templateSpectrum = Spectrum.GenerateSpectrumTemplate(spectrum[0], getMzRange());
        }
        return _templateSpectrum;
    }


    public int getIndexForFeature(Feature f)
    {
        if (f.getTime() > 0)
            return getIndexForTime(f.getTime());

        int index = getIndexForScanNum(f.getScan());
        if (index<0)
            index= -(index+1);
        return index;
    }


    // maps a scan num in this run to an index
    public int getIndexForScanNum(int num)
    {
        return Arrays.binarySearch(_indexMap, num);
    }

    // maps an MS2 scan num in this run to an index
    public int getIndexForMS2ScanNum(int num)
    {
        return Arrays.binarySearch(_ms2IndexMap, num);
    }

    public int getScanNumForIndex(int index)
    {
        return _indexMap[index];
    }

    public int getScanNumForMS2Index(int index)
    {
        return _ms2IndexMap[index];
    }

    /*
	 * maps a scan num in this run to an index. If returnPrecursor is true, returns the
	 * index of the *previous* scan when the exact scan number is not present; useful for
	 * overlaying MS2 features
	 */
    public int getIndexForScanNum(int num, boolean returnPrecursor)
    {
        int i = Arrays.binarySearch(_indexMap, num);
        if (returnPrecursor && i < 0)
            i = -(i + 1);
        return i;
    }

    // map a time to a location in the scan array
    // useful for overlaying feature sets
    public int getIndexForTime(double time)
    {
        int i = Arrays.binarySearch(_timeMap, time);
        if (i < 0)
            i = -(i+1);
        i = Math.min(i, _scans.length-1);
        if (i == 0)
            return 0;
        double prev = _scans[i-1].getDoubleRetentionTime();
        double next = _scans[i].getDoubleRetentionTime();
        assert prev <= time && time <= next;
        return time-prev < next-time ? i-1 : i;
    }


    public class MSScan implements Serializable, Scan
    {
        //private static final long serialVersionUID = 3460298621454357258L;
        private static final long serialVersionUID = 3460298621454357259L;

        //Was the precursor m/z of this scan corrected by some program (namely msInspect)?
        //dhmay adding 3/31/2008
        protected boolean precursorMzCorrected = false;

        public MSScan(org.systemsbiology.jrap.stax.Scan s)
        {
            _scan = s.getHeader();
            float[][] spectrum = convertSpectrumToFloatArray(s.getMassIntensityList());
            if (null != spectrum)
                _spectrumRef = new SoftReference(spectrum);
        }


        public MSScan(ScanHeader s)
        {
            _scan = s;
        }


        public int getIndex()
        {
            return getIndexForScanNum(getNum());
        }


        /** For performance, use this method for caching of resampled spectrum
         public synchronized float[][] getResampledSpectrum()
         {
         FloatRange r = new FloatRange(getLowMz(), getHighMz());
         float[] s = _resampleRef != null ? (float[])_resampleRef.get() : null;
         if (null == s)
         {
         s = Spectrum.Resample(getSpectrum(), r, 36);
         _resampleRef = new SoftReference(s);
         }
         // CONSIDER: cache shared domain array
         float[] domain = new float[s.length];
         for (int i=0 ; i<domain.length ; i++)
         domain[i] = r.min + i/36.0F;
         return new float[][] {domain, s};
         }
         */

        /**
         * Return spectrum for this scan. Sychronize on run, so that different threads
         * can request spectra. Optionally cache spectra for future use.
         */
        public float[][] getSpectrum()
        {
            synchronized(MSRun.this)
            {
                return _getSpectrumInternal();
            }
        }                                    

        /**
         * can return NULL if thread is interrupted
         */
        private synchronized float[][] _getSpectrumInternal()
        {
//System.err.println("_getSpectrumInternal 1");
            float[][] spectrum = null == _spectrumRef ? null : (float[][])_spectrumRef.get();
            if (null != spectrum)
                return spectrum;
//  System.err.println("_getSpectrumInternal 1, null, will read");        `

            if (null == _fileChannel)
                _initIO();

            Throwable throwable = null;

            for (int ttry = 0 ; null == spectrum && ttry<2 ; ttry++)
            {
                try
                {
//System.err.println("Calling _getSpectrum()");
                    //First try without JRAP.  However,
                    //if it's not an mzXML file (i.e., it's an mzML file) we have to use JRAP, because
                    //the local code here doesn't know how to handle mzML
                    spectrum = _getSpectrum(!_filename.toUpperCase().contains(".MZXML"));
                    break;
                }
                catch (ClosedByInterruptException x)
                {
                }
                catch (ClosedChannelException x)
                {
                }
                catch (IOException x)
                {
                    _log.error(x);
                }
                try
                {
                    _checkIO();
                    continue;
                }
                catch (IOException x)
                {
                    throwable = x;
                }
                if (Thread.currentThread().isInterrupted())
                    return null;
            }
            // end of retry loop
            if (null != throwable)
                _log.error(null, throwable);

            if (null == spectrum)
            {
                //Second try.  Whether or not we tried JRAP before, we're trying it now
                try
                {
                    spectrum = _getSpectrum(true);
                    assert null != spectrum;
                    if (null == spectrum)
                        throw new NullPointerException();
                }
                catch (Throwable x)
                {
                    _log.fatal(x);
                   // System.exit(1);
                }
            }

            // Some files have zero MZ values at beginning and/or end!
            int length = spectrum[0].length;
            if (length > 0 && (spectrum[0][0] == 0.0 || spectrum[0][length-1] == 0.0))
            {
                while (length > 0 && spectrum[0][length-1] == 0.0)
                    length--;
                int start;
                for (start=0 ; start<length && spectrum[0][start] == 0.0 ; start++)
                    ;

                float[][] copy = new float[2][length-start];
                System.arraycopy(spectrum[0], start, copy[0], 0, copy[0].length);
                System.arraycopy(spectrum[1], start, copy[1], 0, copy[1].length);
                spectrum = copy;
            }

            _spectrumRef = new SoftReference(spectrum);
            return spectrum;
        }

        /**
         * use getSpectrum(false) when scanning large portions of the file serially
         * to first try avoiding the XML parser.
         * 
         * dhmay changing 20100727.  This had been checking useSequentialParser and passing it into JRAP,
         * but that's not appropriate here.  Sequential parser is only appropriate for full-file read, and this
         * method is only for reads of partial files.
         */
        private synchronized float[][] _getSpectrum(boolean jrap) throws IOException
        {
            float[][] spectrum = null;

            _log.debug("_getSpectrum, scan " + _scan.getNum());
            long offset = _scan.getScanOffset();            
            _log.debug("_getSpectrum, got scan offset, " + offset + ", jrap? " + jrap);
            if (!jrap && offset > 0)
            {
                // PERF HACK: avoid XML parser at all costs
                int count = _scan.getPeaksCount();
                long len = count * 2 * FLOATBYTES;
                long lenEnc = len / 3 * 4;
                lenEnc += 2048; // room for header
                lenEnc = Math.min(lenEnc, _fileChannel.size() - offset);

                ByteBuffer fileBuf = _fileChannel.map(FileChannel.MapMode.READ_ONLY, offset, lenEnc);
                fileBuf.position(0);
                byte[] buf = new byte[2048];
                fileBuf.get(buf);

                int i;
                findPeakList:
                for (i = 0; i < buf.length-6; i++)
                {
                    if (buf[i] != '<')
                        continue;
                    if (buf[i + 1] != 'p' || buf[i + 2] != 'e' || buf[i + 3] != 'a' || buf[i + 4] != 'k' || buf[i + 5] != 's')
                        continue;
                    for (i += 6; i < buf.length; i++)
                        if (buf[i] == '>')
                            break findPeakList;
                }
                if (i < buf.length)
                {
                    fileBuf.position(i + 1);
                    spectrum = _parseSpectrum(fileBuf, count);
                    if (null != spectrum)
                        return spectrum;
                }
            }
            _log.debug("About to try using JRAP, scan " + _scan.getNum() );

            // retry using JRAP parser
            for (int retry=0 ; retry<2 && null == spectrum; retry++)
            {
                //dhmay adding 20091028.  This null-check was missing, so would get NPE every time we reopened an
                //already-indexed file and encountered a scan with spectrum length 2048 (forcing us to hit the parser).
                //Not sure how this bug lasted this long... possibly introduced recently somehow?
                if (parser == null)
                    parser = new MSXMLParser(_file.getAbsolutePath(), false);
                org.systemsbiology.jrap.stax.Scan tmp = parser.rap(_scan.getNum());

                if (null != tmp)
                {
                    spectrum = convertSpectrumToFloatArray(tmp.getMassIntensityList());
                }
                else
                    parser = new MSXMLParser(_file.getAbsolutePath(), false);
            }
//System.err.println("end");
            return spectrum;
        }

        /**
         * JRAP used to deal with spectra as float[][], now it uses double[][].  This method converts one to
         * the other, which is a temporary waste of space and of CPU time.
         *
         * It might be better to move to storing spectra as double[][].  But this would take a lot of work, and
         * would increase the storage requirements for spectra significantly.
         *
         * TODO: move to storing spectra as double[][]?
         * @param spectrumDouble
         * @return
         */
        protected float[][] convertSpectrumToFloatArray(double[][] spectrumDouble)
        {
            if (spectrumDouble == null)
                return null;
            
            float[][] spectrum = new float[spectrumDouble.length][spectrumDouble[0].length];
            for (int i=0; i<spectrumDouble.length; i++)
            {
                for (int j=0; j<spectrumDouble[0].length; j++)
                    spectrum[i][j] = (float) spectrumDouble[i][j];
            }
            return spectrum;
        }


        private ScanHeader _scan;
        transient SoftReference _spectrumRef = null;
        transient SoftReference _resampleRef = null;

        public int getNum()
        {
            return _scan.getNum();
        }

        public int getMsLevel()
        {
            return _scan.getMsLevel();
        }

        public int getPeaksCount()
        {
            return _scan.getPeaksCount();
        }

        public String getPolarity()
        {
            return _scan.getPolarity();
        }

        public String getScanType()
        {
            return _scan.getScanType();
        }

        public int getCentroided()
        {
            return _scan.getCentroided();
        }

        public int getDeisotoped()
        {
            return _scan.getDeisotoped();
        }

        public int getChargeDeconvoluted()
        {
            return _scan.getChargeDeconvoluted();
        }

        public String getRetentionTime()
        {
            return _scan.getRetentionTime();
        }

        public float getStartMz()
        {
            return _scan.getStartMz();
        }

        public float getEndMz()
        {
            return _scan.getEndMz();
        }

        public float getLowMz()
        {
            return _scan.getLowMz();
        }

        public float getHighMz()
        {
            return _scan.getHighMz();
        }

        public float getBasePeakMz()
        {
            return _scan.getBasePeakMz();
        }

        public float getBasePeakIntensity()
        {
            return _scan.getBasePeakIntensity();
        }

        public float getTotIonCurrent()
        {
            return _scan.getTotIonCurrent();
        }

        public float getPrecursorMz()
        {
            return _scan.getPrecursorMz();
        }

        public void setPrecursorMz(float precursorMz)
        {
            _scan.setPrecursorMz(precursorMz);
        }

        public int getPrecursorScanNum()
        {
            return _scan.getPrecursorScanNum();
        }

        public int getPrecursorCharge()
        {
            return _scan.getPrecursorCharge();
        }

        public void setPrecursorCharge(int charge)
        {
            _scan.setPrecursorCharge(charge);
        }

        public float getCollisionEnergy()
        {
            return _scan.getCollisionEnergy();
        }

        public float getIonisationEnergy()
        {
            return _scan.getIonisationEnergy();
        }

        public int getPrecision()
        {
            return _scan.getPrecision();
        }

        public double getDoubleRetentionTime()
        {
            return _scan.getDoubleRetentionTime();
        }

        public String getFilterLine()
        {
            return _scan.getFilterLine();
        }


        public boolean isPrecursorMzCorrected()
        {
            return precursorMzCorrected;
        }

        public void setPrecursorMzCorrected(boolean precursorMzCorrected)
        {
            this.precursorMzCorrected = precursorMzCorrected;
        }

        public String toString()
        {
            return "MSScan(" + MSRun.this._filename + "," + getNum() + ")";
        }

        private float[][] _parseSpectrum(ByteBuffer buf, int count)
        {
            int len = count * 2 * FLOATBYTES;
            int lenEnc = len / 3 * 4;
            if (0 != len % 3)
                lenEnc += 4;

            // validate buffer size and position
            if (buf.position() + lenEnc + 1 >= buf.limit())
                return null;

            byte trailingByte = buf.get(buf.position() + lenEnc);
            // we expect '<' here, if not, then revert to safe parse
            if ('<' != trailingByte)
                return null;

            // seems to be faster to get byte[] than use buf
            assert readTimer.start();
            if (encodedData == null || encodedData.length < lenEnc)
                encodedData = new byte[Math.max(lenEnc,encodedData==null?128*1024:encodedData.length+4*1024)];
            buf.get(encodedData, 0, lenEnc);
            assert readTimer.stop();

            assert decodeTimer.start();
            byte[] byteData = encodedData; // HACK, we can use same array for output (decoded data is always shorter)
            int lenDecode = Base64.decode(encodedData, 0, lenEnc, byteData);
            assert decodeTimer.stop();

            if (null == byteData) // IO error or bad encoding
                return null;

            if (lenDecode % (FLOATBYTES * 2) != 0)
                return null;

            if (lenDecode / FLOATBYTES / 2 != count)
                return null;
            assert toFloatTimer.start();

            float[][] peakList = new float[2][count];
            for (int i = 0, p = 0, intBits = 0; p < count; p++)
            {
                intBits =
                        (((int)byteData[i]) << 24) |
                                ((((int)byteData[i + 1]) & 0xff) << 16) |
                                ((((int)byteData[i + 2]) & 0xff) << 8) |
                                (((int)byteData[i + 3]) & 0xff);
                i += FLOATBYTES;
                peakList[0][p] = Float.intBitsToFloat(intBits);
                intBits =
                        (((int)byteData[i]) << 24) |
                                ((((int)byteData[i + 1]) & 0xff) << 16) |
                                ((((int)byteData[i + 2]) & 0xff) << 8) |
                                (((int)byteData[i + 3]) & 0xff);
                i += FLOATBYTES;
                peakList[1][p] = Float.intBitsToFloat(intBits);
            }
            assert toFloatTimer.stop();

            return peakList;
        }
    }


    public File getFile()
    {
        return _file;
    }


    public void setFile(File file)
    {
        _file = file;
    }


    public String toString()
    {
        return "MSRun(" + _filename + ")";
    }


    public MZXMLFileInfo getHeaderInfo()
    {
        if (null == _fileInfo)
        {
            try
            {
                MSXMLParser parser = new MSXMLParser(_file.getPath(), useSequentialParser);
                _fileInfo = parser.rapFileHeader();
            }
            catch (IOException e)
            {
                throw new RuntimeException("Failed to parse mzXML file", e);
            }
        }
        return _fileInfo;
    }


    public FeatureSet getTandemFeatureSet(int level)
    {
        if (level != 2 && level != 3)
            throw new IllegalArgumentException();

        MSScan[] scans = level == 2 ? _scans2 : _scans3;
        if (null == scans || 0 == scans.length)
            return null;

        ArrayList list = new ArrayList(scans.length);
        for (int i = 0; i < scans.length; i++)
        {
            MSScan msScan = scans[i];
            Feature f = new Feature(msScan.getNum(), msScan.getPrecursorMz(), msScan.getBasePeakIntensity());
            f.setTime((float)msScan.getDoubleRetentionTime());
            if (msScan.getPrecursorCharge() > 0)
            {
                // If we know precursor charge we can fix up the mass as well
                f.setCharge(msScan.getPrecursorCharge());
                f.updateMass();
            }
            list.add(f);
        }
        Feature[] features = (Feature[])list.toArray(new Feature[0]);
        return new FeatureSet(features);
    }



    /**
     * Helper method to return a subset of the scans of this run
     * @param scanStart
     * @param scanEnd
     * @return the partial array of scans
     */
    public MSRun.MSScan[] getPartialScanArray(int scanStart, int scanEnd)
    {
        if (scanEnd<scanStart)
            return null;
        MSRun.MSScan[] result = new MSRun.MSScan[scanEnd - scanStart + 1];
        for (int scannum = scanStart; scannum <= scanEnd; scannum++)
        {
            int index = this.getIndexForScanNum(scannum);
//TODO: need to add handling for pepXml, negative indexes

            if (index < 0)
                result[scannum - scanStart] = this.getScan((-index) - 1);
            else
                result[scannum - scanStart] = this.getScan(index);
        }
        return result;
    }


}
