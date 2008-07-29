package org.fhcrc.cpl.viewer.commandline.modules;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.gui.util.ChartDialog;
import org.fhcrc.cpl.viewer.gui.util.PanelWithBarChart;
import org.fhcrc.cpl.viewer.gui.util.PanelWithBlindImageChart;
import org.fhcrc.cpl.viewer.gui.util.PanelWithChart;
import org.fhcrc.cpl.viewer.gui.util.PanelWithHistogram;
import org.fhcrc.cpl.viewer.gui.util.PanelWithLineChart;
import org.jfree.chart.ChartFactory;
import org.jfree.data.statistics.HistogramBin;
import org.jfree.data.statistics.SimpleHistogramBin;
import org.jfree.data.statistics.SimpleHistogramDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.labkey.common.tools.ApplicationContext;

public class MSeFeatureSetCombinerCLM extends BaseCommandLineModuleImpl
{
	File mzXmlFile, outFile, savePrecursors, saveFragments;
	Map<Feature, Feature[]> pre_to_frag = new HashMap<Feature, Feature[]>();
	FeatureSet precursors, fragments;
	int inBothCounter, inNoneCounter, multipleParentsCounter;
	int skippedFragmentsCounter, skippedFragmentScans, noFragmentsCounter;
	int precursorsCount;
	Map<Integer, Integer> fragsInSkippedScansHistogramData = new HashMap<Integer, Integer>();

	/**
	 * Use the fragments from both neighboring scans of a precursor,
	 * instead of only the fragment scans that contain the precursor. 
	 */
	boolean useFragsFromBoth = true;
	protected boolean plotStatistics = false;
	/**
	 * The maker class that knows how to create a format.
	 * "MGF" or "PKL" for instance.
	 * See FormatMaker.makers for a list.
	 */
	String maker = "MGF";

	public MSeFeatureSetCombinerCLM()
	{
		init();
	}

	protected void init()
	{
		mCommandName = "msefeaturesetcombiner";

		mHelpMessage = "Correlate precursor and fragment featuresets from a MS^E experiment in order to get a 'normal' DDA like peak list.";

		mShortDescription = "Combine precursor and fragment featuresets from a MS^E experiment.";

		CommandLineArgumentDefinition[] argDefs =
			{
				createFileToReadArgumentDefinition("mzXMLfile", true, 
						"The mzXML file that the featuresets has been created from"),
				createFeatureFileArgumentDefinition("precursors", true,
						"featureset of the precursors"),
				createFeatureFileArgumentDefinition("fragments", true,
						"featureset of the fragments"),
				createFileToWriteArgumentDefinition("outfile", true, "Output file"),
				createBooleanArgumentDefinition("plotstats", false,
                     "Plot statistics",
                     plotStatistics),
				createEnumeratedArgumentDefinition("formatMaker", false,
						"Which format should be output?",
						FormatMaker.makers, maker),
			};
		addArgumentDefinitions(argDefs);
	}

	@Override
	public String toString()
	{
		return mCommandName + " " + mzXmlFile.getName();
	}
	
	@Override
	public void assignArgumentValues() throws ArgumentValidationException
	{
		mzXmlFile = getFileArgumentValue("mzXMLfile");
		precursors = getFeatureSetArgumentValue("precursors");
		fragments = getFeatureSetArgumentValue("fragments");
		outFile = getFileArgumentValue("outfile");
		if (!hasArgumentValue("formatMaker"))
		{
			if (outFile.getName().toLowerCase().endsWith("pkl"))
			{
				maker = PKLMaker.makerName;
			}
			else if(outFile.getName().toLowerCase().endsWith("mgf"))
			{
				maker = MGFMaker.makerName;
			}
		}
		else
		{
			maker = getStringArgumentValue("formatMaker");
		}
		plotStatistics = getBooleanArgumentValue("plotstats");
        pre_to_frag.clear();
        inBothCounter = 0;
        inNoneCounter = 0;
        multipleParentsCounter = 0;
        noFragmentsCounter = 0;
        skippedFragmentsCounter = 0;
        skippedFragmentScans = 0;
        precursorsCount = 0;
    	fragsInSkippedScansHistogramData.clear();
	}

	@Override
	public void execute() throws CommandLineModuleExecutionException
	{
		PrintWriter outPW;
		MSRun run;

		try
		{
			outPW = getPrintWriter(outFile);
		}
		catch (FileNotFoundException e)
		{
			throw new CommandLineModuleExecutionException(e);
		}

		try
		{
			run = MSRun.load(mzXmlFile.getAbsolutePath());
			if (run == null)
			{
				throw new CommandLineModuleExecutionException(
						"Error opening run from file "
								+ mzXmlFile.getAbsolutePath());
			}
		}
		catch (IOException e)
		{
			outPW.close();
			throw new CommandLineModuleExecutionException(e);
		}

		Feature[] fragArray = fragments.getFeatures();
		Arrays.sort(fragArray, new Feature.ScanAscComparator());
		Feature[] precArray = precursors.getFeatures();
		Arrays.sort(precArray, new Feature.ScanAscComparator());
		int i = 0;
		int prevScan = 0, scan = 0;
		// Beware, I'm reusing this List Object, it changes with each
		// precursor.
		List<Feature> prevScanFrags = new ArrayList<Feature>();
		List<Feature> nextScanFrags = new ArrayList<Feature>();
		boolean inPrev = false, inNext = false;
		boolean prevAssigned = false, nextAssigned = false;
		int prevFragScanNum = 0, nextFragScanNum = 0;
		/**
		 * use for finding fragmentation scans. run.getIndexForScanNum(scanNum,
		 * true) where scanNum is a precursor scan number will give the scan
		 * index of the previous fragmentation scan.
		 */
		run.setMsLevel(2);
		for (Feature precursor : precArray) {
			prevScan = scan;
			scan = precursor.getScan();
			inPrev = inNext = false;
			if (prevScan != scan) {
				int prevFragScanIndex = run.getIndexForScanNum(scan, true);
				prevFragScanNum = run.getScan(prevFragScanIndex).getNum();
				if (prevFragScanNum == nextFragScanNum) {
					// Precursor only moved one scan number forward.
					prevAssigned = nextAssigned;
					nextAssigned = false;

					// Put the fragments from the next list into the previous
					// list.
					List<Feature> tmp = prevScanFrags;
					prevScanFrags = nextScanFrags;
					nextScanFrags = tmp;
					nextScanFrags.clear();
					if (prevScanFrags.size() > 0) {
						for (Feature frag : prevScanFrags) {
							if (frag.getMass() == precursor.getMass()) {
								inPrev = true;
							}
						}
					}
				} else {
					prevAssigned = nextAssigned = false;
					prevScanFrags.clear();
					nextScanFrags.clear();
					int lastSkippedScan = -2; // no scan numbers should match this.
					int beforeSkippedFrags = skippedFragmentsCounter;
					while (i < fragArray.length
							&& fragArray[i].getScan() < prevFragScanNum) {
						// We are skipping some fragment features, throwing them
						// away, because no precursor was found for them.
						// XXX: Might want to record this, so that it's possible
						// to target a search for precursors for the fragments.
						if (fragArray[i].getScan() != lastSkippedScan) {
							// First skipped fragScan should get in here.
							skippedFragmentScans++;
							
							int fragsInLastSkippedScan = skippedFragmentsCounter - beforeSkippedFrags;
							if (fragsInLastSkippedScan != 0) {
								Integer count = fragsInSkippedScansHistogramData.get(fragsInLastSkippedScan);
								if (count == null) {
									count = 0;
								}
								count += 1;
								fragsInSkippedScansHistogramData.put(fragsInLastSkippedScan, count);
							} else {
								assert lastSkippedScan == -2;
							}
							beforeSkippedFrags = skippedFragmentsCounter;
							// Only get in here again when the skipped fragments are from another scan.
							lastSkippedScan = fragArray[i].getScan();
						}
						skippedFragmentsCounter++;
						
						i++;
					}
					int fragsInLastSkippedScan = skippedFragmentsCounter - beforeSkippedFrags;
					if (fragsInLastSkippedScan != 0) {
						Integer count = fragsInSkippedScansHistogramData.get(fragsInSkippedScansHistogramData);
						if (count == null) {
							count = 0;
						}
						count += 1;
						fragsInSkippedScansHistogramData.put(fragsInLastSkippedScan, count);
					}

					while (i < fragArray.length
							&& fragArray[i].getScan() == prevFragScanNum) {
						if (fragArray[i].getMass() == precursor.getMass()) {
							inPrev = true;
						}
						prevScanFrags.add(fragArray[i]);
						i++;
					}
				}
				nextFragScanNum = run.getScan(prevFragScanIndex + 1).getNum();

				while (i < fragArray.length
						&& fragArray[i].getScan() == nextFragScanNum) {
					if (fragArray[i].getMass() == precursor.getMass()) {
						inNext = true;
					}
					nextScanFrags.add(fragArray[i]);
					i++;
				}
			} else {
				for (Feature frag : prevScanFrags) {
					if (frag.getMass() == precursor.getMass()) {
						inPrev = true;
					}
				}
				for (Feature frag : nextScanFrags) {
					if (frag.getMass() == precursor.getMass()) {
						inNext = true;
					}
				}
			}
			if (prevScanFrags.size() == 0 && nextScanFrags.size() == 0) {
				noFragmentsCounter++;
			}
			if (useFragsFromBoth || inNext == inPrev) {
				// Weird. inNext && inPrev could be true, but it's not good.
				// Means the precursor was found among the fragments in both
				// the previous and the next scan.
				// Would probably mean that the fragments contain the
				// precursor in multiple charge states and
				// the charge states peaked in different scans.
				if (inNext && inPrev) {
					inBothCounter++;
				} else if (!inNext && !inPrev) {
					inNoneCounter++;
				}
				
				if (prevAssigned || nextAssigned) {
					multipleParentsCounter++;
				}
				nextAssigned = prevAssigned = true;
				Feature[] tmp = new Feature[prevScanFrags.size()
						+ nextScanFrags.size()];
				int j = 0;
				for (Feature feature : prevScanFrags) {
					tmp[j++] = feature;
				}
				for (Feature feature : nextScanFrags) {
					tmp[j++] = feature;
				}
				pre_to_frag.put(precursor, tmp);
			} else if (inPrev) {
				if (prevAssigned) {
					multipleParentsCounter++;
				}
				prevAssigned = true;
				pre_to_frag.put(precursor, prevScanFrags
						.toArray(new Feature[0]));
			} else {
				if (nextAssigned) {
					multipleParentsCounter++;
				}
				nextAssigned = true;
				pre_to_frag.put(precursor, nextScanFrags
						.toArray(new Feature[0]));
			}
		}
		fragArray = precArray = null;
		prevScanFrags = nextScanFrags = null;

		Set<Map.Entry<Feature, Feature[]>> entries = pre_to_frag.entrySet();
		precursorsCount = entries.size();
		
		FormatMaker formatMaker;
		try
		{
			Class<? extends FormatMaker> makerClass = (Class<? extends FormatMaker>) Class.forName(getClass().getName() + "$" + maker + FormatMaker.classSuffix);
			Constructor<? extends FormatMaker> makerConstructor = makerClass.getDeclaredConstructor(new Class[] {getClass()});
			formatMaker = makerConstructor.newInstance(this);
		}
		catch (Exception e) 
		{
			throw new CommandLineModuleExecutionException("Unable to use the formatMaker " + maker, e);
		}
		
		formatMaker.init(outPW);
		for (Map.Entry<Feature, Feature[]> entry : entries) {
			Feature pre = entry.getKey();
			Feature[] frags = entry.getValue();
			formatMaker.query(pre, frags);
		}
		formatMaker.terminate();
		outPW.close();
		
		if (plotStatistics) {
			String[] xval = new String[fragsInSkippedScansHistogramData.size()];
			int[] yval = new int[fragsInSkippedScansHistogramData.size()];
			int x = 0;
			for (Entry<Integer, Integer> entry : fragsInSkippedScansHistogramData.entrySet()) {
				xval[x] = entry.getKey().toString();
				yval[x] = entry.getValue();
				x++;
			}
			PanelWithChart pwh = new PanelWithBarChart(xval, yval, "frags in skipped scans");
			ChartDialog chartdlg = new ChartDialog(pwh);
			chartdlg.setVisible(true);
		}
	}
	
	private static boolean massEquals(float mass1, float mass2) {
		return Math.abs(mass1 - mass2) / mass1 * 1E6 < 1;  
	}
	
	interface FormatMaker {
		void init(PrintWriter pw);
		void query(Feature pre, Feature[] frags);
		void terminate();
		public String getMakerName();
		public static final String[] makers = {MGFMaker.makerName, PKLMaker.makerName}; 
		public static final String classSuffix = "Maker";
	}
	
	class PKLMaker implements FormatMaker {
		PrintWriter outPW;
		public static final String makerName = "PKL";
		
		public String getMakerName()
		{
			return makerName;
		}
		
		public void init(PrintWriter pw) {
			outPW = pw;
		}
		
		public void query(Feature pre, Feature[] frags) {
			if (frags.length == 0)
			{
				return;
			}
			outPW.println(pre.getMz() + " " + pre.getTotalIntensity() + " " + pre.getCharge());
			for (Feature frag : frags) {
				//search engines do not expect the precursor in the fragment scan
				if (massEquals(frag.getMass(), pre.getMass())) {
					continue;
				}
				outPW.println(frag.getMz() + " " + frag.getTotalIntensity());
			}
			outPW.println();
		}
		
		public void terminate() { }
	}
	
	class MGFMaker implements FormatMaker {
		PrintWriter outPW;
		public static final String makerName = "MGF";
		
		public String getMakerName() {
			return makerName;
		}
		
		public void init(PrintWriter pw) {
			outPW = pw;
			outPW.println("#inBothCounter = " + inBothCounter);
			outPW.println("#inNoneCounter = " + inNoneCounter);
			outPW.println("#precursorsCount = " + precursorsCount);
			outPW.println("#multipleParentsCounter = " + multipleParentsCounter);
			outPW.println("#skippedFragmentsCounter = " + skippedFragmentsCounter);
			outPW.println("#skippedFragmentScans = " + skippedFragmentScans);
			outPW.println("#noFragmentsCounter = " + noFragmentsCounter);
		}

		public void query(Feature pre, Feature[] frags) {
			outPW.println("BEGIN IONS");
			outPW.println("TITLE=" + pre.toString());
			outPW.println("PEPMASS=" + pre.getMz() + " " + pre.getTotalIntensity());
			if (pre.getCharge() != 0) {
				outPW.println("CHARGE=" + pre.getCharge() + "+");
			}
			outPW.println("SCANS=" + pre.getScanFirst() + "-" + pre.getScanLast());
			outPW.println("RTINSECONDS=" + pre.getTime());
			for (Feature frag : frags) {
				//search engines do not expect the precursor in the fragment scan
				if (frag.getMass() == pre.getMass()) {
					continue;
				}
				outPW.println(frag.getMz() + " " + frag.getTotalIntensity());
			}
			outPW.println("END IONS");
			outPW.println();
		}
		
		public void terminate() { }
	}
}
