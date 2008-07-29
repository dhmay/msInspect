package org.fhcrc.cpl.viewer.commandline.modules;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.MSRun.MSScan;
import org.fhcrc.cpl.viewer.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.feature.Smooth2D;
import org.fhcrc.cpl.viewer.feature.Spectrum;
import org.fhcrc.cpl.viewer.feature.extraction.DefaultPeakCombiner;
import org.fhcrc.cpl.viewer.feature.extraction.PeakCombiner;
import org.fhcrc.cpl.viewer.gui.util.PanelWithScatterPlot;
import org.fhcrc.cpl.viewer.gui.util.ScatterPlotDialog;
import org.jfree.chart.renderer.AbstractRenderer;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.labkey.common.tools.ApplicationContext;
import org.labkey.common.tools.FloatRange;

public class CalibrationUsingLockspray extends BaseCommandLineModuleImpl
		implements CommandLineModule {

	File mzXmlFile, outFile;
	FeatureSet featureSetFile;
	float lockmassMz = 785.8426F;
	int lockmassCharge = 2;
	float massWindow = 0.2F;
	boolean usePPM = false;
	boolean plotCalibration = false;

	public CalibrationUsingLockspray() {
		init();
	}

	protected void init() {
		mCommandName = "CalibrationUsingLockspray";

		mHelpMessage = "Calibrate a FeatureSet or file using the scans with type=calibration in a mzXML raw file.";

		mShortDescription = "Calibrate a FeatureSet file using type=calibration scans.";

		CommandLineArgumentDefinition[] argDefs = {
				createFileToReadArgumentDefinition("mzXmlFile", true,
						"mzXML File"),
				createFeatureFileArgumentDefinition("featureSetFile", false,
						"FeatureSet File"),
				createFileToWriteArgumentDefinition("featureout", false,
						"Output file (featureSet)"),
				createBooleanArgumentDefinition("plotCalibration", false,
						"Plot the calibration", plotCalibration), };
		addArgumentDefinitions(argDefs);
	}

	@Override
	public void assignArgumentValues() throws ArgumentValidationException {
		mzXmlFile = getFileArgumentValue("mzXmlFile");
		featureSetFile = getFeatureSetArgumentValue("featureSetFile");
		outFile = getFileArgumentValue("featureout");
		plotCalibration = getBooleanArgumentValue("plotCalibration");
	}

	@Override
	public void execute() throws CommandLineModuleExecutionException {
		MSRun run;
		try {
			run = MSRun.load(mzXmlFile.toString());
		} catch (IOException e) {
			throw new CommandLineModuleExecutionException(
					"Error processing file", e);
		}
		MSScan[] calibScans = run.getCalibrationScans();
		if (calibScans == null || calibScans.length == 0) {
			throw new CommandLineModuleExecutionException(
					"No calibration scans in " + run.getFileName());
		}

		Feature[] lockmassFeatures = new Feature[calibScans.length];
		for (int i = 0; i < calibScans.length; i++) {
			float[][] spectrum = calibScans[i].getSpectrum();

			// other processing
			spectrum = Spectrum.ResampleSpectrum(spectrum, new FloatRange(
					lockmassMz - 1, lockmassMz + 4), 36, false);
			int window = 72;
			float x[] = spectrum[1];
			float bg[] = Spectrum.MinimaWindow(x, spectrum[0].length, window,
					null);
			for (int j = 0; j < bg.length; j++) {
				bg[j] = Math.max(0, x[j] - bg[j]);
			}
			spectrum = new float[][] { spectrum[0], bg };
			double s = Smooth2D.smoothYfactor;
			spectrum[1] = Spectrum.FFTsmooth(spectrum[1], s, false);
			Spectrum.Peak[] peaks = Spectrum.WaveletPeaksD3(spectrum);

			PeakCombiner combiner = new DefaultPeakCombiner();
			Feature[] features = combiner.createFeaturesFromPeaks(run, peaks);

			float minDiff = massWindow;
			Feature lockmassFeature = null;
			List<Feature> thrownFeatures = new ArrayList<Feature>(
					features.length);
			for (Feature feature : features) {
				if (feature.getCharge() != lockmassCharge) {
					thrownFeatures.add(feature);
					continue;
				}
				float diff = Math.abs(feature.getMz() - lockmassMz);
				if (diff < minDiff) {
					if (lockmassFeature != null) {
						thrownFeatures.add(lockmassFeature);
					}
					lockmassFeature = feature;
					minDiff = diff;
				} else {
					thrownFeatures.add(lockmassFeature);
				}
			}
			lockmassFeatures[i] = lockmassFeature;
		}

		if (plotCalibration) {
			ScatterPlotDialog spd = new ScatterPlotDialog();
			XYSeries series = new XYSeries("Wavelet peaks, within "
					+ massWindow + " Da window");
			double maxDiff = 0;
			for (int index = 0; index < lockmassFeatures.length; index++) {
				if (lockmassFeatures[index] != null) {
					double time = calibScans[index].getDoubleRetentionTime();
					double diff = lockmassFeatures[index].getMz() - lockmassMz;
					if (Math.abs(diff) > maxDiff) {
						maxDiff = Math.abs(diff);
					}
					series.add(time, diff);
				} else {
					series
							.add(calibScans[index].getDoubleRetentionTime(),
									null);
				}
			}
			PanelWithScatterPlot panelWScatterPlot = spd
					.getPanelWithScatterPlot();
			panelWScatterPlot.addSeries(series);
			StandardXYItemRenderer renderer = panelWScatterPlot.getRenderer();
			renderer.setPlotLines(true);
			NumberFormat decFormat = NumberFormat.getInstance();
			panelWScatterPlot.setAxisLabels(
					"Time",
					"Mass Deviation (Da), " + decFormat.format(maxDiff)
							+ "Da=="
							+ decFormat.format(maxDiff / lockmassMz * 1000000)
							+ "ppm");
			spd
					.setTitle("Lockmass peaks within a " + massWindow
							+ " Da window");
			spd.setVisible(true);
		}

		boolean cutEveryOther = false;
		if (cutEveryOther) {
			for (int i = 0; i < lockmassFeatures.length; i++) {
				if (i % 2 == 1) {
					lockmassFeatures[i] = null;
				}
			}
		}

		if (featureSetFile != null) {
			Feature[] features = featureSetFile.getFeatures().clone();
			Arrays.sort(features, new Feature.ScanAscComparator());
			int i = 0;
			while (i < lockmassFeatures.length && lockmassFeatures[i] == null) {
				i++;
			}
			if (i == lockmassFeatures.length - 1) {
				ApplicationContext.infoMessage("No lockmass features found");
				return;
			}
			Feature before = lockmassFeatures[i];
			// if feature is before the first found lockmass feature, then
			// only use the first lockmass feature for correction.
			Feature after = before;
			for (Feature feature : features) {
				int scanNum = feature.getScan();
				// check if the feature has moved past the 'after' feature.
				while (scanNum > after.getScan()) {
					if (i < lockmassFeatures.length - 2) {
						i++;
						if (lockmassFeatures[i + 1] == null) {
							continue;
						}
						before = after;
						after = lockmassFeatures[i + 1];
					} else {
						// since feature is after the last found lockmass
						// feature, then
						// only use the last lockmass feature for correction.
						before = after;
						break;
					}
				}
				float correction = lockmassMz
						- (before.getMz() + after.getMz()) / 2;
				if (usePPM) {
					correction = correction / lockmassMz * feature.getMz();
				}
				feature.setMz(feature.getMz() + correction);
				feature.updateMass();
			}
			try {
				featureSetFile.save(outFile);
			} catch (IOException e) {
				ApplicationContext
						.errorMessage("Error while trying to write file '"
								+ outFile + "'", e);
			}
		}
	}
}
