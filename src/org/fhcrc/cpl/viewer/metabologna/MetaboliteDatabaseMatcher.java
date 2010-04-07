package org.fhcrc.cpl.viewer.metabologna;

import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.matching.FeatureSetMatcher;
import org.fhcrc.cpl.toolbox.proteomics.MassUtilities;
import org.fhcrc.cpl.toolbox.chem.ChemicalCompound;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithScatterPlot;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;
import org.fhcrc.cpl.toolbox.statistics.BasicStatistics;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;

import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;
import java.io.IOException;

public class MetaboliteDatabaseMatcher
{
    public static final float DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM = 5;
    public static final float DEFAULT_CALIBRATION_MAX_STUD_RES = 1;
    public static final float DEFAULT_CALIBRATION_MAX_LEVERAGE = 6;
    public static final int DEFAULT_CALIBRATION_REGRESSION_DEGREES_FREEDOM = 1;



    protected float calMassTolerancePPM = DEFAULT_CALIBRATION_LOOSE_TOLERANCE_PPM;
    protected float calMaxLeverage = DEFAULT_CALIBRATION_MAX_LEVERAGE;
    protected double calMaxStudRes = DEFAULT_CALIBRATION_MAX_STUD_RES;
    protected int calRegressionDoF = DEFAULT_CALIBRATION_REGRESSION_DEGREES_FREEDOM;


    protected List<ChemicalCompound> databaseCompoundsByMass = null;
    protected boolean showCharts = false;


    public void calibrateFeatureMasses(Feature[] features) throws IOException
    {
        {
            //First by time, then by mass
            Map<Feature, List<ChemicalCompound>> featureMassMatchMap =
                massMatch(features, databaseCompoundsByMass, calMassTolerancePPM);
            if (showCharts)
            {
                List<Double> masses = new ArrayList<Double>();
                List<Float> deltaMassesPPM = new ArrayList<Float>();
                for (Feature feature : featureMassMatchMap.keySet())
                {
                    for (ChemicalCompound compound : featureMassMatchMap.get(feature))
                    {
                         masses.add(compound.getMass());
                         deltaMassesPPM.add(MassUtilities.convertDaToPPM(feature.getMass() - (float) compound.getMass(),
                                 (float) compound.getMass()));
                    }
                }
                new PanelWithScatterPlot(masses, deltaMassesPPM, "DeltaMassVsMass_beforeRTCal").displayInTab();
            }
            double[] calibrationCoefficients = calcMassCalibrationUsingMatches(featureMassMatchMap, true);
            for (Feature feature : features)
            {
                feature.setMass(feature.getMass() -
                        (float) RegressionUtilities.mapValueUsingCoefficients(calibrationCoefficients,
                                feature.getTime()));
                if (feature.comprised != null)
                {
                    for (Spectrum.Peak peak : feature.comprised)
                    {
                        if (peak == null)
                            continue;
                        float oldPeakMass = Feature.convertMzToMass(peak.mz, feature.charge);
                        float newPeakMass = oldPeakMass -
                                (float) RegressionUtilities.mapValueUsingCoefficients(calibrationCoefficients,
                                feature.getTime());
                        peak.setMz(Feature.convertMassToMz(newPeakMass, feature.charge));
//System.err.println("old=" + oldPeakMass + ", new=" + newPeakMass + ", ppmdiff: " + MassUtilities.convertDaToPPM(newPeakMass-oldPeakMass, oldPeakMass));
                    }
                }
            }

            featureMassMatchMap =
                massMatch(features, databaseCompoundsByMass, calMassTolerancePPM);
//                if (showCharts)
//                {
//                    List<Float> masses = new ArrayList<Float>();
//                    List<Float> deltaMassesPPM = new ArrayList<Float>();
//                    for (Feature feature : featureMassMatchMap.keySet())
//                    {
//                        for (ChemicalCompound thing : featureMassMatchMap.get(feature))
//                        {
//                             masses.add(thing.getMass());
//                             deltaMassesPPM.add(MassUtilities.convertDaToPPM(feature.getMass() - thing.getMass(), thing.getMass()));
//                        }
//                    }
//                    new PanelWithScatterPlot(masses, deltaMassesPPM, "DeltaMassVsMass_beforeMassCal").displayInTab();
//                }
            calibrationCoefficients = calcMassCalibrationUsingMatches(featureMassMatchMap, false);
            for (Feature feature : features)
            {
                feature.setMass(feature.getMass() -
                        (float) RegressionUtilities.mapValueUsingCoefficients(calibrationCoefficients,
                                feature.getMass()));
                if (feature.comprised != null)
                {
                    for (Spectrum.Peak peak : feature.comprised)
                    {
                        if (peak == null)
                            continue;
                        float oldPeakMass = Feature.convertMzToMass(peak.mz, feature.charge);
                        float newPeakMass = oldPeakMass -
                                (float) RegressionUtilities.mapValueUsingCoefficients(calibrationCoefficients,
                                oldPeakMass);
                        peak.setMz(Feature.convertMassToMz(newPeakMass, feature.charge));
//System.err.println("old=" + oldPeakMass + ", new=" + newPeakMass + ", ppmdiff: " + MassUtilities.convertDaToPPM(newPeakMass-oldPeakMass, oldPeakMass));
                    }
                }
            }
        }//end if shouldCalibrateMasses

    }

    protected double[] calcMassCalibrationUsingMatches(Map<Feature, List<ChemicalCompound>> featureMassMatchMap,
                                                       boolean byTime) throws IOException
    {
        List<Double> massDiffsList = new ArrayList<Double>();
        List<Float> deltaMassPPMListPreCal = new ArrayList<Float>();
        List<Double> featureXvalListPreCal = new ArrayList<Double>();

        String chartSuffix = byTime ? "_RT" : "_DA";

        for (Feature feature : featureMassMatchMap.keySet())
        {
            for (ChemicalCompound massAndName : featureMassMatchMap.get(feature))
            {
                double deltaMassDa = feature.getMass() - massAndName.getMass();
                massDiffsList.add(deltaMassDa);
                deltaMassPPMListPreCal.add(
                        MassUtilities.convertDaToPPM((float)deltaMassDa, (float) massAndName.getMass()));
                featureXvalListPreCal.add(byTime ? (double) feature.getTime() : feature.getMass());
            }
        }
        double[] featureXvals = new double[featureXvalListPreCal.size()];
        double[] massDiffs = new double[massDiffsList.size()];
        for (int i=0; i<featureXvals.length; i++)
        {
            featureXvals[i] = featureXvalListPreCal.get(i);
            massDiffs[i] = massDiffsList.get(i);
        }

        int[] indexesLowStuff = RegressionUtilities.selectIndexesWithLowLeverageAndStudentizedResidual(
                featureXvals, massDiffs, calMaxLeverage,
                calMaxStudRes, false, 1, false, false);

        double[] featureXvalsForRegression = new double[indexesLowStuff.length];
        double[] massDiffsForRegression = new double[indexesLowStuff.length];
        for (int i=0; i<indexesLowStuff.length; i++)
        {
            featureXvalsForRegression[i] = featureXvals[indexesLowStuff[i]];
            massDiffsForRegression[i] = massDiffs[indexesLowStuff[i]];
        }
        double[] resultCoefficients = null;
        resultCoefficients = RegressionUtilities.modalRegression(featureXvalsForRegression,
                massDiffsForRegression, calRegressionDoF);

        if (showCharts)
        {
            double minXval = BasicStatistics.min(featureXvals);
            double maxXval = BasicStatistics.max(featureXvals);

            PanelWithScatterPlot psp = new PanelWithScatterPlot();
            String chartName = "MassCal" + chartSuffix;
            psp.setName(chartName);
            psp.addLineOrCurve(resultCoefficients, minXval, maxXval);

            psp.addData(featureXvalsForRegression, massDiffsForRegression,
                    "Matches used in regression");
            psp.addData(featureXvals, massDiffs, "all mass matches");

            psp.setAxisLabels("mass","mass difference");
            psp.displayInTab();
        }
        return resultCoefficients;
    }



    public Map<Feature, List<ChemicalCompound>> massMatch(Feature[] features,
                                                           List<ChemicalCompound> massesListByMassAsc,
                                                           float massTolPPM)
    {
        int minPossibleMassIndex = 0;
        Map<Feature, List<ChemicalCompound>> result = new HashMap<Feature, List<ChemicalCompound>>();
        boolean newFeature = true;
        for (Feature feature :features)
        {

            float featureMass = feature.getMass();
            float massToleranceDa = MassUtilities.calculateAbsoluteDeltaMass(featureMass, massTolPPM,
                            FeatureSetMatcher.DELTA_MASS_TYPE_PPM);
            float minMatchMass = featureMass - massToleranceDa;
            float maxMatchMass = featureMass + massToleranceDa;


            while (massesListByMassAsc.get(minPossibleMassIndex).getMass() < minMatchMass &&
                   minPossibleMassIndex < massesListByMassAsc.size()-1)
            {
                minPossibleMassIndex++;
            }
            List<ChemicalCompound> massesMatchedThisFeature = new ArrayList<ChemicalCompound>();
            for (int i=minPossibleMassIndex ; i< massesListByMassAsc.size(); i++)
            {
                ChemicalCompound massAndName = massesListByMassAsc.get(i);
                float mass = (float) massAndName.getMass();

                if (mass < minMatchMass)
                    continue;
                else if (maxMatchMass < mass)
                    break;
                else
                {
                    float deltaMassPPM = (feature.getMass() - mass) * 1000000 / mass;
                    if (newFeature)
                    {
                        newFeature = false;
                    }

                    massesMatchedThisFeature.add(massAndName);
                }
            }
            if (massesMatchedThisFeature.size() >= 1)
            {
                result.put(feature, new ArrayList<ChemicalCompound>(massesMatchedThisFeature));
            }
        }
        return result;
    }

    public float getCalMassTolerancePPM()
    {
        return calMassTolerancePPM;
    }

    public void setCalMassTolerancePPM(float calMassTolerancePPM)
    {
        this.calMassTolerancePPM = calMassTolerancePPM;
    }

    public float getCalMaxLeverage()
    {
        return calMaxLeverage;
    }

    public void setCalMaxLeverage(float calMaxLeverage)
    {
        this.calMaxLeverage = calMaxLeverage;
    }

    public double getCalMaxStudRes()
    {
        return calMaxStudRes;
    }

    public void setCalMaxStudRes(double calMaxStudRes)
    {
        this.calMaxStudRes = calMaxStudRes;
    }

    public List<ChemicalCompound> getDatabaseCompoundsByMass()
    {
        return databaseCompoundsByMass;
    }

    public void setDatabaseCompoundsByMass(List<ChemicalCompound> databaseCompoundsByMass)
    {
        this.databaseCompoundsByMass = databaseCompoundsByMass;
    }

    public boolean isShowCharts()
    {
        return showCharts;
    }

    public void setShowCharts(boolean showCharts)
    {
        this.showCharts = showCharts;
    }

    public int getCalRegressionDoF() {
        return calRegressionDoF;
    }

    public void setCalRegressionDoF(int calRegressionDoF) {
        this.calRegressionDoF = calRegressionDoF;
    }
}
