package org.fhcrc.cpl.viewer.quant.gui;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.viewer.feature.Feature;
import org.fhcrc.cpl.viewer.feature.extraInfo.IsotopicLabelExtraInfoDef;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithBarChart;

import java.util.List;


/**
 * Displays a bar chart showing the data that make up individual peptide ratios that contribute to
 * a protein ratio
 */
public class PanelWithProteinRatiosChart extends PanelWithBarChart
{
    protected static Logger _log = Logger.getLogger(PanelWithProteinRatiosChart.class);

    protected List<Feature> ratioPeptideFeatures;

    public PanelWithProteinRatiosChart()
    {
        super();
    }

    public PanelWithProteinRatiosChart(List<Feature> ratioPeptideFeatures)
    {
        this.ratioPeptideFeatures = ratioPeptideFeatures;

        double[] scaledIntensities = new double[ratioPeptideFeatures.size() * 2];
        String[] labels = new String[ratioPeptideFeatures.size() * 2];

        for (int i=0; i<ratioPeptideFeatures.size(); i++)
        {
            Feature feature = ratioPeptideFeatures.get(i);
            double lightArea = IsotopicLabelExtraInfoDef.getLightIntensity(feature);
            double heavyArea = IsotopicLabelExtraInfoDef.getHeavyIntensity(feature);

            double lightHeavySum = lightArea + heavyArea;
            double scaledLightArea = lightArea / lightHeavySum;
            double scaledHeavyArea = heavyArea / lightHeavySum;

            labels[i] = i + "_light";
            scaledIntensities[i] = scaledLightArea;

            labels[ratioPeptideFeatures.size() + i] = i + "_light";
            scaledIntensities[ratioPeptideFeatures.size() + i] = scaledLightArea;
        }
        addData(labels, scaledIntensities, "Ratios");
    }

}
