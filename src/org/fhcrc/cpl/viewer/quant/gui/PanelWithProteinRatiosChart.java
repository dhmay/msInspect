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

package org.fhcrc.cpl.viewer.quant.gui;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.IsotopicLabelExtraInfoDef;
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
