package org.fhcrc.cpl.viewer.feature.extraction;

import org.fhcrc.cpl.viewer.feature.Spectrum;

/**
 * Removes background noise from spectra
 */
public class BackgroundRemover
{
    protected int _resamplingFrequency = 36;

    public float[][] removeBackground(float[][] spectra)
    {
        return Spectrum.RemoveBackground(spectra);
    }

    /**
     * Calculate median intensity at each point on the grid
     * @param spectra
     * @return
     */
    public float[][] calculateMedian(float[][] spectra)
    {
        int numSpectra = spectra.length;
        int spectrumHeight = spectra[0].length;

        float[][] median = new float[numSpectra][];
        for (int i = 0; i < numSpectra; i++)
            median[i] = Spectrum.MedianWindow(spectra[i], spectrumHeight, 2 * _resamplingFrequency, false);
        float[] row = null;
        for (int r = 0; r < spectrumHeight; r++)
        {
            row = Spectrum.getRow(spectra, r, row);
            float[] m = Spectrum.MedianWindow(row, numSpectra, _resamplingFrequency, false);
            for (int s = 0; s < m.length; s++)
                median[s][r] = Math.max(median[s][r], m[s]);
        }
        return median;
    }


    public int getResamplingFrequency()
    {
        return _resamplingFrequency;
    }

    public void setResamplingFrequency(int resamplingFrequency)
    {
        this._resamplingFrequency = resamplingFrequency;
    }
}