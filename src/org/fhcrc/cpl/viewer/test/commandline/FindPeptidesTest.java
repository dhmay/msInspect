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

package org.fhcrc.cpl.viewer.test.commandline;

import junit.framework.Test;
import junit.framework.TestResult;
import junit.framework.AssertionFailedError;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.test.TestUtilities;
import org.fhcrc.cpl.toolbox.test.commandline.BaseCommandLineTest;
import org.fhcrc.cpl.viewer.commandline.modules.FindPeptidesCommandLineModule;

import java.io.File;

public class FindPeptidesTest extends BaseCommandLineTest implements Test
{
    //MD5 sums of files known to be a good result of a test run
    protected static final String sumStart50Count50Accurate = "b246515ec84ef03caedcd0b00e2894c8";
    protected static final String sumStart50Count50NoAccurate = "dea54286edac89d5f8da54a20b7bde70";
    protected static final String sumStrategyCentroidedStart200Count50 = "563e9afc2b024294946271ddc55561cb";
    protected static final String sumStrategyWaveletStart100Count5MinMz700MaxMz800 = "e6b19c1790383f60644697501ac9626f";

    /**
     * Return an instance of the FindPeptidesCommandLineModule
     * @return
     */
    protected CommandLineModule createCommandLineModule()
    {
        return new FindPeptidesCommandLineModule();
    }

    protected void doCleanup()
    {
        //no special cleanup for this test yet
    }

    /**
     * Strategy is to filter a set of features in various ways and then compare MD5 sums of
     * the resulting feature files (with whitespace removed) with sums computed when the tests
     * were written
     * @param result
     */
    protected void doRun(TestResult result)
    {
        try
        {
            //default strategy is PeakClusters
            String outFilePath = constructTempFilePath("findpeptidestest.out");
            setArgumentValue("out", outFilePath);
            setArgumentValue("start", "50");
            setArgumentValue("count", "100");
            setUnnamedSeriesArgumentValue(TestUtilities.getSampleDataDirFilePath("sample.mzXML"));

            log("Testing default strategy, 50 scans...");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumStart50Count50Accurate));

            renewCommandLineModule();

            setArgumentValue("start", "50");
            setArgumentValue("count", "100");
            setArgumentValue("noaccuratemass", "true");
            setUnnamedSeriesArgumentValue(TestUtilities.getSampleDataDirFilePath("sample.mzXML"));

            log("Testing default strategy without accurate mass adjustment, 50 scans...");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumStart50Count50NoAccurate));

            renewCommandLineModule();

            log("Testing Centroided strategy, 50 scans...");
            setArgumentValue("strategy", "FeatureStrategyCentroided");
            setArgumentValue("start", "200");
            setArgumentValue("count", "50");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumStrategyCentroidedStart200Count50));

            renewCommandLineModule();

            log("Testing Wavelet strategy, 5 scans, mz 700-800");
            setArgumentValue("strategy", "FeatureStrategyWavelet");
            setArgumentValue("start", "100");
            setArgumentValue("count", "5");
            setArgumentValue("minmz", "700");
            setArgumentValue("maxmz", "800");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumStrategyWaveletStart100Count5MinMz700MaxMz800));

            //...
        }
        catch (AssertionFailedError afe)
        {
            result.addError(this, afe);
        }
        catch (Exception e)
        {
            result.addError(this, e);
        }
    }

}
