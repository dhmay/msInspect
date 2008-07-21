/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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
import org.fhcrc.cpl.viewer.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.modules.FilterFeaturesCommandLineModule;
import org.fhcrc.cpl.viewer.test.TestUtilities;

import java.io.File;

public class FilterFeaturesTest extends BaseCommandLineTest implements Test
{
    //MD5 sums of files known to be a good result of a test run
    protected static final String sumMinCharge1MaxCharge1MaxKl2 = "d19e82f4aed74395a6ecb4591b2b6b8a";
    protected static final String sumScanFirst100ScanLast200 = "7f23cd68af1f36b16aaa98b9ad21eb33";
    protected static final String sumMinTime1500MaxTime1550 = "474510eb884decc6acdf58f97927b15a";
    protected static final String sumMinIntensity15 = "fb23cc65fb51d63a788e2e5e7d93ff43";
    protected static final String sumMinTotalIntensity1000 = "3f9216a2dbd00ee0cc20740be5c1fff4";
    protected static final String sumMinPeaks3 = "2855bdae748b26923a89faf75103370e";

    /**
     * Return an instance of the FilterFeaturesCommandLineModule
     * @return
     */
    protected CommandLineModule createCommandLineModule()
    {
        return new FilterFeaturesCommandLineModule();
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
            String outFilePath = constructTempFilePath("filtertest.out");
            setArgumentValue("out", outFilePath);
            setUnnamedSeriesArgumentValue(TestUtilities.getSampleDataDirFilePath("sample.peptides.tsv"));


            log("Testing charge and kl filtering...");
            setArgumentValue("mincharge", "1");
            setArgumentValue("maxcharge", "1");
            setArgumentValue("maxkl","2");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumMinCharge1MaxCharge1MaxKl2));

            unsetArgumentValue("mincharge");
            unsetArgumentValue("maxcharge");
            unsetArgumentValue("maxkl");
            renewCommandLineModule();


            log("Testing first/last scan filtering...");
            setArgumentValue("scanfirst", "100");
            setArgumentValue("scanlast", "200");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumScanFirst100ScanLast200));

            unsetArgumentValue("scanfirst");
            unsetArgumentValue("scanlast");
            renewCommandLineModule();


            log("Testing min/max time filtering...");
            setArgumentValue("mintime", "1500");
            setArgumentValue("maxtime", "1550");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                        sumMinTime1500MaxTime1550));

            unsetArgumentValue("mintime");
            unsetArgumentValue("maxtime");
            renewCommandLineModule();

            log("Testing intensity filtering...");
            setArgumentValue("minintensity", "15");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                       sumMinIntensity15));

            unsetArgumentValue("minintensity");
            renewCommandLineModule();

            log("Testing total intensity filtering...");
            setArgumentValue("mintotalintensity", "1000");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                       sumMinTotalIntensity1000));

            unsetArgumentValue("mintotalintensity");
            renewCommandLineModule();

            log("Testing min peaks filtering...");
            setArgumentValue("minpeaks", "3");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(new File(outFilePath),
                       sumMinPeaks3));

//            unsetArgumentValue("minpeaks");
//            renewCommandLineModule();

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