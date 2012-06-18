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
import org.fhcrc.cpl.viewer.commandline.modules.FilterFeaturesCommandLineModule;
import org.fhcrc.cpl.toolbox.test.TestUtilities;
import org.fhcrc.cpl.toolbox.test.commandline.BaseCommandLineTest;

import java.io.File;

public class FilterFeaturesTest extends BaseCommandLineTest implements Test
{
    //MD5 sums of files known to be a good result of a test run
    protected static final String sumMinCharge1MaxCharge1MaxKl2 = "57f9ae1925ba11d9fd14ba72576bdef8";
    protected static final String sumScanFirst100ScanLast200 = "824e427d2eeec3e2b2710f8403f9f616";
    protected static final String sumMinTime1500MaxTime1550 = "4430f8e8460f67d640fe3f2f55121f68";
    protected static final String sumMinIntensity15 = "b0c7e567a0e499eff74b3a96de775cd4";
    protected static final String sumMinTotalIntensity1000 = "b0ccdbece4905be3ffa5bdb9be2b35a4";
    protected static final String sumMinPeaks3 = "4742db6aefc57223907022818a2d75e5";

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
