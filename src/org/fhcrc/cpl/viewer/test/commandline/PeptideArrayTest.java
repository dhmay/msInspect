/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.test.TestUtilities;
import org.fhcrc.cpl.toolbox.test.commandline.BaseCommandLineTest;
import org.fhcrc.cpl.viewer.align.commandline.PeptideArrayCommandLineModule;

import java.io.File;

/**
 * Note:  since this test calls out to R, it is dependent on the version of R on the path.
 * Currently, the MD5 sums are calculated using R 2.0.1
 */
public class PeptideArrayTest extends BaseCommandLineTest implements Test
{
    //MD5 sums of files known to be a good result of a test run

    //msinspect --peptidearray --normalize --masswindow=.1 --scanwindow=75 --out=/home/dhmay/temp/pepArrayTest.out ms1features_1.tsv ms1features_2.tsv
    protected static final String sumMassPoint1Scan75Normalize =
            "079ca865079c098b85e0f8db04bbcb8c";

    //msinspect --peptidearray --normalize --masswindow=.2 --scanwindow=50 --out=/home/dhmay/temp/pepArrayTest.out ms1features_1.tsv ms1features_2.tsv ms1features_3.tsv
    protected static final String sumMassPoint2Scan50 =
            "479ed449eeabd7902ce05ad1bfb8a507";
    protected static final String sumMassPoint2Scan50Details =
            "b447a41af5259e162c60ad46b3f27a4f";



    /**
     * Return an instance of the FindPeptidesCommandLineModule
     * @return
     */
    protected CommandLineModule createCommandLineModule()
    {
        return new PeptideArrayCommandLineModule();
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
            File outFile = new File(constructTempFilePath("peptidearraytest.out"));
            //details filename is automatically generated by the code
            File outDetailsFile = new File(constructTempFilePath("peptidearraytest.details.out"));
            String outFilePath = outFile.getAbsolutePath();

//System.err.println("This test temporarily disabled until we can get Hertz onto R 2.1");
            setArgumentValue("out", outFilePath);
            setArgumentValue("masswindow", ".1");
            setArgumentValue("scanwindow", "75");
            setArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                    TestUtilities.getSampleDataDirFilePath("ms1features_1.tsv") +
                            CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR +
                    TestUtilities.getSampleDataDirFilePath("ms1features_2.tsv"));
            setArgumentValue("normalize","true");
            log("Testing aligning 2 files, masswindow .1, scanwindow 75, normalizing...");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(outFile,
                       sumMassPoint1Scan75Normalize));

            renewCommandLineModule();

            setArgumentValue("out", outFilePath);
            setArgumentValue("masswindow", ".2");
            setArgumentValue("scanwindow", "50");
            String filesArgValue =  TestUtilities.getSampleDataDirFilePath("ms1features_1.tsv") +
                            CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR +
                    TestUtilities.getSampleDataDirFilePath("ms1features_2.tsv") +
                            CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR +
                    TestUtilities.getSampleDataDirFilePath("ms1features_3.tsv");
            setArgumentValue(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                    filesArgValue);

            log("Testing aligning 3 files, masswindow .2, scanwindow 50...");
            assertTrue(executeCommandLineModule(result));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(outFile,
                       sumMassPoint2Scan50));
            assertTrue(TestUtilities.compareHexMD5SumsNoWhiteSpaceNoComments(outDetailsFile,
                       sumMassPoint2Scan50Details));
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
