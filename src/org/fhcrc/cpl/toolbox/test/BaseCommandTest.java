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

package org.fhcrc.cpl.toolbox.test;

import junit.framework.TestCase;
import junit.framework.TestResult;

import java.io.File;

import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;

/**
 * Abstract base class for all msInspect tests.  Lots of utility methods
 */
public abstract class BaseCommandTest extends TestCase
{
    public BaseCommandTest()
    {

    }

    /**
     * Call the doCleanup method in the test itself to do any necessary cleanup.
     * Then recursively delete all files in the test's temp dir, which is assigned
     * by class name.
     *
     * @throws Exception
     */
    public void cleanup() throws Exception
    {
        doCleanup();
        File testTempDir = new File(getTempDirName());
        if (testTempDir.exists())
            TestUtilities.recursiveFileDelete(testTempDir);
    }

    /**
     * Simple covers to get temp directory info for this class
     */

    protected String getTempDirName()
    {
        return TestUtilities.getCurrentTestTempDirName(getClass().getSimpleName());
    }
    protected String constructTempFilePath(String fileName)
    {
        return TestUtilities.constructTempFilePath(getClass().getSimpleName(), fileName);
    }

    /**
     * Each implementing class can do special cleanup here
     * @throws Exception
     */
    protected abstract void doCleanup() throws Exception;

    /**
     * Create a temp dir for the class to work in, then run the test
     * @param result
     */
    public void run(TestResult result)
    {
        if (!TestUtilities.createTestTempDir(getClass().getSimpleName()))
        {
            result.addError(this, new CommandLineModuleExecutionException("Failed to create temp directory " + getTempDirName() + ", quitting"));
            return;
        }
        doRun(result);
    }

    /**
     * Do the actual guts of the test
     * @param result
     */
    protected abstract void doRun(TestResult result);


    protected void log(String message)
    {
        TestUtilities.log(message);
    }
}
