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

package org.fhcrc.cpl.viewer.test;

import junit.framework.TestFailure;
import junit.framework.TestResult;
import junit.framework.TestSuite;

import java.util.*;
import java.io.File;

import org.fhcrc.cpl.viewer.util.ConvertHelper;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.test.BaseCommandTest;
import org.fhcrc.cpl.toolbox.test.TestUtilities;

public class TestRunner extends TestSuite
{
    protected static String _testPackage = "org.fhcrc.cpl.viewer.test";
    protected static String _testBundleName = "tests";



    // Fire the static converter initializers needed by the tab loader
    static
    {
        ConvertHelper.registerHelpers();
    }

    /**
     * Run tests.  If the "test" system property is specified, try to figure out which test
     * it refers to and run that one.  Otherwise run all known tests
     * @param args
     */
    public static void main(String args[])
    {
        //This is a helper utility to calculate md5 sums of files, for comparison
        if (args.length > 0 && "md5sum".equalsIgnoreCase(args[0]))
        {
            File fileToSum = new File(args[1]);
            if (!fileToSum.exists() || !fileToSum.canRead())
                log("Error reading file " + fileToSum.getAbsolutePath());
            else
            {
                ApplicationContext.infoMessage("Calculating MD5 sum...");
                try
                {
                    log(TestUtilities.calculateHexMD5SumNoWhiteSpaceNoComments(fileToSum));
                }
                catch(Exception e)
                {
                    log("Failed to calculate MD5 sum: " + e.getMessage() );
                    e.printStackTrace(System.err);
                }
            }
            return;
        }

        //check for definition of msInspect source root directory
        String toolsRootDir = System.getProperty("viewer.root");
        if (toolsRootDir == null || toolsRootDir.equals(""))
        {
            log("Property viewer.root must be set in order to run tests");
            return;
        }

        Class[] standardTestClasses = loadStandardTestClasses();

        ArrayList<BaseCommandTest> tests = new ArrayList<BaseCommandTest>();
        String testName = System.getProperty("test");

        if (testName != null && testName.length() > 0)
        {
            //Test was specified, try to instantiate it
            for (Class knownClass : standardTestClasses)
            {
                if (knownClass.getSimpleName().equalsIgnoreCase(testName) ||
                    knownClass.getSimpleName().equalsIgnoreCase(testName + "Test"))
                {
                try
                {
                    Object testObject = knownClass.newInstance();
                    if (testObject instanceof BaseCommandTest)
                        tests.add((BaseCommandTest) testObject);
                    else
                        log("Failure adding test class " + knownClass.getName() + ": not an instance of BaseViwerTest");
                }
                catch (Exception e)
                {
                    log("Failure adding test class " + knownClass.getName());
                }
                    break;
                }
            }
            if (tests.size() == 0)
                log("Couldn't recognize test " + testName);
        }
        else //no test specified, run them all
        {
            for (Class testClass : standardTestClasses)
            {
                try
                {
                    Object testObject = testClass.newInstance();
                    if (testObject instanceof BaseCommandTest)
                        tests.add((BaseCommandTest) testObject);
                    else
                        log("Failure adding test class " + testClass.getName() + ": not an instance of BaseViwerTest");
                }
                catch (Exception e)
                {
                    log("Failure adding test class " + testClass.getName());
                }
            }
        }

        if (tests.size() == 0)
        {
            log("No tests specified");
            return;
        }

        if (!cleanUpTests(tests))
            return;
        log("========= Done cleaning up ========= \n");

        if("true".equals(System.getProperty("cleanonly")))
            return;

        ArrayList<String> successClassNames = new ArrayList<String>();
        log("========= Running tests... ========= ");
        for (BaseCommandTest test : tests)
        {
            log("========= Running Test " + test.getClass().getSimpleName() + " =========");
            TestResult result = test.run();
            if (result.errorCount() == 0 && result.failureCount() == 0)
            {
                log("========= " + test.getClass().getSimpleName() + " test success =========");
                successClassNames.add(test.getClass().getSimpleName());
            }
            else
            {
                Enumeration errors = result.errors();
                while (errors.hasMoreElements())
                {
                    System.err.println("Error:");
                    TestFailure failure = (TestFailure) errors.nextElement();
                    if (failure.exceptionMessage() != null)
                        System.err.println(failure.exceptionMessage());
                    System.err.println(failure.trace());
                }

                log("Test " + test.getClass().getSimpleName() + " failed, " + result.errorCount() + " errors.");
                log("==================================================================");
                log("Test Suite Failed.  ");
                log(successClassNames.size() + " tests succeeded");
                for (String successTestName : successClassNames)
                    log("\t" + successTestName);
                log("Test " + test.getClass().getSimpleName() + " failed.");
                System.exit(result.errorCount() + result.failureCount());
            }
        }
        log("All " + successClassNames.size() + " tests passed");
        closeLog();
    }

    /**
     * Clean up after all tests
     * @param tests
     * @return
     */
    public static boolean cleanUpTests(ArrayList<BaseCommandTest> tests)
    {
        //Clean up after all tests.  At this level, we know we want to depopulate the
        //scratch directories
        log("========= Cleaning up after all tests... =========");
        for (BaseCommandTest test : tests)
        {
            log("Cleaning up test " + test.getClass().getSimpleName());
            try
            {
                test.cleanup();
            }
            catch (Exception e)
            {
                log("Exception when cleaning up prior to executing test:");
                log(e.getMessage());
                e.printStackTrace(System.err);
                log("Please clean up manually before trying again.");
                return false;
            }
        }
        return true;
    }

    /**
     * Load the standard test classes specified in the resource bundle
     * @return
     */
    protected static Class[] loadStandardTestClasses()
    {
        List<Class> testClassList = new ArrayList<Class>();

        ResourceBundle testResourceBundle =
                ResourceBundle.getBundle(_testPackage + "." +
                        _testBundleName);
        Enumeration propertyNames = testResourceBundle.getKeys();

        String propertyName = null;
        while(propertyNames.hasMoreElements())
        {

            propertyName = (String) propertyNames.nextElement();
            if ("true".equalsIgnoreCase(testResourceBundle.getString(propertyName)))
            {
                Class testClass = null;
                try
                {
                    testClass = Class.forName(_testPackage + "." + propertyName);
                }
                catch (Exception e)
                {
                    try
                    {
                        //todo: what if someone wants a test in another package?
                        testClass = Class.forName(_testPackage + ".commandline." + propertyName);
                    }
                    catch (Exception x)
                    {
                        ApplicationContext.errorMessage("Exception loading test module with name " +
                                                        propertyName, x);
                    }
                }

                if (testClass != null)
                {
                    
                    testClassList.add(testClass);
                }
            }
        }
        return testClassList.toArray(new Class[testClassList.size()]);
    }

    /**
     * Write a lot message
     * @param mesg
     */
    protected static void log(String mesg)
    {
        TestUtilities.log(mesg);
    }

    /**
     * Close the log file
     */
    protected static void closeLog()
    {
        TestUtilities.closeLog();
    }
}
