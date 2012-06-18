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

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Date;

import org.fhcrc.cpl.toolbox.ApplicationContext;

/**
 * Abstract base class for all msInspect tests.  Lots of utility methods
 */
public class TestUtilities
{
    //Strings to store paths for log, temp files, etc
    protected static String mViewerRootDirName = null;
    protected static String mSampleDataDirName = null;
    protected static String mRootTestTempfileDirName = null;
    protected static String mRootTestDirName = null;
    protected static String mTestLogDirName = null;

    //Variables related to the log file
    protected static File mLogFile = null;
    protected static PrintWriter mLogPrintWriter = null;
    protected static boolean mTriedToOpenLog = false;

    public TestUtilities()
    {

    }

    /**
     * Recursively delete a file and any children from the filesystem
     * @param file
     */
    public static void recursiveFileDelete(File file)
    {
        if (file.isDirectory())
        {
            for (File childFile : file.listFiles())
            {
                recursiveFileDelete(childFile);
            }
        }
        file.delete();
    }

    /*
     * Methods to deal with temp files
     */


    /**
     * Create a temp dir for this test
     * @return
     */
    public static boolean createTestTempDir(String className)
    {
        try
        {
            File testTempDir = new File(getCurrentTestTempDirName(className));
            if (!testTempDir.exists())
            {
                if (!testTempDir.mkdirs())
                    return false;
            }
        }
        catch (Exception e)
        {
            return false;
        }
        return true;
    }

    /**
     * Get the directory name of the sampledata directory, where all files that the test might need
     * to work with should live
     * @return
     */
    public static String getSampleDataDirName()
    {
        if (mSampleDataDirName == null)
            mSampleDataDirName = getViewerRootDirName() + System.getProperty("file.separator") + "sampledata";
        return mSampleDataDirName;
    }

    /**
     * Convenience method to get the filepath of a named file in the sampledata dir
     * @param fileName
     * @return
     */
    public static String getSampleDataDirFilePath(String fileName)
    {
        return getSampleDataDirName() + System.getProperty("file.separator") + fileName;
    }

    /**
     * Get the viewer root directory (webapps/tools)
     * @return
     */
    public static String getViewerRootDirName()
    {
        if (mViewerRootDirName == null)
            mViewerRootDirName = System.getProperty("viewer.root");
        return mViewerRootDirName;
    }

    /**
     * webapps/tools/test.  Contains logs and temp_files subdirs
     * @return
     */
    public static String getRootTestDirName()
    {
        if (mRootTestDirName == null)
            mRootTestDirName = getViewerRootDirName() + System.getProperty("file.separator") + "test";
        return mRootTestDirName;
    }

    /**
     *
     * @return directory name of the test log directory
     */
    public static String getTestLogDirName()
    {
        if (mTestLogDirName == null)
            mTestLogDirName = getRootTestDirName() + System.getProperty("file.separator") + "logs";
        return mTestLogDirName;
    }

    /**
     *
     * @return dir name of the root test temp file directory
     */
    public static String getRootTestTempfileDirName()
    {
        if (mRootTestTempfileDirName == null)
            mRootTestTempfileDirName = getRootTestDirName() + System.getProperty("file.separator") + "temp_files";
        return mRootTestTempfileDirName;
    }

    /**
     *
     * @return dir name of the temp file dir for this test
     */
    public static String getCurrentTestTempDirName(String className)
    {
        return getRootTestTempfileDirName() + System.getProperty("file.separator") + className;
    }

    /**
     * Convenience method to construct a filepath in the temp file dir
     * @param relativeFilePath
     * @return
     */
    public static String constructTempFilePath(String className, String relativeFilePath)
    {
        return getCurrentTestTempDirName(className) + System.getProperty("file.separator") + relativeFilePath;
    }



    /*
     * Methods related to MD5 hashes.  We use MD5 hashes for results verification when we can be
     * absolutely sure that succeeding in a test requires the file to be exactly the same as a
     * file generated when the test is created, MODULO WHITESPACE.  We strip whitespace before
     * comparing the files
     *
     */


    /**
     * Convert a byte array to a hex string
     * @param v
     * @return
     */
    public static String toHexString (byte [] v)
    {
        StringBuffer	sb = new StringBuffer ();
        byte		n1, n2;

        for (int c = 0; c < v.length; c++)
        {
            n1 = (byte)((v[c] < 0 ? -v[c] + 127 : v[c]) / 0x10);
            n2 = (byte)((v[c] < 0 ? -v[c] + 127 : v[c]) % 0x10);

            sb.append (n1 >= 0xA ? (char)(n1 - 0xA + 'a') : (char)(n1 + '0'));
            sb.append (n2 >= 0xA ? (char)(n2 - 0xA + 'a') : (char)(n2 + '0'));
        }

        return sb.toString();
    }

    /**
     * Strip # style comments from a file.  Useful for creating md5 sums for files in which the
     * comments might vary but everything else must be the same
     * @param file
     * @throws Exception
     */
    public static void removeHashCommentLinesFromFile(File file) throws Exception
    {
        ArrayList<Byte> fileBytesList = new ArrayList<Byte>();

        FileInputStream fis = new FileInputStream(file);

        boolean isLineStart = true;
        boolean inComment = false;
        while (true)
        {
            byte nextChar = (byte) fis.read();
            if (nextChar == -1)
                break;
            if (inComment)
            {
                if (nextChar == '\n')
                {
                    inComment = false;
                    isLineStart = true;
                }
                continue;
            }
            if (isLineStart)
            {
                isLineStart = false;
                if (nextChar == '#')
                {
                    inComment = true;
                    continue;
                }
            }
            fileBytesList.add(nextChar);
        }
        fis.close();

        file.delete();
        FileOutputStream fos = new FileOutputStream(file);
        for (int i=0; i<fileBytesList.size(); i++)
        {
            fos.write(fileBytesList.get(i));
        }
        fos.flush();
        fos.close();
    }

    public static String calculateHexMD5SumNoWhiteSpaceNoComments(File file) throws Exception
    {
        removeHashCommentLinesFromFile(file);
        return calculateHexMD5SumNoWhiteSpace(file);
    }

    /**
     * Calculate an MD5 sum for a file, first stripping out all whitespace, then convert it to hex
     * @param file
     * @return
     * @throws Exception
     */
    public static String calculateHexMD5SumNoWhiteSpace(File file) throws Exception
    {
        ArrayList<Byte> fileBytesList = new ArrayList<Byte>();

        FileInputStream fis = new FileInputStream(file);

int realchars=0;
        while (true)
        {
            byte nextChar = (byte) fis.read();
            if (nextChar == -1)
                break;
            if (!Character.isWhitespace(nextChar))
            {
//if (Character.isLetterOrDigit(nextChar))
//{
//    System.err.println(nextChar);
//    realchars++;
//}
                fileBytesList.add(nextChar);
            }
        }
        MessageDigest md5Digest = MessageDigest.getInstance("MD5");
        byte[] fileBytes = new byte[fileBytesList.size()];
//System.err.println("stripped file is " + fileBytes.length + " bytes, with " + realchars + " letters and digits");
        for (int i=0; i<fileBytesList.size(); i++)
            fileBytes[i] = fileBytesList.get(i);
        return toHexString(md5Digest.digest(fileBytes));
    }

    /**
     * Compare the MD5 sum of a file with a previously stored sum
     * @param file
     * @param targetSum
     * @return
     * @throws Exception
     */
    public static boolean compareHexMD5SumsNoWhiteSpace(File file, String targetSum) throws Exception
    {
        String fileSum = calculateHexMD5SumNoWhiteSpace(file);
        if (fileSum.equals(targetSum))
            return true;
        //if we got here, failed
        log("Failed MD5 sum comparison for file " + file.getAbsolutePath());
        log("File exists? " + file.exists());
        log("Target sum: " + targetSum);
        log("Calculated sum: " + fileSum);
        return false;
    }

    public static boolean compareHexMD5SumsNoWhiteSpaceNoComments(File file, String targetSum) throws Exception
    {
        removeHashCommentLinesFromFile(file);
        return compareHexMD5SumsNoWhiteSpace(file,targetSum);
    }

    /**
     * Write a lot message to System.err and also to a log file
     * @param mesg
     */
    public static void log(String mesg)
    {
        Date date = new Date();
        //yeah, yeah, these methods are deprecated.  They'll be around forever.
        String minutes = "" + date.getMinutes();
        if (date.getMinutes() < 10) minutes = "0" + minutes;
        String seconds = "" + date.getSeconds();
        if (date.getSeconds() < 10) seconds = "0" + seconds;
        mesg = "[" + date.getHours() + ":" + minutes + ":" + seconds + "] " + mesg;
        ApplicationContext.infoMessage(mesg);
        if (mLogPrintWriter != null)
        {
            mLogPrintWriter.println(mesg);
            mLogPrintWriter.flush();
        }
        else
        {
            if (!mTriedToOpenLog)
            {
                openLogPrintWriter();
                if (mLogPrintWriter != null)
                    mLogPrintWriter.println(mesg);
            }
        }

    }

    /**
     * Manage the log file
     * @return
     */
    public static File getLogFile()
    {
        if (mLogFile == null)
        {
            File logDirectory = new File(getTestLogDirName());
            if (!logDirectory.exists())
                logDirectory.mkdirs();
            mLogFile = new File(getTestLogDirName() + System.getProperty("file.separator") +
                                "msinspect_test.log");
        }
        return mLogFile;
    }

    /**
     * Open a PrintWriter on the log file
     */
    public static void openLogPrintWriter()
    {
        if (mLogPrintWriter == null && !mTriedToOpenLog)
        {
            try
            {
                mTriedToOpenLog = true;
                mLogPrintWriter = new PrintWriter(getLogFile());
            }
            catch (Exception e)
            {
                System.err.println("Error opening log file for writing, logs will only be written to screen.  Error message: " +
                                   e.getMessage());
            }
        }
    }

    /**
     * Close the PrintWriter on the log file, if it's open
     */
    public static void closeLog()
    {
        if (mLogPrintWriter != null)
            mLogPrintWriter.close();
    }

}
