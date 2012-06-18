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
package org.fhcrc.cpl.toolbox.filehandler;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;

/**
 * manages creation and destruction of temp files.
 *
 * Different classes request the creation of temp files, and TempFileManager records the class that
 * requested the files.  When the calling class is done with all temp files, it can call
 * deleteTempFiles to remove all those files.
 *
 * If a calling class is using static methods to create temp files, it will need to create a fake caller
 * object to associate the temp files with
 */
public class TempFileManager
{
    protected static File _tmpDir;
    private static Logger _log = Logger.getLogger(TempFileManager.class);
    protected static Map<Object, List<File>> objectTempFileMap = new HashMap<Object, List<File>>();

    /**
     * find the temp dir.  Null if errors
     * @return
     */
    public static File getTmpDir()
            throws RuntimeException
    {
        if (null != _tmpDir)
            return _tmpDir;

        File tmpDir = null;
        try
        {
            File tmpFile = File.createTempFile("msInspect", "");
            tmpFile.deleteOnExit();
            tmpDir = new File(tmpFile.getAbsolutePath() + "dir");
            if (!tmpDir.exists())
                tmpDir.mkdir();
        }
        catch (IOException x)
        {
            throw new RuntimeException("Couldn't create temp dir " + (tmpDir !=null ? tmpDir.getAbsolutePath() : ""));
        }

        tmpDir.deleteOnExit();
        _tmpDir = tmpDir;
        return _tmpDir;
    }

    /**
     * return a file in the temp dir.  Throw an exception if doesn't exist
     * @param fileName
     * @return
     */
    public static File getTmpFile(String fileName)
            throws RuntimeException
    {
        File tmpDir = getTmpDir();
        File tmpFile = null;

        tmpFile = new File(tmpDir, fileName);
        if (!tmpFile.exists())
            throw new RuntimeException("Requested temporary file " + tmpFile.getAbsolutePath() + " does not exist");


        return tmpFile;
    }

    /**
     * create a file in the temp dir.  If it already exists, so be it, return it.
     * Register the file with the list of files to be cleaned up when the caller asks for its
     * files to be cleaned up.
     * @param fileName
     * @return
     */
    public static File createTempFile(String fileName, Object caller)
            throws RuntimeException
    {
        File tempFile = new File(getTmpDir(), fileName);
        markFileForDeletion(tempFile, caller);
        return tempFile;
    }

    /**
     * mark this file to be deleted when the caller requests its filet to be deleted
     * @param file
     * @param caller
     */
    public static void markFileForDeletion(File file, Object caller)
    {
        List<File> tempFilesForCaller = getTempFiles(caller);
        if (tempFilesForCaller == null)
        {
            tempFilesForCaller = new ArrayList<File>();
            objectTempFileMap.put(caller, tempFilesForCaller);
        }
        tempFilesForCaller.add(file);
    }

    public static void unmarkFileForDeletion(File file, Object caller)
    {
        if (file == null)
            return;
        List<File> tempFilesForCaller = getTempFiles(caller);
        for (File tempFile : tempFilesForCaller)
        {
            try
            {
                if (tempFile.getCanonicalPath().equals(file.getCanonicalPath()))
                {
                    tempFilesForCaller.remove(tempFile);
                    break;                                                          
                }
            }
            catch (IOException e)
            {
                _log.debug("IOException comparing paths in unmarkFileForDeletion");
            }
        }
    }

    /**
     * Get the list of temp files registered by this caller
     * @param caller
     * @return
     */
    protected static List<File> getTempFiles(Object caller)
    {
        return objectTempFileMap.get(caller);
    }

    /**
     * Delete all registered temp files for this caller
     * @param caller
     */
    public static void deleteTempFiles(Object caller)
    {
        List<File> tempFilesForCaller = getTempFiles(caller);
        if (tempFilesForCaller != null)
        {
            for (File file : tempFilesForCaller)
            {
                file.delete();
            }
        }
        if (objectTempFileMap.containsKey(caller))
            objectTempFileMap.remove(caller);


//        //if there are no more temp files in the temp dir, clean it up.  It'll
//        //be recreated later if necessary
//        File[] tempFiles = getTmpDir().listFiles();
//        if (tempFiles != null && tempFiles.length == 0)
//        {
//            getTmpDir().delete();
//            _tmpDir = null;
//        }
    }

    /**
     * Get a filepath for a file in the tmp dir.  Null if problems
     * @param fileName
     * @return
     */
    public static String getTmpFilepath(String fileName)
    {
        File tmpFile = getTmpFile(fileName);
        if (tmpFile == null)
            return null;
        return tmpFile.getAbsolutePath();
    }

}
