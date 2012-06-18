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
package org.fhcrc.cpl.toolbox;

import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;

import java.io.IOException;
import java.io.File;
import java.io.FileOutputStream;
import java.net.URL;
import java.lang.reflect.Method;

/**
 * This class exists to make it easier to open browser windows,
 * and to provide various kinds of URL construction
 */
public class BrowserController
{
    /**
     * @param url
     * @throws IOException if the attempt to open the browser fails
     */
    public static void navigate(URL url) throws IOException
    {
       navigate(url.toString());
    }

    /**
     * throws IOException 
     * @param urlString
     * @throws IOException if the attempt to open the browser fails
     */
    public static void navigate(String urlString) throws IOException
    {
        if (urlString == null)
            return;

        String osName = String.valueOf(String.valueOf(System.getProperty("os.name")));
        //If we're on Windows, great, there's a standard way of opening a browser window.
        //Otherwise, take a shot in the dark
        if (osName.toLowerCase().contains("windows"))
            Runtime.getRuntime().exec(new String []{"rundll32", "url.dll,FileProtocolHandler", urlString});
        else if (osName.startsWith("Mac OS"))
        {
            try
            {
                Class fileMgr = Class.forName("com.apple.eio.FileManager");
                Method openURL = fileMgr.getDeclaredMethod("openURL",
                        new Class[] {String.class});
                openURL.invoke(null, new Object[] {urlString});
            }
            catch (Exception e)
            {
                throw new IOException("Failed to open Mac browser.  Message: " + e.getMessage());
            }
        }
        else Runtime.getRuntime().exec(new String []{"firefox", urlString}); // might as well try
    }

    /**
     *
     * @param file
     * @throws IOException if the attempt to open the browser fails
     */
    public static void navigate(File file) throws IOException
    {
        navigate("file://" + file.getAbsolutePath());
    }
    
    /**
     * Write the specified contents to a temp file and point the browser at it.  Needs
     * an object to tie the temp file to, so we know when to clean it up
     * @param contents
     * @param caller
     * @throws IOException
     */
    public static void openTempFileWithContents(String contents, String tempFileName,
                                                Object caller) throws IOException
    {
        File upregulatedPeptidesHtmlFile =
                TempFileManager.createTempFile(tempFileName, caller);
        FileOutputStream fos = new FileOutputStream(upregulatedPeptidesHtmlFile);
        fos.write(contents.getBytes());

        navigate(upregulatedPeptidesHtmlFile);
    }

    /**
     * Attempts to point the browser to a temp file filled with "contents".  If that fails,
     * opens an HtmlViewerPanel on the file, with title tempFileName.  Output a warning if browser failed
     * @param contents
     * @param tempFileName
     * @param caller
     * @throws IOException only if we can't write the temp file
     */
    public static void navigateOrPanelTempFileWithContents(String contents, String tempFileName,
                                                Object caller) throws IOException
    {
        File upregulatedPeptidesHtmlFile =
                TempFileManager.createTempFile(tempFileName, caller);
        FileOutputStream fos = new FileOutputStream(upregulatedPeptidesHtmlFile);
        fos.write(contents.getBytes());
        try
        {
            navigate(upregulatedPeptidesHtmlFile);
        }
        catch (IOException e)
        {
            ApplicationContext.infoMessage("WARNING: Failed to open browser, using Swing panel.  Error message: " +
                    e.getMessage());
            HtmlViewerPanel.showFileInDialog(upregulatedPeptidesHtmlFile, tempFileName);
        }
    }


}

