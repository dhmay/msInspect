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
package org.fhcrc.cpl.toolbox;

import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;

import java.io.IOException;
import java.io.File;
import java.io.FileOutputStream;
import java.net.URL;

/**
 * This class exists to make it easier to open browser windows,
 * and to provide various kinds of URL construction
 */
public class BrowserController
{
    public static void navigate(URL url) throws IOException
    {
       navigate(url.toString());
    }

    public static void navigate(String urlString) throws IOException
    {
        if (urlString == null)
            return;

        String os = String.valueOf(String.valueOf(System.getProperty("os.name")));
        //If we're on Windows, great, there's a standard way of opening a browser window.
        //Otherwise, take a shot in the dark
        if (os.toLowerCase().contains("windows"))
            Runtime.getRuntime().exec(new String []{"rundll32", "url.dll,FileProtocolHandler", urlString});
        else
            Runtime.getRuntime().exec(new String []{"firefox", urlString}); // might as well try
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

        navigate("file://" + upregulatedPeptidesHtmlFile.getAbsolutePath());

    }


}

