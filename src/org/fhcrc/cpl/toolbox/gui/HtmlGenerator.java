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

package org.fhcrc.cpl.toolbox.gui;

import javax.swing.*;

/**
 * Helpful HTML generation.  Keep this simple, no need to reinvent wheels -- this stuff
 * is purely for convenience.  Anything that gets too complex should be done in Swing
 */
public class HtmlGenerator extends JPanel
{
    protected static final String HTML_OPEN_TAG = "<html>";
    protected static final String HTML_CLOSE_TAG = "</html>";

    public static String createDocumentHeader(String title)
    {
        return "<head><title>" + title + "</head></title>";
    }

    public static String createDocumentBeginning(String title)
    {
        return HTML_OPEN_TAG + createDocumentHeader(title) + "\n<body>";
    }

    public static String createDocumentEnd()
    {
        return "</body>\n" + HTML_CLOSE_TAG + "\n";
    }

    /**
     * Cover method for createHtmlDocumentForTable, for a table with just one column
     * @param tableDataSingleColumn
     * @param title
     * @return
     */
    public static String createHtmlDocumentForTableOneColumn(
            Object[] tableDataSingleColumn, String title)
    {
        Object[][] tableData = new Object[tableDataSingleColumn.length][1];
        for (int i=0; i<tableDataSingleColumn.length; i++)
            tableData[i][0] = tableDataSingleColumn[i];
        return createHtmlDocumentForTable(tableData, title);
    }

    /**
     * Create an HTML document with the given body and title
     * @param body
     * @param title
     * @return
     */
    public static String createHtmlDocument(String body, String title)
    {
        return createDocumentBeginning(title) + body +
                createDocumentEnd();
    }

    /**
     * Create an entire HTML document containing a table created from the object matrix
     * @param tableData
     * @param title
     * @return
     */
    public static String createHtmlDocumentForTable(Object[][] tableData, String title)
    {
        return createHtmlDocument(createHtmlStringForTable(tableData),title);
    }

    /**
     * Create an HTML table from the object matrix, using toString() to populate the cells
     * TODO: table headers
     * @param tableData
     * @param title
     * @return
     */
    public static String createHtmlStringForTable(Object[][] tableData)
    {
        StringBuffer htmlToWrite = new StringBuffer();
        htmlToWrite.append("<table>");
        for (Object[] row : tableData)
        {
            htmlToWrite.append("<tr>");
            for (Object cell : row)
                htmlToWrite.append("<td>" + cell + "</td>");
            htmlToWrite.append("</tr>");
        }
        htmlToWrite.append("</table>");
        return htmlToWrite.toString();
    }
}
