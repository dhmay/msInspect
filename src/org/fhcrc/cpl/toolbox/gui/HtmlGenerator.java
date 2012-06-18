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

package org.fhcrc.cpl.toolbox.gui;

import javax.swing.*;
import java.util.List;
import java.util.ArrayList;

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
     * @return
     */
    public static String createHtmlStringForTable(String[] columnHeaders, Object[][] tableData)
    {
        StringBuffer htmlToWrite = new StringBuffer();
        htmlToWrite.append("<table>");
        if (columnHeaders != null)
        {
            htmlToWrite.append("<tr>");
            for (String columnHeader : columnHeaders)
                htmlToWrite.append("<th>" + columnHeader + "</td>");
            htmlToWrite.append("</tr>");

        }
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

    public static String createTableRow(List<String> rowValues)
    {
        StringBuffer htmlToWrite = new StringBuffer();
        htmlToWrite.append("<tr>");

        for (String rowValue : rowValues)
        {
            htmlToWrite.append("<td>" + rowValue + "</td>");
        }
        htmlToWrite.append("</tr>");

        return htmlToWrite.toString();
    }

    /**
     * Create an HTML table from the object matrix, using toString() to populate the cells
     * @param tableData
     * @return
     */
    public static String createHtmlStringForTable(Object[][] tableData)
    {
        return createHtmlStringForTable(null, tableData);
    }

    public static String createLink(String url, String text)
    {
        return "<a href=\"" + url + "\">" + text + "</a>";
    }

    /**
     *
     * @param ordered
     * @param listTexts
     * @param urls
     * @return
     */
    public static String createListOfLinks(boolean ordered, List<String> listTexts, List<String> urls)
    {
        List<String> linksAndValues = new ArrayList<String>(listTexts.size());
        for (int i=0; i<listTexts.size(); i++)
        {
            linksAndValues.add(createLink(urls.get(i), listTexts.get(i)));
        }
        return createList(ordered, linksAndValues);
    }

    public static String createList(boolean ordered, List<String> listValues)
    {
        String listTag = ordered ? "ol" : "ul";

        StringBuffer result = new StringBuffer();
        result.append("<" + listTag + ">\n");
        for (String listValue : listValues)
            result.append("\t" + "<li>" + listValue + "</li>\n");
        result.append("</" + listTag + ">\n");

        return result.toString();
    }
}
