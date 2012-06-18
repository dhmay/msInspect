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

import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.text.html.HTMLFrameHyperlinkEvent;
import javax.swing.text.html.HTMLDocument;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.HyperlinkEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.awt.*;

/**
 * Display HTML in a scrollpane
 */
public class HtmlViewerPanel extends JScrollPane
{
    protected static Logger _log = Logger.getLogger(HtmlViewerPanel.class);


    protected JEditorPane editorPane = null;

    public static final int DEFAULT_WIDTH=1000;

    public static final int DEFAULT_HEIGHT=800;

    public HtmlViewerPanel()
    {
        editorPane = new JEditorPane();
        editorPane.setContentType("text/html");
        editorPane.setEditable(false);
        editorPane.addHyperlinkListener(new MyHyperLinkListener());

        setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        setViewportView(editorPane);
    }

    public HtmlViewerPanel(File file)
            throws IOException
    {
        this();
        navigate(file);
    }

    public HtmlViewerPanel(URL url)
            throws IOException
    {
        this();
        navigate(url);
    }

    public static JDialog showURLInDialog(String url, String dialogTitle)
            throws IOException
    {
        return showURLInDialog(new URL(url), dialogTitle);
    }

    public static JDialog showURLInDialog(URL url, String dialogTitle)
            throws IOException
    {
        JDialog dialog = new JDialog();
        dialog.setSize(DEFAULT_WIDTH,DEFAULT_HEIGHT);
        dialog.add(new HtmlViewerPanel(url));
        dialog.setTitle(dialogTitle);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        dialog.setVisible(true);
        return dialog;
    }

    public static JDialog showHTMLInDialog(String html, String dialogTitle)
            throws IOException
    {
        JDialog dialog = new JDialog();
        dialog.setSize(DEFAULT_WIDTH,DEFAULT_HEIGHT);
        HtmlViewerPanel htmlPanel = new HtmlViewerPanel();
        htmlPanel.displayHTML(html);
        dialog.add(htmlPanel);
        dialog.setTitle(dialogTitle);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        dialog.setVisible(true);
        return dialog;
    }    

    public static JDialog showResourceInDialog(String resourceLocation, String dialogTitle)
            throws IOException
    {
        ClassLoader classLoader = HtmlViewerPanel.class.getClassLoader();
        URL url = classLoader.getResource(resourceLocation);
        return showURLInDialog(url, dialogTitle);
    }

    public static JDialog showFileInDialog(File file, String dialogTitle)
            throws IOException
    {
        return showURLInDialog(file.toURL(), dialogTitle);
    }


    public HtmlViewerPanel(String urlString)
            throws IOException
    {
        this();
        navigate(urlString);
    }


    public void navigate(String urlString)
            throws IOException
    {
        navigate(new URL(urlString));
    }

    public void navigate(File file)
            throws IOException
    {
        navigate(file.toURL());
    }

    public void navigate(URL url)
            throws IOException
    {
        editorPane.setPage(url);
    }

    public void displayHTML(String html)
    {
        editorPane.setText(html);
    }



    protected class MyHyperLinkListener implements HyperlinkListener
    {
        public void hyperlinkUpdate(HyperlinkEvent e) {
            if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                JEditorPane pane = (JEditorPane) e.getSource();
                if (e instanceof HTMLFrameHyperlinkEvent) {
                    HTMLFrameHyperlinkEvent evt = (HTMLFrameHyperlinkEvent)e;
                    HTMLDocument doc = (HTMLDocument)pane.getDocument();
                    doc.processHTMLFrameHyperlinkEvent(evt);
                } else {
                    try {
                        pane.setPage(e.getURL());
                    } catch (Throwable t) {
                        t.printStackTrace();
                    }
                }
            }
        }
    }


    public JEditorPane getEditorPane()
    {
        return editorPane;
    }

    public void setEditorPane(JEditorPane editorPane)
    {
        this.editorPane = editorPane;
    }
}
