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
package org.fhcrc.cpl.toolbox.gui.widget;


import org.fhcrc.cpl.toolbox.gui.ImagePanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.net.URL;


/**
 * This is a JFrame that shows a splash image that disappears if the user presses a key or clicks on it.
 * To display this, call the setVisible() method.  To get rid of it, call dispose()
 */
public class SplashFrame extends JFrame
{
    protected Image img;

    public SplashFrame(Image img)
    {
        super();
        this.img = img;
        init();
    }

    /**
     * Read an image from a URL
     * @param url
     */
    public SplashFrame(URL url)
    {
        this(Toolkit.getDefaultToolkit().getImage(url));
        init();
    }

    protected void init()
    {
        ImageIcon icon = new ImageIcon(img);
        ImagePanel panel = new ImagePanel(img);
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        setResizable(false);
        setUndecorated(true);
        setContentPane(panel);
        setBounds((int)(screenSize.getWidth() - icon.getIconWidth())/2,
                (int) (screenSize.getHeight() - icon.getIconHeight())/2, icon.getIconWidth() + 2, icon.getIconHeight() + 2 );
        addKeyListener(new KeyAdapter()
        {
            public void keyPressed(KeyEvent e)
            {
                dispose();
            }
        });
        addMouseListener(new MouseAdapter()
        {
            public void mousePressed(MouseEvent e)
            {
                dispose();
            }
        });
    }
}
