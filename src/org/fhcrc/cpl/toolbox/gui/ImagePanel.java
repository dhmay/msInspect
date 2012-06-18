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
import java.awt.*;


/**
 * User: mbellew
 * Date: May 21, 2004
 * Time: 5:28:36 PM
 */
public class ImagePanel extends JPanel//implements Scrollable
    {
    protected Image img = null;
    protected int _x;
    protected int _y;


    public ImagePanel()
    {

    }

    public ImagePanel(Image image)
        {
        img = image;
        }


    public void paintComponent(Graphics g)
        {
        super.paintComponent(g);
        _x = Math.max(0, getSize().width - img.getWidth(this)) / 2;
        _y = Math.max(0, getSize().height - img.getHeight(this)) / 2;
        g.drawImage(img, _x, _y, this);
        }
    }
