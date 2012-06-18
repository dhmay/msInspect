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

package org.fhcrc.cpl.toolbox.gui.chart;

import javax.swing.*;
import java.awt.*;

/**
 * Render a Venn diagram in a panel
 */
public class PanelWithVenn extends JPanel
{
    int set1Size = 0;
    int set2Size = 0;
    int overlap = 0;

    public PanelWithVenn()
    {
        super();
    }

    public PanelWithVenn(int set1Size, int set2Size, int overlap)
    {
        init(set1Size, set2Size, overlap);
    }

    protected void init(int set1Size, int set2Size, int overlap)
    {
        this.set1Size = set1Size;
        this.set2Size = set2Size;
        this.overlap = overlap;
    }

    /**
     * All the shape calculations are done directly in paint().
     * None of this is exact.  Eyeballing it, it looks pretty good.
     * There will be plenty of corner cases that look bad, though.
     * @param graphics
     */
    public void paint(Graphics graphics)
    {
        super.paint(graphics);

        graphics.setColor(Color.WHITE);
        graphics.fillRect(0,0,getWidth(), getHeight());
        graphics.setColor(Color.BLACK);

        double totalSize = set1Size + set2Size - overlap;


        int set1Diam = (int) (Math.sqrt(set1Size));
        int set2Diam = (int) (Math.sqrt(set2Size));

        double size1TotalRatio = (double) set1Size / (double) (set1Size + set2Size - overlap);
        double diam1TotalRatio = set1Diam / Math.sqrt(totalSize);
        double diam2TotalRatio = set2Diam / Math.sqrt(totalSize);
        double overlapTotalRatio = overlap / totalSize;

        double widthHeightRatio =
                (set1Size + set2Size - overlap) / Math.max(set1Size,set2Size);

        int totalWidth = (int) Math.round(0.9 * this.getWidth());

        totalWidth = (int) Math.round(Math.min(totalWidth, widthHeightRatio * ((int) Math.round(0.9 * this.getHeight()))));

        int set1Diameter = (int) Math.round(diam1TotalRatio * totalWidth);
        int set2Diameter = (int) Math.round(diam2TotalRatio * totalWidth);
        int verticalCenter = getHeight() / 2;
        int horizCenter = (int) (((getWidth() * .5) + (getWidth() * size1TotalRatio))  / 2);

        int set1Top = verticalCenter - (set1Diameter / 2);
        int set2Top = verticalCenter - (set2Diameter / 2);

        int set1Left = horizCenter + (int)Math.round(overlapTotalRatio * totalWidth) -
                       set1Diameter;
        int set2Left = horizCenter - (int)Math.round(overlapTotalRatio * totalWidth);


        graphics.drawOval(set1Left, set1Top, set1Diameter, set1Diameter);
        graphics.drawOval(set2Left, set2Top, set2Diameter, set2Diameter);

        graphics.setFont(new Font("dialog", Font.BOLD, (int) (.07 * totalWidth)));

        String set1OnlyText = ("" + (set1Size - overlap));
        graphics.drawChars(set1OnlyText.toCharArray(),
                0, set1OnlyText.length(), set1Left + 5, verticalCenter-2);

        String overlapText = ("" + (overlap));
        graphics.drawChars(overlapText.toCharArray(),
                0, overlapText.length(), horizCenter - 5, verticalCenter-2);

        String set2OnlyText = ("" + (set2Size - overlap));
        graphics.drawChars(set2OnlyText.toCharArray(),
                0, set2OnlyText.length(), set1Left + set1Diameter + 5, verticalCenter-2);
    }


}
