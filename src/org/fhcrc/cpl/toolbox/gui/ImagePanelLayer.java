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

import java.awt.*;

/**
 * Interface to be implemented by anything that draws to the image panel.
 * Implementing classes should play nice w.r.t. the opacity
 */
public interface ImagePanelLayer
{
    public String getName();

    public void draw(Graphics graphics,
                     int imageWidth, int imageHeight,
                     int minVisibleScan, int maxVisibleScan,
                     double minVisibleMz, double maxVisibleMz);

    public void setTranslucence(int translucence);

    public int getTranslucence();
}
