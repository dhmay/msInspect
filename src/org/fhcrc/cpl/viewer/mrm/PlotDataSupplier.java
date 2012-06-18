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
package org.fhcrc.cpl.viewer.mrm;

import org.jfree.data.xy.XYSeries;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 19, 2007
 * Time: 12:05:16 PM
 * To change this template use File | Settings | File Templates.
 */
public interface PlotDataSupplier {
    public XYSeries getGraphData();
    public void setGraphData(XYSeries xys);
    public Paint getGraphColor();
    public void setGraphColor(Paint p);
    public boolean isContinuous();
}
