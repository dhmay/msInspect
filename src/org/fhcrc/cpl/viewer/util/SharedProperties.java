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
package org.fhcrc.cpl.viewer.util;

/**
 * User: migra
 * Date: Oct 4, 2004
 * Time: 2:41:12 PM
 */
public interface SharedProperties
{
    public static final String MS_RUN = "MSRun";
    public static final String MS_SCAN = "MSScan";
    public static final String FEATURE_SETS = "featureSets"; //Raw features
    public static final String FEATURE_RANGES = "featureRanges"; //Filtered & merged features
    public static final String HIGHLIGHT_FEATURES = "highlightFeatures"; //Features to circle
    public static final String HIGHLIGHT_FEATURES2 = "highlightFeatures2"; //More features to circle, taking precedence over previous
    public static final String HIGHLIGHT_FEATURES3 = "highlightFeatures3"; //More features to circle, taking precedence over previous        
    public static final String SELECTED_POINT = "selectedPoint";
    public static final String ZOOM_REGION = "zoomRegion";
	public static final String AUTO_ZOOM = "autoZoom";
	public static final String LOCK_Y_AXIS = "lockYAxis";
	public static final String CHART_RANGE = "chartRange";
	public static final String SELECTED = "selectedObject";
    public static final String FEATURE_VIEWER = "featureViewer";
}
