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
package org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;


/**
 * Contains column name and datatype information about each column.
 * Each column datatype must be a class with a constructor that accepts
 * one String argument.
 */
public class TimeExtraInfoDef extends FeatureExtraInformationDef
{
    static Logger _log = Logger.getLogger(TimeExtraInfoDef.class);


    public TimeExtraInfoDef()
    {
        super(
                    "TIME",
                    new String[]{"starttime","endtime"},
                    new Class[]{Double.class,Double.class}
            );

    }

    public static double getStartTime(Feature feature)
    {
        return feature.getDoubleProperty("starttime",-1);
    }

    public static void setStartTime(Feature feature, double startTime)
    {
        feature.setProperty("starttime", startTime);
    }

    public static double getEndTime(Feature feature)
    {
        return feature.getDoubleProperty("endtime",-1);
    }

    public static void setEndTime(Feature feature, double endTime)
    {
        feature.setProperty("endtime", endTime);
    }

    protected static TimeExtraInfoDef singletonInstance;

    public static TimeExtraInfoDef getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new TimeExtraInfoDef();
        return singletonInstance;
    }
}
