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

package org.fhcrc.cpl.viewer.commandline;

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleDiscoverer;
import org.apache.log4j.Logger;

import java.util.*;
import java.net.JarURLConnection;
import java.net.URL;
import java.net.URLDecoder;
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Modifier;

/**
 * Utility class for finding all command-line modules.
 * See comments on method discoverAllCommandLineModules() for important details.
 */
public class ViewerCommandLineModuleDiscoverer  extends CommandLineModuleDiscoverer
{
    protected static Logger _log = Logger.getLogger(ViewerCommandLineModuleDiscoverer.class);


    protected static Map<String,CommandLineModule> allCommandLineModuleMap = null;

    protected static ViewerCommandLineModuleDiscoverer singletonInstance = null;


    public static ViewerCommandLineModuleDiscoverer getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new ViewerCommandLineModuleDiscoverer();
        return singletonInstance;
    }

    public ViewerCommandLineModuleDiscoverer()
    {

        //TODO: stick all this stuff into a resource file
        //packages in which all modules must reside
        modulePackageNames =
                new String[]
                        {
                                "org.fhcrc.cpl.viewer.commandline.modules",
                                "org.fhcrc.cpl.viewer.quant.commandline",
                                "org.fhcrc.cpl.viewer.align.commandline",
                                "org.fhcrc.cpl.viewer.ms2.commandline",
                                "org.fhcrc.cpl.viewer.amt.commandline",
                                "org.fhcrc.cpl.viewer.mrm.commandline",
                                "org.fhcrc.cpl.viewer.qa.commandline",
                                "org.fhcrc.cpl.viewer.metabologna.commandline",
                        };

        //an identifier string unique to this package, with only alphanumeric characters
        modulePackageIdentifiers =
                new String[]
                        {
                                "general",
                                "quant",
                                "align",
                                "ms2",
                                "amt",
                                "mrm",
                                "qa",
                                "metabologna",
                        };

        //Short (1-2 words) user-readable descriptions of each commandline module package -- what's it for?
        modulePackageDescriptionsShort =
                new String[]
                        {
                                "General",
                                "Quantitation",
                                "Alignment",
                                "MS/MS",
                                "AMT",
                                "MRM",
                                "Quality Assurance",
                                "Metabolite Utilities"
                        };

        //Longer (1-2 sentences) user-readable descriptions of each commandline module package -- what's it for?
        modulePackageDescriptionsLong =
                new String[]
                        {
                                "General msInspect tools",
                                "Tools related to labeled and unlabeled quantitation",
                                "Tools for aligning multiple runs, including Peptide Arrays",
                                "Tools related to tandem mass spectrometry and the file formats used for tandem MS",
                                "Accurate Mass and Time analysis tools",
                                "The MRMer tools for Multiple Reaction Monitoring",
                                "Quality Assurance tools",
                                "Metabolite identification and quantitation",
                        };
    }


}
