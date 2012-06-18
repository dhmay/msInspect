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
package org.fhcrc.cpl.viewer.quant.commandline;

import org.fhcrc.cpl.viewer.commandline.modules.DeconvoluteCommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.BooleanArgumentDefinition;
import org.apache.log4j.Logger;


/**
 * Quantitation.  Same as DeconvoluteCommandLineModule, but with the "quant"
 * argument set to true, and allowing you to tweak the "deconvolute" argument 
 */
public class QuantCommandLineModule extends DeconvoluteCommandLineModule
        implements CommandLineModule
{
    protected static Logger _log = Logger.getLogger(QuantCommandLineModule.class);

    protected void init()
    {
        mCommandName = "quant";
        mShortDescription = "Quantitate";
        mHelpMessage = "Quantitate";

        addArgumentDefinitions(createCommonArgDefs());
        addArgumentDefinition(new BooleanArgumentDefinition("deconvolute",false,
                "Deconvolute",false));
    }

    public void assignArgumentValues()
            throws ArgumentValidationException
    {
        quant = true;        
        super.assignArgumentValues();

        assertArgumentPresent("labeledresidue");
        assertArgumentPresent("lighttagweight");
        assertArgumentPresent("heavytagweight");
        deconvolute = getBooleanArgumentValue("deconvolute");
    }
}
