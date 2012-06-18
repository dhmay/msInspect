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

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Dec 12, 2007
 * Time: 12:48:08 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TransitionDefinitionParser {

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return this.version;
    }

    protected String version;

    public String getTransitionDefFile() {
        return transitionDefFile;
    }

    public void setTransitionDefFile(String transitionDefFile) {
        this.transitionDefFile = transitionDefFile;
    }

    protected String transitionDefFile;

    public String getMzXMLFile() {
        return mzXMLFile;
    }

    public void setMzXMLFile(String mzXMLFile) {
        this.mzXMLFile = mzXMLFile;
    }

    protected String mzXMLFile;

    public String getComment() {
        return comment;
    }

    public void setComment(String comment) {
        this.comment = comment;
    }

    protected String comment;

    public ArrayList getTransitionDefs() {
        return transitionDefs;
    }

    public void setTransitionDefs(ArrayList transitionDefs) {
        this.transitionDefs = transitionDefs;
    }

    protected ArrayList<TransitionDefinition> transitionDefs;

    public abstract void parse(String file) throws Exception;

    public void parse() throws Exception {
        parse(getTransitionDefFile());
    }

}
