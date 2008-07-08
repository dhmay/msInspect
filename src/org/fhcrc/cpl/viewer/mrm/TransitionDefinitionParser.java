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
