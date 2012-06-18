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

import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Dec 12, 2007
 * Time: 12:26:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransitionDefinitionHeader {
    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
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

    public ArrayList<TransitionDefinition> getTransitionDefs() {
        return transitionDefs;
    }

    public void setTransitionDefs(ArrayList<TransitionDefinition> transitionDefs) {
        this.transitionDefs = transitionDefs;
    }

    protected ArrayList<TransitionDefinition> transitionDefs;

    public TransitionDefinitionParser getParser() {
        return parser;
    }

    public void setParser(TransitionDefinitionParser parser) {
        this.parser = parser;
    }

    protected TransitionDefinitionParser parser;

    public void doParse() throws Exception {
       if(this.getParser() == null) {
           ApplicationContext.errorMessage(TextProvider.getText("NO_TRANSITION_DEFINITION_PARSER_SET"), new Exception("Must set a transition definition parser before trying to read a transition definition file."));
       } else {
           this.getParser().setTransitionDefFile(this.getTransitionDefFile());
           this.getParser().parse();
           this.setVersion(parser.getVersion());
           this.setMzXMLFile(parser.getMzXMLFile());
           this.setComment(parser.getComment());
           this.setTransitionDefs(parser.getTransitionDefs());
       }
    }

    public TransitionDefinitionHeader(String transitionDefFile, TransitionDefinitionParser parser) {
        this.transitionDefFile = transitionDefFile;
        this.parser = parser;
    }

    public TransitionDefinitionHeader(String transitionDefFile, TransitionDefinitionParser parser, boolean shouldParse) {
        this.transitionDefFile = transitionDefFile;
        this.parser = parser;
        if(shouldParse)
           try {doParse();}
           catch (Exception e){};
    }

    public TransitionDefinitionHeader() {
    }

    public HashMap<MRMDaughter, TransitionDefinition> getDToTD() {
        return dToTD;
    }

    public void setDToTD(HashMap<MRMDaughter, TransitionDefinition> dToTD) {
        this.dToTD = dToTD;
    }

    protected HashMap<MRMDaughter,TransitionDefinition> dToTD;

    public void linkUpToTransitionList(MRMTransition mrta[],float ptol, float dtol) {
        setDToTD(new HashMap<MRMDaughter,TransitionDefinition>());
        for(MRMTransition curTrans: mrta) {
            for(MRMDaughter mrmd: curTrans.getDaughters().values()){
                for(TransitionDefinition curDef: getTransitionDefs()){
                    if(curDef.getPeptideMZ()-ptol <= curTrans.getPrecursorMz() &&
                       curDef.getPeptideMZ()+ptol >= curTrans.getPrecursorMz() &&
                       curDef.getProductMZ()-dtol <= mrmd.getMeanMz()          &&
                       curDef.getProductMZ()+dtol >= mrmd.getMeanMz()
                    ) {
                        curDef.setAssociatedProduct(mrmd);
                        getDToTD().put(mrmd,curDef);
                    }
                }
            }
        }
    }

    public class AQUApair {

        public TransitionDefinition getLightMember() {
            return lightMember;
        }

        public void setLightMember(TransitionDefinition lightMember) {
            this.lightMember = lightMember;
        }

        public TransitionDefinition getHeavyMember() {
            return heavyMember;
        }

        public void setHeavyMember(TransitionDefinition heavyMember) {
            this.heavyMember = heavyMember;
        }

        TransitionDefinition lightMember;
        TransitionDefinition heavyMember;


        public AQUApair(TransitionDefinition lowMember, TransitionDefinition highMember) {
            this.lightMember = lowMember;
            this.heavyMember = highMember;
        }
    }


    public HashMap<Integer,AQUApair> getAQUApairs() {
        return AQUApairs;
    }

    public void setAQUApairs(HashMap<Integer,AQUApair> AQUApairs) {
        this.AQUApairs = AQUApairs;
    }

    protected HashMap<Integer,AQUApair> AQUApairs = null;

    public void determinePairs() {
        if(this.getTransitionDefs() == null || this.getTransitionDefs().isEmpty()) return;
        HashMap<Integer,Vector<TransitionDefinition>> pairs =
           new HashMap<Integer,Vector<TransitionDefinition>>();

        //scan potential aquapairs
        for(TransitionDefinition td: getTransitionDefs()) {
            if(td.getAQUAcode() <= 0) continue;
            Integer aqc = new Integer(td.getAQUAcode());
            Vector<TransitionDefinition> tds = pairs.get(aqc);
            if(tds == null) {
                tds = new Vector<TransitionDefinition>();
                pairs.put(aqc,tds);
            }
            tds.add(td);
        }
        //Check for sanity, load into AQUApairs vector if sane
        setAQUApairs(new HashMap<Integer,AQUApair>());
        for(Integer code: pairs.keySet()) {
            Vector<TransitionDefinition> tdefs = pairs.get(code);
            if(tdefs.size() != 2) {
                ApplicationContext.errorMessage(TextProvider.getText("BAD_AQUA_PAIR"," There must be exactly two members for each code, code:"+code+" has "+tdefs.size()+" members"),null);
                continue;
            }
            TransitionDefinition td1 = tdefs.get(0);
            TransitionDefinition td2 = tdefs.get(1);
            if(!td1.getPeptide().equalsIgnoreCase(td2.getPeptide())){
                ApplicationContext.errorMessage(TextProvider.getText("BAD_AQUA_PAIR"," Both members of an AQUA pair must relate to the same peptide, code:"+code+" has "+td1.getPeptide()+" and "+td2.getPeptide()),null);
                continue;
            }
            TransitionDefinition high=null, low=null;
            if(td1.isHigh()) high = td1;
            if(td1.isLow()) low = td1;
            if(td2.isHigh()) high = td2;
            if(td2.isLow()) low = td2;
            if(high == null || low == null || high == low ){
                ApplicationContext.errorMessage(TextProvider.getText("BAD_AQUA_PAIR"," There must be one heavy and one light member of each pair, code:"+code+" does not obey this rule"),null);
                continue;
            }
            getAQUApairs().put(code,new AQUApair(low,high));
        }
    }

    public MRMDaughter getMatchingDaughter(MRMDaughter mrmd) {
        MRMDaughter retVal = null;
        if(getAQUApairs() == null) return retVal;
        if(mrmd == null) return retVal;
        TransitionDefinition curTD = getDToTD().get(mrmd);
        if(curTD == null)return retVal;
        if(getAQUApairs() == null || getAQUApairs().size() == 0) return retVal;
        AQUApair curAQUA = getAQUApairs().get(curTD.getAQUAcode());
        if(curAQUA == null) return retVal;
        if(curTD.isHigh()){
            retVal = curAQUA.getLightMember().getAssociatedProduct();
        } else {
            if(curTD.isLow()) {
                retVal = curAQUA.getHeavyMember().getAssociatedProduct();
            }
        }
        return retVal; 
    }
}
