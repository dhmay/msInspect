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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Dec 12, 2007
 * Time: 1:14:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class TSVTransitionDefinitionParser extends TransitionDefinitionParser{
    public void parse(String f) throws Exception {
       String line;
       try {
           if(!(new File(f)).exists()) throw new Exception("Transition file '"+f+"' does not exist.");
           BufferedReader bf = new BufferedReader(new FileReader(f));
           setTransitionDefFile(f);
           line = bf.readLine();
           String tokens[] = line.split("\\t");
           if(tokens.length < 2 || !tokens[0].equalsIgnoreCase("Transition definition file"))
              throw new Exception(f+" is not a transition definition file");
           setVersion(tokens[1]);
           setMzXMLFile(bf.readLine());
           setComment(bf.readLine());
           setTransitionDefs(new ArrayList<TransitionDefinition>());
           while((line = bf.readLine())!= null) {
               tokens = line.split("\\t");
               if(tokens.length < 3) {
                   ApplicationContext.errorMessage(TextProvider.getText("CANNOT_PARSE_TRANSITION_DEFINITION"),new Exception("Definition lines must have at least 3 tab-delimited tokens"));
                   continue;
               }
               String peptide = tokens[0];
               float PepMZ, ProductMZ;
               try {
                  PepMZ = Float.parseFloat(tokens[1]);
                  ProductMZ = Float.parseFloat(tokens[2]);
               } catch (Exception e) {
                  ApplicationContext.errorMessage(TextProvider.getText("CANNOT_PARSE_TRANSITION_DEFINITION"),new Exception("Precursor and Product MZ fields must contain valid floating point numbers"));
                  continue;
               }
               TransitionDefinition curTD = new TransitionDefinition(tokens[0],PepMZ,ProductMZ);
               getTransitionDefs().add(curTD);
               if(tokens.length >= 4) {
                  if("LH".contains(tokens[3].toUpperCase())) {
                      curTD.setLowOrHigh(tokens[3].charAt(0));
                  } else {
                      if(tokens[3].length()>0)
                         ApplicationContext.errorMessage(TextProvider.getText("CANNOT_PARSE_TRANSITION_DEFINITION"),new Exception("LowOrHigh field must contain 'L' or 'H', not '"+tokens[3]+"'"));
                  }
               }
               if(tokens.length >= 5) {
                 try {
                    int AQUAcode = Integer.parseInt(tokens[4]);
                    curTD.setAQUAcode(AQUAcode);
                 } catch (Exception e) {
                     ApplicationContext.errorMessage(TextProvider.getText("CANNOT_PARSE_TRANSITION_DEFINITION"),new Exception("AQUAPair field must contain an integer, not '"+tokens[4]+"'"));
                 }
               }
           }
       } catch (Exception e) {
           ApplicationContext.errorMessage(TextProvider.getText("CANNOT_PARSE_TRANSITION_DEFINITION_FILE"),e);
           throw e;
       }
    }
}
