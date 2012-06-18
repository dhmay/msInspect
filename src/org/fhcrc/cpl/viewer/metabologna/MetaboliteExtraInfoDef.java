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
package org.fhcrc.cpl.viewer.metabologna;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.Feature;
import org.fhcrc.cpl.toolbox.proteomics.feature.extraInfo.FeatureExtraInformationDef;
import org.fhcrc.cpl.toolbox.chem.ChemicalFormula;

import java.util.*;
import java.lang.reflect.InvocationTargetException;

/**
 * FeatureExtraInformationDef related to metabolite analysis
 */
public class MetaboliteExtraInfoDef extends FeatureExtraInformationDef
{            
    static Logger _log = Logger.getLogger(MetaboliteExtraInfoDef.class);

    /**
     * initialize the name of the extra information, the column list, and the datatypes
     * of the columns
     */
    public MetaboliteExtraInfoDef()
    {
        super(
                    "METABOLITE",
                    new String[]{
                            "formulas",
                    },
                    new Class[]{
                            List.class,
                    },
                    new String[]{
                    }

            );
    }

    protected static MetaboliteExtraInfoDef singletonInstance = null;

    public static MetaboliteExtraInfoDef getSingletonInstance()
    {
        if (singletonInstance == null)
            singletonInstance = new MetaboliteExtraInfoDef();
        return singletonInstance;
    }
    public Object convertStringValue(String columnName, String value)
            throws InstantiationException, IllegalAccessException,
            InvocationTargetException, NoSuchMethodException
    {
        if ("formulas".equalsIgnoreCase(columnName))
        {
            List<String> formulaStrings = parseStringListString(value);
            if (formulaStrings == null || formulaStrings.isEmpty())
                return null;
            List<ChemicalFormula> formulas = new ArrayList<ChemicalFormula>();
            for (String formulaString : formulaStrings)
                formulas.add(new ChemicalFormula(formulaString));
            return formulas;
        }
        else
            return super.convertStringValue(columnName, value);
    }

    /**
     * Special handling for peptide, protein columns
     *
     * TODO: special handling for modifiedaminoacids
     * @param columnName
     * @param value
     * @return
     */
    public String convertToString(String columnName, Object value)
    {
        if ("formulas".equalsIgnoreCase(columnName))
        {
            List<ChemicalFormula> formulas = (List<ChemicalFormula>) value;
            if (formulas == null || formulas.isEmpty())
                return "";
            List<String> resultList = new ArrayList<String>();
            for (ChemicalFormula formula : formulas)
                resultList.add(formula.toString());

            return convertStringListToString(resultList);
        }
        else return super.convertToString(columnName, value);
    }


    public static List<ChemicalFormula> getFormulaList(Feature feature)
    {
        return (List<ChemicalFormula>) feature.getProperty("formulas");
    }

    public static void setFormulaList(Feature feature,
                                      List<ChemicalFormula> formulaList)
    {
        feature.setProperty("formulas", formulaList);
    }

    public static void removeAllFormulas(Feature feature)
    {
        feature.setProperty("formulas", null);
    }

    /**
     * Adds a peptide with an associated Protein
     * @param feature
     * @param formula
     */
    public static void addFormula(Feature feature, ChemicalFormula formula)
    {
        List<ChemicalFormula> formulaList = getFormulaList(feature);
        if (formulaList == null)
        {
            formulaList = new ArrayList<ChemicalFormula>();
            setFormulaList(feature, formulaList);
        }
        formulaList.add(formula);
    }

}
