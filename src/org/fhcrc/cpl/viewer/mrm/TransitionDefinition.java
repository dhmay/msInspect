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

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Dec 12, 2007
 * Time: 11:54:29 AM
 * To change this template use File | Settings | File Templates.
 */
public class TransitionDefinition {

    public String getPeptide() {
        return peptide;
    }

    public void setPeptide(String peptide) {
        this.peptide = peptide;
    }

    protected String peptide;

    public float getPeptideMZ() {
        return peptideMZ;
    }

    public void setPeptideMZ(float peptideMZ) {
        this.peptideMZ = peptideMZ;
    }

    protected float peptideMZ;

    public float getProductMZ() {
        return productMZ;
    }

    public void setProductMZ(float productMZ) {
        this.productMZ = productMZ;
    }

    protected float productMZ;

    public char getLowOrHigh() {
        return lowOrHigh;
    }

    public void setLowOrHigh(char lowOrHigh) {
        this.lowOrHigh = lowOrHigh;
    }

    protected char lowOrHigh;

    public boolean isLow() {
        return (lowOrHigh == 'l') || (lowOrHigh == 'L');
    }

    public boolean isHigh() {
        return (lowOrHigh == 'h') || (lowOrHigh == 'H');
    }

    public int getAQUAcode() {
        return AQUAcode;
    }

    public void setAQUAcode(int AQUAcode) {
        this.AQUAcode = AQUAcode;
    }

    protected int AQUAcode = -1;

    public MRMDaughter getAssociatedProduct() {
        return associatedProduct;
    }

    public void setAssociatedProduct(MRMDaughter associatedProduct) {
        this.associatedProduct = associatedProduct;
    }

    protected MRMDaughter associatedProduct;

    public TransitionDefinition(String peptide, float peptideMZ, float productMZ) {
        this.peptide = peptide;
        this.peptideMZ = peptideMZ;
        this.productMZ = productMZ;
    }


    public TransitionDefinition(String peptide, float peptideMZ, float productMZ, char lowOrHigh, int AQUAcode) {
        this.peptide = peptide;
        this.peptideMZ = peptideMZ;
        this.productMZ = productMZ;
        this.lowOrHigh = lowOrHigh;
        this.AQUAcode = AQUAcode;
    }

    public static class TransitionDefinitionPeptideMZAscComparator implements Comparator<TransitionDefinition>
    {
        public int compare(TransitionDefinition o1, TransitionDefinition o2)
        {
             return o1.peptideMZ == o2.peptideMZ ? 0 : o1.peptideMZ < o2.peptideMZ ? -1 : 1;
        }
    }
}
