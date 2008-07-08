package org.fhcrc.cpl.viewer.mrm;

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
}
