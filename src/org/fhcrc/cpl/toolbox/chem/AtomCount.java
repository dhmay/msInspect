package org.fhcrc.cpl.toolbox.chem;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 26, 2010
 * Time: 3:30:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class AtomCount {
    Double mass;
    Element el;
    int count;

    public Element getEl() {
        return el;
    }

    public void setEl(Element el) {
        this.el = el;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public Double getMass() {
        return mass;
    }

    public void setMass(Double mass) {
        this.mass = mass;
    }

    public AtomCount(Double mass, int count) {
        this.mass = mass;
        this.count = count;
    }

    public AtomCount(Double mass, int count, Element el) {
        this.mass = mass;
        this.count = count;
        this.el = el;
    }

    Double contrib() {
        return mass * count;
    }
}

