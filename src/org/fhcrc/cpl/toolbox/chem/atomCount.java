package org.fhcrc.cpl.toolbox.chem;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 26, 2010
 * Time: 3:30:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class atomCount {
   Double mass;
   Element el;
   int count;

   public Double getMass() {
      return mass;
   }

   public void setMass(Double mass) {
      this.mass = mass;
   }

   public atomCount(Double mass, int count) {
       this.mass=mass;
       this.count=count;
    }

   public atomCount(Double mass, int count, Element el) {
       this.mass=mass;
       this.count=count;
       this.el=el;
    }

   Double contrib() {
      return mass*count;
   }
}

