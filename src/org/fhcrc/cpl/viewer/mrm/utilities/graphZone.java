package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.mrm.MRMTransition;
import org.jfree.data.Range;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Sep 19, 2007
 * Time: 2:38:47 PM
 * To change this template use File | Settings | File Templates.
 */
 public interface graphZone {
        Range getRange(MRMTransition mrmt);
 }