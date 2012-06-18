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
package org.fhcrc.cpl.viewer.amt;

import java.util.*;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.MS2Modification;
import org.fhcrc.cpl.toolbox.statistics.RegressionUtilities;

/**
 * Encapsulates everything we need to retrieve from / store into an AMT database related to a run.
 */
public class AmtRunEntry
{
    static Logger _log = Logger.getLogger(AmtPeptideEntry.class);



    //required stuff

    protected Date mTimeAdded = null;

    // coefficients for the function mapping time to Hydrophobicity in this run
    protected double[] mTimeHydMapCoefficients;

    //known modifications
    MS2Modification[] mModifications = null;


    //optional stuff

    protected Date mTimeAnalyzed = null;
    protected String mPepXmlFilename = null;
    protected String mMzXmlFilename = null;
    protected String mLSID = null;
    //minimum peptideprophet score required for inclusion in this run
    protected double mMinPeptideProphet = 0.0;

    /**
     * Created a new run entry dated right now
     */
    public AmtRunEntry(double[] timeHydMapCoefficients,
                       MS2Modification[] modifications)
    {
        this(timeHydMapCoefficients, modifications, new Date());
    }

    public AmtRunEntry(double[] timeHydMapCoefficients,
                       MS2Modification[] modifications,
                       Date timeAdded)
    {
        mTimeAdded = timeAdded;
        mTimeHydMapCoefficients = timeHydMapCoefficients;
        mModifications = modifications;
    }

    public String toString()
    {
        StringBuffer resultBuf =
            new StringBuffer("AmtRunEntry: TimeHydMapFunction=");
        for (int i = mTimeHydMapCoefficients.length-1; i>=0; i--)
        {
            double coeff = mTimeHydMapCoefficients[i];
            if (i < mTimeHydMapCoefficients.length-1 && coeff > 0)
                resultBuf.append("+");
            resultBuf.append(coeff + "t^" + i);
        }
        return resultBuf.toString();
    }

    /**
     * Recover the time, within this run, that something with a given hydrophobicity
     * value must have occurred.  Work the regression line backward.
     *
     * Can ONLY do this for linear mapping
     * @param hydrophobicity
     * @return
     */
    public double recoverTimeForHydrophobicity(double hydrophobicity)
    {
        if (mTimeHydMapCoefficients.length == 2)
            return RegressionUtilities.predictXFromY(mTimeHydMapCoefficients[1], mTimeHydMapCoefficients[0], hydrophobicity);
        throw new IllegalArgumentException("Can't recover time from hydrophobicity unless mapping is linear");
    }


    /**
     * convert time to hydrophobicity for this run
     * @param time
     * @return
     */
    public double convertTimeToHydrophobicity(double time)
    {
        return RegressionUtilities.mapValueUsingCoefficients(mTimeHydMapCoefficients, time);
    }














    //getters and setters

    public Date getTimeAdded()
    {
        return mTimeAdded;
    }

    public void setTimeAdded(Date mTimeAdded)
    {
        this.mTimeAdded = mTimeAdded;
    }

    public Date getTimeAnalyzed()
    {
        return mTimeAnalyzed;
    }

    public void setTimeAnalyzed(Date mTimeAnalyzed)
    {
        this.mTimeAnalyzed = mTimeAnalyzed;
    }

    public String getPepXmlFilename()
    {
        return mPepXmlFilename;
    }

    public void setPepXmlFilename(String mPepXmlFilename)
    {
        this.mPepXmlFilename = mPepXmlFilename;
    }

    public String getMzXmlFilename()
    {
        return mMzXmlFilename;
    }

    public void setMzXmlFilename(String mMzXmlFilename)
    {
        this.mMzXmlFilename = mMzXmlFilename;
    }

    public double[] getTimeHydMapCoefficients()
    {
        return mTimeHydMapCoefficients;
    }

    public void setTimeHydMapCoefficients(double[] timeHydMapCoefficients)
    {
        mTimeHydMapCoefficients = timeHydMapCoefficients;
    }

    public String getLSID()
    {
        return mLSID;
    }

    public void setLSID(String lsid)
    {
        this.mLSID = lsid;
    }

    public double getMinPeptideProphet()
    {
        return mMinPeptideProphet;
    }

    public void setMinPeptideProphet(double minPeptideProphet)
    {
        this.mMinPeptideProphet = minPeptideProphet;
    }

    public MS2Modification[] getModifications()
    {
        return mModifications;
    }

    //TODO: these two methods should REALLY be cached somewhere
    public MS2Modification[] getStaticModifications()
    {
        List<MS2Modification> resultList = new ArrayList<MS2Modification>();
        if (mModifications != null)
            for (MS2Modification mod : mModifications)
            {
                if (!mod.getVariable())
                    resultList.add(mod);
            }
        if (resultList.size() == 0)
            return null;
        return resultList.toArray(new MS2Modification[resultList.size()]);
    }

    public MS2Modification[] getVariableModifications()
    {
        List<MS2Modification> resultList = new ArrayList<MS2Modification>();
        for (MS2Modification mod : mModifications)
        {
            if (mod.getVariable())
                resultList.add(mod);
        }
        if (resultList.size() == 0)
            return null;
        return resultList.toArray(new MS2Modification[resultList.size()]);
    }


    public void setModifications(MS2Modification[] modifications)
    {
        mModifications = modifications;
    }

    public MS2Modification getStaticMod(String residue)
    {
        for (MS2Modification mod : getStaticModifications())
        {
            if (mod.getAminoAcid().equals(residue))
                return mod;
        }
        return null;
    }

    /**
     * when a run is folded into a new database, that database's mods take precedence.
     * They're equivalent, but this makes the pointers point in the right places.
     * @param amtDB
     */
    public Map<MS2Modification, MS2Modification> overrideDuplicateModifications(AmtDatabase amtDB)
    {

        Map<MS2Modification, MS2Modification> oldNewModMap =
                    new HashMap<MS2Modification, MS2Modification>();
        if (mModifications == null)
            return oldNewModMap;        
        for (int i=0; i<mModifications.length; i++)
        {
            MS2Modification dbMod =
                amtDB.findExistingEquivalentModification(mModifications[i]);
            if (dbMod != null)
            {
                _log.debug("  overriding mod: " + dbMod.getAminoAcid() + ", " + dbMod.getMassDiff());
                oldNewModMap.put(mModifications[i],dbMod);

                mModifications[i] = dbMod;
            }
        }
        return oldNewModMap;
    }



}
