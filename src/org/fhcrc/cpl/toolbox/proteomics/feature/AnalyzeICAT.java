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
package org.fhcrc.cpl.toolbox.proteomics.feature;

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.toolbox.datastructure.Tree2D;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Collections;
import java.io.File;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: mbellew
 * Date: May 23, 2005
 * Time: 11:30:21 AM
 */
public class AnalyzeICAT
	{
	public static int DEFAULT_MAX_LABEL_COUNT = 3;

    protected static ArrayList prePopulatedTags = null;

    public static class IsotopicLabel
		{
        //an identifier for this label.  Human-readable
        protected String name = "";
        protected float light;    // weight of base tag
		protected float heavy;    // delta weight (e.g. 9.03 for icat)
		protected char residue = ' '; // resiude labeled; ' ' means label affects non-residue (e.g. C or N terminus)
		protected int maxLabelCount;  // maximum number of labels to consider per peptide

		public IsotopicLabel(String name, float light, float heavy, char residue, int maxLabelCount)
		    {
            this(light, heavy, residue, maxLabelCount);
		    setName(name);
            }

        public IsotopicLabel(float light, float heavy, char residue, int maxLabelCount)
		    {
		    setLight(light);
		    setHeavy(heavy);
		    setResidue(residue);
		    setMaxLabelCount(maxLabelCount);
		    }

		/**
		 * Parse an isotopic label string. An example string representation is "100.0+4.05#3@C" to
		 * be parsed as floating point light tag weight, plus symbol, delta weight for the heavy
		 * label, number symbol, maximum number of labels to consider, at sign, residue to be
		 * labeled.
		 *
		 * If the maximum number of labels to consider is left off, it defaults to 3. If supplied,
		 * this must be positive.
		 * If the residue to be labeled is left off, it defaults to a blank to indicate that the
		 * labeling scheme a part of the peptide other than a specific residue (e.g. the C or N terminus).
		 */
		public IsotopicLabel(String s)
		    {
		    float light, heavy;
		    char residue = ' ';
		    int maxLabelCount = DEFAULT_MAX_LABEL_COUNT;
		    
		    int plus = s.indexOf("+");
		    int hash = s.indexOf("#");
		    int at = s.indexOf("@");
		    
		    if (plus == -1)
			throw new IllegalArgumentException("Malformed label string; expected light tag weight + delta mass: " + s);
		    
		    if (hash != -1 && hash < plus)
			throw new IllegalArgumentException("Malformed label string; '#' must follow '+': " + s);
		    
		    if (at != -1 && (at < plus || at < hash))
			throw new IllegalArgumentException("Malformed label string; '@' must follow '+' and '#': " + s);
		    
		    light = Float.parseFloat(s.substring(0,plus));
		    
		    if ( -1 != hash )
			heavy = Float.parseFloat(s.substring(plus+1,hash));
		    else if ( -1 != at )
			heavy = Float.parseFloat(s.substring(plus+1,at));
		    else
			heavy = Float.parseFloat(s.substring(plus+1));
		    
		    if ( -1 != hash )
			if ( -1 != at )
			    maxLabelCount = Integer.parseInt(s.substring(hash+1,at));
			else
			    maxLabelCount = Integer.parseInt(s.substring(hash+1));
		    
		    if ( -1 != at )
			residue = s.charAt(at+1);
		    
		    setLight(light);
		    setHeavy(heavy);
		    setResidue(residue);
		    setMaxLabelCount(maxLabelCount);
		    }
		
		public String toString()
		    {
		    return "" + getLight() + "+" + getHeavy()
			+ (getMaxLabelCount() == DEFAULT_MAX_LABEL_COUNT ? "" : "#" + getMaxLabelCount())
			+ (getResidue() == ' ' ? "" : "@" + getResidue());
		    }

		public String getName()
		    {
		    return name;
		    }

		private void setName(String name)
		    {
		    this.name = name;
		    }

        public float getLight()
		    {
		    return light;
		    }

		private void setLight(float light) 
		    {
		    if (light < 0.f)
			throw new IllegalArgumentException("Light tag weight must be non-negative");
		    this.light = light;
		    }

		public float getHeavy()
		    {
		    return heavy;
		    }

		private void setHeavy(float heavy) 
		    {
		    if (heavy <= 0.f)
			throw new IllegalArgumentException("Heavy tag weight must be positive");
		    this.heavy = heavy;
		    }

		public char getResidue()
		    {
		    return residue;
		    }

		private void setResidue(char residue) 
		    {
		    if ( residue != ' ' && (residue < 'A' || residue > 'Y' ) )
			throw new IllegalArgumentException("Illegal residue character '" + residue + "'");
		    this.residue = residue;
		    }

		public int getMaxLabelCount()
		    {
		    return maxLabelCount;
		    }

		private void setMaxLabelCount(int maxLabelCount)
		    {
		    if ( maxLabelCount <= 0 )
			throw new IllegalArgumentException("maxLabelCount must be at least 1");
		    this.maxLabelCount = maxLabelCount;
		    }

	    }
	
//dhmay changing, consolidating default into the one arraylist
//	public static IsotopicLabel icatLabel = new IsotopicLabel( 227.1263F, 9.0297F, 'C', DEFAULT_MAX_LABEL_COUNT); // NOTE: I computed these myself (mbellew)
    public static IsotopicLabel icatLabel = (IsotopicLabel) getPrePopulatedLabels().get(0);

    //constants used to specify how mass tolerance should be calculated
    public static final int DELTA_MASS_ABSOLUTE = 0;
    public static final int DELTA_MASS_PPM = 1;

    //defaults for mass and time tolerances
    public static final float DEFAULT_DELTA_MASS = 0.2F;
	public static final float DEFAULT_DELTA_TIME = 10.0F;
    public static final int DEFAULT_DELTA_MASS_TYPE = DELTA_MASS_ABSOLUTE;

        /**
         * Utility method to calculate the absolute mass tolerance, given a mass tolerance
         * parameter that may be absolute or relative
         * @param centerMass
         * @param deltaMass
         * @param deltaMassType
         * @return
         */
    protected static float calculateAbsoluteDeltaMass(float centerMass,
                                               float deltaMass,
                                               int deltaMassType)
    {
        if (deltaMassType == DELTA_MASS_ABSOLUTE)
            return deltaMass;
        //deltamass must be in ppm
        return (deltaMass * centerMass) / 1000000;
    }

        /**
         * Uses default values for mass and time tolerance
         * @param featuresIN
         * @param label
         * @return
         */
    public static ArrayList analyze1(Feature[] featuresIN, IsotopicLabel label)
    {
        return analyze1(featuresIN, label, DEFAULT_DELTA_MASS, DELTA_MASS_ABSOLUTE,
                        DEFAULT_DELTA_TIME);
    }

        /**
         * can be run on deconvoluted or non-deconvoluted feature set
         * @param featuresIN
         * @param label
         * @param deltaMass
         * @param deltaMassType absolute or ppm
         * @param deltaTime
         * @return
         */
    public static ArrayList analyze1(Feature[] featuresIN, IsotopicLabel label,
                                     float deltaMass, int deltaMassType, float deltaTime)
		{
		Feature[] features = (Feature[])featuresIN.clone();
        Arrays.sort(features, Spectrum.comparePeakMzAsc);

		Tree2D tree = new Tree2D();
		for (int i = 0; i < features.length; i++)
			{
			Feature feature = features[i];
			feature.excluded = false;
			tree.add(feature.mass, feature.time, feature);
			}

		ArrayList pairs = new ArrayList();
		ArrayList heavies = new ArrayList();
		for (int i = 0; i < features.length; i++)
			{
			Feature light = features[i];
			if (light.excluded)
				continue;
            float absoluteDeltaMass = calculateAbsoluteDeltaMass(light.mass,
                                                         deltaMass,
                                                         deltaMassType);
            tree.getPoints(light.mass + label.getHeavy()-absoluteDeltaMass,
                           light.time-deltaTime,
			               light.mass+label.getHeavy()+absoluteDeltaMass,
                           light.time+deltaTime,
                           heavies);
			if (heavies.isEmpty())
				continue;
			Feature heavy = null;
			for (Iterator it = heavies.iterator(); it.hasNext();)
				{
				Feature next = (Feature)it.next();
				if (next.charge != light.charge)
					continue;
                heavy = next;
				}
			if (null == heavy)
				continue;
			heavy.excluded = true;
			Pair pair = new Pair(light, heavy);
			pairs.add(pair);
			}
		return pairs; // (Pair[])pairs.toArray(new Pair[0]);
		}

        /**
         * Uses default values for mass and time tolerance
         * @param featuresIN
         * @param label
         * @return
         */
    public static ArrayList analyze(Feature[] featuresIN, IsotopicLabel label)
    {
        return analyze(featuresIN, label, DEFAULT_DELTA_MASS, DELTA_MASS_ABSOLUTE,
                        DEFAULT_DELTA_TIME);
    }


        /**In some way this seems a lot like deconvolute (FeatureSet.deconvolute()).
	     * However, in that case we can group features by mass.  In this case, we don't
	     * know the adjusted mass until after we align (we need to know the number of
	     * modifications).
         *
         * @param featuresIN
         * @param label
         * @param deltaMass
         * @param deltaMassType absolute or ppm
         * @param deltaTime
         * @return
         */
    public static ArrayList analyze(Feature[] featuresIN, IsotopicLabel label,
                                    float deltaMass, int deltaMassType, float deltaTime)
		{
        Feature[] features = (Feature[])featuresIN.clone();
        Arrays.sort(features, Spectrum.comparePeakIntensityDesc);

		Tree2D tree = new Tree2D();

        for (int i = 0; i < features.length; i++)
			{
			Feature feature = features[i];
			feature.excluded = false;
			tree.add(feature.mass, feature.time, feature);
			}

		ArrayList pairs = new ArrayList();
		ArrayList candidates = new ArrayList();
		ArrayList list = new ArrayList();
		for (int i = 0; i < features.length; i++)
			{
			Feature anchor = features[i];
			if (anchor.excluded)
				continue;

			float mass = anchor.mass;
			float time = anchor.time;
			int charge = anchor.charge; // match by charge, if data is deconvoluted then all charge==1

			// we don't know if we have heavy or light
			candidates.clear();
			for (int count=-label.maxLabelCount ; count<=label.maxLabelCount ; count++)
				{
                if (count == 0) continue;
                //deltaMass may be specified in daltons or in ppm.  If ppm, need to calculate
                //the absolute deltaMass to check around this mass
                float absoluteDeltaMass = calculateAbsoluteDeltaMass(mass,
                                                         deltaMass,
                                                         deltaMassType);
                tree.getPoints(mass + count*label.getHeavy()-absoluteDeltaMass,
                           time-deltaTime,
			               mass+count*label.getHeavy()+absoluteDeltaMass,
                           time+deltaTime,
                           list);
                if (list.isEmpty())
	                continue;
				Feature candidate = null;
				for (Iterator it = list.iterator(); it.hasNext();)
					{
                    Feature next = (Feature)it.next();
					if (next.charge != charge)
						continue;
					if (next.excluded)
						continue;
					if (null == candidate || next.totalIntensity > candidate.totalIntensity)
						candidate = next;
					}
				if (null != candidate)
					candidates.add(candidate);
				}
			if (candidates.isEmpty())
				continue;
			candidates.add(anchor);
			Collections.sort(candidates, Spectrum.comparePeakMassAsc);

			// 99% of the time there is exactly one pair
			Pair pair = null;
			if (candidates.size() == 2)
				{
				pair = new Pair(candidates.get(0), candidates.get(1));
				}
			else
				{
				// what range of values do we have
				Feature first = (Feature)candidates.get(0);
				Feature last = (Feature)candidates.get(candidates.size()-1);
				int range = Math.round((last.mass-first.mass)/label.getHeavy());
				assert(range+1 >= candidates.size());

				// slot them into position
				int anchorIndex=0;
				Feature[] slots = new Feature[range+1];
				for (int f=0 ; f<candidates.size() ; f++)
					{
					Feature feature = (Feature)candidates.get(f);
//					System.out.println("" + (feature.mass-anchor.mass) + "\t" + feature.toString());
					int s = Math.round((feature.mass-first.mass)/label.getHeavy());
					slots[s] = feature;
					if (feature == anchor)
						anchorIndex = f;
					}

				// If we had peptide identifications, this would be different,
				// we're just guessing here...

				// OK, now what.  Maybe they pair up nicely
                if (range+1 == candidates.size() && candidates.size()%2==0)
	                {
	                int i1 = (anchorIndex/2)*2;
	                int i2 = i1+1;
	                pair = new Pair(slots[i1],slots[i2]);
	                assert anchor == pair.first || anchor == pair.second;
	                }
				else
	                {
					// OK they don't pair up nicely, there's no guarantee that our anchor peak is the best match here.
	                // I'm going to use a completely arbitrary distance score to find the best match
	                // if it doesn't contain the anchor throw them out and try again
	                int remaining = candidates.size();
                    while (remaining > 1)
	                    {
	                    // find best pair, with some hokey distance measure
	                    // This case doesn't happen often, therefore this is not tested very much either
	                    int bestA = -1, bestB = -1;
	                    double bestDist = Double.MAX_VALUE;
	                    for (int a=0 ; a<range ; a++)
		                    {
		                    if (slots[a] == null) continue;
		                    for (int b=a+1 ; b<= range; b++)
			                    {
			                    if (slots[b] == null) continue;
			                    Feature f = slots[a];
			                    Feature g = slots[b];
			                    double dist =   Math.pow(a-b,2) +
			                                    Math.pow(f.time-g.time,2) +
			                                    Math.pow(Math.log(f.totalIntensity)-Math.log(g.totalIntensity),2);
			                    if (dist < bestDist)
				                    {
				                    bestDist = dist;
				                    bestA = a;
				                    bestB = b;
				                    }
			                    }
			                }
	                    if (slots[bestA] == anchor || slots[bestB] == anchor)
		                    {
		                    pair = new Pair(slots[bestA],slots[bestB]);
		                    break;
		                    }
	                    slots[bestA] = null;
	                    slots[bestB] = null;
	                    remaining-=2;
	                    }
	                }
				}

			if (null == pair)
				{
				anchor.excluded = true;
				}
			else
				{
				Feature f = (Feature)pair.first;
				Feature g = (Feature)pair.second;
				f.excluded = true;
				g.excluded = true;
				if (g.mass < f.mass)
					{
					pair.first = g;
					pair.second = f;
					}
				pairs.add(pair);
                }
			}
		return pairs;
		}


	public static void main(String[] args) throws Exception
		{
		//MSRun run = MSRun.load("E:\\mzxml\\3protmix\\032505KD013_ICAT_3mix.mzXML");
		FeatureSet fs = new FeatureSet(new File("E:\\mzxml\\3protmix\\032505KD013_ICAT_3mix.peptides.tsv"), Color.BLACK);
		ArrayList list = analyze(fs.getFeatures(), icatLabel);
		System.out.println("mz\tscan\tmass\ttotalIntensity\tdeltaMasss\ta/a+b");
		for (int i = 0; i < list.size(); i++)
			{
			Pair pair = (Pair)list.get(i);
			Feature a = (Feature)pair.first;
			Feature b = (Feature)pair.second;
			System.out.println("" + a.mz + "\t" + a.scan + "\t" + a.mass + "\t" + (a.totalIntensity) + "\t" + (b.mass-a.mass) + "\t" + (a.totalIntensity / b.totalIntensity));
			}
		}

    /**
     * Returns an ArrayList of pre-populated IsotopicLabels.
     * If the list isn't yet populated, populate it here
     *
     * @return an arraylist of pre-populated IsotopicLabels
     */
    public static ArrayList getPrePopulatedLabels()
        {
        if (prePopulatedTags == null)
            {
            prePopulatedTags = new ArrayList();
            prePopulatedTags.add(new IsotopicLabel("Cleavable ICAT",227.1263F,9.0297F,'C',
                    DEFAULT_MAX_LABEL_COUNT));
            prePopulatedTags.add(new IsotopicLabel("O18",0.0F,4.0085F,' ',1));
            prePopulatedTags.add(new IsotopicLabel("Silac",0.0F,6.0201F,'K',
                    DEFAULT_MAX_LABEL_COUNT));
            prePopulatedTags.add(new IsotopicLabel("N-terminal",100.4F,4.0313F,
                    'K',DEFAULT_MAX_LABEL_COUNT));
            prePopulatedTags.add(new IsotopicLabel("Acrylamide (D0/D3)",71.03657F,3.0188F,'C',
                    DEFAULT_MAX_LABEL_COUNT));
            }
        return prePopulatedTags;
        }
    }
