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

import org.fhcrc.cpl.toolbox.chem.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.reaction.type.AdductionProtonPBReaction;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;


/**
 * Represents the reduction of any double bond between two Carbons and addition of hydrogens.
 * Checks if formula has a double bond somewhere.
 *
 * Triple, quadruple bonds not supported.
 *
*/
public class ReduceDoubleBondAddWaterMod implements ChemicalModification
{
    public ReduceDoubleBondAddWaterMod()
    {
    }


    /**
     * 
     * @param adduct
     * @return
     */
    public void perform(Adduct adduct)
    {
        IMolecule cdkMolecule = adduct.getMolecule();
        try
        {
            List<List<RMap>> bondMappings = UniversalIsomorphismTester.getSubgraphMaps(cdkMolecule,
                    ReduceDoubleBondAdd2HMod.cDoubleBondC);
            int bondIdInMolecule = bondMappings.get(0).get(0).getId1();
            IBond bondInMolecule = cdkMolecule.getBond(bondIdInMolecule);

            bondInMolecule.setOrder(IBond.Order.SINGLE);

            for (int i=0; i<=1; i++)
            {


                ReduceDoubleBondAdd2HMod.unsetAtomProperties(bondInMolecule.getAtom(i));
            }

            //Add OH to one end
            IAtom atom1 = bondInMolecule.getAtom(0);
            IAtom newOxygen = cdkMolecule.getBuilder().newAtom("O");
            cdkMolecule.addAtom(newOxygen);
            IAtom newHydrogen1 = cdkMolecule.getBuilder().newAtom("H");
            cdkMolecule.addAtom(newHydrogen1);
            cdkMolecule.addBond(cdkMolecule.getBuilder().newBond(atom1, newOxygen, IBond.Order.SINGLE));
            cdkMolecule.addBond(cdkMolecule.getBuilder().newBond(newOxygen, newHydrogen1, IBond.Order.SINGLE));

            //Add H to other end
            IAtom atom2 = bondInMolecule.getAtom(1);
            IAtom newHydrogen2 = cdkMolecule.getBuilder().newAtom("H");
            cdkMolecule.addAtom(newHydrogen2);
            cdkMolecule.addBond(cdkMolecule.getBuilder().newBond(atom2, newHydrogen2, IBond.Order.SINGLE));            


            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cdkMolecule);
            adduct.updateFormula();
            adduct.getModifications().add(this);            
        }
        catch (CDKException e)
        {
            ApplicationContext.errorMessage("Reaction failure",e);
        }

    }

    /**
     * Just checks whether there's a double bond to reduce.
     * @param adduct
     * @return
     */
    public boolean canPerform(Adduct adduct)
    {
        //if we don't have a molecule, fail
        IMolecule molecule = adduct.getMolecule();
        if (molecule == null)
            return false;

        boolean canPerform = false;

        try
        {
            //if this fails in some way, we can't perform the modification
            canPerform = UniversalIsomorphismTester.isSubgraph(molecule, ReduceDoubleBondAdd2HMod.cDoubleBondC);
        }
        catch (CDKException e)
        {}


        return canPerform;
    }

    public String getSymbol()
    {
        return "+ H2O";
    }

    public String getName()
    {
        return "DoubleBondReductionWaterAddition";
    }
}
