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
public class ReduceDoubleBondAdd2HMod implements ChemicalModification
{
    //C=C
    public static IMolecule cDoubleBondC;

    static
    {
        try
        {
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            cDoubleBondC = sp.parseSmiles("C=C");
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cDoubleBondC);

            LonePairElectronChecker lpcheck = new LonePairElectronChecker();
            lpcheck.saturate(cDoubleBondC);
        }
        catch (CDKException e)
        {
            e.printStackTrace(System.err);
        }
    }

    public ReduceDoubleBondAdd2HMod()
    {
    }

    public static void unsetAtomProperties(IAtom atom)
    {
        //Remove properties from the atom so that they can be inferred
        //Not all of these are necessary
        atom.setAtomTypeName((String) CDKConstants.UNSET);
        atom.setMaxBondOrder((IBond.Order) CDKConstants.UNSET);
        atom.setBondOrderSum((Double) CDKConstants.UNSET);
        atom.setCovalentRadius((Double) CDKConstants.UNSET);
        atom.setValency((Integer) CDKConstants.UNSET);
        atom.setFormalCharge((Integer) CDKConstants.UNSET);
        atom.setHybridization((IAtomType.Hybridization) CDKConstants.UNSET);
        atom.setFormalNeighbourCount((Integer) CDKConstants.UNSET);
        atom.setFlag(CDKConstants.IS_HYDROGENBOND_ACCEPTOR, false);
        atom.setFlag(CDKConstants.IS_HYDROGENBOND_DONOR, false);
        atom.setFlag(CDKConstants.ISAROMATIC, false);
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
            List<List<RMap>> bondMappings = UniversalIsomorphismTester.getSubgraphMaps(cdkMolecule, cDoubleBondC);
            int bondIdInMolecule = bondMappings.get(0).get(0).getId1();
            IBond bondInMolecule = cdkMolecule.getBond(bondIdInMolecule);

            bondInMolecule.setOrder(IBond.Order.SINGLE);

            for (int i=0; i<=1; i++)
            {
                IAtom atom = bondInMolecule.getAtom(i);              
                IAtom newHydrogen = cdkMolecule.getBuilder().newAtom("H");
                cdkMolecule.addAtom(newHydrogen);
                IBond newBond = cdkMolecule.getBuilder().newBond(atom, newHydrogen, IBond.Order.SINGLE);
                cdkMolecule.addBond(newBond);

                unsetAtomProperties(atom);
            }
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cdkMolecule);
            adduct.updateFormula();
            adduct.getModifications().add(this);            
        }
        catch (CDKException e)
        {
            ApplicationContext.errorMessage("Reaction failure",e);
        }

/*
        //assume cdkMolecule not null, i.e., canPerform has returned true
        setOfReactants.addMolecule(cdkMolecule);
        IReactionProcess rectionProcess = new AdductionProtonPBReaction();



        try
        {
            IReactionSet setOfReactions = rectionProcess.initiate(setOfReactants, null);
            IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
            //todo: fix this.  Applying all reactions
            for (IReaction reaction : setOfReactions.reactions())
            {
                IMoleculeSet products = reaction.getProducts();

                //only one product for this reaction, guaranteed
                IMolecule product = products.getMolecule(0);

                for (int i=0; i<product.getAtomCount(); i++)
                {                                  
                    IAtom atom = product.getAtom(i);
                    if (atom.getFlag(CDKConstants.REACTIVE_CENTER))
                    {
                        IAtom newHydrogen = builder.newAtom("H");
                        product.addAtom(newHydrogen);
                        IBond newBond = product.getBuilder().newBond(atom, newHydrogen, IBond.Order.SINGLE);
                        product.addBond(newBond);

                        //Remove properties from the atom so that they can be inferred
                        //Not all of these are necessary
                        atom.setAtomTypeName((String) CDKConstants.UNSET);
                        atom.setMaxBondOrder((IBond.Order) CDKConstants.UNSET);
                        atom.setBondOrderSum((Double) CDKConstants.UNSET);
                        atom.setCovalentRadius((Double) CDKConstants.UNSET);
                        atom.setValency((Integer) CDKConstants.UNSET);
                        atom.setFormalCharge((Integer) CDKConstants.UNSET);
                        atom.setHybridization((IAtomType.Hybridization) CDKConstants.UNSET);
                        atom.setFormalNeighbourCount((Integer) CDKConstants.UNSET);
                        atom.setFlag(CDKConstants.IS_HYDROGENBOND_ACCEPTOR, false);
                        atom.setFlag(CDKConstants.IS_HYDROGENBOND_DONOR, false);
                        atom.setFlag(CDKConstants.ISAROMATIC, false);
                        break;
                    }
                }

//                AtomContainerManipulator.clearAtomConfigurations(product);
                
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(product);
                adduct.setMolecule(product);

            }
        }
        catch (CDKException e)
        {
            ApplicationContext.errorMessage("Reaction failure",e);
        }
//        adduct.setFormula(adduct.getFormula().createFormulaWithAddition(formulaToAdd));
//        adduct.getModifications().add(this);
*/
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
            canPerform = UniversalIsomorphismTester.isSubgraph(molecule, cDoubleBondC);
        }
        catch (CDKException e)
        {}


        return canPerform;
    }

    public String getSymbol()
    {
        return "+ H2";
    }

    public String getName()
    {
        return "DoubleBondReduction2HAddition";
    }
}
