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

/* Rajarshi Guha <rajarshi@presidency.com>
 * 1/11/2004
 */
package org.fhcrc.cpl.viewer.metabologna;

import java.awt.image.BufferedImage;
import java.awt.*;
import java.util.List;
import java.util.ArrayList;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.Renderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.IAtomContainerGenerator;

public class MoleculeRenderer2D
{
    protected int width = 400;
    protected int height = 400;

    protected boolean shouldShowHydrogens = false;
    protected boolean shouldShowCarbons = false;

    public Image renderMolecule(IMolecule molecule)
    {
        // the draw area and the image should be the same size
        Rectangle drawArea = new Rectangle(width, height);
        Image image = new BufferedImage(
                width, height, BufferedImage.TYPE_INT_RGB);

        StructureDiagramGenerator sdg = new StructureDiagramGenerator();       
        sdg.setMolecule(molecule);
        try
        {
            sdg.generateCoordinates();
        }
        catch (Exception e)
        {}
        molecule = sdg.getMolecule();
        
        // generators make the image elements
        List<IAtomContainerGenerator> generators = new ArrayList<IAtomContainerGenerator>();
        generators.add(new BasicBondGenerator());
        MyAtomGenerator atomGenerator = new MyAtomGenerator(shouldShowHydrogens, shouldShowCarbons);
        generators.add(atomGenerator);

        // the renderer needs to have a toolkit-specific font manager
        Renderer renderer = new Renderer(generators, new AWTFontManager());

        // the call to 'setup' only needs to be done on the first paint
        renderer.setup(molecule, drawArea);

        // paint the background
        Graphics2D g2 = (Graphics2D)image.getGraphics();
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, width, height);

        // the paint method also needs a toolkit-specific renderer

        renderer.paintMolecule(molecule, new AWTDrawVisitor(g2), new Rectangle(0, 0, width, height), false);
        return image;
    }

    public int getWidth() {
        return width;
    }

    public void setWidth(int width) {
        this.width = width;
    }

    public int getHeight() {
        return height;
    }

    public void setHeight(int height) {
        this.height = height;
    }

    public boolean isShouldShowHydrogens() {
        return shouldShowHydrogens;
    }

    /**
     * Must provide a molecule with explicit hydrogens in order for hydrogens to be shown
     * @param shouldShowHydrogens
     */
    public void setShouldShowHydrogens(boolean shouldShowHydrogens) {
        this.shouldShowHydrogens = shouldShowHydrogens;
    }

    public void setShouldShowCarbons(boolean shouldShowCarbons) {
        this.shouldShowCarbons = shouldShowCarbons;
    }

    protected class MyAtomGenerator extends BasicAtomGenerator
    {
        protected boolean shouldShowHydrogens = false;
        protected boolean shouldShowCarbons = false;

        public MyAtomGenerator(boolean h, boolean c)
        {
            super();
            shouldShowHydrogens = h;
            shouldShowCarbons = c;
        }

        public boolean invisibleCarbon(IAtom atom, IAtomContainer container, RendererModel model)
        {
            return (isCarbon(atom) && !shouldShowCarbons);
        }

        public boolean invisibleHydrogen(IAtom atom, RendererModel model)
        {
            return (isHydrogen(atom) && !shouldShowHydrogens);
        }
    }


}

