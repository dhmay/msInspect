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
package org.fhcrc.cpl.viewer.gui;

import com.sun.j3d.utils.universe.SimpleUniverse;

import javax.media.j3d.*;
import javax.swing.*;
import javax.vecmath.Color3f;
import javax.vecmath.Matrix4f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3f;
import java.awt.*;

public class SurfaceFrame extends JPanel
	{
	private float[][] _matrix = null;
	private SimpleUniverse u = null;


	// Define colors
	Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
	Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
	Color3f red = new Color3f(0.80f, 0.20f, 0.2f);
//	Color3f ambient = new Color3f(0.25f, 0.25f, 0.25f);
	Color3f ambient = new Color3f(0.3f, 0.25f, 0.25f);
	Color3f diffuse = new Color3f(0.7f, 0.7f, 0.7f);
	Color3f specular = new Color3f(0.9f, 0.9f, 0.9f);
	Color3f ambientRed = new Color3f(0.2f, 0.05f, 0.0f);
//	Color3f bgColor = new Color3f(0.05f, 0.05f, 0.3f);
	Color3f bgColor = new Color3f(0.8f, 0.8f, 1f);


	public BranchGroup createSceneGraph()
		{
		// Create the root of the branch graph
		BranchGroup scene = new BranchGroup();

		// Create the bounding leaf node
		BoundingSphere bounds =
				new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
		BoundingLeaf boundingLeaf = new BoundingLeaf(bounds);
		scene.addChild(boundingLeaf);

		// Create the background
		Background bg = new Background(bgColor);
		bg.setApplicationBounds(bounds);
		scene.addChild(bg);

		// Create the ambient light
		AmbientLight ambLight = new AmbientLight(white);
		ambLight.setInfluencingBounds(bounds);
//		scene.addChild(ambLight);

		// Create the directional light
		Vector3f dir = new Vector3f(-1.0f, -1.0f, -1.0f);
		DirectionalLight dirLight = new DirectionalLight(white, dir);
		dirLight.setInfluencingBounds(bounds);
		scene.addChild(dirLight);


		// Create the TransformGroup node and initialize it to the
		// identity. Enable the TRANSFORM_WRITE capability so that
		// our behavior code can modify it at run time. Add it to
		// the root of the subgraph.
		TransformGroup objTrans = new TransformGroup();
		objTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		scene.addChild(objTrans);

		// Create a simple Shape3D node; add it to the scene graph.
//		objTrans.addChild(new ColorCube(0.4));

		Material surfaceMaterial =
				//new Material(ambient, black, diffuse, specular, 110.f);
				new Material();
		surfaceMaterial.setLightingEnable(true);
		PolygonAttributes attr = new PolygonAttributes(); //appearance.getPolygonAttributes();
		attr.setCullFace(PolygonAttributes.CULL_NONE);
		Appearance appearance = new Appearance();
		appearance.setPolygonAttributes(attr);
//		appearance.setMaterial(surfaceMaterial);
		Surface s = new Surface(_matrix, appearance);
		Transform3D xyz = new Transform3D();
		// swap around axes
//		xyz.set(new Matrix3f(0F,1F,0F, 0F,0F,1F, 1F,0F,0F));
		xyz.set(new Matrix4f(0F, 1F, 0F, 0, 0F, 0F, 1F, -.5F, 1F, 0F, 0F, 0F, 0, 0, 0, 1));
		TransformGroup surface = new TransformGroup(xyz);
		surface.addChild(s.getShape());


		objTrans.addChild(surface);


		// Create a new Behavior object that will perform the
		// desired operation on the specified transform and add
		// it into the scene graph.
		if (true)
			{
			Transform3D yAxis = new Transform3D();
			Alpha rotationAlpha = new Alpha(-1, 4000);

			RotationInterpolator rotator =
					new RotationInterpolator(rotationAlpha, objTrans, yAxis,
							0.0f, (float)Math.PI * 2.0f);
			rotator.setSchedulingBounds(bounds);
			scene.addChild(rotator);
			}

		// Have Java 3D perform optimizations on this scene graph.
		scene.compile();

		return scene;
		}


	public SurfaceFrame(float[][] matrix)
		{
		_matrix = matrix;
		}


	public void init()
		{
		setLayout(new BorderLayout());
		GraphicsConfiguration config =
				SimpleUniverse.getPreferredConfiguration();

		Canvas3D c = new Canvas3D(config);
		add("Center", c);

		// Create a simple scene and attach it to the virtual universe
		BranchGroup scene = createSceneGraph();
		u = new SimpleUniverse(c);

		// This will move the ViewPlatform back a bit so the
		// objects in the scene can be viewed.
		u.getViewingPlatform().setNominalViewingTransform();

		u.addBranchGraph(scene);
		}


	public void destroy()
		{
		u.cleanup();
		}


	public static JFrame ShowSurfaceFrame(float[][] matrix)
		{
		JFrame frame = new JFrame("3D");
		frame.setSize(640,400);
		SurfaceFrame s = new SurfaceFrame(matrix);
		s.init();
		//new MainFrame(new SurfaceFrame(matrix), 640, 400);
		frame.getContentPane().add(s);
		return frame;
		}


	//
	// The following allows HelloUniverse to be run as an application
	// as well as an applet
	//
	public static void main(String[] args)
		{
		float[][] m = new float[][]
			{
				{1, 2, 3, 4},
				{1, 2, 1, 2},
				{2, 3, 2, 3},
				{4, 5, 5, 4}
			};
		JFrame frame = ShowSurfaceFrame(m);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
		}
	}
