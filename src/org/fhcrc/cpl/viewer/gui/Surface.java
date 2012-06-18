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

import javax.media.j3d.QuadArray;
import javax.media.j3d.Shape3D;
import javax.media.j3d.Appearance;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;
import javax.vecmath.Color3f;

/**
 * User: mbellew
 * Date: Feb 12, 2005
 * Time: 9:52:40 AM
 */
public class Surface
	{
	Shape3D _shape;

	static int[] topoRGB = {
		0x4C00FF, 0x4900FF, 0x4500FF, 0x4200FF, 0x3E00FF, 0x3B00FF, 0x3700FF,
		0x3300FF, 0x3000FF, 0x2C00FF, 0x2800FF, 0x2500FF, 0x2100FF, 0x1E00FF,
		0x1A00FF, 0x1600FF, 0x1300FF, 0x0F00FF, 0x0C00FF, 0x0800FF, 0x0400FF,
		0x0100FF, 0x0003FF, 0x0006FF, 0x000AFF, 0x000DFF, 0x0011FF, 0x0015FF,
		0x0018FF, 0x001CFF, 0x001FFF, 0x0023FF, 0x0027FF, 0x002AFF, 0x002EFF,
		0x0031FF, 0x0035FF, 0x0039FF, 0x003CFF, 0x0040FF, 0x0043FF, 0x0047FF,
		0x004BFF, 0x004EFF, 0x0052FF, 0x0055FF, 0x0059FF, 0x005DFF, 0x0060FF,
		0x0064FF, 0x0067FF, 0x006BFF, 0x006FFF, 0x0072FF, 0x0076FF, 0x007AFF,
		0x007DFF, 0x0081FF, 0x0084FF, 0x0088FF, 0x008BFF, 0x008FFF, 0x0093FF,
		0x0096FF, 0x009AFF, 0x009EFF, 0x00A1FF, 0x00A5FF, 0x00A8FF, 0x00ACFF,
		0x00AFFF, 0x00B3FF, 0x00B7FF, 0x00BAFF, 0x00BEFF, 0x00C1FF, 0x00C5FF,
		0x00C9FF, 0x00CCFF, 0x00D0FF, 0x00D3FF, 0x00D7FF, 0x00DBFF, 0x00DEFF,
		0x00E2FF, 0x00E5FF, 0x00FF4D, 0x00FF49, 0x00FF45, 0x00FF42, 0x00FF3E,
		0x00FF3A, 0x00FF37, 0x00FF33, 0x00FF2F, 0x00FF2C, 0x00FF28, 0x00FF24,
		0x00FF21, 0x00FF1D, 0x00FF1A, 0x00FF16, 0x00FF12, 0x00FF0F, 0x00FF0B,
		0x00FF07, 0x00FF04, 0x00FF00, 0x04FF00, 0x07FF00, 0x0BFF00, 0x0FFF00,
		0x12FF00, 0x16FF00, 0x1AFF00, 0x1DFF00, 0x21FF00, 0x24FF00, 0x28FF00,
		0x2CFF00, 0x2FFF00, 0x33FF00, 0x37FF00, 0x3AFF00, 0x3EFF00, 0x42FF00,
		0x45FF00, 0x49FF00, 0x4DFF00, 0x50FF00, 0x54FF00, 0x57FF00, 0x5BFF00,
		0x5FFF00, 0x62FF00, 0x66FF00, 0x6AFF00, 0x6DFF00, 0x71FF00, 0x75FF00,
		0x78FF00, 0x7CFF00, 0x80FF00, 0x83FF00, 0x87FF00, 0x8AFF00, 0x8EFF00,
		0x92FF00, 0x95FF00, 0x99FF00, 0x9DFF00, 0xA0FF00, 0xA4FF00, 0xA8FF00,
		0xABFF00, 0xAFFF00, 0xB3FF00, 0xB6FF00, 0xBAFF00, 0xBDFF00, 0xC1FF00,
		0xC5FF00, 0xC8FF00, 0xCCFF00, 0xD0FF00, 0xD3FF00, 0xD7FF00, 0xDBFF00,
		0xDEFF00, 0xE2FF00, 0xE6FF00, 0xFFFF00, 0xFFFE02, 0xFFFD04, 0xFFFB06,
		0xFFFA08, 0xFFF90B, 0xFFF80D, 0xFFF70F, 0xFFF611, 0xFFF513, 0xFFF415,
		0xFFF317, 0xFFF219, 0xFFF11C, 0xFFF01E, 0xFFEF20, 0xFFEE22, 0xFFED24,
		0xFFEC26, 0xFFEC28, 0xFFEB2A, 0xFFEA2D, 0xFFE92F, 0xFFE831, 0xFFE833,
		0xFFE735, 0xFFE637, 0xFFE639, 0xFFE53C, 0xFFE43E, 0xFFE440, 0xFFE342,
		0xFFE344, 0xFFE246, 0xFFE148, 0xFFE14A, 0xFFE04D, 0xFFE04F, 0xFFDF51,
		0xFFDF53, 0xFFDF55, 0xFFDE57, 0xFFDE59, 0xFFDD5B, 0xFFDD5E, 0xFFDD60,
		0xFFDD62, 0xFFDC64, 0xFFDC66, 0xFFDC68, 0xFFDC6A, 0xFFDB6C, 0xFFDB6F,
		0xFFDB71, 0xFFDB73, 0xFFDB75, 0xFFDB77, 0xFFDB79, 0xFFDB7B, 0xFFDB7D,
		0xFFDB80, 0xFFDB82, 0xFFDB84, 0xFFDB86, 0xFFDB88, 0xFFDB8A, 0xFFDB8C,
		0xFFDB8E, 0xFFDB90, 0xFFDB93, 0xFFDC95, 0xFFDC97, 0xFFDC99, 0xFFDC9B,
		0xFFDD9D, 0xFFDD9F, 0xFFDDA1, 0xFFDDA4, 0xFFDEA6, 0xFFDEA8, 0xFFDFAA,
		0xFFDFAC, 0xFFDFAE, 0xFFE0B0, 0xFFE0B2};
	static float[][] colorMap = new float[topoRGB.length][];

	static
		{
	    for (int i=0 ; i<topoRGB.length ; i++)
			{
			int c = topoRGB[i];
			float r = (0xff & (c>>16)) / 255F;
			float g = (0xff & (c>>8)) / 255F;
			float b = (0xff & (c)) / 255F;
			colorMap[i] = new float[] {r,g,b};
			}
		}


	float scaleX(float x, float width)
		{
        return x/width - .5F;
		}

	float scaleY(float y, float height)
		{
        return y/height - .5F;
		}

	float scaleZ(float z, float max)
		{
		return z/max;
		}


	public Surface(float[][] m, Appearance appearance)
		{
		// compute vertices
		int width = m.length;
		int height = m[0].length;
        float maxZ = 1;

		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++)
				maxZ = Math.max(maxZ, m[x][y]);

		// compute a 3d point for each point in the matrix
		Point3f[][] points = new Point3f[width][height];
		for (int x = 0; x < width; x++)
			{
			for (int y = 0; y < height; y++)
				{
				points[x][y] = new Point3f(
						scaleX(x,width),
						scaleY(y,height),
						scaleZ(m[x][y],maxZ));
				}
			}
		Vector3f[][] vectors = new Vector3f[width][height];
		for (int x = 0; x < width; x++)
			{
			for (int y = 0; y < height; y++)
				{
				Point3f x1 = x>0 ? points[x-1][y] : points[x][y];
				Point3f x2 = x<width-1 ? points[x+1][y] : points[x][y];
				Point3f y1 = y>0 ? points[x][y-1] : points[x][y];
				Point3f y2 = y<height-1 ? points[x][y+1] : points[x][y];
				float dzX = x2.z - x1.z;
				if (x>0 && x<width-1)
					dzX /= 2;
				float dzY = y2.z - y1.z;
				if (y>0 && y<height-1)
					dzY /= 2;

				// a = (1, 0, dzX)
				// b = (0, 1, dzY)
				// n.x = a.y * b.z - a.z * b.y = -a.z
				// n.y = a.z * b.x - a.x * b.z = -b.z
				// n.z = a.x * b.y - a.y * b.x = 1
			    float nX = -dzX;
				float nY = -dzY;
				float nZ = 1;
				double nLen = Math.sqrt(nX*nX+nY*nY+1);
				nX /= nLen;
				nY /= nLen;
				nZ /= nLen;
				vectors[x][y] = new Vector3f((float)(nX/nLen), (float)(nY/nLen), (float)(nZ/nLen));
				}
			}
		Color3f[][] colors = new Color3f[width][height];
		for (int x = 0; x < width; x++)
			{
			for (int y = 0; y < height; y++)
				{
				double z = points[x][y].z;
//				z = Math.log1p(z) / Math.log1p(1);
//				z = Math.max(0F,Math.min(1F,z));
//				colors[x][y] = new Color3f((float)z,(float)z,(float)z);
				z = (float)(Math.log(z+1000)/Math.log(1001) * (colors.length-2));
				z = Math.max(0, Math.min(z, colorMap.length-2));
				int i = (int)Math.floor(z);
				double r = z-i;
				colors[x][y] = new Color3f(
						(float)(colorMap[i][0] * (1-r) + colorMap[i+1][0] * r),
						(float)(colorMap[i][1] * (1-r) + colorMap[i+1][1] * r),
						(float)(colorMap[i][2] * (1-r) + colorMap[i+1][2] * r));
				}
			}

		Point3f[] verts = new Point3f[4 * (width-1) * (height-1)];
		Vector3f[] normals = new Vector3f[4 * (width-1) * (height-1)];
		Color3f[] clrs = new Color3f[4 * (width-1) * (height-1)];
		int v = 0, n=0, c=0;
		for (int x = 0; x < width - 1; x++)
			{
			for (int y = 0; y < height - 1; y++)
				{
				verts[v++] = points[x][y];
				verts[v++] = points[x + 1][y];
				verts[v++] = points[x + 1][y + 1];
				verts[v++] = points[x][y + 1];
				normals[n++] = vectors[x][y];
				normals[n++] = vectors[x+1][y+1];
				normals[n++] = vectors[x+1][y];
				normals[n++] = vectors[x][y+1];
				clrs[c++] = colors[x][y];
				clrs[c++] = colors[x+1][y];
				clrs[c++] = colors[x+1][y+1];
				clrs[c++] = colors[x][y+1];
				}
			}

		QuadArray quad = new QuadArray(verts.length, QuadArray.COORDINATES|QuadArray.NORMALS|QuadArray.COLOR_3);
		quad.setCoordinates(0, verts);
		quad.setNormals(0, normals);
		quad.setColors(0, clrs);

		_shape = new Shape3D(quad, appearance);
		}


	public Shape3D getShape()
		{
		return _shape;
		}
	}
