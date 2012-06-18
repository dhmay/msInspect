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
package org.fhcrc.cpl.toolbox.gui;

/**
 * User: mbellew
 * Date: May 24, 2004
 * Time: 10:33:45 AM
 */

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;


public class ScrollableImage extends JComponent
		implements Scrollable,
		MouseMotionListener
	{
	String name = ""; // for debugging only
	protected int maxUnitIncrement = 1;
	protected Image image = null;

	public ScrollableImage(String name)
		{
		this.name = name;
		setDoubleBuffered(false);
		setAutoscrolls(true);
		setOpaque(true);
		setBackground(Color.white);
		addMouseMotionListener(this);
		}


	public ScrollableImage(String name, Image i, int m)
		{
		this(name);
		setMaxUnitIncrement(m);
		setImage(i);
		}



	public void setImage(Image i)
		{
		image = i;
		if (null == i)
			setPreferredSize(new Dimension(0,0));
		else
			setPreferredSize(new Dimension(image.getWidth(this), image.getHeight(this)));

		invalidate();
		repaint();
		}


	public Image getImage()
		{
		return image;
		}


	public void paint(Graphics graphics)
		{
		if (null == image)
			{
			graphics.clearRect(0, 0, getWidth(), getHeight());
			return;
			}

		// UNDONE: Smart scrolling is broken w/ multiple monitors (JDK 1.4)
		if (getParent() instanceof JViewport)
			{
			if (JViewport.SIMPLE_SCROLL_MODE != ((JViewport) getParent()).getScrollMode())
				{
				((JViewport) getParent()).setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
				repaint();
				}
			}
		graphics.drawImage(image, 0, 0, this);
		}


	public void mouseMoved(MouseEvent e)
		{
		}


	public void mouseDragged(MouseEvent e)
		{
		Rectangle r = new Rectangle(e.getX(), e.getY(), 1, 1);
		scrollRectToVisible(r);
		}


	/*
	public Dimension getPreferredSize()
		{
		if (null == image)
			{
			return new Dimension(-1, -1);
			}
		else
			{
			return new Dimension(image.getWidth(this), image.getHeight(this));
			}
		}
    */


	public Dimension getPreferredScrollableViewportSize()
		{
		return getPreferredSize();
		}


	public int getScrollableUnitIncrement(Rectangle visibleRect,
			int orientation,
			int direction)
		{
		int currentPosition = 0;
		if (orientation == SwingConstants.HORIZONTAL)
			{
			currentPosition = visibleRect.x;
			}
		else
			{
			currentPosition = visibleRect.y;
			}

		if (direction < 0)
			{
			int newPosition = currentPosition -
					(currentPosition / maxUnitIncrement)
					* maxUnitIncrement;
			return (newPosition == 0) ? maxUnitIncrement : newPosition;
			}
		else
			{
			return ((currentPosition / maxUnitIncrement) + 1)
					* maxUnitIncrement
					- currentPosition;
			}
		}


	public int getScrollableBlockIncrement(Rectangle visibleRect,
			int orientation,
			int direction)
		{
		if (orientation == SwingConstants.HORIZONTAL)
			{
			return visibleRect.width - maxUnitIncrement;
			}
		else
			{
			return visibleRect.height - maxUnitIncrement;
			}
		}


	public boolean getScrollableTracksViewportWidth()
		{
		return false;
		}


	public boolean getScrollableTracksViewportHeight()
		{
		return false;
		}


	public void setMaxUnitIncrement(int pixels)
		{
		maxUnitIncrement = pixels;
		}
	}

