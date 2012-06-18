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

import javax.swing.event.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * User: mbellew
 * Date: Sep 17, 2004
 * Time: 7:07:00 PM
 */
public class ListenerHelper
	{
	Object controller;
	Class controllerClass;
	ArrayList listeners;

	protected static final int CHANGE = 1;
	protected static final int PROPERTYCHANGE = 2;
	protected static final int ACTION = 3;
	protected static final int COMPONENT = 4;
	protected static final int MOUSE = 5;
	protected static final int MOUSEMOTION = 6;
	protected static final int KEY = 7;
	protected static final int DOCUMENT = 8;
    protected static final int LIST_SELECTION = 9;

    static class ListenerDescription
		{
		String name;
		Class listenerClass;
		Class eventClass;


		ListenerDescription(String name, Class listenerClass, Class eventClass)
			{
			this.name = name;
			this.listenerClass = listenerClass;
			this.eventClass = eventClass;
			}
		}

	ListenerDescription[] desciptions = new ListenerDescription[]
		{
			null,
			new ListenerDescription("ChangeListener", ChangeListener.class, ChangeEvent.class),
			new ListenerDescription("PropertyChangeListener", PropertyChangeListener.class, PropertyChangeEvent.class),
			new ListenerDescription("ActionListener", ActionListener.class, ActionEvent.class),
			new ListenerDescription("ComponentListener", ComponentListener.class, ComponentEvent.class),
			new ListenerDescription("MouseListener", MouseListener.class, MouseEvent.class),
			new ListenerDescription("MouseMotionListener", MouseMotionListener.class, MouseEvent.class),
			new ListenerDescription("KeyListener", KeyListener.class, KeyEvent.class),
			new ListenerDescription("DocumentListener", DocumentListener.class, DocumentEvent.class),
            new ListenerDescription("ListSelectionListener", ListSelectionListener.class, ListSelectionEvent.class)
        };


	protected static HashMap mapListenMethods = new HashMap();


	static
	{
	mapListenMethods.put("stateChanged", new Integer(CHANGE));
	mapListenMethods.put("propertyChange", new Integer(PROPERTYCHANGE));

	mapListenMethods.put("actionPerformed", new Integer(ACTION));

	mapListenMethods.put("componentResized", new Integer(COMPONENT));
	mapListenMethods.put("componentMoved", new Integer(COMPONENT));
	mapListenMethods.put("componentShown", new Integer(COMPONENT));
	mapListenMethods.put("componentHidden", new Integer(COMPONENT));

	mapListenMethods.put("mouseClicked", new Integer(MOUSE));
	mapListenMethods.put("mousePressed", new Integer(MOUSE));
	mapListenMethods.put("mouseReleased", new Integer(MOUSE));
	mapListenMethods.put("mouseEntered", new Integer(MOUSE));
	mapListenMethods.put("mouseExited", new Integer(MOUSE));

	mapListenMethods.put("mouseDragged", new Integer(MOUSEMOTION));
	mapListenMethods.put("mouseMoved", new Integer(MOUSEMOTION));

	mapListenMethods.put("keyPressed", new Integer(KEY));
	mapListenMethods.put("keyReleased", new Integer(KEY));
	mapListenMethods.put("keyTyped", new Integer(KEY));

	mapListenMethods.put("insertUpdate", new Integer(DOCUMENT));
	mapListenMethods.put("removeUpdate", new Integer(DOCUMENT));
	mapListenMethods.put("changedUpdate", new Integer(DOCUMENT));

    mapListenMethods.put("valueChanged", new Integer(LIST_SELECTION));
    }


	public ListenerHelper(Object controller)
		{
		this.controllerClass = controller.getClass();
		this.controller = controller;
		listeners = new ArrayList();
		}


	public Object addListener(Object component, String method, String propName)
		{
		try
			{
			Class componentClass = component.getClass();
			String eventName = method.substring(method.lastIndexOf('_') + 1);
			int type = ((Integer)mapListenMethods.get(eventName)).intValue();
			Method wrapMethod = null;
			Method addMethod = null;

			ListenerDescription desc = desciptions[type];

			wrapMethod = controllerClass.getMethod(method, new Class[]{desc.eventClass});
			Object listener = newAdapter(component, type, wrapMethod);

			if (type == PROPERTYCHANGE && null != propName)
				{
				addMethod = componentClass.getMethod("addPropertyChangeListener",
						new Class[]{String.class, PropertyChangeListener.class});
				addMethod.invoke(component, new Object[]{propName, listener});
				}
			else
				{
				addMethod = componentClass.getMethod("add" + desc.name, new Class[]{desc.listenerClass});
				addMethod.invoke(component, new Object[]{listener});
				}

			listeners.add(listener);
			return listener;
			}
		catch (NoSuchMethodException x)
			{
			throw new RuntimeException(x);
			}
		catch (IllegalAccessException x)
			{
			throw new RuntimeException(x);
			}
		catch (InvocationTargetException x)
			{
			throw new RuntimeException(x);
			}
		}


	public void addListener(Object component, String method)
		{
		addListener(component, method, null);
		}


	public void removeListeners()
		{
		for (int i = 0; i < listeners.size(); i++)
			{
			Proxy proxy = (Proxy)listeners.get(i);
			ListenerAdapter adapter = (ListenerAdapter)Proxy.getInvocationHandler(proxy);
			Class componentClass = adapter.component.getClass();
			ListenerDescription desc = desciptions[adapter.type];
			try
				{
				Method removeMethod = componentClass.getMethod("remove" + desc.name, new Class[]{desc.listenerClass});
				removeMethod.invoke(adapter.component, new Object[]{adapter});
				}
			catch (NoSuchMethodException x)
				{

				}
			catch (InvocationTargetException x)
				{

				}
			catch (IllegalAccessException x)
				{

				}
			}
		}


	Object newAdapter(Object c, int type, Method m)
		{
		ListenerDescription desc = desciptions[type];
		return java.lang.reflect.Proxy.newProxyInstance(this.getClass().getClassLoader(),
				new Class[]{desc.listenerClass},
				new ListenerAdapter(c, type, m));
		}


	class ListenerAdapter implements
			java.lang.reflect.InvocationHandler
		{
		// so we can remove the listener
		int type;
		Object component;

		// method
		Method wrapMethod;


		private ListenerAdapter(Object c, int t, Method m)
			{
			type = t;
			component = c;
			this.wrapMethod = m;
			}


		public Object invoke(Object o, Method method, Object[] objects) throws Throwable
			{
			if (!wrapMethod.getName().endsWith(method.getName()))
            {
                return null;
            }
			try
				{
				return wrapMethod.invoke(controller, objects);
				}
			catch (Throwable x)
				{
				if (x instanceof InvocationTargetException)
					x = ((InvocationTargetException)x).getTargetException();
				x.printStackTrace();
				throw x;
				}
			}
		}
	}
