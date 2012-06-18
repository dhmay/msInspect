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
package org.fhcrc.cpl.viewer.database;

import org.apache.log4j.Logger;

import org.hibernate.*;
import org.hibernate.cfg.*;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.viewer.gui.WorkbenchFileChooser;

import javax.swing.*;
import java.util.prefs.Preferences;
import java.io.File;


/**
 * Command linemodule for plotting the mass calibration of a feature file
 */
public class HibernateManager
{
    protected static Logger _log = Logger.getLogger(HibernateManager.class);

    protected static HibernateManager singletonInstance = null;

    protected SessionFactory sessionFactory;

    protected static File configFile = null;

    protected static boolean dirtySessionFactory = false;

    public static final String HIBERNATE_CONFIG_FILE_PREF_NAME = "hibernateConfigFile";


    public HibernateManager()
    {
        try
        {
            // Create the SessionFactory from hibernate.cfg.xml
            Configuration configuration = new Configuration();
            if (configFile == null)
                retrieveOrPromptConfigFilePreference();
            if (configFile != null)
                configuration.configure(configFile);
            sessionFactory = configuration.buildSessionFactory();
        } catch (Throwable ex) {
            // Make sure you log the exception, as it might be swallowed
            System.err.println("SessionFactory creation failed." + ex);
            throw new ExceptionInInitializerError(ex);
        }
    }

    protected void retrieveOrPromptConfigFilePreference()
    {
        Preferences prefs = Preferences.userNodeForPackage(Application.class);
        String configFilePath = prefs.get(HIBERNATE_CONFIG_FILE_PREF_NAME, null);
        if (configFilePath == null)
            promptConfigFilePreference();
        else
        {
            try
            {
                setConfigurationFile(new File(configFilePath));
            }
            catch (Exception e)
            {
                promptConfigFilePreference();                
            }
        }
    }

    protected void promptConfigFilePreference()
    {
        WorkbenchFileChooser wfc = new WorkbenchFileChooser();
        wfc.setDialogTitle("Choose Hibernate Configuration File");
        int chooserStatus = wfc.showOpenDialog(null);
        if (chooserStatus != JFileChooser.APPROVE_OPTION)
        {
            throw new RuntimeException("You must choose a Hibernate configuration file to continue");
        }

        final File file = wfc.getSelectedFile();
        if (null != file)
        {
            setConfigurationFile(file);
            Preferences prefs = Preferences.userNodeForPackage(Application.class);
            prefs.put(HIBERNATE_CONFIG_FILE_PREF_NAME, file.getAbsolutePath());
        }
        else
        {
            throw new RuntimeException("You must choose a Hibernate configuration file to continue");
        }
    }

    public void close()
    {
        if (sessionFactory != null)
            sessionFactory.close();
    }

    public static void setConfigurationFile(File newConfigFile)
    {
        if (newConfigFile == null || !newConfigFile.exists())
            throw new IllegalArgumentException("Bad config file");
        configFile = newConfigFile;
        dirtySessionFactory = true;
    }

    public static HibernateManager getInstance()
    {
        if (singletonInstance == null || dirtySessionFactory)
        {
            if (dirtySessionFactory && singletonInstance != null)
                singletonInstance.close();
            singletonInstance = new HibernateManager();
            dirtySessionFactory = false;
        }
        return singletonInstance;
    }

    public SessionFactory getSessionFactory()
    {
        return sessionFactory;
    }

    public Session openSession()
    {
        return sessionFactory.getCurrentSession();
    }




}
