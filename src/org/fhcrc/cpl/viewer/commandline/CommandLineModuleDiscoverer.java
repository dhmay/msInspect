/*
 * Copyright (c) 2003-2007 Fred Hutchinson Cancer Research Center
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

package org.fhcrc.cpl.viewer.commandline;

import org.labkey.common.tools.ApplicationContext;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.jar.JarFile;
import java.net.JarURLConnection;
import java.net.URL;
import java.net.URLDecoder;
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Modifier;

/**
 * Utility class for finding all command-line modules.
 * All _standard_ commandlinemodules should be declared in the commandline_modules
 * resource bundle.
 * See comments on method discoverAllCommandLineModules() for important details.
 */
public class CommandLineModuleDiscoverer
{
    protected static Logger _log = Logger.getLogger(CommandLineModuleDiscoverer.class);


    protected static Map<String,CommandLineModule> allCommandLineModuleMap = null;
    protected static Map<String,CommandLineModule> standardCommandLineModuleMap = null;

    //default values for the resourcebundle package and name
    protected static String _moduleBundlePackage = "org.fhcrc.cpl.viewer.commandline";
    protected static String _moduleBundleName = "commandline_modules";

    //package in which all modules must reside
    protected static final String modulePackageName =
            "org.fhcrc.cpl.viewer.commandline.modules";


    /**
     * Return a map of just the _standard_ command line modules.  This is for performance:
     * only go fishing if we need to
     * @return
     */
    public static Map<String,CommandLineModule> getStandardCommandLineModules()
    {
        if (standardCommandLineModuleMap == null)
        {
            //load standard commandline modules from the resource bundle

            standardCommandLineModuleMap = new HashMap<String,CommandLineModule>();

            ResourceBundle moduleResourceBundle =
                    ResourceBundle.getBundle(_moduleBundlePackage + "." +
                                             _moduleBundleName);
            Enumeration<String> propertyNames = moduleResourceBundle.getKeys();

            String propertyName = null;
            while(propertyNames.hasMoreElements())
            {
                try
                {
                    propertyName = (String) propertyNames.nextElement();
                    if ("true".equalsIgnoreCase(moduleResourceBundle.getString(propertyName)))
                    {
                        String fullClassName = propertyName;
                        if (!fullClassName.contains("."))
                            fullClassName = modulePackageName + "." + propertyName;
                        Class<?> moduleClass = Class.forName(fullClassName);
                        CommandLineModule module =
                                (CommandLineModule) moduleClass.newInstance();
                        standardCommandLineModuleMap.put(module.getCommandName(),
                                module);
                    }
                }
                catch (Exception e)
                {
                    ApplicationContext.setMessage("Failed to load command line module with name " + propertyName);
                    if (_log.isDebugEnabled())
                    {
                        ApplicationContext.setMessage(e.getMessage());
                        e.printStackTrace(System.err);
                    }
                }
            }
        }
        return standardCommandLineModuleMap;
    }

    /**
     * Return all known commandline modules, initializing if necessary
     * @return
     */
    public static Map<String,CommandLineModule> findAllCommandLineModules()
    {
        if (allCommandLineModuleMap == null)
            allCommandLineModuleMap = discoverAllCommandLineModules();
        return allCommandLineModuleMap;
    }

    /**
     * Attempts to discover all the classes in the package org.fhcrc.cpl.viewer.commandline.modules
     * that implement the CommandLineModule interface, as determined by the context class loader
     *
     * Note special handling for One-Jar.  Here's how it works:
     * one-jar has a Very Special classloader called JarClassLoader.  It's Very Special in that
     * you can't interrogate it to find the jar files that the classes are coming from, as resources.
     * You can only locate them as entries on the main jar itself.
     * So, if we're using that classloader, we have to follow a very special convention, which is this:
     *  -Each extension jar contains exactly one extension
     *  -All extension jars go in the lib directory of the one jar
     *  -Each extension jar has the name "extension_ClassName.jar", where ClassName is the name of the
     *   CommandLineModule class contained in the jar
     *
     * @return a map of CommandLineModule-implementing classes that exist within the correct package
     */
    protected static Map<String, CommandLineModule> discoverAllCommandLineModules()
    {
        HashMap<String, CommandLineModule> result = new HashMap<String, CommandLineModule>();


        //populate the standard modules, so they don't get blown away by anything custom.
        //You want custom behavior, create a custom module.  This will keep bug-tracking
        //slightly simpler
        for (String key : getStandardCommandLineModules().keySet())
        {
            result.put(key, getStandardCommandLineModules().get(key));
        }

        //This holds the list of discovered classes
        ArrayList<Class> classesInPackage = new ArrayList<Class>();

        ClassLoader cld = Thread.currentThread().getContextClassLoader();

        if (cld == null)
        {
            ApplicationContext.infoMessage("Failure finding commandline modules: couldn't get classloader");
            return result;
        }

        //Special handling for one-jar classloader.
        if ("com.simontuffs.onejar.JarClassLoader".equals(cld.getClass().getName()))
        {
            try
            {
                //this URL is for the one-jar
                URL jarUrl = cld.getResource("lib");
                JarURLConnection conn = (JarURLConnection) jarUrl.openConnection();
                JarFile jfile = conn.getJarFile();
                java.util.Enumeration e = jfile.entries();
                while (e.hasMoreElements())
                {
                    java.util.zip.ZipEntry entry = (java.util.zip.ZipEntry) e.nextElement();
                    String entryname = entry.getName();
                    if (entryname.startsWith("lib") &&
                        entryname.endsWith(".jar"))
                    {
                        String jarFilename = entryname.substring(entryname.lastIndexOf(System.getProperty("file.separator")) + 1,
                                                                 entryname.lastIndexOf(".jar"));
                        if (jarFilename.startsWith("extension_"))
                        {
                            String extensionClassName = jarFilename.substring("extension_".length());

                            classesInPackage.add(Class.forName(modulePackageName + "." +
                                                               extensionClassName));
                        }
                    }

                }
            }
            catch (Exception e)
            {
                ApplicationContext.infoMessage("Failed while trying to discover custom commandline modules in one-jar.\n" +
                                               "Only standard modules will work.");
            }
        }
        else
        {
        //Regular handling for the non-one-jar situation, in which extensions can be anywhere
        //on the classpath and we can actually discover them

        //This will hold a list of directories matching the modulePackageName. There may be more than
        //one if a package is split over multiple jars/paths
        ArrayList<File> directories = new ArrayList<File>();
        // this will hold a list of jars matching the pckgname.  May be more than one
        ArrayList<URL> jarUrls = new ArrayList<URL>();

        try
        {
            String path = modulePackageName.replace('.', '/');
            // Ask for all resources for the path
            Enumeration<URL> resources = cld.getResources(path);
            while (resources.hasMoreElements())
            {
                URL urlResource = resources.nextElement();
                if (urlResource.getProtocol().equals("jar"))
                {
                    jarUrls.add(urlResource);
                }
                else
                    directories.add(new File(URLDecoder.decode(urlResource.getPath(), "UTF-8")));
            }
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("Failure finding commandline modules", x);
            return result;
        }

        //First, add all the classes in the package to a list
        try
        {
            // For every directory identified capture all the matching .class files
            for (File directory : directories)
            {
                if (directory.exists())
                {
                    // Get the list of the files contained in the package
                    String[] files = directory.list();
                    for (String file : files)
                    {
                        // we are only interested in .class files
                        if (file.endsWith(".class") && !file.contains("$"))
                        {
                            // removes the .class extension
                            classesInPackage.add(Class.forName(modulePackageName + '.' +
                                    file.substring(0, file.length() - ".class".length())));
                        }
                    }
                }
            }
            // For every jar identified capture all the matching .class files
            for (URL jarUrl : jarUrls)
            {
                JarURLConnection conn = (JarURLConnection) jarUrl.openConnection();
                java.util.jar.JarFile jfile = conn.getJarFile();
                String starts = conn.getEntryName();
                java.util.Enumeration e = jfile.entries();
                while (e.hasMoreElements())
                {
                    java.util.zip.ZipEntry entry =
                            (java.util.zip.ZipEntry) e.nextElement();
                    String entryname = entry.getName();
                    if (entryname.startsWith(starts)
                            && (entryname.lastIndexOf('/') <= starts.length())
                            && entryname.endsWith(".class") && !entryname.contains("$"))
                    {
                        String classname =
                                entryname.substring(0, entryname.length() - 6);
                        if (classname.startsWith("/"))
                            classname = classname.substring(1);
                        classname = classname.replace('/', '.');
                        try
                        {
                            // Try to create an instance of the object
                            classesInPackage.add(Class.forName(classname));
                        }
                        catch (Throwable t)
                        {
                        }
                    }
                }
            }
        }
        catch (Exception e)
        {
            ApplicationContext.errorMessage("Failure finding commandline modules, unable to search all directories", e);
            return result;
        }
        }

        //Now, for each class identified, check to make sure it implements the CommandLineModule interface
        for (Class candidateClass : classesInPackage)
        {
            try
            {
                //don't try to instantiate abstract classes
                if (Modifier.isAbstract(candidateClass.getModifiers()))
                    continue;

                Object candidateInstance = candidateClass.newInstance();
                //only add the first discovered instance of a particular command
                if (candidateInstance instanceof CommandLineModule &&
                        !result.containsKey(((CommandLineModule) candidateInstance).getCommandName()))
                {
                    result.put(((CommandLineModule) candidateInstance).getCommandName(),
                            (CommandLineModule) candidateInstance);
                }
            }
            catch (Exception e)
            {
                ApplicationContext.errorMessage("Failure finding commandline module " + candidateClass.getName(), e);
            }
        }
        return result;
    }

    /**
     * Return the commandline module for the specified command.
     * If that module doesn't exist, throw a FileNotFoundExcption.
     * That's a bit of a hack, but it's necessary for an exception to be thrown
     * that's not a RunTimeException, to make the calling code handle it explicitly
     *
     * First try standard modules, then custom.  For performance reasons.
     * @param command
     * @return
     * @throws FileNotFoundException if no module exists for the specified command
     */
    public static CommandLineModule getCommandLineModule(String command)
            throws FileNotFoundException
    {
        CommandLineModule module = null;
        try
        {
            module = getStandardCommandLineModule(command);
        }
        catch (FileNotFoundException e)
        {
            module = findAllCommandLineModules().get(command.toLowerCase());
        }

        if (module == null)
            throw new FileNotFoundException("No commandline Module found for command " + command);
        return module;
    }

    /**
     * Return the commandline module for the specified command.
     * If that module doesn't exist, throw a FileNotFoundExcption.
     * That's a bit of a hack, but it's necessary for an exception to be thrown
     * that's not a RunTimeException, to make the calling code handle it explicitly
     * @param command
     * @return
     * @throws FileNotFoundException if no module exists for the specified command
     */
    public static CommandLineModule getStandardCommandLineModule(String command)
            throws FileNotFoundException
    {
        CommandLineModule module = getStandardCommandLineModules().get(command.toLowerCase());
        if (module == null)
            throw new FileNotFoundException("No standard commandline Module found for command " + command);
        return module;
    }
}
