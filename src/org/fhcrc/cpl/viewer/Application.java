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
package org.fhcrc.cpl.viewer;

import org.fhcrc.cpl.viewer.util.ConvertHelper;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.viewer.gui.OpenFileAction;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.fhcrc.cpl.viewer.gui.SpecifyModuleArgumentsFrame;
import org.fhcrc.cpl.toolbox.gui.AwtPropertyBag;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.viewer.commandline.*;
import org.fhcrc.cpl.viewer.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.viewer.commandline.arguments.BooleanArgumentDefinition;
import org.fhcrc.cpl.viewer.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.viewer.feature.FeatureSet;
import org.fhcrc.cpl.viewer.quant.Q3;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;

import org.apache.log4j.*;

import java.awt.*;
import java.beans.PropertyChangeListener;
import java.io.*;
import java.util.*;
import java.util.List;
import javax.swing.*;

/**
 * If the "log" arg is specified, all ApplicationContext messages and all Log4J messages will all be
 * appended to a log file.  The user can specify a file for this purpose, or one can be written in the
 * temp dir.  If this logging is on and a CommandLineModule fails execution, the full log will be
 * appended to the failure log.  If a CLM fails execution and logging is NOT on, the user is prompted to
 * re-run the command with logging on.
 */
public class Application implements ApplicationContext.ApplicationContextProvider
{
    static Application _instance = new Application();

    private static Logger _log = Logger.getLogger(Application.class);

    // Fire the static converter initializers needed by the tab loader
    static
    {
        ConvertHelper.registerHelpers();
    }

    WorkbenchFrame _frame;
    protected boolean _isWorkbenchInitialized = false;
    AwtPropertyBag _properties;

    //dhmay adding 5/22/08 for logging
    protected static boolean enableLog = false;
    protected static File logFile = null;
    protected static PrintWriter logPrintWriter = null;
    protected static final String LOGFILE_DUMMY_CALLER = "//APPLICATION_LOGFILE_DUMMY_CALLER//";

    //dhmay adding 12/15/2005
    //Two different types of messages that can be displayed to the user
    protected static String MESSAGE_TYPE_INFO = "Information";
    protected static String MESSAGE_TYPE_ERROR = "Error";


    private static final String[] hardcodedCommandNames =
            {
                    "revision",
                    "q3",
                    "help",
                    "interactive"
            };

    private static final String[] hardcodedCommandDescriptions =
            {
                    "Revision information",
                    "Perform Q3 quantitation",
                    "List available commands and get usage information for individual commands",
                    "Show an interactive window for entering command arguments"
            };

    private static final String[] hardcodedCommandUsages =
            {
                    "",
                    "--labeledResidue=<residue> --massDiff=<mass> [--minPeptideProphet=<decimal>] [--maxFracDeltaMass=<decimal>] [--forceOutput] [--mimicXpress] [--noSentinels] [--stripAnalysisResults] --out=<output file> <input file>",
                    "[command] [--html]",
                    "command"
            };

    private Application()
    {
        _properties = new AwtPropertyBag(this);
    }


    public static Application getInstance()
    {
        return _instance;
    }


    public static void setMessage(String message)
    {
        if (null != _instance._frame)
            _instance._frame.setMessage(message);
    }




    private void start(String[] args)
    {
        JFrame splashFrame = WorkbenchFrame.ShowSplashScreen();
        ApplicationContext.setImpl(this);

        //set the appropriate locale and font
        Localizer.setInitialLocaleProperties();

        _frame = new WorkbenchFrame();
        setProperty("frame", _frame);
        if (args.length > 0)
        {
            OpenFile(args[0]);
        }
        else
        {
            final Action a = new OpenFileAction(null);
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    a.actionPerformed(null);
                }
            });
        }

        splashFrame.dispose();
        _frame.setVisible(true);
        _isWorkbenchInitialized = true;
    }

    /**
     * start from scratch with a new Workbench.  For language changing
     */
    public void redrawWorkbench()
    {
        if (!_isWorkbenchInitialized)
            return;
//If there are open files in the old workbenchframe, hopefully they'll get GC'ed when there are
//no more references.  Anyway, that's not likely -- language-switching is likely to be done before
//file loading
        _frame.disposeNoExit();

        _frame = new WorkbenchFrame();
        _frame.setVisible(true);
    }

    public void status(String message)
    {
        setMessage(message); 
    }

    public void errorMessage(String message, Throwable t)
    {
        String text = message == null ? "" : message;
        if (null != t)
        {
            StringWriter sw = new StringWriter();
            PrintWriter w = new PrintWriter(sw);
            t.printStackTrace(w);
            w.flush();
            if (text.length() > 0)
                text += "\n";
            text += sw.toString();
        }

        if (t instanceof RuntimeException)
            t.printStackTrace(System.err);
        showMessage(text,MESSAGE_TYPE_ERROR);
    }

/*
 * dhmay adding 12/15/2005
 * Just shows an informational message, no stack trace
 */
    public void infoMessage(String message)
    {
        String text = message == null ? "" : message;

        showMessage(text,MESSAGE_TYPE_INFO);
    }

/*
 * dhmay adding 12/15/2005
 * common code between errorMessage and infoMessage, not public
 */
    protected void showMessage(String text, String messageType)
    {

        if (EventQueue.isDispatchThread())
        {
            JOptionPane.showMessageDialog(getFrame(), text, messageType, JOptionPane.INFORMATION_MESSAGE);
        }
        else
        {
            final String finalText = text;
            final String finalMessageType = messageType;
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    JOptionPane.showMessageDialog(getFrame(), finalText, finalMessageType, JOptionPane.INFORMATION_MESSAGE);
                }
            });
        }

        if (enableLog)
        {
            logPrintWriter.println(text);
            logPrintWriter.flush();
        }
    }


    public JFrame getFrame()
    {
        return _instance._frame;
    }


    public void setProperty(String name, Object value)
    {
        _properties.put(name, value);
    }


    public Object getProperty(String name)
    {
        return _properties.get(name);
    }


    public void addPropertyChangeListener(PropertyChangeListener listener)
    {
        _properties.addPropertyChangeListener(listener);
    }


    public void addPropertyChangeListener(String name, PropertyChangeListener listener)
    {
        _properties.addPropertyChangeListener(name, listener);
    }

    public void OpenFile(String filename)
    {
        final File file = new File(filename);
        boolean okFile = true;
        if (!file.exists())
        {
            okFile = false;
            ApplicationContext.infoMessage(
                    TextProvider.getText("COULD_NOT_OPEN_FILE_FILENAME", filename));
        }
        else if (!file.canRead())
        {
            okFile = false;
            ApplicationContext.infoMessage(
                    TextProvider.getText("COULD_NOT_OPEN_FILE_FILENAME", filename));            
        }

        if (!okFile)
            return;

        Thread t = new Thread(new Runnable()
        {
            public void run()
            {
                Thread.yield();
                try
                {
                    ApplicationContext.setProperty(SharedProperties.MS_SCAN, null);
                    ApplicationContext.setProperty(SharedProperties.MS_RUN, null);

                    ApplicationContext.setProperty(SharedProperties.FEATURE_RANGES, null);
                    ApplicationContext.setProperty("featureSelector", null);
                    MSRun run = MSRun.load(file.getPath());
                    ApplicationContext.setProperty(SharedProperties.MS_RUN, run);
                    ApplicationContext.setProperty(SharedProperties.MS_SCAN, run.getScan(0));
                    ApplicationContext.setProperty(SharedProperties.SELECTED,  run);
                }
                catch (IOException x)
                {
                }
            }
        });
        t.setPriority(Thread.MIN_PRIORITY);
        t.start();
    }

    /**
     * Handle arguments and run a command
     * @param module
     * @param args
     */
    protected static void runCommand(CommandLineModule module, String[] args)
    {

        if (!processRawArguments(module, args))
        {
            ApplicationContext.infoMessage(module.getUsage());
            quit();
        }

        try
        {
            module.execute();
        }
        catch (Exception e)
        {
            if (e instanceof CommandLineModuleExecutionException &&
                    ((CommandLineModuleExecutionException) e).shouldShowStackTrace())
                ApplicationContext.errorMessage(e.getMessage(), e);
            else
                ApplicationContext.infoMessage(e.getMessage());
            String message =
                    TextProvider.getText("ERROR_RUNNING_COMMAND_COMMAND", module.getCommandName());

            message = message + "\n" + e.getMessage() + "\n";

            StringWriter sw = new StringWriter();
            PrintWriter w = new PrintWriter(sw);
            e.printStackTrace(w);
            w.flush();
            message += "\n";
            message += sw.toString();
            System.err.println(message);
            System.err.println(CommandLineModuleUtilities.createFailureReportAndPrompt(module, e));
        }
    }


    /**
     * Starts the application. Application can be started in either interactive mode or
     * specifically to find peptides in a file.
     * Interactive Mode: Application [mzxmlfile]
     * Find Peptides: Application --findPeptides [--out=outfile] mzxmlfile
     * If outfile is not specified, outfile is mzxmlfile.peptides.tsv
     *
     * @param args
     */
    public static void main(String[] args) throws Exception
    {
        try
        {
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
        }
        catch (Exception x)
        {
            //x.printStackTrace();
        }

        String revision = null;
        try
        {
            InputStream is = Application.class.getClassLoader().getResourceAsStream("revision.properties");
            Properties props = new Properties();
            props.load(is);

            revision = props.getProperty("SVNRevision");
            if (null == revision)
                revision = props.getProperty("revision");
            if ("".equals(revision))
                revision = null;
            ApplicationContext.setProperty("REVISION", revision);
        }
        catch (Exception x)
        {
        }

        //Special handling for the "debug" argument.  This enables global debug logging
        //for all classes under the package org.fhcrc.cpl.  This can be very verbose.
        String debugArg = null;
        boolean enableDebug = false;
        for (int i=0; i<args.length && !enableDebug; i++)
        {
            String thisArg = args[i];
            if (thisArg.startsWith("--"))
                thisArg = thisArg.substring(2);

            if ("debug".equalsIgnoreCase(thisArg) || thisArg.toLowerCase().startsWith("debug="))
            {
                debugArg = thisArg;
                //reshuffle arguments, get rid of debug arg
                String[] newArgs = new String[args.length-1];
                int newIndex = 0;
                for (int j=0; j<args.length; j++)
                {
                    if (j != i)
                       newArgs[newIndex++] = args[j];
                }
                args = newArgs;

                enableDebug = true;
                break;
            }
        }

        //Special handling for the "log" argument.  This enables global logging to a file
        //both ApplicationContext.infoMessage/setMessage/errorMessage messages and log4j messages
        for (int i=0; i<args.length && !enableLog; i++)
        {
            String thisArg = args[i];
            if (thisArg.startsWith("--"))
                thisArg = thisArg.substring(2);

            if (thisArg.toLowerCase().startsWith("log"))
            {
                if (thisArg.length() > 3 && !(thisArg.charAt(3) == '='))
                    continue;

                if (thisArg.length() == 4)
                    continue;

                if (thisArg.contains("="))
                {
                    logFile = new File(thisArg.substring(4));
                }
                else
                    logFile = TempFileManager.createTempFile("msInspect.log", LOGFILE_DUMMY_CALLER);

                try
                {
                    logPrintWriter = new PrintWriter(logFile);
                    logPrintWriter.println("msInspect log starting " + new Date());
                    logPrintWriter.flush();
                }
                catch (Exception e)
                {
                    throw new IllegalArgumentException("Can't write to log file " + logFile.getAbsolutePath());
                }

                FileAppender logFileAppender = new FileAppender();
                logFileAppender.setName("Logfile_appender");
                logFileAppender.setFile(logFile.getAbsolutePath(), true, false, 10);
                logFileAppender.setLayout(new SimpleLayout());

                Logger cplRootLogger = Logger.getLogger("org.fhcrc.cpl");
                cplRootLogger.setAdditivity(true);
                cplRootLogger.addAppender(logFileAppender);

                Logger.getLogger(ApplicationContext.class).addAppender(logFileAppender);

                enableLog = true;


                ApplicationContext.infoMessage("Logging this session to file " + logFile.getAbsolutePath());

                String[] newArgs = new String[args.length-1];
                int newIndex = 0;                
                for (int j=0; j<args.length; j++)
                {
                    if (j != i)
                       newArgs[newIndex++] = args[j];
                }
                args = newArgs;


            }
        }

        if (enableDebug)
        {
            if ("debug".equalsIgnoreCase(debugArg) || "debug=true".equalsIgnoreCase(debugArg))
            {
                ApplicationContext.setMessage("Enabling all debug logging");
                Logger cplRootLogger = Logger.getLogger("org.fhcrc.cpl");
                cplRootLogger.setAdditivity(true);
                cplRootLogger.setLevel(Level.DEBUG);
            }
            else
            {
                String[] debugClassNames = debugArg.substring("debug=".length()).split(",");
                Class[] debugClasses = new Class[debugClassNames.length];

                for (int i=0; i<debugClassNames.length; i++)
                {
                    String debugClassName = debugClassNames[i];
                    try
                    {
                        Class debugClass = Class.forName(debugClassName);
                        debugClasses[i] = debugClass;
                    }
                    catch (ClassNotFoundException e)
                    {
                        ApplicationContext.setMessage("Unable to find class with name " + debugClassName + ". Quitting");
                        return;
                    }
                }
                for (Class debugClass : debugClasses)
                {
                    Logger classLogger = Logger.getLogger(debugClass);
                    classLogger.addAppender(Logger.getRootLogger().getAppender("CONSOLE"));
                    classLogger.setAdditivity(false);
                    classLogger.setLevel(Level.DEBUG);
                }
                ApplicationContext.setMessage("Debug logging is enabled for the following classes:");
                for (Class debugClass : debugClasses)
                    ApplicationContext.setMessage(debugClass.getName());
            }


        }

        String command = args.length > 0 ? args[0] : "";
        command = command.toLowerCase();

        //runmacro is handled specially since it invokes the GUI
        if (command.startsWith("--"))
            command = command.substring(2);
        else
        {
            //No command identified, start GUI
            Application app = Application.getInstance();
            app.setProperty("REVISION", revision); // new ApplicationContext, need to set again
            app.start(args);
            return;
        }

        //If we get here, we've got a command starting with '--', so let's figure out what to do with it
        //First, check if the command is one recognized by one of the standard CommandLineModules.
        //Only if we don't find it there do we go to the hardcoded handlers in this file.
        try
        {
            CommandLineModule standardModule = CommandLineModuleDiscoverer.getStandardCommandLineModule(command);
            runCommand(standardModule, args);
            closeLog();
            return;
        }
        catch (FileNotFoundException e)
        {
            //do nothing -- execution beyond this point assumes we didn't find it
        }

        //No standard module for this.  Try the hardcoded commands
        //TODO: get rid of the hardcoded handlers in favor of CommandLineModule implementations
        if (command.equalsIgnoreCase("ms2Correct"))
        {
            ms2Correct(args);
            quit();
        }
        else if (command.equalsIgnoreCase("q3"))
        {
            Q3.run(args);
            quit();
        }
        else if (command.equalsIgnoreCase("version"))
        {
            ApplicationContext.infoMessage(revision);
            quit();
        }
        else if (command.equalsIgnoreCase("interactive"))
        {
            //a mode that gets rid of the commandline-related buttons and
            //the user notification that the command is complete
            boolean guiOnly = false;
            String commandForInteract = "";

            int beginModuleArgIndex = 2;

            if (args.length > 1)
            {
                commandForInteract = args[1].toLowerCase();

                if (commandForInteract.equalsIgnoreCase("--guionly"))
                {
                    guiOnly = true;
                    if (args.length > 1)
                    {
                        commandForInteract = args[2].toLowerCase();
                        beginModuleArgIndex = 3;
                    }
                }
            }

            if (commandForInteract.startsWith("--"))
                commandForInteract = commandForInteract.substring(2);

            CommandLineModule moduleForInteract = null;
            try
            {
                //this will fail if commandForInteract is null
                moduleForInteract = CommandLineModuleDiscoverer.getCommandLineModule(commandForInteract);
            }
            catch (Exception e)
            {
                if (commandForInteract != null && commandForInteract.length()>0)
                      JOptionPane.showMessageDialog(ApplicationContext.getFrame(),
                                "Unknown command " + commandForInteract,
                                "Information", JOptionPane.INFORMATION_MESSAGE);

                //we'll get a FileNotFoundException if the command isn't found,
                //and an ArrayIndexOutOfBoundsException if the user didn't specify one.
                //Either way, do this
                SpecifyModuleArgumentsFrame.ChooseCommandDialog chooseCommandDialog =
                        new SpecifyModuleArgumentsFrame.ChooseCommandDialog();

                moduleForInteract = chooseCommandDialog.chooseCommand();
            }
            if (moduleForInteract == null)
                quit();

            Map<String,String> argMap = null;
            if (args.length > beginModuleArgIndex)
            {
                //leave a dummy arg in there representing the command
                String[] moduleArguments = new String[args.length - beginModuleArgIndex + 1];
                System.arraycopy(args, beginModuleArgIndex, moduleArguments, 1, moduleArguments.length-1);

                argMap = parseRawArguments(moduleForInteract, moduleArguments);
            }


            SpecifyModuleArgumentsFrame interactFrame =
                    new SpecifyModuleArgumentsFrame(moduleForInteract, guiOnly, argMap);
            boolean shouldExecute = interactFrame.collectArguments();

            if (shouldExecute)
            {

                try
                {
                    moduleForInteract.execute();
                    if (!guiOnly)
                    {
                        JOptionPane.showMessageDialog(ApplicationContext.getFrame(),
                                TextProvider.getText("COMMAND_COMPLETE",moduleForInteract.getCommandName()) + "\n" +
                                        TextProvider.getText("CHECK_COMMAND_WINDOW_FOR_DETAILS"),
                                "Information", JOptionPane.INFORMATION_MESSAGE);
                    }
                }
                catch (Exception e)
                {
                    
                    String message =
                            TextProvider.getText("ERROR_RUNNING_COMMAND_COMMAND", moduleForInteract.getCommandName());

                    message = message + "\n" + e.getMessage() + "\n";

                    StringWriter sw = new StringWriter();
                    PrintWriter w = new PrintWriter(sw);
                    e.printStackTrace(w);
                    w.flush();
                    message += "\n";
                    message += sw.toString();
                    System.err.println(message);
                    message += CommandLineModuleUtilities.createFailureReportAndPrompt(moduleForInteract, e);
                    JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
                }
            }

            closeLog();
            return;
        }
        else if (command.equalsIgnoreCase("help") || command.equals("?"))
        {
            boolean isHtml = false;
            if (args.length > 1)
            {
                for (int i=1; i<args.length; i++)
                {
                    if (args[i].equalsIgnoreCase("--html"))
                        isHtml = true;
                }
            }

            if (args.length > (isHtml ? 2 : 1))
            {
                String commandForHelp = args[1].toLowerCase();
                if (isHtml && commandForHelp.equalsIgnoreCase("--html"))
                    commandForHelp = args[2].toLowerCase();

                try
                {
                    CommandLineModule module =
                            CommandLineModuleDiscoverer.getCommandLineModule(commandForHelp);
                    if (isHtml)
                    {
                        String dummyCaller = "dummy_help_caller";
                        File tempHelpFile =
                                TempFileManager.createTempFile("help_" + module.getCommandName(), dummyCaller);
                        PrintWriter outPW = new PrintWriter(tempHelpFile);
                        CLMUserManualGenerator.generateCommandManualEntry(module, outPW);
                        outPW.flush();
                        outPW.close();
                        HtmlViewerPanel.showFileInDialog(tempHelpFile, "Help for command " + module.getCommandName());
                    }
                    else
                        showHelp(module);
                }
                catch (FileNotFoundException e)
                {
                    boolean foundIt = false;
                    for (int i = 0; i < hardcodedCommandNames.length; i++)
                    {
                        if (commandForHelp.equalsIgnoreCase(hardcodedCommandNames[i]))
                        {
                            ApplicationContext.infoMessage("\nUsage:");
                            ApplicationContext.infoMessage("--" + hardcodedCommandNames[i] +
                                    " " + hardcodedCommandUsages[i]);
                            ApplicationContext.infoMessage("Details:");
                            ApplicationContext.infoMessage(hardcodedCommandDescriptions[i]);
                            ApplicationContext.infoMessage("\n");
                            foundIt = true;
                            break;
                        }
                    }
                    if (!foundIt)
                        showUsage();
                }
            }
            else
            {
                if (isHtml)
                {
                    String dummyCaller = "dummy_help_caller";
                    File tempHelpFile =
                            TempFileManager.createTempFile("help", dummyCaller);
                    PrintWriter outPW = new PrintWriter(tempHelpFile);
                    CLMUserManualGenerator.generateFullManual(outPW);
                    outPW.flush();
                    outPW.close();
                    HtmlViewerPanel.showFileInDialog(tempHelpFile, "Commandline Help");
                }
                else
                {
                    showUsage();
                    quit();
                }
            }
        }
        else
        {
            //Didn't find it in the standard modules, look for custom
            try
            {
                CommandLineModule customModule =
                        CommandLineModuleDiscoverer.getCommandLineModule(command);
                runCommand(customModule, args);
                return;
            }
            catch (FileNotFoundException e)
            {
                //completely unknown command
                ApplicationContext.infoMessage("Unknown command " + command);
                showUsage();
                quit();
            }
        }
    }


    

    /**
     * Show help for a given module
     * @param module
     */
    protected static void showHelp(CommandLineModule module)
    {
        ApplicationContext.infoMessage("\nUsage:");
        ApplicationContext.infoMessage(module.getUsage());
        ApplicationContext.infoMessage("Details:");
        ApplicationContext.infoMessage(module.getFullHelp());
    }

    /**
     * Does this set of argument definitions allow the unnamed (single) argument?
     * @return
     */
    protected static boolean containsUnnamedArgumentDefinition(CommandLineModule module)
    {
        return containsArgumentDefinition(module,
                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT);
    }

    protected static boolean containsUnnamedArgumentSeriesDefinition(CommandLineModule module)
    {
        return containsArgumentDefinition(module,
                CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT);
    }

    /**
     * Does this set of argument definitions allow a particular argument?
     * @return
     */
    protected static boolean containsArgumentDefinition(CommandLineModule module,
                                              String argument)
    {
        for (CommandLineArgumentDefinition def : module.getArgumentDefinitions())
        {
            if (argument.equalsIgnoreCase(def.getArgumentName()))
                return true;
        }
        return false;
    }

    /**
     * For CommandLineModule commands, parse the argument strings, digest the relevant arguments
     * according to the definitions provided by the module.  If any undefined arguments are
     * specified, die
     *
     * @param args
     * @return
     */
    public static boolean processRawArguments(CommandLineModule module, String[] args)
    {
        try
        {
            Map<String,String> argNameValueMap = parseRawArguments(module, args);
            module.digestArguments(argNameValueMap);
        }
        catch (ArgumentValidationException e)
        {
            ApplicationContext.infoMessage("Failure while parsing arguments:");
            if (e.shouldShowStackTrace())
                ApplicationContext.errorMessage(e.getMessage(), e);
            else
                ApplicationContext.infoMessage(e.getMessage());
            return false;
        }
        catch (IllegalArgumentException iae)
        {
            ApplicationContext.infoMessage(iae.getMessage());
            return false;
        }

        return true;
    }

    protected static Map<String,String> parseRawArguments(CommandLineModule module, String[] args)
            throws IllegalArgumentException
    {
        HashMap<String, String> argNameValueMap = new HashMap<String, String>();
        boolean alreadyProcessedNakedArgument = false;
        List<String> unnamedArgumentSeriesList = new ArrayList<String>();

        //Start with the second argument, because the first is the command
        for (int i = 1; i < args.length; i++)
        {
            _log.debug("Processing argument: " + args[i]);
            boolean ignoreThisArg = false;
            String[] param = args[i].split("=");

            //if there's something really funky with the argument, punt
            if (param.length < 1 || param.length > 2)
            {
                throw new IllegalArgumentException("Unable to parse argument " + args[i]);
            }
            if (param.length == 1)
            {
                //an argument with no = sign.  First check for any BooleanArgumentDefinitions
                //that match the name
                if (param[0].startsWith("--"))
                {
                    CommandLineArgumentDefinition argDef =
                            module.getArgumentDefinition(param[0].substring(2).toLowerCase());
                    if (argDef != null)
                    {
                        //if it's a boolean argument definition, set value to true.  Otherwise, bad
                        if (argDef instanceof BooleanArgumentDefinition)
                        {
                            param = new String[2];
                            param[0] = argDef.getArgumentName();
                            param[1] = "true";
                        }
                        else
                        {
                            //make translatable?
                            ApplicationContext.infoMessage("No argument value specified for argument " + argDef.getArgumentName());
                        }
                    }
                }
                // Check for short args (e.g. "-m0.75")
                else if (param[0].startsWith("-"))
                {
                    String paramName = param[0].substring(1,2);
                    String paramVal = param[0].substring(2);

                    // Command-line implementation explicitly ignores case. This can be deadly for
                    // short args, so currently disallow all upper-case args.
                    if (!paramName.equalsIgnoreCase(paramName))
                    {
                        ApplicationContext.infoMessage("Upper case short-form arguments are disallowed; found '" +
                                paramName + "'");
                    }

                    if ("".equals(paramVal))
                        paramVal = checkBooleanDefault(module, paramName);

                    if (null != paramVal)
                    {
                        param = new String[2];
                        param[0] = paramName;
                        param[1] = paramVal;
                    }
                }
                //Not an unnamed BooleanArgumentDefinition... try the naked argument
                else
                {
                    //if we've got an argument with no = sign, then check to see if this
                    //command allows that.  If so, handle it that way.
                    if (containsUnnamedArgumentDefinition(module))
                    {

                        //two naked arguments, die
                        if (alreadyProcessedNakedArgument)
                        {
                            throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                                    args[i]));
                        }
                        param = new String[2];
                        param[0] = CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_ARGUMENT;
                        param[1] = args[i];
                        alreadyProcessedNakedArgument = true;
                    }
                    //ok, maybe this module allows a series of unnamed arguments....
                    else if (containsUnnamedArgumentSeriesDefinition(module))
                    {
                        unnamedArgumentSeriesList.add(args[i]);
                        ignoreThisArg = true;
                    }
                    else
                    {
                        throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                                args[i]));
                    }
                }
            }
            if (ignoreThisArg)
                continue;

            String paramName = param[0];
            if (paramName.startsWith("--"))
                paramName = paramName.substring(2);
            String paramVal = null;
            if (param.length > 1)
                paramVal = param[1];

            _log.debug("Got arg " + paramName + " = '" + paramVal + "'");

            if (containsArgumentDefinition(module, paramName))
            {
                if (argNameValueMap.containsKey(paramName.toLowerCase()))
                {
                    throw new IllegalArgumentException("Argument " + paramName +
                            " is specified more than once.  Quitting.\n");
                }
                argNameValueMap.put(paramName.toLowerCase(), paramVal);
            }
            else
            {
                throw new IllegalArgumentException(TextProvider.getText("UNKNOWN_PARAMETER_COLON_PARAM",
                        paramName));
            }
        } // End of arg loop

        if (containsUnnamedArgumentSeriesDefinition(module) &&
                !unnamedArgumentSeriesList.isEmpty())
        {
            StringBuffer combinedUnnamedArguments = new StringBuffer();
            boolean firstArg=true;
            for (String unnamedArgumentValue : unnamedArgumentSeriesList)
            {
                //use '*ARGSEP*' as arg separator
                if (!firstArg)
                    combinedUnnamedArguments.append(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                combinedUnnamedArguments.append(unnamedArgumentValue);
                firstArg = false;
            }
            argNameValueMap.put(CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT,
                    combinedUnnamedArguments.toString());
        }

        return argNameValueMap;
    }

    private static String checkBooleanDefault(CommandLineModule module, String paramName)
    {
        CommandLineArgumentDefinition argDef =
            module.getArgumentDefinition(paramName);

        // Argument not defined
        if (argDef == null)
            return null;

        // Not a boolean argument
        if (!(argDef instanceof BooleanArgumentDefinition))
        {
            ApplicationContext.infoMessage("No argument value specified for argument " + argDef.getArgumentName());
            return null;
        }

        return "true";
    }

    public static void closeLog()
    {
        if (enableLog)
        {
            ApplicationContext.setMessage("Session log written to " + logFile.getAbsolutePath());
            logPrintWriter.close();
        }
    }

    public static void quit()
    {
        closeLog();
        System.exit(0);
    }


    public static void quit(String err)
    {
        System.err.println(err);
        showUsage();
        System.exit(1);
    }


    private static final String helpStr = "USAGE:\n" + "java -jar -Xmx=memSize viewerApp.jar [mzxmlfile]\n" +
            "java -jar -Xmx=memSize viewerApp.jar <command> [options]\n" +
            "    memSize should be 512M minimum for interactive use or peptide finding. 1G is better.\n\n"    +
            "COMMANDS:\n";
    /**
     "--revision\n" +
     "--ms2Correct --features=ms1FeatureFile [--out=outfilename] [filterOptions]  runFile" + "\n" +
     "--optimize [--out=outfilename] [--scanWindows=int,int,...] [--massWindows=float,float,float] [filterOptions]  featureFiles..." + "\n" +
     "--deconvolute --quant [--out=outfilename] [--scanWindow=int] [--massWindow=float] [--lightTagWeight=float] [--heavyTagWeight=float] [--maxLabelCount=int] [--labeledResidue=char] [--intensityType={total,max,recalculated}] [--msFile=filepath] --deltaTime=[float] --deltaMass=[float][da|ppm] featureFile"  + "\n" +
     "    --deconvolute and --quant may be used together or separately\n" +
     "\n" +
     "[filterOptions] = \n" +
     "[--minMz=float] [--maxMz=float] [--minMass=float] [--maxMass=float]\n" +
     "[--minPeaks=int] [--minCharge=int] [--maxCharge=int] [--maxKL=float]\n" +
     "[--minIntensity=float] [--minTotalIntensity=float]\n" +
     "[--minTime=float] [--maxTime=float]\n" +
     "[--scanFirst=int] [scanLast=int] [--minScans=int]\n" +
     "scanFirst and scanLast filter the scan column. minScans filters the scanCount column. Other filter parameters filter the column matching their name.";
     */


    public static void showUsage()
    {
        ApplicationContext.infoMessage(helpStr);

        Map<String, String> commandHelpTextMap = new HashMap<String, String>();

        for (int i=0; i<hardcodedCommandNames.length; i++)
        {
            commandHelpTextMap.put(hardcodedCommandNames[i],
                    hardcodedCommandDescriptions[i]);
//            System.out.println("--" + hardcodedCommandNames[i] + ":\t" +
//                    hardcodedCommandDescriptions[i]);
        }

        for (CommandLineModule module : CommandLineModuleDiscoverer.findAllCommandLineModules().values())
        {
            commandHelpTextMap.put(module.getCommandName(),
                                   module.getShortDescription());

//            System.out.println("--" + module.getCommandName() + ":\t" +
//                    module.getShortDescription());
        }

        String[] commandArray =
                commandHelpTextMap.keySet().toArray(new String[commandHelpTextMap.size()]);
        Arrays.sort(commandArray);

        for (String command : commandArray)
            ApplicationContext.infoMessage("--" + command + ":\t" + commandHelpTextMap.get(command));

        ApplicationContext.infoMessage(
                "\nFor more information on individual commands, use the '--help' command followed by the command name\n");
    }

    /*
     * dhmay adding to centralize PrintWriter creation
     * Utility method: given an output filename, checks if null.  If not, returns a PrintWriter 
     * on that file. If so, returns a PrintWriter on System.err
     *
     * @param outFileName the output filename specified on the command line
     * @return an appropriate printwriter
     */
    protected static PrintWriter getPrintWriter(String outFileName) throws java.io.FileNotFoundException
    {
        PrintWriter pw = null;
        if (null != outFileName)
        {
            try
            {
                pw = new PrintWriter(new FileOutputStream(outFileName));
            }
            catch (java.io.FileNotFoundException e)
            {
                System.err.println("Error creating PrintWriter from file " + outFileName + ", file not found");
                throw e;
            }
        }
        else
            pw = new PrintWriter(System.out);

        return pw;
    }









    //Everything after this point should be removed or moved somewhere else








    /**
     * TODO: we should remove this, probably
     * @param args
     * @throws Exception
     */
    public static void ms2Correct(String[] args) throws Exception
    {
        String outFileName = null;
        File runFile = null;
        String featureFileName = null;
        File featureFile = null;
        FeatureSet.FeatureSelector sel = new FeatureSet.FeatureSelector();
        sel.setMinCharge(1);

        for (int i = 1; i < args.length; i++)
        {
            if (args[i].startsWith("--"))
            {
                String[] param = args[i].split("=");
                if (param.length != 2)
                    quit("Unknown param: " + args[i]);

                String paramName = param[0];
                String paramVal = param[1];

                if ("--out".equals(paramName))
                    outFileName = paramVal;
                else if ("--features".equals(paramName))
                    featureFileName = paramVal;
                else if (!sel.setFilterParam(paramName, paramVal))
                {
                    quit("Unknown parameter: " + paramName);
                }
            }
            else
            {
                if (null != runFile)
                    quit("Can only correct one run at a time. Found second run: " + args[i]);

                runFile = new File(args[i]);
                if (!runFile.exists())
                    quit("Couldn't find file: " + args[i]);
            }
        }

        if (null == runFile)
            quit("Must specify run to correct");

        if (null == featureFileName)
            quit("Must specify ms1 feature file name using --features=fileName parameter");

        featureFile = new File(featureFileName);
        if (!featureFile.exists())
            quit("Could not find feature file named " + featureFile.getPath());

        if (null == outFileName)
            outFileName = runFile.getPath() + ".ms2Correct.tsv";

        MS2Correct.correct(runFile, featureFile, sel, new File(outFileName));
    }

    public static boolean isLogEnabled()
    {
        return enableLog;
    }

    public static File getLogFile()
    {
        return logFile;
    }

}