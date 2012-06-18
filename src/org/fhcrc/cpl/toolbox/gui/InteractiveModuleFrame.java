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

import org.fhcrc.cpl.toolbox.commandline.arguments.*;
import org.fhcrc.cpl.toolbox.commandline.CLMUserManualGenerator;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.Map;
import java.util.HashMap;
import java.util.prefs.Preferences;
import java.util.prefs.BackingStoreException;
import java.io.*;

/**
 * JFrame class for specifying command-line arguments for a given module
 *
 * Override the method createArgumentFieldPanel() in a subclass to add handling for more arg types
 *
 * TODO: use swiXML?  Layout is variable, but some bits are constant, like the buttons
 */
public class InteractiveModuleFrame extends JDialog
{
    static Logger _log = Logger.getLogger(InteractiveModuleFrame.class);

    //module to invoke
    protected CommandLineModule module;

    //map from arguments to the components that will hold their values
    protected Map<CommandLineArgumentDefinition,JComponent> argComponentMap;

    //for storing previously specified values
    protected Preferences prefs = Preferences.userNodeForPackage(InteractiveModuleFrame.class);

    protected static final int MAX_FIELDPANE_HEIGHT = 600;
    protected static final int ARG_HELP_PANEL_HEIGHT = 90;

    public static final int DEFAULT_MAX_ARG_LABEL_LENGTH = 50;
    protected int maxArgLabelLength = DEFAULT_MAX_ARG_LABEL_LENGTH;

    //help with a specific argument
//    protected JDialog argHelpDialog;
    protected JTextArea argHelpTextArea;
//    protected JPanel argHelpPanel;
    protected JScrollPane argHelpScrollPane;

    //dialog box providing full help for the command
    protected JDialog helpDialog;


    protected Map<String, String> moduleArgMap;

    protected CLMUserManualGenerator userManualGenerator = new CLMUserManualGenerator();

    protected ListenerHelper helper = new ListenerHelper(this);

    protected JPanel buttonPanel = null;
    protected GridBagConstraints buttonGBC = null;

    //arbitrary
    public int width = 650;

    //will be overridden
    public int height = 300;
    public int fieldViewportHeight = 300;
    public int fieldPaneHeight = 300;
    public int fieldPanelWidth = width;


    //are we done with the frame?
    public boolean done = false;
    //has the user specified arguments?
    public boolean argsSpecified = false;

    //to hang listeners off of, to notify things that we're done.  not displayed
    protected JButton fakeButton;

    protected boolean shouldStoreAlgValues = true;


    public InteractiveModuleFrame(CommandLineModule module, Map<String, String> moduleArgMap)
    {
        this(module, moduleArgMap, DEFAULT_MAX_ARG_LABEL_LENGTH);
    }

    /**
     *
     * @param module
     * @param moduleArgMap initial values for specified args
     */
    public InteractiveModuleFrame(CommandLineModule module, Map<String, String> moduleArgMap, int maxArgLabelLength)
    {
        super();
        this.maxArgLabelLength = maxArgLabelLength;
        setTitle(TextProvider.getText("ARGUMENTS_FOR_COMMAND_COMMAND",module.getCommandName()));
        this.setModalityType(ModalityType.APPLICATION_MODAL);
        fakeButton = new JButton("fake");

        this.module = module;
        this.moduleArgMap = moduleArgMap;

        //buttons
        JButton buttonExecute = new JButton(TextProvider.getText("EXECUTE"));
        JButton buttonCancel = new JButton(TextProvider.getText("CANCEL"));
        JButton buttonShowHelp = new JButton(TextProvider.getText("HELP"));
        helper.addListener(buttonExecute, "buttonExecute_actionPerformed");
        helper.addListener(buttonCancel, "buttonCancel_actionPerformed");
        helper.addListener(buttonShowHelp, "buttonShowHelp_actionPerformed");
        //default action is execute
        getRootPane().setDefaultButton(buttonExecute);

        //set up the panel
        JPanel contentPanel = new JPanel();
        GridBagConstraints contentPanelGBC = new GridBagConstraints();
        contentPanelGBC.gridwidth = GridBagConstraints.REMAINDER;
        contentPanel.setLayout(new GridBagLayout());
        add(contentPanel);

        //add all argument fields
        addFieldsForArguments(contentPanel, helper);

        //add buttons
        buttonPanel = new JPanel();
        GridBagConstraints buttonPanelGBC = new GridBagConstraints();
        buttonPanelGBC.gridwidth = GridBagConstraints.REMAINDER;
        contentPanel.add(buttonPanel, buttonPanelGBC);

        buttonGBC = new GridBagConstraints();
        buttonGBC.insets = new Insets(10, 0, 10, 0);
        buttonPanel.add(buttonExecute, buttonGBC);

        GridBagConstraints secondToLastButtonGBC = new GridBagConstraints();
        secondToLastButtonGBC.insets = new Insets(10, 0, 10, 0);
        secondToLastButtonGBC.gridwidth = GridBagConstraints.RELATIVE;
        buttonPanel.add(buttonCancel, secondToLastButtonGBC);


        GridBagConstraints lastButtonGBC = new GridBagConstraints();
        lastButtonGBC.insets = new Insets(10, 0, 10, 0);
        lastButtonGBC.gridwidth = GridBagConstraints.REMAINDER;
        buttonPanel.add(buttonShowHelp, lastButtonGBC);

        //add help
        argHelpTextArea = new JTextArea();
        argHelpTextArea.setEditable(false);
        argHelpTextArea.setOpaque(false);
        argHelpTextArea.setBounds(0,0,500,200);
        argHelpTextArea.setLineWrap(true);
        argHelpTextArea.setWrapStyleWord(true);

        argHelpScrollPane = new JScrollPane();
        argHelpScrollPane.setBorder(BorderFactory.createTitledBorder("Argument Help"));
        argHelpScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        argHelpScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        argHelpScrollPane.setPreferredSize(new Dimension(width, ARG_HELP_PANEL_HEIGHT));
        argHelpScrollPane.setViewportView(argHelpTextArea);
        contentPanel.add(argHelpScrollPane, buttonPanelGBC);

        //20 * the number of text fields, plus the height of the button area,
        //plus some padding
        height = fieldPaneHeight + 100 + 70 + 15;

        setPreferredSize(new Dimension(width+30, height));
        setSize(new Dimension(width+30, height));
        setMinimumSize(new Dimension(width+30, height));

        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int centerH = screenSize.width / 2;
        int centerV = screenSize.height / 2;
        this.setLocation(centerH - width / 2, centerV - height / 2);

        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    }

    public void addDoneListener(ActionListener listener)
    {
        fakeButton.addActionListener(listener);
    }

    /**
     * Wait for argument specification
     * @return
     */
    public boolean collectArguments()
    {
        setVisible(true);
        return argsSpecified;
    }

    /**
     * Notify any listeners that we're done
     * @param event
     */
    protected void notifyDone(ActionEvent event)
    {
        ActionListener[] fakeButtonListeners = fakeButton.getActionListeners();

        if (fakeButtonListeners != null)
        {
            for (ActionListener listener : fakeButtonListeners)
                listener.actionPerformed(event);
        }
    }

    /**
     * Add fields representing each of the module's arguments.  Each arg handled separately by another method
     * TODO: move field creation into the argument classes?  That would be more modular. Also more work
     * @param contentPanel
     * @param helper
     */
    protected void addFieldsForArguments(JPanel contentPanel, ListenerHelper helper)
    {
        //set up scrollpane UI
        boolean firstArg = true;
        JScrollPane fieldScrollPane = new JScrollPane();
        fieldScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        fieldScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);          
        GridBagConstraints fieldPaneGBC = new GridBagConstraints();
        fieldPaneGBC.gridwidth=GridBagConstraints.REMAINDER;

        JPanel panelInViewport = new JPanel();
        panelInViewport.setLayout(new GridBagLayout());
        fieldScrollPane.setViewportView(panelInViewport);

        JPanel allFieldsPanel = new JPanel();
        allFieldsPanel.setLayout(new GridBagLayout());
        GridBagConstraints allFieldsPanelGBC = new GridBagConstraints();
        allFieldsPanelGBC.gridwidth=GridBagConstraints.REMAINDER;

        argComponentMap = new HashMap<CommandLineArgumentDefinition,JComponent>();

        //For each argument, create the GUI components that will capture the value.
        //Basic args first
        for (CommandLineArgumentDefinition argDef :
                module.getBasicArgumentDefinitions())
        {
            addComponentsForArg(argDef, firstArg, allFieldsPanel);
            firstArg = false;
        }

        //add advanced args, after a separator
        int extraAdvancedArgsHeight = 0; //extra height, not from the args themselves, but from the section title
        CommandLineArgumentDefinition[] advancedArgs = module.getAdvancedArgumentDefinitions();
        if (advancedArgs != null && advancedArgs.length > 0)
        {
            extraAdvancedArgsHeight = 25;
            JSeparator visibleSeparator = new JSeparator();
            visibleSeparator.setPreferredSize(new Dimension(fieldPanelWidth-10, 1));
            allFieldsPanel.add(visibleSeparator, allFieldsPanelGBC);
            JLabel advancedArgsLabel = new JLabel("Advanced Arguments");
            advancedArgsLabel.setToolTipText("The default values of these arguments are appropriate for " +
                    "most use cases.  Do not alter these values unless you know what you're doing.");
            allFieldsPanel.add(advancedArgsLabel, allFieldsPanelGBC);
            JSeparator visibleSeparator2 = new JSeparator();
            visibleSeparator2.setPreferredSize(new Dimension(fieldPanelWidth-10, 1));
            allFieldsPanel.add(visibleSeparator2, allFieldsPanelGBC);

            for (CommandLineArgumentDefinition argDef : advancedArgs)
                addComponentsForArg(argDef, false, allFieldsPanel);
        }

        //20 * the number of text fields, plus the height of the button area, plus whatever from the advanced area,
        //plus some padding
        fieldViewportHeight =  41 * argComponentMap.size() + 25 + extraAdvancedArgsHeight;
        fieldPaneHeight = Math.min(fieldViewportHeight, MAX_FIELDPANE_HEIGHT);
        fieldPaneHeight = Math.max(300, fieldPaneHeight);

        allFieldsPanel.setMinimumSize(new Dimension(fieldPanelWidth, fieldViewportHeight));
        allFieldsPanel.setPreferredSize(new Dimension(fieldPanelWidth, fieldViewportHeight));

        fieldScrollPane.setMinimumSize(new Dimension(fieldPanelWidth, fieldPaneHeight));
        fieldScrollPane.setPreferredSize(new Dimension(fieldPanelWidth, fieldPaneHeight));

        contentPanel.add(fieldScrollPane, fieldPaneGBC);
        panelInViewport.add(allFieldsPanel, allFieldsPanelGBC);
    }

    /**
     * Add the components that represent a particular argument.  This will include labels.
     * Inidividual handling for each arg type is done in yet another method, createArgumentFieldPanel
     * @param argDef
     * @param firstArg
     * @param allFieldsPanel
     */
    protected void addComponentsForArg(CommandLineArgumentDefinition argDef, boolean firstArg,
                                       JPanel allFieldsPanel)
    {
        String labelText = argDef.getArgumentDisplayName();
        if (CommandLineModuleUtilities.isUnnamedSeriesArgument(argDef) ||
                CommandLineModuleUtilities.isUnnamedArgument(argDef))
            labelText = argDef.getHelpText();
        if (labelText.length() > maxArgLabelLength)
            labelText = labelText.substring(0, maxArgLabelLength-3) + "...";

        JLabel argLabel = new JLabel(labelText);
        String toolTipText = argDef.getHelpText();
        if (toolTipText == null)
            toolTipText = "";
        if (toolTipText.length() > 50)
            toolTipText = toolTipText.substring(0, 46) + "....";
        argLabel.setToolTipText(toolTipText);

        JComponent argComponent;

        String fieldValue = null;

        //if module argument map is provided, don't use prefs for defaulting
        if (moduleArgMap != null)
        {
            if (moduleArgMap.containsKey(argDef.getArgumentName()))
            {
                fieldValue = moduleArgMap.get(argDef.getArgumentName());
            }
        }
        else
        {
            fieldValue = prefs.get(module.getCommandName() + ":" + argDef.getArgumentName(), null);
        }

        if (fieldValue == null || fieldValue.length() == 0)
        {
            if (argDef.hasDefaultValue())
            {
                //TODO: this may bomb on more complicated custom datatypes
                fieldValue = argDef.getDefaultValueAsString();
            }
        }

        JPanel fieldPanel = new JPanel();
//System.err.println("Handling " + argDef.getArgumentName());
        //treat unnamed series parameters differently
        if (CommandLineModuleUtilities.isUnnamedSeriesArgument(argDef))
            argComponent = argDef.addComponentsForGUISeries(fieldPanel, this, fieldValue);
        else
            argComponent = argDef.addComponentsForGUI(fieldPanel, this, fieldValue);
        argComponentMap.put(argDef, argComponent);

        JPanel labelPanel = new JPanel();
        GridBagConstraints labelPanelGBC = new GridBagConstraints();
        labelPanelGBC.gridwidth=GridBagConstraints.RELATIVE;
        labelPanelGBC.anchor = GridBagConstraints.LINE_END;
        labelPanelGBC.insets = new Insets(0,0,0,0);

        GridBagConstraints labelGBC = new GridBagConstraints();
        labelGBC.anchor = GridBagConstraints.LINE_END;
        labelGBC.gridwidth=GridBagConstraints.RELATIVE;

        GridBagConstraints helpGBC = new GridBagConstraints();
        helpGBC.anchor = GridBagConstraints.LINE_START;

        labelGBC.insets = new Insets(0, 0, 0, 0);

        if (argDef.isRequired())
        {
            JLabel requiredLabel = new JLabel("*");
            requiredLabel.setForeground(Color.BLUE);
            requiredLabel.setToolTipText("This argument is required");
            requiredLabel.addMouseListener(new ArgRequiredListener());
            labelPanel.add(requiredLabel);
        }
        labelPanel.add(argLabel, labelGBC);
        //only add help link if there's help text to show
        if (argDef.getHelpText() != null && argDef.getHelpText().length() > 0)
        {
            JLabel helpLabel = new JLabel("?");
            helpLabel.setToolTipText("Click for argument help");
            helpLabel.setForeground(Color.BLUE);
            helpLabel.addMouseListener(new ArgHelpListener(argDef));

            labelPanel.add(helpLabel, helpGBC);
        }


        allFieldsPanel.add(labelPanel, labelPanelGBC);
        GridBagConstraints fieldGBC = new GridBagConstraints();
        fieldGBC.gridwidth=GridBagConstraints.REMAINDER;
        fieldGBC.anchor = GridBagConstraints.LINE_START;
        fieldGBC.insets = new Insets(0,0,0,0);
        allFieldsPanel.add(fieldPanel, fieldGBC);
    }


    /**
     * Tells the user that an arg is required
     */
    protected class ArgRequiredListener implements MouseListener
    {
        public void mouseClicked(MouseEvent e)
        {
            infoMessage("This argument is required");
        }

        public void mousePressed(MouseEvent e) {}
        public void mouseReleased(MouseEvent e) {}
        public void mouseEntered(MouseEvent e) {}
        public void mouseExited(MouseEvent e) {}
    }

    /**
     * Display arg help in a box
     */
    protected class ArgHelpListener implements MouseListener
    {
        protected CommandLineArgumentDefinition argDef;

        public ArgHelpListener(CommandLineArgumentDefinition argDef)
        {
            this.argDef = argDef;
        }


        public void mouseClicked(MouseEvent e)
        {
/*
//dhmay getting rid of separate argument help dialog, replacing with an always-there text area
            if (argHelpDialog != null && argHelpDialog.isVisible())
            {
                //help already initialized and visible
            }
            else
            {
                argHelpTextArea = new JTextArea();
                argHelpTextArea.setEditable(false);
                argHelpTextArea.setOpaque(false);
                argHelpTextArea.setBounds(0,0,500,200);
                argHelpTextArea.setLineWrap(true);
                argHelpTextArea.setWrapStyleWord(true);

                argHelpDialog = new JDialog();
                argHelpDialog.add(argHelpTextArea);
                argHelpDialog.setVisible(true);
                Point thisLocation = getLocation();
                argHelpDialog.setLocation((int) thisLocation.getX() + 10, (int) thisLocation.getY() + 5);
                argHelpDialog.setSize(600,150);
                argHelpDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                argHelpDialog.setModalityType(ModalityType.DOCUMENT_MODAL);
//                argHelpDialog.setAlwaysOnTop(true);
            }
            argHelpDialog.setTitle("Help for argument '" + argName + "'");

            argHelpTextArea.setVisible(true);
            argHelpDialog.setVisible(true);
*/
            argHelpTextArea.setText(argDef.getHelpText());
            argHelpScrollPane.setBorder(BorderFactory.createTitledBorder("Argument Help for '" +
                    argDef.getArgumentDisplayName() + "'"));
            argHelpScrollPane.getVerticalScrollBar().setValue(argHelpScrollPane.getVerticalScrollBar().getMinimum());
        }

        public void mousePressed(MouseEvent e) {}
        public void mouseReleased(MouseEvent e) {}
        public void mouseEntered(MouseEvent e) {}
        public void mouseExited(MouseEvent e) {}
    }

    /**
     * grab all argument field values
     * @return
     */
    protected Map<String,String> parseFieldArguments()
    {
        Map<String,String> argNameValueMap = new HashMap<String,String>();
        String PROTECTED_SPACE_SEP = "MYSEP_PROTECTEDSPACE_MYSEP";
        for (CommandLineArgumentDefinition argDef : argComponentMap.keySet())
        {
            JComponent argComponent = argComponentMap.get(argDef);

            String argValue = argDef.getValueFromGUIComponent(argComponent);

            if (argValue != null && argValue.length() > 0)
            {
                if (argDef.getArgumentName().equals(
                        CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
                {
                    _log.debug("Unnamed series.  Value before:\n" + argValue);

                    argValue = argValue.trim();

                    if (FileArgumentDefinition.class.isAssignableFrom(argDef.getClass()) ||
                            FileToReadListArgumentDefinition.class.isAssignableFrom(argDef.getClass()))
                    {
                        //First, look for " symbols.  Protect any unprotected spaces between them with a \
                        StringBuffer alteredArg = new StringBuffer();
                        boolean insideQuote = false;
                        boolean afterBackslash = false;
                        for (byte curByte : argValue.getBytes())
                        {
                            if (curByte == '\\')
                                afterBackslash = true;
                            else
                            {
                                if (curByte == '\"')
                                {
                                    insideQuote = !insideQuote;
                                }
                                else if (curByte == ' ' && !afterBackslash && insideQuote)
                                    alteredArg.append('\\');
                                afterBackslash = false;
                            }
                            alteredArg.append((char) curByte);
                        }
                        argValue = alteredArg.toString();
                        _log.debug("After in-quote space protection: " + argValue);

                        //dhmay adding protection for a backslash followed by a space, convention on Linux for escaping
                        //a space
                        argValue = argValue.replaceAll("\\\\ ", PROTECTED_SPACE_SEP);
                        _log.debug("After protect:\n" + argValue);
                        argValue = argValue.replaceAll(" ", CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                        argValue = argValue.replaceAll(PROTECTED_SPACE_SEP, " ");
                    }
                    else
                    {
                        argValue = argValue.replaceAll(" ", CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                    }
                    _log.debug("After:\n" + argValue);                    
                }

                argNameValueMap.put(argDef.getArgumentName(), argValue);
                if (shouldStoreAlgValues)
                {
                    prefs.put(module.getCommandName() + ":" + argDef.getArgumentName(), argValue);
                    try
                    {
                        prefs.flush();
                    }
                    catch (BackingStoreException e)
                    {
                        _log.debug("BackingStoreException saving prefs for " + module.getCommandName() +
                                ":" + argDef.getArgumentName(), e);
                    }
                }
            }
            else
            {
                prefs.remove(module.getCommandName() + ":" + argDef.getArgumentName());
            }
        }

        return argNameValueMap;
    }

    /**
     * Show a dialog box with help for this command
     * @param event
     */
    public void buttonShowHelp_actionPerformed(ActionEvent event)
    {
        File tempFile = null;
        try
        {
            tempFile = TempFileManager.createTempFile("help_" + module.getCommandName() + ".html", this);
            PrintWriter outPW = new PrintWriter(tempFile);
            userManualGenerator.generateCommandManualEntry(module, outPW);
            outPW.flush();
            HtmlViewerPanel.showFileInDialog(tempFile, "Manual for commmand '" + module.getCommandName() + "'");
        }
        catch (IOException e)
        {
            ApplicationContext.infoMessage("Failed to open browser, HTML is in " +
                    tempFile.getAbsolutePath());
        }        
    }


    /**
     * Try to execute command
     * @param event
     */
    public void buttonExecute_actionPerformed(ActionEvent event)
    {
        Map<String,String> argNameValueMap = parseFieldArguments();

        try
        {
            _log.debug("Executing action.  Arguments:");
            for (String argName : argNameValueMap.keySet())
            {
                _log.debug("\t" + argName + "=" + argNameValueMap.get(argName));
            }
            module.digestArguments(argNameValueMap);
        }
        catch (ArgumentValidationException e)
        {
            infoMessage(TextProvider.getText("FAILED_ARGUMENT_VALIDATION") + "\n" + 
                        TextProvider.getText("ERROR") + ": " + e.getMessage());
            return;
        }
        _log.debug("Done digesting arguments");

        argsSpecified = true;

        notifyDone(event);
        disposeAllComponents();
    }

    /**
     * Get rid of everything
     */
    protected void disposeAllComponents()
    {
//        if (argHelpDialog != null)
//        {
//            argHelpDialog.setVisible(false);
//            argHelpDialog.dispose();
//        }
        

        if (helpDialog != null)
        {
            helpDialog.setVisible(false);
            helpDialog.dispose();
        }

        this.setVisible(false);
        this.dispose();
    }

    protected void infoMessage(String message)
    {
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
    }

    protected void errorMessage(String message, Throwable t)
    {
        if (null != t)
        {
            message = message + "\n" + t.getMessage() + "\n";

            StringWriter sw = new StringWriter();
            PrintWriter w = new PrintWriter(sw);
            t.printStackTrace(w);
            w.flush();
            message += "\n";
            message += sw.toString();
        }
        JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information", JOptionPane.INFORMATION_MESSAGE);
    }

    /**
     * Cancel, close window
     * @param e
     */
    public void buttonCancel_actionPerformed(ActionEvent e)
    {
        disposeAllComponents();
        done=true;
        notifyDone(e);
    }

    public CLMUserManualGenerator getUserManualGenerator()
    {
        return userManualGenerator;
    }

    public void setUserManualGenerator(CLMUserManualGenerator userManualGenerator)
    {
        this.userManualGenerator = userManualGenerator;
    }

    public int getMaxArgLabelLength()
    {
        return maxArgLabelLength;
    }

    public void setMaxArgLabelLength(int maxArgLabelLength)
    {
        this.maxArgLabelLength = maxArgLabelLength;
    }

    public boolean isShouldStoreAlgValues()
    {
        return shouldStoreAlgValues;
    }

    public void setShouldStoreAlgValues(boolean shouldStoreAlgValues)
    {
        this.shouldStoreAlgValues = shouldStoreAlgValues;
    }
}

