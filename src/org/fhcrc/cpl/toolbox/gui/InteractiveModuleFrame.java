/*
 * Copyright (c) 2003-2008 Fred Hutchinson Cancer Research Center
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

import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentValidationException;
import org.fhcrc.cpl.toolbox.commandline.arguments.CommandLineArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.EnumeratedValuesArgumentDefinition;
import org.fhcrc.cpl.toolbox.commandline.arguments.ArgumentDefinitionFactory;
import org.fhcrc.cpl.toolbox.commandline.CLMUserManualGenerator;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TempFileManager;
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
public class InteractiveModuleFrame extends JFrame
{
    static Logger _log = Logger.getLogger(InteractiveModuleFrame.class);

    //module to invoke
    protected CommandLineModule module;

    //map from arguments to the components that will hold their values
    protected Map<CommandLineArgumentDefinition,JComponent> argComponentMap;

    //for storing previously specified values
    protected Preferences prefs = Preferences.userNodeForPackage(InteractiveModuleFrame.class);

    protected static final int MAX_FIELDPANE_HEIGHT = 600;

    //dialog box containing help with a specific argument
    protected JDialog argHelpDialog;
    protected JTextArea argHelpTextArea;

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





    /**
     *
     * @param module
     * @param moduleArgMap initial values for specified args
     */
    public InteractiveModuleFrame(CommandLineModule module, Map<String, String> moduleArgMap)
    {
        super(TextProvider.getText("ARGUMENTS_FOR_COMMAND_COMMAND",module.getCommandName()));

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


        //20 * the number of text fields, plus the height of the button area,
        //plus some padding
        height = fieldPaneHeight + 70 + 15;

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
        while (!done)
        {
            try
            {
                Thread.sleep(500);
            }
            catch (InterruptedException e)
            {

            }
        }

        disposeAllComponents();

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
            addComponentsForArg(argDef, firstArg, helper, allFieldsPanel);
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
                addComponentsForArg(argDef, false, helper, allFieldsPanel);
        }

        //20 * the number of text fields, plus the height of the button area, plus whatever from the advanced area,
        //plus some padding
        fieldViewportHeight =  33 * argComponentMap.size() + 25 + extraAdvancedArgsHeight;
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
     * @param helper
     * @param allFieldsPanel
     */
    protected void addComponentsForArg(CommandLineArgumentDefinition argDef, boolean firstArg,
                                       ListenerHelper helper, JPanel allFieldsPanel)
    {
        String labelText = argDef.getArgumentDisplayName();
        if (CommandLineModuleUtilities.isUnnamedSeriesArgument(argDef) ||
                CommandLineModuleUtilities.isUnnamedArgument(argDef))
            labelText = argDef.getHelpText();
        if (labelText.length() > 50)
            labelText = labelText.substring(0, 47) + "...";

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
                //TODO: this will probably bomb on some of the more complicated datatypes
                fieldValue = argDef.getDefaultValueAsString();
            }
        }

        JPanel fieldPanel = null;

        //treat unnamed series parameters differently
        if (CommandLineModuleUtilities.isUnnamedSeriesArgument(argDef))
        {
            //if a series of files to read (most likely case), handle super-specially.  Else, a big text field
            JTextField argTextField = new JTextField();
            argTextField.setPreferredSize(new Dimension(220, 20));
            argTextField.setMinimumSize(new Dimension(220, 20));

            boolean fieldHasValue = (fieldValue != null && fieldValue.length() > 0);

            if (fieldHasValue)
            {
                fieldValue = fieldValue.replaceAll(CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR, " ");
                argTextField.setText(fieldValue);
            }
            argComponent = argTextField;


            fieldPanel = new JPanel();


            GridBagConstraints textFieldGBC = new GridBagConstraints();
            textFieldGBC.anchor = GridBagConstraints.LINE_START;
            textFieldGBC.gridwidth=GridBagConstraints.RELATIVE;
            if (firstArg)
            {
                textFieldGBC.insets = new Insets(10, 0, 0, 0);
            }            
            fieldPanel.add(argComponent, textFieldGBC);

            if (argDef.getDataType() == ArgumentDefinitionFactory.FILE_TO_READ)
            {
                addFileChooser(argDef.getArgumentName(), true, helper, fieldPanel, false);
            }

            argComponentMap.put(argDef, argComponent);
        }
        else
        {
            fieldPanel = createArgumentFieldPanel(argDef, firstArg, helper, fieldValue);
        }

        JPanel labelPanel = new JPanel();
        GridBagConstraints labelPanelGBC = new GridBagConstraints();
        labelPanelGBC.gridwidth=GridBagConstraints.RELATIVE;
        labelPanelGBC.anchor = GridBagConstraints.LINE_END;

        GridBagConstraints labelGBC = new GridBagConstraints();
        labelGBC.anchor = GridBagConstraints.LINE_END;
        labelGBC.gridwidth=GridBagConstraints.RELATIVE;

        GridBagConstraints helpGBC = new GridBagConstraints();
        helpGBC.anchor = GridBagConstraints.LINE_START;



        labelGBC.insets = new Insets(0, 0, 0, 5);
        if (firstArg)
        {
            labelGBC.insets = new Insets(10, 0, 0, 5);
            firstArg = false;
        }

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
        allFieldsPanel.add(fieldPanel, fieldGBC);
    }

    /**
     * Add a filechooser to a field
     * @param argName
     * @param isMulti
     * @param helper
     * @param fieldPanel
     */
    protected void addFileChooser(String argName, boolean isMulti, ListenerHelper helper, JPanel fieldPanel,
                                  boolean choosesDir)
    {
        JButton chooserButton = new JButton(TextProvider.getText("BROWSE_DOTDOTDOT"));
        GridBagConstraints buttonGBC = new GridBagConstraints();
        buttonGBC.gridwidth = GridBagConstraints.REMAINDER;
        chooserButton.setActionCommand(argName);
        if (choosesDir)
        {
            if (isMulti)
                helper.addListener(chooserButton, "buttonChooseMultiDir_actionPerformed");
            else
                helper.addListener(chooserButton, "buttonChooseSingleDir_actionPerformed");
        }
        else
        {
            if (isMulti)
                helper.addListener(chooserButton, "buttonChooseMultiFile_actionPerformed");
            else
                helper.addListener(chooserButton, "buttonChooseSingleFile_actionPerformed");
        }

        fieldPanel.add(chooserButton, buttonGBC);
    }

    /**
     * This is the method that knows how to handle each argument type differently.  This is the one that
     * should be overridden by a child class that is aware of more arg types
     * @param argDef
     * @param firstArg
     * @param helper
     * @param fieldValue
     * @return
     */
    protected JPanel createArgumentFieldPanel(CommandLineArgumentDefinition argDef,
                                                  boolean firstArg,
                                                  ListenerHelper helper,
                                                  String fieldValue)
    {
        JPanel fieldPanel = new JPanel();
        GridBagConstraints fieldGBC = new GridBagConstraints();
        fieldGBC.gridwidth=GridBagConstraints.REMAINDER;
        fieldGBC.anchor = GridBagConstraints.LINE_START;

        boolean shouldAddFileChooser = false;
        boolean fileChooserChoosesDir = false;
        JComponent argComponent = null;

        boolean fieldHasValue = (fieldValue != null && fieldValue.length() > 0);
        Object defaultValue = argDef.getDefaultValue();

        //special handling for each data type
        switch (argDef.getDataType())
        {
            case ArgumentDefinitionFactory.ENUMERATED:
                JComboBox enumeratedComboBox = new JComboBox();
                for (String allowedValue :
                        ((EnumeratedValuesArgumentDefinition) argDef).getEnumeratedValues())
                {
                    enumeratedComboBox.addItem(allowedValue);
                }

                if (fieldHasValue)
                    enumeratedComboBox.setSelectedItem(fieldValue);
                argComponent = enumeratedComboBox;
                break;
            case ArgumentDefinitionFactory.BOOLEAN:
                JComboBox booleanComboBox = new JComboBox();
                booleanComboBox.addItem("true");
                booleanComboBox.addItem("false");

                if (fieldHasValue)
                    booleanComboBox.setSelectedItem(fieldValue);
                argComponent = booleanComboBox;
                break;
            case ArgumentDefinitionFactory.DECIMAL:
                JTextField decimalTextField = new JTextField();
                decimalTextField.setPreferredSize(new Dimension(70, 20));
                decimalTextField.setMinimumSize(new Dimension(70, 20));

                if (fieldHasValue && (defaultValue == null || !fieldValue.equals(defaultValue.toString())))
                    decimalTextField.setText(fieldValue);
                argComponent = decimalTextField;
                break;
            case ArgumentDefinitionFactory.INTEGER:
                JTextField intTextField = new JTextField();
                intTextField.setPreferredSize(new Dimension(50, 20));
                intTextField.setMinimumSize(new Dimension(50, 20));

                if (fieldHasValue && (defaultValue == null || !fieldValue.equals(defaultValue.toString())))
                    intTextField.setText(fieldValue);
                argComponent = intTextField;
                break;
            case ArgumentDefinitionFactory.FILE_TO_READ:
            case ArgumentDefinitionFactory.FILE_TO_WRITE:
            case ArgumentDefinitionFactory.DIRECTORY_TO_READ:
                JTextField fileTextField = new JTextField();
                fileTextField.setPreferredSize(new Dimension(225, 20));
                fileTextField.setMinimumSize(new Dimension(225, 20));

                if (fieldHasValue)
                    fileTextField.setText(fieldValue);
                argComponent = fileTextField;
                shouldAddFileChooser = true;
                if (argDef.getDataType() == ArgumentDefinitionFactory.DIRECTORY_TO_READ)
                    fileChooserChoosesDir = true;
                break;
            default:
                JTextField argTextField = new JTextField();
                argTextField.setPreferredSize(new Dimension(120, 20));
                argTextField.setMinimumSize(new Dimension(120, 20));

                if (fieldHasValue)
                    argTextField.setText(fieldValue);
                argComponent = argTextField;
                break;
        }

        GridBagConstraints argComponentGBC = new GridBagConstraints();
        argComponentGBC.anchor = GridBagConstraints.LINE_START;
        if (shouldAddFileChooser)
            argComponentGBC.gridwidth=GridBagConstraints.RELATIVE;
        else
            argComponentGBC.gridwidth=GridBagConstraints.REMAINDER;

        if (firstArg)
        {
            argComponentGBC.insets = new Insets(10, 0, 0, 0);
        }
        fieldPanel.add(argComponent, argComponentGBC);
        //add a file chooser that drives off of and populates the text field,
        //if this is a file data type
        if (shouldAddFileChooser)
            addFileChooser(argDef.getArgumentName(), false, helper, fieldPanel, fileChooserChoosesDir);

        argComponentMap.put(argDef, argComponent);


        return fieldPanel;
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
            String argName = argDef.getArgumentName();


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
                argHelpDialog.setModal(true);
                argHelpDialog.setAlwaysOnTop(true);
            }

            argHelpTextArea.setText(argDef.getHelpText());
            argHelpDialog.setTitle("Help for argument '" + argName + "'");

            argHelpTextArea.setVisible(true);
            argHelpDialog.setVisible(true);
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

        for (CommandLineArgumentDefinition argDef : argComponentMap.keySet())
        {
            JComponent argComponent = argComponentMap.get(argDef);

            String argValue;
            switch (argDef.getDataType())
            {
                case ArgumentDefinitionFactory.BOOLEAN:
                case ArgumentDefinitionFactory.ENUMERATED:
                    argValue = (String) ((JComboBox) argComponent).getSelectedItem();
                    break;
                default:
                    argValue = ((JTextField) argComponent).getText();
            }

            if (argValue != null && argValue.length() > 0)
            {
                if (argDef.getArgumentName().equals(
                        CommandLineArgumentDefinition.UNNAMED_PARAMETER_VALUE_SERIES_ARGUMENT))
                {
                    _log.debug("Unnamed series.  Before:\n" + argValue);
                    argValue = argValue.trim();
                    argValue = argValue.replaceAll(" ", CommandLineModule.UNNAMED_ARG_SERIES_SEPARATOR);
                    _log.debug("After:\n" + argValue);                    

                }

                argNameValueMap.put(argDef.getArgumentName(), argValue);
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
     * get the chosen file from the appropriate file chooser.  If single file, that becomes the text field
     * value.  If multi, add this to the existing text field value (if there is one)
     * @param event
     * @param isMulti
     */
    public void chooseSingleOrMultiFileOrDir(ActionEvent event, boolean isMulti, boolean isDir)
    {
        String argName = event.getActionCommand();
        for (CommandLineArgumentDefinition argDef : argComponentMap.keySet())
        {
            if (argDef.getArgumentName().equals(argName))
            {
                JTextField fileTextField = (JTextField) argComponentMap.get(argDef);
                JFileChooser fc = new JFileChooser();
                fc.setMultiSelectionEnabled(isMulti);
                String currentFieldValue = fileTextField.getText();
                File directory = null;
                if (currentFieldValue != null &&
                    currentFieldValue.length() > 0)
                {
                    //if multiple files selected, get the first one
                    if (currentFieldValue.contains(" "))
                        currentFieldValue = currentFieldValue.substring(0, currentFieldValue.indexOf(" "));
                    File currentFile = new File(currentFieldValue);
                    fc.setSelectedFile(currentFile);
                    directory = currentFile.getParentFile();
                }
                else
                {
                    directory = new File (".");
                }
                if (directory != null && directory.exists())
                    fc.setCurrentDirectory(directory);
                if (isDir)
                    fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                
                int chooserStatus = fc.showOpenDialog(this);
                //if user didn't hit OK, ignore
                if (chooserStatus != JFileChooser.APPROVE_OPTION)
                    break;
                File[] files = fc.getSelectedFiles();
                if (!isMulti)
                    files = new File[] { fc.getSelectedFile() };
                if (null != files && files.length > 0)
                {
                    StringBuffer newFileTextBuf = new StringBuffer();
                    for (int i=0; i<files.length; i++)
                    {
                        if (i>0)
                            newFileTextBuf.append(" ");
                        newFileTextBuf.append(files[i].getAbsolutePath());
                    }
                    if (isMulti && fileTextField.getText() != null && fileTextField.getText().length() > 0)
                        fileTextField.setText(fileTextField.getText() + " " + newFileTextBuf.toString());
                    else
                        fileTextField.setText(newFileTextBuf.toString());
                }
                break;
            }
        }
    }

    /**
     * action method for choosing multiple dirs
     * @param event
     */
    public void buttonChooseMultiDir_actionPerformed(ActionEvent event)
    {
        chooseSingleOrMultiFileOrDir(event, true, true);
    }

    /**
     * action method for choosing a single dir
     * @param event
     */
    public void buttonChooseSingleDir_actionPerformed(ActionEvent event)
    {
        chooseSingleOrMultiFileOrDir(event, false, true);
    }

    /**
     * action method for choosing multiple files
     * @param event
     */
    public void buttonChooseMultiFile_actionPerformed(ActionEvent event)
    {
        chooseSingleOrMultiFileOrDir(event, true, false);
    }

    /**
     * action method for choosing a single file
     * @param event
     */
    public void buttonChooseSingleFile_actionPerformed(ActionEvent event)
    {
        chooseSingleOrMultiFileOrDir(event, false, false);
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


        argsSpecified=true;
        done=true;

        notifyDone(event);
    }

    /**
     * Get rid of everything
     */
    protected void disposeAllComponents()
    {
        if (argHelpDialog != null)
        {
            argHelpDialog.setVisible(false);
            argHelpDialog.dispose();
        }

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
}

