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

import org.fhcrc.cpl.toolbox.commandline.CommandLineModule;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleUtilities;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.gui.ListenerHelper;
import org.fhcrc.cpl.viewer.commandline.ViewerCommandLineModuleDiscoverer;
import org.fhcrc.cpl.viewer.Application;

import javax.swing.*;
import java.util.Map;
import java.util.Arrays;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.StringWriter;
import java.io.PrintWriter;

/**
     * Simple dialog for choosing a command to specify arguments for
 */
public class ChooseCommandDialog extends JDialog
{
    public JComboBox commandBox;
    public Map<String, CommandLineModule> moduleMap;
    public JTextArea descriptionPanel;

    public static final int height = 175;
    public static final int width = 300;

    protected boolean done = false;
    protected CommandLineModule chosenModule = null;

    public JButton buttonGo;
    public JButton buttonCancel;

    public JButton fakeButton;

    public ChooseCommandDialog()
    {
        super();

        setTitle(TextProvider.getText("CHOOSE_COMMAND"));

        JPanel contentPanel = new JPanel();
        GridBagConstraints contentPanelGBC = new GridBagConstraints();
        contentPanelGBC.gridwidth = GridBagConstraints.REMAINDER;
        contentPanel.setLayout(new GridBagLayout());
        add(contentPanel);

        setPreferredSize(new Dimension(width, height));
        setSize(new Dimension(width, height));
        setMinimumSize(new Dimension(width, height));

        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int centerH = screenSize.width / 2;
        int centerV = screenSize.height / 2;
        this.setLocation(centerH - width / 2, centerV - height / 2);

        commandBox = new JComboBox();
        GridBagConstraints commandBoxGBC = new GridBagConstraints();
        commandBoxGBC.gridwidth = GridBagConstraints.REMAINDER;
        commandBoxGBC.insets = new Insets(10, 0, 10, 0);
        contentPanel.add(commandBox, commandBoxGBC);

        moduleMap = ViewerCommandLineModuleDiscoverer.getSingletonInstance().findAllCommandLineModules();
        String[] commandsArray = moduleMap.keySet().toArray(new String[moduleMap.size()]);
        Arrays.sort(commandsArray);
        for (String command : commandsArray)
        {
            commandBox.addItem(command);
        }

        fakeButton = new JButton("fake");

        buttonGo = new JButton(TextProvider.getText("OK"));
        GridBagConstraints buttonGoGBC = new GridBagConstraints();
        buttonGoGBC.insets = new Insets(0, 30, 0, 0);
        buttonGoGBC.gridwidth = GridBagConstraints.RELATIVE;
        getRootPane().setDefaultButton(buttonGo);

        buttonCancel = new JButton(TextProvider.getText("CANCEL"));
        GridBagConstraints buttonCancelGBC = new GridBagConstraints();
        buttonCancelGBC.gridwidth = GridBagConstraints.REMAINDER;

        contentPanel.add(buttonGo, buttonGoGBC);
        contentPanel.add(buttonCancel, buttonCancelGBC);

        descriptionPanel = new JTextArea(4, 20);
//            JScrollPane scrollPane = new JScrollPane(descriptionPanel);
        descriptionPanel.setOpaque(false);
        descriptionPanel.setEditable(false);
        descriptionPanel.setLineWrap(true);
        descriptionPanel.setWrapStyleWord(true);
        GridBagConstraints descriptionGBC = new GridBagConstraints();
        descriptionGBC.gridwidth = GridBagConstraints.REMAINDER;

        CommandLineModule firstCommand = moduleMap.get(commandsArray[0]);
        descriptionPanel.setText(firstCommand.getShortDescription());

        contentPanel.add(descriptionPanel, descriptionGBC);

        ListenerHelper helper = new ListenerHelper(this);

        helper.addListener(commandBox, "commandBox_actionPerformed");
        helper.addListener(buttonGo, "buttonGo_actionPerformed");
        helper.addListener(buttonCancel, "buttonCancel_actionPerformed");

        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    }

    public CommandLineModule chooseCommand()
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

        setVisible(false);
        dispose();

        return chosenModule;
    }

    public void addDoneListener(ActionListener listener)
    {
        fakeButton.addActionListener(listener);
    }

    public void commandBox_actionPerformed(ActionEvent event)
    {
        String command = (String) commandBox.getSelectedItem();
        CommandLineModule clm = moduleMap.get(command);
        descriptionPanel.setText(clm.getShortDescription());
    }

    public void buttonGo_actionPerformed(ActionEvent event)
    {
        String command = (String) commandBox.getSelectedItem();

        chosenModule = moduleMap.get(command);
        done=true;
        notifyDone(event);

    }

    public void buttonCancel_actionPerformed(ActionEvent event)
    {
        done=true;
        notifyDone(event);
    }

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
     * action for the menu item that kicks macro-running off
     */
    public static class RunCommandAction extends AbstractAction
    {
        protected ChooseCommandDialog chooseCommandDialog;
        protected ViewerInteractiveModuleFrame interactFrame;
        protected CommandLineModule module;
        public void actionPerformed(ActionEvent event)
        {
            chooseCommandDialog =
                    new ChooseCommandDialog();
//            CommandLineModule module = chooseCommandDialog.chooseCommand();
            chooseCommandDialog.addDoneListener(new ChooseModuleListener());
            chooseCommandDialog.setVisible(true);

        }

        protected class ChooseModuleListener implements ActionListener
        {
            public void actionPerformed(ActionEvent event)
            {
                chooseCommandDialog.setVisible(false);
                module = chooseCommandDialog.chosenModule;

                chooseCommandDialog.dispose();

                if (module == null)
                    return;

                interactFrame = new ViewerInteractiveModuleFrame(module, true, null);
                interactFrame.addDoneListener(new ExecuteModuleListener());
                interactFrame.setVisible(true);
            }
        }

        protected class ExecuteModuleListener implements ActionListener
        {
            public void actionPerformed(ActionEvent event)
            {
                interactFrame.setVisible(false);
                if (!interactFrame.argsSpecified)
                {
                    interactFrame.dispose();
                    return;
                }

                interactFrame.dispose();

                try
                {
                    module.execute();
                }
                catch (CommandLineModuleExecutionException ex)
                {
                    String message = TextProvider.getText("ERROR_RUNNING_COMMAND_COMMAND", module.getCommandName());

                    message = message + "\n" + ex.getMessage() + "\n";

                    StringWriter sw = new StringWriter();
                    PrintWriter w = new PrintWriter(sw);
                    ex.printStackTrace(w);
                    w.flush();
                    message += "\n";
                    message += sw.toString();
                    message += CommandLineModuleUtilities.createFailureReportAndPrompt(module, ex,
                            Application.isLogEnabled(), Application.getLogFile(),
                            Application.FAILURE_REPORT_ERRORMESSAGE_TEXT, Application.FAILURE_REPORT_HEADER_TEXT);
                    JOptionPane.showMessageDialog(ApplicationContext.getFrame(), message, "Information",
                            JOptionPane.INFORMATION_MESSAGE);
                }
            }
        }
    }
}
