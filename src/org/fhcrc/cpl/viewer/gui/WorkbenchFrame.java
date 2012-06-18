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

import org.fhcrc.cpl.toolbox.datastructure.Pair;
import org.fhcrc.cpl.viewer.feature.*;
import org.fhcrc.cpl.viewer.feature.extraction.FeatureFinder;
import org.fhcrc.cpl.viewer.util.SharedProperties;
import org.fhcrc.cpl.toolbox.filehandler.TempFileManager;
import org.fhcrc.cpl.toolbox.proteomics.MSRun;
import org.fhcrc.cpl.viewer.gui.MSImageComponent;
import org.fhcrc.cpl.viewer.CommandFileRunner;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.Application;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.toolbox.gui.ImagePanel;
import org.fhcrc.cpl.toolbox.gui.widget.SplashFrame;
import org.fhcrc.cpl.viewer.ViewerUserManualGenerator;
import org.fhcrc.cpl.viewer.quant.gui.QuantitationReviewer;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.commandline.CommandLineModuleExecutionException;
import org.fhcrc.cpl.toolbox.proteomics.feature.Spectrum;
import org.fhcrc.cpl.toolbox.proteomics.feature.FeatureSet;
import org.systemsbiology.jrap.stax.MZXMLFileInfo;
import org.systemsbiology.jrap.stax.MSInstrumentInfo;
import org.systemsbiology.jrap.stax.DataProcessingInfo;
import org.systemsbiology.jrap.stax.SoftwareInfo;

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Clipboard;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.io.StringWriter;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import javax.imageio.ImageIO;


/**
 * User: mbellew
 * Date: May 20, 2004
 * Time: 3:01:47 PM
 */
public class WorkbenchFrame extends JFrame implements PropertyChangeListener
{
    //hardcoded URLs for external support and help links
    public static final String SUPPORT_URL = "https://www.labkey.org/Project/msInspect/begin.view";
    public static final String MSINSPECT_WEBSITE_URL = "https://proteomics.fhcrc.org/CPL/msinspect.html";
    protected static String helpURL = MSINSPECT_WEBSITE_URL;

    private static String startup_appTitle = "msInspect";
    private String appTitle = startup_appTitle;
    JMenuBar menuBar = null;
    JPanel panel2D = null;
    SpectrumComponent panelSpectrum = null;

    public JPanel contentPanel;
    public JLabel messageLabel;
    public JPanel statusPanel;
    public JScrollPane leftScrollPane;
    public JPanel bottomPane;
    public JSplitPane outerSplitPane;
    public JPanel rightPane;
    public GridBagConstraints rightPaneGBC;
    //Inner split pane, splits between properties pane and main image
    public JSplitPane topLeftSplitPane;
    //outer split pane, splits between topLeftSplitPane and detail pane
    public JSplitPane topRightSplitPane;

    public MSImageComponent imageComponent;

    public JPanel topPanel;
    public JTable propertiesTable;
    public JTree filesTree;
    protected DefaultMutableTreeNode filesRoot = new DefaultMutableTreeNode("NYI");
    protected DefaultTreeModel filesModel = new DefaultTreeModel(filesRoot);
    public JTabbedPane infoTabbedPane;
    public JScrollPane propertiesScrollPane;
    protected PropertiesPane propertiesPane = new PropertiesPane();

    //
    // Menu actions
    //

    // file
    public Action openFileAction = new OpenFileAction(null);
    public Action openFastaAction = new OpenFastaDialog.OpenFastaAction(null);

    public Action runInfoAction = new MSRunInfoAction("File Properties");
    public Action saveImageAction = new MSImageComponent.SaveImageAction();
    public Action exitAction = new ExitAction();

    // tools
    public Action copyScanAction = new CopyScanAction();
    public Action extractFeaturesAction = new ExtractFeatureRangesAction();
    public FeatureSelectionFrame selectFeaturesAction = new FeatureSelectionFrame();
    public Action displayHeatMapAction = new HeatMapAction();
    public Action coverageCalculatorAction = new CoverageCalculatorAction();
    public Action autoZoomAction = new AutoZoomAction();
    public Action lockYAxisAction = new LockYAxisAction();

    public Action singleScanExtractorAction = new ChooseFeatureExtractorAction("Single scan", FeatureStrategyWavelet.class);
    public Action centroidedScanExtractorAction = new ChooseFeatureExtractorAction("Centroided scan", FeatureStrategyCentroided.class)
    {
        {
            ApplicationContext.addPropertyChangeListener(SharedProperties.MS_RUN, new PropertyChangeListener()
            {
                public void propertyChange(PropertyChangeEvent evt)
                {
                    MSRun run = (MSRun)evt.getNewValue();
                    boolean enabled = null == run || run.getHeaderInfo().getDataProcessing().getCentroided() != 1;
                    setEnabled(enabled);
                }
            });
        }
    };
    public Action combinedExtractorAction = new ChooseFeatureExtractorAction("Combined", FeatureStrategyCombined.class);
    public Action grossFeaturesAction = new ChooseFeatureExtractorAction("  gross features (intermediate result)", FeatureStrategyGrossFeatures.class);
    public Action peakClustersExtractorAction = new ChooseFeatureExtractorAction("2D Peak Alignment", FeatureFinder.DEFAULT_FEATURE_FINDING_CLASS);
    public Action peaksExtractorAction = new ChooseFeatureExtractorAction("  peaks (intermediate result)", FeatureStrategyPeaks2.class);

    public Action detailOptionsDialogAction = new MSDetailPanel.ShowDialogAction();
    public Action show3DWindowAction = new MSDetailPanel.Show3DAction();
    public Action amtAction = new ProteinMatcherFrame.ProteinMatcherAction();
    public Action runCommandFileAction = new CommandFileRunner.CommandFileRunnerAction();
    public Action runCommandAction = new ChooseCommandDialog.RunCommandAction();
    public Action showSelectedForCIDAction = new SelectedForCIDFinder.SelectedForCIDFinderAction();
    public Action qurateAction = new QurateAction();
    public Action featureViewerAction = new FeatureViewerAction();




    // window
    public Action showPropertiesAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent e)
        {
            if (topLeftSplitPane.getDividerSize() > 0)
                hidePropertiesPane();
            else
                showPropertiesPane();
        }
    };

    public Action showLayerTransparencyAction = new AbstractAction()
    {
        public void actionPerformed(ActionEvent e)
        {
            imageComponent.showLayerTransparencyDialog();
        }
    };

    // Languages
    public JMenu languageMenu;
    // Colors
    public JMenu colorMenu;
    // Fonts
    public JMenu fontMenu;

    // help
    public Action helpAction = new HelpAction();
    public Action supportAction = new SupportAction();
    public Action aboutAction = new AbstractAction("About")
    {
        public void actionPerformed(ActionEvent e)
        {
            WorkbenchFrame.ShowSplashScreen();
        }
    };


    public JMenuItem autoZoomMenuItem;
    public JMenuItem lockYAxisMenuItem;
    public JMenuItem defaultExtratorMenuItem;

    private static Image iconImage = null;



    public WorkbenchFrame()
    {
        super(startup_appTitle);

        try
        {
            Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/WorkbenchFrame.xml",this);
            assert null != contentPanel;
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage("error creating dialog", x);
            throw new RuntimeException(x);
        }

        String revision = (String) Application.getInstance().getProperty("REVISION");
        if (null != revision)
        {
            if (revision.indexOf(' ') >= 0)
                revision = revision.substring(revision.indexOf(' ')).trim();
            appTitle = appTitle + " (build " + revision + ")";
            setTitle(appTitle);
        }

        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setContentPane(contentPanel);
        setSize(Localizer.adjustContainerDimension(new Dimension(800, 600)));
        //adjust the location of the divider between panes for the current locale    
        topLeftSplitPane.setDividerLocation(Localizer.adjustContainerWidth(topLeftSplitPane.getDividerLocation()));

        try
        {
            iconImage = ImageIO.read(WorkbenchFrame.class.getResourceAsStream("icon.gif"));
            setIconImage(iconImage);
        }
        catch (Exception e)
        {
        }

        try
        {
            JMenuBar jmenu = (JMenuBar)Localizer.getSwingEngine(this).render("org/fhcrc/cpl/viewer/gui/WorkbenchFrameMenu.xml");
            for (int i=0 ; i<jmenu.getMenuCount() ; i++)
                jmenu.getMenu(i).getPopupMenu().setLightWeightPopupEnabled(false);
            this.setJMenuBar(jmenu);
            autoZoomMenuItem.doClick();
            defaultExtratorMenuItem.doClick();
        }
        catch (Exception x)
        {
            ApplicationContext.errorMessage(TextProvider.getText("ERROR_LOADING_MENUS"), x);
            throw new RuntimeException(x);
        }

        //Localizer needs to populate all the language and font menu items
        Localizer.initializeLanguageMenu(languageMenu);
        Localizer.initializeFontMenu(fontMenu);

        //populate the available color schemes
        MSImageComponent.initializeColorMenu(colorMenu);

        imageComponent = new MSImageComponent();
        ApplicationContext.setProperty("MSImageComponent", imageComponent);
        imageComponent.setScrollPane(leftScrollPane);


        MSDetailPanel detail = new MSDetailPanel();
        // replace rightPane
        detail.setPreferredSize(rightPane.getPreferredSize());
        detail.setMinimumSize(rightPane.getMinimumSize());
        detail.setMaximumSize(rightPane.getMaximumSize());
        topPanel.remove(rightPane);
        topPanel.add(detail, rightPaneGBC);

        //propertiesScrollPane.add(propertiesPane);
        infoTabbedPane.addTab(TextProvider.getText("PROPERTIES"), propertiesPane.getComponent());

        panelSpectrum = new SpectrumComponent();
        outerSplitPane.setBottomComponent(panelSpectrum.getComponent());
        topRightSplitPane.setRightComponent(detail);
        topRightSplitPane.setLeftComponent(topLeftSplitPane);

        // UNDONE PLAF ??
        messageLabel.setBackground(Color.WHITE);
        messageLabel.setFont(Font.decode("verdana plain 12"));
        messageLabel.setText(" ");

        ApplicationContext.addPropertyChangeListener(SharedProperties.MS_RUN, this);

        this.addComponentListener(new ComponentAdapter()
        {
            public void componentResized(ComponentEvent event)
            {
                _resizeDivider();
            }
            public void componentShown(ComponentEvent event)
            {
                super.componentShown(event);
                _resizeDivider();
            }
        });

        this.invalidate();
        this.doLayout();

        ApplicationContext.addPropertyChangeListener(SharedProperties.SELECTED, new PropertyChangeListener()
        {
            public void propertyChange(PropertyChangeEvent evt)
            {
                updateSelectedObject(evt.getNewValue());
            }
        });
    }




    public void hidePropertiesPane()
    {
        topLeftSplitPane.setDividerSize(0);
        topLeftSplitPane.setDividerLocation(0);
        topLeftSplitPane.setResizeWeight(0);
        infoTabbedPane.setVisible(false);
    }


    public void showPropertiesPane()
    {
        topLeftSplitPane.setDividerSize(5);
        int location = 200;
        if (topLeftSplitPane.getLastDividerLocation() > 20)
            location = topLeftSplitPane.getLastDividerLocation();
        if (getWidth() > 200)
            location = Math.min(getWidth()-100, location);
        topLeftSplitPane.setDividerLocation(location);
        topLeftSplitPane.setResizeWeight(0.5);
        infoTabbedPane.setVisible(true);
    }


    protected void updateSelectedObject(Object o)
    {
        if (null == o)
        {

        }
        else if (o instanceof MSRun)
        {
            MSRun run = (MSRun)o;
            MZXMLFileInfo fileInfo = run.getHeaderInfo();
            if (null == fileInfo)
                return;

            //HashMap m = new HashMap(30);
            ArrayList<Pair<String, String>> list = new ArrayList<Pair<String, String>>();

            list.add(new Pair<String, String>(TextProvider.getText("FILE"), run.getFileName()));

            MSInstrumentInfo inst = fileInfo.getInstrumentInfo();
            if (null != inst)
            {
                list.add(new Pair<String, String>("<html><b>" + TextProvider.getText("INSTRUMENT") + "<b></html>",""));
                list.add(new Pair<String, String>(TextProvider.getText("MANUFACTURER"), inst.getManufacturer()));
                list.add(new Pair<String, String>(TextProvider.getText("DETECTOR"), inst.getDetector()));
                list.add(new Pair<String, String>(TextProvider.getText("MASSANALYSER"), inst.getMassAnalyzer()));
                list.add(new Pair<String, String>(TextProvider.getText("MODEL"), inst.getModel()));
                list.add(new Pair<String, String>(TextProvider.getText("OPERATOR"), valueOf(inst.getOperator())));
                if (inst.getSoftwareInfo() != null) {
                	list.add(new Pair<String, String>(TextProvider.getText("SOFTWARE"), valueOf(inst.getSoftwareInfo().name)));
                }
                list.add(new Pair<String, String>(TextProvider.getText("IONIZATION"), inst.getIonization()));
            }

            DataProcessingInfo data = fileInfo.getDataProcessing();
            if (null != data) //  && (data.getCentroided() != -1 || data.getChargeDeconvoluted() != -1 || data.getDeisotoped() != -1))
            {
                list.add(new Pair<String, String>("<html><b>" + TextProvider.getText("DATA_PROCESSING") + "<b></html>",""));
                list.add(new Pair<String, String>(TextProvider.getText("CENTROIDED"), String.valueOf(data.getCentroided())));
                list.add(new Pair<String, String>(TextProvider.getText("CHARGE_DECONVOLUTED"), String.valueOf(data.getChargeDeconvoluted())));
                list.add(new Pair<String, String>(TextProvider.getText("DEISOTOPED"), String.valueOf(data.getDeisotoped())));
                list.add(new Pair<String, String>(TextProvider.getText("INTENSITY_CUTOFF"), String.valueOf(data.getIntensityCutoff())));
                list.add(new Pair<String, String>(TextProvider.getText("SPOT_INTEGRATION"), String.valueOf(data.getSpotIntegration())));
                java.util.List<SoftwareInfo> soft = data.getSoftwareUsed();
                if (null != soft)
                {
                    String key = TextProvider.getText("SOFTWARE");
                    for (int i=0 ; i<soft.size() ; i++)
                    {
                        list.add(new Pair<String, String>(key, valueOf(soft.get(i).name)));
                        key = "";
                    }
                }
            }
            propertiesPane.setProperties(list);
        }
        else if (o instanceof FeatureSet)
        {
            FeatureSet fs = (FeatureSet)o;
            Map<Object,Object> tm = new TreeMap<Object,Object>(fs.getProperties());
            java.util.List<Map.Entry<Object,Object>> list = new ArrayList<Map.Entry<Object,Object>>();
            // add File to beginning of list
            if (fs.getSourceFile()!=null)
                list.add(new Pair<Object,Object>("File", fs.getSourceFile().getName()));
            list.addAll(tm.entrySet());
            propertiesPane.setProperties(list);
        }
        else
        {
            propertiesPane.setProperties(o);
        }
    }


    static String valueOf(Object o)
    {
        return null == o ? "" : String.valueOf(o);
    }

    public static String getAppName()
    {
        return startup_appTitle;
    }

    private void _resizeDivider()
    {
        if (topLeftSplitPane.getDividerSize() == 0) // hidden
            topLeftSplitPane.setDividerLocation(0);
/*		int width = getContentPane().getWidth();

		if (null != topLeftSplitPane)
			{
			int preferredWidth = 123;
			int dividerSize = topLeftSplitPane.getDividerSize();
			int extra = 5;// = 10;
			int location = width-preferredWidth-dividerSize-extra;
			location = Math.max(100, location);
			topLeftSplitPane.setDividerLocation(location);
			}
		else
			{
//			rightPane.setSize(123, rightPane.getHeight());
			} */
    }

    /**
     * Handles all the Tools actions that select peptide finder types
     */
    class ChooseFeatureExtractorAction extends AbstractAction
    {
        Class<?> _class;


        public ChooseFeatureExtractorAction(String text, Class<?> c)
        {
            super(null == text ? c.getName() : text);
            _class = c;
        }


        public void actionPerformed(ActionEvent event)
        {
            FeatureExtractor.setDefault(_class);
            Object source = event.getSource();
            if (source instanceof JCheckBoxMenuItem)
                ((JCheckBoxMenuItem) source).setState(true);
            else if (source instanceof JRadioButtonMenuItem)
                ((JRadioButtonMenuItem) source).setSelected(true);
        }
    }

    /**
     * handles Tools->Copy Scan
     */
    class CopyScanAction extends AbstractAction
    {
        public CopyScanAction()
        {
            super("Copy Scan");
        }

        public void actionPerformed(ActionEvent event)
        {
            try
            {
                MSRun.MSScan scan = (MSRun.MSScan) ApplicationContext.getProperty(SharedProperties.MS_SCAN);
                if (null == scan)
                    return;
                float[][] spectrum = scan.getSpectrum();

                StringWriter sw = new StringWriter();
                Spectrum.CopyToTSV(spectrum, sw, true);
                StringSelection sel = new StringSelection(sw.toString());
                Clipboard clip = Toolkit.getDefaultToolkit().getSystemClipboard();
                clip.setContents(sel, sel);
            }
            catch (IOException x)
            {

            }
        }
    }

    /**
     * Handles the File->Exit action
     */
    public class ExitAction extends AbstractAction
    {
        public ExitAction()
        {
            super("Exit");
        }

        public void actionPerformed(ActionEvent event)
        {
            Application.quit();
        }
    }

    public static class UrlAction extends AbstractAction
    {
        String url;

        UrlAction(String name, String url)
        {
            super(name);
            this.url = url;
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                HtmlViewerPanel.showURLInDialog(url, "");
                return;
            }
            catch (IOException x)
            {;
            }
            JOptionPane.showConfirmDialog(null, url);
        }
    }


    /**
     * Handle selection of &quot;Display Heat Map...&quot; menu item
     */
    public class HeatMapAction extends AbstractAction implements PropertyChangeListener
    {
        public HeatMapAction()
        {
            super("Display Heat Map...");
            ApplicationContext.addPropertyChangeListener(SharedProperties.FEATURE_RANGES, this);
            setEnabled(null != HeatMapFrame.getDisplayedFeatureSet());
        }

        public void actionPerformed(ActionEvent evt)
        {
            if (null == HeatMapFrame.getDisplayedFeatureSet())
            {
                ApplicationContext.errorMessage("No feature set currently displayed.", null);
                return;
            }
            HeatMapFrame.getInstance().display();
            return;
        }

        public void propertyChange(PropertyChangeEvent event)
        {
            setEnabled(null != HeatMapFrame.getDisplayedFeatureSet());
        }
    }

    /**
     * Handles the Auto Zoom checkbox
     */
    public class AutoZoomAction extends AbstractAction
    {
        public AutoZoomAction()
        {
            super("Auto zoom");
        }

        public void actionPerformed(ActionEvent event)
        {
            Boolean zoomProperty = (Boolean) ApplicationContext.getProperty(SharedProperties.AUTO_ZOOM);
            boolean zoom = null != zoomProperty ? !zoomProperty.booleanValue() : true;
            Object source = event.getSource();
            if (source instanceof JCheckBoxMenuItem)
            {
                JCheckBoxMenuItem menu = (JCheckBoxMenuItem) event.getSource();
                menu.setState(zoom);
            }
            ApplicationContext.setProperty(SharedProperties.AUTO_ZOOM, new Boolean(zoom));
        }
    }

    /**
     * Handles the "Lock Y Scale" checkbox
     */
    public class LockYAxisAction extends AbstractAction
    {
        public LockYAxisAction()
        {
            super("Lock Y Axis");
        }

        public void actionPerformed(ActionEvent event)
        {
            Boolean lockProperty = (Boolean) ApplicationContext.getProperty(SharedProperties.LOCK_Y_AXIS);
            boolean lock = null != lockProperty ? !lockProperty.booleanValue() : true;
            Object source = event.getSource();
            if (source instanceof JCheckBoxMenuItem)
            {
                JCheckBoxMenuItem menu = (JCheckBoxMenuItem) event.getSource();
                menu.setState(lock);
            }
            ApplicationContext.setProperty(SharedProperties.LOCK_Y_AXIS, new Boolean(lock));
        }
    }

    public static class HelpAction extends UrlAction
    {
        public HelpAction()
        {
            super("Help Topics", helpURL);
        }

        public void actionPerformed(ActionEvent e)
        {
            File helpFile = TempFileManager.createTempFile("help.html", this);
            if (!helpFile.exists())
            {
                try
                {
                    PrintWriter helpPW = new PrintWriter(helpFile);
                    new ViewerUserManualGenerator().generateFullManual(helpPW);
                    helpURL = "file://" + helpFile.getAbsolutePath();
                    this.url = helpURL;
                }
                catch (Exception ex)
                {
                    ApplicationContext.setMessage("Failed to generate user manual, defaulting to msInspect website");
                }
            }
            super.actionPerformed(e);
        }
    }


    public static class SupportAction extends UrlAction
    {
        public SupportAction()
        {
            super("Support Forum", SUPPORT_URL);
        }
    }


    public void setMessage(String message)
    {
        if (EventQueue.isDispatchThread())
        {
            if (null == message || 0 == message.length())
                message = " ";
            messageLabel.setText(message);
        }
        else
        {
            final String msg = message;
            EventQueue.invokeLater(new Runnable()
            {
                public void run()
                {
                    setMessage(msg);
                }
            });
        }
    }

    public void disposeNoExit()
    {
        super.dispose();
    }

    public void dispose()
    {
        System.exit(0);
    }


    public void propertyChange(PropertyChangeEvent event)
    {
        if (SharedProperties.MS_RUN.equals(event.getPropertyName()))
        {
            MSRun run = (MSRun)event.getNewValue();
            if (null == run)
                setTitle(appTitle);
            else
            {
                String name = run.getFileName();
                setTitle(appTitle + " -- " + name);
            }
        }
    }


    public static JFrame ShowSplashScreen()
    {
        SplashFrame splashFrame = new SplashFrame(WorkbenchFrame.class.getResource("splash.gif"));
        splashFrame.setVisible(true);
        return splashFrame;
    }

    public static class QurateAction extends AbstractAction
    {
        public void actionPerformed(ActionEvent event)
        {
            QuantitationReviewer quantReviewer = new QuantitationReviewer();
            quantReviewer.setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
            quantReviewer.setVisible(true);
        }
    }
}
