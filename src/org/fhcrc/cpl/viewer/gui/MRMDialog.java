package org.fhcrc.cpl.viewer.gui;
import org.fhcrc.cpl.viewer.Localizer;
import org.fhcrc.cpl.viewer.MSRun;
import org.fhcrc.cpl.viewer.commandline.CommandLineModule;
import org.fhcrc.cpl.viewer.commandline.CommandLineModuleDiscoverer;
import org.fhcrc.cpl.viewer.commandline.CLMUserManualGenerator;
import org.fhcrc.cpl.viewer.mrm.MRMCommandLineModule;
import org.fhcrc.cpl.viewer.util.Clusterer2D;
import org.fhcrc.cpl.viewer.util.ElutionDataPoint;
import org.fhcrc.cpl.toolbox.TempFileManager;
import org.fhcrc.cpl.toolbox.gui.chart.PanelWithChart;
import org.fhcrc.cpl.toolbox.gui.HtmlViewerPanel;
import org.fhcrc.cpl.viewer.mrm.*;
import org.fhcrc.cpl.viewer.mrm.utilities.*;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.toolbox.TextProvider;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.AxisChangeEvent;
import org.jfree.chart.event.AxisChangeListener;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.*;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableModel;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Line2D;
import java.io.File;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
     * GUI
 */
public class MRMDialog extends JFrame {
    protected JLabel titleText;
    protected JLabel elutionTableLabel;
    protected MRMTransition[] _mrmTransitions;
    protected int currentTransitionIndex = 0;
    protected MSRun _run;
    protected float _precursorDiscoveryMzTolerance;
    protected float _daughterMzTolerance;
    protected float _precursorChromatogramWindow;
    protected File _mzXMLFile;
    protected JFileChooser _mzXMLFileChooser;
    protected JFileChooser _outputFileChooser;
    protected boolean _traceAllFragments;
    protected ListSelectionListener _ptmlsl;
    protected transitionListSelectionListener _tlsl;
    protected boolean _synclh;

    //Widgets

    public JList listTransition;
    public JScrollPane listScrollPane;
    public JScrollPane peaksScrollPane;
    public static JTable peaksTable;

    public JLabel fileNameLabel;
    public JPanel precursorContainerPanel;
    public JPanel precursorContainerContainerPanel;
    public JPanel daughterContainerContainerPanel;
    public JPanel daughterContainerPanel;
    public JLabel topGraphLabel;

    public JButton buttonZC;
    public JButton buttonPrev;
    public JButton buttonNext;
    public JButton buttonSave;
    public JButton buttonAccept;
    public JButton buttonReject;
    public JButton buttonZoom;
    public JButton buttonFindMate;
    public JButton buttonBoost;
    public JButton buttonTiming;
    public JButton buttonRejectGroup;
    public JMenuBar menuBarMain;
    public JMenu    menuFile;
    public JMenuItem menuItemQuit;
    public JMenuItem menuItemOpen;
    public JMenu menuOptions;
    public JMenu menuHelp;
    public JMenuItem menuItemArguments;
    public JMenuItem menuItemOptions;
    public JMenuItem menuItemDefinitions;
    public JMenuItem menuItemTips;
    public JMenuItem menuItemSICtolerance;
    public JMenuItem menuItemPrecDiscTol;
    public JMenuItem menuItemDaughterDiscTol;
    public JMenuItem menuItemTraceAllFragments;
    public JMenuItem menuItemChangeStrategy;
    public JMenuItem menuItemSyncLH;
    public JMenuItem menuItemPMin;
    public JMenuItem menuItemAMin;

    ListenerHelper helper = null;

    public boolean _sim = false;
    public boolean _conn = false;
    public float _sictol = 0.0f;
    public Class _ecurveclass = null;
    public Container contentPanel = null;
    public MRMTransition transitionOnPlot = null;
    public TransitionDefinitionHeader transDefHeader = null;
    public float _minPeakCutoff;
    public float _minAreaCutoff;

    /* todo: these are totally hokey and needs to be re-thought */
    private static enum whichGraph{Precursor,Daughter};
    private static enum whichAxis{Domain,Range};

    private void menuBarInitializations() {
        menuBarMain = new JMenuBar();
        this.setJMenuBar(menuBarMain);
        menuFile = new JMenu("File");
        menuBarMain.add(menuFile);
        menuItemOpen = new JMenuItem("Open");
        menuFile.add(menuItemOpen);
        menuItemQuit = new JMenuItem("Quit");
        menuFile.add(menuItemQuit);
        menuOptions = new JMenu("Options");
        menuItemSICtolerance = new JMenuItem("SIC Tolerance");
        menuOptions.add(menuItemSICtolerance);
        menuItemPrecDiscTol = new JMenuItem("Precursor Tolerance");
        menuOptions.add(menuItemPrecDiscTol);
        menuItemDaughterDiscTol = new JMenuItem("Product Tolerance");
        menuOptions.add(menuItemDaughterDiscTol);
        menuItemTraceAllFragments = new JMenuItem("Trace Chromats");
        menuOptions.add(menuItemTraceAllFragments);
        menuItemChangeStrategy = new JMenuItem("Curve Strategy");
        menuOptions.add(menuItemChangeStrategy);
        menuItemSyncLH = new JMenuItem("Synchronize L/H label elution regions");
        menuOptions.add(menuItemSyncLH);
        menuItemPMin = new JMenuItem("Refilter Min Peak Height");
        menuOptions.add(menuItemPMin);
        menuItemAMin = new JMenuItem("Refilter Min Curve Area");
        menuOptions.add(menuItemAMin);
        menuBarMain.add(menuOptions);
        menuHelp = new JMenu("Help");
        menuItemArguments = new JMenuItem("Java commandline options");
        menuItemOptions = new JMenuItem("Options help");
        menuHelp.add(menuItemOptions);
        menuItemDefinitions = new JMenuItem("Definitions");
        menuHelp.add(menuItemDefinitions);
        menuItemTips = new JMenuItem("User Tips");
        menuHelp.add(menuItemTips);
        menuHelp.add(menuItemArguments);
        menuBarMain.add(menuHelp);
    }

    private void listenerHelperInitializations() {
        helper.addListener(menuItemSICtolerance,"buttonSICUpdate_actionPerformed");
        helper.addListener(menuItemPrecDiscTol,"buttonPDT_actionPerformed");
        helper.addListener(menuItemDaughterDiscTol,"buttonDTOL_actionPerformed");
        helper.addListener(menuItemTraceAllFragments,"menuItemTraceAllFragments_actionPerformed");
        helper.addListener(menuItemChangeStrategy,"menuItemChangeStrategy_actionPerformed");
        helper.addListener(menuItemSyncLH,"menuItemSyncLH_actionPerformed");
        helper.addListener(buttonZC, "buttonZC_actionPerformed");
        helper.addListener(menuItemPMin,"menuItemPMin_actionPerformed");
        helper.addListener(menuItemAMin,"menuItemAMin_actionPerformed");
        helper.addListener(buttonNext, "buttonNext_actionPerformed");
        helper.addListener(buttonPrev, "buttonPrev_actionPerformed");
        helper.addListener(buttonSave, "buttonSave_actionPerformed");
        helper.addListener(buttonAccept,"buttonAccept_actionPerformed");
        helper.addListener(buttonBoost,"buttonBoost_actionPerformed");
        helper.addListener(buttonTiming,"buttonTiming_actionPerformed");
        helper.addListener(buttonReject,"buttonReject_actionPerformed");
        helper.addListener(buttonRejectGroup,"buttonRejectGroup_actionPerformed");
        helper.addListener(contentPanel,"contentPanel_mouseClicked");
        helper.addListener(menuItemQuit,"menuItemQuit_actionPerformed");
        helper.addListener(menuItemTips,"menuItemTips_actionPerformed");
        helper.addListener(menuItemDefinitions,"menuItemDefinitions_actionPerformed");
        helper.addListener(menuItemOptions,"menuItemOptions_actionPerformed");
        helper.addListener(menuItemArguments,"menuItemArguments_actionPerformed");
        helper.addListener(menuItemOpen,"menuItemOpen_actionPerformed");
        helper.addListener(buttonZoom,"buttonZoom_actionPerformed");
        helper.addListener(buttonFindMate,"buttonFindMate_actionPerformed");
    }

    /**
      * Initialize the GUI and show the first transition
      * @param mzXMLFile
      */

      public MRMDialog(File mzXMLFile, float precDiscTol, float daughterTol, float chromTol, Class ECurveStrat, boolean traceAllFragments, boolean synclh,float minP, float minA)
     {
         super(TextProvider.getText("MRMer")+" v. "+TextProvider.getText("MRMER_VERSION")+" (build "+ (String) ApplicationContext.getProperty("REVISION")+")");
         frameInit();
         _precursorDiscoveryMzTolerance =precDiscTol;
         _daughterMzTolerance = daughterTol;
         _precursorChromatogramWindow = chromTol;
         _mzXMLFile = mzXMLFile;
         _ecurveclass = ECurveStrat;
         currentTransitionIndex = 0;
         _traceAllFragments = traceAllFragments;
         _synclh = synclh;
         _minPeakCutoff = minP;
         _minAreaCutoff = minA;
         //graphical stuff
         contentPanel = null;
         helper = new ListenerHelper(this);
         try
         {
              contentPanel =
                     Localizer.renderSwixml("org/fhcrc/cpl/viewer/gui/MRMDialog.xml",this);
              contentPanel.setBackground(new Color(255,255,153));

              menuBarInitializations();
              listenerHelperInitializations();
              _mzXMLFileChooser  = new JFileChooser();
              _mzXMLFileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
              _mzXMLFileChooser.setMultiSelectionEnabled(false);
              _mzXMLFileChooser.setFileFilter(
                 new javax.swing.filechooser.FileFilter() {
                    public boolean accept(File f) {
                       if(f.isDirectory()) return true;
                       String fname = f.getName();
                       String parts[] = fname.split("\\.");
                       if(parts.length != 2) return false;
                       return parts[1].equalsIgnoreCase("mzXML");
                    }
                    public String getDescription() {
                       return "mzXML files";
                    }
                 }
              );

              _outputFileChooser = new JFileChooser();
              _outputFileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
              _outputFileChooser.setMultiSelectionEnabled(false);
              _outputFileChooser.setFileFilter(
                     new javax.swing.filechooser.FileFilter() {
                        public boolean accept(File f)
                           {return f.toString().endsWith(".tsv") || f.isDirectory();}
                        public String getDescription()
                           {return "TSV files";}
                      }
              );

               this.add(contentPanel);
               setResizable(true);
         }
         catch (Exception x)
         {
             ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), x);
             throw new RuntimeException(x);
         }

         initStuff();
     }


      private void transitionListInitializations() {
          listTransition = new JList();
          listTransition.setCellRenderer(new coloredMRMListRenderer());
          listTransition.setModel(new DefaultListModel());
          listTransition.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
          _tlsl = new transitionListSelectionListener();
          listTransition.addListSelectionListener(_tlsl);
          for (MRMTransition curTran : _mrmTransitions)
              ((DefaultListModel)listTransition.getModel()).addElement(curTran);
          listTransition.setVisible(true);
          listScrollPane.getViewport().setView(listTransition);
          listScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
          listScrollPane.getVerticalScrollBar().addAdjustmentListener(
            new AdjustmentListener() {
                public void adjustmentValueChanged(AdjustmentEvent ae) {
                   if(!ae.getValueIsAdjusting()) {
                       listTransition.repaint();
                   }
                }
            }
          );
      }

//ArrayList<ReactionClusterable> allDaughters = new ArrayList<ReactionClusterable>();
//Clusterer2D reactionClusters = null;

      private void peaksDataInitializations() {
          String classElements[] = _ecurveclass.getName().split("\\.");
          elutionTableLabel.setText("<html><body><center>Elution Data<br><font size='-1'>("+classElements[classElements.length-1]+")</font></center></body></html>");
          elutionTableLabel.setHorizontalAlignment(JLabel.CENTER);
          elutionTableLabel.setHorizontalTextPosition(JLabel.CENTER);

          if(transDefHeader == null || transDefHeader.getAQUApairs() == null || transDefHeader.getAQUApairs().size() == 0){
             buttonFindMate.setVisible(false);
          } else {
             buttonFindMate.setVisible(true);
          }
          peaksScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
          peaksScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
          peaksTable = new JTable(new PeaksTableModel() );
          //peaksTable.setPreferredScrollableViewportSize(new Dimension(500, 700));
          peaksTable.setSelectionModel(new peaksTableSelectionModel());
          peaksTable.setAutoscrolls(true);
          peaksTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

          int totalTableRows = 1;
          for(MRMTransition curTrans: _mrmTransitions) totalTableRows += (1+ curTrans.getDaughters().size());
          ((PeaksTableModel)(peaksTable.getModel())).data = new Object[totalTableRows-1][peaksData.values().length];
          for(peaksData pd: EnumSet.allOf(peaksData.class)) {
             peaksTable.getColumnModel().getColumn(pd.colno).setPreferredWidth(pd.colWidth);
          }
          peaksTable.doLayout();
          ((DefaultCellEditor)peaksTable.getDefaultEditor(peaksData.Accept.colClass)).setClickCountToStart(1);
          int i = 0;
          for(MRMTransition curTrans: _mrmTransitions) {
             curTrans.setGraphData(makeParentSeries(curTrans));
             int curPrecursorIndex = i;
             curTrans.setTableRow(curPrecursorIndex);
             for(peaksData pd: EnumSet.allOf(peaksData.class)) {
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][pd.colno] = null;
                 pd.makeVisible(true);
             }

             ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Accept.colno] = null;
             ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Precursor.colno] = curTrans;
             for(MRMDaughter d: curTrans.getDaughters().values()) {

                 i++;
                 d.setGraphData(d.makeDaughterSeries());
                 d.setContinDaughterData(d.makeDaughterSeries(d,true));
                 d.setElutionDataTableRow(i);
                 ElutionCurveStrategy bes = ElutionCurveStrategy.getInstance(curTrans,d,_ecurveclass);
                 bes.calculateParentElutionCurves(null);
                 bes.calculateDaughterElutionCurves(null);
                 bes.calculateBestCurves();
                 d.calculateQuality();
                 curTrans.getElutionCurves().put(d, bes);

                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Accept.colno] = new Boolean(!Utils.allYsAre0(d));
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Peptide.colno] ="";
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Precursor.colno] = curTrans;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Daughter.colno] = d;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.CoStart.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.CoEnd.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.CoDelta.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.AUC.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.MaxPeak.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Quality.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Label.colno] = "";
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Code.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.LHRatio.colno] = null;
                 ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Comment.colno] = "";
                 if(transDefHeader != null && transDefHeader.getDToTD() != null &&  transDefHeader.getDToTD().get(d) != null) {
                     TransitionDefinition td = transDefHeader.getDToTD().get(d);
                     ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Peptide.colno] =td.getPeptide();
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.Peptide.colno] =td.getPeptide();
                     ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Code.colno] = new Integer(td.getAQUAcode());
                 }
                 ElutionCurve bestPrecursorCurve = bes.getBestParentCurve();
                 if(bestPrecursorCurve == null  || bestPrecursorCurve.getMinElutionTimeSecs() <= 0.0){
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.AUC.colno] = new Float(-1);
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.MaxPeak.colno] = new Float(-1);
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.Quality.colno] = new Float(-1);
                 }  else {
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.AUC.colno] = new Float(bestPrecursorCurve.getAUC());
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.MaxPeak.colno] = new Float(bestPrecursorCurve.getHighestPointY());
                     ((PeaksTableModel)(peaksTable.getModel())).data[curPrecursorIndex][peaksData.Quality.colno] = new Float(curTrans.getQuality());
                 }

                 ElutionCurve bestDaughterCurve = bes.getBestDaughterCurve();
                 if(bestDaughterCurve == null  || bestDaughterCurve.getMinElutionTimeSecs() <= 0.0) {
                    ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Accept.colno] = new Boolean(false);
                    ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.AUC.colno] = new Float(-1);
                    ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.MaxPeak.colno] = new Float(-1);
                    ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Quality.colno] = new Float(-1);
                 }  else {
                     ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.AUC.colno] = new Float(bestDaughterCurve.getAUC());
                     ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.MaxPeak.colno] = new Float(bestDaughterCurve.getHighestPointY());
                     d.setBestElutionCurve(bestDaughterCurve);
                     ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Quality.colno] = new Float(d.getQuality());
                     if(_minPeakCutoff > 0 && bestDaughterCurve.getHighestPointY() < _minPeakCutoff) ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Accept.colno] = new Boolean(false);
                     if(_minAreaCutoff > 0 && bestDaughterCurve.getAUC() < _minAreaCutoff) ((PeaksTableModel)(peaksTable.getModel())).data[i][peaksData.Accept.colno] = new Boolean(false);
                 }
              }
              curTrans.setElutionRegionStart(curTrans.calculateMinOfAllBestDaughterCurves());
              curTrans.setElutionRegionEnd(curTrans.calculateMaxOfAllBestDaughterCurves()); 
              for(int j=curPrecursorIndex; j<=i; j++){
                ((PeaksTableModel)(peaksTable.getModel())).data[j][peaksData.CoStart.colno] = new Float(curTrans.getElutionRegionStart());
                ((PeaksTableModel)(peaksTable.getModel())).data[j][peaksData.CoEnd.colno] = new Float(curTrans.getElutionRegionEnd());
                ((PeaksTableModel)(peaksTable.getModel())).data[j][peaksData.CoDelta.colno] = new Float(curTrans.getElutionRegionEnd()-curTrans.getElutionRegionStart());
              }
              i++;
          }

          peaksTable.setDefaultRenderer(MRMTransition.class, new MRMTransitionTableRenderer(false));
          peaksTable.setDefaultRenderer(MRMDaughter.class, new MRMDaughterTableRenderer(false));
          peaksTable.setDefaultRenderer(Number.class, new MRMNumberTableRenderer());
          peaksTable.setDefaultRenderer(Integer.class, new MRMNumberTableRenderer());
          peaksTable.setDefaultRenderer(Boolean.class, new MRMBooleanRenderer());
          peaksTable.getColumnModel().getColumn(peaksData.CoStart.colno).setCellEditor(new NumberTableCellEditor());
          peaksTable.getColumnModel().getColumn(peaksData.CoEnd.colno).setCellEditor(new NumberTableCellEditor());
          peaksTable.getColumnModel().getColumn(peaksData.Code.colno).setCellEditor(new NumberTableCellEditor());
          peaksTable.getColumnModel().getColumn(peaksData.LHRatio.colno).setCellEditor(new NumberTableCellEditor());

          peaksScrollPane.getViewport().setView(peaksTable);

          if(transDefHeader == null) {
              peaksData.Peptide.makeVisible(false);
              peaksData.Label.makeVisible(false);
              peaksData.LHRatio.makeVisible(false);
              peaksData.Code.makeVisible(false);
          } else {
              if(transDefHeader.getAQUApairs() == null || transDefHeader.getAQUApairs().isEmpty()) {
                  peaksData.Label.makeVisible(false);
                  peaksData.LHRatio.makeVisible(false);
                  peaksData.Code.makeVisible(false);
              }
          }

//  "Quality" column, currently unused, is invisible unless one or more of its
//  values is not -1          
          peaksData.Quality.makeVisible(!Utils.qualColIsEmpty());

          peaksScrollPane.getVerticalScrollBar().addAdjustmentListener(
            new AdjustmentListener() {
                public void adjustmentValueChanged(AdjustmentEvent ae) {
                   if(!ae.getValueIsAdjusting()) {
                       peaksTable.repaint();
                   }
                }
            }
          );
          peaksTable.getModel().addTableModelListener(new peaksTableListener());
          _ptmlsl = new PeaksTableListSelectionListener();
          peaksTable.getSelectionModel().addListSelectionListener(_ptmlsl);

//reactionClusters = new Clusterer2D(allDaughters.toArray(new ReactionClusterable[1]));
//reactionClusters.setDimensionSplitCalculator(
//        new Clusterer2D.ClusterDimensionSplitCalculator()
//        {
//             public double calculateDimension1ForSplit(double referenceValue, double dimensionValue)
//                {
//                    return (double) _precursorDiscoveryMzTolerance;
//                }
//
//                public double calculateDimension2ForSplit(double referenceValue, double dimensionValue)
//                {
//                    return _daughterMzTolerance;
//                }
//        }
//);
//Pair<Double,Double> bestBucketSizes = reactionClusters.calculateBestBuckets(parentBuckets, daughterBuckets);
//          clusterer.split2D(bestBucketSizes.first, bestBucketSizes.second);
//reactionClusters.split2D(_precursorDiscoveryMzTolerance,_daughterMzTolerance);
//Clusterer2D.BucketSummary buckSums[] = reactionClusters.summarize();
//for(int iii = 0; iii<buckSums.length; iii++)System.out.println(buckSums[iii]);
      }

      public void AQUAinitializations() {
         if(transDefHeader == null || transDefHeader.getAQUApairs() == null) return;
         for(TransitionDefinitionHeader.AQUApair aqp:  transDefHeader.getAQUApairs().values()) {
             if(aqp.getHeavyMember().getAssociatedProduct() == null || aqp.getLightMember().getAssociatedProduct() == null) continue;
             MRMDaughter light = aqp.getLightMember().getAssociatedProduct();
             MRMDaughter heavy = aqp.getHeavyMember().getAssociatedProduct();
             if(light.getBestElutionCurve() == null || heavy.getBestElutionCurve() == null) continue;
             ElutionCurve lec = light.getBestElutionCurve();
             ElutionCurve hec = heavy.getBestElutionCurve();
             Float lhratio = new Float(lec.getAUC()/hec.getAUC());
             ((PeaksTableModel)(peaksTable.getModel())).setValueAt(new Character(aqp.getLightMember().getLowOrHigh()).toString(),light.getElutionDataTableRow(),peaksData.Label.colno);
             ((PeaksTableModel)(peaksTable.getModel())).setValueAt(new Character(aqp.getHeavyMember().getLowOrHigh()).toString(),heavy.getElutionDataTableRow(),peaksData.Label.colno);
             ((PeaksTableModel)(peaksTable.getModel())).setValueAt(lhratio,light.getElutionDataTableRow(),peaksData.LHRatio.colno);
             ((PeaksTableModel)(peaksTable.getModel())).setValueAt(lhratio,heavy.getElutionDataTableRow(),peaksData.LHRatio.colno);
         }      
      }

      public void initStuff()
      {
         transDefHeader = null;
         transitionOnPlot = null;

         try {
            _run = MSRun.load(_mzXMLFile.getAbsolutePath());
            _mrmTransitions = loadMRMTransitions(_run);
         } catch (Exception e) {
             ApplicationContext.errorMessage(TextProvider.getText("ERROR_CREATING_DIALOG"), e);
         }

         String _transitionDefFilePath = _mzXMLFile.getAbsolutePath().replaceAll("\\.mzXML$",".transition.tsv");
         if((new File(_transitionDefFilePath)).exists()){
            transDefHeader = new TransitionDefinitionHeader(_transitionDefFilePath,new TSVTransitionDefinitionParser());
            try {
               transDefHeader.getParser().setTransitionDefFile(_transitionDefFilePath);
               transDefHeader.doParse();
               if(transDefHeader==null || transDefHeader.getTransitionDefs() == null || transDefHeader.getTransitionDefs().size()==0) {
                  transDefHeader = null;
               } else {
                  transDefHeader.linkUpToTransitionList(_mrmTransitions,_precursorDiscoveryMzTolerance,_daughterMzTolerance);
                  transDefHeader.determinePairs();
               }
            } catch (Exception e) {
               System.err.println(e);
            }
         }
         if(transDefHeader != null && transDefHeader.getComment() != null && !transDefHeader.getComment().equals("")) {
             fileNameLabel.setText("<html><body><center>"+_mzXMLFile.getName()+"<br>"+transDefHeader.getComment()+"</center></body></html>");
             fileNameLabel.setHorizontalAlignment(JLabel.CENTER);
             fileNameLabel.setHorizontalTextPosition(JLabel.CENTER);
             fileNameLabel.setMinimumSize(new Dimension(1200,30));
         }
         else {
             fileNameLabel.setText(_mzXMLFile.getName());
         }
         if(transDefHeader != null && transDefHeader.getMzXMLFile() != null && !transDefHeader.getMzXMLFile().equals("")) {
            if(!transDefHeader.getMzXMLFile().trim().equals(_mzXMLFile.getName().trim())) {
                ApplicationContext.infoMessage("Note: current mzXML ("+_mzXMLFile.getName()+") does not match mzXML file name in transition.tsv file ("+transDefHeader.getMzXMLFile()+")");
            }
         }
         helper.addListener(this, "contentPanel_componentResized");

         transitionListInitializations();
         peaksDataInitializations();
         AQUAinitializations();
         String mzToleranceText =  "" + (double) _precursorChromatogramWindow;
         //this is kind of arbitrary: one digit after decimal.  It also doesn't support
         //comma-radix.  Whatever.
         if (mzToleranceText.contains("."))
             mzToleranceText = mzToleranceText.substring(0,mzToleranceText.indexOf(".") + 2);

         listTransition.setSelectedIndex(0);

         //updateChartsAndFields(false);
         peaksScrollPane.updateUI();

         listScrollPane.updateUI();

         pack();

         //debug
//         System.out.println(Runtime.getRuntime().freeMemory()+" bytes available before freeing _run");
         _run.close();
         _run = null;
         for(MRMTransition mrmt: _mrmTransitions ) mrmt.setRun(null);
         System.gc();
//         System.out.println(Runtime.getRuntime().freeMemory()+" bytes available after freeing _run");

     }

  // Listeners
      private volatile boolean _kludge_listener_removal = false;
      public class transitionListSelectionListener implements ListSelectionListener {
          public void valueChanged(ListSelectionEvent e)
          {
              if(e.getValueIsAdjusting()) return;
              if(_kludge_listener_removal) return;
              MRMTransition curTrans = null;
              try
              {
                 listTransition.getSelectionModel().removeListSelectionListener(this);
                 int indices[]=listTransition.getSelectedIndices();
                 DefaultListModel dlm = (DefaultListModel) listTransition.getModel();
                 for(int i=0; i<indices.length;i++) {
                   curTrans = (MRMTransition)dlm.get(indices[i]);
                   curTrans.setCurrentDaughterIndex(0);
                   int dataIndex = -1;
                   Object[][] tableData = ((PeaksTableModel)(peaksTable.getModel())).data;
                   for(int j = 0; j<tableData.length; j++)
                   {
                       if(curTrans == (MRMTransition)(tableData[j][peaksData.Precursor.colno]))
                       {
                           dataIndex = j;
                       }
                   }
                   if(dataIndex > -1)
                   {
                       ListSelectionListener curlsl = null;
                       ListSelectionListener larr[] = ((peaksTableSelectionModel)peaksTable.getSelectionModel()).getListSelectionListeners();
                       for(int ii = 0; ii<larr.length; ii++) {
                           if(larr[ii] instanceof PeaksTableListSelectionListener) {
                              curlsl = larr[ii];
                              peaksTable.getSelectionModel().removeListSelectionListener(curlsl);
                              //break;
                           }
                       }
                       scrollPeakTableToRow(curTrans.getCurrentDaughter().getElutionDataTableRow());
                       if(curlsl != null) peaksTable.getSelectionModel().addListSelectionListener(curlsl);
                   }
                 }
                 transitionOnPlot = curTrans;
                 updateChartsAndFields(false);
                 listTransition.getSelectionModel().addListSelectionListener(this);

              } catch (Exception ee) {
                 ApplicationContext.infoMessage("Can't plot selected index: "+ee);
              }
          }
      }

      public class PeaksTableListSelectionListener implements ListSelectionListener {
         public void valueChanged(ListSelectionEvent event) {
            if(event.getValueIsAdjusting())return;
            try{
                   ((ListSelectionModel)(event.getSource())).removeListSelectionListener(this);
                   // Switch transitions table to appropriate transition
                   int oldIndex = transitionOnPlot.getCurrentDaughter().getElutionDataTableRow();
                   int index = -1;
                   if(event.getLastIndex()!=oldIndex) {
                       index = event.getLastIndex();
                   } else {
                       if(event.getFirstIndex() != oldIndex) {
                           index = event.getFirstIndex();
                       } else {
                           ApplicationContext.errorMessage("Problem in PeaksTableListener: can't find correct row index", null);
                           ((ListSelectionModel)(event.getSource())).addListSelectionListener(this);
                           return;
                       }
                   }
                   // The transition we selected
                   MRMTransition curTrans =(MRMTransition)((PeaksTableModel)(peaksTable.getModel())).data[index][peaksData.Precursor.colno];

                   //Get currently selected daughter.  If we aren't positioned on a row
                   //that contains a daughter, assume the first daughter of the transition
                   MRMDaughter curD =(MRMDaughter)((PeaksTableModel)(peaksTable.getModel())).data[index][peaksData.Daughter.colno];
                   if(curD == null) curD = curTrans.getDaughters().values().iterator().next();

                   //Set transition's current daughter to the daughter on the table row
                   Object dArr[] = curTrans.getDaughters().values().toArray();
                   for(int i = 0; i<dArr.length; i++) {
                       if(dArr[i] == (Object) curD) {
                           curTrans.setCurrentDaughterIndex(i);
                           break;
                       }
                   }

                   //Move the transition list (left table on screen) to the correct index
                   //if we have moved from one transition to another
                   //Make sure to remove and replace listener, so it won't do any
                   //mischief.  Then replace it afterward *only if it had been there before*!

                   //find listTransition listener of the proper variety, if there is one
                   ListSelectionListener curlsl = null;
                   ListSelectionListener larr[] = listTransition.getListSelectionListeners();
                   for(int i = 0; i<larr.length; i++) {
                       if(larr[i] instanceof transitionListSelectionListener) {
                          curlsl = larr[i];
                          listTransition.removeListSelectionListener(curlsl);
                          //break;
                       }
                   }
                   //todo:  find out why removal of transitionListSelectionListener isn't working
                   _kludge_listener_removal = true;
                   if(curTrans != listTransition.getSelectedValues()[0]){
                     listTransition.setSelectedValue(curTrans,true);
                   }
                   _kludge_listener_removal = false;
                   //replace listener if it had been there before
                   if(transDefHeader != null) {
                       buttonFindMate.setEnabled(transDefHeader.getMatchingDaughter(curD) != null);
                   }
                   if(curlsl != null) listTransition.addListSelectionListener(curlsl);

                   // redisplay graph
                   transitionOnPlot = curTrans;
                   updateChartsAndFields(false);
                   //put current selection model back
                   ((ListSelectionModel)(event.getSource())).addListSelectionListener(this);

            } catch (Exception e) {
                ApplicationContext.infoMessage("Can't set MRMer to selected row: "+e+"\nPlease contact maintainers.");
                e.printStackTrace();
                ((ListSelectionModel)(event.getSource())).addListSelectionListener(this);
            }
         }
     }

    public class domainAxisZoomCoordinator implements AxisChangeListener
    {
        public CenterZoomNumberAxis getMyAxis() {
            return myAxis;
        }

        public void setMyAxis(CenterZoomNumberAxis myAxis) {
            this.myAxis = myAxis;
        }

        CenterZoomNumberAxis myAxis;

        public domainAxisZoomCoordinator(CenterZoomNumberAxis na)
        {
            myAxis = na;
        }

        public void axisChanged(AxisChangeEvent ace)
        {
            if(precursorChromatogramEmpty(transitionOnPlot)) return;
            CenterZoomNumberAxis theOtherAxis = null;
            domainAxisZoomCoordinator theOtherCoordinator = null;
            if(getMyAxis() == findAxisOfPanel(precursorContainerPanel,whichAxis.Domain)){
                theOtherAxis = (CenterZoomNumberAxis) findAxisOfPanel(daughterContainerPanel, whichAxis.Domain);
            } else {
                theOtherAxis = (CenterZoomNumberAxis) findAxisOfPanel(precursorContainerPanel,whichAxis.Domain);
            }
            theOtherCoordinator = theOtherAxis.getTheListener();
            theOtherAxis.removeChangeListener(theOtherCoordinator);
            getMyAxis().removeChangeListener(this);
            theOtherAxis.setLowerBound(getMyAxis().getLowerBound());
            theOtherAxis.setUpperBound(getMyAxis().getUpperBound());
            getMyAxis().addChangeListener(this);
            theOtherAxis.addChangeListener(theOtherCoordinator);
//            theOtherCoordinator.axisChanged(new AxisChangeEvent(theOtherAxis));
        }
    }

    public class daughterRangeAxisChangeListener implements AxisChangeListener {
        public void axisChanged(AxisChangeEvent ace) {
            XYPlot xyp = (XYPlot) ace.getAxis().getPlot();
            if(xyp != null && ace.getAxis() == xyp.getRangeAxis())
            {
               ValueAxis precursorRangeAxis = (ValueAxis)findAxisOfPanel(precursorContainerPanel,whichAxis.Range);
               ValueAxis daughterRangeAxis = (ValueAxis) ace.getAxis();
               precursorRangeAxis.setLowerBound(daughterRangeAxis.getLowerBound());
               precursorRangeAxis.setUpperBound(daughterRangeAxis.getUpperBound());
            }
        }
    }

    public class peaksTableSelectionModel extends DefaultListSelectionModel {

        public ListSelectionListener[] getListSelectionListeners() {
            return (ListSelectionListener[])listenerList.getListeners(
                    ListSelectionListener.class);
        }
    }

    public class peaksTableListener implements TableModelListener {
       public void tableChanged(TableModelEvent e) {
           int row = e.getFirstRow();
           int column = e.getColumn();
           TableModel model = (TableModel)e.getSource();
           String columnName = model.getColumnName(column);

           if(!columnName.equalsIgnoreCase("Start") &&
              !columnName.equalsIgnoreCase("End")   &&
              !columnName.equalsIgnoreCase("Accept")
           ) return;
           Object data = model.getValueAt(row, column);
           if(columnName.equalsIgnoreCase("Accept")){
              if(data != null) {
                 boolean accepted = (Boolean) data;
                 if(!accepted) {
                    //model.setValueAt(new Float(-1),row,peaksData.AUC.colno);
                 } else {
                    MRMDaughter curD = (MRMDaughter)model.getValueAt(row,peaksData.Daughter.colno);
                    ElutionCurve curEC = null;
                    if(curD != null) curEC = curD.getBestElutionCurve();
                    if(curEC != null) {
                       model.setValueAt(
                          new Float(curEC.getAUC()),
                          row,
                          peaksData.AUC.colno
                       );
                    }
                 }
              }
              updateDaughterCharts(false); 
              return; 
           }

           float colVal = Float.parseFloat(data.toString());
           float otherVal = 0.0f;
           float start = 0f;
           float end = 0f;
           // Is this a daughter or parent?
           // assume daughter first
           MRMTransition mrt = null;
           MRMDaughter mrd = (MRMDaughter)model.getValueAt(row,peaksData.Daughter.colno);
           if(mrd == null){
               mrt = (MRMTransition)model.getValueAt(row,peaksData.Precursor.colno);
           }  else {
               mrt = mrd.getPrecursor();
           }          

           if(columnName.equalsIgnoreCase("Start")) {
               start = colVal;
               otherVal = Float.parseFloat(model.getValueAt(row,peaksData.CoEnd.colno).toString());
               end = otherVal;
               if(colVal >= otherVal) {
                   JOptionPane.showMessageDialog(null,"Starting elution time must be less than ending elution time","Elution Range Error",JOptionPane.ERROR_MESSAGE);
                   model.removeTableModelListener(this);
                   model.setValueAt(new Float(0.0f),row, column);
                   model.addTableModelListener(this);
                   return;
               }
               mrt.setElutionRegionStart(start);
           } else {
               end = colVal;
               otherVal = Float.parseFloat(model.getValueAt(row,peaksData.CoStart.colno).toString());
               start = otherVal;
               if(colVal <= otherVal){
                   JOptionPane.showMessageDialog(null,"End elution time must be greater than starting elution time","Elution Range Error",JOptionPane.ERROR_MESSAGE);
                   model.removeTableModelListener(this);
                   model.setValueAt(new Float(100000000.0f),row, column);
                   model.addTableModelListener(this);
                   return;
               }
               mrt.setElutionRegionEnd(end);
           }
           // elution boundaries have been changed manually
           // change all daughter columns
           model.removeTableModelListener(this);
           model.setValueAt(new Float(start),mrt.getTableRow(),peaksData.CoStart.colno);
           model.setValueAt(new Float(end),mrt.getTableRow(),peaksData.CoEnd.colno);
           float delta = end-start;
           model.setValueAt(new Float(delta),mrt.getTableRow(),peaksData.CoDelta.colno);
           for(MRMDaughter mrmd: mrt.getDaughters().values()){
              model.setValueAt(new Float(start),mrmd.getElutionDataTableRow(),peaksData.CoStart.colno);
              model.setValueAt(new Float(end),mrmd.getElutionDataTableRow(),peaksData.CoEnd.colno);
              model.setValueAt(new Float(delta),mrmd.getElutionDataTableRow(),peaksData.CoDelta.colno);
           }
           Vector<MRMDaughter> matchingDaughters = null;
           if(transDefHeader != null && _synclh) {
              matchingDaughters = new Vector<MRMDaughter>();
              for(MRMDaughter mrmd: mrt.getDaughters().values()){
                 MRMDaughter matchingDaughter = transDefHeader.getMatchingDaughter(mrmd);
                 if(matchingDaughter != null) {
                    model.setValueAt(new Float(start),matchingDaughter.getElutionDataTableRow(),peaksData.CoStart.colno);
                    model.setValueAt(new Float(end),matchingDaughter.getElutionDataTableRow(),peaksData.CoEnd.colno);
                    model.setValueAt(new Float(delta),matchingDaughter.getElutionDataTableRow(),peaksData.CoDelta.colno);
                    matchingDaughters.add(matchingDaughter);
                 }
              }
           }
           // recalculate all curves??
           //     or reselect region as elution curve?
           //     seems like a dumb idea
           mrt.getElutionCurves().clear();
           for(MRMDaughter d: mrt.getDaughters().values()){
              d.setBestElutionCurve(null);
              ElutionCurveStrategy ecs =
                 ElutionCurveStrategy.getInstance(mrt,d,_ecurveclass);
              ecs.calculateParentElutionCurves(null);
              ecs.calculateDaughterElutionCurves(null);
              ecs.calculateBestCurves();
              d.setBestElutionCurve(ecs.getBestDaughterCurve());
              mrt.getElutionCurves().put(d,ecs);
              ElutionCurve bestPrecursorCurve = ecs.getBestParentCurve();
              ElutionCurve bestDaughterCurve = ecs.getBestDaughterCurve();
              if(bestPrecursorCurve == null ||
                 bestPrecursorCurve.getMinElutionTimeSecs()<=0.0) {
                 model.setValueAt(new Float(-1),mrt.getTableRow(),peaksData.AUC.colno);
              } else {
                  model.setValueAt(new Float(bestPrecursorCurve.getAUC()),mrt.getTableRow(),peaksData.AUC.colno);
              }
              if(bestDaughterCurve == null ||
                 bestDaughterCurve.getMinElutionTimeSecs()<=0.0
              ) {
                    model.setValueAt(new Float(-1),d.getElutionDataTableRow(),peaksData.AUC.colno);
              } else {
                    model.setValueAt(new Float(bestDaughterCurve.getAUC()),d.getElutionDataTableRow(),peaksData.AUC.colno);
              }
           }
           if(matchingDaughters != null && !matchingDaughters.isEmpty()) {
               MRMTransition otherMRT = matchingDaughters.get(0).getPrecursor();
               otherMRT.getElutionCurves().clear();
               model.setValueAt(new Float(start),otherMRT.getTableRow(),peaksData.CoStart.colno);
               model.setValueAt(new Float(end),otherMRT.getTableRow(),peaksData.CoEnd.colno);
               model.setValueAt(new Float(delta),otherMRT.getTableRow(),peaksData.CoDelta.colno);
               otherMRT.setElutionRegionStart(mrt.getElutionRegionStart());
               otherMRT.setElutionRegionEnd(mrt.getElutionRegionEnd());
               for(MRMDaughter md: matchingDaughters) {
                  md.setBestElutionCurve(null);
                  ElutionCurveStrategy mecs = ElutionCurveStrategy.getInstance(otherMRT,md,_ecurveclass);
                  mecs.calculateParentElutionCurves(null);
                  mecs.calculateDaughterElutionCurves(null);
                  mecs.calculateBestCurves();
                  md.setBestElutionCurve(mecs.getBestDaughterCurve());
                  otherMRT.getElutionCurves().put(md,mecs);
                  ElutionCurve bestPrecursorCurve = mecs.getBestParentCurve();
                  ElutionCurve bestDaughterCurve = mecs.getBestDaughterCurve();
                  if(bestPrecursorCurve == null || bestPrecursorCurve.getMinElutionTimeSecs()<=0.0) {
                     model.setValueAt(new Float(-1),otherMRT.getTableRow(),peaksData.AUC.colno);
                  } else {
                     model.setValueAt(new Float(bestPrecursorCurve.getAUC()),otherMRT.getTableRow(),peaksData.AUC.colno);
                  }
                  if(bestDaughterCurve == null || bestDaughterCurve.getMinElutionTimeSecs()<=0.0) {
                     model.setValueAt(new Float(-1),md.getElutionDataTableRow(),peaksData.AUC.colno);
                  } else {
                     model.setValueAt(new Float(bestDaughterCurve.getAUC()),md.getElutionDataTableRow(),peaksData.AUC.colno);
                  }
               }
           }
           // redisplay AUCs
           if(mrt == transitionOnPlot){
               updateChartsAndFields(new defaultGraphZone());
           }
           if(transDefHeader != null) AQUAinitializations();
           model.addTableModelListener(this);
       }
    }

    public void buttonPostPressTasks() {
       listTransition.requestFocus();
    }

    public void buttonSICUpdate_actionPerformed(ActionEvent event)
    {
        try
        {
            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);
            String answer = nf.format(_precursorChromatogramWindow);
            String newAnswer = JOptionPane.showInputDialog("New value for MS1 SIC tolerance:",answer);
            _precursorChromatogramWindow = (float) Double.parseDouble(newAnswer);

            for(MRMTransition t : _mrmTransitions)
            {
               t.setGraphData(null);
               for(MRMDaughter d : t.getDaughters().values())
               {
                  d.setGraphData(null);
               }
            }
            updateChartsAndFields(false);
        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Bad value for MZ tolerance, please try again");
        }
    }

    public void buttonPDT_actionPerformed(ActionEvent event)
     {
         try
         {
             NumberFormat nf = NumberFormat.getNumberInstance();
             nf.setMaximumFractionDigits(2);
             nf.setMinimumFractionDigits(2);
             String answer = nf.format(_precursorDiscoveryMzTolerance);
             String newAnswer = JOptionPane.showInputDialog("New value for precursor tolerance:",answer);
             _precursorDiscoveryMzTolerance = (float) Double.parseDouble(newAnswer);
             initStuff();
         }
         catch (Exception e)
         {
             ApplicationContext.infoMessage("Bad value for precursor tolerance.");
         }
     }

    public void buttonDTOL_actionPerformed(ActionEvent event)
        {
            try
            {
                NumberFormat nf = NumberFormat.getNumberInstance();
                nf.setMaximumFractionDigits(4);
                nf.setMinimumFractionDigits(4);
                String answer = nf.format(_daughterMzTolerance);
                String newAnswer = JOptionPane.showInputDialog("New value for product tolerance:",answer);
                _daughterMzTolerance = (float) Double.parseDouble(newAnswer);
                initStuff();
            }
            catch (Exception e)
            {
                ApplicationContext.infoMessage("Bad value for product tolerance.");
            }
        }

    public void menuItemPMin_actionPerformed(ActionEvent event) {
        try {
            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);
            nf.setGroupingUsed(false);
            String answer = nf.format(_minPeakCutoff);
            String newAnswer = JOptionPane.showInputDialog("New value for min peak cutoff:",answer);
            if(newAnswer == null) return;
            float newCutoff = Float.parseFloat(newAnswer);
            if(newCutoff > 0.0f) {
               for(MRMTransition mrt: _mrmTransitions) {
                   for(MRMDaughter curd: mrt.getDaughters().values()) {
                       if(curd.getBestElutionCurve() != null && curd.getBestElutionCurve().getHighestPointY() < newCutoff) {
                           ((PeaksTableModel)(peaksTable.getModel())).setValueAt(new Boolean(false),curd.getElutionDataTableRow(),peaksData.Accept.colno);
                       }
                   }
               }
            }
            _minPeakCutoff = newCutoff;
        } catch (Exception e) {
            ApplicationContext.infoMessage("Failed to change min acceptable peak height: "+e);
        }
    }

    public void menuItemAMin_actionPerformed(ActionEvent event) {
        try {
            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);
            nf.setGroupingUsed(false);
            String answer = nf.format(_minAreaCutoff);
            String newAnswer = JOptionPane.showInputDialog("New value for min AUC cutoff:",answer);
            if(newAnswer == null) return;
            float newCutoff = Float.parseFloat(newAnswer);
            if(newCutoff > 0.0f) {
               for(MRMTransition mrt: _mrmTransitions) {
                   for(MRMDaughter curd: mrt.getDaughters().values()) {
                       if(curd.getBestElutionCurve() != null && curd.getBestElutionCurve().getAUC() < newCutoff) {
                           ((PeaksTableModel)(peaksTable.getModel())).setValueAt(new Boolean(false),curd.getElutionDataTableRow(),peaksData.Accept.colno);
                       }
                   }
               }
            }
            _minAreaCutoff = newCutoff;
        } catch (Exception e) {
            ApplicationContext.infoMessage("Failed to change min acceptable AUC: "+e);
        }
    }

    public void menuItemQuit_actionPerformed(ActionEvent event)
    {
       try {
           System.exit(1);
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot quit: "+e);
       }
    }

    public void menuItemTraceAllFragments_actionPerformed(ActionEvent event) {
       try {
          int newTraceAllCode = JOptionPane.showConfirmDialog(this,"Trace over elution curves of ALL fragments?","Trace All Fragments",JOptionPane.YES_NO_CANCEL_OPTION);
          if(newTraceAllCode != JOptionPane.YES_OPTION && newTraceAllCode != JOptionPane.NO_OPTION) return;
          boolean newTraceAll = (newTraceAllCode == JOptionPane.YES_OPTION);
          if(newTraceAll != _traceAllFragments) {
             _traceAllFragments = newTraceAll;
             updateChartsAndFields(false);
          }
          _traceAllFragments = newTraceAll;
       }  catch (Exception e) {
           ApplicationContext.infoMessage("Cannot change 'Trace All' mode: "+e);
       }
    }
    
    public void menuItemSyncLH_actionPerformed(ActionEvent event) {
       try {
          int newsync = JOptionPane.showConfirmDialog(this,"Synchronize L/H Elution Regions","Sync Label Elution Regions",JOptionPane.YES_NO_CANCEL_OPTION);
          if(newsync != JOptionPane.YES_OPTION && newsync != JOptionPane.NO_OPTION) return;
          boolean newsyncbool = (newsync == JOptionPane.YES_OPTION);
          _synclh = newsyncbool;
       }  catch (Exception e) {
           ApplicationContext.infoMessage("Cannot change 'Synchronize L/H' mode: "+e);
       }
    }

    public void menuItemChangeStrategy_actionPerformed(ActionEvent event){
       try {
          String newStrategyName =
             (String) JOptionPane.showInputDialog(
                this,
                "Which elution curve strategy would you like to use?",
                "Elution Strategy",
                JOptionPane.PLAIN_MESSAGE,
                null,
                MRMCommandLineModule.strategies,
                MRMCommandLineModule.strategies[0]
            );
          newStrategyName = "org.fhcrc.cpl.viewer.mrm." + newStrategyName;
          if(newStrategyName != _ecurveclass.getName()) {
              Class newClass = Class.forName(newStrategyName);
              _ecurveclass = newClass;
              initStuff();
          }
       } catch (Exception e) {
           ApplicationContext.infoMessage("Cannot change elution curve strategy: "+e);
       }
    }

    public void menuItemOpen_actionPerformed(ActionEvent event)
    {
       try {

          int returnVal = _mzXMLFileChooser.showOpenDialog(this);
          if(returnVal == JFileChooser.APPROVE_OPTION) {
             _mzXMLFile =_mzXMLFileChooser.getSelectedFile();
             initStuff();
          }
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot open file '"+_mzXMLFileChooser.getName()+"': "+e);
       }
    }

    public void menuItemTips_actionPerformed(ActionEvent event)
    {
       try {
           HtmlViewerPanel.showURLInDialog("http://proteomics.fhcrc.org/CPL/mrm_user_tips.html","Tips for Using MRMer");
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot display TIPS help:"+e);
       }
    }
    public void menuItemDefinitions_actionPerformed(ActionEvent event)
    {
       try {
           HtmlViewerPanel.showURLInDialog("http://proteomics.fhcrc.org/CPL/MRMer_help.html#mrmcomopt","OPTION menu commands");
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot display Definitions help:"+e);
       }
    }
    public void menuItemOptions_actionPerformed(ActionEvent event)
    {
       try {
           HtmlViewerPanel.showURLInDialog("http://proteomics.fhcrc.org/CPL/MRMer_help.html#mrmmenopt","OPTION menu commands");
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot display OPTIONS help:"+e);
       }
    }
    public void menuItemArguments_actionPerformed(ActionEvent event)
    {
       try {
           CommandLineModule module = CommandLineModuleDiscoverer.getCommandLineModule("MRM");
           String dummyCaller = "dummy_help_caller";
           File tempHelpFile = TempFileManager.createTempFile("help_mrm", dummyCaller);
           PrintWriter outPW = new PrintWriter(tempHelpFile);
           CLMUserManualGenerator.generateCommandManualEntry(module, outPW);
           outPW.flush();
           outPW.close();
           HtmlViewerPanel.showFileInDialog(tempHelpFile,"Java command line options");
       } catch (Exception e)
       {
           ApplicationContext.infoMessage("Cannot show ARGUMENTS help:"+e);
       }
    }

    public void buttonSave_actionPerformed(ActionEvent event)
    {
        int returnVal = _outputFileChooser.showSaveDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
            PeaksTableModel ptm = (PeaksTableModel) peaksTable.getModel();
            if(!ptm.saveSaveModelAsTSV(_outputFileChooser.getSelectedFile()))
               ApplicationContext.infoMessage("Cannot save data to '"+_outputFileChooser.getName()+"'");
        } else {
            ApplicationContext.infoMessage("Cannot save data to '"+_outputFileChooser.getName()+"', returnVal="+returnVal);    
        }
        buttonPostPressTasks();
    }

    public void buttonZC_actionPerformed(ActionEvent event)
     {
         try
         {
            if(((PeaksTableModel)(peaksTable.getModel())).data[peaksTable.getSelectedRow()][peaksData.Daughter.colno] == null)return; 
            updateChartsAndFields(new zoomToCurveGraphZone());
         }
         catch (Exception e)
         {
             ApplicationContext.infoMessage("Can't zoom to current curve");
         }
         buttonPostPressTasks();
     }

    public void buttonZoom_actionPerformed(ActionEvent event) {
        try
         {
           updateChartsAndFields(new zoomedGraphZone());
         }
         catch (Exception e)
         {
             ApplicationContext.errorMessage("Can't zoom:"+e,e);
         }
        buttonPostPressTasks();
    }

    public void buttonNext_actionPerformed(ActionEvent event) {
        try {
            int newIndex = transitionOnPlot.getCurrentDaughter().getElutionDataTableRow() + 1;
            //what if we're positioned on a precursor line?
            if(((PeaksTableModel)(peaksTable.getModel())).data[peaksTable.getSelectedRow()][peaksData.Daughter.colno] == null) newIndex--;
            //what if we've moved out of the current group?
            if(((PeaksTableModel) peaksTable.getModel()).data[newIndex][peaksData.Daughter.colno] == null) {
                newIndex++;
            }
            if (newIndex < ((PeaksTableModel) peaksTable.getModel()).data.length) {
                scrollPeakTableToRow(newIndex);
            }
            //updateChartsAndFields(false);
        }
        catch (ArrayIndexOutOfBoundsException aioobe) {
           System.gc();
        }
        catch (Exception e) {
            ApplicationContext.infoMessage("Can't display NEXT plot");
        }
        buttonPostPressTasks();
    }

    public void buttonPrev_actionPerformed(ActionEvent event)
    {
        try
        {
               int newIndex = transitionOnPlot.getCurrentDaughter().getElutionDataTableRow()-1;
               if(((PeaksTableModel)peaksTable.getModel()).data[newIndex][peaksData.Daughter.colno] == null) {
                  newIndex--;
               }
               if(newIndex > 0) {
                  scrollPeakTableToRow(newIndex);
               }  else {
                  System.gc();
               }
            //updateChartsAndFields(false);
        }
        catch(Exception e)
        {
           ApplicationContext.infoMessage("Can't display PREVIOUS plot");
        }
        buttonPostPressTasks();
    }

    public void buttonReject_actionPerformed(ActionEvent event)
    {
        try
        {
           int selectedIndex = peaksTable.getSelectedRow();
           if(selectedIndex >= 0) {
            if(((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Daughter.colno] == null) return;
            ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Accept.colno] = new Boolean(false);
            ((PeaksTableModel)(peaksTable.getModel())).fireTableCellUpdated(selectedIndex,peaksData.Accept.colno);
            buttonNext.doClick();
           }
        } catch (Exception e) {
            ApplicationContext.infoMessage("Can't REJECT peak: "+e);
        }
        buttonPostPressTasks();
    }
    
    public void buttonRejectGroup_actionPerformed(ActionEvent event)
    {
        try
        {
           int selectedIndex = peaksTable.getSelectedRow();
           if(selectedIndex >= 0) {
              int maxIndex = -1;
              MRMTransition curTrans =
                 ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Daughter.colno] != null
                   ?  ((MRMDaughter)((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Daughter.colno]).getPrecursor()
                   :  ((MRMTransition)((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Precursor.colno]);

              for(MRMDaughter curd: curTrans.getDaughters().values()) {
                 int curIndex = curd.getElutionDataTableRow();
                 ((PeaksTableModel)(peaksTable.getModel())).data[curIndex][peaksData.Accept.colno] = new Boolean(false);
                 ((PeaksTableModel)(peaksTable.getModel())).fireTableCellUpdated(curIndex,peaksData.Accept.colno);
                 if(curIndex > maxIndex) maxIndex = curIndex;
              }
              peaksTable.getSelectionModel().setSelectionInterval(maxIndex,maxIndex);
              buttonNext.doClick();
           }
        } catch (Exception e) {
            ApplicationContext.infoMessage("Can't REJECT group: "+e);
        }
        buttonPostPressTasks();
    }

    public void buttonAccept_actionPerformed(ActionEvent event)
    {
        try
        {
           int selectedIndex = peaksTable.getSelectedRow();
           if(selectedIndex >= 0) {
            if(((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Daughter.colno] == null) return;
            ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Accept.colno] = new Boolean(true);
            ((PeaksTableModel)(peaksTable.getModel())).fireTableCellUpdated(selectedIndex,peaksData.Accept.colno);
            buttonNext.doClick();
           }
        } catch (Exception e) {
            ApplicationContext.infoMessage("Can't ACCEPT peak: "+e);
        }
        buttonPostPressTasks();
    }

    public void buttonBoost_actionPerformed(ActionEvent event)
    {
        try
        {
           int selectedIndex = peaksTable.getSelectedRow();
           if(selectedIndex >= 0) {
             String comVal = (String) ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Comment.colno];
             if(comVal.equals("D ")) return;
             ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Comment.colno] = comVal+"D ";
             ((PeaksTableModel)(peaksTable.getModel())).fireTableCellUpdated(selectedIndex,peaksData.Comment.colno);
             buttonAccept.doClick();
           }
        } catch (Exception e) {
            ApplicationContext.infoMessage("Can't add Dwell comment to product: "+e);
        }
        buttonPostPressTasks();
    }

    public void buttonTiming_actionPerformed(ActionEvent event)
    {
        try
        {
           int selectedIndex = peaksTable.getSelectedRow();
           if(selectedIndex >= 0) {
            String comVal = (String) ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Comment.colno];
            if(comVal.equals("T ")) return;
            ((PeaksTableModel)(peaksTable.getModel())).data[selectedIndex][peaksData.Comment.colno] = comVal + "T ";
            ((PeaksTableModel)(peaksTable.getModel())).fireTableCellUpdated(selectedIndex,peaksData.Comment.colno);
            buttonAccept.doClick();
           }
        } catch (Exception e) {
            ApplicationContext.infoMessage("Can't add TIMING comment to product: "+e);
        }
        buttonPostPressTasks();
    }

    public void buttonFindMate_actionPerformed(ActionEvent event)
    {
        try
        {
          if(transDefHeader == null || transDefHeader.getAQUApairs() == null)
             throw new Exception("No transitions defined.");
          MRMDaughter curMRMD = transitionOnPlot.getCurrentDaughter();
          if(curMRMD == null) throw new Exception("Can't find current daughter");
          TransitionDefinition curTD = transDefHeader.getDToTD().get(curMRMD);
          if(curTD == null) throw new Exception("Can't find transition definition for "+curMRMD.toString());
          TransitionDefinitionHeader.AQUApair curAQUA = transDefHeader.getAQUApairs().get(curTD.getAQUAcode());
          if(curAQUA == null) throw new Exception("Can't find AQUA or SILAC pair");
          MRMDaughter matchingDaughter = null;
          if(curTD.isHigh()){
              matchingDaughter = curAQUA.getLightMember().getAssociatedProduct();
          } else {
              if(curTD.isLow()) {
                  matchingDaughter = curAQUA.getHeavyMember().getAssociatedProduct();
              }
          }
          if(matchingDaughter == null) throw new Exception("Can't find opposite labeled product");
          scrollPeakTableToRow(matchingDaughter.getElutionDataTableRow());
        } catch (Exception e) {
            ApplicationContext.errorMessage("Can't find matching label: ",e);
        }
        buttonPostPressTasks();
    }

    public void contentPanel_componentResized(ComponentEvent event)
    {
        try
        {
            updateChartsAndFields(false);
        }
        catch (Exception e)
        {
            ApplicationContext.infoMessage("Cannot resize chart");
        }
    }

    public void contentPanel_mouseClicked(MouseEvent event) {
        if(event.getButton() == MouseEvent.BUTTON3) {
           contentPanel.setBackground(
              JColorChooser.showDialog(this,"Choose Background Color",contentPanel.getBackground())
           );
        }
    }

// Renderers

    public class coloredMRMListRenderer extends JPanel implements ListCellRenderer
    {
       public coloredMRMListRenderer() {
          super();
       }

       public Component getListCellRendererComponent(
                    JList list,
                    Object value,
                    int index,
                    boolean isSelected,
                    boolean cellHasFocus)
      {
         this.removeAll();
         this.setLayout(new BoxLayout(this,BoxLayout.Y_AXIS));
         NumberFormat nf = NumberFormat.getNumberInstance();
         nf.setMaximumFractionDigits(3);
         nf.setMinimumFractionDigits(3);

         JLabel precursor = new JLabel(((MRMTransition)value).getName()+"\u00B1"+nf.format(_precursorDiscoveryMzTolerance));
         Font tmpFont = precursor.getFont();
         tmpFont=tmpFont.deriveFont(tmpFont.getSize()+4.0F);
         tmpFont=tmpFont.deriveFont(Font.BOLD);
         precursor.setFont(tmpFont);
         this.add(precursor);
         for(MRMDaughter d : ((MRMTransition)value).getDaughters().values()){
            JLabel dname = new JLabel("    "+d.getName());
            dname.setFont(tmpFont);
            Color dcolor = (Color)d.getGraphColor();
            dname.setForeground(dcolor);
            this.add(dname);
         }
         setBackground(isSelected ? Color.LIGHT_GRAY : Color.WHITE);
         this.setVisible(true);

         return this;
      }
    }

    public class MRMTransitionTableRenderer extends JLabel
                                            implements TableCellRenderer
    {
       Border unselectedBorder = null;
       Border selectedBorder = null;
       boolean isBordered = true;

       public MRMTransitionTableRenderer(boolean isBordered)
       {
          super();
          this.isBordered = isBordered;
          setOpaque(true); //MUST do this for background to show up.
       }

       public Component getTableCellRendererComponent(
                                JTable table, Object transition,
                                boolean isSelected, boolean hasFocus,
                                int row, int column)
       {
          setText("");
          this.setHorizontalAlignment(JLabel.RIGHT);
          Color newColor = Color.BLACK;
          setForeground(newColor);
          if (isBordered) {
             if (isSelected) {
                 if (selectedBorder == null) {
                    selectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                                  table.getSelectionBackground());
                    }
                    setBorder(selectedBorder);
                } else {
                    if (unselectedBorder == null) {
                        unselectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                                  table.getBackground());
                    }
                    setBorder(unselectedBorder);
                }
            }
            if(transition != null){
                this.setText(((MRMTransition)transition).getName());
                Font tmpFont = this.getFont();
                tmpFont = tmpFont.deriveFont(tmpFont.getSize()+4.0f);
            }
            super.setBackground(UIManager.getColor("Table.focusCellBackground"));
            return this;
        }
    }

    protected DrawingSupplier getDrawingSupplierForPlot(MRMTransition plot)
    {
        Paint colorseq[];
        if(plot != null)
        {
           int nColors  = 0;
           nColors += plot.getDaughters().size();
           colorseq = new Paint[nColors];
           int i=0;
           for(MRMDaughter d : plot.getDaughters().values()) {
              if(d == plot.getCurrentDaughter())
              {
                 colorseq[i]= Utils.paleColor((Color)plot.getCurrentDaughter().getGraphColor());
              } else {
                 colorseq[i] = Color.LIGHT_GRAY;
              }
              i++;
           }
        } else {
           colorseq = new Paint[1];
           colorseq[0] = Color.WHITE;
        }
        return new DefaultDrawingSupplier(
           colorseq,
           DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
           DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
           DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
           DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE
        );
    }

   public enum peaksData {
       Accept,
       Peptide,
       Precursor,
       Daughter {public String toString(){return "Product";}},
       CoStart {public String toString(){return "Start";}},
       CoEnd {public String toString(){return "End";}},
       CoDelta {public String toString(){return "Width";}},
       AUC,
       MaxPeak,
       Quality,
       Label,
       Code,
       LHRatio,
       Comment;
       public static NumberFormat floatFormat = NumberFormat.getNumberInstance();
       public int colno = -1;
       public Class colClass = null;
       public int colWidth = -1;
       protected boolean visible = true;
       protected int oldMinWidth = 0;
       protected int oldMaxWidth = 0;

       public void makeVisible(boolean viz) {
           if(viz == visible) return;
           if(viz) {
               peaksTable.getColumnModel().getColumn(this.colno).setMinWidth(this.oldMinWidth);
               peaksTable.getColumnModel().getColumn(this.colno).setMaxWidth(this.oldMaxWidth);
               this.visible=true;
           } else {
               this.oldMinWidth = peaksTable.getColumnModel().getColumn(this.colno).getMinWidth();
               this.oldMaxWidth = peaksTable.getColumnModel().getColumn(this.colno).getMaxWidth();
               peaksTable.getColumnModel().getColumn(this.colno).setMinWidth(0);
               peaksTable.getColumnModel().getColumn(this.colno).setMaxWidth(0);
               this.visible = false;
           }
           peaksTable.updateUI();
       }

       static {
           floatFormat.setMaximumFractionDigits(2);
           floatFormat.setMinimumFractionDigits(2);
           floatFormat.setGroupingUsed(false);

           for(peaksData pd: peaksData.values()) pd.colno = pd.ordinal();

           Accept.colClass = Boolean.class;
           Peptide.colClass = String.class;
           Precursor.colClass = MRMTransition.class;
           Daughter.colClass = MRMDaughter.class;
           CoStart.colClass = Number.class;
           CoEnd.colClass = Number.class;
           CoDelta.colClass = Number.class;
           AUC.colClass = Number.class;
           MaxPeak.colClass = Number.class;
           Quality.colClass = Number.class;
           Label.colClass = String.class;
           Code.colClass = Integer.class;
           LHRatio.colClass = Number.class;
           Comment.colClass = String.class;

           Accept.colWidth = 50;
           Peptide.colWidth = 75;
           Precursor.colWidth = 60;
           Daughter.colWidth = 50;
           CoStart.colWidth = 50;
           CoEnd.colWidth = 50;
           CoDelta.colWidth = 50;
           AUC.colWidth = 90;
           MaxPeak.colWidth = 90;
           Quality.colWidth = 50;
           Label.colWidth = 50;
           Code.colWidth = 35;
           LHRatio.colWidth = 45;
           Comment.colWidth = 150;
       }
        public String SaveFileFormat(Object thing) {
              NumberFormat nf = NumberFormat.getNumberInstance();
              nf.setGroupingUsed(false);
              if(thing == null) return "";
               switch(this) {
                   case Accept: return ((Boolean)thing).toString();
                   case Peptide: return (String) thing;
                   case Code: return ((Integer)thing).toString();
                   case Precursor:
                       nf.setMaximumFractionDigits(3);
                       nf.setMinimumFractionDigits(3);
                       return nf.format(((MRMTransition)thing).getPrecursorMz());
                   case Daughter:
                       nf.setMaximumFractionDigits(3);
                       nf.setMinimumFractionDigits(3);
                       return nf.format(((MRMDaughter)thing).getMeanMz());
                   case CoStart:
                       nf.setMaximumFractionDigits(1);
                       nf.setMinimumFractionDigits(1);
                       return nf.format(((Number)thing).doubleValue());
                   case CoEnd:
                       nf.setMaximumFractionDigits(1);
                       nf.setMinimumFractionDigits(1);
                       return nf.format(((Number)thing).doubleValue());
                   case CoDelta:
                       nf.setMaximumFractionDigits(1);
                       nf.setMinimumFractionDigits(1);
                       return nf.format(((Number)thing).doubleValue());
                   case LHRatio:
                   case AUC:
                   case MaxPeak:
                   case Quality:
                       nf.setMaximumFractionDigits(2);
                       nf.setMinimumFractionDigits(2);
                       return nf.format(((Number)thing).doubleValue());
                   case Label:
                   case Comment:
                       return (String) thing;
                   default: return "";
               }
           }
   }

    protected void scrollPeakTableToRow(int row)
    {

        peaksTable.getSelectionModel().setSelectionInterval(row,row);
        int rowHeight = peaksTable.getHeight()/((PeaksTableModel)peaksTable.getModel()).getRowCount();
        int scrollVal = (row-5)*rowHeight;
        if(scrollVal < 0) scrollVal = 0;
        peaksScrollPane.getVerticalScrollBar().setValue(scrollVal);
        peaksTable.repaint();
    }

    protected void scrollPeakTableToTransition(MRMTransition trans)
    {
        int dataIndex = -1;
        Object[][] tableData = ((PeaksTableModel)(peaksTable.getModel())).data;
        for(int j = 0; j<tableData.length; j++)
        {
            if(trans == (MRMTransition)(tableData[j][peaksData.Precursor.colno]))
            {
                dataIndex = j;
                break;
            }
        }
        if(dataIndex > -1)
        {
            scrollPeakTableToRow(dataIndex+trans.getCurrentDaughterIndex()+1);
        }
        peaksTable.repaint();
    }

    Axis findAxisOfPanel(JPanel panel, whichAxis wa)
    {
       for(Component c: panel.getComponents()) {
          if(c instanceof PanelWithChart) {
             PanelWithChart pwc = (PanelWithChart)c;
             JFreeChart jfc = pwc.getChart();
             XYPlot xyp = (XYPlot)jfc.getPlot();
             switch(wa)
             {
                case Domain: return xyp.getDomainAxis();
                case Range: return xyp.getRangeAxis();
                default: return null;
             }
          }
       }
       return null;
    }

    protected void updateDaughterCharts(boolean clear, graphZone gz)
    {
        XYSeriesCollection daughterDatasets = new XYSeriesCollection();
        if(clear)
        {
            transitionOnPlot = null;
        }
        else
        {
           daughterDatasets  = makeDaughterCollection();
        }

        //Center daughter data on graph by finding the time associated
        //with the most data (a centroid or center of mass), centering
        //this time on the plot, and scaling to assure the extreme
        //timepoints are included
        /*
        double weightedAverageTime = getWeightedAverageDaughtersTime(daughterDatasets);
        double extremeHalf = Math.max(weightedAverageTime-getMinRetentionTimeForPlot(transitionOnPlot),getMaxRetentionTimeForPlot(transitionsOnPlot)-weightedAverageTime);
        double left = weightedAverageTime-extremeHalf;
        double right = weightedAverageTime+extremeHalf;
        */
        Range r = null;
        if(transitionOnPlot != null) r = gz.getRange(transitionOnPlot);
        if(r==null) {
            r=new Range(transitionOnPlot.getMinTimeOfAllDaughters(),transitionOnPlot.getMaxTimeOfAllDaughters());
        }    
        createChartInPanel(daughterContainerPanel, daughterDatasets,
                           r.getLowerBound(),r.getUpperBound(),
                           getDrawingSupplierForPlot(transitionOnPlot),
                           whichGraph.Daughter
        );
    }

    protected void updateDaughterCharts(boolean clear)
    {
       updateDaughterCharts(clear,new defaultGraphZone());
    }

    protected void updateDaughterCharts(graphZone gz)
    {
       updateDaughterCharts(false,gz);
    }


    protected void updatePrecursorChart(boolean clear, graphZone gz)
    {
        XYSeriesCollection precursorDatasets = new XYSeriesCollection();
        if(clear)
        {
           transitionOnPlot = null;
        }
        else
        {
           if(transitionOnPlot.getGraphData() == null)
              //transitionOnPlot.setGraphData(makeParentSeries(transitionOnPlot));
              ApplicationContext.infoMessage("Precursor data for '"+transitionOnPlot.toString()+"' wasn't initialized");
              precursorDatasets.addSeries(transitionOnPlot.getGraphData());
        }
        //Center daughter data on graph by finding the time associated
        //with the most data (a centroid or center of mass), centering
        //this time on the plot, and scaling to assure the extreme
        //timepoints are included
        /*
        double weightedAverageTime = getWeightedAverageDaughtersTime(makeDaughterCollection());
        double extremeHalf = Math.max(weightedAverageTime-getMinRetentionTimeForPlot(transitionsOnPlot),getMaxRetentionTimeForPlot(transitionsOnPlot)-weightedAverageTime);
        double left = weightedAverageTime-extremeHalf;
        double right = weightedAverageTime+extremeHalf;
        */

        Range r = null;

        if(transitionOnPlot != null)r = gz.getRange(transitionOnPlot);
        if(r == null) r = new Range(0,100);
        createChartInPanel(precursorContainerPanel, precursorDatasets,
                              r.getLowerBound(),r.getUpperBound(),
                              null,
                              whichGraph.Precursor
        );
    }

    protected void updatePrecursorChart(boolean clear) {
        updatePrecursorChart(clear,new defaultGraphZone());
    }

    protected void updatePrecursorChart(graphZone gz) {
        updatePrecursorChart(false,gz);
    }


    protected void updateChartsAndFields(boolean clear)
    {
        updatePrecursorChart(clear);
        updateDaughterCharts(clear);
    }

    protected void updateChartsAndFields(graphZone gz) {
        updatePrecursorChart(gz);
        updateDaughterCharts(gz);
        //System.gc();
    }

    protected double getMinRetentionTimeForPlot(MRMTransition plot)
    {
       if(plot == null) return 0;
       double retVal = 1000000d;
       for(MRMDaughter mrmd : plot.getDaughters().values()) {
           //Note: scanVals are sorted so we only have to look at the first one
           Double testVal = mrmd.getScanVals().values().iterator().next().getTime();
           if(testVal < retVal) retVal = testVal;
       }
       return retVal;
    }


    protected double getMaxRetentionTimeForPlot(MRMTransition plot)
    {
       if(plot == null) return 5000d;
       double retVal = 0d;
       for(MRMDaughter mrmd : plot.getDaughters().values()) {
            //Note: scanVals are sorted so we only have to look at the first one
            ElutionDataPoint edpvals[] = mrmd.getScanVals().values().toArray(new ElutionDataPoint[0]);
            Double testVal = edpvals[edpvals.length-1].getTime();
            if(testVal > retVal) retVal = testVal;
        }
        return retVal;
    }    

    protected XYSeries makeParentSeries(MRMTransition parent)
    {

       float minMz = parent.getPrecursorMz() - _precursorDiscoveryMzTolerance - _precursorChromatogramWindow;
       float maxMz = parent.getPrecursorMz() + _precursorDiscoveryMzTolerance + _precursorChromatogramWindow;
       NumberFormat nf = NumberFormat.getNumberInstance();
       nf.setMaximumFractionDigits(1);
       nf.setMinimumFractionDigits(1);

       XYSeries result = new XYSeries(parent.getName()+"\u00B1"+nf.format( _precursorDiscoveryMzTolerance + _precursorChromatogramWindow));
       //precursor scans
       for (int i=parent.getMinScanOfDaughters(); i<=parent.getMaxScanOfDaughters(); i++)
       {
          int scanNum = _run.getIndexForScanNum(i);
          if (scanNum <= 0) continue;
          MSRun.MSScan ms1Scan = _run.getScan(scanNum);
          boolean sim_only = _sim;
          String scanType = ms1Scan.getScanType();
          if(!sim_only || (sim_only && scanType.equalsIgnoreCase("SIM")))
             result.add(ms1Scan.getDoubleRetentionTime(),Utils.getMaxIntensityForScan(ms1Scan, minMz, maxMz));
       }
       return result;
    }

    public boolean precursorChromatogramEmpty(MRMTransition parent){
       if(parent == null) return true;
       XYSeries pData;
       if((pData =parent.getGraphData()) == null) return true;
       for(Object xydio: pData.getItems()){
           XYDataItem xydi = (XYDataItem)xydio;
           if(xydi.getY().doubleValue() != 0.0) return false;
       }
       return true;
    }

    private XYSeriesCollection makeDaughterCollection() {
        XYSeriesCollection retVal = new XYSeriesCollection();

        for(MRMDaughter d: transitionOnPlot.getDaughters().values()){
           if(d.getGraphData()== null)
              d.setGraphData(d.makeDaughterSeries());
//if ((Boolean)(((PeaksTableModel)(peaksTable.getModel())).data[d.getElutionDataTableRow()][peaksData.Accept.colno]))
              retVal.addSeries(d.getGraphData());
           }
         return retVal;
    }

    public class defaultGraphZone implements graphZone {
        public Range getRange(MRMTransition mrmt) {
            double left = 0;
            double right = 0;
            double elutionWidth = 0;
            if(mrmt.getElutionRegionStart() != -1 && mrmt.getElutionRegionEnd() != -1) {
                elutionWidth = mrmt.getElutionRegionEnd()-mrmt.getElutionRegionStart();
                left = Math.max(0.0,mrmt.getElutionRegionStart()-(elutionWidth*0.01));
                right =  mrmt.getElutionRegionEnd() + (elutionWidth*0.01);
            } else {
                left = mrmt.calculateMinOfAllBestDaughterCurves();
                right = mrmt.calculateMaxOfAllBestDaughterCurves();
            }
            if(right <= 0) {
                if(mrmt.getParentData() != null  &&  mrmt.getParentData().getItemCount()>0) {
                    left = mrmt.getParentData().getDataItem(0).getX().doubleValue();
                    right = mrmt.getParentData().getDataItem(mrmt.getParentData().getItemCount()-1).getX().doubleValue();
                    left = Math.max(0,left-((right-left)*0.05));
                    right = right + ((right-left)*0.05);
                }
            }
            if(left>=right)return null;
            return new Range(left,right);
        }
    }

    public class zoomedGraphZone implements graphZone {
        public Range getRange(MRMTransition mrmt) {
           double left = getMinRetentionTimeForPlot(transitionOnPlot);
           double right = getMaxRetentionTimeForPlot(transitionOnPlot);
           if(left >= right)return null;  
           return new Range(left, right);
        }
    }
    // zoom in on best elution curve
    public class zoomToCurveGraphZone implements graphZone {
        public Range getRange(MRMTransition mrmt) {
           MRMDaughter mrmd = mrmt.getCurrentDaughter();
           zoomedGraphZone zgd = new zoomedGraphZone();
           if(mrmd == null) return zgd.getRange(mrmt);
           ElutionCurve ec = mrmd.getBestElutionCurve();
           if(ec == null) return zgd.getRange(mrmt);

           double left = ec.getMinElutionTimeSecs();
           double right = ec.getMaxElutionTimeSecs();
           if(left >= right)return null;
           return new Range(left, right);
        }
    }

    protected void createChartInPanelPrecursorTasksOnly(XYPlot xyp) {
           boolean shouldDisplay = false;

           if(!precursorChromatogramEmpty(transitionOnPlot)) {
              shouldDisplay = true;

              if(!transitionOnPlot.getElutionCurves().isEmpty()) {
                 ElutionCurveStrategy ecs = transitionOnPlot.getElutionCurves().values().iterator().next();
                 List<ElutionCurve> ecl = ecs.getParentCurves();
                 if(ecl != null) {
                    for(ElutionCurve e: ecl) {
                       List<Line2D.Double> ll2dd = e.getSegments();
                       for(Line2D.Double l2dd: ll2dd) {
                          xyp.addAnnotation(Utils.line2Annotation(l2dd,new BasicStroke(2.0f),(ecs.isBestParentCurve(e)) ? Color.GREEN :Color.LIGHT_GRAY));
                       }
                    }
                 }
              }
           }

           precursorContainerContainerPanel.setVisible(shouldDisplay);
//           updateDaughterCharts(false);
           if(this.getWidth() >= 100 && this.getHeight() >= 100)
               this.setPreferredSize(new Dimension(this.getWidth(),this.getHeight()));
           pack();
    }

    protected void createChartInPanelDaughterTasksOnly(XYPlot xyp) {
       XYSeries coloredDataset = transitionOnPlot.getCurrentDaughter().getGraphData();
       Paint daughterColor = Utils.paleColor((Color) transitionOnPlot.getCurrentDaughter().getGraphColor());
       ArrayList<XYLineAnnotation> coloredDaughters = new ArrayList<XYLineAnnotation>();
       //Trace calculated elution curves over data spikes
       if(transitionOnPlot.getElutionCurves() != null && !transitionOnPlot.getElutionCurves().isEmpty()) {
          MRMDaughter curDaughter = transitionOnPlot.getCurrentDaughter();
          ElutionCurveStrategy ecs = transitionOnPlot.getElutionCurves().get(curDaughter);
          //Is current daughter rejected?
          Boolean accepted =
             (Boolean)((PeaksTableModel)peaksTable.getModel()).data[curDaughter.getElutionDataTableRow()][peaksData.Accept.colno];
          if(accepted == null || !accepted) {
             xyp.setBackgroundPaint(
                new Color(255,230,230)
             );
          }
          List<ElutionCurve> ecl = ecs.getDaughterCurves();
          if(ecl != null) {
             for(ElutionCurve e: ecl) {
                List<Line2D.Double> ll2dd = e.getSegments();
                for(Line2D.Double l2dd: ll2dd) {
                   xyp.addAnnotation(Utils.line2Annotation(l2dd,new BasicStroke(2.0f),ecs.isBestDaughterCurve(e) ? Color.BLACK : Color.LIGHT_GRAY));
                }
             }
          }
       }

       // If there is a valid "current" daughter draw the spikes in the daughter's color
       // as annotations (sensu JFree)

       if(coloredDataset != null) {
          int nOfPoints = coloredDataset.getItemCount();
          for(int i = 0; i<(nOfPoints-1); i++){
             XYDataItem p1 = coloredDataset.getDataItem(i);
             XYDataItem p2 = coloredDataset.getDataItem(i+1);
             coloredDaughters.add(
                new XYLineAnnotation(p1.getX().doubleValue(),p1.getY().doubleValue(),p2.getX().doubleValue(),p2.getY().doubleValue(),new BasicStroke(1.5f),daughterColor)
             );
          }
       }
       if(_traceAllFragments) {
          for(MRMDaughter d: transitionOnPlot.getDaughters().values()) {
             if(d == transitionOnPlot.getCurrentDaughter()) continue;
             XYSeries curXYSeries = d.getContinDaughterData();
             if(curXYSeries == null || curXYSeries.getItemCount() == 0) continue;
             if(d.getBestElutionCurve() == null) continue;
             int nOfPoints = curXYSeries.getItemCount();
             for(int i = 0; i<(nOfPoints-1); i++){
                XYDataItem p1 = curXYSeries.getDataItem(i);
                XYDataItem p2 = curXYSeries.getDataItem(i+1);
                coloredDaughters.add(
                   new XYLineAnnotation(p1.getX().doubleValue(),p1.getY().doubleValue(),p2.getX().doubleValue(),p2.getY().doubleValue(),new BasicStroke(1f),Utils.paleColor((Color)d.getGraphColor()))
                );
             }
          }
       }
       if(coloredDaughters != null){
          for(XYLineAnnotation xyla: coloredDaughters) {
             xyp.addAnnotation(xyla);
           }
       }
       coloredDaughters.clear();
        
       //Display L or H label in upper left hand corner
       Range xRange = xyp.getDomainAxis().getRange();
       Range yRange = xyp.getRangeAxis().getRange();
       XYTextAnnotation lab =
           new XYTextAnnotation(
             (String)((PeaksTableModel)peaksTable.getModel()).data[transitionOnPlot.getCurrentDaughter().getElutionDataTableRow()][peaksData.Label.colno],
             xRange.getUpperBound()-(0.05*xRange.getLength()),
             yRange.getUpperBound()-(0.05*yRange.getLength())
           );
       lab.setFont(lab.getFont().deriveFont(Font.BOLD,40.0F));
       xyp.addAnnotation(lab);
   }

    public XYPlot oldPrecursorChart = null;
    public XYPlot oldProductChart = null;

    protected void clearPreviousChartJunk(XYPlot xyp) {
        if(xyp != null) {
            xyp.clearAnnotations();
            xyp.clearDomainAxes();
            xyp.clearRangeAxes();
            xyp.setDataset(null);
        }
    }

   /**
     * Draw a chart in a panel.  Good for precursors and daughters
     * @param parentPanel
     * @param dataset
     * @param domainMin
     * @param domainMax
     * @param supplier
     * @param wg
     */

   
    protected void createChartInPanel(JPanel parentPanel,
                                       XYSeriesCollection dataset,
                                       Double domainMin, Double domainMax,
                                       DrawingSupplier supplier,
                                       whichGraph wg
                                       )
    {
        if(precursorChromatogramEmpty(transitionOnPlot) && wg == whichGraph.Precursor){
            precursorContainerContainerPanel.setVisible(false);
            if(this.getWidth() >= 100 && this.getHeight() >= 100)
                this.setPreferredSize(new Dimension(this.getWidth(),this.getHeight()));
            pack();
            return;
        }
        switch(wg) {
           case Precursor:
               clearPreviousChartJunk(oldPrecursorChart);
               oldPrecursorChart = null;
               break;
           case Daughter:
               clearPreviousChartJunk(oldProductChart);
               oldProductChart = null;
               break;
        }

        JFreeChart chart =
                ChartFactory.createXYLineChart(null,"seconds", null, dataset,
                        PlotOrientation.VERTICAL, true, false, false);

        XYPlot xyp = (XYPlot)(chart.getPlot());
        if(supplier != null)
        {
           xyp.setDrawingSupplier(supplier);
        } else {
           xyp.setDrawingSupplier(Utils.plainDrawingSupplier(Color.LIGHT_GRAY));
        }
        xyp.setSeriesRenderingOrder(SeriesRenderingOrder.REVERSE);

        CenterZoomNumberAxis axisDomain = new CenterZoomNumberAxis("seconds");
        axisDomain.setAutoRangeIncludesZero(false);
        axisDomain.setRange(Math.max(0.0,domainMin), domainMax);
        axisDomain.addChangeListener(new domainAxisZoomCoordinator(axisDomain));

        xyp.clearAnnotations();
        xyp.setDomainAxis(axisDomain);
        xyp.getDomainAxis().setAutoRange(false);
        xyp.getRangeAxis().setAutoRange(false);
        XYLineAndShapeRenderer xylsr = (XYLineAndShapeRenderer) xyp.getRenderer();
        xylsr.setLegendLine(Utils.legendThing(16,6));

        xylsr.setShapesFilled(true);
        xylsr.setBaseShapesFilled(true);
        PanelWithChart panelWithChart = new PanelWithChart(chart);
        ChartPanel cp = panelWithChart.getChartPanel();
        cp.removeMouseListener(cp);
        cp.removeMouseMotionListener(cp);
        if(peaksTable != null) {
           MRMerMouseListener mml = new MRMerMouseListener(cp,(PeaksTableModel)peaksTable.getModel());
           cp.addMouseListener(mml);
           cp.addMouseMotionListener(mml);
        }
        cp.setPreferredSize(
                new Dimension(parentPanel.getWidth(), parentPanel.getHeight() - 10));
        cp.setDomainZoomable(true);
        cp.setRangeZoomable(false);
        cp.setPopupMenu(null);
        parentPanel.removeAll();
        parentPanel.add(panelWithChart);

        switch(wg) {
           case Precursor:
              createChartInPanelPrecursorTasksOnly(xyp);
              oldPrecursorChart = xyp;
              break;
           case Daughter:
              createChartInPanelDaughterTasksOnly(xyp);
              oldProductChart = xyp;
              break;
        }        
        parentPanel.updateUI();
        listTransition.requestFocus();
    }

    private int indexOfMaxY(XYSeries s)
    {
        double _maxy = -10000000000.0;
        int retVal = -1;
        for(int i = 0; i<s.getItemCount(); i++) {
           XYDataItem xydi = s.getDataItem(i);
           if(_maxy < xydi.getY().doubleValue()) {
              _maxy = xydi.getY().doubleValue();
              retVal = i;
           }
        }
        return retVal;
    }

    private static NumberFormat _transitionNumberFormat =
       NumberFormat.getNumberInstance();

    static {
        _transitionNumberFormat.setMaximumFractionDigits(4);
        _transitionNumberFormat.setMinimumFractionDigits(4);
    }

    public interface mzXML2TransitionArray {
        public MRMTransition[] reParse(MSRun run);
    }

    protected class originalReParser implements mzXML2TransitionArray {

           
        public MRMTransition[] reParse(MSRun run) {
           Vector<MRMTransition> transitions = new Vector<MRMTransition>();
           Map<Float, MRMDaughter> meanDaughterTransitionMap = new HashMap<Float, MRMDaughter>();

           for (MSRun.MSScan ms2Scan : run.getMS2Scans())
           {
               if ("SRM".equalsIgnoreCase(ms2Scan.getScanType()) || "MRM".equalsIgnoreCase(ms2Scan.getScanType()))
               {
                // See if a parent MRMTransition within the mz tolerance exists
                //    If not.  Create one.
                float testPrecursor = ms2Scan.getPrecursorMz();
                MRMTransition curTrans = null;
                for(MRMTransition testTrans: transitions) {
                    if(testTrans.getPrecursorMz() >= (testPrecursor-_precursorDiscoveryMzTolerance) &&
                       testTrans.getPrecursorMz() <= (testPrecursor+_precursorDiscoveryMzTolerance)) {
                          curTrans = testTrans;
                          break;
                    }
                }
                if(curTrans == null) {
                    curTrans = new MRMTransition(testPrecursor,run);
                    curTrans.setName(_transitionNumberFormat.format(curTrans.getPrecursorMz()));
                    transitions.add(curTrans);
                }
                // See if a daughter mean mz exists yet for that parent
                //    If not. Create one.  Attach it to parent.
                float meanDaughter = (ms2Scan.getLowMz() + ms2Scan.getHighMz()) / 2;
                MRMDaughter curDaughter = null;
                for(MRMDaughter testDaughter: curTrans.getDaughters().values()) {
                    if(meanDaughter >= (testDaughter.getMeanMz() - _daughterMzTolerance) &&
                       meanDaughter <= (testDaughter.getMeanMz() + _daughterMzTolerance)
                    ) {
                        curDaughter = testDaughter;
                        break;
                    }
                }
                if(curDaughter == null)
                {
                    curDaughter = new MRMDaughter(
                            meanDaughter,
                            ms2Scan.getLowMz(),
                            ms2Scan.getHighMz(),
                            ms2Scan.getNum(),
                            ms2Scan.getNum(),
                            curTrans
                     );
                     curDaughter.setGraphColor(MRMTransition.COLOR_SERIES[curTrans.getDaughters().size() % (MRMTransition.COLOR_SERIES.length)]);
                     curDaughter.setName(_transitionNumberFormat.format(curDaughter.getPrecursor().getPrecursorMz())+"/"+_transitionNumberFormat.format(curDaughter.getMeanMz()));
                     meanDaughterTransitionMap.put(curDaughter.getMeanMz(),curDaughter);
                     curTrans.getDaughters().put(curDaughter.getMeanMz(),curDaughter);
                     //make sure this line doesn't duplicate initial scan number
                     curDaughter.addScanVal(ms2Scan.getNum(),new ElutionDataPoint(ms2Scan.getDoubleRetentionTime(),Utils.getMaxIntensityForScan(ms2Scan, curDaughter.getLowMz(), curDaughter.getHighMz())));
//allDaughters.add(new ReactionClusterable(curTrans.getPrecursorMz(),curDaughter.getMeanMz()));
                } else {
                // Add current scan to daughter.
                  curDaughter.addScanVal(ms2Scan.getNum(),new ElutionDataPoint(ms2Scan.getDoubleRetentionTime(),Utils.getMaxIntensityForScan(ms2Scan, curDaughter.getLowMz(), curDaughter.getHighMz())));
//allDaughters.add(new ReactionClusterable(curTrans.getPrecursorMz(),curDaughter.getMeanMz()));

                }
            }
        }
        ArrayList<MRMDaughter> mrmTransitionLists = new ArrayList<MRMDaughter>();

        for (MRMDaughter daughter : meanDaughterTransitionMap.values())
        {
            mrmTransitionLists.add(daughter);
        }
        MRMTransition[] result = transitions.toArray(new MRMTransition[0]);
        Arrays.sort(result, new MRMTransition.PrecursorMzComparator());
        return result;
      }
    }

    protected class thermoReParser implements mzXML2TransitionArray {
        protected Range[] parseFilterLine(String fline)
        {
           String rangeStrings[] = fline.split("\\[")[1].split(",");
           rangeStrings[rangeStrings.length-1] = rangeStrings[rangeStrings.length-1].replace("]","");
           Range retVal[] = new Range[rangeStrings.length];
           for(int i = 0; i<rangeStrings.length; i++) {
               String ends[] = rangeStrings[i].split("-");
               retVal[i] = new Range(Double.parseDouble(ends[0]),Double.parseDouble(ends[1]));
           }
           return retVal;
        }

        public MRMTransition[] reParse(MSRun run) {
            // New (Post Jan '08) thermo files have "filterLine" attributes
            // on the <run> element.  These define transitions.  There will
            // be one "peak" in the <peaks> element for each daughter.
            // Considerations:  each daughter must (should) be assigned a
            // distinct time.  So the entire ms2/MRM scan should be divided
            // into n retention time segments for n daughters.  This means
            // we have to know the value of the NEXT scan retention time
            // in order to get some estimate for the length of the ms2
            // scan.  It also means we'll keep around average times so that
            // we can process the last scan.

            //for final mean scan time
            double sumRetentionTimeDeltas = 0.0d;
            int countMRMScans = 0;

            //places to keep precursors and products
            Vector<MRMTransition> transitions = new Vector<MRMTransition>();
            Map<Float, MRMDaughter> meanProductTransitionMap = new HashMap<Float, MRMDaughter>();

            for (MSRun.MSScan ms2Scan : run.getMS2Scans())
            {
                if ("SRM".equalsIgnoreCase(ms2Scan.getScanType()) || "MRM".equalsIgnoreCase(ms2Scan.getScanType()))
                {
                    try{
                       if(ms2Scan.getFilterLine() == null) {
                          throw new Exception("No filter line in Thermo MRM scan '"+ms2Scan.toString()+"'");
                       }
                        float testPrecursor = ms2Scan.getPrecursorMz();
                        MRMTransition curTrans = null;
                        // Have we seen a transition with this precursor value yet?
                        for(MRMTransition testTrans: transitions) {
                            if(testTrans.getPrecursorMz() >= (testPrecursor-_precursorDiscoveryMzTolerance) &&
                               testTrans.getPrecursorMz() <= (testPrecursor+_precursorDiscoveryMzTolerance)) {
                                  curTrans = testTrans;
                                  break;
                            }
                        }
                        if(curTrans == null) {
                            curTrans = new MRMTransition(testPrecursor,run);
                            curTrans.setName(_transitionNumberFormat.format(curTrans.getPrecursorMz()));
                            transitions.add(curTrans);
                        }
                        //Determine length of scan
                        double scanLen = 0d;
                        MSRun.MSScan nextScan = _run.getScanByNum(ms2Scan.getNum()+1);
                        if(nextScan != null) {
                           scanLen = nextScan.getDoubleRetentionTime()-ms2Scan.getDoubleRetentionTime();
                           sumRetentionTimeDeltas += scanLen;
                           countMRMScans++;
                        } else {
                           scanLen = sumRetentionTimeDeltas/countMRMScans; //mean daughter scan time
                        }
                        //Create a daughter (if necessary) or just a daughter datapoint
                        //Assume there are as many points in the "spectrum" as there are
                        //daughters.  Assume they come in the same order as they are described
                        //in the filterLine.  We can change the assumptions and complicate the
                        //code if necessary.
                        Range productRanges[] = parseFilterLine(ms2Scan.getFilterLine());
                        Vector<MRMDaughter> productRangeDaughters = new Vector<MRMDaughter>();
                        for(int i = 0; i<productRanges.length; i++){
                           Range daughterCenter = productRanges[i];
                           float meanDaughter = (float)daughterCenter.getCentralValue();
                           MRMDaughter curDaughter = null;
                           for(MRMDaughter testDaughter: curTrans.getDaughters().values()) {
                              if(
                                 meanDaughter >= (testDaughter.getMeanMz() - _daughterMzTolerance) &&
                                 meanDaughter <= (testDaughter.getMeanMz() + _daughterMzTolerance)
                              ) {
                                 curDaughter = testDaughter;
                                 break;
                              }
                           }
                           if(curDaughter == null){
                              curDaughter =
                                 new MRMDaughter(
                                    meanDaughter,
                                    (float)daughterCenter.getLowerBound(),
                                    (float)daughterCenter.getUpperBound(),
                                    ms2Scan.getNum(),
                                    ms2Scan.getNum(),
                                    curTrans
                                 );
                              curDaughter.setGraphColor(MRMTransition.COLOR_SERIES[curTrans.getDaughters().size() % (MRMTransition.COLOR_SERIES.length)]);
                              curDaughter.setName(_transitionNumberFormat.format(curDaughter.getPrecursor().getPrecursorMz())+"/"+_transitionNumberFormat.format(curDaughter.getMeanMz()));
                              meanProductTransitionMap.put(curDaughter.getMeanMz(),curDaughter);
                              curTrans.getDaughters().put(curDaughter.getMeanMz(),curDaughter);
                           }
                           productRangeDaughters.add(curDaughter);
                        }
                        //Now we know a daughter defined by the filterLine is associated with
                        //with the precursor.  We have to find which spectrum datapoint goes
                        //with which daughter.
                        int dcount = -1;
                        for(MRMDaughter mrmd: productRangeDaughters) {
                            dcount++;
                            int spectrumPointCount = 0;
                            double curMZ = -1d;
                            double curIntensity = -1d;
                            for(int j = 0; j<ms2Scan.getSpectrum()[0].length; j++) {
                               if((ms2Scan.getSpectrum()[0][j] >= (mrmd.getMeanMz()-_daughterMzTolerance)) && (ms2Scan.getSpectrum()[0][j] <= (mrmd.getMeanMz()+_daughterMzTolerance))) {
                                   curMZ = ms2Scan.getSpectrum()[0][j];
                                   curIntensity = ms2Scan.getSpectrum()[1][j];
                                   spectrumPointCount++;
                               }
                            }
                            if(spectrumPointCount == 1) {
                               mrmd.addScanVal(
                                 ms2Scan.getNum(),
                                 new ElutionDataPoint(
                                     ms2Scan.getDoubleRetentionTime()+((scanLen/productRanges.length)*dcount),
                                     curIntensity
                                 )
                               );
                            } else {
                                if(spectrumPointCount > 1) {
                                    ApplicationContext.infoMessage("More than one spectrum point can belong to same daughter in same scan");
                                } else {
                                    if(spectrumPointCount == 0) {
                                       mrmd.addScanVal(
                                          ms2Scan.getNum(),
                                          new ElutionDataPoint(
                                             ms2Scan.getDoubleRetentionTime()+((scanLen/productRanges.length)*dcount),
                                             0d
                                          )
                                     );
                                     ApplicationContext.infoMessage("No datapoint found matching daughter "+mrmd+" for scan "+ms2Scan.getNum());
                                    }
                                }
                            }
                        }
                    } catch (Exception e) {
                       ApplicationContext.infoMessage("Can't parse scan "+ms2Scan.getNum()+" in thermo reparser: "+e);
                    }
                }
            }
            ArrayList<MRMDaughter> mrmTransitionLists = new ArrayList<MRMDaughter>();

            for (MRMDaughter daughter : meanProductTransitionMap.values())
            {
                mrmTransitionLists.add(daughter);
            }
            MRMTransition[] result = transitions.toArray(new MRMTransition[0]);
            Arrays.sort(result, new MRMTransition.PrecursorMzComparator());
            return result;
         }
    }

    /**
       /**
     * Detect all the transitions and load into an array
     * @param run
     * @return
     *
     * There is a lot of legacy code here to be cleaned up later
     *
     */
    private MRMTransition[] loadMRMTransitions(MSRun run)
    {
        for (MSRun.MSScan ms2Scan : run.getMS2Scans())
        {
            if ("SRM".equalsIgnoreCase(ms2Scan.getScanType()) || "MRM".equalsIgnoreCase(ms2Scan.getScanType()))
            {
                if(ms2Scan.getFilterLine() != null) {
                    return(new thermoReParser()).reParse(run);
                } else return (new originalReParser()).reParse(run);
            }
            break;
        }
        return(new originalReParser()).reParse(run);
    }

    public class ReactionClusterable implements Clusterer2D.Clusterable {
        public float getPrecursorMZ() {
            return precursorMZ;
        }

        public void setPrecursorMZ(float pMZ) {
            this.precursorMZ = pMZ;
        }

        private float precursorMZ;

        public float getDaughterMZ() {
            return daughterMZ;
        }

        public void setDaughterMZ(float daughterMZ) {
            this.daughterMZ = daughterMZ;
        }

        private float daughterMZ;

        public double getDimension1Value() {return (double) getPrecursorMZ();}
        public double getDimension2Value() {return (double) getDaughterMZ();}

        private ReactionClusterable(){};

        public ReactionClusterable(float p, float d) {
            this.precursorMZ = p;
            this.daughterMZ = d;
        }

        public ReactionClusterable(MRMDaughter mrmd) {
            if(mrmd != null && mrmd.getPrecursor() != null) this.precursorMZ = mrmd.getPrecursor().getPrecursorMz();
            if(mrmd != null) this.daughterMZ = mrmd.getMeanMz();
        }
    }


}