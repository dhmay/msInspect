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

package org.fhcrc.cpl.toolbox.gui.chart;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.gui.ImagePanel;
import org.fhcrc.cpl.toolbox.ApplicationContext;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * PanelWithChart implementation to make it easy to put out Line Charts
 * If you want to do anything super-serious with the chart, use
 * getChart() and getRenderer()
 */
public class PanelWithBlindImageChart extends PanelWithChart
{
    static Logger _log = Logger.getLogger(PanelWithBlindImageChart.class);

    protected List<BufferedImage> allImages = null;
    protected BufferedImage image = null;
    protected ImagePanel imagePanel = null;
    protected JPopupMenu popupMenu = null;

    protected JScrollPane scrollPane = null;

    int currentDisplayedImageIndex = 0;


    public PanelWithBlindImageChart()
    {
        super();
        setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.insets = new Insets(0, 0, 0, 0);

        scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        add(scrollPane, gbc);
        addComponentListener(new ResizeListener());
        allImages = new ArrayList<BufferedImage>();
    }

    public PanelWithBlindImageChart(String dataSetName)
    {
        this();
        setName(dataSetName);
    }

    public void setPreferredSize(Dimension dimension)
    {
        scrollPane.setPreferredSize(dimension);
        updateUI();
//System.err.println("Setting preferred scrollpane size of " + getName() + " to " + dimension.getWidth() + ", " + dimension.getHeight());
    }


    public PanelWithBlindImageChart(File imageFile, String name)
    {
        this(name);
        try
        {
            setImage(imageFile);
        }
        catch (IOException e)
        {
            ApplicationContext.errorMessage("Error while displaying chart image " +
                                            imageFile.getAbsolutePath(),e);
        }
    }

    public PanelWithBlindImageChart(BufferedImage image, String name)
    {
        this(name);
        setImage(image);
    }

    /**
     * hack.  Won't actually resize
     * @param width
     * @param height
     * @return
     */
    public BufferedImage createImage(int width, int height)
    {
        return image;
    }

    public void setImages(List<BufferedImage> images)
    {
        allImages = images;
        setImage(allImages.get(0));
        currentDisplayedImageIndex = 0;
    }

    public void cycleImageNext()
    {
        currentDisplayedImageIndex++;
        if (currentDisplayedImageIndex > allImages.size()-1)
            currentDisplayedImageIndex = 0;
        displayImage(currentDisplayedImageIndex);
    }

    public void cycleImagePrevious()
    {
        currentDisplayedImageIndex--;
        if (currentDisplayedImageIndex < 0)
            currentDisplayedImageIndex = allImages.size()-1;
        displayImage(currentDisplayedImageIndex);
    }

    public void displayImage(int imageIndex)
    {
        setImage(allImages.get(imageIndex));
    }

    /**
     * Resize all chart panels when this panel is resized
     */
    protected class ResizeListener implements ComponentListener
    {
        public void componentResized(ComponentEvent event)
        {
//System.err.println("PWBIC componentResized");            
//System.err.println("Setting scrollpane to " + getWidth() + ", " + getHeight());            
//              scrollPane.setPreferredSize(new Dimension(getWidth(), getHeight()));
        }

        public void componentMoved(ComponentEvent event)  {}
        public void componentShown(ComponentEvent event)  {}
        public void componentHidden(ComponentEvent event)  {}
    }

    public void setImage(BufferedImage image)
    {
//        setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));

        this.image = image;
        imagePanel = new ImagePanel(image);
        imagePanel.setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));
        scrollPane.setViewportView(imagePanel);
//scrollPane.setPreferredSize(new Dimension(300, 300));

        popupMenu = new JPopupMenu();
        JMenuItem saveImageMenuItem = new JMenuItem("Save Image");
        saveImageMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JFileChooser fc = new JFileChooser();
                    int chooserStatus = fc.showOpenDialog(PanelWithBlindImageChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = fc.getSelectedFile();
                    try
                    {
                        saveChartToImageFile(outFile);
                    }
                    catch (IOException e)
                    {
                        ApplicationContext.errorMessage("Failed to save image to file " + outFile.getAbsolutePath(), e);
                    }
                }
            }
        );
        popupMenu.add(saveImageMenuItem);

        PopupListener myPopupListener = new PopupListener(popupMenu);
        scrollPane.addMouseListener(myPopupListener);
        imagePanel.addMouseListener(myPopupListener);
        this.addMouseListener(myPopupListener);
    }

    protected class PopupListener extends MouseAdapter
    {
        protected JPopupMenu popup = null;

        public PopupListener(JPopupMenu popup)
        {
            this.popup = popup;
        }

        public void mousePressed(MouseEvent e) {
            maybeShowPopup(e);
        }

        public void mouseReleased(MouseEvent e) {
            maybeShowPopup(e);
        }

        private void maybeShowPopup(MouseEvent e) {
            if (e.isPopupTrigger()) {
                popup.show(e.getComponent(),
                        e.getX(), e.getY());
            }
        }
    }

    public void setImage(File imageFile)
            throws IOException
    {
        image = ImageIO.read(imageFile);
        _log.debug("Loaded image " + imageFile.getAbsolutePath() + ", size: " +
                image.getWidth() + ", " + image.getHeight());
        allImages = new ArrayList<BufferedImage>();
        allImages.add(image);
        setImage(image);        
    }

    public void setImageFiles(List<File> imageFiles)
            throws IOException
    {
        List<BufferedImage> newImages = new ArrayList<BufferedImage>();
        for (File imageFile : imageFiles)
        {
            BufferedImage currentImage = ImageIO.read(imageFile);
            newImages.add(currentImage);
            _log.debug("Loaded image " + imageFile.getAbsolutePath() + ", size: " + currentImage.getWidth() + ", " +
                    currentImage.getHeight());
//newImages.add(new BufferedImage(100, 100, BufferedImage.TYPE_INT_RGB));
        }
        setImages(newImages);
    }




    public void addItemToPopupMenu(JMenuItem item)
    {
    }

    public void addSeparatorToPopupMenu()
    {
    }



    /**
     * Don't add the usual "save data" items to the popup
     */
    protected void initPopupMenu()
    {
        JMenuItem saveImageMenuItem = new JMenuItem("Save Chart");
        saveImageMenuItem.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JFileChooser fc = new JFileChooser();
                    int chooserStatus = fc.showOpenDialog(PanelWithBlindImageChart.this);
                    //if user didn't hit OK, ignore
                    if (chooserStatus != JFileChooser.APPROVE_OPTION)
                        return;

                    File outFile = fc.getSelectedFile();
                    try
                    {
                        saveChartToImageFile(outFile);
                    }
                    catch (Exception e)
                    {
                        ApplicationContext.errorMessage("Error saving image file " +
                                outFile.getAbsolutePath(),e);
                    }
                }
            }
        );
        addItemToPopupMenu(saveImageMenuItem);
    }

    public void saveChartToImageFile(File outFile) throws IOException
    {
        ImageIO.write(image,"png",outFile);
    }

    public void saveAllImagesToFiles(File outDir) throws IOException
    {
        for (int i=0; i<allImages.size(); i++)
        {
            BufferedImage currentImage = allImages.get(i);
            String imageNumString = "" + (i+1);
            if ((i+1)<10)
                imageNumString = "0" + imageNumString;
            ImageIO.write(currentImage,"png",new File(outDir, "image" + imageNumString + ".png"));
        }

    }

    /**
     * Do nothing
     * @param outFile
     * @param delimiter
     */
    protected void saveChartDataToFile(File outFile, String delimiter)
    {

    }

    public BufferedImage getImage()
    {
        return image;
    }
}
