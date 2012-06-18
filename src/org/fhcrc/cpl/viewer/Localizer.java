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
package org.fhcrc.cpl.viewer;

import org.apache.log4j.Logger;

import javax.swing.plaf.FontUIResource;
import javax.swing.*;
import java.util.Locale;
import java.util.Vector;
import java.util.HashMap;
import java.util.prefs.Preferences;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.fhcrc.cpl.toolbox.TextProvider;
import org.fhcrc.cpl.toolbox.ApplicationContext;
import org.fhcrc.cpl.viewer.gui.WorkbenchFrame;
import org.swixml.SwingEngine;

/**
 * Provides methods for defining aspects of behavior for different Locales, and for
 * switching between Locales
 */
public class Localizer
{
    private static Logger _log = Logger.getLogger(Localizer.class);

    protected static HashMap _languageProperties = null;

    public static final String DEFAULT_LOCALE_STRING = "en";
    public static final Locale DEFAULT_LOCALE = new Locale(DEFAULT_LOCALE_STRING);

    //Index of the language and font menus within the menubar.
    //It would be better if we could reference this by name somehow
    protected static int LANGUAGE_MENU_INDEX = 3;
    protected static int FONT_MENU_INDEX = 4;

    protected static Locale _locale = DEFAULT_LOCALE;
    protected static Font _font;

    //for storing language and font preference
    protected static final Preferences _prefs = Preferences.userNodeForPackage(Application.class);
    protected static final String _localePrefKey = "lastLocale";
    protected static final String _fontPrefKeyPrefix = "lastFont_";

    protected static Locale[] supportedLocales = {new Locale("en")
                                                , new Locale("zh")
                                                 };

    //constants for defining properties for each language
    protected static final String FONT_TEST_STRING="FONT_TEST_STRING";
    protected static final String LABEL_SCALE="LABEL_SCALE";
    protected static final String CONTAINER_WIDTH_SCALE="CONTAINER_WIDTH_SCALE";
    protected static final String CONTAINER_HEIGHT_SCALE="CONTAINER_HEIGHT_SCALE";
    protected static final String FONT_SIZE="FONT_SIZE";
    protected static final String FONT_NAME="FONT_NAME";
    protected static final String FONT_STYLE="FONT_STYLE";
    protected static final String REQUIRES_FONT_CHOICE="REQUIRES_FONT_CHOICE";




    public Localizer()
    {

    }


    /**
     * Initialize properties for each language
     */
    protected static void initLanguageProperties()
    {
        _languageProperties = new HashMap();

        HashMap enProperties = new HashMap();
        enProperties.put(FONT_TEST_STRING,"a");
        enProperties.put(LABEL_SCALE,"1");
        enProperties.put(CONTAINER_WIDTH_SCALE,"1");
        enProperties.put(CONTAINER_HEIGHT_SCALE,"1");
        enProperties.put(FONT_SIZE,"12");
        enProperties.put(FONT_NAME,"Dialog");
        enProperties.put(FONT_STYLE,Integer.toString(Font.PLAIN));
        enProperties.put(REQUIRES_FONT_CHOICE,"FALSE");
        _languageProperties.put("en",enProperties);

        HashMap zhProperties = new HashMap();
        zhProperties.put(FONT_TEST_STRING,"\u4e00");
        zhProperties.put(LABEL_SCALE,"1");
        zhProperties.put(CONTAINER_WIDTH_SCALE,"1.1");
        zhProperties.put(CONTAINER_HEIGHT_SCALE,"1.1");
        zhProperties.put(FONT_SIZE,"16");
        zhProperties.put(FONT_NAME,"Dialog");
        zhProperties.put(FONT_STYLE,Integer.toString(Font.BOLD));
        zhProperties.put(REQUIRES_FONT_CHOICE,"TRUE");
        _languageProperties.put("zh",zhProperties);
    }

    /**
     *
     * @return a HashMap full of HashMaps
     */
    protected static HashMap getLanguageProperties()
    {
        if (_languageProperties == null)
            initLanguageProperties();
        return _languageProperties;
    }

    /**
     * return a single language property
     * @param language
     * @param propertyName
     * @return
     */
    protected static String getLanguageProperty(String language, String propertyName)
    {
        return (String) ((HashMap)(getLanguageProperties().get(language))).get(propertyName);
    }

    /**
     * return a single language property
     * @param propertyName
     * @return
     */
    protected static String getLanguageProperty(String propertyName)
    {
        return getLanguageProperty(getLanguage(), propertyName);
    }

    /**
     * return a single languge property, as an int
     * @param language
     * @param propertyName
     * @return
     */
    protected static int getLanguagePropertyInt(String language, String propertyName)
    {
        int result = 0;
        try
        {
            result = Integer.parseInt(getLanguageProperty(language,propertyName));
        }
        catch(Exception e) {}
        return result;
    }

    protected static int getLanguagePropertyInt(String propertyName)
    {
        return getLanguagePropertyInt(getLanguage(), propertyName);
    }


    /**
     * return a single language property, as a double
     * @param language
     * @param propertyName
     * @return
     */
    protected static double getLanguagePropertyDouble(String language, String propertyName)
    {
        double result = 0;
        try
        {
            result = Double.parseDouble(getLanguageProperty(language,propertyName));
        }
        catch(Exception e) {}
        return result;
    }

    protected static double getLanguagePropertyDouble(String propertyName)
    {
        return getLanguagePropertyDouble(getLanguage(), propertyName);
    }

    /**
     * Create a locale from a string.  Used for loading preferences
     *
     * @param localeString of the format "<language_code>" or "<language_code>_<country_code>"
     * @return the appropriate locale
     */
    protected static Locale createLocale(String localeString)
    {
        String[] stringPieces = localeString.split("_");

        if (stringPieces.length == 1)
        {
            return new Locale(stringPieces[0]);
        }
        return new Locale(stringPieces[0], stringPieces[1]);
    }

    /**
     * Turn a Locale into a string that we can store as a preference
     * @param locale
     * @return a string that we can store as a preference
     */
    protected static String createLocaleString(Locale locale)
    {
        String result = locale.getLanguage();
        if (locale.getCountry() != null && locale.getCountry().length() > 0)
            result = result + "_" + locale.getCountry();
        return result;
    }



    /**
     * set the initial locale and font for the session, from the stored preference if one is available
     */
    public static void setInitialLocaleProperties()
    {
        String localeString = _prefs.get(_localePrefKey, DEFAULT_LOCALE_STRING);
        changeLocale(createLocale(localeString));
        String fontString = _prefs.get(getFontPrefKey(), null);
        if (fontString == null)
        {
            fontString = getWorkingFont().getName();
        }
        changeFont(fontString);
    }

    /**
     * Set the locale for the application.  Sets the default Java Locale, determines the
     * appropriate font, and updates all Swing components to use that font.  Also updates
     * the menu item text in the Language menu
     * @param locale
     */
    public static void changeLocale(Locale locale)
    {
        if (locale.equals(_locale))
            return;

        //if no fonts supporting this language, give the user a warning, don't change
        if (!doesLanguageHaveFonts(locale.getLanguage()))
        {
            ApplicationContext.infoMessage(TextProvider.getText("NO_FONTS_FOR_LANGUAGE_CHOSEN"));
            return;
        }

        _locale = locale;
        _prefs.put(_localePrefKey, createLocaleString(locale));

        Locale.setDefault(locale);
        TextProvider.reloadTextBundle();

        String fontString = _prefs.get(getFontPrefKey(), null);
        changeFont(fontString);

        setSwingDefaultFont(_font);

        //redraw the entire application, losing any data
        Application application = (Application) ApplicationContext.getImpl();
        if (application != null)
        {
            application.redrawWorkbench();
        }
/* it sure would be nice if we could update everything this way. However, this doesn't
 //give us a way to change all the text on the components to the right text for the new
 //Locale.
        WorkbenchFrame frame = (WorkbenchFrame) application.getFrame();
        if (frame != null)
        {
            frame.getSwingEngine().setLocale(getLocale());
            setFontForComponentTree(frame, font);
            recreateLanguageMenuItems(frame);
            //update all current UI components
            frame.update(frame.getGraphics());
        }
*/
    }



        /**
         * returns a swing engine with the appropriate locale set
          * @return a swing engine with the appropriate locale set
         */
    public static SwingEngine getSwingEngine(Object client)
    {
        SwingEngine swingEngine = new SwingEngine(client);
        swingEngine.setLocale(Localizer.getLocale());
        return swingEngine;
    }

    /**
     * Utility method to scale up the width of a dimension by a certain factor
     * @param dimension
     * @param widthFactor
     * @param heightFactor
     * @return
     */
    protected static Dimension adjustDimension(Dimension dimension, double widthFactor, double heightFactor)
    {
        int newHeight = (int) (dimension.getHeight() * heightFactor);
        int newWidth = (int) (dimension.getWidth() * widthFactor);

        return new Dimension(newWidth,newHeight);
    }

    /**
     * Special method used to resize the main window
     * @param dimension
     * @return
     */
    public static Dimension adjustContainerDimension(Dimension dimension)
    {
        return adjustDimension(dimension,getLanguagePropertyDouble(CONTAINER_WIDTH_SCALE),
                               getLanguagePropertyDouble(CONTAINER_HEIGHT_SCALE));
    }

    public static int adjustContainerWidth(int width)
    {
        return (int) ((double) width * getLanguagePropertyDouble(CONTAINER_WIDTH_SCALE));
    }
    public static int adjustContainerHeight(int height)
    {
        return (int) ((double) height * getLanguagePropertyDouble(CONTAINER_HEIGHT_SCALE));
    }

    /**
     * renders a swixml file
     * @param swixmlPath
     * @param client
     * @return the Container that holds everything created by the swixml file
     * @throws Exception
     */
    public static Container renderSwixml(String swixmlPath, Object client) throws Exception
    {
        Container content = getSwingEngine(client).render(swixmlPath);
        localizeComponentTree(content);
        return content;
    }

    /**
     * Apply localization rules to the component, depending on its type
     * @param theComponent
     */
    public static void localizeComponent(Component theComponent)
    {
        double labelScale = getLanguagePropertyDouble(LABEL_SCALE);
        double containerWidthScale = getLanguagePropertyDouble(CONTAINER_WIDTH_SCALE);
        double containerHeightScale = getLanguagePropertyDouble(CONTAINER_HEIGHT_SCALE);

        //if this language is at the same scale as English, do nothing
        if (labelScale == 1.0 && containerWidthScale == 1.0 && containerHeightScale == 1.0)
            return;
        //for labels
        if (theComponent instanceof JLabel)
        {
            theComponent.setPreferredSize(adjustDimension(theComponent.getPreferredSize(),
                                    labelScale, 1));
            theComponent.setMinimumSize(adjustDimension(theComponent.getMinimumSize(),
                                  labelScale, 1));
        }
        //for containers
        else if (theComponent instanceof JPanel || theComponent instanceof JSplitPane ||
//for some reason, vertical resizing of JTables causes their content rows not to render!
//            theComponent instanceof JTable ||
            theComponent instanceof JScrollPane ||
            theComponent instanceof JTabbedPane)
        {
            theComponent.setPreferredSize(adjustDimension(theComponent.getPreferredSize(),
                                    containerWidthScale, containerHeightScale));
            theComponent.setMinimumSize(adjustDimension(theComponent.getMinimumSize(),
                                  containerWidthScale, containerHeightScale));

        }

   }

    /**
     * "Localizes" a component tree by scaling up the dimensions of certain types of objects
     * in the tree in order to try to make things fit despite larger fonts and longer labels
     * @param parent
     */
    public static void localizeComponentTree(Component parent)
    {
        if (parent == null)
            return;

        localizeComponent(parent);

        //recurse
        try
        {
            Component[] children = ((Container) parent).getComponents();
            for (int i = 0; i < children.length; i++)
            {
                if (children[i] != null)
                    localizeComponentTree(children[i]);
            }
        }
        catch (Exception e)
        {
        }
    }



    /**
     * Blow away the existing language menu items and recreate them.
     * This is necessary when switching languages, because the text
     * on the menu items needs to change
     * @param frame
     */
    protected static void recreateLanguageMenuItems(JFrame frame)
    {
        JMenuBar menuBar = frame.getJMenuBar();
        //it would be better if we could reference this by name
        JMenu languageMenu = menuBar.getMenu(LANGUAGE_MENU_INDEX);
        languageMenu.removeAll();
        initializeLanguageMenu(languageMenu);
    }

    /**
     * add menu items for each available language
     *
     * @param languageMenu
     */
    public static void initializeLanguageMenu(JMenu languageMenu)
    {
        //if only one language supported, disable the menu and hide it
        if (supportedLocales.length < 2)
        {
            languageMenu.setEnabled(false);
            languageMenu.setVisible(false);
            return;
        }
        for (int i = 0; i < supportedLocales.length; i++)
        {
            Locale locale = supportedLocales[i];
            JMenuItem menuItem = new JMenuItem(locale.getDisplayName(_locale));
            LanguageActionListener listener = new LanguageActionListener(locale);
            menuItem.addActionListener(listener);
            languageMenu.add(menuItem);
        }
    }

    /**
     * get the current locale
     * @return the current Locale
     */
    public static Locale getLocale()
    {
        return _locale;
    }

    /**
     * return the current language
     * @return the current language
     */
    public static String getLanguage()
    {
        return _locale.getLanguage();
    }




    //Font-related methods

    /**
     * return a font preference key for the current language
     * @return
     */
    protected static String getFontPrefKey()
    {
        return _fontPrefKeyPrefix + getLanguage();
    }    


    /**
     * Set the default font for all Swing components.  This font will be used for all
     * components instantiated after this method is called
     * @param font
     */
    public static void setSwingDefaultFont(Font font)
    {
        _font = font;
        FontUIResource f = new javax.swing.plaf.FontUIResource(font);
        java.util.Enumeration keys = UIManager.getDefaults().keys();
        while (keys.hasMoreElements())
        {
            Object key = keys.nextElement();
            Object value = UIManager.get(key);
            if (value instanceof FontUIResource)
            {
                UIManager.put(key, f);
            }
        }
    }

    /**
     * Set the font for the application.  Sets the default Java Locale, determines the
     * appropriate font, and updates all Swing components to use that font.
     * If, however, the font can't handle the language, find one that can
     * @param fontName
     */
    public static void changeFont(String fontName)
    {
        if (fontName == null)
            fontName = getLanguageProperty(FONT_NAME);

        if ( _font != null && fontName.equals(_font.getName()))
        {
            return;
        }

        _font = createFont(fontName);

        //if this font can't handle the current language, don't use it
        if (!canHandleLanguage(getLanguage(),_font))
        {
            _font = getWorkingFont();
            String newFontName = _font.getName();
            System.err.println("Chosen font " + fontName + " can't handle the current language.  Switching to " + newFontName);
            fontName = newFontName;
        }


        _prefs.put(getFontPrefKey(), fontName);

        setSwingDefaultFont(_font);

        Application application = (Application) ApplicationContext.getImpl();
        WorkbenchFrame frame = (WorkbenchFrame) application.getFrame();
        if (frame != null)
        {
            getSwingEngine(frame).setLocale(getLocale());
            setFontForComponentTree(frame, _font);
            //update all current UI components
            frame.update(frame.getGraphics());
        }

    }

    /**
     * Find a working font for the Locale.  This requires special handling for some languages
     * @return the appropriate font for the locale
     */
    protected static Font getWorkingFont()
    {
        //this will work for languages where font doesn't matter so much
        Font result = createFont(getLanguageProperty(FONT_NAME));
        if ("TRUE".equals(getLanguageProperty(REQUIRES_FONT_CHOICE)))
        {
            //use the previously chosen font if it exists and can handle the language
            //Otherwise, find one that can, if there's one installed.
            if (!canHandleLanguage(getLanguage(),result))
            {
                //Default font can't handle the language.  Looking for other fonts
                Font alternateFont = pickFirstFontForLanguage();
                if (alternateFont != null)
                {
                    result = alternateFont;
                    _log.info("Found fonts for language " + getLanguage() + ". Using font " + alternateFont.getName());
                }
                else
                    ApplicationContext.infoMessage("No fonts found for language " + getLanguage() + ".  Proceeding with default font.");
            }
        }

        return result;
    }

    /**
     * Change the font for a tree of components rooted at parent
     * Not useful, since we don't have a way of changing the text, too
     * @param parent
     * @param font
     */
    protected static void setFontForComponentTree(Component parent, Font font)
    {
        parent.setFont(font);
        try
        {
            Component[] children = ((Container) parent).getComponents();
            for (int i = 0; i < children.length; i++)
            {
                if (children[i] != null)
                    setFontForComponentTree(children[i], font);
            }
        }
        catch (Exception e)
        {
        }
    }



    /**
     * Determines whether a given font can handle a given language
     * @param language
     * @param font
     * @return whether the font can handle the language
     */
    protected static boolean canHandleLanguage(String language, Font font)
    {
        return (font.canDisplayUpTo(getLanguageProperty(language,FONT_TEST_STRING))) == -1;
    }

    /**
     * Get all font names that support this language
     * @param language
     * @return
     */
    public static Vector getFontNamesForLanguage(String language)
    {
        Vector languageFonts = new Vector();

        Font[] allfonts = GraphicsEnvironment.getLocalGraphicsEnvironment().getAllFonts();
        int languageFontCount = 0;
        for (int j = 0; j < allfonts.length; j++)
        {
            if (canHandleLanguage(language,allfonts[j]))
            {
                languageFonts.add(allfonts[j].getName());
            }
        }
        return languageFonts;
    }

    /**
     * Is there at least one supporting font for the language?
     * @param language
     * @return
     */
    public static boolean doesLanguageHaveFonts(String language)
    {
        //always assume english has a supporting font
        if ("en".equals(language))
            return true;
        return (getFontNamesForLanguage(language).size() > 0);
    }

    /**
     * Search through all the fonts and pick the first one that can handle Chinese
     * @return the first font that can handle Chinese, or null if there is none
     */
    public static Font pickFirstFontForLanguage()
    {
        Font result = null;
        Vector languageFontNames = getFontNamesForLanguage(getLanguage());
        if (languageFontNames.size() > 0)
        {
            //pick the first Chinese-capable font, and size it appropriately
            result = new Font((String) languageFontNames.elementAt(0),
                               getLanguagePropertyInt(getLanguage(), FONT_STYLE),
                               getLanguagePropertyInt(getLanguage(), FONT_SIZE));
        }
        return result;
    }

    /**
     * Blow away the existing font menu items and recreate them.
     * This is necessary when switching languages, because the available fonts will change
     * @param frame
     */
    protected static void recreateFontMenuItems(JFrame frame)
    {
        JMenuBar menuBar = frame.getJMenuBar();
        //it would be better if we could reference this by name
        JMenu fontMenu = menuBar.getMenu(FONT_MENU_INDEX);
        fontMenu.removeAll();
        initializeFontMenu(fontMenu);
    }

    /**
     * add menu items for each available font for the current language, if the language is picky
     * about fonts
     *
     * @param fontMenu
     */
    public static void initializeFontMenu(JMenu fontMenu)
    {
        //if only current language does not require font choice, disable the menu and hide it
        if ("FALSE".equals(getLanguageProperty(REQUIRES_FONT_CHOICE)))
        {
            fontMenu.setEnabled(false);
            fontMenu.setVisible(false);
            return;
        }
        Vector fontNames = getFontNamesForLanguage(getLanguage());
        if (fontNames.size() < 1)
        {
            fontMenu.setEnabled(false);
            fontMenu.setVisible(false);
            return;
        }

        for (int i = 0; i < fontNames.size(); i++)
        {
            String fontName = (String) fontNames.elementAt(i);
            JMenuItem menuItem = new JMenuItem(fontName);
            menuItem.setFont(createFont(fontName));
            FontActionListener listener = new FontActionListener(fontName);
            menuItem.addActionListener(listener);
            fontMenu.add(menuItem);
        }
    }

    protected static Font createFont(String fontName)
    {
        return new Font(fontName,
                        getLanguagePropertyInt(getLanguage(), FONT_STYLE),
                        getLanguagePropertyInt(getLanguage(), FONT_SIZE));
    }


    /**
     * actionlistener class for language menu items
     */
    protected static class LanguageActionListener implements ActionListener
    {
        Locale _locale;

        public LanguageActionListener(Locale locale)
        {
            _locale = locale;
        }

        public void actionPerformed(ActionEvent event)
        {
            Localizer.changeLocale(_locale);
        }
    }

    /**
     * actionlistener class for font menu items
     */
    protected static class FontActionListener implements ActionListener
    {
        String _fontName;

        public FontActionListener(String fontName)
        {
            _fontName = fontName;
        }

        public void actionPerformed(ActionEvent event)
        {
            Localizer.changeFont(_fontName);
        }
    }







}
