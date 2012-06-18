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

package org.swixml.converters;

import org.jdom.Attribute;
import org.jdom.DataConversionException;
import org.swixml.Converter;
import org.swixml.Localizer;
import org.swixml.Parser;

import javax.swing.*;
import javax.swing.tree.TreeSelectionModel;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * The <code>PrimitiveConverter</code> class defines a converter that creates primitive objects (wrapper),
 * based on a provided input Class and String.
 *
 * @author <a href="mailto:wolf@paulus.com">Wolf Paulus</a>
 * @version $Revision: 1.1 $
 * @see org.swixml.ConverterLibrary

 */

public class PrimitiveConverter implements Converter, SwingConstants, KeyEvent, InputEvent {

  /** converter's return type */
  public static final Class TEMPLATE = Object.class;
  /** map contains all constant provider types */
  private static Map dictionaries = new HashMap();
  /**
   * Static Initializer, setting up the initial constant providers
   */
  static {

    PrimitiveConverter.addConstantProvider( JTabbedPane.class );
    PrimitiveConverter.addConstantProvider( JScrollPane.class );
    PrimitiveConverter.addConstantProvider( GridBagConstraints.class );
    PrimitiveConverter.addConstantProvider( FlowLayout.class );
    PrimitiveConverter.addConstantProvider( ListSelectionModel.class );
    PrimitiveConverter.addConstantProvider( TreeSelectionModel.class );
    PrimitiveConverter.addConstantProvider( WindowConstants.class );
    PrimitiveConverter.addConstantProvider( JFrame.class );
//dhmay adding, custom for msInspect
    PrimitiveConverter.addConstantProvider( JTextField.class );
    PrimitiveConverter.addConstantProvider( JSplitPane.class );
    PrimitiveConverter.addConstantProvider( SwingConstants.class );
//end dhmay additions


  };
  /**
   * Converts String into java primitive type
   * @param type <code>Class</code> target type
   * @param attr <code>Attribute</code> value field needs to provide convertable String
   * @return <code>Object</code> primitive wrapped into wrapper object
   */
  public static Object conv( final Class type, final Attribute attr,  final Localizer localizer ) throws Exception {
    Attribute a = (Attribute) attr.clone();
    Object obj = null;
    if ( Parser.LOCALIZED_ATTRIBUTES.contains( a.getName().toLowerCase() ))
      if (a.getAttributeType() == Attribute.CDATA_TYPE )
         a.setValue( localizer.getString( a.getValue() ));

    try {
      if (boolean.class.equals( type )) {
        obj = new Boolean( a.getBooleanValue() );
      } else if (int.class.equals( type )) {
        obj = new Integer( a.getIntValue() );
      } else if (long.class.equals( type )) {
        obj = new Long( a.getLongValue() );
      } else if (float.class.equals( type )) {
        obj = new Float( a.getFloatValue() );
      } else if (double.class.equals( type )) {
        obj = new Double( a.getDoubleValue() );
      }
    } catch (DataConversionException e) {
    } finally {
      if (obj==null) {
        try {
          String s = a.getValue();
          int k = s.indexOf( '.' ) - 1;
          Class pp = (Class) dictionaries.get( s.substring( 0, s.indexOf( '.' ) ) );
          obj = pp.getField( s.substring( k + 2 ) ).get( pp );
        } catch (Exception ex) {
          //
          //  Try to find the given value as a Constant in SwingConstants
          //
          obj = PrimitiveConverter.class.getField( a.getValue() ).get( PrimitiveConverter.class );
        }
      }
    }

    return obj;
  }

  /**
   * Converts String into java primitive type
   * @param type <code>Class</code> target type
   * @param attr <code>Attribute</code> value field needs to provide convertable String
   * @return <code>Object</code> primitive wrapped into wrapper object
   * @throws Exception
   */
  public Object convert( final Class type, final Attribute attr, final Localizer localizer ) throws Exception {
    return PrimitiveConverter.conv( type, attr, localizer );
  }

  /**
   * A <code>Converters</code> conversTo method informs about the Class type the converter
   * is returning when its <code>convert</code> method is called
   * @return <code>Class</code> - the Class the converter is returning when its convert method is called
   */
  public Class convertsTo() {
    return TEMPLATE;
  }

  /**
   * Adds a new class or interface to the dictionary of primitive providers.
   * @param type <code>Class</code> providing primitive constants / public (final) fields
   */
  public static void addConstantProvider(final Class type) {
    String shortName = type.getName().substring(type.getName().lastIndexOf( '.' ) + 1 );
    dictionaries.put( shortName, type );
  }
}
