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
package org.fhcrc.cpl.viewer.util;

import org.apache.commons.beanutils.ConversionException;
import org.apache.commons.beanutils.ConvertUtils;
import org.apache.commons.beanutils.Converter;
import org.apache.commons.beanutils.converters.*;
import org.fhcrc.cpl.toolbox.AbstractConvertHelper;

import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

public class ConvertHelper extends AbstractConvertHelper
{
    private static ConvertHelper _myInstance = null;

    public static void registerHelpers()
    {
        // Should register only once
        assert null == _myInstance;

        _myInstance = new ConvertHelper();
        _myInstance.register();
    }


    private ConvertHelper()
    {
    }


    /*
     * Register shared converters, then msInspect converters
     * TODO: Should merge cpl and msInspect converters
     */
    protected void register()
    {
        super.register();

        ConvertUtils.register(new MyBooleanConverter(), Boolean.TYPE);
        ConvertUtils.register(new NullSafeConverter(new MyBooleanConverter()), Boolean.class);
        ConvertUtils.register(new NullSafeConverter(new ByteArrayConverter()), byte[].class);
        ConvertUtils.register(new PercentWrapper(new DoubleConverter()), Double.TYPE);
        ConvertUtils.register(new NullSafeConverter(new PercentWrapper(new DoubleConverter())), Double.class);
        ConvertUtils.register(new PercentWrapper(new FloatConverter()), Float.TYPE);
        ConvertUtils.register(new NullSafeConverter(new PercentWrapper(new FloatConverter())), Float.class);
        ConvertUtils.register(new ChargeWrapper(new IntegerConverter()), Integer.TYPE);
        ConvertUtils.register(new NullSafeConverter(new ChargeWrapper(new IntegerConverter())), Integer.class);
        ConvertUtils.register(new NullSafeConverter(new DateFriendlyStringConverter()), String.class);
        ConvertUtils.register(new LenientTimestampConverter(), java.sql.Timestamp.class);
        ConvertUtils.register(new LenientDateConverter(), java.util.Date.class);
    }


    public static class PercentWrapper implements Converter
    {
        Converter conv;

        PercentWrapper(Converter converter)
        {
            this.conv = converter;
        }

        public Object convert(Class aClass, Object o)
        {
            if (o instanceof String)
            {
                String s = ((String)o).trim();
                if (s.endsWith("%"))
                {
                    s = s.substring(0, s.length()-1);
                    Object v = conv.convert(aClass, s);
                    if (v instanceof Double)
                        return new Double(((Double)v).doubleValue() / 100);
                    else
                        return new Float(((Float)v).floatValue() / 100);
                }
            }
            return conv.convert(aClass, o);
        }
    }


    public static class ChargeWrapper implements Converter
    {
        Converter conv;

        ChargeWrapper(Converter converter)
        {
            this.conv = converter;
        }

        public Object convert(Class aClass, Object o)
        {
            if (o instanceof String)
            {
                String s = ((String)o).trim();
                if (s.endsWith("+"))
                    o = s.substring(0, s.length()-1);
            }
            return conv.convert(aClass, o);
        }
    }


    public static class MyBooleanConverter implements Converter
    {
        static BooleanConverter defaultConverter = new BooleanConverter();

        public Object convert(Class aClass, Object o)
        {
            if (o instanceof Boolean)
                return o;
            if (o instanceof Integer)
                return ((Integer) o).intValue() != 0 ? Boolean.TRUE : Boolean.FALSE;
            if (o instanceof String)
            {
                String s = (String)o;
                if (s.equals("0") || s.equalsIgnoreCase("false"))
                    return Boolean.FALSE;
                if (s.equals("1") || s.equalsIgnoreCase("true"))
                    return Boolean.TRUE;
            }
            return defaultConverter.convert(aClass, o);
        }
    }


    public static class LenientTimestampConverter implements Converter
    {
        private LenientDateConverter _dateConverter = new LenientDateConverter();
        public Object convert(Class clss, Object o)
        {
            if (null == o)
                return null;

            if (o instanceof Timestamp)
                return o;

            return new Timestamp(((Date) _dateConverter.convert(Date.class, o)).getTime());
        }
    }

    /**
     * This format accepts dates in the form MM/dd/yy
     */
    public static class LenientDateConverter implements Converter
    {
        private static SimpleDateFormat _fullReadableFormat = null;
        private static SimpleDateFormat _userFormat = null;
        static
        {
            _fullReadableFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            _fullReadableFormat.setLenient(true);
            _userFormat = new SimpleDateFormat("MM/dd/yy hh:mm a");
            _userFormat.setLenient(true);
        }

        public Object convert(Class clss, Object o)
        {
            if (null == o)
                return null;


            if (o instanceof java.util.Date)
                return o;

            Date dt = null;
            //Try a few different things
            try
            {
                dt = DateFormat.getDateInstance().parse(o.toString());
            }
            catch (Exception x){}

            if (null == dt)
                try
                {
                    dt =(Date)  _fullReadableFormat.parseObject(o.toString());
                }
                catch (Exception x) {}

            if (null == dt)
                try
                {
                    dt = (Date) _userFormat.parseObject(o.toString());
                }
                catch (Exception x) {}

            if (null == dt)
                throw new ConversionException("Couldn't convert date '" + o.toString() + "'");

            return dt;
        }
    }

    public static class DateFriendlyStringConverter implements Converter
    {
        private static Format _dateOnlyFormat = SimpleDateFormat.getDateInstance(DateFormat.SHORT);
        private static Format _fullReadableFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        private static Converter _stringConverter = new StringConverter();

        public Object convert(Class clss, Object o)
        {
            if (o instanceof String)
                return o;

            if (o instanceof Date)
            {
                Calendar cal = Calendar.getInstance();
                cal.setTime((Date) o);
                cal.clear(Calendar.HOUR);
                cal.clear(Calendar.MINUTE);
                cal.clear(Calendar.SECOND);
                cal.clear(Calendar.MILLISECOND);
                if(cal.getTime().equals(o))
                    return _dateOnlyFormat.format(o);
                else
                    return _fullReadableFormat.format(o);
            }

            return _stringConverter.convert(String.class, o);
        }
    }
}
