package org.systemsbiology.jrap;

/**
 * Copyright (C) 2004 Fred Hutchinson Cancer Research Center. All Rights Reserved.
 * User: migra
 * Date: Jul 27, 2004
 * Time: 3:42:06 PM
 *
 *
 * Non-thread-safe buffer meant to behave like a StringBuffer.
 */
public class ByteAppender
{
    byte[] value;
    int count = 0;

    public ByteAppender()
    {
        value = new byte[1024];
    }

    public ByteAppender(int len)
    {
        value = new byte[len];
    }

    public byte[] getBuffer()
    {
        return value;
    }

    public int getCount()
    {
        return count;
    }

    public void reset()
    {
        count = 0;
    }

    public void ensureCapacity(int minimumCapacity)
    {
        if (minimumCapacity > value.length)
        {
            expandCapacity(minimumCapacity);
        }
    }

    private void expandCapacity(int minimumCapacity)
    {
        int newCapacity = (value.length + 1) * 2;
        if (newCapacity < 0)
        {
            newCapacity = Integer.MAX_VALUE;
        }
        else if (minimumCapacity > newCapacity)
        {
            newCapacity = minimumCapacity;
        }

        byte newValue[] = new byte[newCapacity];
        System.arraycopy(value, 0, newValue, 0, count);
        value = newValue;
    }

    public ByteAppender append(byte bytes[], int offset, int len) {
        int newcount = count + len;
        if (newcount > value.length)
            expandCapacity(newcount);
        System.arraycopy(bytes, offset, value, count, len);
        count = newcount;
        return this;
    }

    public ByteAppender appendCharsAsBytes(char chars[], int offset, int len) {
        int newcount = count + len;
        if (newcount > value.length)
            expandCapacity(newcount);

        int end = offset + len;
        int dst = count;
        for (int i = offset; i < end; i++)
            value[dst++] = (byte) chars[i];

        count = newcount;
        return this;
    }

}
