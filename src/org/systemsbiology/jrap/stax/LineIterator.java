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

package org.systemsbiology.jrap.stax;

import java.util.Iterator;
import java.nio.ByteBuffer;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Nov 16, 2009
 * Time: 2:34:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class LineIterator implements Iterator{

    private ByteBufferIterator bbi=null;
    private ByteBuffer bb=null;

    public ByteBufferIterator getBbi() {
        return bbi;
    }

    public void setBbi(ByteBufferIterator bbi) {
        this.bbi = bbi;
    }

    public ByteBuffer getBb() {
        return bb;
    }

    public long getFilePos() {
        return filePos;
    }

    public StringBuilder getCurLine() {
        return curLine;
    }

    public int getLineNum() {
        return lineNum;
    }


    private long filePos = 0;
    private StringBuilder curLine = new StringBuilder();
    private int lineNum=0;

    public LineIterator(ByteBufferIterator bbit) {
       bbi = bbit;
    }

    public boolean hasNext() {
        return bbi.hasNext() || (bb != null && bb.hasRemaining()) ;
    }

    //This code is a little iffy.  If a \r\n pair (or \n\n etc) straddle a ByteBuffer boundary
    // it will not work.  Also, line numbers will be wrong for \n\n.
    public StringBuilder next() {
        curLine.setLength(0);
        if(bb==null && bbi.hasNext()) {
            bb = bbi.next();
        }
        lineNum++;
        filePos = bbi.getFilePos()-bb.remaining();
        for(;;) {
           if(!bb.hasRemaining()) {
              if(!bbi.hasNext()) break;
              bb = bbi.next();
           }
           byte curByte = bb.get();
           if(curByte != '\n' && curByte != '\r') {
               curLine.append((char)curByte);
           } else {
              if(bb.hasRemaining()) {
                 byte nextByte = bb.get();
                 if(nextByte != '\n' && nextByte != '\r') bb.position(bb.position()-1);
              }
              break;
           }
        }
        return curLine;
    }

    public void remove() {}

}
