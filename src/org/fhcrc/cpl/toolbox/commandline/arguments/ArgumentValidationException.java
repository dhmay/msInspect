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

package org.fhcrc.cpl.toolbox.commandline.arguments;

/**
 * This class purely serves as a way to create a known type of Exception from a
 * non-Exception event (like a String not conforming to a required format)
 */
public class ArgumentValidationException extends Exception
{
    protected Exception mNestedException = null;

    //indicates whether stack trace should be shown.  Default, no stack trace
    protected boolean mShouldShowStackTrace = false;

    public ArgumentValidationException(String message)
    {
        super(message);
    }

    public ArgumentValidationException(String message, boolean shouldShowStackTrace)
    {
        super(message);
        mShouldShowStackTrace = shouldShowStackTrace;
    }

    public ArgumentValidationException(String message, Exception e)
    {
        super(message + "\n" + e.getMessage());
        mNestedException = e;
        mShouldShowStackTrace = true;
    }

    public ArgumentValidationException(Exception e)
    {
        super(e.getMessage());
        mNestedException = e;
        mShouldShowStackTrace = true;
    }

    /**
     * Indicates whether stack trace should be shown.  If we're getting an exception
     * we don't understand, best to show stack trace.  If, instead, it's something like
     * a missing file, better to give a meaningful message and leave it at that.
     * @return
     */
    public boolean shouldShowStackTrace()
    {
        return mShouldShowStackTrace;
    }

    public void setShouldShowStackTrace(boolean shouldShowStackTrace)
    {        
        mShouldShowStackTrace = shouldShowStackTrace;
    }
}
