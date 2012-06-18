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

package org.fhcrc.cpl.toolbox.commandline;

import java.io.PrintWriter;
import java.io.PrintStream;

/**
 * This class is used to communicate that execution Did Not Go Well, and 
 * information about why
 */
public class CommandLineModuleExecutionException extends Exception
{
    protected Exception mNestedException = null;

    //indicates whether stack trace should be shown.  Default, no stack trace
    protected boolean mShouldShowStackTrace = false;

    public CommandLineModuleExecutionException(String message)
    {
        super(message);
    }

    public CommandLineModuleExecutionException(Exception e)
    {
        super(e.getMessage());
        mNestedException = e;
        mShouldShowStackTrace = true;
    }

    public CommandLineModuleExecutionException(String message, Exception e)
    {
        super(message + "\n" + e.getMessage());
        mNestedException = e;
        mShouldShowStackTrace = true;
    }

    public CommandLineModuleExecutionException(String message, boolean shouldShowStackTrace)
    {
        super(message);
        mShouldShowStackTrace = shouldShowStackTrace;
    }

    public void printStackTrace(PrintWriter pw)
    {
        if (mNestedException != null)
            mNestedException.printStackTrace(pw);
    }

    public void printStackTrace(PrintStream out)
    {
        if (mNestedException != null)
            mNestedException.printStackTrace(out);
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

    public Exception getNestedException()
    {
        return mNestedException;
    }
}
