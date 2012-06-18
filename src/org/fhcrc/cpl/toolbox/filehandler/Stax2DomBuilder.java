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

package org.fhcrc.cpl.toolbox.filehandler;

import org.w3c.dom.*;

import javax.xml.parsers.*;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.Reader;

/**
 * Builds a DOM {@link org.w3c.dom.Document} using a
 * {@link javax.xml.stream.XMLStreamReader}.
 *
 * @version $Revision: 1.00 $, $Date: 2004/12/11 00:00:00 $
 * @author  Tatu Saloranta
 *
 * Major revisions by Damon May -- this class searches for a start tag with a given name and
 * produces the DOM tree that represents that node.
 */
public class Stax2DomBuilder
{
    // // // Configuration settings:

    /**
     * Whether ignorable white space should be ignored, ie not added
     * in the resulting JDOM tree. If true, it will be ignored; if false,
     * it will be added in the tree. Default value if false.
     */
    protected boolean mCfgIgnoreWs = false;

    protected boolean mNsAware = true;

    // // Trivial caching...

    protected String mLastPrefix = null;
    protected String mLastLocalName = null;
    protected String mLastQName = null;

    //Default document for use in generating DOM nodes
    protected Document _defaultDocument = null;

    protected Reader _fileReader = null;

    // streamreader persisted for multiple sequential checks of the same file
    protected XMLStreamReader _docStreamReader = null;

    /**
     * Creates a new XMLStreamReader for the file and prepares to use it.  This should
     * be the only constructor exposed, since _docStreamReader is required by the parsing methods
     * @param file
     * @throws FileNotFoundException
     * @throws XMLStreamException
     */
    public Stax2DomBuilder(File file) throws FileNotFoundException, XMLStreamException
    {
        _fileReader = new java.io.FileReader(file);
        javax.xml.stream.XMLInputFactory f = javax.xml.stream.XMLInputFactory.newInstance();
        _docStreamReader = f.createXMLStreamReader(_fileReader);
    }

    /**
     * prepares to read from the given streamreader
     * @param xmlStreamReader
     */
    public Stax2DomBuilder(XMLStreamReader xmlStreamReader)
    {
        _docStreamReader = xmlStreamReader;
    }

    /**
     * Method used to change whether the build methods will add ignorable
     * (element) white space in the DOM tree or not.
     *<p>
     * Whether all-whitespace text segment is ignorable white space or
     * not is based on DTD read in, as per XML specifications (white space
     * is only significant in mixed content or pure text elements).
     */
    public void setIgnoreWhitespace(boolean state) {
        mCfgIgnoreWs = state;
    }

    /**
     * Creates a default Document for node creation, if necessary, and returns it
     * @return
     */
    protected Document getDefaultDocument()
    {
        if (_defaultDocument == null)
        {
            try
            {
                _defaultDocument = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
            }
            catch (Exception e) { e.printStackTrace(System.err); }
        }
        return _defaultDocument;
    }

    /**
     * This cover method will create a {@link org.w3c.dom.Document} instance using
     * the default JAXP mechanism and
     * populate using the given StAX stream reader.
     *
     * @return <code>Document</code> - DOM document object.
     * @throws XMLStreamException If the reader threw such exception (to
     *   indicate a parsing or I/O problem)
     */
    public Node findTreeForName(String nodeName, String endNodeNameToStop)
        throws ParserConfigurationException, XMLStreamException
    {
        return findTreeForName(getDefaultDocument(), nodeName, endNodeNameToStop);
//        return findTreeForName(DocumentBuilderFactory.newInstance().newDocumentBuilder(), nodeName);
    }

    /**
     * This cover method uses the passed-in DocumentBuilder to build a new doc
     * @param docbuilder
     * @param nodeName
     * @return
     * @throws XMLStreamException
     */
    public Node findTreeForName(DocumentBuilder docbuilder, String nodeName, String endNodeNameToStop)
        throws XMLStreamException
    {
        Document doc = docbuilder.newDocument();
        return findTreeForName(doc, nodeName, endNodeNameToStop);
    }

    /**
     * This method uses the class <code>XMLStreamReader</code> and builds up
     * a DOM tree for the first occurrence of a tag with name nodeName, using the
     * passed-in <code>Document</code>. Recursion has been eliminated by using nodes'
     * parent/child relationship; this improves performance somewhat
     * (classic recursion-by-iteration-and-explicit stack transformation)
     *
     * @param doc DOM <code>Document</code> for use in constructing Nodes
     * @param nodeName the name of the tag to match
     * @param endElementNameToStop When we see an end tag that matches this name,
     * if we're not currently building a tree we'll stop.  If null, search whole doc
     * @return a Node containing the tree representing the first occurrence
     * of a tag with name nodeName in the document.  If none found, return null
     */
    public Node findTreeForName(Document doc, String nodeName, String endElementNameToStop)
        throws XMLStreamException
    {
        checkReaderSettings(_docStreamReader);

        //Top level node.  In XmlBeans 2.1, they introduced a check to make sure you don't
        //have more than one node at the top level.  Their check is a little strange in that
        //it always seems to stop you from appending a node to a document.  So we can't just
        //add our top node directly to the document, we have to create this one level of indirection.
        Node top = doc.createElement("dummy");
        Node current = top;

        //this will hold the root of the tree we want to keep, if we find it
        Node matchingTreeRoot = null;

        //reflects whether we're inside the tree we actually want to return.
        //If so, great, build everything normally.  If not, refuse to do anything
        //for any element.
        boolean insideTree = false;

        main_loop:

        while (true) {
            int evtType = _docStreamReader.next();
            Node child = null;

            switch (evtType) {
            case XMLStreamConstants.CDATA:
                if (insideTree)
                    child = doc.createCDATASection(_docStreamReader.getText());
                break;

            case XMLStreamConstants.SPACE:
                if (mCfgIgnoreWs) {
                    continue main_loop;
                }
                /* Oh great. DOM is brain-dead in that ignorable white space
                 * can not be added, even though it is legal, and often
                 * reported by StAX/SAX impls...
                 */
                if (current == top) { // better just ignore, thus...
                    continue;
                }
                // fall through

            case XMLStreamConstants.CHARACTERS:
                if (insideTree)
                    child = doc.createTextNode(_docStreamReader.getText());
                break;

            case XMLStreamConstants.COMMENT:
                if (insideTree)
                    child = doc.createComment(_docStreamReader.getText());
                break;

            //probably shouldn't get here
            case XMLStreamConstants.END_DOCUMENT:
            {
                break main_loop;
            }

            //check whether we've found the end tag of the root node of our tree, or
            //the end tag that indicates we should stop.
            //if so, we're done
            case XMLStreamConstants.END_ELEMENT:
                current = current.getParentNode();
                if (current == null) {
                    current = top;
                }
                //rely on the fact that our root tree node is assigned directly to the doc,
                //so stop if the parent is the doc itself.
                //Also stop if we're NOT building a tree right now and we hit our sentinel
                String thisEndNodeName = _docStreamReader.getLocalName();
                if ((insideTree && current == top) ||
                    (!insideTree && thisEndNodeName.equals(endElementNameToStop)))
                {
                    break main_loop;
                }
                continue main_loop;

            case XMLStreamConstants.ENTITY_DECLARATION:
            case XMLStreamConstants.NOTATION_DECLARATION:
                /* Shouldn't really get these, but maybe some stream readers
                 * do provide the info. If so, better ignore it -- DTD event
                 * should have most/all we need.
                 */
                continue main_loop;

            case XMLStreamConstants.ENTITY_REFERENCE:
                if (insideTree)
                    child = doc.createEntityReference(_docStreamReader.getLocalName());
                break;

            case XMLStreamConstants.PROCESSING_INSTRUCTION:
                if (insideTree)
                    child = doc.createProcessingInstruction(_docStreamReader.getPITarget(),
                                                            _docStreamReader.getPIData());
                break;

            case XMLStreamConstants.START_ELEMENT:
                // Ok, need to add a new element...
                {
                    String ln = _docStreamReader.getLocalName();
                    Element newElem;

                    boolean thisIsRootNode = false;

                    if (!insideTree)
                    {
                        //check to see if we should start the tree up now
                        if (nodeName.equals(ln))
                        {
                            insideTree = true;
                            //indicate that, once the node is built, it should be saved
                            //as the root of our tree
                            thisIsRootNode = true;
                        }
                    }

                    //if we're not inside the tree, don't bother
                    if (!insideTree)
                        continue main_loop;
                    if (mNsAware) {

                        String elemPrefix = _docStreamReader.getPrefix();

                        // Doh, DOM requires a silly qualified name...
                        if (elemPrefix != null && elemPrefix.length() > 0) {
                            newElem = doc.createElementNS(_docStreamReader.getNamespaceURI(),
                                                          getQualified(elemPrefix, ln));
                        } else {
                            newElem = doc.createElementNS(_docStreamReader.getNamespaceURI(), ln);
                        }

                    } else { // if non-ns-aware, things are simpler:
                        newElem = doc.createElement(ln);
                    }

                    if (thisIsRootNode)
                    {
                        matchingTreeRoot = newElem;
                    }

                    /* No need to check namespace bindings, unlikes with some
                     * other frameworks (JDOM)
                     */

                    // And then the attributes:
                    for (int i = 0, len = _docStreamReader.getAttributeCount(); i < len; ++i) {
                        ln = _docStreamReader.getAttributeLocalName(i);
                        if (mNsAware) {
                            String prefix = _docStreamReader.getAttributePrefix(i);
                            if (prefix != null && prefix.length() > 0) {
                                ln = getQualified(prefix, ln);
                            }
                            Attr attr = doc.createAttributeNS(_docStreamReader.getAttributeNamespace(i), ln);
                            attr.setValue(_docStreamReader.getAttributeValue(i));
                            newElem.setAttributeNodeNS(attr);
                        } else {
                            Attr attr = doc.createAttribute(ln);
                            attr.setValue(_docStreamReader.getAttributeValue(i));
                            newElem.setAttributeNode(attr);
                        }

                    }
                    // And then 'push' new element...
                    current.appendChild(newElem);
                    current = newElem;

                    continue main_loop;
                }

            case XMLStreamConstants.START_DOCUMENT:
                /* This should only be received at the beginning of document...
                 * so, should we indicate the problem or not?
                 */
                /* For now, let it pass: maybe some (broken) readers pass
                 * that info as first event in beginning of doc?
                 */
                continue main_loop;

            case XMLStreamConstants.DTD:
                /* !!! Note: StAX does not expose enough information about
                 *  doctype declaration (specifically, public and system id!);
                 *  (altough StAX2 would...)
                 *
                 * Worse, DOM1/2 do not specify a way to create the DocType
                 * node, even if StAX provided it. This is pretty silly,
                 * all in all.
                 */
                continue main_loop;

            // Should never get these, from a stream reader:

            /* (commented out entries are just FYI; default catches
            * them all)
            */

            //case XMLStreamConstants.ATTRIBUTE:
            //case XMLStreamConstants.NAMESPACE:
            default:
                throw new XMLStreamException("Unrecognized iterator event type: "+_docStreamReader.getEventType()+"; should not receive such types (broken stream reader?)");
            }

            if (child != null && insideTree) {
                current.appendChild(child);
            }
        }
        return matchingTreeRoot;
    }

    // // // Overridable helper methods:

    protected String getQualified(String prefix, String localName)
    {
        /* This mostly/only helps with empty/text-only elements...
         * might make sense to do 'real' caching...
         */
        if (localName == mLastLocalName &&
            prefix == mLastPrefix) {
            return mLastQName;
        }
        String qn = prefix + ":" + localName;
        mLastQName = qn;
        return qn;
    }

    protected void checkReaderSettings(XMLStreamReader r)
        throws XMLStreamException
    {
        Object o = r.getProperty(XMLInputFactory.IS_NAMESPACE_AWARE);
        /* StAX defaults to namespace aware, so let's use similar
         * logics (although all compliant implementations really should
         * return a valid value)
         */
        if ((o instanceof Boolean) && !((Boolean) o).booleanValue()) {
            mNsAware = false;
        } else {
            mNsAware = true;
        }
    }

    /**
     * Free up resources, in this case just the streamreader
     */
    public void dispose()
    {
        try
        {
            if (_docStreamReader != null)
                _docStreamReader.close();
        }
        catch (Exception e) {}
        try
        {
            if (_fileReader != null)
                _fileReader.close();
        }
        catch (Exception e) {}
    }
}
