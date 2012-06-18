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
package org.fhcrc.cpl.viewer.mrm.utilities;

import org.fhcrc.cpl.viewer.gui.MRMDialog;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Sep 19, 2007
 * Time: 1:30:24 PM
 * To change this template use File | Settings | File Templates.
 */
   public class MRMNumberTableRenderer extends DefaultTableCellRenderer {
       protected Border noFocusBorder = new EmptyBorder(1, 1, 1, 1);
       private final Border SAFE_NO_FOCUS_BORDER = new EmptyBorder(1, 1, 1, 1);
       private Border getNoFocusBorder() {
           if (System.getSecurityManager() != null) {
               return SAFE_NO_FOCUS_BORDER;
           } else {
               return noFocusBorder;
           }
       }

       public Component getTableCellRendererComponent(JTable table, Object number,
                              boolean isSelected, boolean hasFocus, int row, int column) {
        this.setHorizontalAlignment(JLabel.RIGHT);
        if(number != null && number instanceof Float) {
            this.setText(MRMDialog.peaksData.floatFormat.format(number));
        }
        if(number != null && number instanceof Integer) {
            this.setText(((Integer)number).toString());
        }
        if(number == null || !(number instanceof Number)) {
            this.setText("");
        }

        if (isSelected) {
           super.setForeground(table.getSelectionForeground());
           super.setBackground(table.getSelectionBackground());
        }
        else {
            super.setForeground(table.getForeground());
            super.setBackground(table.getBackground());
        }

        setFont(table.getFont());

        if (hasFocus) {
                Border border = null;
                if (isSelected) {
                    border = UIManager.getBorder("Table.focusSelectedCellHighlightBorder");
                }
                if (border == null) {
                    border = UIManager.getBorder("Table.focusCellHighlightBorder");
                }
                setBorder(border);

            if (!isSelected && table.isCellEditable(row, column)) {
                    Color col;
                    col = UIManager.getColor("Table.focusCellForeground");
                    if (col != null) {
                        super.setForeground(col);
                    }
                    col = UIManager.getColor("Table.focusCellBackground");
                    if (col != null) {
                        super.setBackground(col);
                    }
            }
        } else {
                setBorder(getNoFocusBorder());
        }


        return this;
        }
   }