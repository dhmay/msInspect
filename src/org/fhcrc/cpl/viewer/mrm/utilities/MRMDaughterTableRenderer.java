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

import org.fhcrc.cpl.viewer.mrm.MRMDaughter;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Sep 19, 2007
 * Time: 1:26:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class MRMDaughterTableRenderer extends JLabel
    implements TableCellRenderer {
   Border unselectedBorder = null;
   Border selectedBorder = null;
   boolean isBordered = true;

   public MRMDaughterTableRenderer(boolean isBordered) {
      super();
      this.isBordered = isBordered;
      setOpaque(true); //MUST do this for background to show up.
   }

   public Component getTableCellRendererComponent(
                             JTable table, Object daughter,
                             boolean isSelected, boolean hasFocus,
                             int row, int column)
   {
      this.setText("");
      this.setHorizontalAlignment(JLabel.RIGHT);
      Color newColor = Color.LIGHT_GRAY;
      if(daughter != null &&  (Color)((MRMDaughter)daughter).getGraphColor() != null)
         newColor = (Color)((MRMDaughter)daughter).getGraphColor();
      setForeground(newColor);
      if (isBordered) {
         if (isSelected) {
            if (selectedBorder == null) {
               selectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                               table.getSelectionBackground());
            }
            setBorder(selectedBorder);
         } else {
            if (unselectedBorder == null)
            {
               unselectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                      table.getBackground());
            }
            setBorder(unselectedBorder);
         }
      }
      super.setBackground(UIManager.getColor("Table.focusCellBackground"));
      if(daughter != null) this.setText(((MRMDaughter)daughter).getName().split("/")[1]);
      return this;
  }
}
