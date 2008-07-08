package org.fhcrc.cpl.viewer.mrm.utilities;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.plaf.UIResource;
import javax.swing.table.TableCellRenderer;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Sep 19, 2007
 * Time: 11:36:27 AM
 * To change this template use File | Settings | File Templates.
 */

public class MRMBooleanRenderer extends JCheckBox implements TableCellRenderer, UIResource {
        private static final Border noFocusBorder = new EmptyBorder(1, 1, 1, 1);

	public MRMBooleanRenderer() {
	    super();
	    setHorizontalAlignment(JLabel.CENTER);
            setBorderPainted(true);
	}

        public Component getTableCellRendererComponent(JTable table, Object value,
						       boolean isSelected, boolean hasFocus, int row, int column) {

        if(value == null) return new JLabel("");
	    if (isSelected) {
	        setForeground(table.getSelectionForeground());
	        super.setBackground(table.getSelectionBackground());
	    }
	    else {
	        setForeground(table.getForeground());
	        setBackground(table.getBackground());
	    }
            setSelected((value != null && ((Boolean)value).booleanValue()));

            if (hasFocus) {
                setBorder(UIManager.getBorder("Table.focusCellHighlightBorder"));
            } else {
                setBorder(noFocusBorder);
            }

            return this;
        }
    }

