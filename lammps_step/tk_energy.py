# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):
    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
        my_logger=logger,
    ):
        """Initialize a node

        Keyword arguments:
        """
        # Metadata for the properties
        self.property_metadata = {}

        self.results_widgets = []

        # Call the constructor for the energy
        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
            my_logger=my_logger,
        )

    def right_click(self, event):
        """Probably need to add our dialog..."""

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self, title="Edit LAMMPS Energy Step"):
        """Create the dialog!"""
        frame = super().create_dialog(title=title)

        P = self.node.parameters

        # Create the structure-handling widgets
        sframe = self["structure"] = ttk.LabelFrame(
            self["frame"], text="Configuration Handling", labelanchor=tk.N
        )
        row = 0
        widgets = []
        for key in ("structure handling", "system name", "configuration name"):
            self[key] = P[key].widget(sframe)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        sw.align_labels(widgets, sticky=tk.E)

        sframe.grid(row=0, column=0, sticky=tk.N)

        self.setup_results()

        return frame

    def _widget_values(self):
        """The dialog's current values, {name: value}, for the parameters' rules."""
        return widget_values(self)

    def _grid_stress_table(self, frame, headers, rows, weight=1):
        """Grid a table of stress (or damping) widgets into a frame.

        The columns are the components that apply; the parameters' rules decide
        which they are. Only the last widget in each row shows its units.

        Parameters
        ----------
        frame : ttk.Frame
            The frame for the table.
        headers : [str]
            The keys of the column labels, one per column.
        rows : [(str, [str], str)]
            For each row, the key of its label, the keys of its widgets, and the part
            of the widgets to show ("combobox" or "entry") besides the units.
        """
        for slave in frame.grid_slaves():
            slave.grid_forget()
        for column in range(1, 7):
            frame.columnconfigure(column, weight=0, minsize=0)

        for column, key in enumerate(headers, start=1):
            self[key].grid(row=0, column=column)
        row = 1
        for label, keys, part in rows:
            self[label].grid(row=row, column=0, sticky=tk.E)
            for column, key in enumerate(keys, start=1):
                self[key].grid(row=row, column=column, sticky=tk.EW)
                if column == len(keys):
                    self[key].show(part, "units")
                else:
                    self[key].show(part)
            row += 1
        for column in range(1, len(headers) + 1):
            frame.columnconfigure(column, weight=weight, minsize=10)


def widget_values(tk_node):
    """A dialog's current values, {name: value}, for the parameters' rules.

    Parameters
    ----------
    tk_node : seamm.TkNode
        The graphical node with the dialog.

    Returns
    -------
    {str: any}
        The values of the parameters that have widgets, without units.
    """
    values = {}
    for key in tk_node.node.parameters:
        if key == "results" or key not in tk_node:
            continue
        try:
            value = tk_node[key].get()
        except Exception:
            continue
        values[key] = value[0] if isinstance(value, tuple) else value
    return values


def stress_headers(keys, couple):
    """The column labels for the stress components that apply.

    Parameters
    ----------
    keys : [str]
        The stress components shown, as their parameters' names, e.g. 'Sxx,initial'.
    couple : str
        The directions that are coupled, e.g. 'x and y'.

    Returns
    -------
    [str]
        The keys of the column labels, e.g. 'XX+YY'.
    """
    xx = {"x, y and z": "XX+YY+ZZ", "x and y": "XX+YY", "x and z": "XX+ZZ"}
    labels = {
        "xx": xx.get(couple, "XX"),
        "yy": "YY+ZZ" if couple == "y and z" else "YY",
        "zz": "ZZ",
        "xy": "XY",
        "xz": "XZ",
        "yz": "YZ",
    }
    return [labels[key[1:3].lower()] for key in keys]
