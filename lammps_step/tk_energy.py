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
