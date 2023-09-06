# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS velocities step"""

import seamm
import seamm_widgets as sw
import tkinter as tk


class TkVelocities(seamm.TkNode):
    def __init__(
        self, tk_flowchart=None, node=None, canvas=None, x=None, y=None, w=200, h=50
    ):
        """Initialize a node

        Keyword arguments:
        """

        super().__init__(
            tk_flowchart=tk_flowchart, node=node, canvas=canvas, x=x, y=y, w=w, h=h
        )

    def right_click(self, event):
        """Probably need to add our dialog..."""

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self):
        """Create the dialog!"""

        frame = super().create_dialog("Edit LAMMPS Velocity Step")

        # Create all the widgets
        P = self.node.parameters
        for key in P:
            self[key] = P[key].widget(frame)

        self["method"].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
        self["method"].combobox.bind("<Return>", self.reset_dialog)
        self["method"].combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the widgets given the current state"""
        method = self["method"].get()

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        for key in ("method", "T", "remove_momentum"):
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if "scaling" not in method:
            self["seed"].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self["seed"])
            row += 1

            sw.align_labels(widgets, sticky=tk.E)
