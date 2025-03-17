# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS minimization step"""

import lammps_step
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk


class TkMinimization(lammps_step.TkEnergy):
    def __init__(
        self, tk_flowchart=None, node=None, canvas=None, x=None, y=None, w=200, h=50
    ):
        """Initialize a node

        Keyword arguments:
        """

        super().__init__(
            tk_flowchart=tk_flowchart, node=node, canvas=canvas, x=x, y=y, w=w, h=h
        )

    def create_dialog(self):
        """Create the dialog!"""

        # Let parent classes do their thing.
        frame = super().create_dialog(title="Edit Minimization step")

        # Shortcut for parameters
        P = self.node.parameters

        # Frame to isolate widgets
        opt_frame = self["optimization"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Geometry Optimization",
            labelanchor="n",
            padding=10,
        )

        for key in lammps_step.MinimizationParameters.parameters:
            self[key] = P[key].widget(opt_frame)

        # and binding to change as needed
        for key in ("convergence", "minimizer"):
            self[key].bind("<<ComboboxSelected>>", self.reset_optimization)
            self[key].bind("<Return>", self.reset_optimization)
            self[key].bind("<FocusOut>", self.reset_optimization)

        # Top level needs to call reset_dialog
        if self.node.calculation == "minimization":
            self.reset_dialog()

        return frame

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        for key in ("optimization", "structure"):
            self[key].grid(row=row, column=0, sticky="new")
            row += 1

        self.reset_optimization()

        return row

    def reset_optimization(self, widget=None):
        convergence = self["convergence"].get()
        minimizer = self["minimizer"].get()

        frame = self["optimization"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0

        for key in ("minimizer", "convergence"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if convergence == "custom":
            for key in ("etol", "ftol"):
                self[key].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self[key])
                row += 1

        for key in ("nsteps", "nevaluations"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if minimizer in ("Fire", "QuickMin"):
            for key in ("timestep",):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self[key])
                row += 1

        width1 = sw.align_labels(widgets, sticky=tk.E)
        width2 = sw.align_labels(widgets2, sticky=tk.E)
        frame.columnconfigure(0, minsize=width1 - width2 + 30)
