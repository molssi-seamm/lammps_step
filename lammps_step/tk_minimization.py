# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS minimization step"""

import lammps_step
import seamm_widgets as sw
import tkinter.ttk as ttk

from .tk_energy import stress_headers


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

        for key in (
            "convergence",
            "etol",
            "ftol",
            "nsteps",
            "nevaluations",
            "minimizer",
            "timestep",
            "optimize cell",
        ):
            self[key] = P[key].widget(opt_frame)

        # and binding to change as needed
        for key in ("optimize cell",):
            self[key].bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].bind("<Return>", self.reset_dialog)
            self[key].bind("<FocusOut>", self.reset_dialog)
        for key in ("convergence", "minimizer"):
            self[key].bind("<<ComboboxSelected>>", self.reset_optimization)
            self[key].bind("<Return>", self.reset_optimization)
            self[key].bind("<FocusOut>", self.reset_optimization)

        # Pressure/stress controls -- create frame for them
        p_frame = self["pressure frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Pressure",
            labelanchor="n",
            padding=10,
        )

        # Create the widgets
        for key in (
            "system type",
            "P",
            "allow shear",
            "use_stress",
            "couple",
            "nreset",
        ):
            self[key] = P[key].widget(p_frame)

        # The stress/pressure section is quite complicated due to
        # couplings, etc. so put in its own subsection
        s_frame = self["stress_frame"] = ttk.Frame(
            p_frame, borderwidth=4, relief="sunken"
        )

        # The last widgets
        for key in (
            "Sxx",
            "Syy",
            "Szz",
            "Sxy",
            "Sxz",
            "Syz",
        ):
            self[key] = P[key].widget(s_frame, unitswidth=6)
            self[key].units.bind("<<ComboboxSelected>>", self._handle_units)

        # The binding above and below is so that we can keep the units of stress and
        # pressure consistent throughout.
        self["P"].units.bind("<<ComboboxSelected>>", self._handle_units)

        # and labels for the directions, couplings, etc.
        for text in (
            "XX",
            "YY",
            "ZZ",
            "YZ",
            "XZ",
            "XY",
            "XX+YY",
            "XX+ZZ",
            "YY+ZZ",
            "XX+YY+ZZ",
        ):
            self[text] = ttk.Label(s_frame, text=text)

        self["stress"] = ttk.Label(s_frame, text="Stress:")

        # and adding bindings for appropriate widgets
        for key in ("system type", "use_stress"):
            self[key].bind("<<ComboboxSelected>>", self.reset_pressure_frame)
        for key in ("couple", "allow shear"):
            self[key].bind("<<ComboboxSelected>>", self.reset_stress_frame)

        # Top level needs to call reset_dialog
        if self.node.calculation == "minimization":
            self.reset_dialog()

        return frame

    def _handle_units(self, event=None):
        """Callback used to keep the pressure units consistent.

        When the units of a stress or the pressure is changed, sets all the others to
        the same units. This is needed because most of the units are hidden.

        Parameters
        ----------
        event : tkinter event object
        """

        units = event.widget.get()

        for key in (
            "P",
            "Sxx",
            "Syy",
            "Szz",
            "Sxy",
            "Sxz",
            "Syz",
        ):
            self[key].units.set(units)

    # The keys in the pressure frame, which applies only when optimizing the cell
    pressure_frame_keys = ("system type", "P", "allow shear", "use_stress", "couple")
    stress_keys = ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy")

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state.

        Which controls are shown comes from the parameters' rules
        (lammps_step.MinimizationParameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        keys = ["optimization"]
        if any(P.applies(key, values) for key in self.pressure_frame_keys):
            keys.append("pressure frame")
        keys.append("structure")
        for key in keys:
            self[key].grid(row=row, column=0, sticky="new")
            row += 1

        self.reset_optimization()
        self.reset_pressure_frame()

        return row

    def reset_optimization(self, widget=None):
        P = self.node.parameters
        values = self._widget_values()

        frame = self["optimization"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0

        # Full-width controls, with the custom convergence criteria indented
        for key, indented in (
            ("minimizer", False),
            ("convergence", False),
            ("etol", True),
            ("ftol", True),
            ("nsteps", False),
            ("nevaluations", False),
            ("timestep", False),
            ("optimize cell", False),
        ):
            if not P.applies(key, values):
                continue
            if indented:
                self[key].grid(row=row, column=1, sticky="ew")
                widgets2.append(self[key])
            else:
                self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
                widgets.append(self[key])
            row += 1

        width1 = sw.align_labels(widgets, sticky="e")
        width2 = sw.align_labels(widgets2, sticky="e")
        frame.columnconfigure(0, minsize=width1 - width2 + 30)

    def reset_pressure_frame(self, widget=None):
        """Layout the widgets for the pressure/stress control
        as needed for the current state"""
        P = self.node.parameters
        values = self._widget_values()

        # Remove all the current widgets
        p_frame = self["pressure frame"]
        for slave in p_frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        # and place the needed ones back in
        for key in ("system type", "use_stress", "allow shear", "couple"):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky="w")
                widgets.append(self[key])
                row += 1

        # A general stress, or a pressure
        if any(P.applies(key, values) for key in self.stress_keys):
            self["stress_frame"].grid(row=row, column=0, sticky="w")
            row += 1
            # and lay out the stress terms
            self.reset_stress_frame()

        for key in ("P", "nreset"):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky="w")
                widgets.append(self[key])
                row += 1

        sw.align_labels(widgets, sticky="e")

    def reset_stress_frame(self, widget=None):
        """Layout the widgets for the pressure/stress
        as needed for the current state.

        We use labels across the top and left side of the table of stresses, then
        hide the labels and units of all the entries, except for the last entry in
        the row, which displays the units too. The components shown are those that
        apply, given the coupling and whether the cell may shear.
        """
        P = self.node.parameters
        values = self._widget_values()

        keys = [key for key in self.stress_keys if P.applies(key, values)]
        headers = stress_headers(keys, values.get("couple"))
        self._grid_stress_table(
            self["stress_frame"], headers, [("stress", keys, "combobox")]
        )
