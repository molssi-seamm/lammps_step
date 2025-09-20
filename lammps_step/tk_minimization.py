# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS minimization step"""

import lammps_step
import seamm_widgets as sw
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

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        optimize_cell = self["optimize cell"].get() != "no"

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        if optimize_cell:
            keys = ("optimization", "pressure frame", "structure")
        else:
            keys = ("optimization", "structure")
        for key in keys:
            self[key].grid(row=row, column=0, sticky="new")
            row += 1

        self.reset_optimization()
        self.reset_pressure_frame()

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
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        if convergence == "custom":
            for key in ("etol", "ftol"):
                self[key].grid(row=row, column=1, sticky="ew")
                widgets2.append(self[key])
                row += 1

        for key in ("nsteps", "nevaluations"):
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        if minimizer in ("Fire", "QuickMin"):
            for key in ("timestep",):
                self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
                widgets.append(self[key])
                row += 1

        for key in ("optimize cell",):
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        width1 = sw.align_labels(widgets, sticky="e")
        width2 = sw.align_labels(widgets2, sticky="e")
        frame.columnconfigure(0, minsize=width1 - width2 + 30)

    def reset_pressure_frame(self, widget=None):
        """Layout the widgets for the pressure/stress control
        as needed for the current state"""

        # Get the values that control the layout
        system_type = self["system type"].get()
        use_stress = "pressure" not in self["use_stress"].get()

        # Remove all the current widgets
        p_frame = self["pressure frame"]
        for slave in p_frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        # and place the needed ones back in
        self["system type"].grid(row=row, column=0, sticky="w")
        widgets.append(self["system type"])
        row += 1

        if system_type == "fluid":
            self["P"].grid(row=row, column=0, sticky="w")
            widgets.append(self["P"])
            row += 1
        else:
            for key in ("use_stress", "allow shear", "couple"):
                self[key].grid(row=row, column=0, sticky="w")
                widgets.append(self[key])
                row += 1

            if use_stress:
                self["stress_frame"].grid(row=row, column=0, sticky="w")
                row += 1
                # and lay out the stress terms
                self.reset_stress_frame()
            else:
                self["P"].grid(row=row, column=0, sticky="w")
                widgets.append(self["P"])
                row += 1

        self["nreset"].grid(row=row, column=0, sticky="w")
        widgets.append(self["nreset"])
        row += 1

        sw.align_labels(widgets, sticky="e")

    def reset_stress_frame(self, widget=None):
        """Layout the widgets for the pressure/stress
        as needed for the current state.

        We use labels across the top and left side of the
        table of stresses, then hide the labels and units
        of all the entries, except for the last entry in the
        row, which displays the units too.

        It is a bit tricky getting the column labels to lign
        up for the last item in the row because it is wider
        with the units displaying. The weird stuff with the
        columnspan=3 and setting the 'uniform' attribute of
        the columns does this. There might be better ways, but
        this works ... just be careful changing it :-).
        """

        # Get the values that control the layout
        couple = self["couple"].get()
        allow_shear = self["allow shear"].get() != "no"

        # Remove all the current widgets
        frame = self["stress_frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        frame.columnconfigure(1, weight=0)
        frame.columnconfigure(2, weight=0)
        frame.columnconfigure(3, weight=0)
        frame.columnconfigure(4, weight=0)
        frame.columnconfigure(5, weight=0)
        frame.columnconfigure(6, weight=0)

        row = 0
        # and place the needed ones back in
        if couple == "x, y and z":
            # all stresses
            self["XX+YY+ZZ"].grid(row=row, column=1)
            if allow_shear:
                self["YZ"].grid(row=row, column=2)
                self["XZ"].grid(row=row, column=3)
                self["XY"].grid(row=row, column=4)
            row += 1
            self["stress"].grid(row=row, column=0, sticky="e")
            self["Sxx"].grid(row=row, column=1, sticky="ew")
            if not allow_shear:
                self["Sxx"].show("combobox", "units")
                frame.columnconfigure(1, weight=1, uniform="b")
            else:
                self["Sxx"].show("combobox")
                self["Syz"].grid(row=row, column=2, sticky="ew")
                self["Sxz"].grid(row=row, column=3, sticky="ew")
                self["Sxy"].grid(row=row, column=4, sticky="ew")
                self["Syz"].show("combobox")
                self["Sxz"].show("combobox")
                self["Sxy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
                frame.columnconfigure(4, weight=1)
            row += 1
        elif couple == "x and y":
            # couple xx and yy
            self["XX+YY"].grid(row=row, column=1)
            self["ZZ"].grid(row=row, column=2)
            if allow_shear:
                self["YZ"].grid(row=row, column=3)
                self["XZ"].grid(row=row, column=4)
                self["XY"].grid(row=row, column=5)
            row += 1
            self["stress"].grid(row=row, column=0, sticky="e")
            self["Sxx"].grid(row=row, column=1, sticky="ew")
            self["Sxx"].show("combobox")
            self["Szz"].grid(row=row, column=2, sticky="ew")
            if not allow_shear:
                self["Szz"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
            else:
                self["Szz"].show("combobox")
                self["Syz"].grid(row=row, column=3, sticky="ew")
                self["Sxz"].grid(row=row, column=4, sticky="ew")
                self["Sxy"].grid(row=row, column=5, sticky="ew")
                self["Syz"].show("combobox")
                self["Sxz"].show("combobox")
                self["Sxy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
                frame.columnconfigure(4, weight=1)
                frame.columnconfigure(5, weight=1)
            row += 1
        elif couple == "x and z":
            # couple xx and zz
            self["XX+ZZ"].grid(row=row, column=1)
            self["YY"].grid(row=row, column=2)
            if allow_shear:
                self["YZ"].grid(row=row, column=3)
                self["XZ"].grid(row=row, column=4)
                self["XY"].grid(row=row, column=5)
            row += 1
            self["stress"].grid(row=row, column=0, sticky="e")
            self["Sxx"].grid(row=row, column=1, sticky="ew")
            self["Syy"].grid(row=row, column=2, sticky="ew")
            self["Sxx"].show("combobox")
            if not allow_shear:
                self["Syy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
            else:
                self["Syy"].show("combobox")
                self["Syz"].grid(row=row, column=3, sticky="ew")
                self["Sxz"].grid(row=row, column=4, sticky="ew")
                self["Sxy"].grid(row=row, column=5, sticky="ew")
                self["Syz"].show("combobox")
                self["Sxz"].show("combobox", "units")
                self["Sxy"].show("combobox")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
                frame.columnconfigure(4, weight=1)
                frame.columnconfigure(5, weight=1)
            row += 1
        elif couple == "y and z":
            # couple yy and zz
            self["XX"].grid(row=row, column=1)
            self["YY+ZZ"].grid(row=row, column=2)
            if allow_shear:
                self["YZ"].grid(row=row, column=3)
                self["XZ"].grid(row=row, column=4)
                self["XY"].grid(row=row, column=5)
            row += 1
            self["stress"].grid(row=row, column=0, sticky="e")
            self["Sxx"].grid(row=row, column=1, sticky="ew")
            self["Syy"].grid(row=row, column=2, sticky="ew")
            self["Sxx"].show("combobox")
            if not allow_shear:
                self["Syy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
            else:
                self["Syy"].show("combobox")
                self["Syz"].grid(row=row, column=3, sticky="ew")
                self["Sxz"].grid(row=row, column=4, sticky="ew")
                self["Sxy"].grid(row=row, column=5, sticky="ew")
                self["Syz"].show("combobox")
                self["Sxz"].show("combobox")
                self["Sxy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
                frame.columnconfigure(4, weight=1)
                frame.columnconfigure(5, weight=1)
            row += 1
        else:
            # if couple == 'none':
            # all stresses
            self["XX"].grid(row=row, column=1)
            self["YY"].grid(row=row, column=2)
            self["ZZ"].grid(row=row, column=3)
            if allow_shear:
                self["YZ"].grid(row=row, column=4)
                self["XZ"].grid(row=row, column=5)
                self["XY"].grid(row=row, column=6)
            row += 1
            self["stress"].grid(row=row, column=0, sticky="e")
            self["Sxx"].grid(row=row, column=1, sticky="ew")
            self["Sxx"].show("combobox")
            self["Syy"].grid(row=row, column=2, sticky="ew")
            self["Syy"].show("combobox")
            self["Szz"].grid(row=row, column=3, sticky="ew")
            if not allow_shear:
                self["Szz"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
            else:
                self["Szz"].show("combobox")
                self["Syz"].grid(row=row, column=4, sticky="ew")
                self["Sxz"].grid(row=row, column=5, sticky="ew")
                self["Sxy"].grid(row=row, column=6, sticky="ew")
                self["Syz"].show("combobox")
                self["Sxz"].show("combobox")
                self["Sxy"].show("combobox", "units")
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, weight=1)
                frame.columnconfigure(3, weight=1)
                frame.columnconfigure(4, weight=1)
                frame.columnconfigure(5, weight=1)
                frame.columnconfigure(6, weight=1)
            row += 1
