# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import seamm_widgets as sw
import tkinter.ttk as ttk

from .tk_energy import stress_headers


class TkNPT(lammps_step.TkNVT):
    def __init__(
        self, tk_flowchart=None, node=None, canvas=None, x=None, y=None, w=200, h=50
    ):
        """Initialize a node

        Keyword arguments:
        """

        super().__init__(
            tk_flowchart=tk_flowchart, node=node, canvas=canvas, x=x, y=y, w=w, h=h
        )

    def create_dialog(self, title="Edit NPT dynamics parameters"):
        """Create the edit dialog!

        This is reasonably complicated, so a bit of description
        is in order. The superclasses NVT and NVE create the dialog
        along with the basic runtime and timestep (from NVE) and
        the setting the temperature and thermostat (NVT).

        This method adds a third frame with the pressure and barostat.

        The layout is handled in part by the NVT superclass, which
        handles the temperature frame. Our part is handled by three
        methods:

        * reset_dialog does the general layout of the large blocks
        * reset_pressure_frame handles the layout of the pressure
          section except for the detail of the actual pressure/stress
          which is handled by...
        * reset_stress_frame which does the detailed layout of their
          stress or pressure terms, depending on whether we are annealing,
          and the coupling between directions.
        """

        # Let parent classes do their thing.
        super().create_dialog(title=title)

        # Shortcut for parameters
        P = self.node.parameters

        # Pressure/stress controls -- create frame for them
        p_frame = self["pressure_frame"] = ttk.LabelFrame(
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
            "barostat",
            "Panneal",
            "allow shear",
            "use_stress",
            "couple",
            "nreset",
            "mtk",
            "modulus",
        ):
            self[key] = P[key].widget(p_frame)

        # The stress/pressure section is quite complicated due to
        # couplings, etc. so put in its own subsection
        s_frame = self["stress_frame"] = ttk.Frame(
            p_frame, borderwidth=4, relief="sunken"
        )

        # The last widgets
        for key in ("Pinitial", "Pfinal", "Pdamp"):
            self[key] = P[key].widget(s_frame)
        for key in (
            "Sxx,initial",
            "Syy,initial",
            "Szz,initial",
            "Sxy,initial",
            "Sxz,initial",
            "Syz,initial",
            "Sxx,final",
            "Syy,final",
            "Szz,final",
            "Sxy,final",
            "Sxz,final",
            "Syz,final",
        ):
            self[key] = P[key].widget(s_frame, unitswidth=6, width=8)

        for key in (
            "Sxx damp",
            "Syy damp",
            "Szz damp",
            "Sxy damp",
            "Sxz damp",
            "Syz damp",
        ):
            self[key] = P[key].widget(s_frame, unitswidth=6, width=8)

        # and labels for the directions, couplings, etc.
        for text in (
            "XX",
            "YY",
            "ZZ",
            "XY",
            "XZ",
            "YZ",
            "XX+YY",
            "XX+ZZ",
            "YY+ZZ",
            "XX+YY+ZZ",
        ):
            self[text] = ttk.Label(s_frame, text=text)

        self["stress"] = ttk.Label(s_frame, text="Stress:")
        self["initial stress"] = ttk.Label(s_frame, text="Initial stress:")
        self["final stress"] = ttk.Label(s_frame, text="Final stress:")
        self["damping"] = ttk.Label(s_frame, text="Damping:")

        # and adding bindings for appropriate widgets
        for key in ("system type", "barostat"):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_pressure_frame)
        for key in ("Panneal", "use_stress", "allow shear", "couple"):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_stress_frame)
        # Keep the units of stress and pressure consistent throughout.
        for key in (
            "Pinitial",
            "Pfinal",
            "Sxx,initial",
            "Syy,initial",
            "Szz,initial",
            "Sxy,initial",
            "Sxz,initial",
            "Syz,initial",
            "Sxx,final",
            "Syy,final",
            "Szz,final",
            "Sxy,final",
            "Sxz,final",
            "Syz,final",
        ):
            self[key].units.bind("<<ComboboxSelected>>", self._handle_units)

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
            "Pinitial",
            "Pfinal",
            "Sxx,initial",
            "Syy,initial",
            "Szz,initial",
            "Sxy,initial",
            "Sxz,initial",
            "Syz,initial",
            "Sxx,final",
            "Syy,final",
            "Szz,final",
            "Sxy,final",
            "Sxz,final",
            "Syz,final",
        ):
            self[key].units.set(units)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        row = super().reset_dialog()

        self["temperature_frame"].grid(row=0, column=1, sticky="n", padx=10, pady=10)
        self.reset_temperature_frame()

        self["pressure_frame"].grid(
            row=row,
            column=0,
            columnspan=2,
            sticky="n",
            padx=10,
            pady=10,
        )
        self.reset_pressure_frame()

        row += 1

        # And how to handle the structure
        if self.node.calculation == "npt":
            self["structure"].grid(row=row, column=0, columnspan=2)
            row += 1

        frame = self["frame"]
        frame.columnconfigure(1, weight=0, uniform="a")
        frame.columnconfigure(1, weight=1, uniform="a")

        return row

    # The stress components, in the order of the columns of the table
    stress_components = ("Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz")

    def reset_pressure_frame(self, widget=None):
        """Layout the widgets for the pressure/stress control
        as needed for the current state.

        Which controls are shown comes from the parameters' rules
        (lammps_step.NPT_Parameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()

        # Remove all the current widgets
        p_frame = self["pressure_frame"]
        for slave in p_frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []

        # and place the needed ones back in
        for key in (
            "system type",
            "barostat",
            "Panneal",
            "allow shear",
            "use_stress",
            "couple",
        ):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky="ew")
                widgets.append(self[key])
                row += 1

        sw.align_labels(widgets, sticky="e")

        self["stress_frame"].grid(row=row, column=0, sticky="ew")
        row += 1

        # The barostat's own controls
        widgets = []
        for key in ("nreset", "mtk", "modulus"):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky="ew")
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky="e")

        # and lay out the pressure or stress terms
        self.reset_stress_frame()

    def reset_stress_frame(self, widget=None):
        """Layout the widgets for the pressure/stress
        as needed for the current state.

        A table of the stresses that apply -- given the coupling and whether the
        cell may shear -- with rows for the initial and, if annealing, final stress,
        and the damping times for the Nose-Hoover barostat; or the pressure. We use
        labels across the top and left side of the table of stresses, then hide the
        labels and units of all the entries, except for the last entry in the row,
        which displays the units too.
        """
        P = self.node.parameters
        values = self._widget_values()
        anneal = values.get("Panneal") != "no"

        frame = self["stress_frame"]

        columns = [
            key for key in self.stress_components if P.applies(f"{key},initial", values)
        ]
        rows = []
        if len(columns) > 0:
            rows.append(
                (
                    "initial stress" if anneal else "stress",
                    [f"{key},initial" for key in columns],
                    "combobox",
                )
            )
        for label, suffix, part in (
            ("final stress", ",final", "combobox"),
            ("damping", " damp", "entry"),
        ):
            keys = [key + suffix for key in columns if P.applies(key + suffix, values)]
            if len(keys) > 0:
                rows.append((label, keys, part))
        headers = stress_headers(columns, values.get("couple"))
        self._grid_stress_table(frame, headers, rows)

        # The (isotropic) pressure, below any table
        row = len(rows) + 1 if len(rows) > 0 else 0
        widgets = []
        if anneal:
            self["Pinitial"].label.configure(text="Initial pressure:")
        else:
            self["Pinitial"].label.configure(text="Pressure:")
        for key in ("Pinitial", "Pfinal", "Pdamp"):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky="ew")
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky="e")

    def handle_dialog(self, result):
        """Handle when the user clicks a button on the dialog,
        which can either be 'Cancel', 'Help' or 'OK', or they can
        close the dialog with the 'x' button == 'Cancel'"""

        if result == "OK":
            # Shortcut for parameters
            P = self.node.parameters

            # Gather all the parameters. It maybe a bit of an overkill,
            # but it is easy.
            for key in lammps_step.NPT_Parameters.parameters:
                if key not in ("results", "create table"):
                    P[key].set_from_widget()

        # Let base classes reap their parameters
        super().handle_dialog(result)
