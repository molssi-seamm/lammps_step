# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk


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
        for key in lammps_step.NPT_Parameters.parameters:
            if key not in self:
                if key[0] == "S":
                    self[key] = P[key].widget(s_frame, width=8, unitswidth=6)
                else:
                    self[key] = P[key].widget(s_frame)

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
        for key in ("Panneal", "use_stress", "couple"):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_stress_frame)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        row = super().reset_dialog()

        # Reset the trajectory frame to span 2 columns
        self["control_frame"].grid(row=0, column=0, columnspan=2)

        self["temperature_frame"].grid(row=row, column=0, sticky=tk.N, padx=10, pady=10)
        self.reset_temperature_frame()

        self["pressure_frame"].grid(row=row, column=1, sticky=tk.N, padx=10, pady=10)
        self.reset_pressure_frame()

        row += 1

        # And how to handle the structure
        if self.node.calculation == "npt":
            self["structure"].grid(row=row, column=0, columnspan=2)
            row += 1

        return row

    def reset_pressure_frame(self, widget=None):
        """Layout the widgets for the pressure/stress control
        as needed for the current state"""

        # Get the values that control the layout
        barostat = self["barostat"].get()
        system_type = self["system type"].get()

        # Remove all the current widgets
        p_frame = self["pressure_frame"]
        for slave in p_frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # and place the needed ones back in
        self["system type"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        self["barostat"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        self["Panneal"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        if system_type == "fluid":
            sw.align_labels(
                (self["system type"], self["barostat"], self["Panneal"]),
                sticky=tk.E,
            )
        else:
            self["use_stress"].grid(row=row, column=0, sticky=tk.W)
            row += 1

            self["couple"].grid(row=row, column=0, sticky=tk.W)
            row += 1

            sw.align_labels(
                (
                    self["system type"],
                    self["barostat"],
                    self["Panneal"],
                    self["use_stress"],
                    self["couple"],
                ),
                sticky=tk.E,
            )

        self["stress_frame"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        if barostat == "Nose-Hoover":
            self["nreset"].grid(row=row, column=0, sticky=tk.W)
            row += 1
            self["mtk"].grid(row=row, column=0, sticky=tk.W)
            sw.align_labels((self["nreset"], self["mtk"]), sticky=tk.E)
        else:
            self["modulus"].grid(row=row, column=0, sticky=tk.W)

        # and lay out the pressure or stress terms
        self.reset_stress_frame()

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
        system_type = self["system type"].get()
        barostat = self["barostat"].get()
        Panneal = self["Panneal"].get()
        if system_type == "fluid":
            use_stress = "isotropic pressure"
            couple = "x, y and z"
        else:
            use_stress = self["use_stress"].get()
            couple = self["couple"].get()

        # Remove all the current widgets
        frame = self["stress_frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        # and place the needed ones back in

        keep_orthorhombic = self["keep orthorhombic"].get()

        if use_stress != "isotropic pressure":
            # Annealing and stresses
            if couple == "x, y and z":
                # all stresses
                self["XX+YY+ZZ"].grid(row=row, column=1)
                if not keep_orthorhombic:
                    self["XY"].grid(row=row, column=2)
                    self["XZ"].grid(row=row, column=3)
                    self["YZ"].grid(row=row, column=4)

                row += 1

                self["initial stress"].grid(row=row, column=0, sticky=tk.E)
                if keep_orthorhombic:
                    self["Sxx,initial"].grid(
                        row=row, column=1, columnspan=3, sticky=tk.W
                    )
                    self["Sxx,initial"].show("entry", "units")
                    frame.columnconfigure(1, weight=1, uniform="b")
                    frame.columnconfigure(2, weight=1, uniform="b")
                else:
                    self["Sxx,initial"].grid(row=row, column=1)
                    self["Sxx,initial"].show("entry")

                    self["Sxy,initial"].grid(row=row, column=2)
                    self["Sxz,initial"].grid(row=row, column=3)
                    self["Syz,initial"].grid(
                        row=row, column=4, columnspan=3, sticky=tk.W
                    )

                    self["Sxy,initial"].show("entry")
                    self["Sxz,initial"].show("entry")
                    self["Syz,initial"].show("entry", "units")

                    frame.columnconfigure(1, weight=1, uniform="b")
                    frame.columnconfigure(2, weight=1, uniform="b")
                    frame.columnconfigure(3, weight=1, uniform="b")
                    frame.columnconfigure(4, weight=1, uniform="b")
                    frame.columnconfigure(5, weight=1, uniform="b")
                    frame.columnconfigure(6, weight=1, uniform="b")
                row += 1

                if Panneal != "no":
                    self["final stress"].grid(row=row, column=0, sticky=tk.E)
                    if keep_orthorhombic:
                        self["Sxx,final"].grid(
                            row=row, column=1, columnspan=3, sticky=tk.W
                        )
                        self["Sxx,final"].show("entry", "units")
                    else:
                        self["Sxx,final"].grid(row=row, column=1)
                        self["Sxx,final"].show("entry")

                        self["Sxy,final"].grid(row=row, column=2)
                        self["Sxz,final"].grid(row=row, column=3)
                        self["Syz,final"].grid(
                            row=row, column=4, columnspan=3, sticky=tk.W
                        )

                        self["Sxy,final"].show("entry")
                        self["Sxz,final"].show("entry")
                        self["Syz,final"].show("entry", "units")

                    row += 1

                if barostat == "Nose-Hoover":
                    self["damping"].grid(row=row, column=0, sticky=tk.E)
                    if keep_orthorhombic:
                        self["Sxx damp"].grid(
                            row=row, column=1, columnspan=3, sticky=tk.W
                        )
                        self["Sxx damp"].show("entry", "units")
                    else:
                        self["Sxx damp"].grid(row=row, column=1)
                        self["Sxx damp"].show("entry")

                        self["Sxy damp"].grid(row=row, column=2)
                        self["Sxz damp"].grid(row=row, column=3)
                        self["Syz damp"].grid(
                            row=row, column=4, columnspan=3, sticky=tk.W
                        )

                        self["Sxy damp"].show("entry")
                        self["Sxz damp"].show("entry")
                        self["Syz damp"].show("entry", "units")
            elif couple == "x and y":
                # couple xx and yy
                self["XX+YY"].grid(row=row, column=1)
                self["ZZ"].grid(row=row, column=2)
                if not keep_orthorhombic:
                    self["XY"].grid(row=row, column=3)
                    self["XZ"].grid(row=row, column=4)
                    self["YZ"].grid(row=row, column=5)

                row += 1

                if Panneal != "no":
                    self["initial stress"].grid(row=row, column=0, sticky=tk.E)
                else:
                    self["stress"].grid(row=row, column=0, sticky=tk.E)
                self["Sxx,initial"].grid(row=row, column=1)

                self["Sxx,initial"].show("entry")
                if keep_orthorhombic:
                    self["Szz,initial"].grid(
                        row=row, column=2, columnspan=3, sticky=tk.W
                    )
                    self["Szz,initial"].show("entry", "units")
                    frame.columnconfigure(1, weight=1, uniform="c")
                    frame.columnconfigure(2, weight=1, uniform="c")
                    frame.columnconfigure(3, weight=1, uniform="c")
                else:
                    self["Szz,initial"].show("entry")

                    self["Sxy,initial"].grid(row=row, column=3)
                    self["Sxz,initial"].grid(row=row, column=4)
                    self["Syz,initial"].grid(
                        row=row, column=5, columnspan=3, sticky=tk.W
                    )

                    self["Sxy,initial"].show("entry")
                    self["Sxz,initial"].show("entry")
                    self["Syz,initial"].show("entry", "units")

                    frame.columnconfigure(1, weight=1, uniform="c")
                    frame.columnconfigure(2, weight=1, uniform="c")
                    frame.columnconfigure(3, weight=1, uniform="c")
                    frame.columnconfigure(4, weight=1, uniform="c")
                    frame.columnconfigure(5, weight=1, uniform="c")
                    frame.columnconfigure(6, weight=1, uniform="c")
                row += 1

                if Panneal != "no":
                    self["final stress"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx,final"].grid(row=row, column=1)

                    self["Sxx,final"].show("entry")
                    if keep_orthorhombic:
                        self["Szz,final"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Szz,final"].show("entry", "units")
                    else:
                        self["Szz,final"].grid(row=row, column=2)
                        self["Szz,final"].show("entry")

                        self["Sxy,final"].grid(row=row, column=3)
                        self["Sxz,final"].grid(row=row, column=4)
                        self["Syz,final"].grid(
                            row=row, column=5, columnspan=3, sticky=tk.W
                        )

                        self["Sxy,final"].show("entry")
                        self["Sxz,final"].show("entry")
                        self["Syz,final"].show("entry", "units")

                    row += 1

                if barostat == "Nose-Hoover":
                    self["damping"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx damp"].grid(row=row, column=1)

                    self["Sxx damp"].show("entry")
                    if keep_orthorhombic:
                        self["Szz damp"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Szz damp"].show("entry", "units")
                    else:
                        self["Szz damp"].grid(row=row, column=2)
                        self["Szz damp"].show("entry")

                        self["Sxy damp"].grid(row=row, column=3)
                        self["Sxz damp"].grid(row=row, column=4)
                        self["Syz damp"].grid(
                            row=row, column=5, columnspan=3, sticky=tk.W
                        )

                        self["Sxy damp"].show("entry")
                        self["Sxz damp"].show("entry")
                        self["Syz damp"].show("entry", "units")
            elif couple == "x and z":
                # couple xx and zz
                self["XX+ZZ"].grid(row=row, column=1)
                self["YY"].grid(row=row, column=2)
                if not keep_orthorhombic:
                    self["XY"].grid(row=row, column=3)
                    self["XZ"].grid(row=row, column=4)
                    self["YZ"].grid(row=row, column=5)

                row += 1

                if Panneal != "no":
                    self["initial stress"].grid(row=row, column=0, sticky=tk.E)
                else:
                    self["stress"].grid(row=row, column=0, sticky=tk.E)
                self["Sxx,initial"].grid(row=row, column=1)

                self["Sxx,initial"].show("entry")
                if keep_orthorhombic:
                    self["Syy,initial"].grid(
                        row=row, column=2, columnspan=3, sticky=tk.W
                    )
                    self["Syy,initial"].show("entry", "units")
                    frame.columnconfigure(1, weight=1, uniform="d")
                    frame.columnconfigure(2, weight=1, uniform="d")
                    frame.columnconfigure(3, weight=1, uniform="d")
                else:
                    self["Syy,initial"].grid(row=row, column=2)
                    self["Syy,initial"].show("entry")

                    self["Sxy,initial"].grid(row=row, column=3)
                    self["Sxz,initial"].grid(row=row, column=4)
                    self["Syz,initial"].grid(
                        row=row, column=5, columnspan=3, sticky=tk.W
                    )

                    self["Sxy,initial"].show("entry")
                    self["Sxz,initial"].show("entry")
                    self["Syz,initial"].show("entry", "units")

                    frame.columnconfigure(1, weight=1, uniform="d")
                    frame.columnconfigure(2, weight=1, uniform="d")
                    frame.columnconfigure(3, weight=1, uniform="d")
                    frame.columnconfigure(4, weight=1, uniform="d")
                    frame.columnconfigure(5, weight=1, uniform="d")
                    frame.columnconfigure(6, weight=1, uniform="d")
                row += 1

                if Panneal != "no":
                    self["final stress"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx,final"].grid(row=row, column=1)

                    self["Sxx,final"].show("entry")
                    if keep_orthorhombic:
                        self["Syy,final"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Syy,final"].show("entry", "units")
                    else:
                        self["Syy,final"].grid(row=row, column=2)
                        self["Syy,final"].show("entry")

                        self["Sxy,final"].grid(row=row, column=3)
                        self["Sxz,final"].grid(row=row, column=4)
                        self["Syz,final"].grid(
                            row=row, column=5, columnspan=3, sticky=tk.W
                        )

                        self["Sxy,final"].show("entry")
                        self["Sxz,final"].show("entry")
                        self["Syz,final"].show("entry", "units")
                    row += 1

                if barostat == "Nose-Hoover":
                    self["damping"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx damp"].grid(row=row, column=1)

                    self["Sxx damp"].show("entry")
                    if keep_orthorhombic:
                        self["Syy damp"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Syy damp"].show("entry", "units")
                    else:
                        self["Syy damp"].grid(row=row, column=2)
                        self["Syy damp"].show("entry")

                        self["Sxy damp"].grid(row=row, column=3)
                        self["Sxz damp"].grid(row=row, column=4)
                        self["Syz damp"].grid(row=row, column=5)

                        self["Sxy damp"].show("entry")
                        self["Sxz damp"].show("entry")
                        self["Syz damp"].show("entry", "units")
            elif couple == "y and z":
                # couple yy and zz
                self["XX"].grid(row=row, column=1)
                self["YY+ZZ"].grid(row=row, column=2)
                if not keep_orthorhombic:
                    self["XY"].grid(row=row, column=3)
                    self["XZ"].grid(row=row, column=4)
                    self["YZ"].grid(row=row, column=5)

                row += 1

                if Panneal != "no":
                    self["initial stress"].grid(row=row, column=0, sticky=tk.E)
                else:
                    self["stress"].grid(row=row, column=0, sticky=tk.E)
                self["Sxx,initial"].grid(row=row, column=1)

                self["Sxx,initial"].show("entry")
                if keep_orthorhombic:
                    self["Syy,initial"].grid(
                        row=row, column=2, columnspan=3, sticky=tk.W
                    )
                    self["Syy,initial"].show("entry", "units")
                    frame.columnconfigure(1, weight=1, uniform="e")
                    frame.columnconfigure(2, weight=1, uniform="e")
                    frame.columnconfigure(3, weight=1, uniform="e")
                else:
                    self["Syy,initial"].grid(row=row, column=2)
                    self["Syy,initial"].show("entry")

                    self["Sxy,initial"].grid(row=row, column=3)
                    self["Sxz,initial"].grid(row=row, column=4)
                    self["Syz,initial"].grid(row=row, column=5)

                    self["Sxy,initial"].show("entry")
                    self["Sxz,initial"].show("entry")
                    self["Syz,initial"].show("entry", "units")

                    frame.columnconfigure(1, weight=1, uniform="e")
                    frame.columnconfigure(2, weight=1, uniform="e")
                    frame.columnconfigure(3, weight=1, uniform="e")
                    frame.columnconfigure(4, weight=1, uniform="e")
                    frame.columnconfigure(5, weight=1, uniform="e")
                    frame.columnconfigure(6, weight=1, uniform="e")
                row += 1

                if Panneal != "no":
                    self["final stress"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx,final"].grid(row=row, column=1)

                    self["Sxx,final"].show("entry")
                    if keep_orthorhombic:
                        self["Syy,final"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Syy,final"].show("entry", "units")
                    else:
                        self["Syy,final"].grid(row=row, column=2)
                        self["Syy,final"].show("entry")

                        self["Sxy,final"].grid(row=row, column=3)
                        self["Sxz,final"].grid(row=row, column=4)
                        self["Syz,final"].grid(
                            row=row, column=5, columnspan=3, sticky=tk.W
                        )

                        self["Sxy,final"].show("entry")
                        self["Sxz,final"].show("entry")
                        self["Syz,final"].show("entry", "units")

                    row += 1

                if barostat == "Nose-Hoover":
                    self["damping"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx damp"].grid(row=row, column=1)

                    self["Sxx damp"].show("entry")
                    if keep_orthorhombic:
                        self["Syy damp"].grid(
                            row=row, column=2, columnspan=3, sticky=tk.W
                        )
                        self["Syy damp"].show("entry", "units")
                    else:
                        self["Syy damp"].grid(row=row, column=2)
                        self["Syy damp"].show("entry")

                        self["Sxy damp"].grid(row=row, column=3)
                        self["Sxz damp"].grid(row=row, column=4)
                        self["Syz damp"].grid(
                            row=row, column=5, columnspan=3, sticky=tk.W
                        )

                        self["Sxy damp"].show("entry")
                        self["Sxz damp"].show("entry")
                        self["Syz damp"].show("entry", "units")
            else:
                # if couple == 'none':
                # all stresses
                self["XX"].grid(row=row, column=1)
                self["YY"].grid(row=row, column=2)
                self["ZZ"].grid(row=row, column=3)
                if not keep_orthorhombic:
                    self["XY"].grid(row=row, column=4)
                    self["XZ"].grid(row=row, column=5)
                    self["YZ"].grid(row=row, column=6)

                row += 1

                if Panneal != "no":
                    self["initial stress"].grid(row=row, column=0, sticky=tk.E)
                else:
                    self["stress"].grid(row=row, column=0, sticky=tk.E)
                self["Sxx,initial"].grid(row=row, column=1)
                self["Syy,initial"].grid(row=row, column=2)

                self["Sxx,initial"].show("entry")
                self["Syy,initial"].show("entry")
                if keep_orthorhombic:
                    self["Szz,initial"].grid(
                        row=row, column=3, columnspan=3, sticky=tk.W
                    )
                    self["Szz,initial"].show("entry", "units")
                    frame.columnconfigure(1, weight=1, uniform="a")
                    frame.columnconfigure(2, weight=1, uniform="a")
                    frame.columnconfigure(3, weight=1, uniform="a")
                    frame.columnconfigure(4, weight=1, uniform="a")
                else:
                    self["Szz,initial"].grid(row=row, column=3)
                    self["Szz,initial"].show("entry")

                    self["Sxy,initial"].grid(row=row, column=4)
                    self["Sxz,initial"].grid(row=row, column=5)
                    self["Syz,initial"].grid(
                        row=row, column=6, columnspan=3, sticky=tk.W
                    )

                    self["Sxy,initial"].show("entry")
                    self["Sxz,initial"].show("entry")
                    self["Syz,initial"].show("entry", "units")

                    frame.columnconfigure(1, weight=1, uniform="a")
                    frame.columnconfigure(2, weight=1, uniform="a")
                    frame.columnconfigure(3, weight=1, uniform="a")
                    frame.columnconfigure(4, weight=1, uniform="a")
                    frame.columnconfigure(5, weight=1, uniform="a")
                    frame.columnconfigure(6, weight=1, uniform="a")
                    frame.columnconfigure(7, weight=1, uniform="a")
                row += 1

                if Panneal != "no":
                    self["final stress"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx,final"].grid(row=row, column=1)
                    self["Syy,final"].grid(row=row, column=2)

                    self["Sxx,final"].show("entry")
                    self["Syy,final"].show("entry")

                    if keep_orthorhombic:
                        self["Szz,final"].grid(
                            row=row, column=3, columnspan=3, sticky=tk.W
                        )
                        self["Szz,final"].show("entry", "units")
                    else:
                        self["Szz,final"].grid(row=row, column=3)
                        self["Szz,final"].show("entry")

                        self["Sxy,final"].grid(row=row, column=4)
                        self["Sxz,final"].grid(row=row, column=5)
                        self["Syz,final"].grid(
                            row=row, column=6, columnspan=3, sticky=tk.W
                        )

                        self["Sxy,final"].show("entry")
                        self["Sxz,final"].show("entry")
                        self["Syz,final"].show("entry", "units")
                    row += 1

                if barostat == "Nose-Hoover":
                    self["damping"].grid(row=row, column=0, sticky=tk.E)
                    self["Sxx damp"].grid(row=row, column=1)
                    self["Syy damp"].grid(row=row, column=2)
                    self["Szz damp"].grid(row=row, column=3)

                    self["Sxx damp"].show("entry")
                    self["Syy damp"].show("entry")
                    if keep_orthorhombic:
                        self["Szz damp"].grid(
                            row=row, column=3, columnspan=3, sticky=tk.W
                        )
                        self["Szz damp"].show("entry", "units")
                    else:
                        self["Szz damp"].grid(row=row, column=3)
                        self["Szz damp"].show("entry")

                        self["Sxy damp"].grid(row=row, column=4)
                        self["Sxz damp"].grid(row=row, column=5)
                        self["Syz damp"].grid(
                            row=row, column=6, columnspan=3, sticky=tk.W
                        )

                        self["Sxy damp"].show("entry")
                        self["Sxz damp"].show("entry")
                        self["Syz damp"].show("entry", "units")
        else:
            widgets = []
            if Panneal != "no":
                self["Pinitial"].label.configure(text="Initial pressure:")
            else:
                self["Pinitial"].label.configure(text="Pressure:")
            self["Pinitial"].grid(row=row, column=0, sticky=tk.W)
            widgets.append(self["Pinitial"])
            row += 1

            if Panneal != "no":
                self["Pfinal"].grid(row=row, column=0, sticky=tk.W)
                widgets.append(self["Pfinal"])
                row += 1

            if barostat == "Nose-Hoover":
                self["Pdamp"].grid(row=row, column=0, sticky=tk.W)
                widgets.append(self["Pdamp"])
                row += 1
            sw.align_labels(widgets, sticky=tk.E)

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
