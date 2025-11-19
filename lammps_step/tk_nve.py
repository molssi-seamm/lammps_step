# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import logging
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk

logger = logging.getLogger(__name__)


class TkNVE(lammps_step.TkEnergy):
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

    def create_dialog(self, title="Edit NVE dynamics parameters"):
        """Create the dialog!"""

        # Let parent classes do their thing.
        frame = super().create_dialog(title=title)

        # Shortcut for parameters
        P = self.node.parameters

        # Add a tab for the trajectory setup
        notebook = self["notebook"]
        tframe = ttk.Frame(notebook)
        self["trajectory frame"] = tframe
        notebook.insert(
            self["results frame"], tframe, text="Trajectories", sticky=tk.NSEW
        )
        for key in (
            "trajectory",
            "trajectory save",
            "trajectory export",
            "atomic positions",
            "com positions",
            "atomic velocities",
            "com velocities",
            "heat flux",
            "shear stress",
        ):
            if title == "Heat Flux" and (key == "heat flux" or "centroid" in key):
                continue
            self[key] = P[key].widget(tframe)
            self[key].bind("<<ComboboxSelected>>", self.reset_trajectory_frame)
            self[key].bind("<Return>", self.reset_trajectory_frame)
            self[key].bind("<FocusOut>", self.reset_trajectory_frame)

        for key in (
            "trajectory rate",
            "trajectory number of samples",
            "trajectory forces",
            "trajectory velocities",
            "trajectory system name",
            "make current",
            "trajectory extxyz",
            "trajectory extxyz skip frames",
            "trajectory extxyz filename",
            "trajectory extxyz append",
            "atomic positions rate",
            "atomic positions number of samples",
            "com positions rate",
            "com positions number of samples",
            "atomic velocities rate",
            "atomic velocities number of samples",
            "com velocities rate",
            "com velocities number of samples",
            "heat flux rate",
            "heat flux number of samples",
            "shear stress rate",
            "shear stress number of samples",
        ):
            if title == "Heat Flux" and (key == "heat flux" or "centroid" in key):
                continue
            self[key] = P[key].widget(tframe)

        self.reset_trajectory_frame()

        if title == "Heat Flux":
            return frame

        # Frame to isolate widgets
        c_frame = self["control_frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="General Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in lammps_step.NVE_Parameters.parameters:
            if key == "control_properties":
                self[key] = P[key].widget(c_frame, metadata=self.property_metadata)
            else:
                self[key] = P[key].widget(c_frame)

        # make the control combobox wide enough
        self["run_control"].combobox.configure(width=40)

        # and binding to change as needed
        self["run_control"].combobox.bind(
            "<<ComboboxSelected>>", self.reset_control_frame
        )

        return frame

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self["frame"]
        # Clear the dialog
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        # Put in our control frame
        self["control_frame"].grid(row=row, column=0)
        row += 1

        # and the widgets in it
        self.reset_control_frame()

        # And how to handle the structure
        if self.node.calculation == "nve":
            self["structure"].grid(row=row, column=0)
            row += 1

        return row

    def reset_control_frame(self, widget=None):
        """Layout the control widgets as needed for the current state"""

        run_control = self["run_control"].get()

        # Clear out the previous widgets
        c_frame = self["control_frame"]
        for slave in c_frame.grid_slaves():
            slave.grid_forget()

        # And put them back in depending...
        row = 0
        # self["run_control"].grid(row=row, column=0, sticky=tk.W)
        # row += 1

        widgets = []
        if "fixed length" in run_control:
            widgets.append(self["time"])
            self["time"].grid(row=row, column=0, sticky=tk.W)
            row += 1
        else:
            widgets.append(self["maximum_time"])
            self["maximum_time"].grid(row=row, column=0, sticky=tk.W)
            row += 1
        widgets.append(self["timestep"])
        self["timestep"].grid(row=row, column=0, sticky=tk.W)
        row += 1
        if "fixed length" in run_control:
            widgets.append(self["sampling"])
            self["sampling"].grid(row=row, column=0, sticky=tk.W)
            row += 1
        else:
            self["control_properties"].grid(row=row, column=0, sticky=tk.NSEW)
            row += 1

        sw.align_labels(widgets, sticky=tk.E)

    def reset_trajectory_frame(self, widget=None):
        """Layout the trajectory frame according to its contents."""
        trajectory = self["trajectory"].get()
        save = self["trajectory save"].get()
        export = self["trajectory export"].get()
        atomic_positions = self["atomic positions"].get()
        com_positions = self["com positions"].get()
        atomic_velocities = self["atomic velocities"].get()
        com_velocities = self["com velocities"].get()
        heat_flux = self["heat flux"].get()
        shear_stress = self["shear stress"].get()

        # Clear out the previous widgets
        frame = self["trajectory frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        widgets1 = []
        widgets2 = []
        widgets3 = []
        widgets4 = []

        self["trajectory"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["trajectory"])
        row += 1
        if trajectory != "never":
            if self.is_expr(trajectory) or "number" in trajectory:
                self["trajectory number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["trajectory number of samples"])
                row += 1
            if self.is_expr(trajectory) or "interval" in trajectory:
                self["trajectory rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["trajectory rate"])
                row += 1
            for key in (
                "trajectory forces",
                "trajectory velocities",
                "trajectory save",
            ):
                self[key].grid(row=row, column=1, columnspan=4, sticky=tk.EW)
                widgets1.append(self[key])
                row += 1

            if save != "no":
                for key in ("trajectory system name", "make current"):
                    self[key].grid(row=row, column=2, columnspan=3, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1

            for key in ("trajectory export",):
                self[key].grid(row=row, column=1, columnspan=4, sticky=tk.EW)
                widgets1.append(self[key])
                row += 1

            if export != "no":
                for key in ("trajectory extxyz",):
                    self[key].grid(row=row, column=2, columnspan=3, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1
                for key in (
                    "trajectory extxyz skip frames",
                    "trajectory extxyz filename",
                    "trajectory extxyz append",
                ):
                    self[key].grid(row=row, column=3, columnspan=2, sticky=tk.EW)
                    widgets3.append(self[key])
                    row += 1

        self["atomic positions"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["atomic positions"])
        row += 1
        if atomic_positions != "never":
            if self.is_expr(atomic_positions) or "number" in atomic_positions:
                self["atomic positions number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["atomic positions number of samples"])
                row += 1
            if self.is_expr(atomic_positions) or "interval" in atomic_positions:
                self["atomic positions rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["atomic positions rate"])
                row += 1

        self["com positions"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["com positions"])
        row += 1
        if com_positions != "never":
            if self.is_expr(com_positions) or "number" in com_positions:
                self["com positions number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["com positions number of samples"])
                row += 1
            if self.is_expr(com_positions) or "interval" in com_positions:
                self["com positions rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["com positions rate"])
                row += 1

        self["atomic velocities"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["atomic velocities"])
        row += 1
        if atomic_velocities != "never":
            if self.is_expr(atomic_velocities) or "number" in atomic_velocities:
                self["atomic velocities number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["atomic velocities number of samples"])
                row += 1
            if self.is_expr(atomic_velocities) or "interval" in atomic_velocities:
                self["atomic velocities rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["atomic velocities rate"])
                row += 1

        self["com velocities"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["com velocities"])
        row += 1
        if com_velocities != "never":
            if self.is_expr(com_velocities) or "number" in com_velocities:
                self["com velocities number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["com velocities number of samples"])
                row += 1
            if self.is_expr(com_velocities) or "interval" in com_velocities:
                self["com velocities rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["com velocities rate"])
                row += 1

        if self.__class__.__name__ == "TkHeatFlux":
            pass
        else:
            self["heat flux"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
            widgets.append(self["heat flux"])
            row += 1
            if heat_flux != "never":
                if self.is_expr(heat_flux) or "number" in heat_flux:
                    self["heat flux number of samples"].grid(
                        row=row, column=1, columnspan=4, sticky=tk.EW
                    )
                    widgets1.append(self["heat flux number of samples"])
                    row += 1
                if self.is_expr(heat_flux) or "interval" in heat_flux:
                    self["heat flux rate"].grid(
                        row=row, column=1, columnspan=4, sticky=tk.EW
                    )
                    widgets1.append(self["heat flux rate"])
                    row += 1

        self["shear stress"].grid(row=row, column=0, columnspan=5, sticky=tk.EW)
        widgets.append(self["shear stress"])
        row += 1
        if shear_stress != "never":
            if self.is_expr(shear_stress) or "number" in shear_stress:
                self["shear stress number of samples"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["shear stress number of samples"])
                row += 1
            if self.is_expr(shear_stress) or "interval" in shear_stress:
                self["shear stress rate"].grid(
                    row=row, column=1, columnspan=4, sticky=tk.EW
                )
                widgets1.append(self["shear stress rate"])
                row += 1

        width = sw.align_labels(widgets, sticky=tk.E)
        width1 = sw.align_labels(widgets1, sticky=tk.E)
        indent = 75
        size = max(width - width1 + indent, 0)
        frame.columnconfigure(0, minsize=size)
        if len(widgets2) > 0:
            width2 = sw.align_labels(widgets2, sticky=tk.E)
            size = max(width1 - width2 + indent, 0)
            frame.columnconfigure(1, minsize=size)
            if len(widgets3) > 0:
                width3 = sw.align_labels(widgets3, sticky=tk.E)
                size = max(width2 - width3 + indent, 0)
                frame.columnconfigure(2, minsize=size)
                if len(widgets4) > 0:
                    width4 = sw.align_labels(widgets4, sticky=tk.E)
                    size = max(width3 - width4 + indent, 0)
                    frame.columnconfigure(2, minsize=size)
        frame.columnconfigure(4, weight=1)

    def handle_dialog(self, result):
        if result == "OK":
            # Shortcut for parameters
            P = self.node.parameters

            value, units = self["time"].get()
            P["time"].value = value
            P["time"].units = units

            tmp = self["timestep"].get()
            if tmp in P["timestep"].enumeration:
                P["timestep"].value = tmp
            else:
                P["timestep"].value = tmp[0]
                P["timestep"].units = tmp[1]

            tmp = self["sampling"].get()
            if tmp in P["sampling"].enumeration:
                P["sampling"].value = tmp
            else:
                P["sampling"].value = tmp[0]
                P["sampling"].units = tmp[1]

        # Let base classes reap their parameters
        super().handle_dialog(result)
