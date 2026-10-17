# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import logging
import Pmw
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk

from lammps_step.extras import DISPLAY_UNITS, PREFACTORS, TREATMENTS, check

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
        # The rows of the extra properties table: a dict of widgets per row
        self._extra_rows = []

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
            "pressure tensor",
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
            "pressure tensor rate",
            "pressure tensor number of samples",
            "keep sampled files",
        ):
            if title == "Heat Flux" and (key == "heat flux" or "centroid" in key):
                continue
            self[key] = P[key].widget(tframe)

        self.reset_trajectory_frame()

        if title == "Heat Flux":
            return frame

        self._create_extras_tab(notebook, P)

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
        self.reset_trajectory_frame()

        # And how to handle the structure
        if self.node.calculation == "nve":
            self["structure"].grid(row=row, column=0)
            row += 1

        return row

    def reset_control_frame(self, widget=None):
        """Layout the control widgets as needed for the current state.

        Which controls are shown comes from the parameters' rules
        (lammps_step.NVE_Parameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()

        # Clear out the previous widgets
        c_frame = self["control_frame"]
        for slave in c_frame.grid_slaves():
            slave.grid_forget()

        # And put them back in depending...
        row = 0
        widgets = []
        for key in ("time", "maximum_time", "timestep", "sampling"):
            if P.applies(key, values):
                widgets.append(self[key])
                self[key].grid(row=row, column=0, sticky=tk.W)
                row += 1
        if P.applies("control_properties", values):
            self["control_properties"].grid(row=row, column=0, sticky=tk.NSEW)
            row += 1

        sw.align_labels(widgets, sticky=tk.E)

    # The layout of the trajectory tab: each control and its level of indentation
    trajectory_layout = (
        ("trajectory", 0),
        ("trajectory number of samples", 1),
        ("trajectory rate", 1),
        ("trajectory forces", 1),
        ("trajectory velocities", 1),
        ("trajectory save", 1),
        ("trajectory system name", 2),
        ("make current", 2),
        ("trajectory export", 1),
        ("trajectory extxyz", 2),
        ("trajectory extxyz skip frames", 3),
        ("trajectory extxyz filename", 3),
        ("trajectory extxyz append", 3),
        *(
            item
            for what in (
                "atomic positions",
                "com positions",
                "atomic velocities",
                "com velocities",
                "heat flux",
                "pressure tensor",
            )
            for item in (
                (what, 0),
                (f"{what} number of samples", 1),
                (f"{what} rate", 1),
            )
        ),
        ("keep sampled files", 0),
    )

    # Controls not in the trajectory tab, though in the layout above
    trajectory_skip = ()

    def reset_trajectory_frame(self, widget=None):
        """Layout the trajectory frame according to its contents.

        Which controls are shown comes from the parameters' rules
        (lammps_step.NVE_Parameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()

        # Clear out the previous widgets
        frame = self["trajectory frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        # The widgets at each level of indentation
        levels = [[], [], [], []]
        for key, level in self.trajectory_layout:
            if key in self.trajectory_skip or key not in self:
                continue
            if not P.applies(key, values):
                continue
            self[key].grid(row=row, column=level, columnspan=5 - level, sticky=tk.EW)
            levels[level].append(self[key])
            row += 1

        # Indent each level from the one above it
        indent = 75
        width = sw.align_labels(levels[0], sticky=tk.E)
        for level in range(1, 4):
            if level > 1 and len(levels[level]) == 0:
                break
            width_level = sw.align_labels(levels[level], sticky=tk.E)
            frame.columnconfigure(
                level - 1, minsize=max(width - width_level + indent, 0)
            )
            width = width_level
        frame.columnconfigure(4, weight=1)

    # The extra commands and properties (lammps_step.extras)

    def _create_extras_tab(self, notebook, P):
        """The tab for the extra commands and the table of extra properties."""
        eframe = ttk.Frame(notebook)
        self["extras frame"] = eframe
        notebook.insert(self["results frame"], eframe, text="Extras", sticky=tk.NSEW)

        ttk.Label(
            eframe,
            text=(
                "Extra commands: LAMMPS compute, fix, variable and dump commands for "
                "this sub-step only, removed when it ends."
            ),
        ).grid(row=0, column=0, sticky=tk.W)
        text = self["extra commands text"] = Pmw.ScrolledText(
            eframe,
            labelpos="n",
            text_wrap="none",
            text_height=6,
            text_width=80,
        )
        text.insert("1.0", P["extra commands"].value or "")
        text.component("text").bind("<KeyRelease>", self._check_extras)
        text.grid(row=1, column=0, sticky=tk.NSEW)

        ttk.Label(
            eframe,
            text=(
                "Extra properties: LAMMPS expressions sampled and analysed with the "
                "other properties."
            ),
        ).grid(row=2, column=0, sticky=tk.W, pady=(10, 0))
        table = self["extra properties table"] = ttk.Frame(eframe)
        table.grid(row=3, column=0, sticky=tk.NSEW)
        self["extras message"] = ttk.Label(eframe, text="", foreground="red")
        self["extras message"].grid(row=4, column=0, sticky=tk.W)
        eframe.columnconfigure(0, weight=1)
        eframe.rowconfigure(1, weight=1)

        for row in P["extra properties"].value or []:
            self._add_extra_row(row, layout=False)
        self._layout_extra_properties()

    def _add_extra_row(self, data=None, layout=True):
        """Add a row to the extra properties table."""
        data = data or {}
        table = self["extra properties table"]
        row = {
            "name": ttk.Entry(table, width=12),
            "expression": ttk.Entry(table, width=24),
            "units": ttk.Combobox(table, width=14, values=list(DISPLAY_UNITS)),
            "treatment": ttk.Combobox(
                table, width=18, values=TREATMENTS, state="readonly"
            ),
            "prefactor": ttk.Combobox(table, width=12, values=PREFACTORS),
            "store": tk.BooleanVar(value=bool(data.get("store", False))),
        }
        row["store widget"] = ttk.Checkbutton(table, variable=row["store"])
        row["remove"] = ttk.Button(
            table, text="-", width=2, command=lambda r=row: self._remove_extra_row(r)
        )
        row["name"].insert(0, data.get("name", ""))
        row["expression"].insert(0, data.get("expression", ""))
        row["units"].set(data.get("units", "dimensionless") or "dimensionless")
        row["treatment"].set(data.get("treatment", "average") or "average")
        row["prefactor"].set(data.get("prefactor", "") or "V/kT")
        for key in ("name", "expression", "units", "prefactor"):
            row[key].bind("<KeyRelease>", self._check_extras)
        for key in ("units", "treatment", "prefactor"):
            row[key].bind("<<ComboboxSelected>>", self._layout_extra_properties)
        self._extra_rows.append(row)
        if layout:
            self._layout_extra_properties()

    def _remove_extra_row(self, row):
        """Remove a row from the extra properties table."""
        for key, widget in row.items():
            if isinstance(widget, tk.Widget):
                widget.destroy()
        self._extra_rows.remove(row)
        self._layout_extra_properties()

    def _layout_extra_properties(self, event=None):
        """Lay out the table: the prefactor only for a Green-Kubo integral."""
        table = self["extra properties table"]
        for slave in table.grid_slaves():
            slave.grid_forget()
        headers = ("Name", "Expression", "Units", "Treatment", "Prefactor", "Store")
        for column, text in enumerate(headers):
            if "header " + text not in self:
                self["header " + text] = ttk.Label(table, text=text)
            self["header " + text].grid(row=0, column=column, sticky=tk.W)
        for i, row in enumerate(self._extra_rows, start=1):
            row["name"].grid(row=i, column=0, sticky=tk.EW)
            row["expression"].grid(row=i, column=1, sticky=tk.EW)
            row["units"].grid(row=i, column=2, sticky=tk.EW)
            row["treatment"].grid(row=i, column=3, sticky=tk.EW)
            if row["treatment"].get() == "Green-Kubo integral":
                row["prefactor"].grid(row=i, column=4, sticky=tk.EW)
            row["store widget"].grid(row=i, column=5)
            row["remove"].grid(row=i, column=6)
        if "add extra property" not in self:
            self["add extra property"] = ttk.Button(
                table, text="+", width=2, command=self._add_extra_row
            )
        self["add extra property"].grid(
            row=len(self._extra_rows) + 1, column=6, sticky=tk.W
        )
        table.columnconfigure(1, weight=1)
        self._check_extras()

    def _read_extras(self):
        """The extra commands and the rows of the table, as the parameters hold
        them."""
        commands = self["extra commands text"].get(1.0, tk.END).rstrip()
        rows = []
        for row in self._extra_rows:
            treatment = row["treatment"].get()
            rows.append(
                {
                    "name": row["name"].get().strip(),
                    "expression": row["expression"].get().strip(),
                    "units": row["units"].get().strip(),
                    "treatment": treatment,
                    "prefactor": (
                        row["prefactor"].get().strip()
                        if treatment == "Green-Kubo integral"
                        else ""
                    ),
                    "store": bool(row["store"].get()),
                }
            )
        return commands, rows

    def _check_extras(self, event=None):
        """Show any problems with the extra commands and properties."""
        if "extras message" not in self:
            return ""
        commands, rows = self._read_extras()
        problems = check({"extra commands": commands, "extra properties": rows})
        self["extras message"].configure(text=problems)
        return problems

    def handle_dialog(self, result):
        if result == "OK":
            # Shortcut for parameters
            P = self.node.parameters

            if "extras frame" in self:
                commands, rows = self._read_extras()
                P["extra commands"].value = commands
                P["extra properties"].value = [
                    r for r in rows if r["name"] != "" or r["expression"] != ""
                ]

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
