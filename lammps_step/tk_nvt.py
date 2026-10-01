# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS NVT dynamics step"""

import lammps_step
import logging
import seamm_widgets as sw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk

logger = logging.getLogger(__name__)


class TkNVT(lammps_step.TkNVE):
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

        # Set the logging level for this module if requested
        # if 'lammps_tk_nvt_log_level' in self.options:
        #     logger.setLevel(self.options.lammps_tk_nvt_log_level)
        #     logger.critical(
        #         'Set log level to {}'.format(
        #             self.options.lammps_tk_nvt_log_level
        #         )
        #     )

        # Call the constructor for the superclass
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

    def create_dialog(self, title="Edit NVT dynamics parameters"):
        """Create the edit dialog!

        This is reasonably complicated, so a bit of description
        is in order. The superclass NVE creates the dialog
        along with the basic runtime and timestep.

        This method adds a second frame for setting the temperature and
        thermostat (this is the T part of NVT).

        The layout is handled in part by the NVT superclass, which
        handles the temperature frame. Our part is handled by two
        methods:

        * reset_dialog does the general layout of the large blocks
        * reset_temperature_frame handles the layout of the temperature
          section except
        """

        logger.debug("TkNVT.create_dialog")

        # Let parent classes do their thing.
        super().create_dialog(title=title)

        # Shortcut for parameters
        P = self.node.parameters

        logger.debug("Parameters:\n{}".format(pprint.pformat(P.to_dict())))

        # Temperature frame to isolate widgets
        t_frame = self["temperature_frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Temperature",
            labelanchor="n",
            padding=10,
        )

        for key in lammps_step.NVT_Parameters.parameters:
            self[key] = P[key].widget(t_frame)

        # and binding to change as needed
        self["thermostat"].combobox.bind(
            "<<ComboboxSelected>>", self.reset_temperature_frame
        )

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        row = super().reset_dialog()

        self["temperature_frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_temperature_frame()

        # And how to handle the structure
        if self.node.calculation == "nvt":
            self["structure"].grid(row=row, column=0)
            row += 1

        return row

    def reset_temperature_frame(self, widget=None):
        """Layout the widgets in the temperature frame
        as needed for the current state.

        Which controls are shown comes from the parameters' rules
        (lammps_step.NVT_Parameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()

        t_frame = self["temperature_frame"]
        for slave in t_frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # Main controls
        for key in ("T0", "T1", "thermostat"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
            row += 1

        sw.align_labels((self["T0"], self["T1"], self["thermostat"]), sticky=tk.E)

        # and controls for specific thermostats, indented
        widgets = []
        for key in (
            "Tdamp",
            "drag",
            "Tchain",
            "Tloop",
            "seed",
            "frequency",
            "window",
            "fraction",
        ):
            if P.applies(key, values):
                self[key].grid(row=row, column=1, sticky=tk.W)
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky=tk.E)

        t_frame.columnconfigure(0, minsize=50)
