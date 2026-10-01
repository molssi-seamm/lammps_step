# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Initialization step"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import seamm
import seamm_widgets as sw

from .tk_energy import widget_values

logger = logging.getLogger(__name__)


class TkInitialization(seamm.TkNode):
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
        """Setup  the LAMMPS initialization node.

        Keyword arguments:
        """

        # # Set the logging level for this module if requested
        # if 'lammps_tk_initialization_log_level' in self.options:
        #     logger.setLevel(self.options.lammps_tk_initialization_log_level)
        #     logger.critical(
        #         'Set log level to {}'.format(
        #             self.options.lammps_tk_initialization_log_level
        #         )
        #     )

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

    def create_dialog(self):
        """Create the dialog!"""
        frame = super().create_dialog("Edit LAMMPS Initialization Step")

        # Shortcut for parameters
        P = self.node.parameters

        # Frame for options for all systems
        general = self["general"] = ttk.LabelFrame(
            frame,
            text="For all systems",
            relief=tk.SUNKEN,
            borderwidth=5,
            labelanchor=tk.N,
        )

        for key in (
            "cutoff",
            "shift_nonbond",
            "atomic charges",
            "qeq convergence",
            "qeq iterations",
        ):
            self[key] = P[key].widget(general)

        # Frame for the periodic system options, i.e. kspace, etc.
        periodic = self["periodic"] = ttk.LabelFrame(
            frame,
            text="For periodic systems",
            relief=tk.SUNKEN,
            borderwidth=5,
            labelanchor=tk.N,
        )

        for key in (
            "tail_correction",
            "kspace_method",
            "kspace_accuracy",
            "kspace_smallq",
        ):
            self[key] = P[key].widget(periodic)

        self["atomic charges"].bind("<<ComboboxSelected>>", self.general_cb)
        self["kspace_method"].bind("<<ComboboxSelected>>", self.kspace_method_cb)

        # Grid in the static part of the dialog

        self["general"].grid(row=0, column=0, sticky=tk.EW, pady=10)

        self["periodic"].grid(row=1, column=0, sticky=tk.EW, pady=10)
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out both parts of the dialog for the current values."""
        self.general_cb()
        self.kspace_method_cb()

    def _widget_values(self):
        """The dialog's current values, {name: value}, for the parameters' rules."""
        return widget_values(self)

    def general_cb(self, event=None):
        """Grid the widgets into the dialog, depending on the current values
        of key variables. This provides a dyamic presentation to the user.

        Which controls are shown comes from the parameters' rules
        (lammps_step.InitializationParameters), which the flowchart builder uses too.
        """
        P = self.node.parameters
        values = self._widget_values()
        frame = self["general"]

        # Remove any widgets previously packed
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        widgets2 = []
        for key in ("cutoff", "shift_nonbond", "atomic charges"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        # The controls for charge equilibration, indented
        for key in ("qeq convergence", "qeq iterations"):
            if P.applies(key, values):
                self[key].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self[key])
                row += 1

        w1 = sw.align_labels(widgets, sticky=tk.E)
        if len(widgets2) > 0:
            w2 = sw.align_labels(widgets2, sticky=tk.E)
            frame.columnconfigure(0, minsize=w1 - w2 + 30)

    def kspace_method_cb(self, event=None):
        """Grid the widgets into the dialog, depending on the current values
        of key variables. This provides a dyamic presentation to the user.
        """
        P = self.node.parameters
        values = self._widget_values()

        # Remove any widgets previously packed
        for slave in self["periodic"].grid_slaves():
            slave.grid_forget()

        row = 0
        for key in (
            "kspace_method",
            "kspace_accuracy",
            "kspace_smallq",
            "tail_correction",
        ):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky=tk.W)
                row += 1
