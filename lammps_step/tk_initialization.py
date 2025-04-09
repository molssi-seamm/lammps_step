# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Initialization step"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import seamm
import seamm_widgets as sw

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
        self.general_cb()
        self.kspace_method_cb()

    def general_cb(self, event=None):
        """Grid the widgets into the dialog, depending on the current values
        of key variables. This provides a dyamic presentation to the user.
        """
        frame = self["general"]
        charge_method = self["atomic charges"].get()

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

        if charge_method.startswith("charge "):
            for key in ("qeq convergence", "qeq iterations"):
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

        # Remove any widgets previously packed
        for slave in self["periodic"].grid_slaves():
            slave.grid_forget()

        row = 0
        self["kspace_method"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        method = self["kspace_method"].get()
        if method != "none":
            self["kspace_accuracy"].grid(row=row, column=0, sticky=tk.W)
            row += 1

            if "few charged" in method or method[0] == "$":
                self["kspace_smallq"].grid(row=row, column=0, sticky=tk.W)
                row += 1

        if method[0] == "$" or "dispersion" not in method:
            self["tail_correction"].grid(row=row, column=0, sticky=tk.W)
            row += 1
