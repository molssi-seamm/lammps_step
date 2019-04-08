# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS velocities step"""

import molssi_workflow
import molssi_util.molssi_widgets as mw
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkVelocities(molssi_workflow.TkNode):
    def __init__(self, tk_workflow=None, node=None, canvas=None,
                 x=None, y=None, w=200, h=50):
        '''Initialize a node

        Keyword arguments:
        '''

        super().__init__(tk_workflow=tk_workflow, node=node,
                         canvas=canvas, x=x, y=y, w=w, h=h)

    def right_click(self, event):
        """Probably need to add our dialog...
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self):
        """Create the dialog!"""

        self.dialog = Pmw.Dialog(
            self.toplevel,
            buttons=('OK', 'Help', 'Cancel'),
            defaultbutton='OK',
            master=self.toplevel,
            title='Edit velocity step',
            command=self.handle_dialog)
        self.dialog.withdraw()

        frame = ttk.Frame(self.dialog.interior())
        frame.pack(expand=tk.YES, fill=tk.BOTH)
        self['frame'] = frame

        # Create all the widgets
        P = self.node.parameters
        for key in P:
            self[key] = P[key].widget(frame)

        self['method'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_dialog
        )
        self['method'].combobox.bind("<Return>", self.reset_dialog)
        self['method'].combobox.bind("<FocusOut>", self.reset_dialog)

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS velocities"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Lay out the widgets given the current state"""
        method = self['method'].get()

        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        for key in ('method', 'T', 'remove_momentum'):
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if 'scaling' in method:
            self[seed].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[seed])
            row += 1

            mw.alignlabels(widgets)

    def handle_dialog(self, result):
        if result is None or result == 'Cancel':
            self.dialog.deactivate(result)
            return

        if result == 'Help':
            # display help!!!
            return

        if result != "OK":
            self.dialog.deactivate(result)
            raise RuntimeError(
                "Don't recognize dialog result '{}'".format(result))

        self.dialog.deactivate(result)

        # Shortcut for parameters
        P = self.node.parameters

        for key in P:
            P[key].set_from_widget()
        
