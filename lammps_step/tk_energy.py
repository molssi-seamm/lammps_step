# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Energy step"""

import molssi_workflow
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkEnergy(molssi_workflow.TkNode):
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
            title='Edit Energy parameters',
            command=self.handle_dialog)
        self.dialog.withdraw()

        frame = ttk.Frame(self.dialog.interior())
        frame.pack(expand=tk.YES, fill=tk.BOTH)
        self['frame'] = frame

        self['message'] = ttk.Label(
            frame,
            text='The LAMMPS energy step has no parameters\n'
            'All relevant parameters are set in the initialization step.'
        )
        self['message'].grid()

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS energy
        calculation"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Layout the dialog according to the current control
        parameters. For the energy, there are no parameters, so
        do nothing."""

        pass

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
