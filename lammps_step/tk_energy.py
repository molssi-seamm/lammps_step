# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Energy step"""

import chemflowchart
import lammps_step
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkEnergy(chemflowchart.TkNode):
    def __init__(self, node=None, canvas=None, x=None, y=None, w=None, h=None):
        '''Initialize a node

        Keyword arguments:
        '''

        self.dialog = None
        self.tmp = {}

        super().__init__(node=node, canvas=canvas, x=x, y=y, w=w, h=h)

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
            title='Edit Energy step',
            command=self.handle_dialog)
        self.dialog.withdraw()

        self.tmp = {}
        frame = self.dialog.interior()

        # which structure? may need to set default first...
        if self.node.structure is None:
            if isinstance(self.node.previous(), chemflowchart.StartNode):
                self.node.structure = 'initial'
            else:
                self.node.structure = 'current'

        structure_label = ttk.Label(frame, text='Structure:')
        structure = ttk.Combobox(
            frame,
            state='readonly',
            values=list(lammps_step.Energy.structures))
        structure.set(self.node.structure)
        self.tmp['structure'] = structure

        structure_label.grid(row=0, column=0, columnspan=2, sticky=tk.E)
        structure.grid(row=0, column=2, sticky=tk.W)

        frame.grid_columnconfigure(0, minsize=30)

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS energy
        calculation"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        pass

    def handle_dialog(self, result):
        if result == 'Cancel':
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

        self.node.structure = self.tmp['structure'].get()
        self.tmp = None
