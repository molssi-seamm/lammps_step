# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS velocities step"""

import molssi_workflow
import molssi_util.molssi_widgets as mw
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkVelocities(molssi_workflow.TkNode):
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
            title='Edit velocity step',
            command=self.handle_dialog)
        self.dialog.withdraw()

        self.tmp = {}
        frame = self.dialog.interior()

        # how to set the velocities
        method_label = ttk.Label(frame, text='Set the velocities')
        self.tmp['method_label'] = method_label
        method = ttk.Combobox(
            frame,
            state='readonly',
            values=['using a random distribution', 'by scaling them']
        )
        method.set(self.node.method)
        self.tmp['method'] = method

        # the temperature
        temperature_method_label = ttk.Label(frame, text='Temperature')
        self.tmp['temperature_method_label'] = temperature_method_label
        temperature_method = ttk.Combobox(
            frame,
            state='readonly',
            values=['as given', 'from variable']
        )
        temperature_method.set(self.node.temperature_method)
        self.tmp['temperature_method'] = temperature_method
        temperature = mw.UnitEntry(frame, width=15)
        temperature.set(self.node.temperature)
        self.tmp['temperature'] = temperature
        temperature_variable = ttk.Entry(frame, width=15)
        temperature_variable.insert(0, self.node.temperature_variable)
        self.tmp['temperature_variable'] = temperature_variable

        # the random seed
        seed_method_label = ttk.Label(frame, text='Random seed')
        self.tmp['seed_method_label'] = seed_method_label
        seed_method = ttk.Combobox(
            frame,
            state='readonly',
            values=['random', 'as given', 'from variable']
        )
        seed_method.set(self.node.seed_method)
        self.tmp['seed_method'] = seed_method
        seed = ttk.Entry(frame, width=15)
        seed.insert(0, self.node.seed)
        self.tmp['seed'] = seed
        seed_variable = ttk.Entry(frame, width=15)
        seed_variable.insert(0, self.node.seed_variable)
        self.tmp['seed_variable'] = seed_variable

        # zero translational momentum
        self.tmp['remove_translations_value'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            self.tmp['remove_translations_value'].set(1)
        else:
            self.tmp['remove_translations_value'].set(0)
        remove_translations = ttk.Checkbutton(
            frame,
            variable=self.tmp['remove_translations_value'],
            text='remove any translational momentum'
        )
        self.tmp['remove_translations'] = remove_translations

        # zero rotational momentum
        self.tmp['remove_rotations_value'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            self.tmp['remove_rotations_value'].set(1)
        else:
            self.tmp['remove_rotations_value'].set(0)
        remove_rotations = ttk.Checkbutton(
            frame,
            variable=self.tmp['remove_rotations_value'],
            text='remove any rotational momentum'
        )
        self.tmp['remove_rotations'] = remove_rotations

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS velocities"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Lay out the widgets given the current state"""
        w = self.tmp
        method = w['method'].get()
        temperature_method = w['temperature_method'].get()
        seed_method = w['seed_method'].get()

        frame = self.dialog.interior()
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['method_label'].grid(row=row, column=0, sticky=tk.E)
        w['method'].grid(row=row, column=1, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            w['temperature_method_label'].grid(row=row, column=0, sticky=tk.E)
            w['temperature_method'].grid(row=row, column=1, sticky=tk.EW)
            if temperature_method == 'as given':
                w['temperature'].grid(row=row, column=2, sticky=tk.W)
            else:
                w['temperature_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method:
            row += 1
            w['seed_method_label'].grid(row=row, column=0, sticky=tk.E)
            w['seed_method'].grid(row=row, column=1, sticky=tk.EW)
            if seed_method == 'as given':
                w['seed'].grid(row=row, column=2, sticky=tk.W)
            else:
                w['seed_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            w['remove_translations'].grid(row=row, column=1, sticky=tk.W)
            row += 1
            w['remove_rotations'].grid(row=row, column=1, sticky=tk.W)

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

        w = self.tmp
        method = w['method'].get()
        temperature_method = w['temperature_method'].get()
        seed_method = w['seed_method'].get()

        self.node.method = method
        if 'random' in method or 'scaling' in method:
            self.node.temperature_method = temperature_method
            if temperature_method == 'as given':
                self.node.temperature['temperature'].get()
            else:
                self.node.temperature_variable = \
                    w['temperature_variable'].get()

        if 'random' in method:
            self.node.seed_method = seed_method
            if seed_method == 'as given':
                self.node.seed['seed'].get()
            else:
                self.node.seed_variable = \
                    w['seed_variable'].get()

        if 'random' in method or 'scaling' in method:
            self.node.remove_linear_momentum = \
                (w['remove_translations_value'] == 1)
            self.node.remove_angular_momentum = \
                (w['remove_rotations_value'] == 1)
