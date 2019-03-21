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

        # how to set the velocities
        method_label = ttk.Label(frame, text='Set the velocities')
        self['method_label'] = method_label
        method = ttk.Combobox(
            frame,
            state='readonly',
            values=['using a random distribution', 'by scaling them']
        )
        method.set(self.node.method)
        self['method'] = method

        # the temperature
        temperature_method_label = ttk.Label(frame, text='Temperature')
        self['temperature_method_label'] = temperature_method_label
        temperature_method = ttk.Combobox(
            frame, width=10,
            state='readonly',
            values=['is', 'from variable']
        )
        temperature_method.set(self.node.temperature_method)
        self['temperature_method'] = temperature_method
        temperature = mw.UnitEntry(frame, width=10)
        temperature.set(self.node.temperature)
        self['temperature'] = temperature
        temperature_variable = ttk.Entry(frame, width=15)
        temperature_variable.insert(0, self.node.temperature_variable)
        self['temperature_variable'] = temperature_variable

        # the random seed
        seed_method_label = ttk.Label(frame, text='Random seed')
        self['seed_method_label'] = seed_method_label
        seed_method = ttk.Combobox(
            frame, width=10,
            state='readonly',
            values=['random', 'is', 'from variable']
        )
        seed_method.set(self.node.seed_method)
        self['seed_method'] = seed_method
        seed = ttk.Entry(frame, width=10)
        seed.insert(0, self.node.seed)
        self['seed'] = seed
        seed_variable = ttk.Entry(frame, width=15)
        seed_variable.insert(0, self.node.seed_variable)
        self['seed_variable'] = seed_variable

        # zero translational momentum
        self.tk_var['remove_translations'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            self.tk_var['remove_translations'].set(1)
        else:
            self.tk_var['remove_translations'].set(0)
        remove_translations = ttk.Checkbutton(
            frame,
            variable=self.tk_var['remove_translations'],
            text='remove any translational momentum'
        )
        self['remove_translations'] = remove_translations

        # zero rotational momentum
        self.tk_var['remove_rotations'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            self.tk_var['remove_rotations'].set(1)
        else:
            self.tk_var['remove_rotations'].set(0)
        remove_rotations = ttk.Checkbutton(
            frame,
            variable=self.tk_var['remove_rotations'],
            text='remove any rotational momentum'
        )
        self['remove_rotations'] = remove_rotations

        self['method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        self['temperature_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        self['seed_method'].bind("<<ComboboxSelected>>", self.reset_dialog)

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS velocities"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Lay out the widgets given the current state"""
        method = self['method'].get()
        temperature_method = self['temperature_method'].get()
        seed_method = self['seed_method'].get()

        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        self['method_label'].grid(row=row, column=0, sticky=tk.E)
        self['method'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            self['temperature_method_label'].grid(row=row, column=0, sticky=tk.E)
            self['temperature_method'].grid(row=row, column=1, sticky=tk.EW)
            if temperature_method == 'is':
                self['temperature'].grid(row=row, column=2, sticky=tk.W)
            else:
                self['temperature_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method:
            row += 1
            self['seed_method_label'].grid(row=row, column=0, sticky=tk.E)
            self['seed_method'].grid(row=row, column=1, sticky=tk.EW)
            if seed_method == 'is':
                self['seed'].grid(row=row, column=2, sticky=tk.W)
            elif 'variable' in seed_method:
                self['seed_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            self['remove_translations'].grid(row=row, column=1, columnspan=2,
                                          sticky=tk.W)
            row += 1
            self['remove_rotations'].grid(row=row, column=1, columnspan=2,
                                       sticky=tk.W)

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

        method = self['method'].get()
        temperature_method = self['temperature_method'].get()
        seed_method = self['seed_method'].get()

        self.node.method = method
        if 'random' in method or 'scaling' in method:
            self.node.temperature_method = temperature_method
            if temperature_method == 'is':
                self.node.temperature['temperature'].get()
            else:
                self.node.temperature_variable = \
                    self['temperature_variable'].get()

        if 'random' in method:
            self.node.seed_method = seed_method
            if seed_method == 'is':
                self.node.seed['seed'].get()
            else:
                self.node.seed_variable = \
                    self['seed_variable'].get()

        if 'random' in method or 'scaling' in method:
            self.node.remove_linear_momentum = \
                (self.tk_var['remove_translations'] == 1)
            self.node.remove_angular_momentum = \
                (self.tk_var['remove_rotations'] == 1)
