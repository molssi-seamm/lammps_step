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

        w = self._widget

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
        w['frame'] = frame

        # how to set the velocities
        method_label = ttk.Label(frame, text='Set the velocities')
        w['method_label'] = method_label
        method = ttk.Combobox(
            frame,
            state='readonly',
            values=['using a random distribution', 'by scaling them']
        )
        method.set(self.node.method)
        w['method'] = method

        # the temperature
        temperature_method_label = ttk.Label(frame, text='Temperature')
        w['temperature_method_label'] = temperature_method_label
        temperature_method = ttk.Combobox(
            frame, width=10,
            state='readonly',
            values=['is', 'from variable']
        )
        temperature_method.set(self.node.temperature_method)
        w['temperature_method'] = temperature_method
        temperature = mw.UnitEntry(frame, width=10)
        temperature.set(self.node.temperature)
        w['temperature'] = temperature
        temperature_variable = ttk.Entry(frame, width=15)
        temperature_variable.insert(0, self.node.temperature_variable)
        w['temperature_variable'] = temperature_variable

        # the random seed
        seed_method_label = ttk.Label(frame, text='Random seed')
        w['seed_method_label'] = seed_method_label
        seed_method = ttk.Combobox(
            frame, width=10,
            state='readonly',
            values=['random', 'is', 'from variable']
        )
        seed_method.set(self.node.seed_method)
        w['seed_method'] = seed_method
        seed = ttk.Entry(frame, width=10)
        seed.insert(0, self.node.seed)
        w['seed'] = seed
        seed_variable = ttk.Entry(frame, width=15)
        seed_variable.insert(0, self.node.seed_variable)
        w['seed_variable'] = seed_variable

        # zero translational momentum
        w['remove_translations_value'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            w['remove_translations_value'].set(1)
        else:
            w['remove_translations_value'].set(0)
        remove_translations = ttk.Checkbutton(
            frame,
            variable=w['remove_translations_value'],
            text='remove any translational momentum'
        )
        w['remove_translations'] = remove_translations

        # zero rotational momentum
        w['remove_rotations_value'] = tk.IntVar()
        if self.node.remove_linear_momentum:
            w['remove_rotations_value'].set(1)
        else:
            w['remove_rotations_value'].set(0)
        remove_rotations = ttk.Checkbutton(
            frame,
            variable=w['remove_rotations_value'],
            text='remove any rotational momentum'
        )
        w['remove_rotations'] = remove_rotations

        w['method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['temperature_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['seed_method'].bind("<<ComboboxSelected>>", self.reset_dialog)

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS velocities"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Lay out the widgets given the current state"""
        w = self._widget
        method = w['method'].get()
        temperature_method = w['temperature_method'].get()
        seed_method = w['seed_method'].get()

        frame = w['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['method_label'].grid(row=row, column=0, sticky=tk.E)
        w['method'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            w['temperature_method_label'].grid(row=row, column=0, sticky=tk.E)
            w['temperature_method'].grid(row=row, column=1, sticky=tk.EW)
            if temperature_method == 'is':
                w['temperature'].grid(row=row, column=2, sticky=tk.W)
            else:
                w['temperature_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method:
            row += 1
            w['seed_method_label'].grid(row=row, column=0, sticky=tk.E)
            w['seed_method'].grid(row=row, column=1, sticky=tk.EW)
            if seed_method == 'is':
                w['seed'].grid(row=row, column=2, sticky=tk.W)
            elif 'variable' in seed_method:
                w['seed_variable'].grid(row=row, column=2, sticky=tk.W)

        if 'random' in method or 'scaling' in method:
            row += 1
            w['remove_translations'].grid(row=row, column=1, columnspan=2,
                                          sticky=tk.W)
            row += 1
            w['remove_rotations'].grid(row=row, column=1, columnspan=2,
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

        w = self._widget
        method = w['method'].get()
        temperature_method = w['temperature_method'].get()
        seed_method = w['seed_method'].get()

        self.node.method = method
        if 'random' in method or 'scaling' in method:
            self.node.temperature_method = temperature_method
            if temperature_method == 'is':
                self.node.temperature['temperature'].get()
            else:
                self.node.temperature_variable = \
                    w['temperature_variable'].get()

        if 'random' in method:
            self.node.seed_method = seed_method
            if seed_method == 'is':
                self.node.seed['seed'].get()
            else:
                self.node.seed_variable = \
                    w['seed_variable'].get()

        if 'random' in method or 'scaling' in method:
            self.node.remove_linear_momentum = \
                (w['remove_translations_value'] == 1)
            self.node.remove_angular_momentum = \
                (w['remove_rotations_value'] == 1)
