# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import molssi_util.molssi_widgets as mw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk


class TkNVE(lammps_step.TkEnergy):
    def __init__(self, tk_workflow=None, node=None, canvas=None,
                 x=None, y=None, w=None, h=None):
        '''Initialize a node

        Keyword arguments:
        '''

        super().__init__(tk_workflow=tk_workflow, node=node,
                         canvas=canvas, x=x, y=y, w=w, h=h)

    def create_dialog(self):
        """Create the dialog!"""

        super().create_dialog()

        self.dialog.configure(title='Edit NVE dynamics step')
        frame = self._widget['frame']

        # Simulation time...
        time_label = ttk.Label(frame, text='Simulation time:')
        self._widget['time_label'] = time_label

        time = mw.UnitEntry(frame, width=15)
        time.set(self.node.time)
        self._widget['time'] = time

        # Timestep: automatic, a variable or a value...
        timestep_label = ttk.Label(frame, text='Timestep:')
        self._widget['timestep_label'] = timestep_label

        timestep_method = ttk.Combobox(
            frame, state='readonly',
            values=['normal',
                    'accurate but slow',
                    'coarse but fast',
                    'from variable',
                    'is'
                    ],
            justify=tk.LEFT, width=15
        )
        timestep_method.set(self.node.timestep_method)
        self._widget['timestep_method'] = timestep_method

        timestep_variable = ttk.Entry(frame, width=15)
        timestep_variable.insert(0, self.node.timestep_variable)
        self._widget['timestep_variable'] = timestep_variable

        timestep = mw.UnitEntry(frame, width=15)
        timestep.set(self.node.timestep)
        self._widget['timestep'] = timestep

        # Sampling: automatic, a variable or a value...
        sampling_label = ttk.Label(frame, text='Sampling frequency:')
        self._widget['sampling_label'] = sampling_label

        sampling_method = ttk.Combobox(
            frame, state='readonly',
            values=['none',
                    'from variable',
                    'is'
                    ],
            justify=tk.LEFT, width=15
        )
        sampling_method.set(self.node.sampling_method)
        self._widget['sampling_method'] = sampling_method

        sampling_variable = ttk.Entry(frame, width=15)
        sampling_variable.insert(0, self.node.sampling_variable)
        self._widget['sampling_variable'] = sampling_variable

        sampling = mw.UnitEntry(frame, width=15)
        sampling.set(self.node.sampling)
        self._widget['sampling'] = sampling

        timestep_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        sampling_method.bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        w = self._widget
        timestep_method = w['timestep_method'].get()
        sampling_method = w['sampling_method'].get()

        frame = self._widget['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['time_label'].grid(row=row, column=0, sticky=tk.E)
        w['time'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        row += 1
        w['timestep_label'].grid(row=row, column=0, sticky=tk.E)
        w['timestep_method'].grid(row=row, column=1, sticky=tk.W)
        if timestep_method == 'from variable':
            w['timestep_variable'].grid(row=row, column=2, sticky=tk.W)
        elif timestep_method == 'is':
            w['timestep'].grid(row=row, column=2, sticky=tk.W)

        row += 1
        w['sampling_label'].grid(row=row, column=0, sticky=tk.E)
        w['sampling_method'].grid(row=row, column=1, sticky=tk.W)
        if sampling_method == 'from variable':
            w['sampling_variable'].grid(row=row, column=2, sticky=tk.W)
        elif sampling_method == 'is':
            w['sampling'].grid(row=row, column=2, sticky=tk.W)

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
        timestep_method = w['timestep_method'].get()
        sampling_method = w['sampling_method'].get()

        self.node.time = w['time'].get()

        self.node.timestep_method = timestep_method
        if timestep_method == 'from variable':
            self.node.timestep_variable = w['timestep_variable'].get()
        elif timestep_method == 'is':
            self.node.timestep = w['timestep'].get()

        self.node.sampling_method = sampling_method
        if sampling_method == 'from variable':
            self.node.sampling_method = w['sampling_variable'].get()
        elif sampling_method == 'is':
            self.node.sampling = w['sampling'].get()
