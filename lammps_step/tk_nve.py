# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import molssi_util.molssi_widgets as mw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk


class TkNVE(lammps_step.TkEnergy):
    def __init__(self, node=None, canvas=None, x=None, y=None, w=None, h=None):
        '''Initialize a node

        Keyword arguments:
        '''

        super().__init__(node=node, canvas=canvas, x=x, y=y, w=w, h=h)

    def create_dialog(self):
        """Create the dialog!"""

        super().create_dialog()

        frame = self.dialog.interior()

        # Simulation time...
        time_label = ttk.Label(frame, text='Simulation time:')
        self.tmp['time_label'] = time_label

        time = mw.UnitEntry(frame, width=15)
        time.set(self.node.time)
        self.tmp['time'] = time

        # Timestep: automatic, a variable or a value...
        timestep_label = ttk.Label(frame, text='Timestep:')
        self.tmp['timestep_label'] = timestep_label

        timestep_method = ttk.Combobox(
            frame, state='readonly',
            values=['normal',
                    'accurate but slow',
                    'coarse but fast',
                    'from a variable',
                    'give a value'
                    ],
            justify=tk.RIGHT, width=28
        )
        timestep_method.set(self.node.timestep_method)
        self.tmp['timestep_method'] = timestep_method

        timestep_variable = ttk.Entry(frame, width=15)
        timestep_variable.insert(0, self.node.timestep_variable)
        self.tmp['timestep_variable'] = timestep_variable

        timestep = mw.UnitEntry(frame, width=15)
        timestep.set(self.node.timestep)
        self.tmp['timestep'] = timestep

        # Sampling: automatic, a variable or a value...
        sampling_label = ttk.Label(frame, text='Sampling frequency:')
        self.tmp['sampling_label'] = sampling_label

        sampling_method = ttk.Combobox(
            frame, state='readonly',
            values=['none',
                    'from a variable',
                    'give a value'
                    ],
            justify=tk.RIGHT, width=28
        )
        sampling_method.set(self.node.sampling_method)
        self.tmp['sampling_method'] = sampling_method

        sampling_variable = ttk.Entry(frame, width=15)
        sampling_variable.insert(0, self.node.sampling_variable)
        self.tmp['sampling_variable'] = sampling_variable

        sampling = mw.UnitEntry(frame, width=15)
        sampling.set(self.node.sampling)
        self.tmp['sampling'] = sampling

        timestep_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        sampling_method.bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        w = self.tmp
        timestep_method = w['timestep_method'].get()
        sampling_method = w['sampling_method'].get()

        frame = self.dialog.interior()
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['time_label'].grid(row=row, column=0, sticky=tk.E)
        w['time'].grid(row=row, column=1, sticky=tk.W)

        row += 1
        w['timestep_label'].grid(row=row, column=0, sticky=tk.E)
        w['timestep_method'].grid(row=row, column=1, sticky=tk.EW)
        if timestep_method == 'from a variable':
            w['timestep_variable'].grid(row=row, column=2, sticky=tk.W)
        elif timestep_method == 'give a value':
            w['timestep'].grid(row=row, column=2, sticky=tk.W)

        row += 1
        w['sampling_label'].grid(row=row, column=0, sticky=tk.E)
        w['sampling_method'].grid(row=row, column=1, sticky=tk.EW)
        if sampling_method == 'from a variable':
            w['sampling_variable'].grid(row=row, column=2, sticky=tk.W)
        elif sampling_method == 'give a value':
            w['sampling'].grid(row=row, column=2, sticky=tk.W)

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
        timestep_method = w['timestep_method'].get()
        sampling_method = w['sampling_method'].get()

        self.node.time = w['time'].get()

        self.node.timestep_method = timestep_method
        if timestep_method == 'from a variable':
            self.node.timestep_variable = w['timestep_variable'].get()
        elif timestep_method == 'give a value':
            self.node.timestep = w['timestep'].get()

        self.node.sampling_method = sampling_method
        if sampling_method == 'from a variable':
            self.node.sampling_method = w['sampling_variable'].get()
        elif sampling_method == 'give a value':
            self.node.sampling = w['sampling'].get()
