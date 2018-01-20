# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS minimization step"""

import lammps_step
import molssi_util.molssi_widgets as mw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk


class TkMinimization(lammps_step.TkEnergy):
    def __init__(self, node=None, canvas=None, x=None, y=None, w=None, h=None):
        '''Initialize a node

        Keyword arguments:
        '''

        super().__init__(node=node, canvas=canvas, x=x, y=y, w=w, h=h)

    def create_dialog(self):
        """Create the dialog!"""

        super().create_dialog()

        frame = self.dialog.interior()

        # Convergence
        convergence_label = ttk.Label(frame, text='Simulation convergence:')
        self.tmp['convergence_label'] = convergence_label

        convergence = ttk.Combobox(
            frame, state='readonly',
            values=['normal',
                    'tight',
                    'loose',
                    'crude',
                    'as given',
                    'from a variable',
                    ],
            justify=tk.LEFT, width=20
        )
        convergence.set(self.node.convergence)
        self.tmp['convergence'] = convergence

        # etol: a variable or a value...
        etol_label = ttk.Label(frame, text='change in energy:')
        self.tmp['etol_label'] = etol_label

        etol_method = ttk.Combobox(
            frame, state='readonly',
            values=['as given',
                    'from a variable',
                    ],
            justify=tk.LEFT, width=20
        )
        etol_method.set(self.node.etol_method)
        self.tmp['etol_method'] = etol_method

        etol_variable = ttk.Entry(frame, width=15)
        etol_variable.insert(0, self.node.etol_variable)
        self.tmp['etol_variable'] = etol_variable

        etol = ttk.Entry(frame, width=15)
        etol.insert(0, self.node.etol)
        self.tmp['etol'] = etol

        # ftol: a variable or a value...
        ftol_label = ttk.Label(frame, text='change in energy:')
        self.tmp['ftol_label'] = ftol_label

        ftol_method = ttk.Combobox(
            frame, state='readonly',
            values=['as given',
                    'from a variable',
                    ],
            justify=tk.LEFT, width=20
        )
        ftol_method.set(self.node.ftol_method)
        self.tmp['ftol_method'] = ftol_method

        ftol_variable = ttk.Entry(frame, width=15)
        ftol_variable.insert(0, self.node.ftol_variable)
        self.tmp['ftol_variable'] = ftol_variable

        ftol = mw.UnitEntry(frame, width=15)
        ftol.set(self.node.ftol)
        self.tmp['ftol'] = ftol

        # maxiters: a variable or a value...
        maxiters_label = ttk.Label(frame, text='change in energy:')
        self.tmp['maxiters_label'] = maxiters_label

        maxiters_method = ttk.Combobox(
            frame, state='readonly',
            values=['default',
                    'as given',
                    'from a variable',
                    ],
            justify=tk.LEFT, width=20
        )
        maxiters_method.set(self.node.maxiters_method)
        self.tmp['maxiters_method'] = maxiters_method

        maxiters_variable = ttk.Entry(frame, width=15)
        maxiters_variable.insert(0, self.node.maxiters_variable)
        self.tmp['maxiters_variable'] = maxiters_variable

        maxiters = ttk.Entry(frame, width=15)
        maxiters.insert(0, self.node.maxiters)
        self.tmp['maxiters'] = maxiters

        # maxevals: a variable or a value...
        maxevals_label = ttk.Label(frame, text='change in energy:')
        self.tmp['maxevals_label'] = maxevals_label

        maxevals_method = ttk.Combobox(
            frame, state='readonly',
            values=['default',
                    'as given',
                    'from a variable',
                    ],
            justify=tk.LEFT, width=20
        )
        maxevals_method.set(self.node.maxevals_method)
        self.tmp['maxevals_method'] = maxevals_method

        maxevals_variable = ttk.Entry(frame, width=15)
        maxevals_variable.insert(0, self.node.maxevals_variable)
        self.tmp['maxevals_variable'] = maxevals_variable

        maxevals = ttk.Entry(frame, width=15)
        maxevals.insert(0, self.node.maxevals)
        self.tmp['maxevals'] = maxevals

        convergence.bind("<<ComboboxSelected>>", self.reset_dialog)
        etol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        ftol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxiters_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxevals_method.bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        w = self.tmp
        convergence = w['convergence'].get()
        etol_method = w['etol_method'].get()
        ftol_method = w['ftol_method'].get()
        maxiters_method = w['maxiters_method'].get()
        maxevals_method = w['maxevals_method'].get()

        frame = self.dialog.interior()
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['convergence_label'].grid(row=row, column=0, sticky=tk.E)
        w['convergence'].grid(row=row, column=1, sticky=tk.W)

        if 'energy' in convergence:
            row += 1
            w['etol_label'].grid(row=row, column=0, sticky=tk.E)
            w['etol_method'].grid(row=row, column=1, sticky=tk.EW)
            if etol_method == 'from a variable':
                w['etol_variable'].grid(row=row, column=1, sticky=tk.W)
            elif etol_method == 'give a value':
                w['etol'].grid(row=row, column=1, sticky=tk.W)

        if 'forces' in convergence:
            row += 1
            w['ftol_label'].grid(row=row, column=0, sticky=tk.E)
            w['ftol_method'].grid(row=row, column=1, sticky=tk.EW)
            if ftol_method == 'from a variable':
                w['ftol_variable'].grid(row=row, column=1, sticky=tk.W)
            elif ftol_method == 'give a value':
                w['ftol'].grid(row=row, column=1, sticky=tk.W)

        if 'energy' in convergence or 'forces' in convergence:
            row += 1
            w['maxiters_label'].grid(row=row, column=0, sticky=tk.E)
            w['maxiters_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxiters_method == 'from a variable':
                w['maxiters_variable'].grid(row=row, column=1, sticky=tk.W)
            elif maxiters_method == 'give a value':
                w['maxiters'].grid(row=row, column=1, sticky=tk.W)

            row += 1
            w['maxevals_label'].grid(row=row, column=0, sticky=tk.E)
            w['maxevals_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxevals_method == 'from a variable':
                w['maxevals_variable'].grid(row=row, column=1, sticky=tk.W)
            elif maxevals_method == 'give a value':
                w['maxevals'].grid(row=row, column=1, sticky=tk.W)

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
        convergence = w['convergence'].get()
        etol_method = w['etol_method'].get()
        ftol_method = w['ftol_method'].get()
        maxiters_method = w['maxiters_method'].get()
        maxevals_method = w['maxevals_method'].get()

        self.node.convergence = convergence
        if 'energy' in convergence:
            self.node.etol_method = etol_method
            if etol_method == 'from a variable':
                self.node.etol_variable = w['etol_variable'].get()
            elif etol_method == 'give a value':
                self.node.etol = int(w['etol'].get())

        if 'forces' in convergence:
            self.node.ftol_method = ftol_method
            if ftol_method == 'from a variable':
                self.node.ftol_variable = w['ftol_variable'].get()
            elif ftol_method == 'give a value':
                self.node.ftol = w['ftol'].get()

        if 'energy' in convergence or 'forces' in convergence:
            self.node.maxiters_method = maxiters_method
            if maxiters_method == 'from a variable':
                self.node.maxtiters_variable = w['maxiters_variable'].get()
            elif maxiters_method == 'give a value':
                self.node.maxiters = int(w['maxiters'].get())

            self.node.maxevals_method = maxevals_method
            if maxevals_method == 'from a variable':
                self.node.maxevals_variable = w['maxevals_variable'].get()
            elif maxevals_method == 'give a value':
                self.node.maxevals = int(w['maxevals'].get())
