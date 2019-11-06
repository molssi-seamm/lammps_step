# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS minimization step"""

import lammps_step
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk


class TkMinimization(lammps_step.TkEnergy):

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50
    ):
        '''Initialize a node

        Keyword arguments:
        '''

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h
        )

    def create_dialog(self):
        """Create the dialog!"""
        super().create_dialog()
        self.dialog.configure(title='Edit Minimization step')

        frame = self['frame']

        # Convergence
        convergence_label = ttk.Label(frame, text='Simulation convergence:')
        self['convergence_label'] = convergence_label

        convergence = ttk.Combobox(
            frame,
            state='readonly',
            values=[
                'normal',
                'tight',
                'loose',
                'crude',
                'on forces',
                'on energy and forces',
                'on just the energy change',
                'from variable',
            ],
            justify=tk.LEFT,
            width=20
        )
        convergence.set(self.node.convergence)
        self['convergence'] = convergence

        # etol: a variable or a value...
        etol_label = ttk.Label(frame, text='change in energy:')
        self['etol_label'] = etol_label

        etol_method = ttk.Combobox(
            frame,
            state='readonly',
            values=[
                'is',
                'from variable',
            ],
            justify=tk.LEFT,
            width=10
        )
        etol_method.set(self.node.etol_method)
        self['etol_method'] = etol_method

        etol_variable = ttk.Entry(frame, width=15)
        etol_variable.insert(0, self.node.etol_variable)
        self['etol_variable'] = etol_variable

        etol = ttk.Entry(frame, width=15)
        etol.insert(0, self.node.etol)
        self['etol'] = etol

        # ftol: a variable or a value...
        ftol_label = ttk.Label(frame, text='magnitude of the force vector:')
        self['ftol_label'] = ftol_label

        ftol_method = ttk.Combobox(
            frame,
            state='readonly',
            values=[
                'is',
                'from variable',
            ],
            justify=tk.LEFT,
            width=10
        )
        ftol_method.set(self.node.ftol_method)
        self['ftol_method'] = ftol_method

        ftol_variable = ttk.Entry(frame, width=15)
        ftol_variable.insert(0, self.node.ftol_variable)
        self['ftol_variable'] = ftol_variable

        ftol = sw.UnitEntry(frame, width=15)
        ftol.set(self.node.ftol)
        self['ftol'] = ftol

        # maxiters: a variable or a value...
        maxiters_label = ttk.Label(frame, text='maximum iterations:')
        self['maxiters_label'] = maxiters_label

        maxiters_method = ttk.Combobox(
            frame,
            state='readonly',
            values=[
                'default',
                'is',
                'from variable',
            ],
            justify=tk.LEFT,
            width=10
        )
        maxiters_method.set(self.node.maxiters_method)
        self['maxiters_method'] = maxiters_method

        maxiters_variable = ttk.Entry(frame, width=15)
        maxiters_variable.insert(0, self.node.maxiters_variable)
        self['maxiters_variable'] = maxiters_variable

        maxiters = ttk.Entry(frame, width=15)
        maxiters.insert(0, self.node.maxiters)
        self['maxiters'] = maxiters

        # maxevals: a variable or a value...
        maxevals_label = ttk.Label(frame, text='maximum energy evaluations:')
        self['maxevals_label'] = maxevals_label

        maxevals_method = ttk.Combobox(
            frame,
            state='readonly',
            values=[
                'default',
                'is',
                'from variable',
            ],
            justify=tk.LEFT,
            width=10
        )
        maxevals_method.set(self.node.maxevals_method)
        self['maxevals_method'] = maxevals_method

        maxevals_variable = ttk.Entry(frame, width=15)
        maxevals_variable.insert(0, self.node.maxevals_variable)
        self['maxevals_variable'] = maxevals_variable

        maxevals = ttk.Entry(frame, width=15)
        maxevals.insert(0, self.node.maxevals)
        self['maxevals'] = maxevals

        convergence.bind("<<ComboboxSelected>>", self.reset_dialog)
        etol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        ftol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxiters_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxevals_method.bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        convergence = self['convergence'].get()
        etol_method = self['etol_method'].get()
        ftol_method = self['ftol_method'].get()
        maxiters_method = self['maxiters_method'].get()
        maxevals_method = self['maxevals_method'].get()

        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        self['convergence_label'].grid(row=row, column=0, sticky=tk.E)
        self['convergence'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        if 'energy' in convergence:
            row += 1
            self['etol_label'].grid(row=row, column=0, sticky=tk.E)
            self['etol_method'].grid(row=row, column=1, sticky=tk.EW)
            if etol_method == 'from variable':
                self['etol_variable'].grid(row=row, column=2, sticky=tk.W)
            elif etol_method == 'is':
                self['etol'].grid(row=row, column=2, sticky=tk.W)

        if 'forces' in convergence:
            row += 1
            self['ftol_label'].grid(row=row, column=0, sticky=tk.E)
            self['ftol_method'].grid(row=row, column=1, sticky=tk.EW)
            if ftol_method == 'from variable':
                self['ftol_variable'].grid(row=row, column=2, sticky=tk.W)
            elif ftol_method == 'is':
                self['ftol'].grid(row=row, column=2, sticky=tk.W)

        if 'energy' in convergence or 'forces' in convergence:
            row += 1
            self['maxiters_label'].grid(row=row, column=0, sticky=tk.E)
            self['maxiters_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxiters_method == 'from variable':
                self['maxiters_variable'].grid(row=row, column=2, sticky=tk.W)
            elif maxiters_method == 'is':
                self['maxiters'].grid(row=row, column=2, sticky=tk.W)

            row += 1
            self['maxevals_label'].grid(row=row, column=0, sticky=tk.E)
            self['maxevals_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxevals_method == 'from variable':
                self['maxevals_variable'].grid(row=row, column=2, sticky=tk.W)
            elif maxevals_method == 'is':
                self['maxevals'].grid(row=row, column=2, sticky=tk.W)

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
                "Don't recognize dialog result '{}'".format(result)
            )

        self.dialog.deactivate(result)

        convergence = self['convergence'].get()
        etol_method = self['etol_method'].get()
        ftol_method = self['ftol_method'].get()
        maxiters_method = self['maxiters_method'].get()
        maxevals_method = self['maxevals_method'].get()

        self.node.convergence = convergence
        if 'energy' in convergence:
            self.node.etol_method = etol_method
            if etol_method == 'from variable':
                self.node.etol_variable = self['etol_variable'].get()
            elif etol_method == 'is':
                self.node.etol = int(self['etol'].get())

        if 'forces' in convergence:
            self.node.ftol_method = ftol_method
            if ftol_method == 'from variable':
                self.node.ftol_variable = self['ftol_variable'].get()
            elif ftol_method == 'is':
                self.node.ftol = self['ftol'].get()

        if 'energy' in convergence or 'forces' in convergence:
            self.node.maxiters_method = maxiters_method
            if maxiters_method == 'from variable':
                self.node.maxtiters_variable = self['maxiters_variable'].get()
            elif maxiters_method == 'is':
                self.node.maxiters = int(self['maxiters'].get())

            self.node.maxevals_method = maxevals_method
            if maxevals_method == 'from variable':
                self.node.maxevals_variable = self['maxevals_variable'].get()
            elif maxevals_method == 'is':
                self.node.maxevals = int(self['maxevals'].get())
