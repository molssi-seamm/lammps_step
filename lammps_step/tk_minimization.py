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

        w = self._widget
        frame = w['frame']

        # Convergence
        convergence_label = ttk.Label(frame, text='Simulation convergence:')
        w['convergence_label'] = convergence_label

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
        w['convergence'] = convergence

        # etol: a variable or a value...
        etol_label = ttk.Label(frame, text='change in energy:')
        w['etol_label'] = etol_label

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
        w['etol_method'] = etol_method

        etol_variable = ttk.Entry(frame, width=15)
        etol_variable.insert(0, self.node.etol_variable)
        w['etol_variable'] = etol_variable

        etol = ttk.Entry(frame, width=15)
        etol.insert(0, self.node.etol)
        w['etol'] = etol

        # ftol: a variable or a value...
        ftol_label = ttk.Label(frame, text='magnitude of the force vector:')
        w['ftol_label'] = ftol_label

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
        w['ftol_method'] = ftol_method

        ftol_variable = ttk.Entry(frame, width=15)
        ftol_variable.insert(0, self.node.ftol_variable)
        w['ftol_variable'] = ftol_variable

        ftol = sw.UnitEntry(frame, width=15)
        ftol.set(self.node.ftol)
        w['ftol'] = ftol

        # maxiters: a variable or a value...
        maxiters_label = ttk.Label(frame, text='maximum iterations:')
        w['maxiters_label'] = maxiters_label

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
        w['maxiters_method'] = maxiters_method

        maxiters_variable = ttk.Entry(frame, width=15)
        maxiters_variable.insert(0, self.node.maxiters_variable)
        w['maxiters_variable'] = maxiters_variable

        maxiters = ttk.Entry(frame, width=15)
        maxiters.insert(0, self.node.maxiters)
        w['maxiters'] = maxiters

        # maxevals: a variable or a value...
        maxevals_label = ttk.Label(frame, text='maximum energy evaluations:')
        w['maxevals_label'] = maxevals_label

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
        w['maxevals_method'] = maxevals_method

        maxevals_variable = ttk.Entry(frame, width=15)
        maxevals_variable.insert(0, self.node.maxevals_variable)
        w['maxevals_variable'] = maxevals_variable

        maxevals = ttk.Entry(frame, width=15)
        maxevals.insert(0, self.node.maxevals)
        w['maxevals'] = maxevals

        convergence.bind("<<ComboboxSelected>>", self.reset_dialog)
        etol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        ftol_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxiters_method.bind("<<ComboboxSelected>>", self.reset_dialog)
        maxevals_method.bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        w = self._widget
        convergence = w['convergence'].get()
        etol_method = w['etol_method'].get()
        ftol_method = w['ftol_method'].get()
        maxiters_method = w['maxiters_method'].get()
        maxevals_method = w['maxevals_method'].get()

        frame = w['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['convergence_label'].grid(row=row, column=0, sticky=tk.E)
        w['convergence'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        if 'energy' in convergence:
            row += 1
            w['etol_label'].grid(row=row, column=0, sticky=tk.E)
            w['etol_method'].grid(row=row, column=1, sticky=tk.EW)
            if etol_method == 'from variable':
                w['etol_variable'].grid(row=row, column=2, sticky=tk.W)
            elif etol_method == 'is':
                w['etol'].grid(row=row, column=2, sticky=tk.W)

        if 'forces' in convergence:
            row += 1
            w['ftol_label'].grid(row=row, column=0, sticky=tk.E)
            w['ftol_method'].grid(row=row, column=1, sticky=tk.EW)
            if ftol_method == 'from variable':
                w['ftol_variable'].grid(row=row, column=2, sticky=tk.W)
            elif ftol_method == 'is':
                w['ftol'].grid(row=row, column=2, sticky=tk.W)

        if 'energy' in convergence or 'forces' in convergence:
            row += 1
            w['maxiters_label'].grid(row=row, column=0, sticky=tk.E)
            w['maxiters_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxiters_method == 'from variable':
                w['maxiters_variable'].grid(row=row, column=2, sticky=tk.W)
            elif maxiters_method == 'is':
                w['maxiters'].grid(row=row, column=2, sticky=tk.W)

            row += 1
            w['maxevals_label'].grid(row=row, column=0, sticky=tk.E)
            w['maxevals_method'].grid(row=row, column=1, sticky=tk.EW)
            if maxevals_method == 'from variable':
                w['maxevals_variable'].grid(row=row, column=2, sticky=tk.W)
            elif maxevals_method == 'is':
                w['maxevals'].grid(row=row, column=2, sticky=tk.W)

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

        w = self._widget
        convergence = w['convergence'].get()
        etol_method = w['etol_method'].get()
        ftol_method = w['ftol_method'].get()
        maxiters_method = w['maxiters_method'].get()
        maxevals_method = w['maxevals_method'].get()

        self.node.convergence = convergence
        if 'energy' in convergence:
            self.node.etol_method = etol_method
            if etol_method == 'from variable':
                self.node.etol_variable = w['etol_variable'].get()
            elif etol_method == 'is':
                self.node.etol = int(w['etol'].get())

        if 'forces' in convergence:
            self.node.ftol_method = ftol_method
            if ftol_method == 'from variable':
                self.node.ftol_variable = w['ftol_variable'].get()
            elif ftol_method == 'is':
                self.node.ftol = w['ftol'].get()

        if 'energy' in convergence or 'forces' in convergence:
            self.node.maxiters_method = maxiters_method
            if maxiters_method == 'from variable':
                self.node.maxtiters_variable = w['maxiters_variable'].get()
            elif maxiters_method == 'is':
                self.node.maxiters = int(w['maxiters'].get())

            self.node.maxevals_method = maxevals_method
            if maxevals_method == 'from variable':
                self.node.maxevals_variable = w['maxevals_variable'].get()
            elif maxevals_method == 'is':
                self.node.maxevals = int(w['maxevals'].get())
