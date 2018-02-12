# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Energy step"""

import lammps_step
from molssi_workflow import units, Q_, units_class  # nopep8
import molssi_util.molssi_widgets as mw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk


class TkNVT(lammps_step.TkNVE):
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

        frame = self._widget['frame']
        self.dialog.configure(title='Edit NVT dynamics step')
        w = self._widget

        # Simulation time...
        time_label = ttk.Label(frame, text='Simulation time:')
        w['time_label'] = time_label

        time = mw.UnitEntry(frame, width=15)
        time.set(self.node.time)
        w['time'] = time

        # Timestep: automatic, a variable or a value...
        timestep_label = ttk.Label(frame, text='Timestep:')
        w['timestep_label'] = timestep_label

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
        w['timestep_method'] = timestep_method

        timestep_variable = ttk.Entry(frame, width=15)
        timestep_variable.insert(0, self.node.timestep_variable)
        w['timestep_variable'] = timestep_variable

        timestep = mw.UnitEntry(frame, width=15)
        timestep.set(self.node.timestep)
        w['timestep'] = timestep

        # Temperature control
        tc_frame = w['Tcontrol_frame'] = ttk.Frame(
            frame, borderwidth=4, relief='sunken')

        Tcontrol_method_label = ttk.Label(
            tc_frame, text='Temperature control:'
        )
        w['Tcontrol_method_label'] = Tcontrol_method_label

        methods = lammps_step.NVT.methods
        width = 0
        for method in list(methods):
            width = len(method) if len(method) > width else width
        width += 3
        Tcontrol_method = ttk.Combobox(
            tc_frame, state='readonly',
            values=list(methods),
            justify=tk.LEFT, width=width
        )
        Tcontrol_method.set(self.node.Tcontrol_method)
        w['Tcontrol_method'] = Tcontrol_method

        # Temperature control Parameters
        self.create_unit_entry(frame, 'T0', text='Initial temperature:')
        self.create_unit_entry(frame, 'T1', text='Final temperature:')
        self.create_unit_entry(tc_frame, 'Tdamp', text='Damping time:')
        self.create_unit_entry(tc_frame, 'Tchain', text='Thermostat chain:')
        self.create_unit_entry(
            tc_frame, 'Tloop', text='# of thermostat loops:')
        self.create_unit_entry(tc_frame, 'drag', text='Drag:')
        self.create_unit_entry(tc_frame, 'seed', text='Random seed:',
                               methods=['is', 'from variable', 'random'])
        self.create_unit_entry(
            tc_frame, 'frequency', text='Rescaling frequency:')
        self.create_unit_entry(tc_frame, 'window', text='Rescaling window:')
        self.create_unit_entry(tc_frame, 'fraction',
                               text='Fraction to rescale:')

        # Sampling
        self.create_unit_entry(frame, 'sampling', text='Sampling frequency:',
                               methods=['is', 'from variable', 'none'])

        w['timestep_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['Tcontrol_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['T0_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['T1_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['Tdamp_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['Tchain_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['Tloop_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['drag_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['seed_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['frequency_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['window_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['fraction_method'].bind("<<ComboboxSelected>>", self.reset_dialog)
        w['sampling_method'].bind("<<ComboboxSelected>>", self.reset_dialog)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""
        w = self._widget
        timestep_method = w['timestep_method'].get()
        Tcontrol_method = w['Tcontrol_method'].get()
        T0_method = w['T0_method'].get()
        T1_method = w['T1_method'].get()
        Tdamp_method = w['Tdamp_method'].get()
        Tchain_method = w['Tchain_method'].get()
        Tloop_method = w['Tloop_method'].get()
        drag_method = w['drag_method'].get()
        seed_method = w['seed_method'].get()
        frequency_method = w['frequency_method'].get()
        window_method = w['window_method'].get()
        fraction_method = w['fraction_method'].get()
        sampling_method = w['sampling_method'].get()

        frame = self._widget['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()
        tc_frame = w['Tcontrol_frame']
        for slave in tc_frame.grid_slaves():
            slave.grid_forget()

        row = 0
        w['time_label'].grid(row=row, column=0, sticky=tk.E)
        w['time'].grid(row=row, column=1, columnspan=2, sticky=tk.W)

        row += 1
        w['timestep_label'].grid(row=row, column=0, sticky=tk.E)
        w['timestep_method'].grid(row=row, column=1, sticky=tk.EW)
        if timestep_method == 'from variable':
            w['timestep_variable'].grid(row=row, column=2, sticky=tk.W)
        elif timestep_method == 'is':
            w['timestep'].grid(row=row, column=2, sticky=tk.W)

        row += 1
        w['T0_label'].grid(row=row, column=0, sticky=tk.E)
        w['T0_method'].grid(row=row, column=1, sticky=tk.EW)
        if T0_method == 'is':
            w['T0'].grid(row=row, column=2, sticky=tk.W)
        elif T0_method == 'from variable':
            w['T0_variable'].grid(row=row, column=2, sticky=tk.W)

        row += 1
        w['T1_label'].grid(row=row, column=0, sticky=tk.E)
        w['T1_method'].grid(row=row, column=1, sticky=tk.EW)
        if T1_method == 'is':
            w['T1'].grid(row=row, column=2, sticky=tk.W)
        elif T1_method == 'from variable':
            w['T1_variable'].grid(row=row, column=2, sticky=tk.W)

        row += 1
        w['sampling_label'].grid(row=row, column=0, sticky=tk.E)
        w['sampling_method'].grid(row=row, column=1, sticky=tk.EW)
        if sampling_method == 'from variable':
            w['sampling_variable'].grid(row=row, column=2, sticky=tk.W)
        elif sampling_method == 'is':
            w['sampling'].grid(row=row, column=2, sticky=tk.W)

        frame.columnconfigure(2, weight=1)

        row += 1
        tc_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW)

        tc_row = 0
        w['Tcontrol_method_label'].grid(row=tc_row, column=0, columnspan=2,
                                        sticky=tk.W)
        w['Tcontrol_method'].grid(row=tc_row, column=2, columnspan=3,
                                  sticky=tk.W)

        tc_row += 1
        if Tcontrol_method != 'velocity rescaling':
            tc_row += 1
            w['Tdamp_label'].grid(row=tc_row, column=1, columnspan=2,
                                  sticky=tk.E)
            w['Tdamp_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if Tdamp_method == 'is':
                w['Tdamp'].grid(row=tc_row, column=4, sticky=tk.W)
            elif Tdamp_method == 'from variable':
                w['Tdamp_variable'].grid(row=tc_row, column=4, sticky=tk.W)

        if Tcontrol_method == 'Nose-Hoover':
            tc_row += 1
            w['Tchain_label'].grid(row=tc_row, column=1, columnspan=2,
                                   sticky=tk.E)
            w['Tchain_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if Tchain_method == 'is':
                w['Tchain'].grid(row=tc_row, column=4, sticky=tk.W)
            elif Tchain_method == 'from variable':
                w['Tchain_variable'].grid(row=tc_row, column=4, sticky=tk.W)

            tc_row += 1
            w['Tloop_label'].grid(row=tc_row, column=1, columnspan=2,
                                  sticky=tk.E)
            w['Tloop_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if Tloop_method == 'is':
                w['Tloop'].grid(row=tc_row, column=4, sticky=tk.W)
            elif Tloop_method == 'from variable':
                w['Tloop_variable'].grid(row=tc_row, column=4, sticky=tk.W)

            tc_row += 1
            w['drag_label'].grid(row=tc_row, column=1, columnspan=2,
                                 sticky=tk.E)
            w['drag_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if drag_method == 'is':
                w['drag'].grid(row=tc_row, column=4, sticky=tk.W)
            elif drag_method == 'from variable':
                w['drag_variable'].grid(row=tc_row, column=4, sticky=tk.W)
        elif Tcontrol_method == 'Berendsen':
            pass
        elif 'csvr' in Tcontrol_method or 'csld' in Tcontrol_method:
            tc_row += 1
            w['seed_label'].grid(row=tc_row, column=1, columnspan=2,
                                 sticky=tk.E)
            w['seed_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if seed_method == 'random':
                pass
            elif seed_method == 'is':
                w['seed'].grid(row=tc_row, column=4, sticky=tk.W)
            elif seed_method == 'from variable':
                w['seed_variable'].grid(row=tc_row, column=4, sticky=tk.W)
        elif Tcontrol_method == 'velocity rescaling':
            tc_row += 1
            w['frequency_label'].grid(row=tc_row, column=1, columnspan=2,
                                      sticky=tk.E)
            w['frequency_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if frequency_method == 'is':
                w['frequency'].grid(row=tc_row, column=4, sticky=tk.W)
            elif frequency_method == 'from variable':
                w['frequency_variable'].grid(row=tc_row, column=4, sticky=tk.W)

            tc_row += 1
            w['window_label'].grid(row=tc_row, column=1, columnspan=2,
                                   sticky=tk.E)
            w['window_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if window_method == 'is':
                w['window'].grid(row=tc_row, column=4, sticky=tk.W)
            elif window_method == 'from variable':
                w['window_variable'].grid(row=tc_row, column=4, sticky=tk.W)

            tc_row += 1
            w['fraction_label'].grid(row=tc_row, column=1, columnspan=2,
                                     sticky=tk.E)
            w['fraction_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if fraction_method == 'is':
                w['fraction'].grid(row=tc_row, column=4, sticky=tk.W)
            elif fraction_method == 'from variable':
                w['fraction_variable'].grid(row=tc_row, column=4, sticky=tk.W)
        elif Tcontrol_method == 'Langevin':
            tc_row += 1
            w['seed_label'].grid(row=tc_row, column=1, columnspan=2,
                                 sticky=tk.E)
            w['seed_method'].grid(row=tc_row, column=3, sticky=tk.EW)
            if seed_method == 'random':
                pass
            elif seed_method == 'is':
                w['seed'].grid(row=tc_row, column=4, sticky=tk.W)
            elif seed_method == 'from variable':
                w['seed_variable'].grid(row=tc_row, column=4, sticky=tk.W)
        else:
            raise RuntimeError("Don't recognize temperature control " +
                               "'{}'".format(Tcontrol_method))

        tc_frame.columnconfigure(0, minsize=50)

    def handle_dialog(self, result):
        super().handle_dialog(result)

        if result is None or result == 'Cancel':
            return

        if result == 'Help':
            # display help!!!
            return

        if result != "OK":
            raise RuntimeError(
                "Don't recognize dialog result '{}'".format(result))

        w = self._widget
        Tcontrol_method = w['Tcontrol_method'].get()
        T0_method = w['T0_method'].get()
        T1_method = w['T1_method'].get()
        Tdamp_method = w['Tdamp_method'].get()
        Tchain_method = w['Tchain_method'].get()
        Tloop_method = w['Tloop_method'].get()
        drag_method = w['drag_method'].get()
        seed_method = w['seed_method'].get()
        frequency_method = w['frequency_method'].get()
        window_method = w['window_method'].get()
        fraction_method = w['fraction_method'].get()

        self.node.T0_method = T0_method
        if T0_method == 'is':
            self.node.T0 = w['T0'].get()
        elif T0_method == 'from variable':
            self.node.T0_variable = w['T0_variable'].get()

        self.node.T1_method = T1_method
        if T1_method == 'is':
            self.node.T1 = w['T1'].get()
        elif T1_method == 'from variable':
            self.node.T1_variable = w['T1_variable'].get()

        self.node.Tcontrol_method = Tcontrol_method

        if Tcontrol_method != 'velocity rescaling':
            self.node.Tdamp_method = Tdamp_method
            if Tdamp_method == 'is':
                self.node.Tdamp = w['Tdamp'].get()
            elif Tdamp_method == 'from variable':
                self.node.Tdamp_variable = w['Tdamp_variable'].get()

        if Tcontrol_method == 'Nose-Hoover':
            self.node.Tchain_method = Tchain_method
            if Tchain_method == 'is':
                self.node.Tchain = w['Tchain'].get()
            elif Tchain_method == 'from variable':
                self.node.Tchain_variable = w['Tchain_variable'].get()

            self.node.Tloop_method = Tloop_method
            if Tloop_method == 'is':
                self.node.Tloop = w['Tloop'].get()
            elif Tloop_method == 'from variable':
                self.node.Tloop_variable = w['Tloop_variable'].get()

            self.node.drag_method = drag_method
            if drag_method == 'is':
                self.node.drag = w['drag'].get()
            elif drag_method == 'from variable':
                self.node.drag_variable = w['drag_variable'].get()
            elif 'csvr' in Tcontrol_method or 'csld' in Tcontrol_method or \
                 Tcontrol_method == 'Langevin':
                self.node.seed_method = seed_method
                if seed_method == 'is':
                    self.node.seed = w['seed'].get()
                elif seed_method == 'from variable':
                    self.node.seed_variable = w['seed_variable'].get()
            elif Tcontrol_method == 'velocity rescaling':
                self.node.frequency_method = frequency_method
                if frequency_method == 'is':
                    self.node.frequency = w['frequency'].get()
                elif frequency_method == 'from variable':
                    self.node.frequency_variable = \
                        w['frequency_variable'].get()

                self.node.window_method = window_method
                if window_method == 'is':
                    self.node.window = w['window'].get()
                elif window_method == 'from variable':
                    self.node.window_variable = w['window_variable'].get()

                self.node.fraction_method = fraction_method
                if fraction_method == 'is':
                    self.node.fraction = w['fraction'].get()
                elif fraction_method == 'from variable':
                    self.node.fraction_variable = w['fraction_variable'].get()

    def create_unit_entry(self,
                          parent,
                          name,
                          text=None,
                          methods=['is', 'from variable'],
                          ):
        # the label
        widget = name + '_label'
        self._widget[widget] = ttk.Label(parent, text=text)

        # the combobox for how to set the variable
        widget = name + '_method'
        self._widget[widget] = ttk.Combobox(
            parent, state='readonly',
            values=methods,
            justify=tk.LEFT, width=15
        )
        self._widget[widget].set(self.node.__dict__[widget])

        # the entry or unitentry
        widget = name
        value = self.node.__dict__[widget]
        if isinstance(value, units_class):
            self._widget[widget] = mw.UnitEntry(parent, width=15)
            self._widget[widget].set(value)
        else:
            self._widget[widget] = ttk.Entry(parent, width=15)
            self._widget[widget].insert(0, value)

        # an entry to specify the variable, if needed
        widget = name + '_variable'
        value = self.node.__dict__[widget]
        self._widget[widget] = ttk.Entry(parent, width=15)
        self._widget[widget].insert(0, value)
