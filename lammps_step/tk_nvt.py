# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS NVT dynamics step"""

import lammps_step
import logging
import seamm_widgets as sw
import pprint  # nopep8
import tkinter as tk
import tkinter.ttk as ttk

logger = logging.getLogger(__name__)


class TkNVT(lammps_step.TkNVE):

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

        logger.debug('TkNVT.create_dialog')

        # Let parent classes do their thing.
        super().create_dialog()

        self.dialog.configure(title='Edit NVT dynamics parameters')

        # Shortcut for parameters
        P = self.node.parameters

        logger.debug('Parameters:\n{}'.format(pprint.pformat(P.to_dict())))

        # Temperature frame to isolate widgets
        t_frame = self['temperature_frame'] = ttk.LabelFrame(
            self['frame'],
            borderwidth=4,
            relief='sunken',
            text='Temperature',
            labelanchor='n',
            padding=10
        )

        for key in lammps_step.NVT_Parameters.parameters:
            self[key] = P[key].widget(t_frame)

        # and binding to change as needed
        self['thermostat'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_temperature_frame
        )

        self.setup_results('nvt')

    def reset_temperature_frame(self, widget=None):
        """Layout the widgets in the temperature frame
        as needed for the current state"""

        thermostat = self['thermostat'].get()

        t_frame = self['temperature_frame']
        for slave in t_frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # Main controls
        self['T0'].grid(row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1

        self['T1'].grid(row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1

        self['thermostat'].grid(row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1

        sw.align_labels((self['T0'], self['T1'], self['thermostat']))

        # and controls for specific thermostats
        if thermostat != 'velocity rescaling':
            self['Tdamp'].grid(row=row, column=1, sticky=tk.W)
            row += 1

        if thermostat == 'Nose-Hoover':
            self['drag'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            self['Tchain'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            self['Tloop'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            sw.align_labels(
                (self['Tdamp'], self['Tchain'], self['Tloop'], self['drag'])
            )
        elif thermostat == 'Berendsen':
            pass
        elif 'csvr' in thermostat or 'csld' in thermostat:
            self['seed'].grid(row=row, column=1, sticky=tk.W)
            row += 1
        elif thermostat == 'velocity rescaling':
            self['frequency'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            self['window'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            self['fraction'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            sw.align_labels(
                (self['frequency'], self['window'], self['fraction'])
            )
        elif thermostat == 'Langevin':
            self['seed'].grid(row=row, column=1, sticky=tk.W)
            row += 1
            sw.align_labels((self['Tdamp'], self['seed']))
        else:
            raise RuntimeError(
                "Don't recognize thermostat " + "'{}'".format(thermostat)
            )

        t_frame.columnconfigure(0, minsize=50)

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        row = super().reset_dialog()

        self['temperature_frame'].grid(
            row=row, column=0, sticky=tk.EW, pady=10
        )
        row += 1
        self.reset_temperature_frame()

        return row

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

        # Let base classes reap their parameters
        super().handle_dialog(result)

        # Shortcut for parameters
        P = self.node.parameters

        thermostat = self['thermostat'].get()

        P['T0'].set(self['T0'].get())
        P['T1'].set(self['T1'].get())
        P['thermostat'].set(thermostat)

        if thermostat != 'velocity rescaling':
            P['Tdamp'].set(self['Tdamp'].get())

        if thermostat == 'Nose-Hoover':
            P['Tchain'].set(self['Tchain'].get())
            P['Tloop'].set(self['Tloop'].get())
            P['drag'].set(self['drag'].get())
        elif 'csvr' in thermostat or 'csld' in thermostat or \
             thermostat == 'Langevin':
            P['seed'].set(self['seed'].get())
        elif thermostat == 'velocity rescaling':
            P['frequency'].set(self['frequency'].get())
            P['window'].set(self['window'].get())
            P['fraction'].set(self['fraction'].get())
